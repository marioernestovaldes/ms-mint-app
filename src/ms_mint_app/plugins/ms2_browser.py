from dash import html, dcc
from dash_tabulator import DashTabulator
from dash.dependencies import Input, Output
import plotly.graph_objects as go
import numpy as np
import pandas as pd
from ms_mint.io import mzxml_to_df

# Replace with your actual test file
test_file = '/home/mario/Workspace/LSARP/SRMDataForSoren/2024_10_16_Chr_STD12_Rep1.mzXML'


class MS2BrowserPlugin:
    label = "MS2 Browser"
    order = 999
    df_static = None

    @staticmethod
    def layout(**kwargs):
        try:
            df = mzxml_to_df(test_file)
        except Exception as e:
            return html.Div([html.H3("Error loading MS2 file"), html.Pre(str(e))])

        MS2BrowserPlugin.df_static = df

        # Table
        columns = [{"title": col, "field": col, "headerFilter": True} for col in df.columns]
        table = DashTabulator(
            id="tabulator-ms2-browser",
            columns=columns,
            data=df.to_dict("records"),
            options={
                "layout": "fitDataStretch",
                "pagination": "local",
                "paginationSize": 10,
                "movableColumns": True,
                "resizableColumns": True,
            },
        )

        # Dropdown of unique channels (either in "filterLine" or "filterLine_to_ELMAVEN" columns)
        cols = [c for c in ["filterLine", "filterLine_to_ELMAVEN"] if c in df.columns]
        values = df[cols].values.flatten().tolist()
        clean_values = [v for v in values if isinstance(v, str) and pd.notna(v)]
        channel_options = sorted(list(np.unique(clean_values)))

        dropdown = dcc.Dropdown(
            id="channel-selector",
            options=list(channel_options),
            placeholder="Select a channel...",
            style={"width": "400px", "marginBottom": "10px"},
        )

        return html.Div([
            html.Div([
                html.Button("EXPORT", id="btn-export-ms2", n_clicks=0, className="btn btn-primary"),
                html.Button("CLEAR FILTERS", id="btn-clear-ms2", n_clicks=0, className="btn btn-secondary"),
            ], style={"marginBottom": "10px", "display": "flex", "gap": "10px"}),

            table,

            html.Div([
                html.H4("Select Channel"),
                dropdown
            ], style={"marginTop": "20px"}),

            html.Div(id="channel-fragment-plot-container", style={"marginTop": "30px"}),
        ])

    @staticmethod
    def create_channel_timeline_plot(df, channel):
        if "ESI SRM ms2" in channel:
            mask = df["filterLine"] == channel
        else:
            mask = df["filterLine_to_ELMAVEN"] == channel

        filtered = df[mask & (df["intensity"] != 0)].copy()

        if filtered.empty:
            return go.Figure().update_layout(title=f"No MS2 spectra for channel {channel}")

        expanded = []
        for _, row in filtered.iterrows():
            scan_time = row["scan_time"]

            if "ESI SRM ms2" in channel:
                channel = row["filterLine"]
            else:
                channel = row["filterLine_to_ELMAVEN"]

            expanded.extend(
                {
                    "scan_time": scan_time,
                    "channel": channel,
                    "mz": mz_val,
                    "intensity": intensity_val,
                }
                for mz_val, intensity_val in zip(row["mz"], row["intensity"])
            )
        df_expanded = pd.DataFrame(expanded)

        fig = go.Figure()
        for _, row in df_expanded.iterrows():
            fig.add_trace(go.Scatter(
                x=[row["scan_time"], row["scan_time"]],
                y=[0, row["intensity"]],
                mode="lines",
                line=dict(width=1, color="black"),
                hoverinfo="text",
                hovertext=f"Scan time: {row['scan_time']}<br>"
                          f"m/z: {row['mz']:.4f}<br>"
                          f"Intensity: {row['intensity']:.1f}<br>"
                          f"Channel: {row['channel']}"
            ))

        fig.update_layout(
            title=f"Fragmentation Pattern Over Time for Channel {channel}",
            xaxis_title="Scan Time",
            yaxis_title="Fragment Intensity",
            yaxis_title_standoff=20,
            showlegend=False,
            height=350,
            margin=dict(l=60, r=20, t=30, b=40),
        )

        return fig

    @staticmethod
    def outputs():
        return None

    @staticmethod
    def callbacks(app, fsc, cache):
        @app.callback(
            Output("channel-fragment-plot-container", "children"),
            Input("channel-selector", "value"),
            prevent_initial_call=True,
        )
        def update_channel_plot(channel):
            if channel is None:
                return dcc.Markdown("Select a channel to view fragmentation over time.")

            df = MS2BrowserPlugin.df_static.copy()
            df["mz"] = df["mz"].apply(lambda x: x if isinstance(x, list) else [x])
            df["intensity"] = df["intensity"].apply(lambda x: x if isinstance(x, list) else [x])

            fig = MS2BrowserPlugin.create_channel_timeline_plot(df, channel)
            return dcc.Graph(figure=fig)
