from dash import html, dcc
from dash_tabulator import DashTabulator
from dash.dependencies import Input, Output
import plotly.graph_objects as go
import numpy as np
import pandas as pd
from ms_mint.io import mzxml_to_df

# Replace with your actual test file
test_file = ('/media/mario/92bb1dc0-1f5c-4566-bfea-709a64de54f8/mario/Workspace/Proteomics_data/LSARP_data/proteomics'
             '/Papers/ms-MINT/SRMDataForSoren/2024_10_16_Chr_Med_Rep1.mzXML')


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

        # Dropdown of unique precursor m/z values
        precursor_options = sorted(df["mz_precursor"].dropna().unique())
        dropdown = dcc.Dropdown(
            id="precursor-selector",
            options=[{"label": f"{mz:.4f}", "value": mz} for mz in precursor_options],
            placeholder="Select a precursor m/z...",
            style={"width": "300px", "marginBottom": "10px"},
        )

        return html.Div([
            html.Div([
                html.Button("EXPORT", id="btn-export-ms2", n_clicks=0, className="btn btn-primary"),
                html.Button("CLEAR FILTERS", id="btn-clear-ms2", n_clicks=0, className="btn btn-secondary"),
            ], style={"marginBottom": "10px", "display": "flex", "gap": "10px"}),

            table,

            html.Div([
                html.H4("Select Precursor m/z"),
                dropdown
            ], style={"marginTop": "20px"}),

            html.Div([
                html.H4("MS2 Fragmentation Spectrum"),
                dcc.Graph(figure=MS2BrowserPlugin.create_ms2_spectrum_plot(df))
            ], style={"marginTop": "20px"}),

            html.Div(id="precursor-fragment-plot-container", style={"marginTop": "30px"}),
        ])

    @staticmethod
    def create_ms2_spectrum_plot(df):
        filtered = df[df["intensity"] != 0][["mz", "intensity"]]
        grouped = filtered.groupby("mz", as_index=False).sum()

        fig = go.Figure()
        for _, row in grouped.iterrows():
            fig.add_trace(go.Scatter(
                x=[row["mz"], row["mz"]],
                y=[0, row["intensity"]],
                mode="lines",
                line=dict(width=1, color="black"),
                hoverinfo="text",
                hovertext=f"m/z: {row['mz']:.4f}<br>Summed Intensity: {row['intensity']:.1f}"
            ))

        fig.update_layout(
            title="MS2 Spectrum",
            xaxis_title="m/z",
            yaxis_title="Intensity",
            yaxis_title_standoff=20,
            showlegend=False,
            height=300,
            margin=dict(l=60, r=20, t=30, b=40),
        )

        return fig

    @staticmethod
    def create_precursor_timeline_plot(df, mz_precursor, atol=0.001):
        mask = np.isclose(df["mz_precursor"], mz_precursor, atol=atol)
        filtered = df[mask & (df["intensity"] != 0)].copy()

        if filtered.empty:
            return go.Figure().update_layout(title=f"No MS2 spectra for m/z ≈ {mz_precursor:.4f}")

        expanded = []
        for _, row in filtered.iterrows():
            scan_time = row["scan_time"]
            precursor = row["mz_precursor"]
            for mz_val, intensity_val in zip(row["mz"], row["intensity"]):
                expanded.append({
                    "scan_time": scan_time,
                    "mz_precursor": precursor,
                    "mz": mz_val,
                    "intensity": intensity_val
                })

        df_expanded = pd.DataFrame(expanded)

        fig = go.Figure()
        for _, row in df_expanded.iterrows():
            fig.add_trace(go.Scatter(
                x=[row["scan_time"], row["scan_time"]],
                y=[0, row["intensity"]],
                mode="lines",
                line=dict(width=1, color="black"),
                hoverinfo="text",
                hovertext=f"Scan time: {row['scan_time']}<br>m/z: {row['mz']:.4f}<br>Intensity: {row['intensity']:.1f}"
            ))

        fig.update_layout(
            title=f"Fragmentation Pattern Over Time for Precursor m/z ≈ {mz_precursor:.4f}",
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
            Output("precursor-fragment-plot-container", "children"),
            Input("precursor-selector", "value"),
            prevent_initial_call=True,
        )
        def update_precursor_plot(mz_precursor):
            if mz_precursor is None:
                return dcc.Markdown("Select a precursor m/z to view fragmentation over time.")

            df = MS2BrowserPlugin.df_static.copy()
            df["mz"] = df["mz"].apply(lambda x: x if isinstance(x, list) else [x])
            df["intensity"] = df["intensity"].apply(lambda x: x if isinstance(x, list) else [x])

            fig = MS2BrowserPlugin.create_precursor_timeline_plot(df, float(mz_precursor))
            return dcc.Graph(figure=fig)
