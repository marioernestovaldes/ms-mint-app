from dash import html, dash_table
from ms_mint.io import mzxml_to_df
from dash_tabulator import DashTabulator

test_file = ('/media/mario/92bb1dc0-1f5c-4566-bfea-709a64de54f8/mario/Workspace/Proteomics_data/LSARP_data'
             '/proteomics/Papers/ms-MINT/SRMDataForSoren/2024_10_16_Chr_Med_Rep1.mzXML')  # replace with real


class MS2BrowserPlugin:
    label = "MS2 Browser"
    order = 999

    @staticmethod
    def layout(**kwargs):

        try:
            df = mzxml_to_df(test_file)
        except Exception as e:
            return html.Div([html.H3("Error loading MS2 file"), html.Pre(str(e))])

        # Show table
        columns = [
            {"title": col, "field": col, "headerFilter": True}
            for col in df.columns
        ]

        return html.Div([
            html.Div([
                html.Button("EXPORT", id="btn-export-ms2", n_clicks=0, className="btn btn-primary"),
                html.Button("CLEAR FILTERS", id="btn-clear-ms2", n_clicks=0, className="btn btn-secondary"),
            ], style={"marginBottom": "10px", "display": "flex", "gap": "10px"}),

            DashTabulator(
                id="tabulator-ms2-browser",
                columns=columns,
                data=df.to_dict("records"),
                options={
                    "layout": "fitDataStretch",
                    "pagination": "local",
                    "paginationSize": 10,
                    "movableColumns": True,
                    "resizableColumns": True,
                    "selectable": 1,
                    "clipboard": True,
                    "clipboardPasteAction": "replace",
                },
            ),
        ])

    @staticmethod
    def outputs():
        # Not used currently, but required for layout registry
        return None

    @staticmethod
    def callbacks(app, fsc, cache):
        # No interactivity yet — placeholder required
        pass


