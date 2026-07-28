from __future__ import annotations

import base64
import json
import os
import re
import subprocess
import sys
import threading
import traceback
import uuid
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple
from urllib.parse import parse_qs

from dash import Dash, Input, Output, State, dcc, html, no_update

APP_DIR = Path(__file__).resolve().parent
PIPELINE_SCRIPT = APP_DIR / "humanize.py"
ASSETS_DIR = APP_DIR / "assets"

RUNS_DIR = Path(os.getenv("HUMANIZER_RUNS_DIR", str(APP_DIR / "app_runs"))).expanduser()

BANNER_PATH = ASSETS_DIR / "app_banner.png"
IMAGE_PATH = APP_DIR / "banner.png"
DEFAULT_DATABASE = Path(
    os.getenv(
        "HUMANIZER_DB_PATH",
        str(APP_DIR / "imgt_germline_database" / "human_germline_frameworks"),
    )
)
STATUS_POLL_MS = 1500

RUNS_DIR.mkdir(parents=True, exist_ok=True)
ASSETS_DIR.mkdir(parents=True, exist_ok=True)


app = Dash(__name__, suppress_callback_exceptions=True, title="Antibody Humanizer")
server = app.server

app.index_string = """
<!DOCTYPE html>
<html>
    <head>
        {%metas%}
        <title>{%title%}</title>
        {%favicon%}
        {%css%}
        <style>
            :root {
                color-scheme: light;
                font-synthesis: none;
                text-rendering: optimizeLegibility;
            }

            * {
                box-sizing: border-box;
            }

            html {
                background: #eef3f8;
            }

            body {
                margin: 0;
                min-width: 320px;
                background: #eef3f8;
                color: #0f172a;
                font-family: Inter, ui-sans-serif, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
                -webkit-font-smoothing: antialiased;
                -moz-osx-font-smoothing: grayscale;
            }

            button, input, textarea {
                font: inherit;
            }

            ::selection {
                background: #bfdbfe;
                color: #0f172a;
            }

            .content-card {
                transition: box-shadow 180ms ease, border-color 180ms ease;
            }

            .content-card:hover {
                border-color: #d7e0eb !important;
                box-shadow: 0 10px 32px rgba(15, 23, 42, 0.07) !important;
            }

            .sequence-textarea,
            .standard-input {
                transition: border-color 160ms ease, box-shadow 160ms ease, background-color 160ms ease;
            }

            .sequence-textarea:hover,
            .standard-input:hover {
                border-color: #9fb1c6 !important;
            }

            .sequence-textarea:focus,
            .standard-input:focus {
                border-color: #3b82f6 !important;
                background-color: #ffffff !important;
                box-shadow: 0 0 0 4px rgba(59, 130, 246, 0.12) !important;
                outline: none !important;
            }

            .upload-zone {
                transition: border-color 160ms ease, background-color 160ms ease, transform 160ms ease;
            }

            .upload-zone:hover {
                border-color: #4f8fd8 !important;
                background-color: #eff7ff !important;
                transform: translateY(-1px);
            }

            .choice-group input[type="radio"] {
                accent-color: #2563eb;
            }

            .choice-group label {
                transition: background-color 140ms ease, color 140ms ease;
            }

            .choice-group label:hover {
                background: #edf5ff;
                color: #123d73 !important;
            }

            .help-badge,
            .mini-help-badge {
                transition: transform 140ms ease, background-color 140ms ease;
            }

            .help-badge:hover,
            .mini-help-badge:hover {
                transform: scale(1.06);
            }

            .primary-run-button {
                transition: transform 150ms ease, box-shadow 150ms ease, filter 150ms ease;
            }

            .primary-run-button:hover {
                transform: translateY(-1px);
                box-shadow: 0 11px 25px rgba(37, 99, 235, 0.28) !important;
                filter: saturate(1.06);
            }

            .primary-run-button:active {
                transform: translateY(0);
                box-shadow: 0 6px 14px rgba(37, 99, 235, 0.22) !important;
            }

            .primary-run-button:focus-visible,
            a:focus-visible {
                outline: 3px solid rgba(59, 130, 246, 0.35);
                outline-offset: 3px;
            }

            .feedback-box:empty {
                display: none;
            }

            .file-preview {
                scrollbar-color: #64748b #111827;
                scrollbar-width: thin;
            }

            .file-preview::-webkit-scrollbar {
                width: 10px;
                height: 10px;
            }

            .file-preview::-webkit-scrollbar-track {
                background: #111827;
            }

            .file-preview::-webkit-scrollbar-thumb {
                background: #475569;
                border: 2px solid #111827;
                border-radius: 999px;
            }

            .result-table th,
            .mutation-table th {
                padding: 11px 14px;
                text-align: left;
                color: #475569;
                background: #f8fafc;
                border-bottom: 2px solid #dbe4ee;
                font-size: 0.75rem;
                text-transform: uppercase;
                letter-spacing: 0.04em;
            }

            .result-table td,
            .mutation-table td {
                padding: 12px 14px;
                border-bottom: 1px solid #e2e8f0;
                color: #334155;
                line-height: 1.45;
            }

            .result-table tbody tr:hover {
                background: #f8fbff;
            }

            .technical-results > summary::-webkit-details-marker {
                display: none;
            }

            @media (max-width: 700px) {
                .app-banner-image {
                    height: 165px !important;
                }

                .metadata-grid {
                    grid-template-columns: minmax(0, 1fr) !important;
                }

                .metadata-grid > div {
                    grid-column: span 1 !important;
                }
            }

            @media (prefers-reduced-motion: reduce) {
                *, *::before, *::after {
                    scroll-behavior: auto !important;
                    transition-duration: 0.01ms !important;
                    animation-duration: 0.01ms !important;
                    animation-iteration-count: 1 !important;
                }
            }
        </style>
    </head>
    <body>
        {%app_entry%}
        <footer>
            {%config%}
            {%scripts%}
            {%renderer%}
        </footer>
    </body>
</html>
"""


AA_PATTERN = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$", re.IGNORECASE)
SAFE_NAME_PATTERN = re.compile(r"[^a-zA-Z0-9_-]+")

RUN_MODE_LABELS = {
    "complete": "Run complete",
    "structures": "Run structures only",
    "scores": "Run scores only",
    "graft": "Run CDR grafting only",
    "graft_backmutation": "Run CDR grafting only + Backmutation",
}

OPTIMIZATION_LABELS = {
    "1": "Level 1 - Critical validation",
    "2": "Level 2 - Automatic framework correction",
    "3": "Level 3 - Back-mutation analysis",
    "4": "Level 4 - Full scientific assessment (default)",
}

OPTIMIZATION_HELP = {
    "1": "Fixes only critical structural issues.",
    "2": "Applies targeted framework corrections.",
    "3": "Analyzes selective murine residue restoration.",
    "4": "Performs broader therapeutic risk assessment.",
}

# ---------- File and metadata helpers ----------

def sanitize_work_name(value: str) -> str:
    cleaned = SAFE_NAME_PATTERN.sub("_", value.strip())
    cleaned = re.sub(r"_+", "_", cleaned).strip("_")
    return cleaned[:80] or "work"


def parse_upload(contents: str) -> str:
    if not contents:
        raise ValueError("No file content received.")
    try:
        _, encoded = contents.split(",", 1)
        decoded = base64.b64decode(encoded)
        return decoded.decode("utf-8").strip()
    except UnicodeDecodeError as exc:
        raise ValueError("The uploaded file must be UTF-8 text.") from exc
    except Exception as exc:
        raise ValueError("Could not decode the uploaded file.") from exc


def extract_sequences(file_text: str) -> Tuple[str, str]:
    lines = [line.strip() for line in file_text.splitlines() if line.strip()]
    if len(lines) < 2:
        raise ValueError("The input file must contain at least 2 non-empty lines: VH on line 1 and VL on line 2.")

    return validate_sequence_pair(lines[0], lines[1])


def normalize_sequence_text(value: str) -> str:
    """Normalize a sequence pasted with spaces or line wrapping."""
    return re.sub(r"\s+", "", value or "").upper()


def validate_sequence_pair(vh_text: str, vl_text: str) -> Tuple[str, str]:
    """Validate and normalize a VH/VL sequence pair."""
    vh = normalize_sequence_text(vh_text)
    vl = normalize_sequence_text(vl_text)

    if len(vh) < 50 or len(vl) < 50:
        raise ValueError("VH and VL sequences must have at least 50 amino acids.")
    if not AA_PATTERN.fullmatch(vh):
        raise ValueError("VH contains invalid amino acid characters.")
    if not AA_PATTERN.fullmatch(vl):
        raise ValueError("VL contains invalid amino acid characters.")

    return vh, vl


def sequence_pair_to_input_text(vh: str, vl: str) -> str:
    """Create the two-line input format expected by the pipeline."""
    return f"{vh}\n{vl}\n"


def create_run_structure(
    work_name: str,
    uploaded_text: str,
    run_mode: str,
    optimization_level: Optional[str],
    backmutation: bool = False,
) -> Dict[str, Any]:
    safe_name = sanitize_work_name(work_name)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    run_id = f"{safe_name}_{timestamp}_{uuid.uuid4().hex[:8]}"
    run_dir = RUNS_DIR / run_id
    input_dir = run_dir / "input"
    output_dir = run_dir / "output"
    logs_dir = run_dir / "logs"

    input_dir.mkdir(parents=True, exist_ok=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)

    input_file = input_dir / "input_sequences.txt"
    input_file.write_text(uploaded_text.strip() + "\n", encoding="utf-8")

    metadata = {
        "run_id": run_id,
        "work_name": safe_name,
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "status": "queued",
        "input_file": str(input_file),
        "output_dir": str(output_dir),
        "log_file": str(logs_dir / "pipeline.log"),
        "return_code": None,
        "error": None,
        "database_path": str(DEFAULT_DATABASE),
        "run_mode": run_mode,
        "optimization_level": optimization_level,
        "backmutation": backmutation,
    }
    write_metadata(run_id, metadata)
    return metadata


def metadata_path(run_id: str) -> Path:
    return RUNS_DIR / run_id / "metadata.json"


def write_metadata(run_id: str, metadata: Dict[str, Any]) -> None:
    path = metadata_path(run_id)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8")


def read_metadata(run_id: str) -> Optional[Dict[str, Any]]:
    path = metadata_path(run_id)
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def update_metadata(run_id: str, **updates: Any) -> Dict[str, Any]:
    data = read_metadata(run_id) or {"run_id": run_id}
    data.update(updates)
    write_metadata(run_id, data)
    return data


def read_text_file(path: Path) -> str:
    if not path.exists():
        return "File not found."
    try:
        return path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        return path.read_text(encoding="latin-1", errors="replace")


def list_result_files(output_dir: Path) -> List[Path]:
    if not output_dir.exists():
        return []
    return sorted([path for path in output_dir.iterdir() if path.is_file()], key=lambda p: p.name.lower())


def format_status(status: str) -> str:
    labels = {
        "queued": "Queued",
        "running": "Running",
        "completed": "Completed",
        "failed": "Failed",
    }
    return labels.get(status, status.title())


# ---------- Pipeline execution ----------

def run_pipeline_in_background(run_id: str) -> None:
    metadata = read_metadata(run_id)
    if not metadata:
        return

    input_file = Path(metadata["input_file"])
    output_dir = Path(metadata["output_dir"])
    log_file = Path(metadata["log_file"])
    database_path = metadata.get("database_path") or str(DEFAULT_DATABASE)
    run_mode = metadata.get("run_mode", "complete")
    optimization_level = metadata.get("optimization_level")
    backmutation = bool(metadata.get("backmutation", False))

    command = [
        sys.executable,
        str(PIPELINE_SCRIPT),
        str(input_file),
        "-o",
        str(output_dir),
        "--database",
        str(database_path),
    ]

    if run_mode == "complete":
        if optimization_level:
            command.extend(["--optimization", str(optimization_level)])
    elif run_mode == "structures":
        command.append("--structures")
    elif run_mode == "scores":
        command.append("--scores")
    elif run_mode == "graft":
        command.append("--graft")
    elif run_mode == "graft_backmutation":
        command.append("--graft")
        backmutation = True
    else:
        raise ValueError(f"Invalid run mode: {run_mode}")

    if backmutation:
        command.append("--backmutation")

    update_metadata(run_id, status="running", started_at=datetime.now().isoformat(timespec="seconds"))

    try:
        with log_file.open("w", encoding="utf-8") as log_handle:
            log_handle.write(f"RUN MODE: {run_mode}\n")
            log_handle.write(f"OPTIMIZATION LEVEL: {optimization_level or '-'}\n")
            log_handle.write("COMMAND:\n")
            log_handle.write(" ".join(command) + "\n\n")
            log_handle.flush()

            process = subprocess.run(
                command,
                cwd=str(APP_DIR),
                stdout=log_handle,
                stderr=subprocess.STDOUT,
                text=True,
                check=False,
            )

        if process.returncode == 0:
            update_metadata(
                run_id,
                status="completed",
                return_code=process.returncode,
                finished_at=datetime.now().isoformat(timespec="seconds"),
            )
        else:
            update_metadata(
                run_id,
                status="failed",
                return_code=process.returncode,
                finished_at=datetime.now().isoformat(timespec="seconds"),
                error="The pipeline finished with a non-zero return code. Check the log file.",
            )
    except Exception as exc:
        trace = traceback.format_exc()
        log_file.write_text(trace, encoding="utf-8")
        update_metadata(
            run_id,
            status="failed",
            finished_at=datetime.now().isoformat(timespec="seconds"),
            error=str(exc),
        )


def start_background_run(run_id: str) -> None:
    worker = threading.Thread(target=run_pipeline_in_background, args=(run_id,), daemon=True)
    worker.start()


# ---------- Layout helpers ----------

def build_banner() -> html.Div:
    banner = html.Div(
        [
            html.Div(
                [
                    html.Div(
                        "ANTIBODY ENGINEERING",
                        style={
                            "fontSize": "0.72rem",
                            "fontWeight": "800",
                            "letterSpacing": "0.15em",
                            "color": "#315ea8",
                            "marginBottom": "10px",
                        },
                    ),
                    html.H1(
                        "Antibody Humanizer",
                        style={
                            "margin": "0",
                            "fontSize": "clamp(1.9rem, 4vw, 2.65rem)",
                            "lineHeight": "1.05",
                            "letterSpacing": "-0.045em",
                            "color": "#13295b",
                        },
                    ),
                    html.P(
                        "Sequence-based antibody humanization and scientific assessment.",
                        style={
                            "margin": "13px 0 0",
                            "maxWidth": "650px",
                            "fontSize": "1rem",
                            "lineHeight": "1.6",
                            "color": "#60718a",
                        },
                    ),
                ],
                style={
                    "position": "relative",
                    "zIndex": "2",
                },
            ),

            # Decoración suave en el lado derecho
            html.Div(
                style={
                    "position": "absolute",
                    "width": "320px",
                    "height": "320px",
                    "right": "-110px",
                    "top": "-160px",
                    "borderRadius": "50%",
                    "background": (
                        "radial-gradient(circle, "
                        "rgba(49, 94, 168, 0.15) 0%, "
                        "rgba(49, 94, 168, 0) 70%)"
                    ),
                    "pointerEvents": "none",
                }
            ),

            html.Div(
                style={
                    "position": "absolute",
                    "width": "240px",
                    "height": "240px",
                    "right": "90px",
                    "bottom": "-180px",
                    "borderRadius": "50%",
                    "background": (
                        "radial-gradient(circle, "
                        "rgba(34, 211, 238, 0.16) 0%, "
                        "rgba(34, 211, 238, 0) 72%)"
                    ),
                    "pointerEvents": "none",
                }
            ),
        ],
        style={
            "padding": "clamp(22px, 4vw, 30px)",
            "minHeight": "140px",
            "display": "flex",
            "alignItems": "center",
            "borderRadius": "18px",
            "background": (
                "linear-gradient(135deg, "
                "#fbfdff 0%, "
                "#f1f6fd 52%, "
                "#edf8fb 100%)"
            ),
            "border": "1px solid #dbe5f0",
            "boxShadow": "0 8px 28px rgba(15, 23, 42, 0.055)",
            "position": "relative",
            "overflow": "hidden",
        },
    )

    return html.Div(
        banner,
        style={"marginBottom": "18px"},
    )


def build_app_image() -> html.Div:
    if IMAGE_PATH.exists():
        encoded_image = base64.b64encode(
            IMAGE_PATH.read_bytes()
        ).decode("ascii")

        body = html.Img(
            src=f"data:image/png;base64,{encoded_image}",
            className="app-hero-image",
            style={
                "width": "100%",
                "height": "240px",
                "objectFit": "cover",
                "objectPosition": "center",
                "display": "block",
            },
        )
    else:
        body = html.Div(
            [
                html.Div(
                    "Antibody Humanization",
                    style={
                        "fontSize": "1.2rem",
                        "fontWeight": "800",
                    },
                ),
                html.Div(
                    "VH and VL sequence assessment",
                    style={
                        "marginTop": "6px",
                        "color": "#64748b",
                        "fontSize": "0.92rem",
                    },
                ),
            ],
            style={
                "minHeight": "150px",
                "display": "flex",
                "flexDirection": "column",
                "justifyContent": "center",
                "alignItems": "center",
                "background": "linear-gradient(145deg, #f8fbff, #eef5ff)",
                "color": "#16345c",
            },
        )

    return html.Div(
        body,
        style={
            "width": "100%",
            "marginBottom": "18px",
            "border": "1px solid #e2e8f0",
            "borderRadius": "18px",
            "background": "white",
            "boxShadow": "0 6px 20px rgba(15, 23, 42, 0.045)",
            "overflow": "hidden",
        },
    )


def card(children: Any) -> html.Div:
    return html.Div(
        children,
        className="content-card",
        style={
            "background": "rgba(255, 255, 255, 0.96)",
            "padding": "clamp(20px, 4vw, 30px)",
            "borderRadius": "18px",
            "border": "1px solid #e2e8f0",
            "boxShadow": "0 8px 28px rgba(15, 23, 42, 0.055)",
            "marginBottom": "18px",
        },
    )


def page_container(children: Any) -> html.Div:
    return html.Div(
        html.Div(
            children,
            style={
                "width": "100%",
                "maxWidth": "1080px",
                "margin": "0 auto",
                "padding": "clamp(16px, 4vw, 34px)",
                "boxSizing": "border-box",
            },
        ),
        style={
            "minHeight": "100vh",
            "background": "linear-gradient(180deg, #edf3f9 0%, #f5f8fc 42%, #eef3f8 100%)",
            "color": "#0f172a",
        },
    )


def home_page() -> html.Div:
    textarea_style = {
        "width": "100%",
        "minHeight": "178px",
        "padding": "15px 16px",
        "border": "1px solid #cbd5e1",
        "borderRadius": "13px",
        "backgroundColor": "#f8fafc",
        "color": "#0f172a",
        "fontFamily": "'SFMono-Regular', Consolas, 'Liberation Mono', monospace",
        "fontSize": "0.92rem",
        "lineHeight": "1.58",
        "resize": "vertical",
        "boxSizing": "border-box",
        "outline": "none",
    }

    label_style = {
        "display": "block",
        "fontWeight": "750",
        "marginBottom": "9px",
        "color": "#172033",
        "fontSize": "0.94rem",
    }

    choice_label_style = {
        "display": "flex",
        "alignItems": "center",
        "gap": "9px",
        "padding": "9px 10px",
        "marginBottom": "3px",
        "fontSize": "0.95rem",
        "lineHeight": "1.45",
        "color": "#172033",
        "cursor": "pointer",
        "borderRadius": "9px",
    }

    choice_box_style = {
        "width": "100%",
        "padding": "8px",
        "border": "1px solid #dbe3ed",
        "borderRadius": "14px",
        "backgroundColor": "#f8fafc",
        "boxSizing": "border-box",
    }

    return page_container(
        [
            build_banner(),
            build_app_image(),
            card(
                [
                    html.Div(
                        [
                            html.H2(
                                "New humanization run",
                                style={
                                    "margin": "0",
                                    "fontSize": "clamp(1.45rem, 3vw, 1.8rem)",
                                    "lineHeight": "1.2",
                                    "letterSpacing": "-0.03em",
                                    "color": "#0f172a",
                                },
                            ),
                            html.P(
                                "Paste the VH and VL amino acid sequences below, or upload a text file.",
                                style={
                                    "color": "#64748b",
                                    "fontSize": "0.96rem",
                                    "lineHeight": "1.6",
                                    "margin": "8px 0 0",
                                },
                            ),
                        ],
                        style={"marginBottom": "24px"},
                    ),
                    html.Div(
                        [
                            html.Div(
                                [
                                    html.Label(
                                        "VH sequence",
                                        htmlFor="vh-sequence-input",
                                        style=label_style,
                                    ),
                                    dcc.Textarea(
                                        id="vh-sequence-input",
                                        placeholder="Paste the VH amino acid sequence",
                                        spellCheck=False,
                                        className="sequence-textarea",
                                        style=textarea_style,
                                    ),
                                ],
                                style={"flex": "1 1 390px", "minWidth": "0"},
                            ),
                            html.Div(
                                [
                                    html.Label(
                                        "VL sequence",
                                        htmlFor="vl-sequence-input",
                                        style=label_style,
                                    ),
                                    dcc.Textarea(
                                        id="vl-sequence-input",
                                        placeholder="Paste the VL amino acid sequence",
                                        spellCheck=False,
                                        className="sequence-textarea",
                                        style=textarea_style,
                                    ),
                                ],
                                style={"flex": "1 1 390px", "minWidth": "0"},
                            ),
                        ],
                        style={
                            "display": "flex",
                            "gap": "18px",
                            "flexWrap": "wrap",
                            "marginBottom": "22px",
                        },
                    ),
                    html.Div(
                        [
                            html.Div(style={"height": "1px", "background": "#e2e8f0", "flex": "1"}),
                            html.Div(
                                "OR UPLOAD A FILE",
                                style={
                                    "padding": "0 14px",
                                    "textAlign": "center",
                                    "color": "#7b8a9f",
                                    "fontSize": "0.73rem",
                                    "fontWeight": "800",
                                    "letterSpacing": "0.09em",
                                    "whiteSpace": "nowrap",
                                },
                            ),
                            html.Div(style={"height": "1px", "background": "#e2e8f0", "flex": "1"}),
                        ],
                        style={"display": "flex", "alignItems": "center", "margin": "2px 0 16px"},
                    ),
                    dcc.Upload(
                        id="upload-input-file",
                        className="upload-zone",
                        children=html.Div(
                            [
                                html.Span("Drag and drop or "),
                                html.Strong("select a file", style={"color": "#1d4ed8"}),
                                html.Span(" (VH on line 1, VL on line 2)"),
                            ]
                        ),
                        multiple=False,
                        style={
                            "width": "100%",
                            "padding": "27px 22px",
                            "borderWidth": "1.5px",
                            "borderStyle": "dashed",
                            "borderColor": "#93b5df",
                            "borderRadius": "14px",
                            "textAlign": "center",
                            "background": "#f7fbff",
                            "color": "#4b6078",
                            "cursor": "pointer",
                            "margin": "0 auto 10px",
                            "boxSizing": "border-box",
                            "lineHeight": "1.6",
                        },
                    ),
                    html.Div(
                        id="uploaded-file-name",
                        style={
                            "minHeight": "20px",
                            "marginBottom": "20px",
                            "color": "#64748b",
                            "fontSize": "0.88rem",
                        },
                    ),
                    html.Label("Work name", htmlFor="work-name-input", style=label_style),
                    dcc.Input(
                        id="work-name-input",
                        type="text",
                        placeholder="Enter a clear name for this analysis",
                        className="standard-input",
                        style={
                            "width": "100%",
                            "height": "50px",
                            "padding": "0 15px",
                            "borderRadius": "12px",
                            "border": "1px solid #cbd5e1",
                            "backgroundColor": "#ffffff",
                            "marginBottom": "24px",
                            "fontSize": "0.97rem",
                            "lineHeight": "50px",
                            "boxSizing": "border-box",
                            "outline": "none",
                            "color": "#0f172a",
                        },
                    ),
                    html.Div(
                        [
                            html.Div(
                                [
                                    html.Label(
                                        "Select run mode",
                                        style={**label_style, "marginBottom": "10px"},
                                    ),
                                    dcc.RadioItems(
                                        id="run-mode-radio",
                                        className="choice-group",
                                        options=[
                                            {"label": "Run complete", "value": "complete"},
                                            {"label": "Run structures only", "value": "structures"},
                                            {"label": "Run scores only", "value": "scores"},
                                            {"label": "Run CDR grafting only", "value": "graft"},
                                            {"label": "Run CDR grafting only + Backmutation", "value": "graft_backmutation"},
                                        ],
                                        value="complete",
                                        labelStyle=choice_label_style,
                                        inputStyle={
                                            "margin": "0",
                                            "cursor": "pointer",
                                            "width": "16px",
                                            "height": "16px",
                                            "flexShrink": "0",
                                        },
                                        style={**choice_box_style, "marginBottom": "18px"},
                                    ),
                                ]
                            ),
                            html.Div(
                                id="optimization-container",
                                children=[
                                    html.Div(
                                        [
                                            html.Div(
                                                [
                                                    html.Span(
                                                        "Optimization level",
                                                        style={
                                                            "fontWeight": "750",
                                                            "color": "#172033",
                                                            "fontSize": "0.94rem",
                                                        },
                                                    ),
                                                    html.Span(
                                                        "?",
                                                        title="Optimization is used only in Run complete mode. The program supports levels 1 to 4, and the default full pipeline uses level 4 when no level is explicitly provided.",
                                                        className="help-badge",
                                                        style={
                                                            "display": "inline-flex",
                                                            "alignItems": "center",
                                                            "justifyContent": "center",
                                                            "width": "20px",
                                                            "height": "20px",
                                                            "marginLeft": "8px",
                                                            "borderRadius": "999px",
                                                            "background": "#dbeafe",
                                                            "color": "#1d4ed8",
                                                            "fontWeight": "800",
                                                            "cursor": "help",
                                                            "fontSize": "0.78rem",
                                                        },
                                                    ),
                                                ],
                                                style={
                                                    "display": "flex",
                                                    "alignItems": "center",
                                                    "marginBottom": "10px",
                                                },
                                            ),
                                            dcc.RadioItems(
                                                id="optimization-level-radio",
                                                className="choice-group",
                                                options=[
                                                    {
                                                        "label": html.Span(
                                                            [
                                                                html.Span("Level 1 - Critical validation", style={"marginRight": "8px"}),
                                                                html.Span(
                                                                    "?",
                                                                    title=OPTIMIZATION_HELP["1"],
                                                                    className="mini-help-badge",
                                                                    style={
                                                                        "display": "inline-flex",
                                                                        "alignItems": "center",
                                                                        "justifyContent": "center",
                                                                        "width": "18px",
                                                                        "height": "18px",
                                                                        "borderRadius": "999px",
                                                                        "background": "#e2e8f0",
                                                                        "color": "#475569",
                                                                        "fontWeight": "800",
                                                                        "cursor": "help",
                                                                        "fontSize": "0.7rem",
                                                                    },
                                                                ),
                                                            ],
                                                            style={"display": "inline-flex", "alignItems": "center", "flexWrap": "wrap"},
                                                        ),
                                                        "value": "1",
                                                    },
                                                    {
                                                        "label": html.Span(
                                                            [
                                                                html.Span("Level 2 - Automatic framework correction", style={"marginRight": "8px"}),
                                                                html.Span(
                                                                    "?",
                                                                    title=OPTIMIZATION_HELP["2"],
                                                                    className="mini-help-badge",
                                                                    style={
                                                                        "display": "inline-flex",
                                                                        "alignItems": "center",
                                                                        "justifyContent": "center",
                                                                        "width": "18px",
                                                                        "height": "18px",
                                                                        "borderRadius": "999px",
                                                                        "background": "#e2e8f0",
                                                                        "color": "#475569",
                                                                        "fontWeight": "800",
                                                                        "cursor": "help",
                                                                        "fontSize": "0.7rem",
                                                                    },
                                                                ),
                                                            ],
                                                            style={"display": "inline-flex", "alignItems": "center", "flexWrap": "wrap"},
                                                        ),
                                                        "value": "2",
                                                    },
                                                    {
                                                        "label": html.Span(
                                                            [
                                                                html.Span("Level 3 - Back-mutation analysis", style={"marginRight": "8px"}),
                                                                html.Span(
                                                                    "?",
                                                                    title=OPTIMIZATION_HELP["3"],
                                                                    className="mini-help-badge",
                                                                    style={
                                                                        "display": "inline-flex",
                                                                        "alignItems": "center",
                                                                        "justifyContent": "center",
                                                                        "width": "18px",
                                                                        "height": "18px",
                                                                        "borderRadius": "999px",
                                                                        "background": "#e2e8f0",
                                                                        "color": "#475569",
                                                                        "fontWeight": "800",
                                                                        "cursor": "help",
                                                                        "fontSize": "0.7rem",
                                                                    },
                                                                ),
                                                            ],
                                                            style={"display": "inline-flex", "alignItems": "center", "flexWrap": "wrap"},
                                                        ),
                                                        "value": "3",
                                                    },
                                                    {
                                                        "label": html.Span(
                                                            [
                                                                html.Span("Level 4 - Full scientific assessment (default)", style={"marginRight": "8px"}),
                                                                html.Span(
                                                                    "?",
                                                                    title=OPTIMIZATION_HELP["4"],
                                                                    className="mini-help-badge",
                                                                    style={
                                                                        "display": "inline-flex",
                                                                        "alignItems": "center",
                                                                        "justifyContent": "center",
                                                                        "width": "18px",
                                                                        "height": "18px",
                                                                        "borderRadius": "999px",
                                                                        "background": "#e2e8f0",
                                                                        "color": "#475569",
                                                                        "fontWeight": "800",
                                                                        "cursor": "help",
                                                                        "fontSize": "0.7rem",
                                                                    },
                                                                ),
                                                            ],
                                                            style={"display": "inline-flex", "alignItems": "center", "flexWrap": "wrap"},
                                                        ),
                                                        "value": "4",
                                                    },
                                                ],
                                                value="4",
                                                labelStyle=choice_label_style,
                                                inputStyle={
                                                    "margin": "0",
                                                    "cursor": "pointer",
                                                    "width": "16px",
                                                    "height": "16px",
                                                    "flexShrink": "0",
                                                },
                                            ),
                                        ],
                                        style={**choice_box_style, "marginBottom": "22px"},
                                    ),
                                ],
                            ),
                            html.Button(
                                "Run",
                                id="run-button",
                                n_clicks=0,
                                className="primary-run-button",
                                style={
                                    "width": "100%",
                                    "minHeight": "52px",
                                    "padding": "13px 20px",
                                    "border": "none",
                                    "borderRadius": "12px",
                                    "background": "linear-gradient(135deg, #2563eb 0%, #1d4ed8 100%)",
                                    "color": "white",
                                    "fontWeight": "800",
                                    "cursor": "pointer",
                                    "fontSize": "1rem",
                                    "letterSpacing": "0.01em",
                                    "boxShadow": "0 8px 20px rgba(37, 99, 235, 0.22)",
                                },
                            ),
                        ],
                        style={"width": "100%"},
                    ),
                    html.Div(id="run-feedback", className="feedback-box", style={"marginTop": "18px"}),
                    html.Div(id="run-live-status", style={"marginTop": "10px"}),
                    html.Div(id="results-link-container", style={"marginTop": "16px"}),
                ]
            ),
            dcc.Interval(id="run-status-interval", interval=STATUS_POLL_MS, n_intervals=0),
        ]
    )


def render_file_preview(file_path: Path) -> html.Div:
    title = html.H4(
        file_path.name,
        style={
            "margin": "0",
            "fontSize": "1rem",
            "letterSpacing": "-0.01em",
            "color": "#172033",
        },
    )
    content = read_text_file(file_path)
    preview = html.Pre(
        content[:20000] if content else "",
        className="file-preview",
        style={
            "whiteSpace": "pre-wrap",
            "wordBreak": "break-word",
            "maxHeight": "440px",
            "overflowY": "auto",
            "background": "#111827",
            "color": "#dbe7f3",
            "padding": "18px",
            "borderRadius": "13px",
            "fontSize": "0.88rem",
            "lineHeight": "1.55",
            "border": "1px solid #25334a",
            "margin": "14px 0 0",
        },
    )
    return card([title, preview])


def parse_candidate_summary(content: str) -> List[Dict[str, Any]]:
    """Extract candidate rows from humanization_summary.txt for visual presentation."""
    candidates: List[Dict[str, Any]] = []
    current_chain = None
    current_candidate = None
    for line in content.splitlines():
        chain_match = re.match(r"^(VH|VL) Chain:\s+(\d+) candidates", line)
        if chain_match:
            current_chain = chain_match.group(1)
            continue
        candidate_match = re.match(
            r"^\s+(\d+)\.\s+(.+?)\s+\(identity:\s+([\d.]+)%,\s+length:\s+(\d+) AA\)", line
        )
        if candidate_match and current_chain:
            current_candidate = {
                "chain": current_chain,
                "candidate": int(candidate_match.group(1)),
                "source": candidate_match.group(2),
                "identity": float(candidate_match.group(3)),
                "length": int(candidate_match.group(4)),
                "germline": "-",
            }
            candidates.append(current_candidate)
            continue
        germline_match = re.match(r"^\s+IMGT germline:\s+(.+)$", line)
        if germline_match and current_candidate:
            current_candidate["germline"] = germline_match.group(1)
    return candidates


def parse_backmutation_display(content: str) -> Dict[Tuple[str, int], Dict[str, Any]]:
    """Extract effective back-mutation counts and changes from the enriched report."""
    results: Dict[Tuple[str, int], Dict[str, Any]] = {}
    current_key = None
    for line in content.splitlines():
        candidate_match = re.match(r"^(VH|VL) candidate (\d+)$", line.strip())
        if candidate_match:
            current_key = (candidate_match.group(1), int(candidate_match.group(2)))
            results[current_key] = {"count": None, "changes": []}
            continue
        if current_key:
            count_match = re.match(r"^Back-mutations realizadas:\s+(\d+)$", line.strip())
            if count_match:
                results[current_key]["count"] = int(count_match.group(1))
                continue
            change_match = re.match(r"^\s+(\d+):\s+([A-Z])\s+->\s+([A-Z])$", line)
            if change_match:
                results[current_key]["changes"].append(
                    {"position": int(change_match.group(1)), "from": change_match.group(2), "to": change_match.group(3)}
                )
    return results


def result_metric(label: str, value: Any, accent: str = "#2563eb") -> html.Div:
    return html.Div(
        [
            html.Div(label, style={"fontSize": "0.75rem", "fontWeight": "800", "color": "#64748b", "textTransform": "uppercase"}),
            html.Div(str(value), style={"fontSize": "1.8rem", "fontWeight": "850", "color": "#0f172a", "marginTop": "6px"}),
        ],
        style={
            "background": "#ffffff", "border": "1px solid #e2e8f0", "borderTop": f"4px solid {accent}",
            "borderRadius": "13px", "padding": "16px", "minWidth": "0",
        },
    )


def render_candidate_table(chain: str, candidates: List[Dict[str, Any]], mutations: Dict[Tuple[str, int], Dict[str, Any]]) -> html.Div:
    rows = []
    for candidate in candidates:
        mutation_count = mutations.get((chain, candidate["candidate"]), {}).get("count")
        rows.append(
            html.Tr([
                html.Td(candidate["candidate"]),
                html.Td(candidate["source"], style={"fontFamily": "monospace"}),
                html.Td(candidate["germline"]),
                html.Td(f"{candidate['identity']:.1f}%"),
                html.Td("—" if mutation_count is None else mutation_count, style={"fontWeight": "800", "color": "#b45309"}),
            ])
        )
    table = html.Table(
        [
            html.Thead(html.Tr([html.Th(label) for label in ("Candidate", "Source", "IMGT germline", "Identity", "Back-mutations")])),
            html.Tbody(rows),
        ],
        className="result-table",
        style={"width": "100%", "borderCollapse": "collapse", "fontSize": "0.88rem"},
    )
    return card([
        html.H3(f"{chain} candidates", style={"margin": "0 0 14px", "color": "#0f172a"}),
        html.Div(table, style={"overflowX": "auto"}),
    ])


def render_backmutation_details(mutations: Dict[Tuple[str, int], Dict[str, Any]]) -> Optional[html.Div]:
    panels = []
    for (chain, candidate), data in sorted(mutations.items()):
        if data["count"] is None:
            continue
        changes = data["changes"]
        change_rows = [
            html.Tr([
                html.Td(change["position"]),
                html.Td(html.Span(change["from"], style={"background": "#dbeafe", "color": "#1d4ed8", "padding": "3px 9px", "borderRadius": "7px", "fontWeight": "800"})),
                html.Td("→", style={"color": "#94a3b8"}),
                html.Td(html.Span(change["to"], style={"background": "#ffedd5", "color": "#c2410c", "padding": "3px 9px", "borderRadius": "7px", "fontWeight": "800"})),
            ]) for change in changes
        ]
        panels.append(html.Div([
            html.H4(f"{chain} candidate {candidate} · {data['count']} back-mutations", style={"margin": "0 0 10px"}),
            html.Table([
                html.Thead(html.Tr([html.Th("Position"), html.Th("Humanized"), html.Th(""), html.Th("Backmutated")])),
                html.Tbody(change_rows),
            ], className="mutation-table", style={"borderCollapse": "collapse", "minWidth": "390px"})
            if changes else html.P("No effective residue changes.", style={"color": "#64748b", "margin": "0"}),
        ], style={"padding": "16px", "background": "#f8fafc", "border": "1px solid #e2e8f0", "borderRadius": "12px"}))
    if not panels:
        return None
    return card([
        html.H3("Back-mutation details", style={"margin": "0 0 6px"}),
        html.P("Positions use one-based sequence numbering.", style={"color": "#64748b", "margin": "0 0 16px"}),
        html.Div(panels, style={"display": "grid", "gap": "12px"}),
    ])


def results_page(query_string: str) -> html.Div:
    params = parse_qs((query_string or "").lstrip("?"))
    run_id = params.get("run", [""])[0]

    back_link_style = {
        "display": "inline-flex",
        "alignItems": "center",
        "color": "#1d4ed8",
        "fontWeight": "750",
        "textDecoration": "none",
    }

    if not run_id:
        return page_container(
            [
                build_banner(),
                card(
                    [
                        html.H2("Results", style={"margin": "0 0 10px", "letterSpacing": "-0.03em"}),
                        html.P("No run was selected.", style={"color": "#64748b"}),
                        dcc.Link("Back to home", href="/", style=back_link_style),
                    ]
                ),
            ]
        )

    metadata = read_metadata(run_id)
    if not metadata:
        return page_container(
            [
                build_banner(),
                card(
                    [
                        html.H2("Results", style={"margin": "0 0 10px", "letterSpacing": "-0.03em"}),
                        html.P("The selected run was not found.", style={"color": "#64748b"}),
                        dcc.Link("Back to home", href="/", style=back_link_style),
                    ]
                ),
            ]
        )

    output_dir = Path(metadata["output_dir"])
    log_file = Path(metadata["log_file"])
    summary_file = output_dir / "humanization_summary.txt"
    backmutation_file = output_dir / "backmutation.txt"
    backmutation_log_file = output_dir / "backmutation.log"
    vh_file = output_dir / "humanized_vh.fasta"
    vl_file = output_dir / "humanized_vl.fasta"
    files = list_result_files(output_dir)
    summary_content = read_text_file(summary_file) if summary_file.exists() else ""
    backmutation_content = read_text_file(backmutation_file) if backmutation_file.exists() else ""
    candidates = parse_candidate_summary(summary_content)
    mutations = parse_backmutation_display(backmutation_content)
    vh_candidates = [candidate for candidate in candidates if candidate["chain"] == "VH"]
    vl_candidates = [candidate for candidate in candidates if candidate["chain"] == "VL"]
    known_mutation_counts = [data["count"] for data in mutations.values() if data["count"] is not None]
    total_backmutations = sum(known_mutation_counts) if known_mutation_counts else "—"

    status = metadata.get("status", "unknown")
    status_badge = html.Span(
        format_status(status),
        style={
            "display": "inline-flex",
            "alignItems": "center",
            "padding": "7px 12px",
            "borderRadius": "999px",
            "background": "#dcfce7" if status == "completed" else "#fee2e2" if status == "failed" else "#e2e8f0",
            "color": "#166534" if status == "completed" else "#991b1b" if status == "failed" else "#475569",
            "fontWeight": "800",
            "fontSize": "0.85rem",
            "border": "1px solid #bbf7d0" if status == "completed" else "1px solid #fecaca" if status == "failed" else "1px solid #cbd5e1",
        },
    )

    file_list = (
        html.Ul(
            [html.Li(path.name, style={"marginBottom": "8px"}) for path in files],
            style={"margin": "12px 0 0", "paddingLeft": "20px", "color": "#334155"},
        )
        if files
        else html.P("No result files available yet.", style={"color": "#64748b", "marginBottom": "0"})
    )

    metadata_item_style = {
        "padding": "12px 14px",
        "borderRadius": "11px",
        "background": "#f8fafc",
        "border": "1px solid #e2e8f0",
        "minWidth": "0",
    }
    metadata_label_style = {
        "fontSize": "0.72rem",
        "fontWeight": "800",
        "letterSpacing": "0.07em",
        "textTransform": "uppercase",
        "color": "#7b8a9f",
        "marginBottom": "5px",
    }
    metadata_value_style = {
        "fontSize": "0.91rem",
        "lineHeight": "1.5",
        "color": "#243247",
        "overflowWrap": "anywhere",
    }
    metadata_code_style = {
        **metadata_value_style,
        "fontFamily": "'SFMono-Regular', Consolas, 'Liberation Mono', monospace",
        "fontSize": "0.82rem",
    }

    metadata_grid = html.Div(
        [
            html.Div(
                [
                    html.Div("Created at", style=metadata_label_style),
                    html.Div(metadata.get("created_at", "-"), style=metadata_value_style),
                ],
                style=metadata_item_style,
            ),
            html.Div(
                [
                    html.Div("Run mode", style=metadata_label_style),
                    html.Div(
                        RUN_MODE_LABELS.get(metadata.get("run_mode", "complete"), metadata.get("run_mode", "complete")),
                        style=metadata_value_style,
                    ),
                ],
                style=metadata_item_style,
            ),
            html.Div(
                [
                    html.Div("Optimization", style=metadata_label_style),
                    html.Div(
                        OPTIMIZATION_LABELS.get(str(metadata.get("optimization_level")), "-")
                        if metadata.get("optimization_level")
                        else "-",
                        style=metadata_value_style,
                    ),
                ],
                style=metadata_item_style,
            ),
            html.Div(
                [
                    html.Div("Return code", style=metadata_label_style),
                    html.Div(str(metadata.get("return_code")), style=metadata_code_style),
                ],
                style=metadata_item_style,
            ),
            html.Div(
                [
                    html.Div("Output folder", style=metadata_label_style),
                    html.Div(metadata.get("output_dir", "-"), style=metadata_code_style),
                ],
                style={**metadata_item_style, "gridColumn": "span 2"},
            ),
            html.Div(
                [
                    html.Div("Log file", style=metadata_label_style),
                    html.Div(metadata.get("log_file", "-"), style=metadata_code_style),
                ],
                style={**metadata_item_style, "gridColumn": "span 2"},
            ),
            html.Div(
                [
                    html.Div("Error", style=metadata_label_style),
                    html.Div(metadata.get("error") or "-", style=metadata_value_style),
                ],
                style={**metadata_item_style, "gridColumn": "span 2"},
            ),
        ],
        className="metadata-grid",
        style={
            "display": "grid",
            "gridTemplateColumns": "repeat(2, minmax(0, 1fr))",
            "gap": "10px",
            "marginTop": "20px",
        },
    )

    content_blocks = [
        build_banner(),
        card(
            [
                html.Div(
                    [
                        html.Div(
                            [
                                html.Div(
                                    "ANALYSIS RESULTS",
                                    style={
                                        "fontSize": "0.72rem",
                                        "fontWeight": "800",
                                        "letterSpacing": "0.1em",
                                        "color": "#64748b",
                                        "marginBottom": "7px",
                                    },
                                ),
                                html.H2(
                                    f"Results — {metadata.get('work_name', run_id)}",
                                    style={
                                        "margin": "0",
                                        "fontSize": "clamp(1.35rem, 3vw, 1.75rem)",
                                        "lineHeight": "1.2",
                                        "letterSpacing": "-0.035em",
                                        "color": "#0f172a",
                                    },
                                ),
                                html.Div(
                                    f"Run ID: {run_id}",
                                    style={
                                        "color": "#64748b",
                                        "fontFamily": "monospace",
                                        "fontSize": "0.82rem",
                                        "marginTop": "8px",
                                        "overflowWrap": "anywhere",
                                    },
                                ),
                            ],
                            style={"minWidth": "0"},
                        ),
                        html.Div(status_badge),
                    ],
                    style={
                        "display": "flex",
                        "justifyContent": "space-between",
                        "alignItems": "flex-start",
                        "gap": "16px",
                        "flexWrap": "wrap",
                    },
                ),
                metadata_grid,
                html.Div(
                    dcc.Link("Back to home", href="/", style=back_link_style),
                    style={"marginTop": "18px"},
                ),
            ]
        ),
        card(
            [
                html.H3("Results overview", style={"margin": "0 0 6px", "color": "#0f172a"}),
                html.P(
                    "Scientific results are summarized below. Technical files remain available at the end of the page.",
                    style={"color": "#64748b", "margin": "0 0 16px"},
                ),
                html.Div(
                    [
                        result_metric("VH candidates", len(vh_candidates), "#2563eb"),
                        result_metric("VL candidates", len(vl_candidates), "#0891b2"),
                        result_metric("Total candidates", len(candidates), "#7c3aed"),
                        result_metric("Effective back-mutations", total_backmutations, "#ea580c"),
                    ],
                    style={"display": "grid", "gridTemplateColumns": "repeat(auto-fit, minmax(165px, 1fr))", "gap": "12px"},
                ),
            ]
        ),
        card(
            [
                html.H3(
                    "Generated files",
                    style={"margin": "0", "fontSize": "1.08rem", "letterSpacing": "-0.02em"},
                ),
                file_list,
            ]
        ),
    ]

    if vh_candidates:
        content_blocks.append(render_candidate_table("VH", vh_candidates, mutations))
    if vl_candidates:
        content_blocks.append(render_candidate_table("VL", vl_candidates, mutations))
    mutation_details = render_backmutation_details(mutations)
    if mutation_details:
        content_blocks.append(mutation_details)

    technical_blocks = []
    for technical_file in (summary_file, backmutation_file, backmutation_log_file, vh_file, vl_file, log_file):
        if technical_file.exists():
            technical_blocks.append(render_file_preview(technical_file))
    if technical_blocks:
        content_blocks.append(
            html.Details(
                [
                    html.Summary(
                        "Show technical results and logs",
                        style={
                            "cursor": "pointer", "fontWeight": "800", "color": "#1d4ed8",
                            "padding": "18px 20px", "listStyle": "none",
                        },
                    ),
                    html.Div(technical_blocks, style={"display": "grid", "gap": "16px", "padding": "0 0 4px"}),
                ],
                className="technical-results",
                style={
                    "background": "#ffffff", "border": "1px solid #dce5ef", "borderRadius": "14px",
                    "boxShadow": "0 7px 22px rgba(15, 23, 42, 0.04)",
                },
            )
        )

    return page_container(content_blocks)


app.layout = html.Div(
    [
        dcc.Location(id="url", refresh=False),
        dcc.Store(id="current-run-store", storage_type="session"),
        html.Div(id="page-content"),
    ]
)


# ---------- Routing ----------
@app.callback(
    Output("page-content", "children"),
    Input("url", "pathname"),
    State("url", "search"),
)
def render_page(pathname: str, search: str) -> html.Div:
    if pathname == "/results":
        return results_page(search)
    return home_page()


# ---------- Home page interactions ----------
@app.callback(
    Output("uploaded-file-name", "children"),
    Input("upload-input-file", "filename"),
)
def show_uploaded_filename(filename: Optional[str]) -> str:
    if not filename:
        return "No file selected yet."
    return f"Selected file: {filename}"


@app.callback(
    Output("optimization-container", "style"),
    Input("run-mode-radio", "value"),
)
def toggle_optimization_container(run_mode: Optional[str]) -> Dict[str, str]:
    if run_mode == "complete":
        return {"display": "block"}
    return {"display": "none"}

@app.callback(
    Output("run-feedback", "children"),
    Output("current-run-store", "data"),
    Input("run-button", "n_clicks"),
    State("vh-sequence-input", "value"),
    State("vl-sequence-input", "value"),
    State("upload-input-file", "contents"),
    State("upload-input-file", "filename"),
    State("work-name-input", "value"),
    State("run-mode-radio", "value"),
    State("optimization-level-radio", "value"),
    prevent_initial_call=True,
)
def launch_run(
    n_clicks: int,
    vh_text: Optional[str],
    vl_text: Optional[str],
    uploaded_contents: Optional[str],
    filename: Optional[str],
    work_name: Optional[str],
    run_mode: Optional[str],
    optimization_level: Optional[str],
):
    if not n_clicks:
        return no_update, no_update

    if not work_name or not work_name.strip():
        return html.Div("Please write the work name.", style={"color": "#b91c1c", "fontWeight": "700"}), no_update

    run_mode = run_mode or "complete"
    effective_optimization = optimization_level if run_mode == "complete" else None
    backmutation = run_mode == "graft_backmutation"

    try:
        has_vh_text = bool(vh_text and vh_text.strip())
        has_vl_text = bool(vl_text and vl_text.strip())

        if has_vh_text or has_vl_text:
            if not (has_vh_text and has_vl_text):
                raise ValueError("Please fill in both VH and VL sequence fields.")
            vh, vl = validate_sequence_pair(vh_text or "", vl_text or "")
            input_text = sequence_pair_to_input_text(vh, vl)
            input_source = "Pasted VH/VL sequences"
        elif uploaded_contents:
            input_text = parse_upload(uploaded_contents)
            vh, vl = extract_sequences(input_text)
            input_text = sequence_pair_to_input_text(vh, vl)
            input_source = f"Uploaded file: {filename or '-'}"
        else:
            raise ValueError("Paste both VH and VL sequences or select an input file.")

        metadata = create_run_structure(
            work_name,
            input_text,
            run_mode,
            effective_optimization,
            backmutation,
        )
        start_background_run(metadata["run_id"])

        message_children = [
            html.Div(
                "Run started successfully.",
                style={"color": "#166534", "fontWeight": "700", "marginBottom": "6px"},
            ),
            html.Div(f"Selected mode: {RUN_MODE_LABELS.get(run_mode, run_mode)}"),
        ]

        if effective_optimization:
            message_children.append(
                html.Div(f"Optimization level: {OPTIMIZATION_LABELS.get(effective_optimization, effective_optimization)}")
            )
        if backmutation:
            message_children.append(html.Div("Backmutation: enabled"))

        message_children.extend(
            [
                html.Div(f"Input source: {input_source}"),
                html.Div(f"VH length: {len(vh)} aa"),
                html.Div(f"VL length: {len(vl)} aa"),
                html.Div(f"Run ID: {metadata['run_id']}"),
            ]
        )

        return html.Div(message_children), {"run_id": metadata["run_id"]}

    except Exception as exc:
        return html.Div(str(exc), style={"color": "#b91c1c", "fontWeight": "700"}), no_update


@app.callback(
    Output("run-live-status", "children"),
    Output("results-link-container", "children"),
    Input("run-status-interval", "n_intervals"),
    State("current-run-store", "data"),
)
def poll_run_status(_: int, store_data: Optional[Dict[str, Any]]):
    if not store_data or not store_data.get("run_id"):
        return "", ""

    run_id = store_data["run_id"]
    metadata = read_metadata(run_id)
    if not metadata:
        return html.Div("Run metadata not found."), ""

    status = metadata.get("status", "unknown")
    status_text = html.Div(
        f"Current status: {format_status(status)}",
        style={
            "fontWeight": "700",
            "color": "#1d4ed8" if status == "running" else "#166534" if status == "completed" else "#b91c1c" if status == "failed" else "#334155",
        },
    )

    if status == "completed":
        link = dcc.Link(
            "Open results page",
            href=f"/results?run={run_id}",
            style={
                "display": "inline-block",
                "padding": "10px 14px",
                "borderRadius": "12px",
                "background": "#166534",
                "color": "white",
                "fontWeight": "700",
                "textDecoration": "none",
            },
        )
        return status_text, link

    if status == "failed":
        link = dcc.Link(
            "Open results/log page",
            href=f"/results?run={run_id}",
            style={
                "display": "inline-block",
                "padding": "10px 14px",
                "borderRadius": "12px",
                "background": "#b91c1c",
                "color": "white",
                "fontWeight": "700",
                "textDecoration": "none",
            },
        )
        return status_text, link

    return status_text, ""


if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=8050)
