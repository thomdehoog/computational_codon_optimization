"""Local web server behind the codonopt interface.

Standard library only, so the interface works wherever the package does. The
browser (or a ``pywebview`` window) loads a pre-built page from
``codonopt/web`` and talks to this module over a small JSON API:

    GET  /api/meta               tables, motif groups, defaults, limits
    POST /api/jobs               start an optimisation, returns {"id": ...}
    GET  /api/jobs/<id>          progress, and the result once finished
    POST /api/jobs/<id>/cancel   stop early and keep the best sequence so far

Optimisation can take minutes, so a request starts a job on a worker thread and
the page polls for its state rather than holding a connection open.

The server binds to the loopback interface only, and still defends itself,
because anything running in the user's browser can reach ``127.0.0.1``:

* the ``Host`` header must name the loopback address and this server's port,
  which stops DNS rebinding;
* every ``/api`` call must carry a random per-process token that is only ever
  embedded in the page this server serves, which stops other sites from driving
  it;
* the API never reads or writes files -- sequences and tables arrive as text;
* responses carry a strict Content-Security-Policy.
"""

from __future__ import annotations

import hmac
import json
import math
import secrets
import sys
import threading
import time
import traceback
import uuid
import webbrowser
from dataclasses import dataclass, field, fields
from http import HTTPStatus
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from importlib import resources
from urllib.parse import urlsplit

from . import __version__
from .config import DesignConfig, OptimizerSettings, Weights
from .fasta import Record, format_fasta, parse_fasta
from .metrics import gc_profile, motif_hits
from .motifs import (
    MOTIF_GROUPS,
    REGULATORY_MOTIFS,
    RESTRICTION_SITES,
    resolve_motifs,
)
from .optimize import OptimizationResult, optimize
from .report import report_dict
from .sequence import AA_TO_CODONS, translate, validate_dna, validate_protein
from .splice import get_scorer, maxent_available
from .tables import TABLES, get_table, table_from_text

__all__ = [
    "App",
    "JobSpec",
    "RequestError",
    "build_request",
    "create_server",
    "run_gui",
]

# --------------------------------------------------------------------------
# Limits. The server is local, but a runaway request should still be bounded.
# --------------------------------------------------------------------------

MAX_BODY_BYTES = 2_000_000
MAX_RESIDUES = 15_000
MAX_SECONDS = 900
MAX_RESTARTS = 20
MAX_ANNEAL_STEPS = 1_000_000
MAX_RUNNING_JOBS = 2
MAX_KEPT_JOBS = 30

#: How hard each effort level searches.
EFFORT: dict[str, dict[str, int]] = {
    "quick": {"max_sweeps": 3, "restarts": 1, "anneal_steps": 0},
    "standard": {"max_sweeps": 6, "restarts": 1, "anneal_steps": 0},
    "thorough": {"max_sweeps": 8, "restarts": 3, "anneal_steps": 5000},
}

EXAMPLE_NAME = "example"
EXAMPLE_PROTEIN = (
    "MKWVTFISLLLLFSSAYSRGVFRRDTHKSEIAHRFKDLGEENFKALVLIAFAQYLQQCPF"
    "DEHVKLVNELTEFAKTCVADESHAGCEKSLHTLFGDELCKVASLRETYGDMADCCEKQEP"
)

_ASSETS: dict[str, tuple[str, str]] = {
    "/": ("index.html", "text/html; charset=utf-8"),
    "/index.html": ("index.html", "text/html; charset=utf-8"),
    "/app.js": ("app.js", "text/javascript; charset=utf-8"),
    "/app.css": ("app.css", "text/css; charset=utf-8"),
    "/favicon.svg": ("favicon.svg", "image/svg+xml"),
}

_CSP = (
    "default-src 'none'; script-src 'self'; style-src 'self'; "
    "img-src 'self' data:; connect-src 'self'; base-uri 'none'; "
    "form-action 'none'; frame-ancestors 'none'"
)
_TOKEN_PLACEHOLDER = "__CODONOPT_TOKEN__"
_TOKEN_HEADER = "X-Codonopt-Token"


class RequestError(ValueError):
    """The request was understood but cannot be acted on; the text is for the user."""


# --------------------------------------------------------------------------
# Turning a JSON payload into a validated job
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class JobSpec:
    """A validated optimisation request.

    Attributes:
        name: Record name for the output.
        protein: Protein to encode, with a trailing ``*`` if a stop is wanted.
        start: Coding sequence to improve, or None to build one.
        config: Design constraints.
        settings: Keyword arguments for :class:`OptimizerSettings`.
        kind: ``"protein"`` or ``"dna"``, as interpreted.
    """

    name: str
    protein: str
    start: str | None
    config: DesignConfig
    settings: dict
    kind: str


def _number(
    payload: dict,
    key: str,
    default,
    *,
    low: float,
    high: float,
    integer: bool = False,
    label: str | None = None,
):
    """Read an optional numeric field, with a readable error if it is wrong."""
    value = payload.get(key)
    if value is None or value == "":
        return default
    label = label or key.replace("_", " ")
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise RequestError(f"{label} must be a number")
    try:
        # A JSON integer has no width limit, so this is where one too large for
        # a float arrives; it is a bad request, not a server fault.
        if not math.isfinite(value):
            raise RequestError(f"{label} must be a finite number")
    except OverflowError:
        raise RequestError(f"{label} must be between {low:g} and {high:g}") from None
    if integer and int(value) != value:
        raise RequestError(f"{label} must be a whole number")
    if not low <= value <= high:
        raise RequestError(f"{label} must be between {low:g} and {high:g}")
    return int(value) if integer else float(value)


def _parse_sequence(payload: dict) -> tuple[str, str, str]:
    """Return ``(name, sequence, kind)`` from the pasted text."""
    raw = payload.get("sequence")
    if not isinstance(raw, str) or not raw.strip():
        raise RequestError("Paste a protein or DNA sequence first.")
    name = payload.get("name")
    name = name.strip() if isinstance(name, str) and name.strip() else None

    text = raw.strip()
    if text.startswith(">"):
        records = list(parse_fasta(text.splitlines()))
        if len(records) != 1:
            raise RequestError(
                f"Found {len(records)} FASTA records; paste one sequence at a time."
            )
        name = name or records[0].name
        sequence = records[0].sequence
    else:
        sequence = "".join(text.split())
    sequence = sequence.upper()

    kind = payload.get("kind", "auto")
    if kind not in {"auto", "protein", "dna"}:
        raise RequestError("kind must be auto, protein or dna")
    if kind == "auto":
        looks_like_dna = len(sequence) >= 9 and set(sequence) <= set("ACGTU")
        kind = "dna" if looks_like_dna else "protein"
    return name or "sequence", sequence, kind


def _build_table(payload: dict):
    text = payload.get("table_text")
    if isinstance(text, str) and text.strip():
        stripped = text.lstrip()
        fmt = "json" if stripped.startswith("{") else ("tsv" if "\t" in text else "csv")
        return table_from_text(text, fmt, name="custom", source="custom table")
    if text is not None and "table" not in payload:
        # A custom table was asked for and left blank. Falling back to a
        # built-in host here would optimise against the wrong organism.
        raise RequestError("The custom codon table is empty; paste one or choose a host.")
    name = payload.get("table", "human")
    if not isinstance(name, str):
        raise RequestError("table must be a name")
    try:
        return get_table(name)
    except KeyError as exc:
        raise RequestError(str(exc.args[0])) from None


def _build_weights(raw) -> Weights:
    if raw is None:
        return Weights()
    if not isinstance(raw, dict):
        raise RequestError("weights must be an object")
    known = {field.name for field in fields(Weights)}
    values = {}
    for key, value in raw.items():
        if key not in known:
            raise RequestError(f"unknown weight {key!r}")
        values[key] = _number(
            {key: value},
            key,
            getattr(Weights(), key),
            low=0,
            high=1e6,
            label=f"weight {key}",
        )
    return Weights(**values)


def build_request(payload: object) -> JobSpec:
    """Validate a decoded JSON request body into a :class:`JobSpec`.

    Raises:
        RequestError: With a message fit to show the user, for anything wrong
            with the request.
    """
    if not isinstance(payload, dict):
        raise RequestError("The request body must be a JSON object.")
    try:
        name, sequence, kind = _parse_sequence(payload)
        if kind == "dna":
            start = validate_dna(sequence.replace("U", "T"), in_frame=True)
            protein = translate(start)
        else:
            start = None
            protein = validate_protein(sequence)
            if payload.get("add_stop") and not protein.endswith("*"):
                protein += "*"
        if len(protein) > MAX_RESIDUES:
            raise RequestError(
                f"That is {len(protein):,} residues; the interface accepts up to "
                f"{MAX_RESIDUES:,}. Use the command line for longer genes."
            )

        avoid = payload.get("avoid", [])
        if not isinstance(avoid, list) or not all(isinstance(a, str) for a in avoid):
            raise RequestError("avoid must be a list of motif names or patterns")
        if len(avoid) > 100 or any(len(a) > 200 for a in avoid):
            raise RequestError("too many, or too long, motifs to avoid")
        motifs = tuple(resolve_motifs(avoid))

        gc_range = payload.get("gc_range", [0.30, 0.70])
        if not (isinstance(gc_range, list) and len(gc_range) == 2):
            raise RequestError("gc range must be two numbers")
        gc_low = _number(
            {"v": gc_range[0]}, "v", 0.30, low=0, high=1, label="GC range low"
        )
        gc_high = _number(
            {"v": gc_range[1]}, "v", 0.70, low=0, high=1, label="GC range high"
        )

        table = _build_table(payload)
        backend = payload.get("splice_backend", "none")
        if backend not in {"auto", "maxent", "motif", "none"}:
            raise RequestError("splice backend must be auto, maxent, motif or none")
        try:
            scorer = get_scorer(backend, host=table.name)
        except RuntimeError as exc:
            raise RequestError(str(exc).splitlines()[0]) from None

        def whole(key: str, default: int, low: int, high: int) -> int:
            return _number(payload, key, default, low=low, high=high, integer=True)

        def items(key: str, default: list, kind: type) -> tuple:
            value = payload.get(key, default)
            if not isinstance(value, list) or any(type(x) is not kind for x in value):
                raise RequestError(f"{key} must be a list of {kind.__name__} values")
            return tuple(value)

        left_flank = payload.get("left_flank", "")
        right_flank = payload.get("right_flank", "")
        if not isinstance(left_flank, str) or not isinstance(right_flank, str):
            raise RequestError("flanks must be DNA strings")
        # Compare bases with bases: `sequence` is residues for a protein input.
        construct_bases = len(left_flank) + 3 * len(protein) + len(right_flank)
        if construct_bases > 3 * MAX_RESIDUES:
            raise RequestError(
                f"the complete construct is {construct_bases:,} bases; the interface "
                f"accepts up to {3 * MAX_RESIDUES:,}, flanks included"
            )

        config = DesignConfig(
            table=table,
            gc_target=_number(payload, "gc_target", 0.50, low=0, high=1),
            gc_range=(gc_low, gc_high),
            gc_window=whole("gc_window", 30, 6, 300),
            gc_step=whole("gc_step", 3, 1, 300),
            repeat_k=whole("repeat_k", 12, 4, 60),
            hairpin_stem=whole("hairpin_stem", 10, 4, 40),
            hairpin_min_loop=whole("hairpin_min_loop", 3, 0, 100),
            hairpin_max_loop=whole("hairpin_max_loop", 60, 0, 300),
            max_homopolymer=whole("max_homopolymer", 6, 2, 30),
            motifs=motifs,
            splice_scorer=scorer,
            donor_cutoff=_number(payload, "donor_cutoff", None, low=-100, high=100),
            acceptor_cutoff=_number(
                payload, "acceptor_cutoff", None, low=-100, high=100
            ),
            weights=_build_weights(payload.get("weights")),
            required=items("required", ["motif"], str),
            protected_codons=items("protected_codons", [], int),
            excluded_codons=items("excluded_codons", [], str),
            splice_strand=payload.get("splice_strand", "forward"),
            left_flank=left_flank,
            right_flank=right_flank,
        )

        # The optimizer rejects these too, but inside the worker thread, where
        # the caller would see them as a failed job rather than a bad request.
        if config.protected_codons and start is None:
            raise RequestError(
                "protecting codons needs a starting DNA sequence to protect them from; "
                "paste DNA rather than a protein"
            )
        if any(i >= len(protein) for i in config.protected_codons):
            raise RequestError(
                f"protected codon index is outside the {len(protein)}-codon sequence"
            )
        for index, residue in enumerate(protein):
            if index in config.protected_codons:
                continue
            if not set(AA_TO_CODONS[residue]) - set(config.excluded_codons):
                raise RequestError(
                    f"every codon for {residue} is excluded, so position {index} "
                    "cannot be encoded"
                )

        effort = payload.get("effort", "standard")
        if effort not in EFFORT:
            raise RequestError("effort must be quick, standard or thorough")
        settings = dict(EFFORT[effort])
        settings["restarts"] = min(settings["restarts"], MAX_RESTARTS)
        settings["anneal_steps"] = min(settings["anneal_steps"], MAX_ANNEAL_STEPS)
        settings["max_seconds"] = _number(
            payload, "seconds", 120.0, low=1, high=MAX_SECONDS
        )
        seed = _number(payload, "seed", None, low=0, high=2**31 - 1, integer=True)
        if seed is not None:
            settings["seed"] = seed
        # Surface bad settings now rather than inside the worker thread.
        OptimizerSettings(**settings)
    except RequestError:
        raise
    except ValueError as exc:  # SequenceError, FASTA, table and config errors
        raise RequestError(str(exc)) from None

    return JobSpec(
        name=name,
        protein=protein,
        start=start,
        config=config,
        settings=settings,
        kind=kind,
    )


# --------------------------------------------------------------------------
# Jobs
# --------------------------------------------------------------------------


@dataclass
class Job:
    """One optimisation run and everything the page needs to show about it."""

    id: str
    spec: JobSpec
    status: str = "running"  # running | done | cancelled | error
    messages: list[str] = field(default_factory=list)
    restarts_done: int = 0
    started: float = field(default_factory=time.perf_counter)
    finished: float | None = None
    result: dict | None = None
    error: str | None = None
    stop: threading.Event = field(default_factory=threading.Event)

    def snapshot(self) -> dict:
        """Return the JSON-serialisable state of the job."""
        end = self.finished if self.finished is not None else time.perf_counter()
        body = {
            "id": self.id,
            "status": self.status,
            "elapsed": round(end - self.started, 2),
            "max_seconds": self.spec.settings["max_seconds"],
            "restarts_total": self.spec.settings["restarts"],
            "restarts_done": self.restarts_done,
            "messages": self.messages[-12:],
        }
        if self.error is not None:
            body["error"] = self.error
        if self.result is not None:
            body["result"] = self.result
        return body


class JobStore:
    """Thread-safe registry of jobs, bounded in size and concurrency."""

    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._jobs: dict[str, Job] = {}

    def start(self, spec: JobSpec) -> Job:
        """Create a job and run it on a worker thread.

        Raises:
            RequestError: If too many jobs are already running.
        """
        with self._lock:
            running = sum(job.status == "running" for job in self._jobs.values())
            if running >= MAX_RUNNING_JOBS:
                raise RequestError(
                    "Already running as many optimisations as allowed; "
                    "wait for one to finish or cancel it."
                )
            job = Job(id=uuid.uuid4().hex[:12], spec=spec)
            self._jobs[job.id] = job
            finished = [j for j in self._jobs.values() if j.status != "running"]
            for old in finished[: max(0, len(self._jobs) - MAX_KEPT_JOBS)]:
                del self._jobs[old.id]
        threading.Thread(target=_run_job, args=(job,), daemon=True).start()
        return job

    def get(self, job_id: str) -> Job | None:
        with self._lock:
            return self._jobs.get(job_id)


def _run_job(job: Job) -> None:
    spec = job.spec

    def progress(message: str) -> None:
        job.messages.append(message)
        # One "final score" line is printed as each restart finishes.
        if "final score" in message:
            job.restarts_done += 1

    try:
        settings = OptimizerSettings(
            **spec.settings, progress=progress, should_stop=job.stop.is_set
        )
        result = optimize(spec.protein, spec.config, settings, start=spec.start)
        job.result = result_payload(result, spec)
        job.status = "cancelled" if result.cancelled else "done"
    except (ValueError, KeyError) as exc:
        job.error = str(exc.args[0]) if exc.args else str(exc)
        job.status = "error"
    except Exception as exc:  # pragma: no cover - a bug, not a user error
        traceback.print_exc(file=sys.stderr)
        job.error = f"Unexpected error: {type(exc).__name__}: {exc}"
        job.status = "error"
    finally:
        job.finished = time.perf_counter()


def result_payload(result: OptimizationResult, spec: JobSpec) -> dict:
    """Everything the page needs to show a finished run."""
    config = spec.config
    payload = report_dict(result, config)
    profile_before = gc_profile(
        config.left_flank + result.initial_sequence + config.right_flank,
        window=config.gc_window, step=config.gc_step
    )
    profile_after = gc_profile(
        config.left_flank + result.sequence + config.right_flank,
        window=config.gc_window, step=config.gc_step
    )
    lengths = {motif.name: motif.length for motif in config.motifs}
    hits = motif_hits(config.left_flank + result.sequence + config.right_flank, list(config.motifs))
    payload.update(
        name=spec.name,
        kind=spec.kind,
        protein=result.protein,
        initial_sequence=result.initial_sequence,
        cancelled=result.cancelled,
        gc={
            "window": config.gc_window,
            "step": config.gc_step,
            "range": list(config.gc_range),
            "target": config.gc_target,
            "coordinates": "construct",
            "starts": [start for start, _ in profile_after],
            "initial": [round(value, 4) for _, value in profile_before],
            "final": [round(value, 4) for _, value in profile_after],
        },
        motif_hits=[
            {"name": name, "length": lengths[name], "positions": list(positions)}
            for name, positions in hits.items()
        ],
        fasta=format_fasta(
            Record(
                name=spec.name,
                description=(
                    f"{spec.name} codonopt table={config.table.name} "
                    f"score={result.score:.2f} valid={str(result.valid).lower()}"
                ),
                sequence=result.sequence,
            )
        ),
    )
    return payload


# --------------------------------------------------------------------------
# HTTP
# --------------------------------------------------------------------------


def meta_payload() -> dict:
    """Static facts the page needs to build its form."""
    defaults = DesignConfig()
    return {
        "version": __version__,
        "tables": [
            {"name": name, "description": table.description}
            for name, table in sorted(TABLES.items())
        ],
        "default_table": "human",
        "motif_groups": {name: list(members) for name, members in MOTIF_GROUPS.items()},
        "restriction_sites": dict(sorted(RESTRICTION_SITES.items())),
        "regulatory_motifs": dict(sorted(REGULATORY_MOTIFS.items())),
        "splice_backends": ["auto", "maxent", "motif", "none"],
        "maxent_available": maxent_available(),
        "weights": {
            field.name: getattr(Weights(), field.name) for field in fields(Weights)
        },
        "defaults": {
            "gc_target": defaults.gc_target,
            "gc_range": list(defaults.gc_range),
            "gc_window": defaults.gc_window,
            "gc_step": defaults.gc_step,
            "repeat_k": defaults.repeat_k,
            "hairpin_stem": defaults.hairpin_stem,
            "hairpin_min_loop": defaults.hairpin_min_loop,
            "hairpin_max_loop": defaults.hairpin_max_loop,
            "max_homopolymer": defaults.max_homopolymer,
            "seconds": 120,
        },
        "effort": EFFORT,
        "limits": {
            "max_residues": MAX_RESIDUES,
            "max_seconds": MAX_SECONDS,
        },
        "example": {"name": EXAMPLE_NAME, "protein": EXAMPLE_PROTEIN},
    }


class App:
    """State shared by every request: the token, the jobs, the allowed hosts."""

    def __init__(self, port: int = 0) -> None:
        self.token = secrets.token_urlsafe(24)
        self.jobs = JobStore()
        self.port = port

    def host_allowed(self, host: str | None) -> bool:
        return host is not None and host.lower() in {
            f"127.0.0.1:{self.port}",
            f"localhost:{self.port}",
        }

    def origin_allowed(self, origin: str | None) -> bool:
        return origin is None or origin.lower() in {
            f"http://127.0.0.1:{self.port}",
            f"http://localhost:{self.port}",
        }

    def token_ok(self, supplied: str | None) -> bool:
        return supplied is not None and hmac.compare_digest(supplied, self.token)


class _Handler(BaseHTTPRequestHandler):
    # HTTP/1.0: one request per connection, so an early error response can never
    # leave an unread request body to be mistaken for the next request.
    server_version = "codonopt"

    def log_message(self, format: str, *args) -> None:  # noqa: A002
        if getattr(self.server, "verbose", False):
            super().log_message(format, *args)

    @property
    def app(self) -> App:
        return self.server.app  # type: ignore[attr-defined]

    def do_GET(self) -> None:  # noqa: N802
        self._handle("GET")

    def do_POST(self) -> None:  # noqa: N802
        self._handle("POST")

    def _handle(self, method: str) -> None:
        try:
            path = urlsplit(self.path).path
            if not self.app.host_allowed(self.headers.get("Host")):
                return self._error(HTTPStatus.FORBIDDEN, "Unexpected Host header.")
            if path.startswith("/api/"):
                if not self.app.token_ok(self.headers.get(_TOKEN_HEADER)):
                    return self._error(HTTPStatus.FORBIDDEN, "Missing or wrong token.")
                if method == "POST" and not self.app.origin_allowed(
                    self.headers.get("Origin")
                ):
                    return self._error(HTTPStatus.FORBIDDEN, "Unexpected Origin.")
                return self._api(method, path)
            if method != "GET":
                return self._error(HTTPStatus.METHOD_NOT_ALLOWED, "GET only.")
            return self._static(path)
        except (BrokenPipeError, ConnectionResetError):
            return
        except Exception:  # pragma: no cover - a bug, not a user error
            traceback.print_exc(file=sys.stderr)
            self._error(HTTPStatus.INTERNAL_SERVER_ERROR, "Internal server error.")

    # -- static ------------------------------------------------------------

    def _static(self, path: str) -> None:
        asset = _ASSETS.get(path)
        if asset is None:
            return self._error(HTTPStatus.NOT_FOUND, "Not found.")
        filename, content_type = asset
        try:
            data = resources.files("codonopt").joinpath("web", filename).read_bytes()
        except (FileNotFoundError, OSError):
            return self._error(
                HTTPStatus.INTERNAL_SERVER_ERROR,
                "The interface has not been built. From a source checkout run "
                "'npm install && npm run build' inside web/.",
            )
        if filename == "index.html":
            data = data.replace(_TOKEN_PLACEHOLDER.encode(), self.app.token.encode())
        self._send(HTTPStatus.OK, data, content_type)

    # -- api ---------------------------------------------------------------

    def _api(self, method: str, path: str) -> None:
        parts = [part for part in path.split("/") if part]  # ["api", ...]
        if method == "GET" and parts == ["api", "meta"]:
            return self._json(HTTPStatus.OK, meta_payload())
        if method == "POST" and parts == ["api", "jobs"]:
            return self._start_job()
        if len(parts) >= 3 and parts[:2] == ["api", "jobs"]:
            job = self.app.jobs.get(parts[2])
            if job is None:
                return self._error(HTTPStatus.NOT_FOUND, "No such job.")
            if method == "GET" and len(parts) == 3:
                return self._json(HTTPStatus.OK, job.snapshot())
            if method == "POST" and parts[3:] == ["cancel"]:
                job.stop.set()
                return self._json(HTTPStatus.OK, job.snapshot())
        return self._error(HTTPStatus.NOT_FOUND, "Not found.")

    def _start_job(self) -> None:
        try:
            length = int(self.headers.get("Content-Length", ""))
        except ValueError:
            return self._error(HTTPStatus.LENGTH_REQUIRED, "Content-Length required.")
        if length < 0 or length > MAX_BODY_BYTES:
            return self._error(
                HTTPStatus.REQUEST_ENTITY_TOO_LARGE, "That request is too large."
            )
        try:
            # ValueError also covers a number with more digits than int() will
            # parse (CPython's 4300-digit limit).
            payload = json.loads(self.rfile.read(length) or b"null")
        except (UnicodeDecodeError, ValueError):
            return self._error(HTTPStatus.BAD_REQUEST, "The request is not valid JSON.")
        try:
            job = self.app.jobs.start(build_request(payload))
        except RequestError as exc:
            return self._error(HTTPStatus.BAD_REQUEST, str(exc))
        self._json(HTTPStatus.ACCEPTED, job.snapshot())

    # -- responses ---------------------------------------------------------

    def _json(self, status: HTTPStatus, body: object) -> None:
        self._send(
            status,
            json.dumps(body, separators=(",", ":")).encode(),
            "application/json",
        )

    def _error(self, status: HTTPStatus, message: str) -> None:
        self._json(status, {"error": message})

    def _send(self, status: HTTPStatus, data: bytes, content_type: str) -> None:
        self.send_response(status)
        self.send_header("Content-Type", content_type)
        self.send_header("Content-Length", str(len(data)))
        self.send_header("Cache-Control", "no-store")
        self.send_header("X-Content-Type-Options", "nosniff")
        self.send_header("Referrer-Policy", "no-referrer")
        self.send_header("Content-Security-Policy", _CSP)
        self.end_headers()
        self.wfile.write(data)


class _Server(ThreadingHTTPServer):
    daemon_threads = True
    allow_reuse_address = True

    def __init__(self, address, handler, app: App, verbose: bool = False) -> None:
        super().__init__(address, handler)
        self.app = app
        self.verbose = verbose


def create_server(
    host: str = "127.0.0.1", port: int = 0, *, verbose: bool = False
) -> ThreadingHTTPServer:
    """Create (but do not start) the interface server.

    Args:
        host: Interface to bind; only loopback addresses are accepted.
        port: Port to bind, or 0 to let the system choose a free one.
        verbose: Log each request to stderr.

    Raises:
        ValueError: If ``host`` is not a loopback address.
    """
    if host not in {"127.0.0.1", "localhost"}:
        raise ValueError(
            "the interface only binds to the loopback address (127.0.0.1); "
            "it has no authentication and must not be exposed to a network"
        )
    app = App()
    server = _Server(("127.0.0.1", port), _Handler, app, verbose)
    app.port = server.server_address[1]
    return server


def server_url(server: ThreadingHTTPServer) -> str:
    """Return the address a browser should open."""
    return f"http://127.0.0.1:{server.server_address[1]}/"


def run_gui(
    port: int = 0,
    *,
    window: bool = True,
    open_browser: bool = True,
    verbose: bool = False,
) -> int:
    """Serve the interface until it is closed.

    With ``pywebview`` installed (``pip install pywebview``) and ``window``
    true, the page opens in a native window and the server stops when that
    window closes. Otherwise the address is printed, the page is opened in the
    default browser if ``open_browser`` is true, and the server runs until
    interrupted.
    """
    server = create_server(port=port, verbose=verbose)
    url = server_url(server)

    webview = None
    if window:
        try:
            import webview  # type: ignore[import-not-found,no-redef]
        except ImportError:
            print(
                "pywebview is not installed, so the interface opens in your browser "
                "instead.\nInstall it for a native window: pip install pywebview",
                file=sys.stderr,
            )

    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        if webview is not None:
            try:  # let the Download buttons save files from the native window
                webview.settings["ALLOW_DOWNLOADS"] = True
            except Exception:  # older pywebview without this setting
                pass
            webview.create_window(
                "codonopt", url, width=1180, height=860, min_size=(820, 600)
            )
            webview.start()  # blocks until the window is closed
        else:
            print(f"codonopt interface: {url}   (Ctrl+C to stop)")
            if open_browser:
                webbrowser.open(url)
            while thread.is_alive():
                thread.join(0.5)
    except KeyboardInterrupt:
        pass
    finally:
        server.shutdown()
        server.server_close()
    return 0
