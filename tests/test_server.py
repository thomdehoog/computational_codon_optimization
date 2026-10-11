"""Tests for the local web server.

These drive a real server on a free loopback port over HTTP, because the
properties that matter -- the token and Host checks, the job lifecycle, the
status codes -- live in the request handling rather than in any one function.
"""

from __future__ import annotations

import http.client
import json
import threading
import time

import pytest

from codonopt import translate
from codonopt.server import (
    MAX_RESIDUES,
    RequestError,
    build_request,
    create_server,
)

from conftest import random_protein

TOKEN_HEADER = "X-Codonopt-Token"
PROTEIN = random_protein(40, seed=50)
FAST = {"splice_backend": "none", "effort": "quick", "seed": 1, "seconds": 30}


class Client:
    """A small HTTP client for the test server."""

    def __init__(self, server) -> None:
        self.server = server
        self.port = server.server_address[1]
        self.token = server.app.token

    def call(
        self,
        method: str,
        path: str,
        body=None,
        *,
        token: str | None = "default",
        host: str | None = "default",
        origin: str | None = None,
        raw: bytes | None = None,
        content_length: str | None = None,
    ):
        connection = http.client.HTTPConnection("127.0.0.1", self.port, timeout=30)
        connection.putrequest(method, path, skip_host=True, skip_accept_encoding=True)
        if host is not None:
            connection.putheader(
                "Host", f"127.0.0.1:{self.port}" if host == "default" else host
            )
        if token == "default":
            token = self.token
        if token is not None:
            connection.putheader(TOKEN_HEADER, token)
        if origin is not None:
            connection.putheader("Origin", origin)
        payload = raw if raw is not None else (
            json.dumps(body).encode() if body is not None else b""
        )
        if method == "POST":
            connection.putheader(
                "Content-Length",
                content_length if content_length is not None else str(len(payload)),
            )
            connection.putheader("Content-Type", "application/json")
        connection.endheaders(payload if method == "POST" and content_length is None else b"")
        response = connection.getresponse()
        data = response.read()
        headers = {key.lower(): value for key, value in response.getheaders()}
        connection.close()
        try:
            parsed = json.loads(data)
        except (ValueError, UnicodeDecodeError):
            parsed = data
        return response.status, headers, parsed

    def start(self, payload: dict):
        return self.call("POST", "/api/jobs", payload)

    def wait(self, job_id: str, timeout: float = 90.0) -> dict:
        deadline = time.monotonic() + timeout
        while time.monotonic() < deadline:
            status, _, body = self.call("GET", f"/api/jobs/{job_id}")
            assert status == 200, body
            if body["status"] != "running":
                return body
            time.sleep(0.05)
        raise AssertionError(f"job {job_id} did not finish in {timeout}s")


@pytest.fixture
def client():
    server = create_server()
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    yield Client(server)
    server.shutdown()
    server.server_close()
    thread.join(timeout=5)


# -- static serving and request hygiene --------------------------------------


def test_index_is_served_with_the_session_token_and_a_strict_policy(client):
    status, headers, body = client.call("GET", "/", token=None)
    assert status == 200
    assert headers["content-type"].startswith("text/html")
    assert client.token.encode() in body
    assert b"__CODONOPT_TOKEN__" not in body
    assert "default-src 'none'" in headers["content-security-policy"]
    assert "script-src 'self'" in headers["content-security-policy"]
    assert headers["cache-control"] == "no-store"
    assert headers["x-content-type-options"] == "nosniff"


def test_each_server_has_its_own_token():
    first, second = create_server(), create_server()
    try:
        assert first.app.token != second.app.token
        assert len(first.app.token) >= 24
    finally:
        first.server_close()
        second.server_close()


@pytest.mark.parametrize(
    "path",
    ["/nope", "/../pyproject.toml", "/web/../cli.py", "/%2e%2e/cli.py", "/app.js/../x"],
)
def test_unknown_and_traversal_paths_are_not_found(client, path):
    status, _, _ = client.call("GET", path, token=None)
    assert status == 404


def test_static_paths_reject_other_methods(client):
    status, _, _ = client.call("POST", "/", {})
    assert status == 405


def test_foreign_host_header_is_refused(client):
    """A DNS-rebinding page reaches this server under its own hostname."""
    status, _, _ = client.call("GET", "/", token=None, host="evil.example")
    assert status == 403
    status, _, _ = client.call("GET", "/api/meta", host="evil.example")
    assert status == 403
    status, _, _ = client.call("GET", "/", token=None, host=None)
    assert status == 403


def test_localhost_name_is_accepted(client):
    status, _, _ = client.call("GET", "/api/meta", host=f"localhost:{client.port}")
    assert status == 200


def test_wrong_port_in_host_header_is_refused(client):
    status, _, _ = client.call("GET", "/", token=None, host="127.0.0.1:1")
    assert status == 403


def test_api_requires_the_token(client):
    assert client.call("GET", "/api/meta", token=None)[0] == 403
    assert client.call("GET", "/api/meta", token="wrong")[0] == 403
    assert client.call("POST", "/api/jobs", {"sequence": "MKV"}, token=None)[0] == 403


def test_cross_origin_posts_are_refused(client):
    status, _, body = client.call(
        "POST", "/api/jobs", {"sequence": PROTEIN, **FAST}, origin="https://evil.example"
    )
    assert status == 403
    status, _, body = client.call(
        "POST",
        "/api/jobs",
        {"sequence": PROTEIN, **FAST},
        origin=f"http://127.0.0.1:{client.port}",
    )
    assert status == 202


def test_only_loopback_addresses_can_be_bound():
    for host in ("0.0.0.0", "192.168.1.5", "example.com", ""):
        with pytest.raises(ValueError, match="loopback"):
            create_server(host)


# -- meta --------------------------------------------------------------------


def test_meta_describes_the_form(client):
    status, _, meta = client.call("GET", "/api/meta")
    assert status == 200
    assert {table["name"] for table in meta["tables"]} >= {"human", "ecoli", "uniform"}
    assert "BsaI" in meta["motif_groups"]["golden_gate"]
    assert meta["restriction_sites"]["EcoRI"] == "GAATTC"
    assert meta["defaults"]["gc_range"] == [0.3, 0.7]
    assert set(meta["weights"]) >= {"gc_global", "repeat", "motif", "usage"}
    assert set(meta["effort"]) == {"quick", "standard", "thorough"}
    assert isinstance(meta["maxent_available"], bool)
    assert meta["example"]["protein"].startswith("M")
    assert meta["limits"]["max_residues"] == MAX_RESIDUES


# -- job lifecycle -----------------------------------------------------------


def test_a_job_runs_to_completion_and_keeps_the_protein(client):
    status, _, started = client.start(
        {"sequence": PROTEIN, "avoid": ["golden_gate", "EcoRI"], **FAST}
    )
    assert status == 202
    # The worker may finish before the HTTP start response is serialized.
    assert started["status"] in {"running", "done"}

    done = client.wait(started["id"])
    assert done["status"] == "done"
    result = done["result"]
    assert translate(result["sequence"]) == PROTEIN
    assert result["protein"] == PROTEIN
    assert translate(result["initial_sequence"]) == PROTEIN
    assert result["score"]["final"] <= result["score"]["initial"]
    assert result["final"]["forbidden_motif_hits"] == 0
    assert result["fasta"].startswith(">sequence ")
    assert result["fasta"].replace("\n", "").endswith(result["sequence"])
    assert result["cancelled"] is False


def test_the_gc_profile_lines_up_for_plotting(client):
    _, _, started = client.start({"sequence": PROTEIN, **FAST})
    gc = client.wait(started["id"])["result"]["gc"]
    assert len(gc["starts"]) == len(gc["initial"]) == len(gc["final"]) > 0
    assert gc["range"] == [0.3, 0.7]
    assert all(0.0 <= value <= 1.0 for value in gc["final"])
    assert gc["starts"] == sorted(gc["starts"])


def test_progress_is_reported_while_running(client):
    _, _, started = client.start(
        {"sequence": random_protein(120, seed=51), **FAST, "effort": "thorough"}
    )
    snapshot = client.call("GET", f"/api/jobs/{started['id']}")[2]
    assert snapshot["restarts_total"] == 3
    assert snapshot["max_seconds"] == 30
    done = client.wait(started["id"])
    assert done["restarts_done"] == 3
    assert any("final score" in message for message in done["messages"])


def test_dna_is_detected_and_improved_in_place(client):
    dna = "ATGAAATGGGTTACCTTTATCTCTCTGCTGTTG"
    _, _, started = client.start({"sequence": f">gene one\n{dna}\n", **FAST})
    result = client.wait(started["id"])["result"]
    assert result["kind"] == "dna"
    assert result["name"] == "gene"
    assert result["initial_sequence"] == dna
    assert translate(result["sequence"]) == translate(dna)


def test_a_trailing_stop_codon_can_be_requested(client):
    _, _, started = client.start({"sequence": "MKWVTF", "add_stop": True, **FAST})
    result = client.wait(started["id"])["result"]
    assert result["protein"] == "MKWVTF*"
    assert translate(result["sequence"]) == "MKWVTF*"


def test_a_custom_codon_table_can_be_pasted(client):
    table = json.dumps({"CTG": 1, "CTC": 0.1})
    _, _, started = client.start(
        {"sequence": "MLLLLL", "table_text": table, **FAST}
    )
    result = client.wait(started["id"])["result"]
    assert result["codon_usage_table"] == "custom"


def test_remaining_motif_hits_are_reported_for_highlighting(client):
    # Weight the motif term to zero so a site is allowed to survive.
    _, _, started = client.start(
        {
            "sequence": "MKWVTF",
            "avoid": ["x=ATG"],
            "weights": {"motif": 0},
            **FAST,
        }
    )
    hits = client.wait(started["id"])["result"]["motif_hits"]
    assert hits and hits[0]["name"] == "x" and hits[0]["length"] == 3
    assert hits[0]["positions"][0] == 0


def test_cancelling_returns_the_best_sequence_so_far(client):
    protein = random_protein(500, seed=52)
    _, _, started = client.start(
        {"sequence": protein, **FAST, "effort": "thorough", "seconds": 300}
    )
    begun = time.monotonic()
    status, _, _ = client.call("POST", f"/api/jobs/{started['id']}/cancel")
    assert status == 200
    done = client.wait(started["id"], timeout=60)
    assert time.monotonic() - begun < 30  # nowhere near the 300 s budget
    assert done["status"] == "cancelled"
    assert translate(done["result"]["sequence"]) == protein
    assert done["result"]["cancelled"] is True


def test_the_concurrency_limit_is_enforced_and_recovers(client):
    protein = random_protein(600, seed=53)
    slow = {"sequence": protein, **FAST, "effort": "thorough", "seconds": 300}
    ids = []
    for _ in range(2):
        status, _, body = client.start(slow)
        assert status == 202
        ids.append(body["id"])
    status, _, body = client.start(slow)
    assert status == 400
    assert "as many optimisations as allowed" in body["error"]

    for job_id in ids:
        client.call("POST", f"/api/jobs/{job_id}/cancel")
        client.wait(job_id, timeout=60)
    assert client.start({"sequence": PROTEIN, **FAST})[0] == 202


def test_unknown_jobs_are_not_found(client):
    assert client.call("GET", "/api/jobs/nonexistent")[0] == 404
    assert client.call("POST", "/api/jobs/nonexistent/cancel")[0] == 404
    assert client.call("GET", "/api/elsewhere")[0] == 404


# -- bad requests ------------------------------------------------------------


@pytest.mark.parametrize(
    "payload, message",
    [
        ({"sequence": ""}, "Paste a protein or DNA sequence"),
        ({"sequence": "   \n "}, "Paste a protein or DNA sequence"),
        ({"sequence": "MKZ"}, "unsupported residues"),
        ({"sequence": "MK*V"}, "internal stop"),
        ({"sequence": "ATGAAATGGG"}, "multiple of three"),
        ({"sequence": ">a\nMKV\n>b\nMKW\n"}, "2 FASTA records"),
        ({"sequence": "MKV", "gc_range": [30, 70]}, "between 0 and 1"),
        ({"sequence": "MKV", "gc_range": [0.7, 0.3]}, "low <= high"),
        ({"sequence": "MKV", "gc_range": [0.3]}, "two numbers"),
        ({"sequence": "MKV", "gc_target": 2}, "between 0 and 1"),
        ({"sequence": "MKV", "gc_window": 2.5}, "whole number"),
        ({"sequence": "MKV", "gc_window": "wide"}, "must be a number"),
        ({"sequence": "MKV", "table": "mouse"}, "unknown codon usage table"),
        ({"sequence": "MKV", "table_text": "TTT,abc"}, "not a number"),
        ({"sequence": "MKV", "table_text": "nothing useful"}, "expected two columns"),
        ({"sequence": "MKV", "table_text": "codon,count\n"}, "no codon rows"),
        ({"sequence": "MKV", "avoid": ["NotAnEnzyme"]}, "cannot resolve motif"),
        ({"sequence": "MKV", "avoid": "EcoRI"}, "list of motif names"),
        ({"sequence": "MKV", "weights": {"nonsense": 1}}, "unknown weight"),
        ({"sequence": "MKV", "weights": {"usage": -1}}, "between 0 and"),
        ({"sequence": "MKV", "weights": []}, "must be an object"),
        ({"sequence": "MKV", "effort": "extreme"}, "effort must be"),
        ({"sequence": "MKV", "splice_backend": "magic"}, "splice backend must be"),
        ({"sequence": "MKV", "seconds": 0}, "between 1 and"),
        ({"sequence": "MKV", "seconds": 10**6}, "between 1 and"),
        ({"sequence": "MKV", "seed": -1}, "between 0 and"),
        # An integer too wide for a float used to escape as a 500.
        ({"sequence": "MKV", "seconds": 10**400}, "between 1 and"),
        ({"sequence": "MKV", "gc_target": -(10**400)}, "between 0 and"),
        # "Custom table" with nothing in it must not fall back to another host.
        ({"sequence": "MKV", "table_text": ""}, "custom codon table is empty"),
        ({"sequence": "MKV", "table_text": "   "}, "custom codon table is empty"),
        ({"sequence": "MKV", "kind": "rna"}, "kind must be"),
        ({"sequence": "A" * (MAX_RESIDUES + 1), "kind": "protein"}, "accepts up to"),
    ],
)
def test_bad_requests_get_a_readable_400(client, payload, message):
    status, _, body = client.call("POST", "/api/jobs", payload)
    assert status == 400, body
    assert message in body["error"]


def test_malformed_bodies(client):
    assert client.call("POST", "/api/jobs", raw=b"{not json")[0] == 400
    assert client.call("POST", "/api/jobs", raw=b"[1, 2]")[0] == 400
    assert client.call("POST", "/api/jobs", raw=b"\xff\xfe")[0] == 400
    # More digits than CPython will parse into an int: a bad request, not a crash.
    huge = b'{"sequence": "MKV", "seconds": ' + b"9" * 5000 + b"}"
    status, _, body = client.call("POST", "/api/jobs", raw=huge)
    assert status == 400
    assert "not valid JSON" in body["error"]


def test_an_explicit_table_still_wins_over_an_empty_custom_one(client):
    """The page omits `table` for a custom table; other callers may send both."""
    status, _, started = client.start(
        {"sequence": "MKV", "table": "ecoli", "table_text": "", **FAST}
    )
    assert status == 202
    assert client.wait(started["id"])["status"] == "done"
    status, _, body = client.call("POST", "/api/jobs", raw=b"null")
    assert status == 400


def test_oversized_and_unsized_bodies_are_refused(client):
    status, _, _ = client.call(
        "POST", "/api/jobs", raw=b"", content_length=str(3_000_000)
    )
    assert status == 413
    status, _, _ = client.call("POST", "/api/jobs", raw=b"", content_length="lots")
    assert status == 411


# -- build_request, without HTTP --------------------------------------------


def test_build_request_interprets_the_input():
    spec = build_request({"sequence": "MKWVTF", "add_stop": True})
    assert (spec.protein, spec.start, spec.kind) == ("MKWVTF*", None, "protein")

    spec = build_request({"sequence": "atgaaatgggtt"})
    assert spec.kind == "dna" and spec.start == "ATGAAATGGGTT"
    assert spec.protein == "MKWV"

    # RNA spelling is accepted and normalised.
    spec = build_request({"sequence": "AUGAAAUGGGUU"})
    assert spec.start == "ATGAAATGGGTT"


def test_add_stop_does_not_alter_a_dna_input():
    spec = build_request({"sequence": "ATGAAATGGGTT", "add_stop": True})
    assert spec.protein == "MKWV"


def test_kind_overrides_detection():
    # Alanine, cysteine, glycine and threonine spell valid DNA letters, so a
    # sequence of only those is ambiguous and the caller may need to say.
    assert build_request({"sequence": "ACGTACGTA"}).kind == "dna"
    spec = build_request({"sequence": "ACGTACGTA", "kind": "protein"})
    assert spec.kind == "protein" and spec.start is None
    assert spec.protein == "ACGTACGTA"
    # Short strings are read as protein unless told otherwise.
    assert build_request({"sequence": "ATGAAATG"}).kind == "protein"


def test_build_request_names_the_record():
    assert build_request({"sequence": ">albumin sample\nMKV"}).name == "albumin"
    assert build_request({"sequence": "MKV", "name": "mine"}).name == "mine"
    assert build_request({"sequence": "MKV"}).name == "sequence"


def test_build_request_resolves_every_option():
    spec = build_request(
        {
            "sequence": "MKV",
            "table": "ecoli",
            "avoid": ["golden_gate", "AATAAA"],
            "gc_target": 0.6,
            "gc_range": [0.35, 0.65],
            "repeat_k": 9,
            "weights": {"usage": 0, "motif": 10},
            "effort": "thorough",
            "seed": 7,
            "seconds": 45,
            "splice_backend": "motif",
        }
    )
    config = spec.config
    assert config.table.name == "ecoli"
    assert config.gc_target == 0.6 and config.gc_range == (0.35, 0.65)
    assert config.repeat_k == 9
    assert config.weights.usage == 0 and config.weights.motif == 10
    assert len(config.motifs) == 5  # four golden-gate enzymes plus AATAAA
    assert config.splice_scorer.name == "motif"
    assert spec.settings["restarts"] == 3 and spec.settings["seed"] == 7
    assert spec.settings["max_seconds"] == 45


def test_request_errors_are_value_errors():
    assert issubclass(RequestError, ValueError)
    with pytest.raises(RequestError):
        build_request("not an object")


def test_request_host_policy_and_context_options():
    spec = build_request({"sequence": "TTTCTCTTT", "kind": "dna", "table": "ecoli",
                          "splice_backend": "auto", "left_flank": "GGT",
                          "right_flank": "A", "protected_codons": [0],
                          "excluded_codons": ["TTA"], "required": ["motif", "gc_window"]})
    assert spec.config.splice_scorer is None
    assert spec.config.splice_strand == "forward"
    assert spec.config.left_flank == "GGT" and spec.config.right_flank == "A"
    assert spec.config.protected_codons == (0,)
    assert spec.config.required == ("gc_window", "motif")
    assert build_request({"sequence": "MKV"}).config.splice_scorer is None
    for key, value in (("required", "motif"), ("protected_codons", [True]),
                       ("left_flank", []), ("splice_strand", "unknown"), ("splice_strand", [])):
        with pytest.raises(RequestError):
            build_request({"sequence": "MKV", key: value})


# ---------------------------------------------------------------- built assets


def test_built_interface_is_shipped_and_wired_to_the_token():
    """The Node-free install serves committed files; make sure they exist and fit together."""
    from importlib import resources

    web = resources.files("codonopt") / "web"
    names = {p.name for p in web.iterdir()}
    assert {"index.html", "app.js", "app.css", "favicon.svg"} <= names
    html = (web / "index.html").read_text(encoding="utf-8")
    assert "__CODONOPT_TOKEN__" in html
    assert "app.js" in html and "app.css" in html
    assert "<script>" not in html  # the content security policy forbids inline script
