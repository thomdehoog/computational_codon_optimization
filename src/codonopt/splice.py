"""Cryptic splice site scoring.

Two interchangeable backends are provided:

``maxent``
    MaxEntScan scoring (Yeo & Burge 2004, Journal of Computational Biology
    11:377-394) via the optional ``maxentpy`` package. Scoring matrices are
    loaded once per process and reused; the original scripts reloaded both
    matrices on every call, which cost about 50 ms per call and dominated the
    runtime of the whole pipeline.

``motif``
    A dependency-free consensus screen. It counts matches to the donor
    consensus and the pyrimidine content of the acceptor tract. It is a coarse
    filter, not a model of splicing, and its scores are not comparable to
    MaxEntScan scores -- each backend carries its own default cutoffs.

Both backends use the MaxEntScan window conventions, so a reported position
means the same thing either way:

* donor (5') window: 9 bases, 3 exonic then 6 intronic, with ``GT`` at window
  offsets 3 and 4;
* acceptor (3') window: 23 bases, 20 intronic then 3 exonic, with ``AG`` at
  window offsets 18 and 19.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from typing import Literal

from .sequence import reverse_complement

__all__ = [
    "DONOR_WINDOW",
    "ACCEPTOR_WINDOW",
    "SpliceSite",
    "SpliceScorer",
    "MaxEntScorer",
    "MotifScorer",
    "get_scorer",
    "maxent_available",
    "find_splice_sites",
]

DONOR_WINDOW = 9
ACCEPTOR_WINDOW = 23

Strand = Literal["forward", "reverse"]
SiteType = Literal["donor", "acceptor"]


@dataclass(frozen=True)
class SpliceSite:
    """A scored candidate splice site.

    Attributes:
        position: 0-based start of the scoring window in forward-strand
            coordinates.
        score: Backend-specific score; higher means a stronger site.
        site_type: ``"donor"`` or ``"acceptor"``.
        strand: Which strand the window was scored on.
    """

    position: int
    score: float
    site_type: SiteType
    strand: Strand


class SpliceScorer:
    """Base class for splice site scorers."""

    name = "base"
    #: Sensible starting cutoffs for this backend's score scale. A cutoff is
    #: the highest score still treated as acceptable: a window is reported and
    #: penalised only when it scores strictly above it.
    default_donor_cutoff = 0.0
    default_acceptor_cutoff = 0.0

    def score_donor(self, window: str) -> float:
        """Score a 9-base donor window."""
        raise NotImplementedError

    def score_acceptor(self, window: str) -> float:
        """Score a 23-base acceptor window."""
        raise NotImplementedError


class MaxEntScorer(SpliceScorer):
    """MaxEntScan scoring backed by the optional ``maxentpy`` package."""

    name = "maxent"
    default_donor_cutoff = 3.0
    default_acceptor_cutoff = 3.0

    def __init__(self) -> None:
        self._score5, self._score3, self._matrix5, self._matrix3 = _load_maxent()

    def score_donor(self, window: str) -> float:
        return float(self._score5(window, matrix=self._matrix5))

    def score_acceptor(self, window: str) -> float:
        return float(self._score3(window, matrix=self._matrix3))


class MotifScorer(SpliceScorer):
    """Dependency-free consensus screen; see the module docstring."""

    name = "motif"
    # Scores are small integers here, so the cutoff is the highest acceptable
    # one: four or more consensus matches beyond the mandatory GT is a donor
    # worth designing out.
    default_donor_cutoff = 3.0
    # Likewise 14 of 18 pyrimidines is already a credible polypyrimidine tract.
    default_acceptor_cutoff = 13.0

    _NO_SITE = -10.0
    # Informative donor consensus positions, excluding the mandatory GT.
    _DONOR_CONSENSUS: tuple[tuple[int, str], ...] = (
        (1, "A"),
        (2, "G"),
        (5, "AG"),
        (6, "A"),
        (7, "G"),
        (8, "T"),
    )
    _PYRIMIDINES = frozenset("CT")

    def score_donor(self, window: str) -> float:
        if window[3:5] != "GT":
            return self._NO_SITE
        return float(
            sum(1 for offset, allowed in self._DONOR_CONSENSUS if window[offset] in allowed)
        )

    def score_acceptor(self, window: str) -> float:
        if window[18:20] != "AG":
            return self._NO_SITE
        return float(sum(1 for base in window[:18] if base in self._PYRIMIDINES))


class _MissingMaxEnt(RuntimeError):
    """Raised when the maxent backend is requested but unavailable."""


@lru_cache(maxsize=1)
def _load_maxent():
    """Import ``maxentpy`` and load both matrices once per process."""
    try:
        from maxentpy import maxent
        from maxentpy.maxent import load_matrix3, load_matrix5
    except ImportError as exc:  # pragma: no cover - depends on environment
        raise _MissingMaxEnt(
            "the 'maxent' splice backend needs the maxentpy package:\n"
            "    pip install 'git+https://github.com/kepbod/maxentpy.git'\n"
            "or choose --splice-backend motif / none"
        ) from exc
    return maxent.score5, maxent.score3, load_matrix5(), load_matrix3()


def maxent_available() -> bool:
    """Return True if the ``maxent`` backend can be used."""
    try:
        _load_maxent()
    except _MissingMaxEnt:
        return False
    return True


def get_scorer(backend: str) -> SpliceScorer | None:
    """Return a scorer for ``backend``.

    Args:
        backend: ``"maxent"``, ``"motif"``, ``"none"``, or ``"auto"`` to prefer
            MaxEntScan and fall back to the motif screen.

    Returns:
        A scorer, or None when splice scoring is switched off.

    Raises:
        ValueError: On an unknown backend name.
        RuntimeError: If ``maxent`` is requested but ``maxentpy`` is missing.
    """
    key = backend.strip().lower()
    if key in {"none", "off", ""}:
        return None
    if key == "maxent":
        return MaxEntScorer()
    if key == "motif":
        return MotifScorer()
    if key == "auto":
        return MaxEntScorer() if maxent_available() else MotifScorer()
    raise ValueError(
        f"unknown splice backend {backend!r}; expected maxent, motif, none or auto"
    )


def find_splice_sites(
    sequence: str,
    scorer: SpliceScorer,
    *,
    donor_cutoff: float | None = None,
    acceptor_cutoff: float | None = None,
    strand: str = "both",
) -> list[SpliceSite]:
    """Return every window scoring strictly above its cutoff.

    The comparison matches the one the objective uses, so a site this function
    reports is a site the optimizer is paid to remove.

    Reverse-strand windows are scored on the reverse complement of the forward
    window, so a returned position is always a forward-strand offset and stays
    comparable between strands.

    Raises:
        ValueError: On an unknown ``strand``.
    """
    if strand not in {"forward", "reverse", "both"}:
        raise ValueError(f"unknown strand {strand!r}")
    donor_cutoff = (
        scorer.default_donor_cutoff if donor_cutoff is None else donor_cutoff
    )
    acceptor_cutoff = (
        scorer.default_acceptor_cutoff if acceptor_cutoff is None else acceptor_cutoff
    )

    sites: list[SpliceSite] = []
    strands: tuple[Strand, ...] = (
        ("forward", "reverse") if strand == "both" else (strand,)  # type: ignore[assignment]
    )
    for current in strands:
        flip = current == "reverse"
        for start in range(len(sequence) - DONOR_WINDOW + 1):
            window = sequence[start : start + DONOR_WINDOW]
            if flip:
                window = reverse_complement(window)
            score = scorer.score_donor(window)
            if score > donor_cutoff:
                sites.append(SpliceSite(start, score, "donor", current))
        for start in range(len(sequence) - ACCEPTOR_WINDOW + 1):
            window = sequence[start : start + ACCEPTOR_WINDOW]
            if flip:
                window = reverse_complement(window)
            score = scorer.score_acceptor(window)
            if score > acceptor_cutoff:
                sites.append(SpliceSite(start, score, "acceptor", current))
    sites.sort(key=lambda site: (site.position, site.site_type, site.strand))
    return sites
