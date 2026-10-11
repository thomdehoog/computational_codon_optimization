"""Named sequence motifs to keep out of a synthetic construct.

Restriction sites are listed so a construct stays compatible with a cloning
strategy, and the regulatory motifs are the ones that commonly cause trouble
when a coding sequence is expressed from a transgene. Patterns may use IUPAC
ambiguity codes and are matched on both strands.
"""

from __future__ import annotations

import re

__all__ = [
    "RESTRICTION_SITES",
    "REGULATORY_MOTIFS",
    "MOTIF_GROUPS",
    "iupac_to_regex",
    "resolve_motifs",
    "Motif",
]

_IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "[AG]", "Y": "[CT]", "S": "[GC]", "W": "[AT]",
    "K": "[GT]", "M": "[AC]", "B": "[CGT]", "D": "[AGT]",
    "H": "[ACT]", "V": "[ACG]", "N": "[ACGT]",
}

#: Recognition sequences of commonly used restriction enzymes.
RESTRICTION_SITES: dict[str, str] = {
    "AatII": "GACGTC",
    "AgeI": "ACCGGT",
    "ApaI": "GGGCCC",
    "AscI": "GGCGCGCC",
    "AvrII": "CCTAGG",
    "BamHI": "GGATCC",
    "BbsI": "GAAGAC",
    "BglII": "AGATCT",
    "BsaI": "GGTCTC",
    "BsmBI": "CGTCTC",
    "BssHII": "GCGCGC",
    "ClaI": "ATCGAT",
    "DraI": "TTTAAA",
    "EcoRI": "GAATTC",
    "EcoRV": "GATATC",
    "FseI": "GGCCGGCC",
    "HindIII": "AAGCTT",
    "HpaI": "GTTAAC",
    "KpnI": "GGTACC",
    "MluI": "ACGCGT",
    "NarI": "GGCGCC",
    "NcoI": "CCATGG",
    "NdeI": "CATATG",
    "NheI": "GCTAGC",
    "NotI": "GCGGCCGC",
    "PacI": "TTAATTAA",
    "PstI": "CTGCAG",
    "PvuI": "CGATCG",
    "PvuII": "CAGCTG",
    "SacI": "GAGCTC",
    "SalI": "GTCGAC",
    "SapI": "GCTCTTC",
    "SbfI": "CCTGCAGG",
    "ScaI": "AGTACT",
    "SmaI": "CCCGGG",
    "SpeI": "ACTAGT",
    "SphI": "GCATGC",
    "SspI": "AATATT",
    "StuI": "AGGCCT",
    "XbaI": "TCTAGA",
    "XhoI": "CTCGAG",
}

#: Regulatory motifs worth avoiding inside a coding sequence.
REGULATORY_MOTIFS: dict[str, str] = {
    "polyA_AATAAA": "AATAAA",
    "polyA_ATTAAA": "ATTAAA",
    "TATA_box": "TATAAA",
    "Kozak_internal": "GCCACCATGG",
    "chi_site": "GCTGGTGG",
}

#: Convenience bundles selectable by name.
MOTIF_GROUPS: dict[str, tuple[str, ...]] = {
    "golden_gate": ("BsaI", "BsmBI", "BbsI", "SapI"),
    "common_cloning": (
        "EcoRI", "BamHI", "HindIII", "XhoI", "SalI", "XbaI", "SpeI",
        "NotI", "NcoI", "NdeI", "KpnI", "SacI", "PstI",
    ),
    "regulatory": tuple(REGULATORY_MOTIFS),
    "none": (),
}


class Motif:
    """A named pattern, compiled for both-strand matching.

    Attributes:
        name: Display name, e.g. an enzyme name.
        pattern: The IUPAC pattern as given.
    """

    __slots__ = ("name", "pattern", "_regexes", "length")

    def __init__(self, name: str, pattern: str) -> None:
        pattern = pattern.strip().upper().replace("U", "T")
        if not pattern:
            raise ValueError(f"motif {name!r} has an empty pattern")
        self.name = name
        self.pattern = pattern
        self.length = len(pattern)
        forward = iupac_to_regex(pattern)
        reverse = iupac_to_regex(_revcomp_iupac(pattern))
        patterns = [forward] if reverse == forward else [forward, reverse]
        # Overlapping matches matter, so each regex is applied with a lookahead.
        self._regexes = tuple(re.compile(f"(?=({p}))") for p in patterns)

    def find_all(self, sequence: str) -> list[int]:
        """Return the 0-based start positions of every match on either strand.

        A position is reported once even when both strands match there, which
        an ambiguous pattern can do: one position is one site.
        """
        hits: set[int] = set()
        for regex in self._regexes:
            hits.update(match.start() for match in regex.finditer(sequence))
        return sorted(hits)

    def matches_at(self, sequence: str, start: int) -> bool:
        """Return True if a match starts exactly at ``start``."""
        return any(regex.match(sequence, start) for regex in self._regexes)

    def __repr__(self) -> str:  # pragma: no cover - debugging aid
        return f"Motif({self.name!r}, {self.pattern!r})"


def canonicalize_motifs(motifs) -> tuple[Motif, ...]:
    """Deduplicate both-strand patterns and reject ambiguous report names."""
    patterns: dict[str, Motif] = {}
    names: dict[str, str] = {}
    for motif in motifs:
        key = min(motif.pattern, _revcomp_iupac(motif.pattern))
        if motif.name in names and names[motif.name] != key:
            raise ValueError(f"motif name {motif.name!r} refers to different patterns")
        names[motif.name] = key
        patterns.setdefault(key, motif)
    return tuple(patterns.values())


_IUPAC_COMPLEMENT = str.maketrans(
    "ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN"
)


def _revcomp_iupac(pattern: str) -> str:
    """Reverse complement an IUPAC pattern."""
    return pattern.translate(_IUPAC_COMPLEMENT)[::-1]


def iupac_to_regex(pattern: str) -> str:
    """Translate an IUPAC nucleotide pattern into a regular expression.

    Raises:
        ValueError: On a character that is not an IUPAC nucleotide code.
    """
    try:
        return "".join(_IUPAC[base] for base in pattern)
    except KeyError as exc:
        raise ValueError(
            f"{exc.args[0]!r} is not an IUPAC nucleotide code"
        ) from None


def resolve_motifs(specs: list[str] | tuple[str, ...] | None) -> list[Motif]:
    """Resolve motif specifications into :class:`Motif` objects.

    Each specification may be a group name from :data:`MOTIF_GROUPS`, an enzyme
    name from :data:`RESTRICTION_SITES`, a motif name from
    :data:`REGULATORY_MOTIFS`, a bare IUPAC pattern, or ``name=PATTERN``.

    Duplicate patterns are collapsed, and so are two patterns that are reverse
    complements of each other: both strands are searched anyway, so ``GAAGAC``
    and ``GTCTTC`` name one site and are counted once rather than twice.

    Raises:
        ValueError: If a specification cannot be resolved.
    """
    resolved: dict[str, Motif] = {}

    def add(name: str, pattern: str) -> None:
        motif = Motif(name, pattern)
        key = min(motif.pattern, _revcomp_iupac(motif.pattern))
        if any(m.name == name and k != key for k, m in resolved.items()):
            raise ValueError(f"motif name {name!r} refers to different patterns")
        resolved.setdefault(key, motif)

    for raw in specs or ():
        spec = raw.strip()
        if not spec:
            continue
        if "=" in spec:
            name, _, pattern = spec.partition("=")
            add(name.strip() or pattern.strip(), pattern)
            continue
        key = spec.lower()
        if key in MOTIF_GROUPS:
            for member in MOTIF_GROUPS[key]:
                add(member, _lookup_named(member))
            continue
        named = _lookup_named(spec, required=False)
        if named is not None:
            add(_canonical_name(spec), named)
            continue
        try:
            iupac_to_regex(spec.upper())
        except ValueError:
            raise ValueError(
                f"cannot resolve motif {raw!r}: not a group, a known site, "
                "or an IUPAC pattern"
            ) from None
        add(spec.upper(), spec.upper())

    return list(resolved.values())


def _canonical_name(name: str) -> str:
    lowered = name.lower()
    for table in (RESTRICTION_SITES, REGULATORY_MOTIFS):
        for key in table:
            if key.lower() == lowered:
                return key
    return name


def _lookup_named(name: str, *, required: bool = True) -> str | None:
    lowered = name.lower()
    for table in (RESTRICTION_SITES, REGULATORY_MOTIFS):
        for key, pattern in table.items():
            if key.lower() == lowered:
                return pattern
    if required:
        raise ValueError(f"unknown motif name: {name!r}")
    return None
