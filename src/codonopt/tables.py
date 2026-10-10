"""Codon usage tables and codon adaptation scoring.

The original scripts picked synonymous codons uniformly at random, which throws
away the single most useful piece of prior information in codon optimisation:
not all synonymous codons are translated equally well in a given host. A
:class:`CodonUsage` table carries a relative frequency per codon and derives the
relative adaptiveness weights used by the Codon Adaptation Index (CAI)
of Sharp & Li (1987), Nucleic Acids Research 15:1281-1295.

The bundled ``human`` and ``ecoli`` tables are *approximate* relative
synonymous codon frequencies, rounded to two decimals, suitable for ranking
codons. For production design work, supply a table measured on your own
expression host with :func:`load_table`.
"""

from __future__ import annotations

import csv
import io
import json
import math
import random
from dataclasses import dataclass
from pathlib import Path

from .sequence import AA_TO_CODONS, CODON_TO_AA, STOP, split_codons

__all__ = [
    "CodonUsage",
    "TABLES",
    "available_tables",
    "get_table",
    "load_table",
    "table_from_text",
]

# Approximate relative synonymous codon usage, Homo sapiens.
_HUMAN: dict[str, float] = {
    "TTT": 0.45, "TTC": 0.55,
    "TTA": 0.07, "TTG": 0.13, "CTT": 0.13, "CTC": 0.20, "CTA": 0.07, "CTG": 0.40,
    "ATT": 0.36, "ATC": 0.48, "ATA": 0.16,
    "ATG": 1.00,
    "GTT": 0.18, "GTC": 0.24, "GTA": 0.11, "GTG": 0.47,
    "TCT": 0.18, "TCC": 0.22, "TCA": 0.15, "TCG": 0.06, "AGT": 0.15, "AGC": 0.24,
    "CCT": 0.28, "CCC": 0.33, "CCA": 0.27, "CCG": 0.11,
    "ACT": 0.24, "ACC": 0.36, "ACA": 0.28, "ACG": 0.12,
    "GCT": 0.26, "GCC": 0.40, "GCA": 0.23, "GCG": 0.11,
    "TAT": 0.43, "TAC": 0.57,
    "CAT": 0.41, "CAC": 0.59,
    "CAA": 0.25, "CAG": 0.75,
    "AAT": 0.46, "AAC": 0.54,
    "AAA": 0.42, "AAG": 0.58,
    "GAT": 0.46, "GAC": 0.54,
    "GAA": 0.42, "GAG": 0.58,
    "TGT": 0.45, "TGC": 0.55,
    "TGG": 1.00,
    "CGT": 0.08, "CGC": 0.19, "CGA": 0.11, "CGG": 0.21, "AGA": 0.20, "AGG": 0.21,
    "GGT": 0.16, "GGC": 0.34, "GGA": 0.25, "GGG": 0.25,
    "TAA": 0.30, "TAG": 0.24, "TGA": 0.46,
}

# Approximate relative synonymous codon usage, Escherichia coli K-12.
_ECOLI: dict[str, float] = {
    "TTT": 0.57, "TTC": 0.43,
    "TTA": 0.13, "TTG": 0.13, "CTT": 0.10, "CTC": 0.10, "CTA": 0.04, "CTG": 0.50,
    "ATT": 0.51, "ATC": 0.42, "ATA": 0.07,
    "ATG": 1.00,
    "GTT": 0.26, "GTC": 0.22, "GTA": 0.15, "GTG": 0.37,
    "TCT": 0.15, "TCC": 0.15, "TCA": 0.12, "TCG": 0.15, "AGT": 0.15, "AGC": 0.28,
    "CCT": 0.16, "CCC": 0.12, "CCA": 0.19, "CCG": 0.53,
    "ACT": 0.17, "ACC": 0.44, "ACA": 0.13, "ACG": 0.26,
    "GCT": 0.16, "GCC": 0.27, "GCA": 0.21, "GCG": 0.36,
    "TAT": 0.57, "TAC": 0.43,
    "CAT": 0.57, "CAC": 0.43,
    "CAA": 0.35, "CAG": 0.65,
    "AAT": 0.45, "AAC": 0.55,
    "AAA": 0.76, "AAG": 0.24,
    "GAT": 0.63, "GAC": 0.37,
    "GAA": 0.68, "GAG": 0.32,
    "TGT": 0.45, "TGC": 0.55,
    "TGG": 1.00,
    "CGT": 0.38, "CGC": 0.40, "CGA": 0.06, "CGG": 0.10, "AGA": 0.04, "AGG": 0.02,
    "GGT": 0.34, "GGC": 0.40, "GGA": 0.11, "GGG": 0.15,
}


@dataclass(frozen=True)
class CodonUsage:
    """Relative synonymous codon frequencies for one expression host.

    Attributes:
        name: Identifier used in reports and on the command line.
        frequencies: Codon -> relative frequency within its amino acid family.
            Families are normalised to sum to 1 on construction.
        weights: Codon -> relative adaptiveness ``w = f / f_max`` within the
            family, as used by the Codon Adaptation Index.
        description: Free-text provenance note.
    """

    name: str
    frequencies: dict[str, float]
    weights: dict[str, float]
    description: str = ""

    @classmethod
    def from_frequencies(
        cls,
        name: str,
        frequencies: dict[str, float],
        *,
        description: str = "",
    ) -> "CodonUsage":
        """Build a table from raw per-codon values (counts or frequencies).

        Values are normalised within each amino acid family, so either raw
        codon counts or already-normalised frequencies may be supplied. Codons
        missing from ``frequencies`` default to 0 and are never chosen, unless
        every codon in a family is missing, in which case the family falls back
        to a uniform distribution.

        Raises:
            ValueError: If a key is not a codon of the standard genetic code or
                a value is negative.
        """
        raw = {}
        for codon, value in frequencies.items():
            codon = codon.strip().upper().replace("U", "T")
            if codon not in CODON_TO_AA:
                raise ValueError(f"not a codon of the standard genetic code: {codon!r}")
            if value < 0:
                raise ValueError(f"negative frequency for codon {codon}: {value}")
            raw[codon] = float(value)

        normalised: dict[str, float] = {}
        weights: dict[str, float] = {}
        for aa, codons in AA_TO_CODONS.items():
            values = [raw.get(codon, 0.0) for codon in codons]
            total = sum(values)
            if total <= 0:
                values = [1.0] * len(codons)
                total = float(len(codons))
            peak = max(values)
            for codon, value in zip(codons, values):
                normalised[codon] = value / total
                weights[codon] = value / peak
        return cls(
            name=name,
            frequencies=normalised,
            weights=weights,
            description=description,
        )

    def codons(self, amino_acid: str) -> tuple[str, ...]:
        """Return the synonymous codons for ``amino_acid``."""
        return AA_TO_CODONS[amino_acid]

    def frequency(self, codon: str) -> float:
        """Return the within-family relative frequency of ``codon``."""
        return self.frequencies[codon]

    def weight(self, codon: str) -> float:
        """Return the relative adaptiveness of ``codon`` (``f / f_max``)."""
        return self.weights[codon]

    def most_frequent(self, amino_acid: str) -> str:
        """Return the single most frequent codon for ``amino_acid``."""
        return max(AA_TO_CODONS[amino_acid], key=self.frequencies.__getitem__)

    def sample(self, amino_acid: str, rng: random.Random) -> str:
        """Draw a codon for ``amino_acid`` in proportion to its frequency."""
        codons = AA_TO_CODONS[amino_acid]
        if len(codons) == 1:
            return codons[0]
        weights = [self.frequencies[codon] for codon in codons]
        return rng.choices(codons, weights=weights, k=1)[0]

    def cai(self, sequence: str) -> float:
        """Return the Codon Adaptation Index of an in-frame ``sequence``.

        Single-codon families (Met, Trp) and stop codons carry no information
        about codon choice and are excluded, following the usual convention.
        Returns 0.0 if no scorable codon is present.
        """
        logs = []
        for codon in split_codons(sequence.upper()):
            aa = CODON_TO_AA[codon]
            if aa == STOP or len(AA_TO_CODONS[aa]) == 1:
                continue
            weight = self.weights[codon]
            logs.append(math.log(weight) if weight > 0 else math.log(1e-6))
        if not logs:
            return 0.0
        return math.exp(sum(logs) / len(logs))


UNIFORM = CodonUsage.from_frequencies(
    "uniform",
    {codon: 1.0 for codon in CODON_TO_AA},
    description="Every synonymous codon equally likely; matches the original scripts.",
)
HUMAN = CodonUsage.from_frequencies(
    "human",
    _HUMAN,
    description="Approximate Homo sapiens relative synonymous codon usage.",
)
ECOLI = CodonUsage.from_frequencies(
    "ecoli",
    _ECOLI,
    description="Approximate Escherichia coli K-12 relative synonymous codon usage.",
)

TABLES: dict[str, CodonUsage] = {t.name: t for t in (UNIFORM, HUMAN, ECOLI)}


def available_tables() -> list[str]:
    """Return the names of the bundled codon usage tables."""
    return sorted(TABLES)


def get_table(name: str) -> CodonUsage:
    """Return a bundled table by name.

    Raises:
        KeyError: If ``name`` is not a bundled table.
    """
    key = name.strip().lower()
    if key not in TABLES:
        raise KeyError(
            f"unknown codon usage table {name!r}; "
            f"available: {', '.join(available_tables())}"
        )
    return TABLES[key]


def load_table(path: str | Path, *, name: str | None = None) -> CodonUsage:
    """Load a codon usage table from a JSON, CSV or TSV file.

    JSON files hold an object mapping codon to frequency or count. CSV and TSV
    files hold two columns, codon and value, with an optional header row.

    Raises:
        ValueError: If the file extension is not ``.json``, ``.csv`` or
            ``.tsv``, or a row cannot be parsed.
    """
    path = Path(path)
    return table_from_text(
        path.read_text(),
        path.suffix,
        name=name or path.stem,
        source=str(path),
    )


def table_from_text(
    text: str,
    fmt: str,
    *,
    name: str,
    source: str = "table",
) -> CodonUsage:
    """Build a codon usage table from the text of a JSON, CSV or TSV table.

    Args:
        text: The table contents.
        fmt: ``"json"``, ``"csv"`` or ``"tsv"``, with or without a leading dot.
        name: Name to give the table.
        source: Where the text came from, used to prefix error messages.

    Raises:
        ValueError: If ``fmt`` is unsupported, a row cannot be parsed, or no
            codon rows are found.
    """
    fmt = fmt.lower().lstrip(".")
    if fmt == "json":
        payload = json.loads(text)
        if not isinstance(payload, dict):
            raise ValueError(f"{source}: expected a JSON object of codon -> value")
        try:
            frequencies = {str(k): float(v) for k, v in payload.items()}
        except (TypeError, ValueError):
            raise ValueError(f"{source}: every JSON value must be a number") from None
    elif fmt in {"csv", "tsv"}:
        delimiter = "\t" if fmt == "tsv" else ","
        frequencies = {}
        for lineno, row in enumerate(
            csv.reader(io.StringIO(text, newline=""), delimiter=delimiter), 1
        ):
            if not row or row[0].startswith("#"):
                continue
            if len(row) < 2:
                raise ValueError(f"{source}:{lineno}: expected two columns")
            codon = row[0].strip().upper().replace("U", "T")
            if codon not in CODON_TO_AA:
                continue  # header or annotation row
            try:
                frequencies[codon] = float(row[1])
            except ValueError:
                raise ValueError(
                    f"{source}:{lineno}: {row[1]!r} is not a number"
                ) from None
    else:
        raise ValueError(f"{source}: expected a .json, .csv or .tsv table")

    if not frequencies:
        raise ValueError(f"{source}: no codon rows found")
    return CodonUsage.from_frequencies(
        name,
        frequencies,
        description=f"Loaded from {source}",
    )
