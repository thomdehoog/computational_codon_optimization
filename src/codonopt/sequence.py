"""Core nucleotide/protein primitives.

Everything here is pure, side-effect free and cheap: the optimizer calls these
functions in its inner loop, so lookup tables are built once at import time
rather than per call.
"""

from __future__ import annotations

from typing import Iterator

__all__ = [
    "AA_TO_CODONS",
    "CODON_TO_AA",
    "AMINO_ACIDS",
    "STOP",
    "translate",
    "reverse_complement",
    "iter_codons",
    "split_codons",
    "gc_fraction",
    "gc_count",
    "gc_prefix",
    "window_gc_count",
    "validate_dna",
    "validate_protein",
]

STOP = "*"

#: The standard genetic code, grouped by the amino acid each codon encodes.
AA_TO_CODONS: dict[str, tuple[str, ...]] = {
    "A": ("GCT", "GCC", "GCA", "GCG"),
    "C": ("TGT", "TGC"),
    "D": ("GAT", "GAC"),
    "E": ("GAA", "GAG"),
    "F": ("TTT", "TTC"),
    "G": ("GGT", "GGC", "GGA", "GGG"),
    "H": ("CAT", "CAC"),
    "I": ("ATT", "ATC", "ATA"),
    "K": ("AAA", "AAG"),
    "L": ("TTA", "TTG", "CTT", "CTC", "CTA", "CTG"),
    "M": ("ATG",),
    "N": ("AAT", "AAC"),
    "P": ("CCT", "CCC", "CCA", "CCG"),
    "Q": ("CAA", "CAG"),
    "R": ("CGT", "CGC", "CGA", "CGG", "AGA", "AGG"),
    "S": ("TCT", "TCC", "TCA", "TCG", "AGT", "AGC"),
    "T": ("ACT", "ACC", "ACA", "ACG"),
    "V": ("GTT", "GTC", "GTA", "GTG"),
    "W": ("TGG",),
    "Y": ("TAT", "TAC"),
    STOP: ("TAA", "TAG", "TGA"),
}

#: Reverse mapping, codon -> amino acid.
CODON_TO_AA: dict[str, str] = {
    codon: aa for aa, codons in AA_TO_CODONS.items() for codon in codons
}

#: Amino acid one-letter codes, stop excluded.
AMINO_ACIDS: tuple[str, ...] = tuple(sorted(aa for aa in AA_TO_CODONS if aa != STOP))

_COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
_DNA_ALPHABET = frozenset("ACGT")


class SequenceError(ValueError):
    """Raised when an input sequence is not valid DNA or protein."""


def validate_dna(sequence: str, *, in_frame: bool = False) -> str:
    """Return ``sequence`` upper-cased, raising if it is not unambiguous DNA.

    Args:
        sequence: Candidate nucleotide sequence.
        in_frame: Also require the length to be a multiple of three.

    Raises:
        SequenceError: On an empty sequence, a non-ACGT character, or a length
            that is not a multiple of three when ``in_frame`` is set.
    """
    seq = sequence.strip().upper()
    if not seq:
        raise SequenceError("empty nucleotide sequence")
    bad = set(seq) - _DNA_ALPHABET
    if bad:
        raise SequenceError(
            "sequence contains non-ACGT characters: " + ", ".join(sorted(bad))
        )
    if in_frame and len(seq) % 3:
        raise SequenceError(
            f"sequence length {len(seq)} is not a multiple of three; "
            "it cannot be read as codons"
        )
    return seq


def validate_protein(sequence: str, *, allow_internal_stop: bool = False) -> str:
    """Return ``sequence`` upper-cased, raising if it is not a valid protein.

    A single trailing ``*`` is allowed and preserved; internal stops are
    rejected unless ``allow_internal_stop`` is set.
    """
    seq = sequence.strip().upper()
    if not seq:
        raise SequenceError("empty protein sequence")
    bad = set(seq) - set(AA_TO_CODONS)
    if bad:
        raise SequenceError(
            "protein contains unsupported residues: " + ", ".join(sorted(bad))
        )
    if not allow_internal_stop and STOP in seq[:-1]:
        raise SequenceError(
            f"internal stop codon at protein position {seq.index(STOP) + 1}"
        )
    return seq


def iter_codons(sequence: str) -> Iterator[str]:
    """Yield successive in-frame codons of ``sequence``."""
    for i in range(0, len(sequence) - len(sequence) % 3, 3):
        yield sequence[i : i + 3]


def split_codons(sequence: str) -> list[str]:
    """Return ``sequence`` as a list of in-frame codons.

    Raises:
        SequenceError: If the length is not a multiple of three.
    """
    if len(sequence) % 3:
        raise SequenceError(
            f"sequence length {len(sequence)} is not a multiple of three"
        )
    return [sequence[i : i + 3] for i in range(0, len(sequence), 3)]


def translate(sequence: str) -> str:
    """Translate an in-frame nucleotide sequence into one-letter amino acids.

    Raises:
        SequenceError: If the length is not a multiple of three or a codon is
            not part of the standard genetic code.
    """
    codons = split_codons(sequence.upper())
    try:
        return "".join([CODON_TO_AA[codon] for codon in codons])
    except KeyError as exc:  # pragma: no cover - exercised via tests
        raise SequenceError(f"not a valid codon: {exc.args[0]!r}") from None


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of ``sequence``."""
    return sequence.translate(_COMPLEMENT)[::-1]


def gc_count(sequence: str) -> int:
    """Return the number of G and C bases in ``sequence``."""
    return sequence.count("G") + sequence.count("C")


def gc_fraction(sequence: str) -> float:
    """Return the GC content of ``sequence`` as a fraction between 0 and 1."""
    if not sequence:
        return 0.0
    return gc_count(sequence) / len(sequence)


def gc_prefix(sequence: str) -> list[int]:
    """Return cumulative GC counts, so window GC is an O(1) lookup.

    ``prefix[i]`` is the number of G/C bases in ``sequence[:i]``; the returned
    list has ``len(sequence) + 1`` entries.
    """
    prefix = [0] * (len(sequence) + 1)
    total = 0
    for i, base in enumerate(sequence):
        if base == "G" or base == "C":
            total += 1
        prefix[i + 1] = total
    return prefix


def window_gc_count(prefix: list[int], start: int, end: int) -> int:
    """Return the GC count of ``sequence[start:end]`` from a ``gc_prefix`` list."""
    return prefix[end] - prefix[start]
