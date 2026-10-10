"""Shared fixtures and helpers."""

from __future__ import annotations

import random

import pytest

from codonopt import DesignConfig, Weights, get_table, resolve_motifs, seed_sequence
from codonopt.splice import MotifScorer

RESIDUES = "ACDEFGHIKLMNPQRSTVWY"


def random_protein(length: int, seed: int = 0, *, stop: bool = False) -> str:
    """Return a reproducible pseudo-random protein of ``length`` residues."""
    rng = random.Random(seed)
    protein = "M" + "".join(rng.choice(RESIDUES) for _ in range(length - 1))
    return protein + "*" if stop else protein


def duplicated_protein(block: int = 60, seed: int = 0) -> str:
    """Return a protein made of one block repeated twice."""
    rng = random.Random(seed)
    half = "".join(rng.choice(RESIDUES) for _ in range(block))
    return "M" + half + half


def duplicated_pair(block: int = 60, seed: int = 0) -> tuple[str, str]:
    """Return a protein and a coding sequence that repeats itself exactly.

    A protein with a duplicated block does not by itself produce repeats in
    the DNA -- independently chosen codons differ between the two copies --
    so the coding sequence is built by repeating one block's codons, which
    guarantees the repeat term has something to work on.
    """
    rng = random.Random(seed)
    half_protein = "".join(rng.choice(RESIDUES) for _ in range(block))
    half_dna = random_dna(half_protein, seed=seed + 1)
    return "M" + half_protein + half_protein, "ATG" + half_dna + half_dna


def random_dna(protein: str, seed: int = 0) -> str:
    """Return a random coding sequence for ``protein``."""
    return seed_sequence(protein, get_table("uniform"), random.Random(seed), "random")


@pytest.fixture
def config() -> DesignConfig:
    """A configuration exercising every penalty term without a hard dependency."""
    return DesignConfig(
        table=get_table("human"),
        gc_window=18,
        gc_step=3,
        repeat_k=9,
        hairpin_stem=8,
        hairpin_min_loop=3,
        hairpin_max_loop=24,
        max_homopolymer=4,
        motifs=tuple(resolve_motifs(["EcoRI", "BsaI", "AATAAA"])),
        splice_scorer=MotifScorer(),
        weights=Weights(),
    )


@pytest.fixture
def plain_config() -> DesignConfig:
    """A configuration with the costly terms switched off, for fast tests."""
    return DesignConfig(
        table=get_table("human"),
        splice_scorer=None,
        motifs=(),
    )
