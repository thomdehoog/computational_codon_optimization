"""Tests for splice site scoring."""

from __future__ import annotations

import pytest

from codonopt.sequence import reverse_complement
from codonopt.splice import (
    ACCEPTOR_WINDOW,
    DONOR_WINDOW,
    MaxEntScorer,
    MotifScorer,
    find_splice_sites,
    get_scorer,
    maxent_available,
)

CANONICAL_DONOR = "CAGGTAAGT"
CANONICAL_ACCEPTOR = "TTTTTTTTTTTTTTTTTCAGGTC"


def test_window_sizes_match_the_maxentscan_convention():
    assert DONOR_WINDOW == len(CANONICAL_DONOR) == 9
    assert ACCEPTOR_WINDOW == len(CANONICAL_ACCEPTOR) == 23


def test_motif_scorer_ranks_the_canonical_sites_highest():
    scorer = MotifScorer()
    assert scorer.score_donor(CANONICAL_DONOR) == 6.0
    assert scorer.score_donor("CAGGTTTTT") < scorer.score_donor(CANONICAL_DONOR)
    # Without GT at the exon/intron boundary there is no donor at all.
    assert scorer.score_donor("CAGATAAGT") == -10.0
    # Eighteen pyrimidines in the tract, with AG at the intron's 3' end.
    assert scorer.score_acceptor(CANONICAL_ACCEPTOR) == 18.0
    assert scorer.score_acceptor("AAAAAAAAAAAAAAAAAAAAGTC") == -10.0


def test_motif_scorer_cutoffs_are_on_its_own_scale():
    scorer = MotifScorer()
    assert scorer.score_donor(CANONICAL_DONOR) > scorer.default_donor_cutoff
    assert scorer.score_acceptor(CANONICAL_ACCEPTOR) > scorer.default_acceptor_cutoff


def test_find_splice_sites_reports_forward_positions():
    scorer = MotifScorer()
    sequence = "AAACCC" + CANONICAL_DONOR + "AAACCC"
    sites = find_splice_sites(sequence, scorer, strand="forward")
    assert [site.position for site in sites if site.site_type == "donor"] == [6]
    assert all(site.strand == "forward" for site in sites)


def test_reverse_strand_sites_are_found_at_the_mirrored_window():
    scorer = MotifScorer()
    sequence = "AAACCC" + reverse_complement(CANONICAL_DONOR) + "AAACCC"
    forward = find_splice_sites(sequence, scorer, strand="forward")
    reverse = find_splice_sites(
        sequence, scorer, strand="reverse", donor_cutoff=5.0, acceptor_cutoff=99
    )
    assert not [site for site in forward if site.site_type == "donor"]
    assert [site.position for site in reverse if site.site_type == "donor"] == [6]


def test_a_cutoff_is_the_highest_acceptable_score():
    """Reporting uses the same comparison as the penalty: strictly above.

    With a backend whose scores are integers, reporting sites at exactly the
    cutoff while charging ``score - cutoff`` would report sites the optimizer
    has no reason to remove.
    """
    scorer = MotifScorer()
    sequence = "AAACCC" + CANONICAL_DONOR + "AAACCC"
    assert find_splice_sites(sequence, scorer, donor_cutoff=5.0, acceptor_cutoff=99)
    assert not find_splice_sites(
        sequence, scorer, donor_cutoff=6.0, acceptor_cutoff=99
    )


def test_strand_argument_is_validated():
    with pytest.raises(ValueError, match="unknown strand"):
        find_splice_sites("ACGT" * 10, MotifScorer(), strand="both-ish")


def test_get_scorer_selects_backends():
    assert get_scorer("none") is None
    assert isinstance(get_scorer("motif"), MotifScorer)
    assert get_scorer("auto") is not None
    with pytest.raises(ValueError, match="unknown splice backend"):
        get_scorer("bogus")


def test_auto_prefers_maxent_when_available():
    scorer = get_scorer("auto")
    expected = "maxent" if maxent_available() else "motif"
    assert scorer.name == expected


@pytest.mark.skipif(not maxent_available(), reason="maxentpy is not installed")
def test_maxent_scores_the_canonical_sites_highly():
    scorer = MaxEntScorer()
    assert scorer.score_donor(CANONICAL_DONOR) > 8
    assert scorer.score_donor("AAAAAAAAA") < 0
    assert scorer.score_acceptor(CANONICAL_ACCEPTOR) > 8
    assert scorer.score_acceptor("AAAAAAAAAAAAAAAAAAAAGTC") < 0


@pytest.mark.skipif(not maxent_available(), reason="maxentpy is not installed")
def test_maxent_matrices_are_loaded_once():
    """Two scorers must share the cached matrices.

    The original implementation reloaded both matrices on every call, which
    cost about 50 ms per call; the whole point of caching them is that this
    stops being the dominant cost of the run.
    """
    first = MaxEntScorer()
    second = MaxEntScorer()
    assert first._matrix5 is second._matrix5
    assert first._matrix3 is second._matrix3


@pytest.mark.skipif(maxent_available(), reason="maxentpy is installed")
def test_missing_maxent_gives_an_actionable_error():
    with pytest.raises(RuntimeError, match="maxentpy"):
        get_scorer("maxent")


def test_reported_sites_are_exactly_the_penalised_ones():
    """Whatever the report shows, the objective must be paid to remove."""
    from codonopt import DesignConfig, Weights, get_table
    from codonopt.state import SequenceState

    scorer = MotifScorer()
    only_splice = Weights(
        gc_global=0, gc_window=0, repeat=0, hairpin=0, homopolymer=0,
        motif=0, splice=1.0, usage=0,
    )
    config = DesignConfig(
        table=get_table("uniform"), splice_scorer=scorer, weights=only_splice
    )
    for sequence in (
        "CAGGTAAGTATG",            # a strong donor
        "CAGGTAAACATG",            # a donor at exactly the default cutoff
        "ATGAAACCCTTTGGGAAACCCTTT",
        CANONICAL_ACCEPTOR + "ATGG",  # padded to a whole number of codons
    ):
        penalty = SequenceState(sequence, config).breakdown().raw["splice"]
        reported = find_splice_sites(
            sequence,
            scorer,
            donor_cutoff=config.effective_donor_cutoff,
            acceptor_cutoff=config.effective_acceptor_cutoff,
        )
        assert bool(reported) == (penalty > 0), sequence
