"""Shared test fixtures and utilities for VarSim.

Uses G6PD (glucose-6-phosphate dehydrogenase) as the canonical test gene.
All Entrez API calls are cached at session scope to avoid NCBI rate limiting.
"""

import os
import pytest

# --- Test gene ---
GENE = "G6PD"

# Expected accessions for G6PD MANE Select (v2 transcript)
EXPECTED_NM = "NM_001360016.2"
EXPECTED_NP = "NP_001347945.1"  # Wait - protein may be NP_001346945.1, verified dynamically
EXPECTED_NC = "NC_000023.11"


# --- Environment check ---
def has_entrez_credentials():
    """Check if NCBI Entrez EMAIL and API_KEY environment variables are set."""
    return bool(os.environ.get("EMAIL") and os.environ.get("API_KEY"))


requires_entrez = pytest.mark.skipif(
    not has_entrez_credentials(),
    reason="EMAIL and API_KEY environment variables must be set for NCBI Entrez access",
)


# ============================================================================
# Session-scoped cached fixtures — each API endpoint is called exactly once
# ============================================================================

@pytest.fixture(scope="session", autouse=True)
def _clear_entrez_cache():
    """Clear any stale Entrez caches before the test session."""
    from varsim._fetch import nm, np, nc
    nm.cache_clear()
    np.cache_clear()
    nc.cache_clear()


@pytest.fixture(scope="session")
def cached_nm_seqrecord():
    """G6PD MANE Select nucleotide SeqRecord (cached once per test session)."""
    import time
    time.sleep(1)  # Rate-limit courtesy delay
    from varsim._fetch import nm
    return nm(GENE)


@pytest.fixture(scope="session")
def cached_np_seqrecord():
    """G6PD MANE Select protein SeqRecord (cached once per test session)."""
    import time
    time.sleep(1)  # Rate-limit courtesy delay
    from varsim._fetch import np
    return np(GENE)


@pytest.fixture(scope="session")
def cached_nc_accession():
    """G6PD genomic NC_ accession string (cached once per test session)."""
    import time
    time.sleep(1)  # Rate-limit courtesy delay
    from varsim._fetch import nc
    return nc(GENE)


# ============================================================================
# Cached public function results — each function computed once
# ============================================================================

@pytest.fixture(scope="session")
def cds_result():
    """Cached cds() result for G6PD."""
    from varsim.snv import cds
    return cds(GENE)


@pytest.fixture(scope="session")
def utr5_result():
    """Cached utr5() result for G6PD."""
    from varsim.snv import utr5
    return utr5(GENE)


@pytest.fixture(scope="session")
def utr3_result():
    """Cached utr3() result for G6PD."""
    from varsim.snv import utr3
    return utr3(GENE)


@pytest.fixture(scope="session")
def splice_site_result():
    """Cached splice_site() result for G6PD."""
    from varsim.splicing import splice_site
    return splice_site(GENE)


@pytest.fixture(scope="session")
def aa_sub_result():
    """Cached aa_sub() result for G6PD."""
    from varsim.protein import aa_sub
    return aa_sub(GENE)


@pytest.fixture(scope="session")
def codon_sub_result():
    """Cached codon_sub() result for G6PD."""
    from varsim.codon import codon_sub
    return codon_sub(GENE)


@pytest.fixture(scope="session")
def missense_result():
    """Cached missense() result for G6PD."""
    from varsim.missense import missense
    return missense(GENE)


@pytest.fixture(scope="session")
def frameshift_result():
    """Cached frameshift() result for G6PD."""
    from varsim.frameshift import frameshift
    return frameshift(GENE)


@pytest.fixture(scope="session")
def cds_length(cached_nm_seqrecord):
    """G6PD CDS length computed from the cached nucleotide record."""
    for feature in cached_nm_seqrecord.features:
        if feature.type == "CDS":
            return len(feature.extract(cached_nm_seqrecord).seq)
    return 0


# ============================================================================
# Common validation helpers
# ============================================================================

def assert_valid_c_hgvs(hgvs_str: str, prefix_required: bool = True):
    """Assert that a c.HGVS string has the expected format.

    Expected format: NC_ACC(NM_ACC):c.COORDINATE or NM_ACC:c.COORDINATE
    """
    assert isinstance(hgvs_str, str), f"Expected str, got {type(hgvs_str)}"
    assert ":c." in hgvs_str, f"Missing ':c.' in {hgvs_str}"
    if prefix_required:
        # Should have NC_...(NM_...): format
        assert "(" in hgvs_str and ")" in hgvs_str, f"Missing NC_ prefix in {hgvs_str}"


def assert_valid_p_hgvs(hgvs_str: str):
    """Assert that a p.HGVS string has the expected format.

    Expected format: NP_ACC:p.(XxxN...) or NP_ACC:p.(XN...)
    """
    assert isinstance(hgvs_str, str), f"Expected str, got {type(hgvs_str)}"
    assert ":p.(" in hgvs_str, f"Missing ':p.(' in {hgvs_str}"


def assert_valid_variant_tuple(variant, prefix_required: bool = True):
    """Assert a CDS/missense/frameshift variant is a valid 3-tuple.

    Returns (NC_ACC(NM_ACC):c.HGVS, NP_ACC:p.HGVS_1letter, NP_ACC:p.HGVS_3letter)
    """
    assert isinstance(variant, tuple), f"Expected tuple, got {type(variant)}"
    assert len(variant) == 3, f"Expected 3-tuple, got {len(variant)}-tuple"
    assert_valid_c_hgvs(variant[0], prefix_required=prefix_required)
    assert_valid_p_hgvs(variant[1])
    assert_valid_p_hgvs(variant[2])
