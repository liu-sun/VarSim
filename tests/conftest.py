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
# Common validation helpers
# ============================================================================

def assert_valid_p_hgvs(hgvs_str: str):
    """Assert that a p.HGVS string has the expected format.

    Expected format: NP_ACC:p.(XxxN...) or NP_ACC:p.(XN...)
    """
    assert isinstance(hgvs_str, str), f"Expected str, got {type(hgvs_str)}"
    assert ":p.(" in hgvs_str, f"Missing ':p.(' in {hgvs_str}"
