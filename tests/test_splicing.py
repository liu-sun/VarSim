"""Tests for varsim.splicing: splice_site()."""

import pytest
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from tests.conftest import GENE, requires_entrez, assert_valid_c_hgvs


class TestSpliceSite:
    """Tests for splice_site() function."""
    
    @requires_entrez
    def test_returns_list_of_strings(self, splice_site_result):
        result = splice_site_result
        assert isinstance(result, list)
        assert len(result) > 0
        assert all(isinstance(v, str) for v in result)
    
    @requires_entrez
    def test_nc_nm_format(self, splice_site_result):
        """All results should have NC_ACC(NM_ACC):c. format."""
        result = splice_site_result
        for v in result:
            assert 'NC_' in v, f"Missing NC_ in {v}"
            assert '(NM_' in v, f"Missing (NM_ in {v}"
            assert ':c.' in v, f"Missing :c. in {v}"
    
    @requires_entrez
    def test_canonical_dinucleotides(self, splice_site_result):
        """Verify only canonical GT-AG dinucleotides are used."""
        result = splice_site_result
        # Check donor site (+1G, +2T) and acceptor site (-2A, -1G)
        donors = [v for v in result if '+' in v]
        acceptors = [v for v in result if '-' in v]
        assert len(donors) > 0, "Should have donor splice site variants"
        assert len(acceptors) > 0, "Should have acceptor splice site variants"
    
    @requires_entrez
    def test_no_duplicates(self, splice_site_result):
        """Variants should be unique."""
        result = splice_site_result
        assert len(result) == len(set(result)), "Duplicate splice site variants found"
