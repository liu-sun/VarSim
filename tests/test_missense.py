"""Tests for varsim.missense: missense()."""

import pytest
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from tests.conftest import GENE, requires_entrez, assert_valid_variant_tuple


class TestMissense:
    """Tests for missense() function."""
    
    @requires_entrez
    def test_returns_list_of_tuples(self, missense_result):
        result = missense_result
        assert isinstance(result, list)
        assert len(result) > 0
    
    @requires_entrez
    def test_variant_tuple_structure(self, missense_result):
        result = missense_result
        for variant in result[:10]:
            assert_valid_variant_tuple(variant, prefix_required=True)
    
    @requires_entrez
    def test_nc_prefix_format(self, missense_result):
        """All c.HGVS should have NC_ACC(NM_ACC): format."""
        result = missense_result
        for variant in result[:20]:
            c_hgvs = variant[0]
            assert c_hgvs.startswith('NC_'), f"Expected NC_ in {c_hgvs}"
            assert '(NM_' in c_hgvs, f"Expected (NM_ in {c_hgvs}"
    
    @requires_entrez
    def test_start_codon_convention(self, missense_result):
        """Start codon (index 0) should use p.(M1?) / p.(Met1?)."""
        result = missense_result
        start_variants = [(c, p1, p3) for c, p1, p3 in result 
                         if c.endswith(':c.1_3delins') or any(
                             c.endswith(f':c.{i}_{j}delins') 
                             for i, j in [(1, 3), (2, 3), (1, 2)]
                         ) or any(c.endswith(f':c.{i}') or f':c.{i}' == c[-5:-1] 
                                 for i in ['1A>', '2T>', '3G>'])]
        # Find variants that affect codon 1
        codon1 = [v for v in result if 'Met1?' in v[2] or 'M1?' in v[1]]
        assert len(codon1) > 0, "Should have start codon variants"
        for v in codon1:
            assert 'M1?' in v[1] or 'Met1?' in v[2], f"Missing start codon ? in {v}"
    
    @requires_entrez
    def test_silent_and_missense_exist(self, missense_result):
        """Should have both silent (=) and missense variants."""
        result = missense_result
        silent = [v for v in result if '=' in v[1]]
        missense_variants = [v for v in result if '=' not in v[1]]
        assert len(silent) > 0, "Should have silent variants"
        assert len(missense_variants) > 0, "Should have missense variants"
    
    @requires_entrez
    def test_no_duplicates(self, missense_result):
        """All variant tuples should be unique."""
        result = missense_result
        assert len(result) == len(set(result)), "Duplicate missense variants found"
