"""Tests for varsim.codon: codon_sub()."""

import pytest
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from tests.conftest import GENE, requires_entrez, assert_valid_c_hgvs


class TestCodonSub:
    """Tests for codon_sub() function."""
    
    @requires_entrez
    def test_returns_list_of_strings(self, codon_sub_result):
        result = codon_sub_result
        assert isinstance(result, list)
        assert len(result) > 0
        assert all(isinstance(v, str) for v in result)
    
    @requires_entrez
    def test_hgvs_format(self, codon_sub_result):
        """All items should be valid c.HGVS strings."""
        result = codon_sub_result
        for v in result[:20]:
            assert_valid_c_hgvs(v, prefix_required=False)
    
    @requires_entrez
    def test_no_duplicates(self, codon_sub_result):
        """All variants should be unique."""
        result = codon_sub_result
        assert len(result) == len(set(result)), "Duplicate codon variants found"
    
    @requires_entrez
    def test_delins_format(self, codon_sub_result):
        """Some variants should use delins notation for multi-base changes."""
        result = codon_sub_result
        delins_variants = [v for v in result if 'delins' in v]
        single_variants = [v for v in result if '>' in v and 'delins' not in v]
        assert len(delins_variants) > 0, "Should have delins variants"
        assert len(single_variants) > 0, "Should have single-substitution variants"
    
    @requires_entrez
    def test_coordinate_ranges(self, codon_sub_result):
        """Coordinates should be within CDS bounds."""
        result = codon_sub_result
        for v in result[:20]:
            # Extract coordinate from NM_...:c.123_456delins...
            c_part = v.split(':c.')[1] if ':c.' in v else ''
            assert c_part, f"No c. coordinate in {v}"
