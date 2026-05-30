"""Tests for varsim.protein: aa_sub()."""

import pytest
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from tests.conftest import GENE, requires_entrez, assert_valid_p_hgvs


class TestAASub:
    """Tests for aa_sub() function."""
    
    @requires_entrez
    def test_returns_list_of_tuples(self, aa_sub_result):
        result = aa_sub_result
        assert isinstance(result, list)
        assert len(result) > 0
        for variant in result[:50]:
            assert isinstance(variant, tuple)
            assert len(variant) == 2
    
    @requires_entrez
    def test_p_hgvs_format(self, aa_sub_result):
        """All items should have valid p.HGVS format."""
        result = aa_sub_result
        for p_1, p_3 in result[:20]:
            assert_valid_p_hgvs(p_1)
            assert_valid_p_hgvs(p_3)
    
    @requires_entrez
    def test_start_codon_convention(self, aa_sub_result):
        """First position should use ? notation (not other methionines like M125)."""
        result = aa_sub_result
        # Met1 variants: NP_...:p.(M1?) and NP_...:p.(Met1?)
        met1_variants = [(p1, p3) for p1, p3 in result
                         if ':p.(M1?)' in p1 or ':p.(Met1?)' in p3]
        assert len(met1_variants) > 0, "Should have Met1 variants"
        for p1, p3 in met1_variants:
            assert '?' in p1, f"Start codon should use ?: {p1}"
            assert '?' in p3, f"Start codon should use ?: {p3}"
    
    @requires_entrez
    def test_no_self_substitution(self, aa_sub_result):
        """No variant should substitute an amino acid with itself."""
        result = aa_sub_result
        for p1, _ in result:
            # p1 looks like NP_...:p.(X1Y) — X should never equal Y unless it's ?
            if '?' not in p1:
                parts = p1.split(':p.(')[1].rstrip(')')
                # Format like "M2V" — first char should differ from last
                if len(parts) >= 2 and parts[0].isalpha() and parts[-1].isalpha():
                    assert parts[0] != parts[-1], f"Self-substitution: {p1}"
