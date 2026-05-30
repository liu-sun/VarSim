"""Tests for varsim.frameshift: frameshift()."""

import pytest
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from tests.conftest import GENE, requires_entrez, assert_valid_variant_tuple


class TestFrameshift:
    """Tests for frameshift() function."""
    
    @requires_entrez
    def test_returns_list_of_tuples(self, frameshift_result):
        result = frameshift_result
        assert isinstance(result, list)
        assert len(result) > 0
    
    @requires_entrez
    def test_variant_tuple_structure(self, frameshift_result):
        result = frameshift_result
        for variant in result[:10]:
            assert_valid_variant_tuple(variant, prefix_required=False)
    
    @requires_entrez
    def test_deletion_variants(self, frameshift_result, cds_length):
        """Should have deletion variants with 'del' in c.HGVS."""
        result = frameshift_result
        deletions = [v for v in result if 'del' in v[0] and 'ins' not in v[0]]
        assert len(deletions) > 0, "Should have deletion variants"
        assert len(deletions) == cds_length, f"Expected {cds_length} deletions, got {len(deletions)}"
    
    @requires_entrez
    def test_insertion_variants(self, frameshift_result, cds_length):
        """Should have insertion variants with 'ins' in c.HGVS."""
        result = frameshift_result
        insertions = [v for v in result if 'ins' in v[0]]
        assert len(insertions) > 0, "Should have insertion variants"
        assert len(insertions) == 4 * cds_length, f"Expected {4 * cds_length} insertions, got {len(insertions)}"
    
    @requires_entrez
    def test_frameshift_protein_notation(self, frameshift_result):
        """Protein effect should use fs* notation for frameshift."""
        result = frameshift_result
        fs_variants = [v for v in result if 'fs' in v[1]]
        assert len(fs_variants) > 0, "Should have fs variants"
        for c_hgvs, p_1, p_3 in fs_variants[:10]:
            assert 'fs*' in p_1, f"Expected fs* in {p_1}"
            assert 'fs*' in p_3, f"Expected fs* in {p_3}"
    
    @requires_entrez
    def test_start_codon_convention(self, frameshift_result):
        """Start codon variants should use p.(M1?) / p.(Met1?)."""
        result = frameshift_result
        start_variants = [v for v in result if 'Met1?' in v[2] or 'M1?' in v[1]]
        assert len(start_variants) > 0, "Should have start codon variants"
        for v in start_variants:
            assert 'M1?' in v[1] or 'Met1?' in v[2], f"Missing start codon ? in {v}"
    
    @requires_entrez
    def test_no_duplicates(self, frameshift_result):
        """All variant tuples should be unique."""
        result = frameshift_result
        assert len(result) == len(set(result)), "Duplicate frameshift variants found"
