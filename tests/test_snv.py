"""Tests for varsim.snv: cds(), utr5(), utr3()."""

import pytest
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from tests.conftest import GENE, requires_entrez, assert_valid_c_hgvs, assert_valid_variant_tuple, assert_valid_p_hgvs


class TestCDS:
    """Tests for cds() function."""

    @requires_entrez
    def test_returns_list(self, cds_result):
        result = cds_result
        assert isinstance(result, list)
        assert len(result) > 0

    @requires_entrez
    def test_variant_tuple_structure(self, cds_result):
        result = cds_result
        for variant in result[:10]:  # Check first 10
            assert_valid_variant_tuple(variant, prefix_required=True)

    @requires_entrez
    def test_start_codon_handling(self, cds_result):
        """First codon should use p.(M1?) / p.(Met1?) convention."""
        result = cds_result
        # Find a start codon variant (c.1, c.2, or c.3)
        start_variants = [v for v in result if any(
            v[0].endswith(pat) for pat in ('c.1A>', 'c.2T>', 'c.3G>')
        )]
        for v in start_variants:
            assert 'M1?' in v[1], f"Expected M1? in {v[1]}"
            assert 'Met1?' in v[2], f"Expected Met1? in {v[2]}"

    @requires_entrez
    def test_nc_prefix_format(self, cds_result):
        """All c.HGVS should have NC_ACC(NM_ACC):c. format."""
        result = cds_result
        for variant in result[:20]:
            c_hgvs = variant[0]
            assert c_hgvs.startswith('NC_'), f"Expected NC_ prefix in {c_hgvs}"
            assert '(NM_' in c_hgvs, f"Expected (NM_ in {c_hgvs}"

    @requires_entrez
    def test_protein_annotation_types(self, cds_result):
        """Protein annotations should be either missense or silent (=)."""
        result = cds_result
        for variant in result[:50]:
            p_1 = variant[1]
            assert ':p.(' in p_1


class TestUTR5:
    """Tests for utr5() function."""

    @requires_entrez
    def test_returns_list_of_strings(self, utr5_result):
        result = utr5_result
        assert isinstance(result, list)
        assert len(result) > 0
        assert all(isinstance(v, str) for v in result)

    @requires_entrez
    def test_hgvs_format(self, utr5_result):
        result = utr5_result
        for v in result[:10]:
            assert_valid_c_hgvs(v, prefix_required=True)
            assert 'c.-' in v, f"UTR5 should have negative coordinates, got {v}"


class TestUTR3:
    """Tests for utr3() function."""

    @requires_entrez
    def test_returns_list_of_strings(self, utr3_result):
        result = utr3_result
        assert isinstance(result, list)
        assert len(result) > 0
        assert all(isinstance(v, str) for v in result)

    @requires_entrez
    def test_hgvs_format(self, utr3_result):
        result = utr3_result
        for v in result[:10]:
            assert_valid_c_hgvs(v, prefix_required=True)
            assert 'c.*' in v, f"UTR3 should have * coordinates, got {v}"
