"""Tests for varsim.transcription: genome-transcript coordinate mapping."""

import pytest
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from varsim.transcription import (
    _cdna_to_genomic,
    _genomic_to_cdna,
    _get_true_cdna_and_cds_length,
    _hgvs_cdna_coord,
    _format_coord,
    _format_genomic_coord,
    _format_hgvs,
    _get_variant_suffix,
    get_cds_exon_map,
    c_to_g,
    g_to_c,
)
from varsim.parser import parse, HGVSTag

# ---------------------------------------------------------------------------
# Mock exon maps for testing
# ---------------------------------------------------------------------------

EXON_MAP_PLUS = [
    {"exon": 1, "cds_start": 0, "cds_end": 100, "genomic_start": 1000, "genomic_end": 1099, "strand": 1},
    {"exon": 2, "cds_start": 100, "cds_end": 200, "genomic_start": 2000, "genomic_end": 2099, "strand": 1},
]

EXON_MAP_MINUS = [
    {"exon": 1, "cds_start": 0, "cds_end": 100, "genomic_start": 5000, "genomic_end": 5099, "strand": -1},
    {"exon": 2, "cds_start": 100, "cds_end": 200, "genomic_start": 4000, "genomic_end": 4099, "strand": -1},
]

# ---------------------------------------------------------------------------
# _hgvs_cdna_coord
# ---------------------------------------------------------------------------

class TestHgvsCdnaCoord:
    """Tests for _hgvs_cdna_coord helper."""

    def test_coding_position(self):
        assert _hgvs_cdna_coord(1, 1000) == "1"
        assert _hgvs_cdna_coord(500, 1000) == "500"
        assert _hgvs_cdna_coord(1000, 1000) == "1000"

    def test_5utr_position(self):
        assert _hgvs_cdna_coord(-1, 1000) == "-1"
        assert _hgvs_cdna_coord(-5, 1000) == "-5"
        assert _hgvs_cdna_coord(0, 1000) == "0"

    def test_3utr_position(self):
        assert _hgvs_cdna_coord(1001, 1000) == "*1"
        assert _hgvs_cdna_coord(1050, 1000) == "*50"


# ---------------------------------------------------------------------------
# _cdna_to_genomic — forward strand
# ---------------------------------------------------------------------------

class TestCdnaToGenomicPlus:
    """Tests for _cdna_to_genomic on forward strand."""

    def test_first_coding_base(self):
        g, in_exon = _cdna_to_genomic(1, EXON_MAP_PLUS)
        assert g == 1000
        assert in_exon is True

    def test_last_coding_base(self):
        g, in_exon = _cdna_to_genomic(200, EXON_MAP_PLUS)
        assert g == 2099
        assert in_exon is True

    def test_exon_boundary(self):
        g, in_exon = _cdna_to_genomic(100, EXON_MAP_PLUS)
        assert g == 1099
        g, in_exon = _cdna_to_genomic(101, EXON_MAP_PLUS)
        assert g == 2000

    def test_5utr(self):
        g, in_exon = _cdna_to_genomic(-5, EXON_MAP_PLUS)
        assert g == 995
        assert in_exon is False

    def test_3utr(self):
        g, in_exon = _cdna_to_genomic(201, EXON_MAP_PLUS)
        assert g == 2100
        assert in_exon is False


# ---------------------------------------------------------------------------
# _cdna_to_genomic — reverse strand
# ---------------------------------------------------------------------------

class TestCdnaToGenomicMinus:
    """Tests for _cdna_to_genomic on reverse strand."""

    def test_first_coding_base(self):
        g, in_exon = _cdna_to_genomic(1, EXON_MAP_MINUS)
        assert g == 5099
        assert in_exon is True

    def test_last_coding_base(self):
        g, in_exon = _cdna_to_genomic(200, EXON_MAP_MINUS)
        assert g == 4000
        assert in_exon is True

    def test_5utr(self):
        g, in_exon = _cdna_to_genomic(-1, EXON_MAP_MINUS)
        assert g == 5100
        assert in_exon is False

    def test_3utr(self):
        g, in_exon = _cdna_to_genomic(201, EXON_MAP_MINUS)
        assert g == 3999
        assert in_exon is False


# ---------------------------------------------------------------------------
# _genomic_to_cdna — forward strand
# ---------------------------------------------------------------------------

class TestGenomicToCdnaPlus:
    """Tests for _genomic_to_cdna on forward strand."""

    def test_coding_positions(self):
        c, ie = _genomic_to_cdna(1000, EXON_MAP_PLUS)
        assert c == 1
        c, ie = _genomic_to_cdna(2099, EXON_MAP_PLUS)
        assert c == 200

    def test_5utr(self):
        c, ie = _genomic_to_cdna(995, EXON_MAP_PLUS)
        assert c == -5
        assert ie is False

    def test_3utr(self):
        c, ie = _genomic_to_cdna(2100, EXON_MAP_PLUS)
        assert c == 201
        assert ie is False


# ---------------------------------------------------------------------------
# _genomic_to_cdna — reverse strand
# ---------------------------------------------------------------------------

class TestGenomicToCdnaMinus:
    """Tests for _genomic_to_cdna on reverse strand."""

    def test_coding_positions(self):
        c, ie = _genomic_to_cdna(5099, EXON_MAP_MINUS)
        assert c == 1
        c, ie = _genomic_to_cdna(4000, EXON_MAP_MINUS)
        assert c == 200

    def test_5utr(self):
        c, ie = _genomic_to_cdna(5100, EXON_MAP_MINUS)
        assert c == -1
        assert ie is False

    def test_3utr(self):
        c, ie = _genomic_to_cdna(3999, EXON_MAP_MINUS)
        assert c == 201
        assert ie is False


# ---------------------------------------------------------------------------
# Round-trip tests
# ---------------------------------------------------------------------------

class TestRoundTrip:
    """Verify c→g→c and g→c→g round-trip identity."""

    @pytest.mark.parametrize("exon_map", [EXON_MAP_PLUS, EXON_MAP_MINUS])
    @pytest.mark.parametrize("cdna_pos", [-10, -5, -1, 1, 50, 100, 101, 150, 200, 201, 250])
    def test_cdna_round_trip(self, exon_map, cdna_pos):
        """cDNA position maps to genomic and back unchanged."""
        g, _ = _cdna_to_genomic(cdna_pos, exon_map)
        c, _ = _genomic_to_cdna(g, exon_map)
        assert c == cdna_pos, f"c.{cdna_pos} -> g.{g} -> c.{c}"

    @pytest.mark.parametrize("exon_map", [EXON_MAP_PLUS, EXON_MAP_MINUS])
    @pytest.mark.parametrize("genomic_pos", [990, 1000, 1050, 1099, 2000, 2050, 2099, 2100, 2150])
    def test_genomic_round_trip(self, exon_map, genomic_pos):
        """Genomic position maps to cDNA and back unchanged."""
        c, _ = _genomic_to_cdna(genomic_pos, exon_map)
        g, _ = _cdna_to_genomic(c, exon_map)
        assert g == genomic_pos, f"g.{genomic_pos} -> c.{c} -> g.{g}"


# ---------------------------------------------------------------------------
# _format_coord tests
# ---------------------------------------------------------------------------

class TestFormatCoord:
    """Tests for _format_coord."""

    def test_coding_simple(self):
        tag = parse("NM_000207.3:c.1A>G")
        assert _format_coord(1, None, 200, tag) == "1"
        assert _format_coord(150, None, 200, tag) == "150"

    def test_3utr(self):
        tag = parse("NM_000207.3:c.*1A>G")
        assert _format_coord(201, None, 200, tag) == "*1"
        assert _format_coord(250, None, 200, tag) == "*50"

    def test_5utr(self):
        tag = parse("NM_000207.3:c.-5A>G")
        assert _format_coord(-5, None, 200, tag) == "-5"

    def test_offset(self):
        tag = parse("NM_000207.3:c.637+1G>T")
        assert _format_coord(637, None, 200, tag) == "*437+1"

    def test_range(self):
        tag = parse("NM_000207.3:c.123_456del")
        assert _format_coord(123, 256, 200, tag) == "123_*56"


# ---------------------------------------------------------------------------
# _format_genomic_coord tests
# ---------------------------------------------------------------------------

class TestFormatGenomicCoord:
    """Tests for _format_genomic_coord."""

    def test_simple(self):
        tag = parse("NM_000207.3:c.1A>G")
        assert _format_genomic_coord(1000, None, tag) == "1000"

    def test_offset(self):
        tag = parse("NM_000207.3:c.637+1G>T")
        assert _format_genomic_coord(1500, None, tag) == "1500+1"

    def test_range(self):
        tag = parse("NM_000207.3:c.123_456del")
        assert _format_genomic_coord(1000, 1333, tag) == "1000_1333"


# ---------------------------------------------------------------------------
# _format_hgvs tests
# ---------------------------------------------------------------------------

class TestFormatHgvs:
    """Tests for _format_hgvs."""

    def test_substitution_g(self):
        tag = parse("NM_000207.3:c.1A>G")
        suffix = _get_variant_suffix(tag)
        result = _format_hgvs("NC_000023.11:g.", "1000", suffix, tag)
        assert result == "NC_000023.11:g.1000A>G"

    def test_deletion(self):
        tag = parse("NM_000207.3:c.123_456del")
        suffix = _get_variant_suffix(tag)
        result = _format_hgvs("NC_000023.11:g.", "1000_1333", suffix, tag)
        assert result == "NC_000023.11:g.1000_1333del"

    def test_insertion(self):
        tag = parse("NM_000207.3:c.123_124insACGT")
        suffix = _get_variant_suffix(tag)
        result = _format_hgvs("NC_000023.11:g.", "1000_1001", suffix, tag)
        assert result == "NC_000023.11:g.1000_1001insACGT"

    def test_duplication(self):
        tag = parse("NM_000207.3:c.123dup")
        suffix = _get_variant_suffix(tag)
        result = _format_hgvs("NC_000023.11:g.", "1000", suffix, tag)
        assert result == "NC_000023.11:g.1000dup"

    def test_delins(self):
        tag = parse("NM_000207.3:c.123_456delinsACGT")
        suffix = _get_variant_suffix(tag)
        result = _format_hgvs("NC_000023.11:g.", "1000_1333", suffix, tag)
        assert result == "NC_000023.11:g.1000_1333delinsACGT"

    def test_inversion(self):
        tag = parse("NM_000207.3:c.123_456inv")
        suffix = _get_variant_suffix(tag)
        result = _format_hgvs("NC_000023.11:g.", "1000_1333", suffix, tag)
        assert result == "NC_000023.11:g.1000_1333inv"


# ---------------------------------------------------------------------------
# _get_variant_suffix tests
# ---------------------------------------------------------------------------

class TestVariantSuffix:
    """Tests for _get_variant_suffix."""

    def test_substitution(self):
        tag = parse("NM_000207.3:c.1A>G")
        assert _get_variant_suffix(tag) == "A>G"

    def test_deletion(self):
        tag = parse("NM_000207.3:c.123_456del")
        assert _get_variant_suffix(tag) == "del"

    def test_insertion(self):
        tag = parse("NM_000207.3:c.123_124insACGT")
        assert _get_variant_suffix(tag) == "insACGT"

    def test_dup_with_bases(self):
        tag = parse("NM_000207.3:c.123dupA")
        assert _get_variant_suffix(tag) == "dupA"

    def test_offset_substitution(self):
        tag = parse("NM_000207.3:c.637+1G>T")
        assert _get_variant_suffix(tag) == "G>T"

    def test_delins(self):
        tag = parse("NM_000207.3:c.123_456delinsACGT")
        assert _get_variant_suffix(tag) == "delinsACGT"

    def test_inv(self):
        tag = parse("NM_000207.3:c.123_456inv")
        assert _get_variant_suffix(tag) == "inv"


# ---------------------------------------------------------------------------
# _get_true_cdna_and_cds_length tests
# ---------------------------------------------------------------------------

class TestTrueCdna:
    """Tests for _get_true_cdna_and_cds_length."""

    def test_coding(self):
        tag = parse("NM_000207.3:c.1A>G")
        pos, cds_len, is5, is3 = _get_true_cdna_and_cds_length(tag, EXON_MAP_PLUS)
        assert pos == 1
        assert cds_len == 200
        assert is5 is False
        assert is3 is False

    def test_3utr(self):
        tag = parse("NM_000207.3:c.*1A>G")
        pos, cds_len, is5, is3 = _get_true_cdna_and_cds_length(tag, EXON_MAP_PLUS)
        assert pos == 201  # cds_len + 1
        assert is3 is True

    def test_5utr(self):
        tag = parse("NM_000207.3:c.-5A>G")
        pos, cds_len, is5, is3 = _get_true_cdna_and_cds_length(tag, EXON_MAP_PLUS)
        assert pos == -5
        assert is5 is True
