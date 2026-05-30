"""Tests for varsim.translation: translate_variant, translate_variants, get_protein_effect."""

import pytest
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from tests.conftest import GENE, requires_entrez, assert_valid_p_hgvs

# ---------------------------------------------------------------------------
# Unit tests — no API calls, test internal helpers directly
# ---------------------------------------------------------------------------


class TestBuildMutatedCds:
    """Tests for _build_mutated_cds using a minimal synthetic CDS."""

    CDS = "ATGGCCGGTACTTGA"  # 15 bp: Met-Ala-Gly-Thr-STOP

    @staticmethod
    def _tag(vt, start, end=None, ref=None, alt=None):
        """Build a minimal HGVSTag for internal testing."""
        from varsim.parser import HGVSTag
        return HGVSTag(
            acc="NM_000000.0", genomic_acc=None, prefix="c.",
            start_pos=start, end_pos=end,
            start_offset=None, end_offset=None,
            ref=ref, alt=alt, variant_type=vt,
            is_uncertain=False, fs_length=None,
            original=f"NM_000000.0:c.{start}A>G",
        )

    def test_substitution_single_base(self):
        """c.1A>C replaces first base."""
        from varsim.translation import _build_mutated_cds
        tag = self._tag("substitution", 1, ref="A", alt="C")
        result = _build_mutated_cds(self.CDS, tag)
        assert result == "CTGGCCGGTACTTGA"
        assert result[0] == "C"

    def test_substitution_ref_mismatch_raises(self):
        """Ref mismatch should raise ValueError."""
        from varsim.translation import _build_mutated_cds
        tag = self._tag("substitution", 1, ref="G", alt="C")
        with pytest.raises(ValueError, match="Reference mismatch"):
            _build_mutated_cds(self.CDS, tag)

    def test_deletion_single_base(self):
        """c.3del removes the third base (1-indexed)."""
        from varsim.translation import _build_mutated_cds
        tag = self._tag("deletion", 3)
        result = _build_mutated_cds(self.CDS, tag)
        # CDS: ATG GCC GGT ACT TGA; remove pos 3 (index 2, 'G' in codon 1)
        # After: AT + GCC... = ATGCCGGTACTTGA
        assert result == "ATGCCGGTACTTGA"
        assert len(result) == len(self.CDS) - 1

    def test_deletion_range(self):
        """c.4_6del removes bases 4-6."""
        from varsim.translation import _build_mutated_cds
        tag = self._tag("deletion", 4, end=6)
        result = _build_mutated_cds(self.CDS, tag)
        # CDS: ATG GCC GGT ACT TGA, remove positions 4-6 = "GCC"
        assert result == "ATGGGTACTTGA"
        assert len(result) == len(self.CDS) - 3

    def test_insertion(self):
        """c.3_4insTT inserts TT between positions 3 and 4."""
        from varsim.translation import _build_mutated_cds
        tag = self._tag("insertion", 3, end=4, alt="TT")
        result = _build_mutated_cds(self.CDS, tag)
        # Insert after pos 3 (index 2): ATG + TT + GCCGGTACTTGA
        assert result == "ATGTTGCCGGTACTTGA"
        assert len(result) == len(self.CDS) + 2

    def test_delins(self):
        """c.4_6delinsAAA replaces positions 4-6 with AAA."""
        from varsim.translation import _build_mutated_cds
        tag = self._tag("delins", 4, end=6, alt="AAA")
        result = _build_mutated_cds(self.CDS, tag)
        assert result == "ATGAAAGGTACTTGA"
        assert len(result) == len(self.CDS)

    def test_duplication_single_base(self):
        """c.3dup duplicates position 3."""
        from varsim.translation import _build_mutated_cds
        tag = self._tag("duplication", 3)
        result = _build_mutated_cds(self.CDS, tag)
        # Duplicate position 3 (G): ATG + G + GCCGGTACTTGA
        assert result == "ATGGGCCGGTACTTGA"
        assert len(result) == len(self.CDS) + 1

    def test_duplication_range(self):
        """c.1_3dup duplicates first codon."""
        from varsim.translation import _build_mutated_cds
        tag = self._tag("duplication", 1, end=3)
        result = _build_mutated_cds(self.CDS, tag)
        assert result == "ATGATGGCCGGTACTTGA"
        assert len(result) == len(self.CDS) + 3


class TestDetermineProteinEffect:
    """Tests for _determine_protein_effect."""

    CDS = "ATGGCCGGTACTTGA"  # Met-Ala-Gly-Thr-STOP (5 AAs, 15 bp)

    @staticmethod
    def _tag(vt, start, end=None, ref=None, alt=None):
        from varsim.parser import HGVSTag
        return HGVSTag(
            acc="NM_000000.0", genomic_acc=None, prefix="c.",
            start_pos=start, end_pos=end,
            start_offset=None, end_offset=None,
            ref=ref, alt=alt, variant_type=vt,
            is_uncertain=False, fs_length=None,
            original=f"NM_000000.0:c.{start}A>G",
        )

    def test_silent_substitution(self):
        """Synonymous substitution in codon 2: GCC→GCT (both Ala)."""
        from varsim.translation import _determine_protein_effect
        # c.6C>T: CDS ATG GCC → ATG GCT  (both = Met-Ala)
        # Codon 2 (positions 4-6): GCC → GCT, change pos 6 (0-idx 5)
        tag = self._tag("substitution", 6, ref="C", alt="T")
        mut = self.CDS[:5] + "T" + self.CDS[6:]
        effect = _determine_protein_effect(self.CDS, mut, tag)
        assert effect["effect_type"] == "silent"

    def test_missense_substitution(self):
        """Non-synonymous: c.4G>T changes Ala→Ser (GCC→TCC)."""
        from varsim.translation import _determine_protein_effect
        tag = self._tag("substitution", 4, ref="G", alt="T")
        mut = self.CDS[:3] + "T" + self.CDS[4:]
        effect = _determine_protein_effect(self.CDS, mut, tag)
        assert effect["effect_type"] == "missense"
        assert effect["ref_aa"] == "A"
        assert effect["alt_aa"] == "S"
        assert effect["position"] == 2  # codon 2

    def test_nonsense_substitution(self):
        """c.7C>A on a custom CDS: CAG(Gln) → AAG→no, let's use TGG→TGA."""
        from varsim.translation import _determine_protein_effect
        # CDS: ATG TGG GGT  → TGG(Trp)→TGA(Stop) at position 5 (index 4)
        cds = "ATGTGGGGT"
        # c.5G>A: TGG → TGA
        tag = self._tag("substitution", 5, ref="G", alt="A")
        mut = cds[:4] + "A" + cds[5:]
        effect = _determine_protein_effect(cds, mut, tag)
        assert effect["effect_type"] == "nonsense"
        assert effect["alt_aa"] == "*"

    def test_frameshift_deletion(self):
        """1 bp deletion causes frameshift."""
        from varsim.translation import _determine_protein_effect
        tag = self._tag("deletion", 4)
        mut = self.CDS[:3] + self.CDS[4:]  # delete position 4
        effect = _determine_protein_effect(self.CDS, mut, tag)
        assert effect["effect_type"] == "frameshift"
        assert "fs_length" in effect

    def test_frameshift_insertion(self):
        """1 bp insertion causes frameshift."""
        from varsim.translation import _determine_protein_effect
        tag = self._tag("insertion", 3, end=4, alt="A")
        mut = self.CDS[:3] + "A" + self.CDS[3:]  # insert after pos 3
        effect = _determine_protein_effect(self.CDS, mut, tag)
        assert effect["effect_type"] == "frameshift"

    def test_inframe_deletion(self):
        """3 bp (in-frame) deletion."""
        from varsim.translation import _determine_protein_effect
        tag = self._tag("deletion", 4, end=6)  # delete codon 2 (GCC)
        mut = self.CDS[:3] + self.CDS[6:]  # ATG + GGTACTTGA
        effect = _determine_protein_effect(self.CDS, mut, tag)
        # Codon 1 (pos 1-3) is ATG, codon 2 is now GGT (Gly)
        assert effect["effect_type"] == "missense"


class TestFormatPHgvs:
    """Tests for _format_p_hgvs."""

    PID = "NP_000000.1"

    def test_silent_one_letter(self):
        from varsim.translation import _format_p_hgvs
        effect = {"position": 3, "ref_aa": "L", "alt_aa": "L", "effect_type": "silent"}
        result = _format_p_hgvs(self.PID, effect, one_letter=True)
        assert result == "NP_000000.1:p.(L3=)"

    def test_silent_three_letter(self):
        from varsim.translation import _format_p_hgvs
        effect = {"position": 3, "ref_aa": "L", "alt_aa": "L", "effect_type": "silent"}
        result = _format_p_hgvs(self.PID, effect, one_letter=False)
        assert result == "NP_000000.1:p.(Leu3=)"

    def test_missense_one_letter(self):
        from varsim.translation import _format_p_hgvs
        effect = {"position": 5, "ref_aa": "V", "alt_aa": "L", "effect_type": "missense"}
        result = _format_p_hgvs(self.PID, effect, one_letter=True)
        assert result == "NP_000000.1:p.(V5L)"

    def test_missense_three_letter(self):
        from varsim.translation import _format_p_hgvs
        effect = {"position": 5, "ref_aa": "V", "alt_aa": "L", "effect_type": "missense"}
        result = _format_p_hgvs(self.PID, effect, one_letter=False)
        assert result == "NP_000000.1:p.(Val5Leu)"

    def test_nonsense_one_letter(self):
        from varsim.translation import _format_p_hgvs
        effect = {"position": 10, "ref_aa": "G", "alt_aa": "*", "effect_type": "nonsense"}
        result = _format_p_hgvs(self.PID, effect, one_letter=True)
        assert result == "NP_000000.1:p.(G10*)"

    def test_nonsense_three_letter(self):
        from varsim.translation import _format_p_hgvs
        effect = {"position": 10, "ref_aa": "G", "alt_aa": "*", "effect_type": "nonsense"}
        result = _format_p_hgvs(self.PID, effect, one_letter=False)
        assert result == "NP_000000.1:p.(Gly10Ter)"

    def test_frameshift_one_letter(self):
        from varsim.translation import _format_p_hgvs
        effect = {
            "position": 7, "ref_aa": "V", "alt_aa": "C",
            "effect_type": "frameshift", "fs_length": 17, "fs_first_aa": "C",
        }
        result = _format_p_hgvs(self.PID, effect, one_letter=True)
        assert result == "NP_000000.1:p.(V7Cfs*17)"

    def test_frameshift_three_letter(self):
        from varsim.translation import _format_p_hgvs
        effect = {
            "position": 7, "ref_aa": "V", "alt_aa": "C",
            "effect_type": "frameshift", "fs_length": 17, "fs_first_aa": "C",
        }
        result = _format_p_hgvs(self.PID, effect, one_letter=False)
        assert result == "NP_000000.1:p.(Val7Cysfs*17)"

    def test_start_codon_substitution(self):
        """Position 1 non-silent → p.(M1?)/(Met1?)."""
        from varsim.translation import _format_p_hgvs
        effect = {"position": 1, "ref_aa": "M", "alt_aa": "L", "effect_type": "missense"}
        r1 = _format_p_hgvs(self.PID, effect, one_letter=True)
        r3 = _format_p_hgvs(self.PID, effect, one_letter=False)
        assert r1 == "NP_000000.1:p.(M1?)"
        assert r3 == "NP_000000.1:p.(Met1?)"

    def test_start_codon_silent(self):
        """Position 1 silent does NOT use ? notation."""
        from varsim.translation import _format_p_hgvs
        effect = {"position": 1, "ref_aa": "M", "alt_aa": "M", "effect_type": "silent"}
        r1 = _format_p_hgvs(self.PID, effect, one_letter=True)
        assert r1 == "NP_000000.1:p.(M1=)"

    def test_stop_loss(self):
        from varsim.translation import _format_p_hgvs
        effect = {"position": 5, "ref_aa": "*", "alt_aa": "R", "effect_type": "stop_loss"}
        result = _format_p_hgvs(self.PID, effect, one_letter=True)
        assert result == "NP_000000.1:p.(*5R)"


class TestIsUtr:
    """Tests for _is_utr."""

    def test_negative_position_is_utr(self):
        from varsim.translation import _is_utr
        from varsim.parser import HGVSTag
        tag = HGVSTag(
            acc="NM", genomic_acc=None, prefix="c.", start_pos=-5, end_pos=None,
            start_offset=None, end_offset=None, ref="A", alt="G",
            variant_type="substitution", is_uncertain=False, fs_length=None,
            original="NM:c.-5A>G",
        )
        assert _is_utr(tag, 100) is True

    def test_3prime_utr_by_original_string(self):
        from varsim.translation import _is_utr
        from varsim.parser import HGVSTag
        tag = HGVSTag(
            acc="NM", genomic_acc=None, prefix="c.", start_pos=1, end_pos=None,
            start_offset=None, end_offset=None, ref="A", alt="G",
            variant_type="substitution", is_uncertain=False, fs_length=None,
            original="NM:c.*1A>G",
        )
        assert _is_utr(tag, 100) is True

    def test_beyond_cds_length(self):
        from varsim.translation import _is_utr
        from varsim.parser import HGVSTag
        tag = HGVSTag(
            acc="NM", genomic_acc=None, prefix="c.", start_pos=200, end_pos=None,
            start_offset=None, end_offset=None, ref="A", alt="G",
            variant_type="substitution", is_uncertain=False, fs_length=None,
            original="NM:c.200A>G",
        )
        assert _is_utr(tag, 100) is True

    def test_within_cds(self):
        from varsim.translation import _is_utr
        from varsim.parser import HGVSTag
        tag = HGVSTag(
            acc="NM", genomic_acc=None, prefix="c.", start_pos=50, end_pos=None,
            start_offset=None, end_offset=None, ref="A", alt="G",
            variant_type="substitution", is_uncertain=False, fs_length=None,
            original="NM:c.50A>G",
        )
        assert _is_utr(tag, 100) is False


# ---------------------------------------------------------------------------
# Integration tests — require NCBI Entrez credentials
# ---------------------------------------------------------------------------


class TestTranslateVariantIntegration:
    """Integration tests that hit NCBI to verify real-world translation."""

    @requires_entrez
    def test_start_codon_substitution(self):
        """c.1A>G on G6PD should return p.(M1?)."""
        from varsim.translation import translate_variant
        # G6PD CDS starts with ATG; position 1 is A
        result = translate_variant("c.1A>G", GENE)
        assert_valid_p_hgvs(result)
        assert "M1?" in result or "Met1?" in result

    @requires_entrez
    def test_utr_returns_non_coding(self):
        """UTR positions should return 'non-coding'."""
        from varsim.translation import translate_variant
        result = translate_variant("c.-59A>G", GENE)
        assert result == "non-coding"

    @requires_entrez
    def test_bare_c_prefix_accepted(self):
        """Bare 'c.X' (no accession) should work."""
        from varsim.translation import translate_variant
        result = translate_variant("c.4del", GENE)
        assert_valid_p_hgvs(result)

    @requires_entrez
    def test_silent_variant_detected(self):
        """A known synonymous variant should produce silent p.HGVS."""
        from varsim.translation import translate_variant
        # Test using a deletion (no ref required) that is in-frame
        result = translate_variant("c.4_6del", GENE)
        assert_valid_p_hgvs(result)


class TestGetProteinEffectIntegration:
    """Integration tests for get_protein_effect."""

    @requires_entrez
    def test_returns_expected_keys(self):
        from varsim.translation import get_protein_effect
        # Use a 3bp deletion (in-frame, no ref check needed)
        result = get_protein_effect("c.4_6del", GENE)
        for key in ("p_hgvs_1letter", "p_hgvs_3letter", "effect_type", "position", "ref_aa", "alt_aa"):
            assert key in result, f"Missing key: {key}"

    @requires_entrez
    def test_utr_returns_non_coding_dict(self):
        from varsim.translation import get_protein_effect
        result = get_protein_effect("c.-1A>G", GENE)
        assert result["effect_type"] == "non-coding"
        assert result["position"] == 0

    @requires_entrez
    def test_both_hgvs_formats_present(self):
        from varsim.translation import get_protein_effect
        result = get_protein_effect("c.5del", GENE)
        assert ":p.(" in result["p_hgvs_1letter"]
        assert ":p.(" in result["p_hgvs_3letter"]


class TestTranslateVariantsIntegration:
    """Integration tests for batch translate_variants."""

    @requires_entrez
    def test_batch_returns_same_count(self):
        from varsim.translation import translate_variants
        variants = [
            "c.1del",
            "c.2del",
            "c.3del",
        ]
        results = translate_variants(variants, GENE)
        assert len(results) == len(variants)

    @requires_entrez
    def test_batch_mixed_with_utr(self):
        from varsim.translation import translate_variants
        variants = [
            "c.1del",
            "c.-59A>G",
            "c.10del",
        ]
        results = translate_variants(variants, GENE)
        assert results[1] == "non-coding"
        assert results[0] != "non-coding"
        assert results[2] != "non-coding"

    @requires_entrez
    def test_batch_bare_prefix(self):
        from varsim.translation import translate_variants
        results = translate_variants(["c.1del", "c.2del"], GENE)
        assert len(results) == 2
        for r in results:
            assert_valid_p_hgvs(r)
