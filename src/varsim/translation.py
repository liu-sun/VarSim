"""HGVS coding-to-protein translation.

Translates coding (c.) HGVS variant descriptions to their protein (p.)
consequences. Supports substitutions, deletions, insertions, delins,
duplications, and frameshift-causing indels.

Edge cases handled:
- Start-codon mutations → ``p.(M1?)``
- UTR variants → ``"non-coding"``
- Silent / missense / nonsense / frameshift / stop-loss classification
"""

from __future__ import annotations

import re
from typing import Dict, List

from Bio.Seq import Seq
from Bio.SeqUtils import seq3

from . import _fetch
from ._logging import get_logger
from .parser import HGVSTag, parse

logger = get_logger(__name__)

# ---------------------------------------------------------------------------
# One-letter → three-letter amino acid mapping (incl. stop)
# ---------------------------------------------------------------------------

_AA_1TO3: Dict[str, str] = {}
for _aa1 in "ACDEFGHIKLMNPQRSTVWY":
    _AA_1TO3[_aa1] = seq3(_aa1)
_AA_1TO3["*"] = "Ter"


# ====================================================================
# Internal helpers
# ====================================================================


def _extract_cds(seqrecord) -> str:
    """Return the CDS nucleotide string from a GenBank SeqRecord.

    Raises ValueError if no CDS feature exists.
    """
    for feature in seqrecord.features:
        if feature.type == "CDS":
            return str(feature.extract(seqrecord).seq)
    raise ValueError("No CDS feature found in record")


def _is_utr(tag: HGVSTag, cds_length: int) -> bool:
    """Return True if the variant position falls outside the CDS."""
    # 5' UTR: positions are negative (c.-N)
    if tag.start_pos < 0:
        return True
    # 3' UTR: positions are prefixed with * (c.*N)
    if re.search(r":c\.\*", tag.original):
        return True
    # Beyond CDS end
    if tag.end_pos and tag.end_pos > cds_length:
        return True
    if tag.start_pos > cds_length:
        return True
    return False


def _build_mutated_cds(cds_seq: str, tag: HGVSTag) -> str:
    """Apply the coding variant to *cds_seq* and return the mutated CDS.

    Parameters
    ----------
    cds_seq : str
        Original CDS sequence (all uppercase).
    tag : HGVSTag
        Parsed coding variant.

    Returns
    -------
    str
        Mutated CDS sequence.
    """
    vt = tag.variant_type
    pos = tag.start_pos - 1                     # 0-indexed start
    end = (tag.end_pos - 1) if tag.end_pos else pos  # 0-indexed end

    # --- substitution ---
    if vt == "substitution":
        ref = (tag.ref or "").upper()
        alt = (tag.alt or "").upper()
        # verify reference
        span = cds_seq[pos : pos + max(len(ref), 1)]
        if ref and span.upper() != ref:
            raise ValueError(
                f"Reference mismatch at c.{tag.start_pos}: "
                f"expected {ref}, found {span}"
            )
        return cds_seq[:pos] + alt + cds_seq[pos + len(ref or alt or "N"):]

    # --- deletion ---
    if vt == "deletion":
        return cds_seq[:pos] + cds_seq[end + 1:]

    # --- insertion ---
    if vt == "insertion":
        alt = (tag.alt or "").upper()
        # insert *after* position pos (between pos and pos+1 in 1-based)
        return cds_seq[:pos + 1] + alt + cds_seq[pos + 1:]

    # --- deletion-insertion ---
    if vt == "delins":
        alt = (tag.alt or "").upper()
        return cds_seq[:pos] + alt + cds_seq[end + 1:]

    # --- duplication ---
    if vt == "duplication":
        dup_region = cds_seq[pos : end + 1]
        return cds_seq[:end + 1] + dup_region + cds_seq[end + 1:]

    raise ValueError(f"Unsupported variant type for translation: {vt}")


def _determine_protein_effect(
    cds_seq: str,
    mut_cds: str,
    tag: HGVSTag,
) -> dict:
    """Compare original and mutated CDS to classify the protein effect.

    Returns a dictionary with keys:

    * ``position``  — 1-based amino acid position
    * ``ref_aa``    — original amino acid (one-letter)
    * ``alt_aa``    — alternate amino acid / fs marker (one-letter)
    * ``effect_type`` — one of ``"missense"``, ``"silent"``, ``"nonsense"``,
      ``"frameshift"``, ``"stop_loss"``, ``"extension"``
    * ``fs_length`` — (frameshift only) number of novel amino acids until stop
    * ``fs_first_aa`` — (frameshift only) first novel amino acid after the
      affected codon
    """
    codon_idx = (tag.start_pos - 1) // 3      # 0-indexed codon
    aa_pos = codon_idx + 1                     # 1-indexed AA position

    orig_full = str(Seq(cds_seq).translate(to_stop=False))
    mut_full = str(Seq(mut_cds).translate(to_stop=False))

    ref_aa = orig_full[codon_idx] if codon_idx < len(orig_full) else "?"
    mut_at_pos = mut_full[codon_idx] if codon_idx < len(mut_full) else "?"

    # ---- frameshift detection ----
    net_change = len(mut_cds) - len(cds_seq)
    is_fs = net_change % 3 != 0

    if is_fs:
        # Translate from the affected codon to the first in-frame stop
        fs_start = codon_idx * 3
        fs_aa = str(Seq(mut_cds[fs_start:]).translate(to_stop=True))

        if len(fs_aa) == 0:
            # Immediate stop (no novel AAs after mutation point)
            return {
                "position": aa_pos,
                "ref_aa": ref_aa,
                "alt_aa": "*",
                "effect_type": "nonsense",
            }

        # First novel AA  (skip ref_aa if the first translated AA is unchanged)
        first_novel = fs_aa[0]
        if first_novel == ref_aa and len(fs_aa) > 1:
            first_novel = fs_aa[1]

        return {
            "position": aa_pos,
            "ref_aa": ref_aa,
            "alt_aa": first_novel,
            "effect_type": "frameshift",
            "fs_length": len(fs_aa),
            "fs_first_aa": first_novel,
        }

    # ---- in-frame comparison ----
    ref_aa = orig_full[codon_idx] if codon_idx < len(orig_full) else "?"
    alt_aa = mut_full[codon_idx] if codon_idx < len(mut_full) else "?"

    # Walk forward to find the first *actual* difference (handles multi-codon
    # in-frame changes where the first codon might coincidentally be silent).
    scan_idx = codon_idx
    while scan_idx < min(len(orig_full), len(mut_full)):
        if orig_full[scan_idx] != mut_full[scan_idx]:
            ref_aa = orig_full[scan_idx]
            alt_aa = mut_full[scan_idx]
            aa_pos = scan_idx + 1
            break
        scan_idx += 1
    else:
        # All scanned positions identical — length difference only?
        if len(orig_full) != len(mut_full):
            bigger = "mut" if len(mut_full) > len(orig_full) else "orig"
            return {
                "position": aa_pos,
                "ref_aa": "?",
                "alt_aa": "?",
                "effect_type": "extension" if bigger == "mut" else "truncation",
            }
        # Truly silent
        return {
            "position": aa_pos,
            "ref_aa": orig_full[codon_idx] if codon_idx < len(orig_full) else "?",
            "alt_aa": orig_full[codon_idx] if codon_idx < len(orig_full) else "?",
            "effect_type": "silent",
        }

    # Classify
    if ref_aa == "*":
        effect_type = "stop_loss"
    elif alt_aa == "*":
        effect_type = "nonsense"
    else:
        effect_type = "missense"

    return {
        "position": aa_pos,
        "ref_aa": ref_aa,
        "alt_aa": alt_aa,
        "effect_type": effect_type,
    }


def _format_p_hgvs(
    protein_id: str,
    effect: dict,
    *,
    one_letter: bool = True,
) -> str:
    """Build a ``NP_xxxxx.x:p.(…)`` HGVS string from an effect dictionary."""
    aa_pos = effect["position"]
    ref_aa = effect["ref_aa"]
    alt_aa = effect["alt_aa"]
    etype = effect["effect_type"]

    if one_letter:
        fmt = lambda a: a
    else:
        fmt = lambda a: _AA_1TO3.get(a, a)

    # Start-codon mutation → uncertain
    if aa_pos == 1 and etype != "silent":
        if one_letter:
            return f"{protein_id}:p.(M1?)"
        return f"{protein_id}:p.(Met1?)"

    if etype == "silent":
        return f"{protein_id}:p.({fmt(ref_aa)}{aa_pos}=)"

    if etype == "frameshift":
        fs_len = effect.get("fs_length", 0)
        first_novel = effect.get("fs_first_aa", alt_aa)
        return f"{protein_id}:p.({fmt(ref_aa)}{aa_pos}{fmt(first_novel)}fs*{fs_len})"

    if etype == "nonsense":
        t = "Ter" if not one_letter else "*"
        return f"{protein_id}:p.({fmt(ref_aa)}{aa_pos}{t})"

    if etype == "stop_loss":
        return f"{protein_id}:p.({fmt(ref_aa)}{aa_pos}{fmt(alt_aa)})"

    if etype == "extension":
        return f"{protein_id}:p.({fmt(ref_aa)}{aa_pos}{fmt(alt_aa)}ext*)"

    # missense & catch-all
    return f"{protein_id}:p.({fmt(ref_aa)}{aa_pos}{fmt(alt_aa)})"


# ====================================================================
# Public API
# ====================================================================


def translate_variant(c_hgvs: str, gene: str) -> str:
    """Translate a single coding HGVS variant to its protein consequence.

    Parameters
    ----------
    c_hgvs : str
        Coding HGVS variant, e.g. ``"NM_000207.3:c.1A>G"``.
        A bare ``"c.1A>G"`` is also accepted (a placeholder accession is
        prepended so the parser can process it).
    gene : str
        Gene symbol such as ``"INS"`` or ``"G6PD"``.

    Returns
    -------
    str
        Protein HGVS consequence, e.g. ``"NP_000198.1:p.(M1?)"``.
        Returns ``"non-coding"`` for UTR variants.

    Raises
    ------
    ValueError
        If the variant cannot be translated (unsupported type, reference
        mismatch, or non-coding prefix).

    Examples
    --------
    >>> translate_variant("NM_000207.3:c.1A>G", "INS")  # doctest: +SKIP
    'NP_000198.1:p.(M1?)'

    >>> translate_variant("NM_000207.3:c.-59A>G", "INS")  # doctest: +SKIP
    'non-coding'
    """
    # Accept bare "c.X" by prepending a dummy accession for the parser
    original = c_hgvs.strip()
    if not re.search(r"[A-Z]{2}_\d+\.\d+:", original):
        original = "NM_000000.0:" + original

    logger.info("Translating: %s for %s", c_hgvs, gene)
    tag = parse(original)

    if tag.prefix not in ("c.", "n."):
        raise ValueError(
            f"Expected a coding variant (c./n.), got prefix '{tag.prefix}'"
        )

    # Fetch reference data
    seqrecord = _fetch.nm(gene)
    protein_seqrecord = _fetch.np(gene)
    protein_id = protein_seqrecord.id
    cds_seq = _extract_cds(seqrecord)

    # UTR guard
    if _is_utr(tag, len(cds_seq)):
        return "non-coding"

    # Mutate & translate
    mut_cds = _build_mutated_cds(cds_seq, tag)
    effect = _determine_protein_effect(cds_seq, mut_cds, tag)

    return _format_p_hgvs(protein_id, effect, one_letter=True)


def translate_variants(c_hgvs_list: List[str], gene: str) -> List[str]:
    """Batch-translate multiple c.HGVS variants for the same gene.

    More efficient than calling :func:`translate_variant` repeatedly
    because the NM / NP records are fetched only once.

    Parameters
    ----------
    c_hgvs_list : list of str
        Coding HGVS variant strings.
    gene : str
        Gene symbol.

    Returns
    -------
    list of str
        Protein HGVS strings in the same order as the input list.

    Examples
    --------
    >>> translate_variants(["NM_000207.3:c.1A>G", "NM_000207.3:c.3G>T"], "INS")  # doctest: +SKIP
    ['NP_000198.1:p.(M1?)', 'NP_000198.1:p.(M1?)']
    """
    seqrecord = _fetch.nm(gene)
    protein_seqrecord = _fetch.np(gene)
    protein_id = protein_seqrecord.id
    cds_seq = _extract_cds(seqrecord)
    cds_len = len(cds_seq)

    logger.info("Batch translating %d variants for %s", len(c_hgvs_list), gene)
    results: List[str] = []
    for raw in c_hgvs_list:
        original = raw.strip()
        if not re.search(r"[A-Z]{2}_\d+\.\d+:", original):
            original = "NM_000000.0:" + original

        tag = parse(original)
        if tag.prefix not in ("c.", "n."):
            raise ValueError(
                f"Expected coding variant, got prefix '{tag.prefix}': {raw}"
            )
        if _is_utr(tag, cds_len):
            results.append("non-coding")
            continue

        mut_cds = _build_mutated_cds(cds_seq, tag)
        effect = _determine_protein_effect(cds_seq, mut_cds, tag)
        results.append(_format_p_hgvs(protein_id, effect, one_letter=True))

    return results


def get_protein_effect(c_hgvs: str, gene: str) -> dict:
    """Detailed protein effect for a coding variant.

    Parameters
    ----------
    c_hgvs : str
        Coding HGVS variant string.
    gene : str
        Gene symbol.

    Returns
    -------
    dict
        Keys:

        * ``"p_hgvs_1letter"`` — one-letter protein HGVS
        * ``"p_hgvs_3letter"`` — three-letter protein HGVS
        * ``"effect_type"`` — one of ``"missense"``, ``"silent"``,
          ``"frameshift"``, ``"nonsense"``, ``"start_loss"``,
          ``"stop_loss"``, ``"extension"``, or ``"non-coding"``
        * ``"position"`` — 1-indexed amino acid position
        * ``"ref_aa"`` — original amino acid (one-letter)
        * ``"alt_aa"`` — alternate amino acid (one-letter)

    Examples
    --------
    >>> eff = get_protein_effect("NM_000207.3:c.4A>G", "INS")  # doctest: +SKIP
    >>> eff["effect_type"]  # doctest: +SKIP
    'missense'
    >>> "p.(Ala" in eff["p_hgvs_3letter"]  # doctest: +SKIP
    True
    """
    original = c_hgvs.strip()
    if not re.search(r"[A-Z]{2}_\d+\.\d+:", original):
        original = "NM_000000.0:" + original

    logger.info("Getting protein effect: %s for %s", c_hgvs, gene)
    tag = parse(original)

    if tag.prefix not in ("c.", "n."):
        raise ValueError(
            f"Expected a coding variant (c./n.), got prefix '{tag.prefix}'"
        )

    seqrecord = _fetch.nm(gene)
    protein_seqrecord = _fetch.np(gene)
    protein_id = protein_seqrecord.id
    cds_seq = _extract_cds(seqrecord)

    if _is_utr(tag, len(cds_seq)):
        return {
            "p_hgvs_1letter": "non-coding",
            "p_hgvs_3letter": "non-coding",
            "effect_type": "non-coding",
            "position": 0,
            "ref_aa": "",
            "alt_aa": "",
        }

    mut_cds = _build_mutated_cds(cds_seq, tag)
    effect = _determine_protein_effect(cds_seq, mut_cds, tag)

    # Promote position-1 non-silent effects to "start_loss"
    etype = effect["effect_type"]
    if effect["position"] == 1 and etype != "silent":
        etype = "start_loss"

    return {
        "p_hgvs_1letter": _format_p_hgvs(protein_id, effect, one_letter=True),
        "p_hgvs_3letter": _format_p_hgvs(protein_id, effect, one_letter=False),
        "effect_type": etype,
        "position": effect["position"],
        "ref_aa": effect["ref_aa"],
        "alt_aa": effect.get("fs_first_aa", effect["alt_aa"]),
    }
