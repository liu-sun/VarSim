"""Internal module: NCBI Entrez data fetching helpers.

These functions are not part of the public API. They fetch MANE Select/Plus
Clinical records from NCBI's nucleotide and protein databases.

Results are cached via functools.lru_cache. API calls include retry with
exponential backoff to handle rate limiting during batch operations.
"""

import os
import functools
import time

from Bio import Entrez, SeqIO

from ._logging import get_logger

logger = get_logger(__name__)

Entrez.email = os.environ["EMAIL"]
Entrez.api_key = os.environ.get("API_KEY") or None


# Supported MANE transcript types and their canonical NCBI keywords.
_MANE_TYPES = {
    "mane select": "MANE Select",
    "mane plus clinical": "MANE Plus Clinical",
}


def _mane_keyword(mane: str) -> str:
    """Return the canonical NCBI keyword for a MANE transcript type.

    Raises
    ------
    ValueError
        If *mane* is not ``"MANE Select"`` or ``"MANE Plus Clinical"``
        (case-insensitive).
    """
    keyword = _MANE_TYPES.get(mane.strip().lower())
    if keyword is None:
        raise ValueError(
            "mane must be 'MANE Select' or 'MANE Plus Clinical', got {!r}".format(mane)
        )
    return keyword


def _retry_on_error(func, max_retries=5, base_delay=2.0):
    """Decorator: retry an Entrez API call with exponential backoff.

    If the function raises RuntimeError (NCBI backend failure), it will
    retry up to max_retries times with increasing delays.
    Errors are NOT cached — only successful results are cached.
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        last_error = None
        for attempt in range(max_retries):
            try:
                return func(*args, **kwargs)
            except RuntimeError as e:
                last_error = e
                if "Search Backend failed" not in str(e):
                    raise
                if attempt < max_retries - 1:
                    delay = base_delay * (2 ** attempt)
                    logger.warning("NCBI rate-limit retry %d/%d (%.1fs delay)", attempt + 1, max_retries, delay)
                    time.sleep(delay)
        raise last_error
    return wrapper


def _dont_cache_errors(cached_func):
    """Wrapper: if the cached lru_cache function raises, clear cache and re-raise.
    
    Exposes cache_clear() on the wrapper for test fixture cleanup.
    """
    @functools.wraps(cached_func)
    def wrapper(*args, **kwargs):
        try:
            return cached_func(*args, **kwargs)
        except Exception:
            cached_func.cache_clear()
            raise
    wrapper.cache_clear = cached_func.cache_clear
    return wrapper


@_dont_cache_errors
@functools.lru_cache(maxsize=16)
@_retry_on_error
def nm(gene: str, mane: str = "MANE Select"):
    """Fetch the MANE nucleotide record for a gene.

    Parameters
    ----------
    gene : str
        Gene symbol (e.g. "G6PD").
    mane : str
        MANE transcript type: ``"MANE Select"`` (default) or
        ``"MANE Plus Clinical"``. Case-insensitive.

    Returns a Bio.SeqRecord from GenBank format.

    >>> from varsim._fetch import nm
    >>> seq = nm("G6PD")  # doctest: +SKIP
    >>> seq.id  # doctest: +SKIP
    'NM_001360016.2'
    >>> len(seq.seq) > 1000  # doctest: +SKIP
    True
    """
    keyword = _mane_keyword(mane)
    logger.info("Fetching nucleotide record for %s (%s) from NCBI Entrez ...", gene, keyword)
    t0 = time.perf_counter()
    stream = Entrez.esearch(
        db="nucleotide",
        term=f'{gene}[Gene Name] AND "{keyword}"[Keyword]',
    )
    record = Entrez.read(stream)
    stream = Entrez.efetch(
        db="nucleotide", id=record["IdList"], rettype="gb", retmode="text"
    )
    # Some genes match more than one MANE record (e.g. GNAL); use the first.
    seqrecord = next(SeqIO.parse(stream, "genbank"))
    dt = time.perf_counter() - t0
    logger.debug("nm(%s) → %s (%.2fs)", gene, seqrecord.id, dt)
    return seqrecord


@_dont_cache_errors
@functools.lru_cache(maxsize=16)
@_retry_on_error
def np(gene: str, mane: str = "MANE Select"):
    """Fetch the protein record paired with the gene's MANE transcript.

    The protein accession is taken from the transcript's CDS /protein_id
    qualifier, which guarantees correct transcript–protein pairing (a
    gene-level protein search can match several MANE records).

    Parameters
    ----------
    gene : str
        Gene symbol (e.g. "G6PD").
    mane : str
        MANE transcript type: ``"MANE Select"`` (default) or
        ``"MANE Plus Clinical"``. Case-insensitive.

    Returns a Bio.SeqRecord from FASTA format.

    >>> from varsim._fetch import np
    >>> seq = np("G6PD")  # doctest: +SKIP
    >>> seq.id  # doctest: +SKIP
    'NP_001346945.1'
    >>> len(seq.seq) > 100  # doctest: +SKIP
    True
    """
    logger.info("Fetching protein record for %s ...", gene)
    t0 = time.perf_counter()
    transcript = nm(gene, mane=mane)
    protein_id = None
    for feature in transcript.features:
        if feature.type == "CDS":
            ids = feature.qualifiers.get("protein_id")
            if ids:
                protein_id = ids[0]
            break
    if protein_id is None:
        raise ValueError(f"No protein_id qualifier on the CDS feature of {transcript.id}")
    stream = Entrez.efetch(
        db="protein", id=protein_id, rettype="fasta", retmode="text"
    )
    seqrecord = SeqIO.read(stream, "fasta")
    dt = time.perf_counter() - t0
    logger.debug("np(%s) → %s (%.2fs)", gene, seqrecord.id, dt)
    return seqrecord


@_dont_cache_errors
@functools.lru_cache(maxsize=16)
@_retry_on_error
def nc(gene: str) -> str:
    """Fetch the genomic (NC_) accession for a gene's primary assembly.

    Returns the accession string (e.g. "NC_000011.10").

    >>> from varsim._fetch import nc
    >>> acc = nc("G6PD")
    >>> acc  # doctest: +SKIP
    'NC_000023.11'
    >>> acc.startswith("NC_")  # doctest: +SKIP
    True
    """
    logger.info("Fetching genomic accession for %s ...", gene)
    t0 = time.perf_counter()
    stream = Entrez.esearch(
        db="nucleotide",
        term=f'{gene}[Gene Name] AND "Primary Assembly"[Title] AND human[Organism]',
    )
    record = Entrez.read(stream)
    stream = Entrez.efetch(
        db="nucleotide", id=record["IdList"], rettype="acc", retmode="text"
    )
    # The search can match more than one assembly record; the first
    # accession corresponds to the gene's primary assembly.
    acc = stream.read().strip().split("\n")[0]
    dt = time.perf_counter() - t0
    logger.debug("nc(%s) → %s (%.2fs)", gene, acc, dt)
    return acc
