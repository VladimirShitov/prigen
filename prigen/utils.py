import tempfile

from Bio.Blast.Applications import NcbiblastnCommandline as BlastN


def gc_content(sequence: str) -> float:
    """Calculate GC-content of a nucleotide sequence"""
    sequence = sequence.upper()
    return (sequence.count("G") + sequence.count("C")) / len(sequence)


def parse_blast_result(result: str) -> set[str]:
    """Parse BLAST TSV result (format #6). Return IDs of found queries

    Parameters
    ----------
    result : string with BLAST result in TSV format. This string contains the following columns:
        1. query id
        2. subject id
        3. identity %
        4. length
        5. Mismatch
        6. Gapopen
        7. Query start
        8. Query end
        9. Subject start
        10. Subject end
        11. E-value
        12. Bitscore

    Returns
    -------
    found_queries : set[str]
        IDs of queries that are present in result file
    """

    found_queries = set()

    if not result:
        return found_queries

    for row in result.strip().split("\n"):
        found_queries.add(row.split("\t")[0])

    return found_queries


def remove_keys_from_dict(dictionary: dict, excess_keys) -> dict:
    """Remove `excess_keys` from `dictionary`

    Parameters
    ----------
    dictionary: dict
        A dictionary from that you want to filter out keys
    excess_keys: iterable
        Any iterable object (e.g list or set) with keys that you want to remove from `dictionary`

    Returns
    -------
    `dictionary` without keys that were in `excess_keys`
    """
    for key in excess_keys:
        dictionary.pop(key, None)

    return dictionary


def filter_primers_by_blast(
        primers: dict[str, float], blast_db_path=None, remote: bool = False
) -> dict[str, float]:
    """BLAST primers and remove the ones that were found in database

    In the current version, 0.001 e-value threshold is set, and no other statistics are
    supported. If the sequence is found with any score, it will be filtered out

    Parameters
    ----------
    primers: dict[str, float]
        Dictionary, where keys are nucleotide sequences, and values are melting temperatures
    blast_db_path:
        Path to the directory with blast datablase. Ignored if `remote` is set to True
    remote: bool = False
        If True, BLAST search will be done in remote database. If False, a local database
        in `blast_db_path` will be used

    Returns
    -------
    filtered_primers: dict[str, float]
        `primers` without keys that were found in blast database
    """
    with tempfile.NamedTemporaryFile("w") as f:
        for primer in primers.keys():
            f.write(f">{primer}\n")
            f.write(f"{primer}\n")

        if remote:
            cline = BlastN(
                query=f.name,
                evalue=0.001,
                outfmt=6,  # TSV
                remote=True,
                db="nt"
            )
        else:
            cline = BlastN(
                query=f.name,
                evalue=0.001,
                outfmt=6,  # TSV
                db=blast_db_path
            )

        result = cline()
        found_primers: set[str] = parse_blast_result(result[0])

        return remove_keys_from_dict(primers, found_primers)
