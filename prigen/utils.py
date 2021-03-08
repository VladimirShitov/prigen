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

    for row in result.strip().split("\n"):
        found_queries.add(row.split("\t")[0])

    return found_queries
