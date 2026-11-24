#!/usr/bin/env python3
"""Verify that a masked superstring represents exactly the k-mers of a FASTA file.

Usage: verify_ms.py <fasta> <masked superstring> <k> [-u] [-s]
Prints the result with |S|, the number of runs of ones and of zeros in the mask, exits with 1 if the sets differ.
Without -u, k-mers and their reverse complements are treated as identical.
With -s, also every input sequence (its parts between N's of length at least k) has to be a substring of the superstring
(or its reverse complement without -u), as is guaranteed for unitigs; not_substring counts those that are not.
"""
import gzip
import re
import sys

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def canonical(kmer, bidirectional):
    return min(kmer, kmer.translate(COMPLEMENT)[::-1]) if bidirectional else kmer


def read_fasta(path):
    opener = gzip.open if path.endswith(".gz") else open
    sequences, current = [], []
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if current:
                    sequences.append("".join(current))
                current = []
            else:
                current.append(line.strip())
    if current:
        sequences.append("".join(current))
    return sequences


def count_not_substrings(sequences, superstring, k, bidirectional):
    """Number of sequence parts of length >= k that do not occur in the superstring (in either orientation)."""
    parts = [part for sequence in sequences for part in re.split("[^ACGT]+", sequence.upper()) if len(part) >= k]
    by_first_kmer = {}
    for index, part in enumerate(parts):
        by_first_kmer.setdefault(part[:k], []).append((index, part))
        if bidirectional:
            complement = part.translate(COMPLEMENT)[::-1]
            by_first_kmer.setdefault(complement[:k], []).append((index, complement))
    found = set()
    for i in range(len(superstring) - k + 1):
        for index, part in by_first_kmer.get(superstring[i:i + k], ()):
            if index not in found and superstring.startswith(part, i):
                found.add(index)
    return len(parts) - len(found)


def main():
    arguments = [a for a in sys.argv[1:] if a not in ("-u", "-s")]
    fasta, ms_path, k = arguments[0], arguments[1], int(arguments[2])
    bidirectional = "-u" not in sys.argv
    check_substrings = "-s" in sys.argv

    expected = set()
    sequences = read_fasta(fasta)
    for sequence in sequences:
        sequence = sequence.upper()
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i + k]
            if set(kmer) <= set("ACGT"):
                expected.add(canonical(kmer, bidirectional))

    ms = "".join(read_fasta(ms_path))
    superstring = ms.upper()
    mask = [c.isupper() for c in ms]
    represented = {canonical(superstring[i:i + k], bidirectional)
                   for i in range(len(superstring) - k + 1) if mask[i]}

    runs = sum(1 for i, m in enumerate(mask) if m and (i == 0 or not mask[i - 1]))
    zeros = mask.count(False)
    ok = represented == expected
    substrings = ""
    if check_substrings:
        not_substring = count_not_substrings(sequences, superstring, k, bidirectional)
        ok = ok and not_substring == 0
        substrings = f" not_substring={not_substring}"
    print(f"ok={ok} |S|={len(superstring)} runs={runs} zeros={zeros} "
          f"missing={len(expected - represented)} extra={len(represented - expected)} n={len(expected)}{substrings}")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
