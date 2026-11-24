#!/usr/bin/env python3
"""Generate inputs for the pareto optimization tests.

Usage: gen_inputs.py <genome fasta> <output directory>
small.fa, mid.fa -- prefixes of the genome (about 40k and 400k k-mers)
random.fa        -- sparse random 31-mers (no overlaps, many runs, small index type)
prefix_a.fa      -- k-mers sharing a prefix (empty buckets of the search index)
mixed.fa(.gz)    -- multi-line FASTA with N's and lowercase letters, also gzipped
unitigs_*.fa     -- maximal unitigs of small.fa / mid.fa (some reverse complemented, in random order), which have to stay
                    substrings of the result: k31, k63, mid_k31, k10 (even k, self-complementary k-mers),
                    k15_split (unitigs cut into pieces), k31_u (unidirectional)
"""
import gzip
import os
import random
import sys

COMPLEMENT = str.maketrans("ACGT", "TGCA")


def reverse_complement(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def unitigs(genome, k, bidirectional=True):
    """Maximal non-branching paths of the de Bruijn graph of the k-mers of the genome."""
    canonical = (lambda x: min(x, reverse_complement(x))) if bidirectional else (lambda x: x)
    kmers = {canonical(genome[i:i + k]) for i in range(len(genome) - k + 1)}
    successors = lambda x: [x[1:] + b for b in "ACGT" if canonical(x[1:] + b) in kmers]
    predecessors = lambda x: [b + x[:-1] for b in "ACGT" if canonical(b + x[:-1]) in kmers]
    visited, result = set(), []
    for start in sorted(kmers):
        if start in visited:
            continue
        visited.add(start)
        paths = []
        for forward, backward in ((successors, predecessors), (predecessors, successors)):
            path, x = [], start
            while True:
                following = forward(x)
                if len(following) != 1 or len(backward(following[0])) != 1 or canonical(following[0]) in visited:
                    break
                x = following[0]
                visited.add(canonical(x))
                path.append(x)
            paths.append(path)
        path = paths[1][::-1] + [start] + paths[0]
        result.append(path[0] + "".join(x[-1] for x in path[1:]))
    return result


def shuffled_unitigs(genome, k, rng, bidirectional=True, split=False):
    result = []
    for unitig in unitigs(genome, k, bidirectional):
        while split and len(unitig) > k + 5 and rng.random() < 0.5: # Pieces of a unitig are unitigs too
            cut = rng.randrange(1, len(unitig) - k)
            result.append(unitig[:cut + k - 1])
            unitig = unitig[cut:]
        result.append(unitig)
    if bidirectional:
        result = [reverse_complement(u) if rng.random() < 0.5 else u for u in result]
    rng.shuffle(result)
    return result


def write_fasta(path, sequences, width=None):
    with open(path, "w") as f:
        for i, sequence in enumerate(sequences):
            f.write(f">{i}\n")
            if width:
                f.write("\n".join(sequence[j:j + width] for j in range(0, len(sequence), width)) + "\n")
            else:
                f.write(sequence + "\n")


def main():
    genome_path, out = sys.argv[1], sys.argv[2]
    os.makedirs(out, exist_ok=True)
    genome = "".join(line.strip() for line in open(genome_path) if not line.startswith(">"))
    write_fasta(f"{out}/small.fa", [genome[:40000]])
    write_fasta(f"{out}/mid.fa", [genome[:400000]])

    rng = random.Random(1)
    write_fasta(f"{out}/random.fa", ["".join(rng.choice("ACGT") for _ in range(31)) for _ in range(3000)])
    write_fasta(f"{out}/prefix_a.fa", ["A" + "".join(rng.choice("AC") for _ in range(30)) for _ in range(300)])

    mixed = "".join(rng.choice("ACGTacgtN") if rng.random() < 0.05 else rng.choice("ACGT") for _ in range(20000))
    write_fasta(f"{out}/mixed.fa", [mixed, mixed[:5000].lower()], width=60)
    with open(f"{out}/mixed.fa", "rb") as source, gzip.open(f"{out}/mixed.fa.gz", "wb") as target:
        target.write(source.read())

    write_fasta(f"{out}/unitigs_k31.fa", shuffled_unitigs(genome[:40000], 31, rng))
    write_fasta(f"{out}/unitigs_k63.fa", shuffled_unitigs(genome[:40000], 63, rng))
    write_fasta(f"{out}/unitigs_mid_k31.fa", shuffled_unitigs(genome[:400000], 31, rng))
    write_fasta(f"{out}/unitigs_k10.fa", shuffled_unitigs(genome[:40000], 10, rng))
    write_fasta(f"{out}/unitigs_k15_split.fa", shuffled_unitigs(genome[:40000], 15, rng, split=True))
    write_fasta(f"{out}/unitigs_k31_u.fa", shuffled_unitigs(genome[:40000], 31, rng, bidirectional=False))


if __name__ == "__main__":
    main()
