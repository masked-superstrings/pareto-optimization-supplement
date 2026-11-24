# Pareto optimization of masked superstrings

High-level details about of Pareto optimization algorithm are provided in [the paper](https://www.biorxiv.org/content/10.64898/2026.03.18.712440) and its supplement.

More details on the algorithm can be found in [Chapter 3, section 2](https://dspace.cuni.cz/bitstream/handle/20.500.11956/204030/130441237.pdf?sequence=1&isAllowed=y#section.3.2)
of the bachelor thesis of Ján Plachý.
More details on implementation can be found in [Chapter 4](https://dspace.cuni.cz/bitstream/handle/20.500.11956/204030/130441237.pdf?sequence=1&isAllowed=y#chapter.4)
of the thesis.

Pareto optimization was implemented in fork of KmerCamel🐫, see [original README](./KMERCAMEL_README.md) in the time of forking for underlying details of the framework.

## Options of `kmercamel compute -a pareto`

- `-O runs|zeros` -- objective, `-p INT` -- penalty (default 7 for runs, 2 for zeros).
- `-t INT` -- number of threads (default 1). Searches of a batch run in parallel and are committed sequentially;
  with more than one thread the result is not deterministic (it differs from the single-threaded one only negligibly).
- `-m` -- low memory mode: no skip pointers, no failure cache and a 4 times smaller search table
  (about 30 % less memory, 5 bytes + 4 bytes + ~1.5 bytes per k-mer and its reverse complement).
  It never changes the result, only trades memory for speed; the slowdown depends on the data
  (about 12 % on the N. gonorrhoeae pangenome, about 2.2 times on the SARS-CoV-2 pangenome). Combine with `-t`.

## Phases of the computation

0. Chains of leaves start at the starts of the input sequences, in their direction: the first k-mer and the next base
   of every sequence (between N's) are kept while reading, the first step goes to that base and further ones are greedy
   as in phase 1. If the input sequences are unitigs (or their pieces), no chain starts or ends inside a unitig,
   so every unitig (or its reverse complement) is a substring of the result. Other inputs (e.g. simplitigs) are
   not guaranteed to stay whole; sequences that cannot be followed from their start are reported in a warning.
1. Leaves (k-mers) are chained along overlaps of length k - 1.
2. Iterative deepening -- in every iteration, each uncompleted leaf is searched for a connection cheaper than the current budget.
3. When only a few leaves remain but many budgets are left (typical for large penalties), every remaining leaf is searched
   once exactly (best-first), and the connections are made in the order of their cost. With the default penalties
   this phase is usually not reached.

Building: `make` (`make OMPFLAGS=` disables OpenMP and thus `-t`; on x86_64 the `popcnt` instruction is used, override with `ARCHFLAGS=`).

## Tests

`make pareto-test` checks on small inputs (generated from `data/spneumoniae.fa`) that the results represent the k-mer sets
in all modes, that the reported objective is correct, that `-m` does not change the result, large penalties,
input formats, rejection of invalid options, and that unitigs (generated in the test) stay substrings of the result
without worsening the objective. `make pareto-test-sanitize` runs the same with AddressSanitizer
and UndefinedBehaviorSanitizer. Setting `REF_BIN=<binary>` also checks that the results are identical to another build (builds before phase 0 give
different results for inputs with sequences longer than k).
