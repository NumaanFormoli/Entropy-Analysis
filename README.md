# Entropy-Analysis
Conducting entropy analysis to measure the variability between the N and C terminal of sars-cov-2.

## Reproduce

Procedure to reproduce results

1. download the nucleotide sequences from GISAID (`EPI_SET_240302kp`)
2. align sequences and isolate the N protein
3. translate the result to protein (we used MUSCLE alignment in Ugene)
4. convert the resulting fasta file to excel
  - can use `fasta_to_excel.py`
5. run the entropy scripts

## Data Cleaning Scripts

The other scripts here were used in the process of selecting the set of 1800 sequences:

- `remove_low_coverage.py`
  - remove sequences with more than 10 consecutive unknowns ("N") or deletions ("-")
  - only run on the N protein
- `remove_duplicates.py`
  - just remove any sequences which were duplicated by name in the fasta file
- `equalize_variant_frequency.py`
  - pick 150 for each named WHO variant
