# nsLtpFinder
This script aims to identify putative nsLTP protein in plant proteomes.
nsLTPFinder takes as input a directory containing proteomes to analyze (FASTA format).
Analytical steps are:
  1. Perform HMMSEARCH using HMM models (PF14547.6, PF14368.6 and PF00234.22) (on whole proteomes)
  2. Identify proteins with this conserved motif: "C.{6,14}C.{8,19}CC.{9,20}C.C.{12,34}C.{5,14}C" (motifSearch on whole proteomes)
  3. Merge results from HMMSEARCH and motifSearch
  4. Extract protein sequences retained from step 3
  5. Run signalp on retained proteins
  6. Calculate protein properties (gravy, isoelectric point and molecular weight) on protein from step 3
  7. Extract top candidates (define as protein identified through both motifSearch and HMMSEARCH with signal peptide identified with signalP)
  8. run MEME on:
    - all candidate proteins for each species
    - only top candidate proteins for each species
    - only top candiate proteins of all species merged
  
