# nsLtpFinder
This script aims to identify putative nsLTP protein in plant proteomes.
nsLTPFinder takes as input a directory containing proteomes to analyze (FASTA format).
Analytical steps are:
  1. Perform HMMSEARCH using HMM models (PF14547.6, PF14368.6 and PF00234.22) (on whole proteomes)
  2. Identify proteins with this conserved motif: "C.(6,15)C.(6,80)CC.(8,29)C.C.(8,37)C.(4,25)C" (motifSearch on whole proteomes)
  3. Merge results from HMMSEARCH and motifSearch
  4. Extract protein sequences retained from step 3
  5. Run signalp on retained proteins
  6. Calculate protein properties (gravy, isoelectric point and molecular weight) on protein from step 3
  7. Extract top candidates (define as protein identified through both motifSearch and HMMSEARCH with signal peptide identified with signalP)
  8. run MEME on:
    - all candidate proteins for each species
    - only top candidate proteins for each species
    - only top candiate proteins of all species merged
  
Requirements:
Python modules:
  - pandas (tested on version 1.3.0)
  - biopython (tested on version 1.79)
  - argparse (tested on version 1.4.0)

External programs:
  - SignalP5 (won't work with previous versions)
  - HMMER (tested on version 3.3)
  - MEME (tested on version 5.3.3)
