# BarcodeBabel

Python script for generation of peptide barcodes for LC-MS/MS

# Overview

BarcodeBabel is an algorithm to generate libraries of peptide barcodes with user-defined ranges of features for optimal detectability by nanoscale liquid chromatography tandem mass spectrometry (LC-MS/MS). User-selected rules include lack of presence in a user-provided proteome, m/z range, presence of homopolymers, hydrophobicity range, enzyme cleavage sites, residue frequencies, and library size. Sample human proteome with common contaminants is provided. 

# Manifest

NoteBooks: BarcodeBabel.ipynb - Jupiter notebook with BarcodeBabel script

BarcodeBabel: barcodes.py - BarcodeBabel script

Data: Sample human proteome and common contaminants, and example 48 test barcodes

# Dependencies 

matplotlib

numpy

scipy

sci-kitlearn

pandas

biopython

tqdm (progress bars)

xlrd

lxml

pyteomics (preferred way to install is pip Python package manager)

# Usage

1.	Select enzyme cleavage sites: def append_enzyme_cleavage_sites

•	Defaults: GRA on both ends, for trypsin cleavage

•	User may choose to specify separately for beginning and end

•	If the default is changed, edit enzyme function def trypsinize(peptide) to reflect cleavage 

2.	User defines frequencies of amino acids:

•	Defaults: Frequencies will be 0 for: K, R, H, C, M, P, I

•	Rationale: K, R, H may cause trypsin cleavage within sequences. M, C may be oxidized. P skews fragmentation. I and L cannot be distinguished. Q and N easily deamidated. 

3.	User defines m/z range at user-defined pH

•	Defaults: pH = 3 (for standard LC-MS at low pH), m/z range: 550-850

•	Peptide properties calculated for the enzyme-cleaved barcode 

4.	User defines maximum length of allowed repeats

•	Default: 2 repeats

5.	User defines hydrophobicity range

•	Based on the maximally hydrophobic 5-residue in candidate peptides

•	Default: -0.5 to 2.5, determined by examining library of tryptic peptides with favorable range of retention 

•	Hydrophobocity is calculated for the enzyme-cleaved barcode 

6.	User selects reference proteome, to disallow overlap: 

•	Provided proteome: human reference proteome and common contaminants

7.	User defines library size:

•	Default: 2000

•	Generates .csv file

8.	Default: BarcodeBabel also reports 5'-3' DNA sequence for synthesis, 3'-5' reverse complement, no overhangs, out-of-frame start / stop codons in DNA sequence

•	Included in the generated csv file

# Authors

Josh Fass | josh.fass@choderalab.org

Nicole McNeer | mcneern@mskcc.org

Alex Kentsis | kentsisresearchgroup@gmail.com

John Chodera | john.chodera@choderalab.org
