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

# Authors

Josh Fass | josh.fass@choderalab.org
Nicole McNeer | mcneern@mskcc.org
Alex Kentsis | kentsisresearchgroup@gmail.com
John Chodera | john.chodera@choderalab.org
