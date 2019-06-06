# BarcodeBabel

BarcodeBabel is an algorithm to generate libraries of peptide barcodes with user-defined ranges of features for optimal detectability by nanoscale liquid chromatography tandem mass spectrometry (LC-MS/MS). User-selected rules include lack of presence in a user-provided proteome, m/z range, presence of homopolymers, hydrophobicity range, enzyme cleavage sites, residue frequencies, and library size. Sample human proteome with common contaminants is provided. 

Requirements:
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

Josh Fass and Nicole McNeer
