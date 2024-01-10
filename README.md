# BarcodeBabel

Python script for generation of peptide barcodes for electrospray LC-MS/MS.


# Overview

BarcodeBabel is an algorithm to generate libraries of peptide barcodes with user-defined ranges of features for optimal detectability by nanoscale liquid chromatography tandem mass spectrometry (LC-MS/MS). User-selected rules include enzyme cleavage sites, residue frequencies, pH, m/z range, hydrophobicity range, presence of homopolymers, desired library size, and lack of presence in a provided proteome. Sample human proteome with common contaminants is provided. 


# Manifest

NoteBooks: 	

BarcodeBabel_KRG_v8_script.py -- BarcodeBabel python script
		
Barcode_Babel_KRG_v8_script_notebook.ipynb - BarcodeBabel python script as Jupyter notebook 

BarcodeBabel_KRG_v8_notebook.ipynb - BarcodeBabel as Jupyter notebook 

Input File: 	

inputfile.csv - sample .CSV file containing parameters and default values

Data: 		

Sample human proteome and common contaminants, and example 48 test barcodes


# Dependencies 

sys

getopt

matplotlib

pandas 

numpy

tqdm (progress bars)

biopython

csv

sci-kitlearn

xlrd

lxml

pyteomics (preferred way to install is pip Python package manager)


# Usage

A. Create a csv file (e.g., inputfile.csv) that provides definitions for parameters (steps 1-7). Save this file in the same directory as the BarcodeBabel script.

1.	Select enzyme cleavage sites ('cs1' and 'cs2'), provide length of sites ('len': both should be same length), define target residue ('digest'), and define 	on which terminus of the target residue to cleave the peptide ('side': either n or c).

•	Defaults: GRA on both ends, len = 3, digest = R, side = c (c-terminal) for trypsin cleavage.

•	User may choose to specify separately for first and second cleavage sites.

2.	Define frequencies of amino acids (twenty standard amino acids).

•	Defaults: Frequencies are 0 for: K, R, H, C, M, P, I.

•	Rationale: K, R, H may cause trypsin cleavage within sequences. M, C may be oxidized. P skews fragmentation. I and L cannot be distinguished. Q and N 		easily deamidated. 

3.	Define m/z range at user-defined pH.

•	Defaults: pH = 3 (for standard LC-MS at low pH), m/z range: 550-850

•	Peptide properties calculated for the enzyme-cleaved barcode 

4.	Define maximum length of allowed repeats

•	Default: 2 repeats

5.	Define hydrophobicity range

•	Based on the maximally hydrophobic 5-residue in candidate peptides

•	Default: -0.5 to 2.5, determined by examining library of tryptic peptides with favorable range of retention 

•	Hydrophobocity is calculated for the enzyme-cleaved barcode 

6.	Select reference proteome, to disallow overlap: 

•	Provided proteome: human reference proteome and common contaminants.

7.	Define library size:

•	Default: 2000

•	Generates .csv file

B. Run BarcodeBabel.

1.	User can run the Python script (.py) through Mac Terminal. 

•	Navigate to the directory containing the script, input file, and desired proteome.

•	Run the script: python SCRIPT_NAME.py -i <inputfile>

•	Error messages provided to troubleshoot.

2.	User can beta-test/customize the code through two Jupyter notebooks.

•	BarcodeBabel_KRG_v8.ipynb: specify an input file in the notebook to assign values to parameters.

•	BarcodeBabel_KRG_v8_script_notebook: notebook version of the .py script, which can be easily edited to refine the .py script (must save as a new .py file).

C. BarcodeBabel outputs several files into the directory.

1.	output.csv

•	Included in output.csv: date/time; parameters and their assignments; barcode library; protease-digested peptide sequence; m/z of peptide; max hydrophobicity of peptide; 5'-3' DNA sequence and 3'-5' reverse complement for peptide synthesis (generated by using the Genscript table for commonly used codons for expression in E. Coli).

2.	m-z-for-random-trial-peptides.png

•	Distribution of peptide m/z for random trial peptides.
•	Green box highlights desired m/z range.

3.	m-z-vs.-length-for-random-trial-peptides.png

•	m/z vs. length for random trial peptides.
•	Green box highlights desired m/z range.

4.	sample-librar-length.png

•	Library distribution of peptide length.

5.	sample-librar-m-z.png

•	Library distribution of peptide m/z.
•	Green box highlights desired m/z range.

6.	sample-librar-hydrophobicity.png

•	Library distribution of peptide hydrophobicity.


# Authors

Josh Fass | josh.fass@choderalab.org

Nicole McNeer | mcneern@mskcc.org

Elliot Eton | eoe4001@med.cornell.edu

Alex Kentsis | kentsisresearchgroup@gmail.com

