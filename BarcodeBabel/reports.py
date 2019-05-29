# Plots and such

import matplotlib.pyplot as plt
import numpy as np

from BarcodeBabel.barcodes import compute_mz, trypsinize, min_mz, max_mz, compute_hydrophobicity_profile


def mz_hist(peptides):
    """Histogram of m/z ratios"""
    mz = np.array(list(map(lambda p: compute_mz(trypsinize(p)), peptides)))
    hist_heights = plt.hist(mz, bins=50, alpha=0.5)[0];
    plt.ylabel('probability density')
    plt.yticks([])
    plt.xlabel('m/z')
    plt.title('m/z for random trial peptides')

    # plt.vlines([min_mz, max_mz], 0, max(hist_heights)*1.1)
    plt.fill_between([min_mz, max_mz], [0] * 2, [max(hist_heights) * 1.1] * 2, color='green', alpha=0.3)
    # plt.fill_between([min(mz_for_random_peptides), min_mz], [0]*2, [max(hist_heights)]*2, color='grey', alpha=0.3)
    # plt.fill_between([max_mz, max(mz_for_random_peptides)], [0]*2, [max(hist_heights)]*2, color='grey', alpha=0.3)


def lengths_hist(peptides):
    """Histogram of sampled peptide lengths"""
    lengths = list(map(len, peptides))
    plt.bar(np.arange(max(lengths) + 1), np.bincount(lengths))
    plt.xlabel('length')
    plt.ylabel('# occurrences in sampled barcode library')
    plt.xticks(np.arange(max(lengths) + 2)[::5])

    #plt.savefig('sample-librar-length.png', dpi=300)


def mz_vs_length_scatter(peptides):
    """Scatter-plot of peptide length vs m/z"""
    lengths = list(map(len, peptides))
    mz = np.array(list(map(lambda p: compute_mz(trypsinize(p)), peptides)))

    plt.scatter(lengths, mz, s=1)
    plt.title('m/z vs. length for random trial peptides')
    plt.xlabel('length')
    plt.ylabel('m/z')
    # plt.hlines([min_mz, max_mz], min(lengths), max(lengths))
    plt.fill_between([min(lengths), max(lengths)], [max_mz] * 2, [min_mz] * 2, color='green', alpha=0.3)


def hydrophobicity_profile_summary(peptides):
    """Summarize "hydrophobicity profiles" of a collection of peptides"""
    hydrophobicity_profiles = list(map(compute_hydrophobicity_profile, peptides))
    min_hydrophobicities = list(map(np.min, hydrophobicity_profiles))
    max_hydrophobicities = list(map(np.max, hydrophobicity_profiles))

    ax = plt.subplot(111)
    plt.hist(min_hydrophobicities, bins=50, alpha=0.5, label='min(profile)');
    plt.hist(max_hydrophobicities, bins=50, alpha=0.5, label='max(profile)');
    plt.xlabel('hydrophobicity')
    plt.ylabel('# occurrences in sampled barcode library')
    plt.legend(loc='best')
    plt.title('sampled barcode library')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    #plt.savefig('sample-library-hydrophobicity.png', dpi=300)
