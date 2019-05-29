import numpy as np
from Bio.SeqUtils.ProtParam import ProteinAnalysis

from BarcodeBabel.reference_proteome import not_in_reference_proteome

# Construct sampling distributions over amino acids
standard_amino_acids = sorted(list('AGILPVFWYDERHKSTCMNQ'))
index_of_aa = dict(zip(standard_amino_acids, range(len(standard_amino_acids))))

# uniform distribution
uniform_distribution = np.ones(len(standard_amino_acids)) / len(standard_amino_acids)


def sample_random_peptide(length=10, residue_frequencies=uniform_distribution):
    """Sample a random peptide of specified length, with specified AA probabilities"""
    return ''.join(np.random.choice(standard_amino_acids, size=length, p=residue_frequencies))


# proportions for starting random peptide
flycode_dict = {
    'A': 0.1,
    'S': 0.1,
    'T': 0.1,
    'N': 0.00,
    'Q': 0.00,
    'D': 0.1,
    'E': 0.1,
    'V': 0.1,
    'L': 0.1,
    'F': 0.1,
    'Y': 0.00,
    'W': 0.1,
    'G': 0.1,
    'P': 0.00,
}

flycode_distribution = np.zeros(len(standard_amino_acids))
for aa in flycode_dict:
    flycode_distribution[index_of_aa[aa]] = flycode_dict[aa]
assert (np.sum(flycode_distribution) == 1)

residues_to_avoid = ['K', 'R', 'H', 'M', 'C', 'P', 'Q', 'N', 'I']


def avoid_certain_residues(initial_distribution, residues_to_avoid):
    """Take an initial distribution, and set the probability for each element of residues_to_avoid to 0, then renormalize"""
    censored_distribution = np.array(initial_distribution)
    for aa in residues_to_avoid:
        censored_distribution[index_of_aa[aa]] = 0
    assert (np.sum(censored_distribution) > 0)
    return censored_distribution / np.sum(censored_distribution)


residue_frequencies = avoid_certain_residues(flycode_distribution, residues_to_avoid)
print('initial sampling distribution: ')
for i in range(len(standard_amino_acids)):
    print('\t{}: {:.3}%'.format(standard_amino_acids[i], residue_frequencies[i] * 100))


def compute_mz(peptide, pH=1):
    """Return the mass over charge ratio for a peptide at a specified pH"""
    from pyteomics.mass import calculate_mass
    from pyteomics.electrochem import charge
    return calculate_mass(peptide) / charge(sequence=peptide, pH=pH)


def check_longest_repetition(peptide):
    """Return the length of the longest string of repeated characters in the peptide"""
    prev_aa = peptide[0]
    longest_repetition = 1
    current_repetition = 1
    for i in range(1, len(peptide)):
        if peptide[i] == prev_aa:
            current_repetition += 1
        else:
            prev_aa = peptide[i]
            longest_repetition = max(current_repetition, longest_repetition)
            current_repetition = 1
    return longest_repetition


kyte_doolittle_scale = {'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5,
                        'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
                        'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6,
                        'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2}


def compute_hydrophobicity_profile(peptide, window=5):
    """Return the hydrophobicity profile along the peptide using the Kyte-Doolittle scale"""
    analysis = ProteinAnalysis(peptide)
    window = min(window, len(analysis.sequence))
    hydrophobicity_profile = analysis.protein_scale(kyte_doolittle_scale, window=window)
    return hydrophobicity_profile


def acceptable_hydrophobicity_profile(peptide, min_hydrophobicity=-2, max_hydrophobicity=2):
    """Return True if the MAXIMUM along the hydrophobicity profile is >= min_hydrophobicity,
    and the MAXIMUM along the hydrophobicity_profile is <= max_hydrophobicity"""
    hydrophobicity_profile = compute_hydrophobicity_profile(peptide)
    return (np.max(hydrophobicity_profile) >= min_hydrophobicity) and (
                np.max(hydrophobicity_profile) <= max_hydrophobicity)


def no_long_repetitions(peptide, max_n_repeat=2):
    """Return True only if the sequence contains no strings of repeated characters
    longer than threshold."""
    return check_longest_repetition(peptide) <= max_n_repeat


def mz_in_range(peptide, min_mz=550, max_mz=850, pH=1):
    """Return True only if m/z ratio is within desired interval"""
    mz = compute_mz(peptide, pH=pH)
    return (mz >= min_mz) and (mz <= max_mz)


def avoids_certain_residues(peptide, residues_to_avoid):
    """Return True only if the peptide contains none of the residues_to_avoid"""
    return len(set(peptide).intersection(set(residues_to_avoid))) == 0


def append_enzyme_cleavage_sites(peptide, initial='GRA', end='GRA'):
    """Add fixed sequences to the beginning and end of the barcode"""
    return ''.join([initial, peptide, end])


def trypsinize(peptide):
    """GRAxxxGRA --> AxxxGR"""
    split = peptide.split("GRA")
    assert (len(split) == 3)
    return ''.join(['A', split[1], 'GR'])


def sample_random_length(min_length=3, max_length=25):
    """Pick a random integer given inclusive bounds"""
    return np.random.randint(min_length, max_length + 1)


def sample_from_trial_distribution():
    """Sample a random length, then sample a random peptide of that length"""
    return append_enzyme_cleavage_sites(sample_random_peptide(sample_random_length(), residue_frequencies))


max_n_repeat = 2
min_mz = 500
max_mz = 800
pH = 1.0
min_hydrophobicity = -0.5
max_hydrophobicity = 2.5

constraints = [
    lambda peptide: no_long_repetitions(peptide, max_n_repeat=max_n_repeat),
    lambda peptide: mz_in_range(trypsinize(peptide), min_mz=min_mz, max_mz=max_mz, pH=pH),
    lambda peptide: acceptable_hydrophobicity_profile(trypsinize(peptide), min_hydrophobicity=min_hydrophobicity,
                                                      max_hydrophobicity=max_hydrophobicity),
    lambda peptide: not_in_reference_proteome(trypsinize(peptide)),
]


def satisfies_all_constraints(peptide):
    """Return True only if each function in the constraints list returns True"""
    for constraint in constraints:
        if not constraint(peptide):
            return False
    return True


from tqdm import tqdm


def generate_library(library_size, sample_from_trial_distribution, check_constraints, max_n_trials=int(1e9)):
    """Build a library by sampling peptides from the trial distribution, and adding them to the library
    if they satisfy all the constraints, and if they're not already in the library."""
    library = []
    trange = tqdm(range(1, max_n_trials + 1))
    for i in trange:
        trial = sample_from_trial_distribution()
        if trial not in library:
            if check_constraints(trial):
                library.append(trial)
        trange.set_postfix({'current library size': len(library),
                            'current yield percentage': 100.0 * len(library) / i
                            })
        if len(library) >= library_size:
            print('target library size achieved! terminating early')
            break
    if len(library) < library_size:
        print('target library size not achieved, returning the library assembled so far')
    return library


def describe_peptide(peptide):
    barcode = peptide
    tryptic = trypsinize(peptide)
    mz_of_tryptic_digest = compute_mz(trypsinize(peptide))
    maxhyd = np.max(compute_hydrophobicity_profile(trypsinize(peptide)))
    return ','.join([barcode, tryptic, str(mz_of_tryptic_digest), str(maxhyd)])


if __name__ == '__main__':
    # generate library of random peptides that satisfy all the constraints
    np.random.seed(0)
    library = generate_library(library_size=100,
                               sample_from_trial_distribution=sample_from_trial_distribution,
                               check_constraints=satisfies_all_constraints)
    lines = list(map(describe_peptide, sorted(library)))

    with open('barcode-library-{}-{}-descriptors.csv'.format(min_mz, max_mz), 'w') as f:
        f.writelines([
            'barcode,tryptic digest, m/z of tryptic digest,maximum hydrophobocity tryptic digest, TAT-barcode-HIS-fwd,TAT-barcode-HIS-rev\n'])
        f.writelines(['{}\n'.format(line) for line in lines])
