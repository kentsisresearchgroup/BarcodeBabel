import pandas as pd
import numpy as np
import matplotlib.pyplot as plt



from pkg_resources import resource_filename


def test_mz_matches_skyline():
    """Test that BarcodeBabel's m/z calculator agrees with Skyline's"""

    # load the excel file
    path = resource_filename('BarcodeBabel', 'data/48 test barcodes for synthesis.xlsx')
    excel = pd.read_excel(path)

    # get the peptide sequences and the m/z ratios predicted by Skyline
    peptides = excel[:48].peptide
    skyline_mz = excel[:48]['m/z - Skyline']

    # compute m/z ratios using BarcodeBabel
    from BarcodeBabel.barcodes import compute_mz
    barcode_babel_mz = list(map(compute_mz, peptides))

    # assert they are within some tolerance of each other
    assert(np.allclose(barcode_babel_mz, skyline_mz, rtol=1e-5))
