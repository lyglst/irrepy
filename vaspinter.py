import pymatgen
import numpy as np
from pymatgen.io.vasp import outputs
from pymatgen.io.vasp import inputs
import pymatgen.io.vasp.outputs
import math

def phonon_readin(filename='vasprun.xml'):
    """
    read in phonon calculation result from vasp
    
    Args:
    filename: the vasprun file

    Return:
    coord: atomic coordinates (reduced coords)
    freq: phonon mode frequency
    normal_modes: mode vector
    """
    try:
        output = outputs.Vasprun(filename)
    except BaseException:
        print('A converged vasprun.xml file is required')
    normal_modes = output.normalmode_eigenvecs
    freq = -np.sign(output.normalmode_eigenvals)*np.sqrt(np.abs(output.normalmode_eigenvals)) * 15.633302 * 33.35641
    return output.structures[0].frac_coords, freq, normal_modes, output

if __name__=='__main__':
    frac_coords, freq, normal_modes, vasprun = phonon_readin('vasprun.xml')
    print('there are {} modes'.format(len(freq)))
