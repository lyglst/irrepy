import numpy as np
from isointer import *
from vaspinter import *

def tcart(pt, lattice):
    return np.matmul(pt, lattice)

def tred(pt, lattice):
    return np.matmul(pt, np.linalg.inv(lattice))

def close_up(x):
    disp = []
    for i in x:
        for j in i:
            if (j%1)>0.5:
                disp.append((j%1) -1)
            else:
                disp.append(j%1)
    return np.reshape(disp, np.shape(x))
    
def coord_mapping(pt, frac_coords, supercell, origin_shift = np.array([0, 0, 0])):
    mapping = []
    for ipt in pt:
        dist = np.array(list(map(np.linalg.norm, np.matmul(close_up(ipt + origin_shift - frac_coords),supercell))))
        mapping.append(dist.argmin())
    return mapping 

def gen_mode(lpt, lpv, frac_coords, lattice, lname, supercell, origin_shift = np.array([0,0,0])):
    lvec = []
    lnmode = []
    for i, pt in enumerate(lpt, 0):
        mapping = coord_mapping(np.matmul(pt, np.linalg.inv(supercell)), frac_coords, supercell, origin_shift)
        for k in range(len(lpv[i][0])):
            vec = np.zeros(frac_coords.shape)
            for m, im in enumerate(mapping, 0):
                vec[im] = lpv[i][m][k]
            lvec.append(vec/np.linalg.norm(vec))
            lnmode.append(lname[i]+'_'+str(k))
    return lnmode, lvec

if __name__=='__main__':
    print('test mode')
    lname, lpt, lpv = parse_iso('iso.log')
    print('there are {} irreps: {}'.format(len(lname), " ".join(lname)))
    print('each has dimension of : {}'.format("".join([str(len(pv[0]))+' ' for pv in lpv])))
    frac_coords, freq, normal_modes, vasprun = phonon_readin('vasprun.xml')
    print('there are {} modes'.format(len(freq)))
    lnmode, lvec = gen_mode(lpt, lpv, frac_coords, vasprun.structures[0].lattice.matrix, [[1,1,0], [-1,1,0], [0,0,1]])

    for i in range(134, 144):
        print('mode {} is composed of: \n'.format(i))
        for j,vec in enumerate(lvec):
            if np.dot(np.reshape(normal_modes[i], -1), np.reshape(vec, -1))>0.1:
                print('{} (no.{}) : {}'.format(lnmode[j], j, np.dot(np.reshape(normal_modes[i], -1), np.reshape(vec, -1))))
