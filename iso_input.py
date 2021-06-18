import pymatgen
import numpy as np
from pymatgen.io.vasp import outputs
from pymatgen.io.vasp import inputs
import pymatgen.io.vasp.outputs
import math
import spglib
import sympy 
import re
import os
from isointer import *
from vaspinter import *
from decomp import *

lattice_symbols = {
    'P': [[0, 0, 0]],
    'A': [[0, 0, 0], [0, 12, 12]],
    'B': [[0, 0, 0], [12, 0, 12]],
    'C': [[0, 0, 0], [12, 12, 0]],
    'I': [[0, 0, 0], [12, 12, 12]],
    'R': [[0, 0, 0]],
    'H': [[0, 0, 0], [16, 8, 8], [8, 16, 16]],
    'F': [[0, 0, 0], [0, 12, 12], [12, 0, 12], [12, 12, 0]]
}

def read_wyckoff_csv(filename):
    with open(filename) as wyckoff_file:
        return parse_wyckoff_csv(wyckoff_file)


def parse_wyckoff_csv(wyckoff_file):
    """Parse Wyckoff.csv
    There are 530 data sets. For one example:
    9:C 1 2 1:::::::
    ::4:c:1:(x,y,z):(-x,y,-z)::
    ::2:b:2:(0,y,1/2):::
    ::2:a:2:(0,y,0):::
    """

    rowdata = []
    points = []
    hP_nums = [433, 436, 444, 450, 452, 458, 460]
    for i, line in enumerate(wyckoff_file):
        if line.strip() == 'end of data':
            break
        rowdata.append(line.strip().split(':'))

        # 2:P -1  ::::::: <-- store line number if first element is number
        if rowdata[-1][0].isdigit():
            points.append(i)
    points.append(i)

    wyckoff = []
    for i in range(len(points) - 1):  # 0 to 529
        symbol = rowdata[points[i]][1]  # e.g. "C 1 2 1"
        if i + 1 in hP_nums:
            symbol = symbol.replace('R', 'H', 1)
        wyckoff.append({'symbol': symbol.strip()})

    # When the number of positions is larger than 4,
    # the positions are written in the next line.
    # So those positions are connected.
    for i in range(len(points) - 1):
        count = 0
        wyckoff[i]['wyckoff'] = []
        for j in range(points[i] + 1, points[i + 1]):
            # Hook if the third element is a number (multiplicity), e.g.,
            #
            # 232:P 2/b 2/m 2/b:::::::  <- ignored
            # ::8:r:1:(x,y,z):(-x,y,-z):(x,-y+1/2,-z):(-x,-y+1/2,z)
            # :::::(-x,-y,-z):(x,-y,z):(-x,y+1/2,z):(x,y+1/2,-z)  <- ignored
            # ::4:q:..m:(x,0,z):(-x,0,-z):(x,1/2,-z):(-x,1/2,z)
            # ::4:p:..2:(0,y,1/2):(0,-y+1/2,1/2):(0,-y,1/2):(0,y+1/2,1/2)
            # ::4:o:..2:(1/2,y,0):(1/2,-y+1/2,0):(1/2,-y,0):(1/2,y+1/2,0)
            # ...
            if rowdata[j][2].isdigit():
                pos = []
                w = {'letter': rowdata[j][3].strip(),
                     'multiplicity': int(rowdata[j][2]),
                     'site_symmetry': rowdata[j][4].strip(),
                     'positions': pos}
                wyckoff[i]['wyckoff'].append(w)

                for k in range(4):
                    if rowdata[j][k + 5]:  # check if '(x,y,z)' or ''
                        count += 1
                        pos.append(rowdata[j][k + 5])
            else:
                for k in range(4):
                    if rowdata[j][k + 5]:
                        count += 1
                        pos.append(rowdata[j][k + 5])

        # assertion
        for w in wyckoff[i]['wyckoff']:
            n_pos = len(w['positions'])
            n_pos *= len(lattice_symbols[wyckoff[i]['symbol'][0]])
            assert n_pos == w['multiplicity']

    return wyckoff

def iso_input(filename, spacegroup, supercell, kpt, wy, xyz):
    with open(filename+'.in', 'w+') as f:
        f.write('v par {} \n'.format(spacegroup))
        f.write('v cell {} {} {} \n'.format(",".join(list(map(str,map(int, supercell[0])))),
        ",".join(list(map(str,map(int, supercell[1])))),
        ",".join(list(map(str,map(int, supercell[2]))))))
        f.write('v wy {} \n'.format(wy))
        f.write('v wy xyz {} {} {} \n'.format(xyz[0], xyz[1], xyz[2]))
        f.write('sh irrep \nsh micro vec \n')
        for ikpt in kpt:
            f.write('v kp {}\n'.format(ikpt))
            f.write('d dist\n')
            f.write('\n'*30)
        f.write('quit')


def gen_iso_input(filename, supercell, kpt):
    try:
        output = outputs.Vasprun(filename)
    except BaseException:
        print('A converged vasprun.xml file is required')
    input_list = []
    Ninput = 0
    Natom = len(output.structures[0].species)
    species = np.unique(output.structures[0].species)
    cell = (output.structures[0].lattice.matrix,
            output.structures[0].frac_coords,
            np.zeros(Natom))
    for i, element in enumerate(species, 0):
            cell[2][np.array(output.structures[0].species)==element] = i
    dataset = spglib.get_symmetry_dataset(cell, symprec=1e-5)
    print('v par {}'.format(dataset['number']))
    newcell = spglib.standardize_cell(cell)
    newdataset = spglib.get_symmetry_dataset(newcell)
    pattern_coor = re.compile(r"(?<=\().*?(?=\))")
    
    for wy in np.unique(newdataset['equivalent_atoms']):
        #print('v wy {}'.format(newdataset['wyckoffs'][wy]))
        for info_wy in read_wyckoff_csv('Wyckoff.csv')[dataset['hall_number']-1]['wyckoff']:
            if info_wy['letter'] == newdataset['wyckoffs'][wy]:
                wy_coor = pattern_coor.search(info_wy['positions'][0]).group().split(',')
                x, y, z = sympy.symbols('x, y, z', real=True)
                eq1 = sympy.sympify(wy_coor[0], locals={'x': x, 'y': y, 'z':z})
                eq2 = sympy.sympify(wy_coor[1], locals={'x': x, 'y': y, 'z':z})
                eq3 = sympy.sympify(wy_coor[2], locals={'x': x, 'y': y, 'z':z})
                for try_pt in newdataset['std_positions'][newdataset['equivalent_atoms']==wy]:
                    result = sympy.solve([eq1-try_pt[0], eq2-try_pt[1], eq3-try_pt[2]])
                    if len(result)>0:
                        break
        xyz = {x:0, y:0, z:0}
        try:
            for i in result.keys():
                xyz[i]=result[i]
        except BaseException:
            pass
        #print('v wy xyz {} {} {}'.format(xyz[x], xyz[y], xyz[z]))
        Ninput += 1 
        input_list.append('iso_'+str(Ninput)+'_'+species[newdataset['std_types'][wy]].name
                            +'_'+newdataset['wyckoffs'][wy])
        iso_input(input_list[-1], dataset['number'], supercell, kpt, newdataset['wyckoffs'][wy], [xyz[x], xyz[y], xyz[z]])
    return input_list


if __name__ == "__main__":
    
    supercell = [[2,0,0],[0,2,0],[0,0,2]]
    filename = 'vasprun_Ba3GeO.xml'
    kpt = ['GM', 'M', 'X', 'R']

    input_list = gen_iso_input(filename, supercell, kpt)
    print(input_list)

    lname = []
    lpt = []
    lpv = []

    for file in input_list:
        os.system("iso <"+file+".in >"+file+".log")
        tname, tpt, tpv = parse_iso(file+".log")
        lname += list(map(lambda x: x+'_'+file.split('_')[2], tname))
        lpt += tpt
        lpv += tpv
    
    print('there are {} irreps: {}'.format(len(lname), " ".join(lname)))
    print('each has dimension of : {}'.format("".join([str(len(pv[0]))+' ' for pv in lpv])))

    frac_coords, freq, normal_modes, vasprun = phonon_readin(filename)
    lnmode, lvec = gen_mode(lpt, lpv, frac_coords, vasprun.structures[0].lattice.matrix, lname, supercell)
    print('total irreps: {}'.format(len(lvec)))
    print('there are {} modes'.format(len(freq)))


    for i in range(len(normal_modes)):
        print('mode {0} ({1:2.4f}cm-1) is composed of:'.format(i, freq[i]))
        for j,vec in enumerate(lvec):
            if np.abs(np.dot(np.reshape(normal_modes[i], -1), np.reshape(vec, -1)))>0.01:
                print('{0} (no.{1}) : {2:2.4f}'.format(lnmode[j], j, np.dot(np.reshape(normal_modes[i], -1), np.reshape(vec, -1))))
        print('\n')

    
