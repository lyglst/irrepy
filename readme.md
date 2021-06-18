# Irrepy

## Introduction
Irrepy can take into phonon calculation result from vasp (isif=5,6,7,8) and decompose the phonon modes with irrep modes.


## Installation
Please install spglib, pymatgen and isotropy before use this:

	spglib:
	pip install --user spglib
	
	pymatgen:
	conda install --yes numpy scipy matplotlib
	pip install pymatgen
	
	isotropy:
	check https://stokes.byu.edu/iso/isolinux.php
	after installation, add isotropy into your system path. THIS IS IMPORTANT! Can be checked by this way:
	Type 'iso' into your command line, ok if you can see an interactive interface.

## Usage
After finish this installation, you can use:

	irrepy -i vasprun.xml -k kpoints list -s supercell vector -o output directory
	-i : The vasprun filename.
	-k : The kpoints being processed.
	-s : The 3x3 matrix that you use to generate the supercell, please input it as a 9-dimensional vector.
	-o : Output the irrep modes into a directory.

	e.g. irrepy -i vasprun_Ba3GeO.xml -k 'GM','X','R','M' -s 2,0,0,0,2,0,0,0,2 -o irrep_modes
