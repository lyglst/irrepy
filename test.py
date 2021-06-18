import sys, getopt

def main(argv):
	filename = 'vasprun.xml'
	supercell = [[1,0,0], [0,1,0], [0,0,1]]
	kpoints = ['GM']
	try:
		opts, args = getopt.getopt(argv, "hi:k:s:")
	except getopt.GetoptError:
		print("""
			irrepy -i vasprun.xml -k kpoints list -s supercell vector
			e.g. irepy -i vasprun.xml -k 'GM','X','R' -s 1,0,0,0,1,0,0,0,1
			""")
	for opt, arg in opts:
		if opt=='-h':
			print('irrepy -i vasprun.xml -k kpoints list -s supercell vector')
		elif opt=='-i':
			filename = arg
		elif opt=='-k':
			kpoints = arg
		elif opt=='-s':
			supercell = arg
	print(filename)
	print(kpoints)
	print(supercell)	

if __name__=='__main__':
	main(sys.argv[1:])
