import re

def parse_iso(filename='iso.log'):
    """ parse the isotropy output file

    Args:
        filename: the isotropy output file name

    Returns:
        lname: list of irreps
        lpt: list of atom coordinate
        lpv: list of distortion vectors, might be multi-dimensional
    """

    #read in the isotropy output
    try:
        with open(filename,'r') as f:
            read_data = f.read()
    except BaseException:
        print('the output of isotropy is required here')
        return
    #parse the isotropy output
    #pt - atom coordinates (kind of weird definition, pt = original reduced coordinate * supercell matrix)
    #pv - distortion vectors
    #lpt, lpv - list of wy, pt, pv
    #lname - name of modes
    #nmode - number of modes
    nmode = 0
    lname = []
    lpt = []
    lpv = []

    pattern_name = re.compile(r"^[A-Z0-9\+\-]+(?=\s)")
    pattern_coor = re.compile(r"(?<=\().*?(?=\))")
    pattern_vec = re.compile(r"(?<=\()[0-9,\.\-]*(?=\))")

    for line in read_data.split('\n'):
        if pattern_name.search(line):
            if nmode>0:
                lpt.append(pt)
                lpv.append(pv)
            pt = []
            pv = []
            nmode += 1 
            lname.append(pattern_name.search(line).group())
        if nmode==0:
            continue
        if re.search(r"Irrep|Enter", line):
            continue
        find = pattern_coor.findall(line)
        find2 = pattern_vec.findall(line)
        if (len(find)!=len(find2)):
            npv = 0
            for element in find:
                coor = list(map(float, element.split(',')))
                if npv==0:
                    pt.append(coor)
                if npv==1:
                    pv.append([coor])
                if npv>1:
                    pv[-1].append(coor)
                npv += 1
        else:
            for element in find:
                coor = list(map(float, element.split(',')))
                if npv==1:
                    pv.append([coor])
                if npv>1:
                    pv[-1].append(coor)
                npv += 1
        
    lpt.append(pt)
    lpv.append(pv)
    return lname, lpt, lpv

if __name__=='__main__':
    print('test mode')
    lname, lpt, lpv = parse_iso('iso.log')
    print('there are {} irreps: {}'.format(len(lname), " ".join(lname)))
    print('each has dimension of : {}'.format("".join([str(len(pv[0]))+' ' for pv in lpv])))
    