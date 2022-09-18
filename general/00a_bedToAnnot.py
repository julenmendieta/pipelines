# -*- coding: utf-8 -*-

# The input should have the following format:
#############################################
# chr1    724805  725050  id1   
# chr1    725859  725955  id2   

# The output will ensure all IDs are unique

import getopt, sys, os.path
from collections import defaultdict

def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hi:o:", ["help", "input=",
                                                        "output="])
    except getopt.GetoptError as err:
        print(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    infile = None
    outfile = None

    for o, a in opts:
        if o in ("-h","--help"):
            usage()
            sys.exit()
        elif o in ("-i", "--input"):
            if os.path.isfile(a):
                infile = a
        elif o in ("-o", "--output"):
            outfile = a
        else:
            assert False, "Unhandled option"

   
    if infile is not None and outfile is not None:
        run(infile, outfile)
    else:
        usage()
    


def usage():
    print()


def run(infile, outfile):

    inf  = open(infile, 'r')
    outf = open(outfile,'w')

    seen = defaultdict(int)
    for linea in inf:
        linea_split = linea.split()
        chrom = linea_split[0]
        ini_pos = int(linea_split[1])
        fin_pos = int(linea_split[2])
        id1 = linea_split[3]
        

        outf.write(f"{id1}" + "\t" + chrom + "\t" +  str(ini_pos) + 
                    "\t" + str(fin_pos) + f'\t+\t0\n')
        seen[id1] += 1

    inf.close()
    outf.close()


if __name__ == "__main__":
    main()

