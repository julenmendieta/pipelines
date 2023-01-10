# -*- coding: utf-8 -*-

# Author: Pedro Furió Tarí
# Data: 30/05/2013
#

# The input should have the following format:
#############################################
# chr1    724805  725050  id1   ...
# chr1    725859  725955  id2   ...

# The output will be named the following way:
# chr1    pfurio     [feature]     724705  724804   .      +       .      [feature]_id="id1"
# chr1    pfurio     [feature]     724805  725050   .      +       .      [feature]_id="id2"
# chr1    pfurio     [feature]     725051  725150   .      +       .      [feature]_id="id3"

import getopt, sys, os.path
from collections import defaultdict

def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hi:o:f:", ["help", "input=",
                                                        "output=", "feature="])
    except getopt.GetoptError as err:
        print(err) # will print something like "option -a not recognized"
        usage()
        sys.exit(2)
    infile = None
    outfile = None
    feature = None

    for o, a in opts:
        if o in ("-h","--help"):
            usage()
            sys.exit()
        elif o in ("-i", "--input"):
            if os.path.isfile(a):
                infile = a
        elif o in ("-o", "--output"):
            outfile = a
        elif o in ("-f", "--feature"):
            feature = a
        else:
            assert False, "Unhandled option"

    if feature == None:
        feature = 'bedFile'

    if infile is not None and outfile is not None:
        run(infile, outfile, feature)
    else:
        usage()
    


def usage():
    print("\nUsage: python bed2gtf [options] <mandatory>")
    print("Options:")
    print("\t-h, --help:\n\t\t show this help message and exit")
    print("Mandatory:")
    print("\t-i, --input:\n\t\t File with the regions in bed format")
    print("\t-o, --output:\n\t\t Name of the gtf file output file. Directory where the file will be created should exist!")
    print("\t-f, --feature:\n\t\t (Optional) feature type name, e.g. Gene, Variation, Similarity")
    print("\n30/05/2013. Pedro Furió Tarí.\n")
    print("\n07/09/2022. Julen Mendieta Esteban.\n")

def run(infile, outfile, feature):

    inf  = open(infile, 'r')
    outf = open(outfile,'w')

    cont = 1
    seen = defaultdict(int)
    for linea in inf:
        linea_split = linea.split()
        chrom = linea_split[0]
        ini_pos = int(linea_split[1])
        fin_pos = int(linea_split[2])
        id1 = linea_split[3]
        

        #outf.write(chrom + "\tpfurio\tpeak\t" + str(ini_pos) + "\t" + str(fin_pos) + '\t.\t+\t.\tpeak_id "' + peak + '";\n')
        outf.write(chrom + "\tpfurio\t" + str(feature) + "\t" + str(ini_pos) + 
                    "\t" + str(fin_pos) + f'\t.\t+\t.\t{feature}_id "' +
                    f"{id1}-{seen[id1]}" + '";\n')
        seen[id1] += 1
        cont += 1

    inf.close()
    outf.close()


if __name__ == "__main__":
    main()

