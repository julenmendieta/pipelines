reference = '/scratch/julen/chiadrop/CIMAdata/10X_barcodes_white_list_4M-with-alts-february-2016.txt'
otherfile = '/scratch/julen/chiadrop/data/outData/CIMAdata/fastqCheck/10xcode/SLX-20379.SIGAH3.HYHJ2DRXX_16bp_noDup.txt'
refK = set()
with open(reference, 'r') as f:
    for line in f:
        refK.add(line[:-1])
present = 0
with open(otherfile, 'r') as f:
    for line in f:
        if line[:-1] in refK:
            present += 1
print(present)