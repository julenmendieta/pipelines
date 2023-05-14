#GPT translated code
#!/usr/bin/env python3

import sys
import os
import getopt
import subprocess

def main(argv):
    # Initializing variables with default values
    chipBG = None
    inputBG = None
    chipBGstat = None
    inputBGstat = None
    BGout = None
    fai = None
    tmpDir = '.'
    scoretype = 'KmetStat'

    # Parsing command line arguments
    try:
        opts, args = getopt.getopt(argv, "t:c:sT:sC:out:faidx:tmpDir:s:s")
    except getopt.GetoptError:
        print("Error: invalid input argument(s)")
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-t':
            chipBG = arg
        elif opt == '-c':
            inputBG = arg
        elif opt == '-sT':
            chipBGstat = arg
        elif opt == '-sC':
            inputBGstat = arg
        elif opt == '-out':
            BGout = arg
        elif opt == '-faidx':
            fai = arg
        elif opt == '-tmpDir':
            tmpDir = arg
        elif opt == '-s':
            scoretype = arg

    # Checking if required files are present
    if not os.path.isfile(chipBG):
        print("Error: %s file not found" % chipBG)
        sys.exit(1)
    if not os.path.isfile(inputBG):
        print("Error: %s file not found" % inputBG)
        sys.exit(1)
    if os.path.isfile(BGout):
        print("Error: %s already exists" % BGout)
        sys.exit(1)
    if not os.path.isfile(fai):
        print("Error: %s file not found" % fai)
        sys.exit(1)

    # Checking if the required tool is available
    try:
        subprocess.check_call(['which', 'bedGraphToBigWig'])
    except subprocess.CalledProcessError:
        print('Error: bedGraphToBigWig from UCSC not available!!')
        sys.exit(1)

    # Setting up some variables
    BGKmet = os.path.splitext(BGout)[0] + 'KMETChromosomes.bedgraph'

    # Checking if required files are present
    if not os.path.isfile(chipBGstat):
        print("Error: %s file not found" % chipBGstat)
        sys.exit(1)
    if not os.path.isfile(inputBGstat):
        print("Error: %s file not found" % inputBGstat)
        sys.exit(1)

    chipVal = {}
    inputVal = {}
    get_kmet_stat_vals(chipBGstat, chipVal)
    get_kmet_stat_vals(inputBGstat, inputVal)

    H4me3A_IPe = chipVal['KmetStat_H3K4me3_A'] / inputVal['KmetStat_H3K4me3_A']
    H4me3B_IPe = chipVal['KmetStat_H3K4me3_B'] / inputVal['KmetStat_H3K4me3_B']
    H4me2A_IPe = chipVal['KmetStat_H3K4me2_A'] / inputVal['KmetStat_H3K4me2_A']
    H4me2B_IPe = chipVal['KmetStat_H3K4me2_B'] / inputVal['KmetStat_H3K4me2_B']

    IPe_H3K4me3 = (chipVal['KmetStat_H3K4me3_A'] + chipVal['KmetStat_H3K4me3_B']) / (inputVal['KmetStat_H3K4me3_A'] + inputVal['KmetStat_H3K4me3_A'])
    nFI = (inputVal['KmetStat_H3K4me3_A'] + inputVal['KmetStat_H3K4me3_A']) / inputVal['tot']
    nFC = (chipVal['KmetStat_H3K4me3_A'] + chipVal['KmetStat_H3K4me3_A']) / chipVal['tot']
    nF = nFI / nFC

    tmpBGInit = tmpDir + '/init.bg'
    tmpBGCorr = tmpDir + '/corrected.bg'
    tmpBGKM = BGKmet

    with open(tmpBGInit, 'w') as BGOUT, open(tmpBGKM, 'w') as BGKMET:
        HMDcs = {}
        with subprocess.Popen(['paste', chipBG, inputBG], stdout=subprocess.PIPE) as BGPIPE:
            for line in BGPIPE.stdout:
                line = line.decode('utf-8').rstrip()
                cs, frm, to, cVal, mt1, mt2, mt3, iVal = line.split('\t')

                if int(iVal) == 0 or int(cVal) == 0:
                    HMD = 0
                else:
                    HMD = round(100 * ((int(cVal) / int(iVal)) / IPe_H3K4me3), 3)
                    if cs in HMDcs:
                        HMDcs[cs]['val'] += HMD
                        HMDcs[cs]['cnt'] += 1
                    else:
                        HMDcs[cs] = {'val': HMD, 'cnt': 1}

                if re.search(r'Kmet', cs):
                    BGKMET.write('\t'.join([cs, frm, to, str(HMD)]) + '\n')
                else:
                    BGOUT.write('\t'.join([cs, frm, to, str(HMD)]) + '\n')

    normByCSmean(tmpBGInit, HMDcs, BGout)

    #########################################################################
    def getKmetStatVals(b, ret):
        with open(b) as IN:
            for line in IN:
                line = line.rstrip()
                F = line.split('\t')
                ret[F[0]] = F[1]
                if re.search(r'Kmet', line):
                    ret['totKmet'] += int(F[1])
                ret['tot'] = int(F[2])
    #########################################################################
    def normByCSmean(inBG, csMean, outBG):
      with open(inBG, 'r') as INBG, open(outBG, 'w') as OUTBG:
          for line in INBG:
              line = line.strip()
              F = line.split('\t')
              if F[0].startswith('chr'):
                  normHMD = csMean[F[0]]['val'] / csMean[F[0]]['cnt']
                  OUTBG.write('\t'.join(F[:3]) + '\t' + str(F[3]-normHMD if F[3] else 0) + '\n')
              else:
                  OUTBG.write('\t'.join(F[:4]) + '\n')
        


### Original Perl code
use strict; 

use Getopt::Long;
use File::Basename; 
use File::Find;
use File::Which;

GetOptions ('t=s'     		=> \(my $chipBG),
	    	'c=s'      		=> \(my $inputBG),
	    	'sT=s'  		=> \(my $chipBGstat),
	    	'sC=s'  		=> \(my $inputBGstat),
	    	'out=s'    		=> \(my $BGout),
	    	'faidx=s' 		=> \(my $fai),
	    	'tmpDir=s' 		=> \(my $tmpDir = '.'),
	    	's=s'      		=> \(my $scoretype = 'KmetStat'));

die unless (-e $chipBG);
die unless (-e $inputBG);
die if     (-e $BGout);
die unless (-e $fai);
die ('bedGraphToBigWig from UCSC not available!!') unless (which('bedGraphToBigWig'));

my $BGKmet = $BGout; $BGKmet =~ s/(bg|bedgraph)$//; $BGKmet = $BGKmet.'KMETChromosomes.bedgraph';

die unless (-e $chipBGstat);
die unless (-e $inputBGstat);

my (%chipVal,%inputVal);

getKmetStatVals($chipBGstat,\%chipVal);
getKmetStatVals($inputBGstat,\%inputVal);

my $H4me3A_IPe = $chipVal{'KmetStat_H3K4me3_A'} / $inputVal{'KmetStat_H3K4me3_A'};
my $H4me3B_IPe = $chipVal{'KmetStat_H3K4me3_B'} / $inputVal{'KmetStat_H3K4me3_B'};
my $H4me2A_IPe = $chipVal{'KmetStat_H3K4me2_A'} / $inputVal{'KmetStat_H3K4me2_A'};
my $H4me2B_IPe = $chipVal{'KmetStat_H3K4me2_B'} / $inputVal{'KmetStat_H3K4me2_B'};

my $IPe_H3K4me3 = ($chipVal{'KmetStat_H3K4me3_A'} + $chipVal{'KmetStat_H3K4me3_B'}) / ($inputVal{'KmetStat_H3K4me3_A'} + $inputVal{'KmetStat_H3K4me3_A'});
my $nFI         = ($inputVal{'KmetStat_H3K4me3_A'} + $inputVal{'KmetStat_H3K4me3_A'})/($inputVal{'tot'});
my $nFC         = ($chipVal{'KmetStat_H3K4me3_A'} + $chipVal{'KmetStat_H3K4me3_A'})/($chipVal{'tot'});
my $nF          = $nFI/$nFC;

my $tmpBGInit   = $tmpDir.'/init.bg';
my $tmpBGCorr   = $tmpDir.'/corrected.bg';

my $tmpBGKM     = $BGKmet;

open BGOUT,  '>', $tmpBGInit;
open BGKMET, '>', $tmpBGKM;

my %HMDcs;

open my $BGPIPE, '-|', "paste $chipBG $inputBG";
#open BGChIP, $chipBG;
#open BGInput, $inputBG;

#while (defined(my $nChIP=<BGChIP>) && defined(my $nInput=<BGInput>)){
while (<$BGPIPE>){
	chomp ;
	my ($cs,$from,$to,$cVal,$mt1,$mt2,$mt3,$iVal) = split(/\t/,$_);

	my $HMD;

	if ($iVal == 0 || $cVal == 0){
		$HMD=0;
	}else{
		$HMD = sprintf("%4.3f",100 * (($cVal/$iVal)/$IPe_H3K4me3));

		$HMDcs{$cs}->{val} += $HMD;
		$HMDcs{$cs}->{cnt} ++;
	}

	if ($cs =~ /Kmet/){
		print BGKMET join("\t",$cs,$from,$to,$HMD)."\n";
	}else{
		print BGOUT  join("\t",$cs,$from,$to,$HMD)."\n" ;
	}
}

close BGOUT; 
close BGKMET;

normByCSmean($tmpBGInit,\%HMDcs,$BGout);

#########################################################################
sub getKmetStatVals{
	my ($b,$ret) = @_;
	open IN, $b;
	while (<IN>){
		chomp;
		my @F = split(/\t/,$_);
		$$ret{$F[0]}      = $F[1];
		$$ret{'totKmet'} += $F[1] if ($_ =~ /Kmet/);
		$$ret{'tot'}      = $F[2];
	}
	close IN;
}

#########################################################################
sub normByCSmean{
	my ($inBG,$csMean,$outBG) = @_;
	open INBG, $inBG;
	open OUTBG, '>' ,$outBG;
	while (<INBG>){
		chomp;
		my @F = split(/\t/,$_);
		if ($F[0] =~ /^chr\S+$/){
			my $normHMD = $$csMean{$F[0]}->{val}/$$csMean{$F[0]}->{cnt};
			print OUTBG join("\t",@F[0..2],($F[3]?($F[3]-$normHMD):0)."\n");
		}else{
			print OUTBG join("\t",@F[0..3]."\n");
		}
	}
	close INBG;close OUTBG;
}