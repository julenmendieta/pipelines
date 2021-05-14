
/scratch/julen/chiadrop/CIMAdata/CRUK_chiadrop2021

bases=16;

fip=/scratch/julen/chiadrop/CIMAdata/CRUK_chiadrop2021;
fi=SLX-20379.SIGAH9.HYHJ2DRXX.s_2.r_1.fq.gz;

infip=${fip}"/"${fi};
# get first part of name
outfi=`echo $fi | sed "s/\.s.*//g"`;
# define following file names
outfull=${outfi}"_"${bases}"bp_full.txt";
outfilt=${outfi}"_"${bases}"bp_noDup.txt";

zcat ${infip} | 
    # having introduced $bases varaible inside, set field
    #separator to blank space (each character is a field)
    # when NR (line number) divded by 4 has a remainder of 2 
    #we are in second line, were sequence is
    # We iterate over values of i less than or equal to
    #$bases
    awk -v len="$bases" -F "" '{
        if (NR % 4 == 2) 
            {for (i=1; i<=len; i++) {
                    printf $(i)
                }; printf "\n"
            }
        }' > ${outfull};

    # we now filter duplicated IDs
    sort ${outfull} | uniq > ${outfilt}
