# input file
ifile=$1
# output file
ofile=$2
# beta file
bfile=$3


# extract id and weight
grep " W " $ifile | awk '{printf "%6d  %16.10e\n",NR-1,$NF}' > w_$ifile

# extract centers and convert to nm
grep " M " $ifile | awk '{printf "%13.7lf %13.7lf %13.7lf\n",$(NF-2)/10,$(NF-1)/10,$NF/10}' > m_$ifile

# extract covariance and convert to PLUMED format
grep " CovM " $ifile | awk '{printf "%13.10e %13.10e %13.10e ",$(NF-4)/100,$(NF-2)/100,$NF/100; if($(NF-1)=="zz")printf "\n"}' | awk '{print $1,$2,$3,$4,$5,$6}' > c_$ifile

# create header for PLUMED file
echo "#! FIELDS Id Weight Mean_0 Mean_1 Mean_2 Cov_00 Cov_01 Cov_02 Cov_11 Cov_12 Cov_22 Beta" > $ofile

# paste together and format
# in case beta file is provided
if [ ! -z $bfile  ]; then
 paste w_$ifile m_$ifile c_$ifile $bfile | awk '{printf "%6d  %16.10e  %13.7lf %13.7lf %13.7lf  ",$1,$2,$3,$4,$5; for(i=6;i<NF;i++) {printf "%13.10e ",$i}; printf " %d\n",$NF}' >> $ofile
else
# otherwise just one group of components (labelled as zero)
 paste w_$ifile m_$ifile c_$ifile        | awk '{printf "%6d  %16.10e  %13.7lf %13.7lf %13.7lf  ",$1,$2,$3,$4,$5; for(i=6;i<=NF;i++){printf "%13.10e ",$i}; printf " 0\n"}' >> $ofile
fi

# clear stuff
rm  w_$ifile m_$ifile c_$ifile
