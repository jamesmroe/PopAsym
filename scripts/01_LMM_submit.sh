#!/bin/bash

#calls matlab func PopAsym_LMM()
cohorts=" \
DLBS_55minus \
SALD_55minus \
IXI_55minus" #\ 
#INDIVIDUAL-LEVEL DATA FOR THESE SPECIFIC COHORTS IS AVAILABLE WITHOUT RESTRICTIONS
#AND IS THUS PROVIDED HERE

#LCBC_55minus \
#CamCam_55minus \
#HCP_55minus \
#UKB_55minus"


for f in area thickness; do
		
	for cohort in $cohorts; do
		
		#dirs
		wdir=$(echo `dirname \`pwd\``) #portable projectdir
		adir=$wdir/data/$cohort
		sdir=$wdir/subjects
		odir=$wdir/results/$f/$cohort
		if [ ! -e $odir ]; then mkdir -p $odir; fi

		#files
		Y="${adir}/MIXEDHEMI.lh.${f}.CONCAT.40tvs8ref_sym_20.sm8.mgh"
		Qdec="Qdec_${cohort}.dat"
		surf="lh.sphere"
		subj="40tvs8ref_sym_20"
		mask="lh.cortex.label"


echo -e "
### FOR TESTING IN MATLAB ###
wdir='$wdir'
adir='$adir'
sdir='$sdir'
odir='$odir' 
Y='$Y'
Qdec='$Qdec'
subj='$subj'
surf='$surf'
mask='$mask'
cohort='$cohort'
f='$f'

"

		echo "Submitting LME: ${cohort} ${f}"
		#sbatch PopAsym_LMM.sh $wdir $adir $sdir $odir $Y $Qdec $subj $surf $mask $cohort $f
		sleep 1;
	done
done
