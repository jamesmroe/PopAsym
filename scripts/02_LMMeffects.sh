#!/bin/bash

# ------------------------------------------------------------------
# Author: James M Roe
#         Description: compute overlapping effects across cohorts
# ------------------------------------------------------------------

#SCRIPT TAKES SUMMARY-LEVEL DATA AS INPUT AND IS THUS FULLY REPRODUCIBLE WITH THE MAPS PROVIDED HERE


echo "SOURCING FREESURFER"
export FREESURFER_HOME=/cluster/projects/p23/tools/mri/freesurfer/freesurfer.6.0.0
source $FREESURFER_HOME/SetUpFreeSurfer.sh
export SUBJECTS_DIR=$FREESURFER_HOME/subjects
echo "SOURCING FSL"
FSLDIR=/cluster/projects/p23/tools/mri/fsl/current
. ${FSLDIR}/etc/fslconf/fsl.sh
PATH=${FSLDIR}/bin:${PATH}
export FSLDIR PATH


cohorts="\
SALD_55minus \
DLBS_55minus \
IXI_55minus" #\
#LCBC_55minus \
#CamCan_55minus \
#HCP_55minus \
#UKB_55minus"

metric="\
area \
thickness"


wdir=$(echo `dirname \`pwd\``) #portable projectdir


for f in $metric; do
	odir=${wdir}/results/${f}/crossCohorts; if [ ! -d $odir ]; then mkdir $odir; fi
	for c in $cohorts; do

		d=${wdir}/results/${f}/${c}/Beta
		if [ ! -d $d ]; then 
			echo "> $d does not exist. Run LMM first"; exit 1
		else

			if [ ! -e "$odir/${f}.asymPercent.${c}.nii.gz" ]; then
				#compute Asym %
				#HemiBeta/(InterceptBeta+.5*HemiBeta)
				cd $d
				echo "> computing %asym for $c"
				sleep 1
				for i in Beta1 Beta3; do mri_convert $i.mgh $i.nii.gz; done
				fslmaths Beta3.nii.gz -mul 0.5 tmp1.nii.gz
				fslmaths Beta1.nii.gz -add tmp1.nii.gz tmp2.nii.gz
				fslmaths Beta3.nii.gz -div tmp2.nii.gz asymPercent.nii.gz
				rm tmp*

				echo "> copying to results/crossCohorts"
				sleep 1
				cp asymPercent.nii.gz $odir/${f}.asymPercent.${c}.nii.gz
			else
				echo "> $odir/${f}.asymPercent.${c}.nii.gz exists"
				sleep 1
			fi
		fi
	done
done


#COUNT ASYM MAPS
cohorts=" \
SALD_55minus \
DLBS_55minus \
IXI_55minus \
LCBC_55minus \
CamCan_55minus \
HCP_55minus \
UKB_55minus"
areaCount=0
cthCount=0
for f in $metric; do
	pdir=${wdir}/publishedMaps/${f}Asym
	for c in $cohorts; do
		if [ -e $pdir/${f}.asymPercent.${c}.nii.gz ]; then
			if [ $f == "area" ]; then
				areaCount=$((areaCount+1))
			elif [ $f == "thickness" ]; then
				cthCount=$((cthCount+1))
			fi
		fi
	done
done
echo -e "
> areaCount=$areaCount
> cthCount=$cthCount"
sleep 1


###COMPUTE MASKS FOR 6 DATASET OVERLAP###
#mask at 1% pos or 5% pos depending on thck or area
for f in $metric; do
	
	if [ $f == "area" ]; then
		thr=5 #effect threshold (%)
	elif [ $f == "thickness" ]; then
		thr=1 #effect threshold (%)
	fi
	odir=${wdir}/results/${f}/crossCohorts
	
	if [ ! -e ${odir}/sum*Neg.nii.gz ]; then	
		echo "> creating conjunction maps"
		sleep 1
		for c in $cohorts; do
			
			#take from pubMaps
			idir=$wdir/publishedMaps/${f}Asym

			echo "> ${thr}% mask ${c} ${f}"
			#left
			mri_binarize \
			--i $idir/${f}.asymPercent.${c}.nii.gz \
			--min 0.0${thr} \
			--o $odir/m${thr}Pos.${c}.nii.gz

			#right
			fslmaths $idir/${f}.asymPercent.${c}.nii.gz -mul -1 $odir/inv-${f}.asymPercent.${c}.nii.gz
			mri_binarize \
			--i $odir/inv-${f}.asymPercent.${c}.nii.gz \
			--min 0.0${thr} \
			--o $odir/m${thr}Neg.${c}.nii.gz

			rm $odir/inv-${f}.asymPercent.${c}.nii.gz
		done

		echo "> 6 dataset overlap"
		sleep 1
		mod=($cohorts)
		cd ${wdir}/results/${f}/crossCohorts
		if [ $f == "thickness" ]; then

		fslmaths m1Pos.${mod[0]}.nii.gz -add m1Pos.${mod[1]}.nii.gz -add m1Pos.${mod[2]}.nii.gz -add m1Pos.${mod[3]}.nii.gz -add m1Pos.${mod[4]}.nii.gz -add m1Pos.${mod[5]}.nii.gz -add m1Pos.${mod[6]}.nii.gz sum1Pos.nii.gz
		fslmaths m1Neg.${mod[0]}.nii.gz -add m1Neg.${mod[1]}.nii.gz -add m1Neg.${mod[2]}.nii.gz -add m1Neg.${mod[3]}.nii.gz -add m1Neg.${mod[4]}.nii.gz -add m1Neg.${mod[5]}.nii.gz -add m1Neg.${mod[6]}.nii.gz sum1Neg.nii.gz

		elif [ $f == "area" ]; then

		fslmaths m5Pos.${mod[0]}.nii.gz -add m5Pos.${mod[1]}.nii.gz -add m5Pos.${mod[2]}.nii.gz -add m5Pos.${mod[3]}.nii.gz -add m5Pos.${mod[4]}.nii.gz -add m5Pos.${mod[5]}.nii.gz -add m5Pos.${mod[6]}.nii.gz sum5Pos.nii.gz
		fslmaths m5Neg.${mod[0]}.nii.gz -add m5Neg.${mod[1]}.nii.gz -add m5Neg.${mod[2]}.nii.gz -add m5Neg.${mod[3]}.nii.gz -add m5Neg.${mod[4]}.nii.gz -add m5Neg.${mod[5]}.nii.gz -add m5Neg.${mod[6]}.nii.gz sum5Neg.nii.gz

		fi
		#remove tmp masks
		rm m*nii.gz 
	else
		echo "> conjunction maps exist"
	fi
done



###EXTRACT CLUSTERS WITH OVERLAPPING EFFECTS IN 6/7 COHORTS###
overlap=6
for f in $metric; do

	odir=${wdir}/results/${f}/crossCohorts

	if [ ! -e ${odir}/cluster.$f.6N*merge*.label ]; then

		echo "> extracting clusters"
		sleep 1
		clustdir=$odir/surfclust${overlap}
		if [ ! -d $clustdir ]; then mkdir -p $clustdir; fi
		#for left (pos) and right (neg)
		for i in `seq 2`; do

			if [ $i == 1 ]; then
				if [ $f == "thickness" ]; then
					inp=sum1Pos.nii.gz
				elif [ $f == "area" ]; then
					inp=sum5Pos.nii.gz
				fi
				direction=p
			else
				if [ $f == "thickness" ]; then
					inp=sum1Neg.nii.gz
				elif [ $f == "area" ]; then
					inp=sum5Neg.nii.gz
				fi
				direction=n
			fi
				

			mri_surfcluster \
			--in $odir/$inp \
			--subject 40tvs8ref_sym_20 \
			--hemi lh \
			--surf inflated \
			--minarea 200 \
			--annot aparc.DKTatlas40 \
			--sum $clustdir/summary${direction}.txt \
			--olab $clustdir/cluster.$f.${overlap}${direction} \
			--thmin $overlap

			#threshold summed map at 3 for visualization
			mri_surfcluster \
			--in $odir/$inp \
			--subject 40tvs8ref_sym_20 \
			--hemi lh \
			--surf inflated \
			--o $odir/thr3.$inp \
			--thmin 3


			cd $clustdir

			#specific based on summary file (clusters > 200m2)
			if [ $f == "thickness" ] && [ $i == 1 ]; then			
				mri_mergelabels \
				-i cluster.thickness.6p-0001.label \
				-i cluster.thickness.6p-0002.label \
				-i cluster.thickness.6p-0003.label \
				-i cluster.thickness.6p-0004.label \
				-i cluster.thickness.6p-0005.label \
				-i cluster.thickness.6p-0006.label \
				-i cluster.thickness.6p-0007.label \
				-i cluster.thickness.6p-0008.label \
				-i cluster.thickness.6p-0009.label \
				-i cluster.thickness.6p-0010.label \
				-i cluster.thickness.6p-0011.label \
				-o ../cluster.thickness.6P-merge111.label
			elif [ $f == "thickness" ] && [ $i == 2 ]; then
				mri_mergelabels \
				-i cluster.thickness.6n-0001.label \
				-i cluster.thickness.6n-0002.label \
				-i cluster.thickness.6n-0003.label \
				-i cluster.thickness.6n-0004.label \
				-i cluster.thickness.6n-0005.label \
				-i cluster.thickness.6n-0006.label \
				-i cluster.thickness.6n-0007.label \
				-i cluster.thickness.6n-0008.label \
				-i cluster.thickness.6n-0009.label \
				-o ../cluster.thickness.6N-merge19.label
			elif [ $f == "area" ] && [ $i == 1 ]; then
				mri_mergelabels \
				-i cluster.area.6p-0001.label \
				-i cluster.area.6p-0002.label \
				-i cluster.area.6p-0003.label \
				-i cluster.area.6p-0004.label \
				-i cluster.area.6p-0005.label \
				-i cluster.area.6p-0006.label \
				-i cluster.area.6p-0007.label \
				-o ../cluster.area.6P-merge17.label
			elif [ $f == "area" ] && [ $i == 2 ]; then
			 	mri_mergelabels \
				-i cluster.area.6n-0001.label \
				-i cluster.area.6n-0002.label \
				-i cluster.area.6n-0003.label \
				-i cluster.area.6n-0004.label \
				-i cluster.area.6n-0005.label \
				-i cluster.area.6n-0006.label \
				-i cluster.area.6n-0007.label \
				-o ../cluster.area.6N-merge17.label
			fi

		done
	else
		echo "> clusters already extracted"
	fi
done

