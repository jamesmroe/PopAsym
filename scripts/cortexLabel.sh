#/bin/bash

#1 - 40tvs8ref_sym_20 --> $SUBJECTS_DIR (https://www.dropbox.com/s/cdjlm119jypye7b/LH_Sym.zip?dl=0)

#2 - make lh.cortex.label
for h in lh rh; do
	if [ ! -e "$SUBJECTS_DIR/40tvs8ref_sym_20/label/${h}.cortex.label" ]; then
		echo "> making ${h}.cortex.label"
		
		mkdir -p $SUBJECTS_DIR/40tvs8ref_sym_20/label/aparc.DKTatlas40/$h
		mri_annotation2label \
		--annotation aparc.DKTatlas40 \
		--subject 40tvs8ref_sym_20 \
		--hemi ${h} \
		--surf inflated \
		--sd $SUBJECTS_DIR \
		--outdir $SUBJECTS_DIR/40tvs8ref_sym_20/label/aparc.DKTatlas40/${h}

		mri_mergelabels \
		-d $SUBJECTS_DIR/40tvs8ref_sym_20/label/aparc.DKTatlas40/${h} \
		-o $SUBJECTS_DIR/40tvs8ref_sym_20/label/${h}.cortex.label
		
	fi
done

