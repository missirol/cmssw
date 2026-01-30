#!/bin/bash

OUTDIR=tmp_plots

rm -rf "${OUTDIR}"

./fitJetResponse.py -k l1s_run3_jecFits \
  -i out2/harvesting/*root \
  -o "${OUTDIR}" \
  -m "*_pt_GENoverREC_Median_wrt_pt" \
  -e png \
  -l 'QCD-#hat{p}_{T}[15-7000], PU[0-120] (Run3Winter25)'
