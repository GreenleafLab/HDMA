#!/usr/bin/bash

# Purpose: download detailed elements BED files from the ENCODE V3 registry
# of candidate cis-regulatory elements (cCREs).
# via: https://screen.encodeproject.org/index/cversions

mkdir -p output/04-enhancers/07/ENCODE
cd output/04-enhancers/07/ENCODE

## ENCODE V3
# wget https://downloads.wenglab.org/Registry-V3/GRCh38-cCREs.bed
# wget https://downloads.wenglab.org/Registry-V3/GRCh38-cCREs.PLS.bed
# wget https://downloads.wenglab.org/Registry-V3/GRCh38-cCREs.pELS.bed
# wget https://downloads.wenglab.org/Registry-V3/GRCh38-cCREs.dELS.bed
# wget https://downloads.wenglab.org/Registry-V3/GRCh38-cCREs.CTCF-only.bed
# wget https://downloads.wenglab.org/Registry-V3/GRCh38-cCREs.DNase-H3K4me3.bed

# ENCODE V4
wget https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.bed
wget https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.PLS.bed
wget https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.pELS.bed
wget https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.dELS.bed
wget https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.CA-CTCF.bed
wget https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.CA-H3K4me3.bed
wget https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.CA-TF.bed
wget https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.CA.bed
wget https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.TF.bed
wget https://downloads.wenglab.org/GRCh38-cCREs.CTCF-bound.bed