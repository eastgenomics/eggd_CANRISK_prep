#!/bin/bash

set -exo pipefail

main() {

    # Download inputs from DNAnexus in parallel, these will be downloaded to /home/dnanexus/in/
    dx-download-all-inputs --parallel

    # Install packages from the python asset
    pip3 install /pytz-*.whl /numpy-*.whl /pandas-*.whl

    # make the bed file chr, start and end, then remove the first rows
    # which is the alpha value for the PRS risk and headers
    awk -F "," '{print $1"\t"$2-1"\t"$2}' $PRS_variants_path | tail -n +3  > PRS_variants.bed

    # intersect the VCF and PRS variant list to save mem in python
    vcf_filename=${vcf_path##*/}
    PRS_variants_filename=${PRS_variants_path##*/}

    out_filename=${vcf_filename%%.*}_${PRS_variants_filename%.*}
    echo $out_filename
    bedtools intersect -a $vcf_path -b  PRS_variants.bed -header > ${out_filename}.vcf

    head ${out_filename}.vcf

    # filtered_VCF=$(dx upload filtered_VCF --brief)
    # dx-jobutil-add-output filtered_VCF "$filtered_VCF" --class=file
}
