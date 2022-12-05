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

    # make the output filename a combo of sample VCF and PRS cancer type
    # to track the sample and filtering applied
    vcf_filename=${vcf_path##*/}
    PRS_variants_filename=${PRS_variants_path##*/}
    out_filename=${vcf_filename%%.*}_${PRS_variants_filename%.*}

    # intersect the VCF and PRS variant list to save mem in python
    bedtools intersect -a $vcf_path -b  PRS_variants.bed -header > ${out_filename}.vcf

    # run the VCF and PRS variant list in python to get the same
    # variants in the VCF to upload to CANRISK
    python3 vcf_filtering.py -v ${out_filename}.vcf -p $PRS_variants_path
    out_file=$(ls | grep "CANRISK.vcf")

    # python didnt output the header, so lets add that from the original
    # VCF and output it to the same filename
    cat <(zcat $vcf_path | grep ^"#") <(grep -v ^"#" $out_file | sort -k1,1V -k2,2n ) > $(echo $out_file)

    # upload VCF
    filtered_VCF=$(dx upload $out_file --brief)
    dx-jobutil-add-output filtered_VCF "$filtered_VCF" --class=file
}
