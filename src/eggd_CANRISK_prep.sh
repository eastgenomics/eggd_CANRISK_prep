#!/bin/bash

set -exo pipefail

main() {

    # Download inputs from DNAnexus in parallel, these will be downloaded to /home/dnanexus/in/
    dx-download-all-inputs --parallel

    # Install packages from the python asset
    pip3 install /pytz-*.whl /numpy-*.whl /pandas-*.whl

    bedtools intersect -a $vcf_path -b $PRS_variants_path -header > NA12878-NA12878-1-CEN-F-EGG5_markdup_recalibrated_Haplotyper_filteredPOS.vcf

    filtered_VCF=$(dx upload filtered_VCF --brief)
    dx-jobutil-add-output filtered_VCF "$filtered_VCF" --class=file
}
