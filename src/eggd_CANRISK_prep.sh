#!/bin/bash

set -exo pipefail

main() {

    # Download inputs from DNAnexus in parallel, these will be downloaded to /home/dnanexus/in/
    dx-download-all-inputs --parallel

    # Install packages from the python asset
    pip3 install /pytz-*.whl /numpy-*.whl /pandas-*.whl

    # samtools in htslib doesn't work as its missing a library, so
    # will install the missing libraries from the downloaded deb files
    # (so as to not use internet)
    sudo dpkg -i libtinfo5_6.2-0ubuntu2_amd64.deb
    sudo dpkg -i libncurses5_6.2-0ubuntu2_amd64.deb

    # make the bed file chr, start and end, then remove the first rows
    # which is the alpha value for the PRS risk and headers
    awk -F "," '{print $1"\t"$2-1"\t"$2"\t"$3"\t"$4}' $PRS_variants_path | tail -n +3  > PRS_variants.bed

    # make the output filename a combo of sample VCF and PRS cancer type
    # to track the sample and filtering applied
    vcf_filename=${vcf_path##*/}
    PRS_variants_filename=${PRS_variants_path##*/}
    out_filename=${vcf_filename%%.*}_${PRS_variants_filename%.*}

    # normalise variants

    # intersect the VCF and PRS variant list to save mem in python
    # keeps bed regions from the VCF that are missed but annotates with "."
    # first three columns is from -a which is the PRS bed file
    # bedtools intersect -b NA12878-NA12878-1-CEN-F-EGG5_markdup_recalibrated_Haplotyper.vcf.gz -a  PRS_variants.bed -loj
    # we cant intersect with log on gvcf, so filter to PRS sites, then do loj on it
    scp $reference_path .
    gunzip -d $reference_name
    samtools faidx ${reference_prefix}.fa

    bedtools intersect -a $vcf_path  -b PRS_variants.bed -header | bcftools norm -f ${reference_prefix}.fa -m -any --keep-sum AD - > ${out_filename}.vcf

    # run the VCF and PRS variant list in python to get the same
    # variants in the VCF to upload to CANRISK
    python3 vcf_filtering.py -v ${out_filename}.vcf -p $PRS_variants_path
    out_file=$(ls | grep "CANRISK.vcf")
    grep -v ^"#" $out_file | sort -k1,1V -k2,2n > ${out_file}_filtered.vcf
    rm ${out_filename}.vcf

    # python didnt output the header, so lets add that from the original
    # VCF and output it to the same filename
    cat <(zcat $vcf_path | grep ^"#") <( cat ${out_file}_filtered.vcf) > ${out_filename}.vcf

    # upload VCF
    filtered_VCF=$(dx upload $out_file --brief)
    dx-jobutil-add-output filtered_VCF "$filtered_VCF" --class=file
}