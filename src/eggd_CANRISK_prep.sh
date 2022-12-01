#!/bin/bash

set -exo pipefail

main() {

    dx download "$VCF" -o VCF

    dx download "$PRS_variants" -o PRS_variants

    filtered_VCF=$(dx upload filtered_VCF --brief)

    dx-jobutil-add-output filtered_VCF "$filtered_VCF" --class=file
}
