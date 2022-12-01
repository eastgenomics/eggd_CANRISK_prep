#!/bin/bash

set -exo pipefail

main() {

    # Download inputs from DNAnexus in parallel, these will be downloaded to /home/dnanexus/in/
    dx-download-all-inputs --parallel

    # Install packages from the python asset
    pip3 install /pytz-*.whl /numpy-*.whl /pandas-*.whl

    

    filtered_VCF=$(dx upload filtered_VCF --brief)
    dx-jobutil-add-output filtered_VCF "$filtered_VCF" --class=file
}
