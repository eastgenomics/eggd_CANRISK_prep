# eggd_CANRISK_prep

## What does this app do?

This app filters and formats VCFs to calculate PRS scores using the online CanRisk tool (https://www.canrisk.org/canrisk_tool/).

## What inputs are required for this app to run?

- gVCF (VCF also works but this is not accurate as homozygous reference variants are missed)
- PRS file from CanRisk (https://canrisk.atlassian.net/wiki/spaces/FAQS/pages/35979266/What+variants+are+used+in+the+PRS)
- Reference genome

## How does this app work?

The app filters the gVCF to the PRS variants. Then it normalises multiallelic variants. This intermediate file is run through the custom python script to reformat the VCF and add any missing variants with null genotype "./."


## What does this app output?

A single VCF with the PRS name appended to the VCF name.

## What limitations does this app have?

This has been tested for germline gVCFs outputted from Sentieon DNAnexus app.

## This app was created by East Genomics GLH
