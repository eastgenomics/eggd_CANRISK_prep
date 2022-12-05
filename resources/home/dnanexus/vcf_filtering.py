import argparse
from re import A, M
from sre_constants import BRANCH
import pandas as pd

def parse_args():
    """Parse through arguments
    Returns:
        args: Variable that you can extract relevant
        arguments inputs needed
    """
    # Read in arguments
    parser = argparse.ArgumentParser()

    parser.add_argument(
        '-v', '--vcf',
        help='VCF filtered for the PRS regions',
        required=True
        )

    parser.add_argument(
        '-p', '--PRS_variant_file',
        help='PRS variant file',
        required=True
        )

    args = parser.parse_args()

    return args

def main():
    args = parse_args()
    # get the sample name from the VCF
    SAMPLE = args.vcf.split("_")[0]

    # read in the filtered VCF
    cols = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", SAMPLE]
    # read vcf records into df
    vcf_df = pd.read_csv(args.vcf, sep="\t", comment='#',
                        names=cols, compression='infer')
    # From the FORMAT column, we are only interested in the GT info
    vcf_df['FORMAT'] = vcf_df['FORMAT'].str.split(':', expand=True)[0]
    # From the sample column, we are only interested in the GT info
    vcf_df[SAMPLE] = vcf_df[SAMPLE].str.split(':', expand=True)[0]

    # Remove some data from columns
    vcf_df['ID'] = "."
    vcf_df['QUAL'] = "."
    vcf_df['INFO'] = "."

    # add empty column to state whether we have the right genotype
    vcf_df["right_Ref2Alt"] = ""
    print(vcf_df)

    # Read in the BCAC_313 PRS and filter to match REF>ALT change
    BCAC_PRS = pd.read_csv(args.PRS_variant_file, sep=",", skiprows=1,
                        compression='infer')


    # Now loop through our VCF, get the chr & pos of each variant and check
    # that the REF>ALT matches whats in BCAC_313 PRS

    for variant in range(0, len(vcf_df.index)):
        chr = vcf_df.loc[variant,"CHROM"]
        pos = vcf_df.loc[variant,"POS"]
        ref = vcf_df.loc[variant,"REF"]
        alt = vcf_df.loc[variant,"ALT"]
        BCAC_PRS_subset = BCAC_PRS[(BCAC_PRS['Chromosome'] == chr) & (BCAC_PRS['Position'] == pos)]
        if BCAC_PRS_subset.empty:
            # print("empty dataframe")
            next
        else:
            if BCAC_PRS_subset['Reference_Allele'].values[0] == ref and BCAC_PRS_subset['Effect_Allele'].values[0] == alt:
                vcf_df.loc[variant, "right_Ref2Alt"] = "Yes"
            else:
                vcf_df.loc[variant, "right_Ref2Alt"] = "No"

    print(vcf_df)
    print(vcf_df['right_Ref2Alt'].value_counts())

    vcf_df2 = vcf_df.loc[vcf_df['right_Ref2Alt'] == 'Yes']
    vcf_df2 = vcf_df2[["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", SAMPLE]]

    # CANRISK is sassy and needs all sites. To deal with it, we will
    # check if position and REF>ALT dont match, then we can add this is
    # but genotype 0/0

    print(vcf_df2)

    new_vcf_df = pd.DataFrame()

    for variant in range(0, len(BCAC_PRS.index)):
        chr = BCAC_PRS.loc[variant,"Chromosome"]
        pos = BCAC_PRS.loc[variant,"Position"]
        ref = BCAC_PRS.loc[variant,"Reference_Allele"]
        alt = BCAC_PRS.loc[variant,"Effect_Allele"]
        vcf_df2_subset = vcf_df2[(vcf_df2['CHROM'] == chr) & (vcf_df2['POS'] == pos) & (vcf_df2['REF'] == ref) & (vcf_df2['ALT'] == alt)]
        # if the chr,pos,ref and alt doesnt exist in our VCF, we should include it,
        if vcf_df2_subset.empty:
            new_row = pd.Series({'CHROM': chr,
                    'POS': pos,
                    'ID': '.',
                    'REF': ref,
                    'ALT': alt,
                    'QUAL': '.',
                    'FILTER': '.',
                    'INFO': '.',
                    'FORMAT': 'GT',
                    SAMPLE:'0/0'
                    })
            #append row to the dataframe
            # if you concate a dict it makes it into a single col so transpose
            new_vcf_df = pd.concat([new_vcf_df, new_row.to_frame().T], ignore_index=True)



    # this contains all the GT's that don't exist/match what CANRISK wants
    # let's combine vcf_df3 and new_vcf_df
    final_vcf =  pd.concat([vcf_df2, new_vcf_df], ignore_index=True)
    final_vcf['CHROM'] = 'chr' + final_vcf['CHROM'].astype(str)
    final_vcf = final_vcf.rename(columns={"CHROM": "#CHROM"})
    print(final_vcf)
    final_vcf.to_csv(
            args.vcf.split(".")[0] + "_CANRISK.vcf",
            sep="\t", index=False, header=True
            )
    # we are missing headers so lets add that back in
    #command = "grep ^'##' NA12878-NA12878-1-CEN-F-EGG5_markdup_recalibrated_Haplotyper.vcf > header.txt"
    #command2 = "cat header.txt NA12878-NA12878-1-CEN-F-EGG5_markdup_recalibrated_Haplotyper_filteredPOS_REF_ALT.vcf > NA12878-NA12878-1-CEN-F-EGG5_markdup_recalibrated_Haplotyper_filteredPOS_REF_ALT_final.vcf"
    #subprocess.check_output(command, stderr=subprocess.STDOUT, shell = True)
    #subprocess.check_output(command2, stderr=subprocess.STDOUT, shell = True)

if __name__ == "__main__":

    main()