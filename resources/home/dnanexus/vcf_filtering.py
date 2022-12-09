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


def read_vcf_in(args, SAMPLE):
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

    # g.vcf has noref info in it, the first is usually the alt variant
    vcf_df['ALT'] = vcf_df['ALT'].str.split(',', expand=True)[0]

    return vcf_df

def main():
    args = parse_args()

    # get the sample name from the VCF
    SAMPLE = args.vcf.split("_")[0]

    # Read in VCF
    vcf_df = read_vcf_in(args, SAMPLE)
    PRS = pd.read_csv(args.PRS_variant_file, sep=",", skiprows=1,
                        compression='infer')

    # CANRISK is sassy and needs all sites. To deal with it, we will
    # check if position and REF>ALT dont match, then we can add this is
    # but genotype 0/0
    match = 0
    no_match = 0

    new_vcf_df = pd.DataFrame()

    for variant in range(0, len(PRS.index)):
        chr = PRS.loc[variant,"Chromosome"]
        pos = PRS.loc[variant,"Position"]
        ref = PRS.loc[variant,"Reference_Allele"]
        alt = PRS.loc[variant,"Effect_Allele"]
        vcf_df_subset = vcf_df[(vcf_df['CHROM'] == chr) & (vcf_df['POS'] == pos) & (vcf_df['REF'] == ref) & (vcf_df['ALT'] == alt)]
        # if the chr,pos,ref and alt from PRS doesnt exist in our VCF,
        # we need to append this to a new df with details from PRS with
        # GT = "./". If there is data, then append this to the empty df
        if vcf_df_subset.empty:
            no_match = no_match +1
            new_row = pd.Series({'CHROM': chr,
                    'POS': pos,
                    'ID': '.',
                    'REF': ref,
                    'ALT': alt,
                    'QUAL': '.',
                    'FILTER': '.',
                    'INFO': '.',
                    'FORMAT': 'GT',
                    SAMPLE:'./.'
                    })
            #append row to the dataframe
            # if you concate a dict it makes it into a single col so transpose
            new_vcf_df = pd.concat([new_vcf_df, new_row.to_frame().T], ignore_index=True)
        else:
            #if there is a match then add this to the new_df
            match = match + 1
            new_vcf_df = pd.concat([new_vcf_df, vcf_df_subset], ignore_index=True)

    print(f"Number of variants in VCF as expected : {match}")
    print(f"Number of variants in VCF modified as it was not expected genotype : {no_match}")

    new_vcf_df = new_vcf_df.rename(columns={"CHROM": "#CHROM"})
    new_vcf_df = new_vcf_df[["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", SAMPLE]]
    print(new_vcf_df)
    new_vcf_df.to_csv(
            args.vcf.split(".")[0] + "_CANRISK.vcf",
            sep="\t", index=False, header=True
            )


if __name__ == "__main__":

    main()