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


    return vcf_df

def main():
    args = parse_args()

    # get the sample name from the VCF
    SAMPLE = args.vcf.split("_")[0]

    # Read in VCF
    vcf_df = read_vcf_in(args, SAMPLE)
    # Read in the PRS varaints file
    PRS = pd.read_csv(args.PRS_variant_file, sep=",", skiprows=1,
                        compression='infer')


    # Read in the BCAC_313 PRS
    PRS = pd.read_csv(args.PRS_variant_file, sep=",", skiprows=1,
                        compression='infer')

    # Now loop through our VCF, get the chr & pos of each variant and check
    # that the REF>ALT matches whats in PRS variants file. If there is
    # a match keep, else delete row
    variants_to_drop = []
    for variant in range(0, len(vcf_df.index)):
        chr = vcf_df.loc[variant,"CHROM"]
        pos = vcf_df.loc[variant,"POS"]
        ref = vcf_df.loc[variant,"REF"]
        alt = vcf_df.loc[variant,"ALT"]
        PRS_subset = PRS[(PRS['Chromosome'] == chr) & (PRS['Position'] == pos)]
        if PRS_subset.empty:
            next
        else:
            if PRS_subset['Reference_Allele'].values[0] == ref and PRS_subset['Effect_Allele'].values[0] == alt:
                next
            else:
                print("Incorrect Ref>Alt change, removing this variant")
                variants_to_drop.append(variant)

    # remove variants that have unexpected Ref>Alt change
    vcf_df = vcf_df.drop(variants_to_drop)
    # CANRISK is sassy and needs all sites. To deal with it, we will
    # check if position and REF>ALT dont match, then we can add this is
    # but genotype 0/0

    print(vcf_df)

    new_vcf_df = pd.DataFrame()

    for variant in range(0, len(PRS.index)):
        chr = PRS.loc[variant,"Chromosome"]
        pos = PRS.loc[variant,"Position"]
        ref = PRS.loc[variant,"Reference_Allele"]
        alt = PRS.loc[variant,"Effect_Allele"]
        vcf_df2_subset = vcf_df[(vcf_df['CHROM'] == chr) & (vcf_df['POS'] == pos) & (vcf_df['REF'] == ref) & (vcf_df['ALT'] == alt)]
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
    # let's combine vcf_df and new_vcf_df
    final_vcf =  pd.concat([vcf_df, new_vcf_df], ignore_index=True)
    final_vcf['CHROM'] = 'chr' + final_vcf['CHROM'].astype(str)
    final_vcf = final_vcf.rename(columns={"CHROM": "#CHROM"})
    print(final_vcf)
    final_vcf.to_csv(
            args.vcf.split(".")[0] + "_CANRISK.vcf",
            sep="\t", index=False, header=True
            )

if __name__ == "__main__":

    main()