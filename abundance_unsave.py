import argparse
import math
import pandas as pd
import os


def main():
    # parse arguments
    parser = argparse.ArgumentParser(description='Process data files.')
    parser.add_argument('-f', required=True, help='Path to tsv file')
    parser.add_argument('-c', required=True, help='Path to count file')
    parser.add_argument('-o', required=True, help='Path to output directory')

    args = parser.parse_args()

    # get input file name
    input_file_name = args.f.split('/')[-1].split('.')[0]

    print(f"Processing {input_file_name} file")

    # read files & merge
    tsv_file = pd.read_csv(args.f, sep='\t', low_memory=False)
    count_file = pd.read_csv(args.c, sep='\t')
    merged_df = pd.merge(tsv_file, count_file, left_on='qseqid', right_on='Name', how='inner')

    # count total reads
    total_reads = merged_df.NumReads.sum()

    # col of kindom, family, genus, NumReads
    selected_columns = ['family', 'genus', 'NumReads', 'phylum', 'class', 'order']

    # viruses dataframe
    vi = merged_df[merged_df.kindom == 'k__Viruses'][selected_columns]

    # rename NumReads column
    vi.rename(columns={'NumReads': f'{input_file_name}:reads'}, inplace=True)

    # Family data processing
    f_res = vi.groupby(['family', 'phylum', 'class', 'order'])[f'{input_file_name}:reads'].sum().reset_index()

    # calculate rpm & logRPM
    f_res[f'{input_file_name}:rpm'] = f_res.iloc[:, 4].apply(lambda x: x * 1e6 / total_reads)
    f_res[f'{input_file_name}:logRPM'] = f_res[f_res.iloc[:, 5] > 1].iloc[:, 5].apply(math.log10)

    # select rows that logRPM is not NaN
    f_final = f_res[~f_res[f'{input_file_name}:logRPM'].isna()].iloc[:, list(range(4)) + [6]]

    # construct output path
    family_output_file = os.path.join(args.o, f'F{input_file_name}.csv')

    # write to file
    f_final.to_csv(family_output_file, index=False)

    # Genus data processing
    g_res = vi.groupby(['genus', 'phylum', 'class', 'order', 'family'])[f'{input_file_name}:reads'].sum().reset_index()

    # calculate rpm & logRPM
    g_res[f'{input_file_name}:rpm'] = g_res.iloc[:, 5].apply(lambda x: x * 1e6 / total_reads)
    g_res[f'{input_file_name}:logRPM'] = g_res[g_res.iloc[:, 6] > 1].iloc[:, 6].apply(math.log10)

    # select rows that logRPM is not NaN
    g_final = g_res[~g_res[f'{input_file_name}:logRPM'].isna()].iloc[:, list(range(5)) + [7]]

    # construct output path
    genus_output_file = os.path.join(args.o, f'G{input_file_name}.csv')

    # write to file
    g_final.to_csv(genus_output_file, index=False)


if __name__ == "__main__":
    main()
