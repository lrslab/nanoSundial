import pandas as pd
import numpy as np
import argparse

def analyze_data(input_file, output_file, coverage_cutoff, mean_cutoff, dwell_cutoff, pvalue_cutoff, shift_size):
    df = pd.read_csv(input_file)
    df = df[(df['sample_coverage'] >= coverage_cutoff) & (df['control_coverage'] >= coverage_cutoff)]

    df['-log10(fdr)'] = -np.log10(df['adj_p'])
    df = df[(df['mean_differ'].abs() >= mean_cutoff) & (df['dwell_differ'].abs() >= dwell_cutoff)]
    df = df[(df['-log10(fdr)'] >= pvalue_cutoff)]
    df['start'] = df['position']
    df['end'] = df['position'] + 1

    Start = None
    End = None
    Strand = None
    result_list = []
    grouped_df = df.groupby(['chrom', 'strand'])

    for group_id, temp_df in grouped_df:
        chrome, strand = group_id
        temp_df.sort_values(by='position', inplace=True)
        temp_df['differ'] = temp_df['position'].diff().fillna(-1).astype(int)

        for idx, row in temp_df.iterrows():
            if row['differ'] > shift_size or row['differ'] == -1:
                if Start is not None:
                    result_list.append([chrome, Start - shift_size, End + shift_size, '.', '.', Strand])
                Start = row['start']
                End = row['end']
                Strand = row['strand']
            elif row['differ'] == 0:
                raise RuntimeError("table error")
            else:
                End = row['end']

        result_list.append([chrome, Start - shift_size, End + shift_size, '.', '.', Strand])

    result_df = pd.DataFrame(result_list)
    result_df.to_csv(output_file, sep='\t', index=False, header=False)
    print("Merging finished! Results written to {}".format(output_file))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merging the raw result of nanoSundial")
    parser.add_argument("--input", required=True, help="Input CSV file")
    parser.add_argument("--output", default='nanoSundial_merged_region.bed', help="Output BED file")
    parser.add_argument("--coverage_cutoff", type=int, default=50, help="Coverage cutoff value")
    parser.add_argument("--mean_cutoff", type=float, default=0.18, help="cutoff of absolute mean difference ")
    parser.add_argument("--dwell_cutoff", type=float, default=1, help="cutoff of absolute dwell time difference ")
    parser.add_argument("--pvalue_cutoff", type=float, default=3, help="cutoff of -log10(adj pvalue)")
    parser.add_argument("--shift_size", type=int, default=4, help="shift size pf merging adjacent position")

    args = parser.parse_args()
    analyze_data(args.input, args.output, args.coverage_cutoff, args.mean_cutoff, args.dwell_cutoff, args.pvalue_cutoff, args.shift_size)
