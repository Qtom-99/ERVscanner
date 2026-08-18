import pandas as pd
import argparse


def create_position_sample_matrix(data_file, samples_file, output_file):

    # Read position IDs and sample names
    df = pd.read_csv(
        data_file,
        sep=r'\s+',
        header=None,
        names=['Position', 'Sample']
    )

    # Read the sample list
    with open(samples_file, 'r') as f:
        samples = f.read().splitlines()

    # Sort position IDs in natural numerical order
    # while preserving the original POS IDs
    position_order = sorted(
        df['Position'].unique(),
        key=lambda x: int(x[3:])
    )

    df['Position'] = pd.Categorical(
        df['Position'],
        categories=position_order,
        ordered=True
    )

    # Create a sample-by-position matrix
    pivot_df = df.pivot_table(
        index='Sample',
        columns='Position',
        aggfunc='size',
        fill_value=0,
        observed=False
    )

    # Convert counts to binary presence/absence values
    pivot_df = (pivot_df > 0).astype(int)

    # Ensure that all samples in the sample list are retained
    samples_df = pd.DataFrame(
        samples,
        columns=['Sample']
    )

    pivot_df = (
        samples_df
        .merge(
            pivot_df,
            on='Sample',
            how='left'
        )
        .fillna(0)
        .set_index('Sample')
        .astype(int)
    )

    # Write the matrix as a tab-delimited file
    pivot_df.to_csv(
        output_file,
        sep='\t'
    )

    return pivot_df


def main():

    parser = argparse.ArgumentParser(
        description="Create a position-by-sample presence/absence matrix."
    )

    parser.add_argument(
        'data_file',
        help='Input file containing position IDs and sample names.'
    )

    parser.add_argument(
        'samples_file',
        help='File containing the list of samples.'
    )

    parser.add_argument(
        'output_file',
        help='Output TSV file.'
    )

    args = parser.parse_args()

    create_position_sample_matrix(
        args.data_file,
        args.samples_file,
        args.output_file
    )


if __name__ == "__main__":
    main()
