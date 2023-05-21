#!/usr/bin/env python3

import pandas as pd
import argparse

def load_data(input_file, metadata_file):
    # Load input_file data
    blast_results = pd.read_csv(input_file, sep='\t', header=None, names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])

    # Load metadata_file data
    card_metadata = pd.read_csv(metadata_file, sep='\t')

    print("----Blast Results----")
    print(blast_results.head(10))

    print("----CARD Metadata----")
    print(card_metadata.head(10))

    return blast_results, card_metadata

def process_data(blast_results, card_metadata):
    # Extract CARD Short Name from sseqid column in blast_results
    blast_results['CARD Short Name'] = blast_results['sseqid'].apply(lambda x: x.split('|')[-1])

    # Merge the blast_results and card_metadata dataframes on CARD Short Name
    merged_data = pd.merge(blast_results, card_metadata, on='CARD Short Name')

    # Keep only the required columns
    merged_data = merged_data[['qseqid', 'Drug Class', 'Resistance Mechanism']]

    print("----Merged Data----")
    print(merged_data.head(10))

    return merged_data

def save_data(merged_data, output_file):
    # Save the merged_data dataframe to the output_file
    merged_data.to_csv(output_file, sep='\t', index=False)

def main():
    parser = argparse.ArgumentParser(description='Annotate the diamond blast result with a CARD_metadata file')
    parser.add_argument('input_file', help='The input file containing the diamond blast result (e.g., A_genome_CARD_results.txt)')
    parser.add_argument('metadata_file', help='The metadata file (e.g., CARD_metadata.tsv)')
    parser.add_argument('output_file', help='The output file to save the annotated result (e.g., A_genome_summary.tsv)')

    args = parser.parse_args()

    # Load data
    blast_results, card_metadata = load_data(args.input_file, args.metadata_file)

    # Process data
    merged_data = process_data(blast_results, card_metadata)

    # Save data
    save_data(merged_data, args.output_file)

if __name__ == "__main__":
    main()
