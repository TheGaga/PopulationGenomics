# Author:  DINDIN Meryll
# Date:    30 April 2019
# Project: Population Genomics

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from io import StringIO
from itertools import islice

def vcf_to_pandas(lines):
    """
    Read a .vcf file as pandas DataFrame
    Input:
        lines: list of strings
    Returns:
        data: pd.DataFrame
    """
    # Filter lines
    lines = [l for l in lines if not l.startswith('##')]
    # Transform into dataframe
    frame = pd.read_csv(StringIO(''.join(lines)), dtype={'POS': int}, sep='\t')
    frame.set_index('POS', inplace=True)
    frame.drop(columns=['#CHROM', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'], inplace=True)
    # Memory efficiency
    del lines

    return frame

def chunk_transformer(filename, size=25000):

    writer = None
    chromosome = filename.split('/')[-1].split('.')[0]
    output = '/'.join(['data', chromosome + '.pq'])

    with open(filename, 'r') as fle:

        while True:

            lines = list(islice(fle, size))
            if not lines: break

            lines = pa.Table.from_pandas(vcf_to_pandas(lines), preserve_index=True)
            if writer is None: writer = pq.ParquetWriter(output, lines.schema)
            writer.write_table(lines)

        writer.close()   

if __name__ == '__main__':

    chunk_transformer('data/chr20.vcf')