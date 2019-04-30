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
    # Dictionnary to map values to integers
    mapper = {'0|0': 0, '0|1': 1, '0|2': 2, '1|0': 3, '1|1': 4, '1|2': 5, '2|0': 6, '2|1': 7, '2|2': 8,
              '0|3': 9, '1|3': 10, '2|3': 11, '3|0': 12, '3|1': 13, '3|2': 14, '3|3': 15}
    # Filter lines
    lines = [l for l in lines if not l.startswith('##')]
    # Transform into dataframe
    frame = pd.read_csv(StringIO(''.join(lines)), dtype={'POS': int}, sep='\t')
    frame.set_index('POS', inplace=True)
    frame.drop(columns=['#CHROM', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'], inplace=True)
    frame.replace(to_replace=mapper, inplace=True)
    # Memory efficiency
    del lines

    return frame.astype('int8')

def chunk_transformer(filename, size=25000):

    writer = None
    chromosome = filename.split('/')[-1].split('.')[0]
    output = '/'.join(['data', chromosome, '.pq'])

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