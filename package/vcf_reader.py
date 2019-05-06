# Author:  DINDIN Meryll
# Date:    30 April 2019
# Project: Population Genomics

import os
import tqdm
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from io import StringIO
from itertools import islice

def count_lines(filename):

    return int(os.popen('wc -l {}'.format(filename)).read().split()[0])

def collect_columns(filename):

    with open(filename, 'r') as fle: lines = list(islice(fle, 500))
    lines = np.asarray([l[:-1].split('\t') for l in lines if not l.startswith('##')])

    return lines[0]

def chunk_transformer(filename, chunk=1000):

    writer, ind, size = None, 0, count_lines(filename)
    columns = collect_columns(filename)
    chromosome = filename.split('/')[-1].split('.')[0]
    output = '/'.join(['data', chromosome + '.pq'])

    print('Columns to be extracted: {}'.format(len(columns)))
    print('Number of lines to be read: {}'.format(size))

    with open(filename, 'r') as fle:

        for _ in tqdm.tqdm(np.arange((size // chunk)+1)):

            lines = np.asarray([l[:-1].split('\t') for l in islice(fle, chunk) if not l.startswith('##')])
            if (ind == 0): lines = lines[1:]
            lines = pd.DataFrame(lines, columns=columns)  
            lines.set_index('POS', inplace=True)
            lines.drop(columns=['#CHROM', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'], inplace=True)
            lines = pa.Table.from_pandas(lines, preserve_index=True)
            if writer is None: writer = pq.ParquetWriter(output, lines.schema)
            writer.write_table(lines)
            # Increase the counters
            ind += 1
            # Memory efficiency
            del lines

            writer.close()
            break

if __name__ == '__main__':

    chunk_transformer('data/chr20.vcf')