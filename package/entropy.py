# Author:  GALEON Thomas
# Date:    06 May 2019
# Project: Population Genomics

try: from package.vcf_reader import *
except: from vcf_reader import *

def entropy_line(array):

    uni, cnt = np.unique(array, return_counts=True)

    return entropy(cnt/sum(cnt))

def chunk_entropy(filename, chunk=1000):

    writer, ind, size = None, 0, count_lines(filename)
    columns = collect_columns(filename)
    chromosome = filename.split('/')[-1].split('.')[0]
    output = '/'.join(['data', chromosome + '_entropy.pq'])

    print('Columns to be extracted: {}'.format(len(columns)))
    print('Number of lines to be read: {}'.format(size))

    with open(filename, 'r') as fle:

        for _ in tqdm.tqdm(np.arange(ceil(size / chunk))):

            lines = np.asarray([l[:-1].split('\t') for l in islice(fle, chunk) if not l.startswith('##')])
            if (ind == 0): lines = lines[1:]
            lines = pd.DataFrame(lines, columns=columns)  
            lines.set_index('POS', inplace=True)
            lines.drop(columns=['#CHROM', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'], inplace=True)
            lines = pd.DataFrame(lines.apply(entropy_line, axis=1))
            lines = pa.Table.from_pandas(lines, preserve_index=True)
            if writer is None: writer = pq.ParquetWriter(output, lines.schema)
            writer.write_table(lines)
            # Increase the counters
            ind += 1
            # Memory efficiency
            del lines

        writer.close()
        
if __name__ == '__main__':

    # Initialize the arguments
    prs = argparse.ArgumentParser()
    prs.add_argument('-f', '--file', help='Filename', type=str)
    prs.add_argument('-s', '--size', help='Chunk Size', type=int, default=1000)
    prs = prs.parse_args()

    # Build the dataframe
    chunk_entropy(prs.file, chunk=prs.size)