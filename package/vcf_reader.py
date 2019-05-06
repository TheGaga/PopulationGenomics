# Author:  DINDIN Meryll
# Date:    30 April 2019
# Project: Population Genomics

try: from package.imports import *
except: from imports import *

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

if __name__ == '__main__':

    # Initialize the arguments
    prs = argparse.ArgumentParser()
    prs.add_argument('-f', '--file', help='Filename', type=str)
    prs.add_argument('-s', '--size', help='Chunk Size', type=int, default=1000)
    prs = prs.parse_args()

    # Serialize the list of individuals
    out = '.'.join(prs.file.split('.')[:-1]) + '.jb'
    joblib.dump(collect_columns(prs.file), out)

    # Build the pyarrow structure
    chunk_transformer(prs.file, prs.size)