try: from package.imports import *
except: from imports import *

def count_lines(filename):

    return int(os.popen('wc -l {}'.format(filename)).read().split()[0])

def collect_columns(filename):

    with open(filename, 'r') as fle: lines = list(islice(fle, 500))
    lines = np.asarray([l[:-1].split('\t') for l in lines if not l.startswith('##')])

    return lines[0]

def entropy_line(array):

    uni, cnt = np.unique(array, return_counts=True)
    return entropy(cnt/sum(cnt))


def chunk_transformer(filename, chunk=1000):

    data, ind, size = None, 0, count_lines(filename)
    columns = collect_columns(filename)
    chromosome = filename.split('/')[-1].split('.')[0]
    output = '/'.join(['data', chromosome + '_entropy.csv'])

    print('Columns to be extracted: {}'.format(len(columns)))
    print('Number of lines to be read: {}'.format(size))

    with open(filename, 'r') as fle:

        for _ in tqdm.tqdm(np.arange(ceil(size / chunk))):

            lines = np.asarray([l[:-1].split('\t') for l in islice(fle, chunk) if not l.startswith('##')])
            if (ind == 0): lines = lines[1:]
            lines = pd.DataFrame(lines, columns=columns)  
            lines.set_index('POS', inplace=True)
            lines.drop(columns=['#CHROM', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT'], inplace=True)
            lines = lines.apply(entropy_line, axis=1)
            if data is None: data = lines
            else: data = pd.concat([data, lines])
            # Increase the counters
            ind += 1
            # Memory efficiency
            del lines
            time.sleep(1)

        data.to_csv(output)
        

if __name__ == '__main__':

    # Initialize the arguments
    prs = argparse.ArgumentParser()
    prs.add_argument('-f', '--file', help='Filename', type=str)
    prs.add_argument('-s', '--size', help='Chunk Size', type=int, default=1000)
    prs = prs.parse_args()

    # Build the dataframe
    chunk_transformer(prs.file, prs.size)