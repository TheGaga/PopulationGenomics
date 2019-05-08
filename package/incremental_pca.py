# Author:  DINDIN Meryll
# Date:    05 May 2019
# Project: PopulationGenomics

try: from package.imports import *
except: from imports import *

def build_estimator(filename, n_components=3, chunk=50):

    warnings.simplefilter('ignore')

    nme = filename.split('/')[-1].split('.')[0]
    lst = joblib.load('data/{}_patients.jb'.format(nme))
    pca = IncrementalPCA(n_components=n_components, copy=False)

    for index in tqdm.tqdm(range(len(lst)//chunk + 1)):

        vec = pd.read_parquet(filename, columns=lst[chunk*index:chunk*(index+1)]).values
        # Temporary file modification
        vec[vec == '0|0'] = 0
        vec[vec != 0] = 1
        vec = vec.astype('int8').transpose()
        # Run partial fit with approximated components
        pca.partial_fit(vec)
        # Memory efficiency
        del vec

    joblib.dump(pca, 'embedding/{}_pca_{}.jb'.format(nme, n_components))

def embed_chromosome(filename, n_components=3, chunk=50):

    col, res, nme = [], [], filename.split('/')[-1].split('.')[0]
    lst = joblib.load('data/{}_patients.jb'.format(nme))
    msk = ['#CHROM', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
    lst = [l for l in lst if l not in msk]
    pca = joblib.load('embedding/{}_pca_{}.jb'.format(nme, n_components))

    for index in tqdm.tqdm(range(len(lst)//chunk + 1)):

        beg, end = max(0, chunk*index), min(len(lst), chunk*(index+1))
        vec = pd.read_parquet(filename, columns=lst[beg:end]).values
        # Temporary file modification
        vec[vec == '0|0'] = 0
        vec[vec != 0] = 1
        vec = vec.astype('int8').transpose()
        print(len(lst[beg:end]), len(vec))
        # Transform the data based on approximated components
        res.append(pca.transform(vec))
        # Memory efficiency
        del vec

    res = pd.DataFrame(np.vstack(tuple(res)), index=col)
    res.to_pickle('embedding/{}_dimension_{}.df'.format(nme, n_components))

if __name__ == '__main__':

    # Initialize the arguments
    prs = argparse.ArgumentParser()
    prs.add_argument('-f', '--file', help='Filename', type=str)
    prs.add_argument('-s', '--size', help='Chunk Size', type=int, default=50)
    prs.add_argument('-c', '--npca', help='Number of components', type=int, default=3)
    prs = prs.parse_args()

    # Build the PCA estimator
    print('> Building the PCA estimator (n_components {}) for {} ...'.format(prs.npca, prs.file))
    build_estimator(prs.file, n_components=prs.npca, chunk=prs.size)

    # Build the pyarrow structure
    print('> Embed {} to {} ...'.format(prs.file, prs.npca))
    embed_chromosome(prs.file, n_components=prs.npca, chunk=prs.size)