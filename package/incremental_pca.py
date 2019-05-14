# Author:  DINDIN Meryll
# Date:    05 May 2019
# Project: PopulationGenomics

try: from package.vcf_reader import *
except: from vcf_reader import *

def continent_individuals(continent):

    phe = pd.read_csv('data/phenotypes.ped', sep='\t')[['Individual ID', 'Population']]
    pop = pd.read_csv('data/populations.tsv', sep='\t')[['Population Code', 'Super Population']]
    pop = phe.merge(pop, how='left', left_on='Population', right_on='Population Code')
    pop.drop(['Population', 'Population Code'], axis=1, inplace=True)
    
    return list(pop[pop['Super Population'] == continent]['Individual ID'].values.ravel())

def build_estimator(filename, n_components=3, chunk=50, trim=None):

    warnings.simplefilter('ignore')

    nme = filename.split('/')[-1].split('.')[0]
    pca = IncrementalPCA(n_components=n_components, copy=False)

    if trim is None: 
        out = list_patients()
        pop = out.Population.values.ravel()
        out = out.index.values.ravel()
        ser = 'embedding/{}_pca_{}.jb'.format(nme, n_components)
    else: 
        out = np.asarray(continent_individuals(trim))
        pop = [trim for _ in range(len(out))]
        ser = 'embedding/{}_{}_pca_{}.jb'.format(nme, trim, n_components)

    skf = StratifiedKFold(n_splits=len(out)//chunk, shuffle=True, random_state=42)
    for _, lst in tqdm.tqdm(skf.split(out, pop)):

        vec = pd.read_parquet(filename, columns=out[lst]).values
        # Temporary file modification
        vec[vec == '0|0'] = 0
        vec[vec != 0] = 1
        vec = vec.astype('int8').transpose()
        # Run partial fit with approximated components
        pca.partial_fit(vec)
        # Memory efficiency
        del vec

    joblib.dump(pca, ser)

def embed_chromosome(filename, n_components=3, chunk=50, trim=None):

    res, nme = [], filename.split('/')[-1].split('.')[0]

    # Select the right patients
    if trim is None: 
        lst = joblib.load('data/{}_patients.jb'.format(nme))
        pca = joblib.load('embedding/{}_pca_{}.jb'.format(nme, n_components))
        out = 'embedding/{}_dimension_{}.df'.format(nme, n_components)
    else: 
        lst = continent_individuals(trim)
        pca = joblib.load('embedding/{}_{}_pca_{}.jb'.format(nme, trim, n_components))
        out = 'embedding/{}_{}_dimension_{}.df'.format(nme, trim, n_components)

    # Apply column filtering
    msk = ['#CHROM', 'ID', 'POS', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT']
    lst = [l for l in lst if l not in msk]

    for index in tqdm.tqdm(range(len(lst)//chunk + 1)):

        beg, end = max(0, chunk*index), min(len(lst), chunk*(index+1))
        vec = pd.read_parquet(filename, columns=lst[beg:end]).values
        # Temporary file modification
        vec[vec == '0|0'] = 0
        vec[vec != 0] = 1
        vec = vec.astype('int8').transpose()
        # Transform the data based on approximated components
        res.append(pca.transform(vec))
        # Memory efficiency
        del vec

    res = pd.DataFrame(np.vstack(tuple(res)), index=lst)
    res.to_pickle(out)

if __name__ == '__main__':

    # Initialize the arguments
    prs = argparse.ArgumentParser()
    prs.add_argument('-f', '--file', help='Filename', type=str)
    prs.add_argument('-s', '--size', help='Chunk size', type=int, default=50)
    prs.add_argument('-t', '--trim', help='Continent specific', default=None)
    prs.add_argument('-c', '--npca', help='Number of components', type=int, default=3)
    prs = prs.parse_args()

    # Build the PCA estimator
    print('> Building the PCA estimator (n_components {}) for {} ...'.format(prs.npca, prs.file))
    build_estimator(prs.file, n_components=prs.npca, chunk=prs.size)

    # Build the pyarrow structure
    print('> Embed {} to {} ...'.format(prs.file, prs.npca))
    embed_chromosome(prs.file, n_components=prs.npca, chunk=prs.size, trim=prs.trim)