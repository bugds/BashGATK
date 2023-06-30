import os
import pandas as pd
from add_omim import get_vardict, gene_style, genotype_style, omim_style

phylongs = os.environ['phylongs']
wd = os.environ['outputFolder']

if __name__ == "__main__":
    with open(phylongs, 'r') as inp:
        phylongs_ids = set([l.strip() for l in inp.readlines()])

    df = pd.read_csv(os.path.join(wd, 'xl_results', 'results.tsv'), sep = '\t')

    samples = list(df['Проба'].unique())
    samples = [i + '.xlsx' for i in samples]

    for s in samples:
        print('Adding PhyloNGS info to', s)
        df_map = pd.read_excel(os.path.join(wd, 'xl_results', s), sheet_name=None)
        with pd.ExcelWriter(
            os.path.join(wd, 'xl_results', 'plus' + s),
            mode='w'
        ) as writer:
            for k in df_map:
                ids = set(df_map[k]['ID'])
                phylo_ids = ids.intersection(phylongs_ids)
                print(str(len(phylo_ids)), 'phylogenetic polymorphisms in list', k, "(out of {})".format(str(len(ids))))
                mapdict = dict()
                for i in ids:
                    mapdict[i] = '.'
                    if i in phylo_ids:
                        mapdict[i] = '1'
                df_map[k]['phyloNGS'] = df_map[k]['ID'].map(mapdict)
                vardict = get_vardict(df_map[k])
                if k != 'all':
                    df_map[k] = df_map[k].style \
                        .applymap(gene_style, vard = vardict, subset=pd.IndexSlice[:, ['Ген']])\
                        .applymap(genotype_style, subset=pd.IndexSlice[:, ['Генотип']])\
                        .applymap(omim_style, subset=pd.IndexSlice[:, ['OMIM']])
                df_map[k].to_excel(writer, sheet_name = k, index = None)
