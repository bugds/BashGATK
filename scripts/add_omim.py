import pandas as pd
import os

order = [
    'Проба',
    'Хромосома',
    'Позиция',
    'Реф_аллель',
    'Альт_аллель',
    'ID',
    'Аллельная_частота',
    'Глубина_прочтения',
    'Генотип',
    'Ген',
    'Последствие',
    'Кодирующее_последствие',
    'Полная_запись',
    'Макс_попул_ч-та',
    'In_silico_прогноз',
    'COSMIC_кодир',
    'COSMIC_некодир',
    'Сумма_COSMIC',
    'Клин_знач',
    'Клин_диаг',
    'rsID',
    'Коллер',
    'MIM_ID',
    'OMIM',
    'Missense_Z',
    'pLI',
    'SpliceAI',
    'Маска',
    'Число_проб_с_вариантом',
    'Доля_проб_с_вариантом',
    'Добро',
    'Доверие',
    'Target_region',
    'AlphaMissense',
    # 'PVS1',
    # 'PS1',
    # 'PS2',
    # 'PS3',
    # 'PS4',
    # 'PM1',
    # 'PM2',
    # 'PM3',
    # 'PM4',
    # 'PM5',
    # 'PM6',
    # 'PP1',
    # 'PP2',
    # 'PP3',
    # 'PP4',
    # 'PP5',
    'Баллы'
]

omim_file = os.environ['omim']
regions_file = os.environ['regions']
wd = os.environ['outputFolder']

def get_vardict(df):
    vardict = dict()
    grouped = df.groupby('Ген')
    for name, group in grouped:
        vardict[name] = len(group['ID'].unique())
    return vardict

def gene_style(v, vard):
    if vard[str(v)] > 1:
        return 'font-weight:bold;'
    return None

def genotype_style(v):
    if '1/1' in str(v):
        return 'font-weight:bold;'
    if '1|1' in str(v):
        return 'font-weight:bold;'
    return None

def omim_style(v):
    if 'dominant' in str(v).lower():
        return 'font-weight:bold;'
    return None

def add_omim(df):
    df['MIM_ID'] = ''
    df['OMIM'] = ''
    omim = pd.read_csv(omim_file, sep = '\t')
    for _, row in omim.iterrows():
        chr = row['Chromosome']
        start = row['Genomic Position Start']
        end = row['Genomic Position End']
        df.loc[
            (df['Хромосома'] == chr)
            & (df['Позиция'] <= end)
            & (df['Позиция'] >= start),
        'MIM_ID'] += str(row['MIM Number']) + ';'
        df.loc[
            (df['Хромосома'] == chr)
            & (df['Позиция'] <= end)
            & (df['Позиция'] >= start),
        'OMIM'] += row['Phenotypes'] + ';'
    df['MIM_ID'] = df['MIM_ID'].replace('', '.')
    df['OMIM'] = df['OMIM'].replace('', '.')
    return df

def get_region_info(df):
    df['Target_region'] = '.'
    regions = pd.read_csv(regions_file, sep = '\t', header = None)
    for _, row in regions.iterrows():
        chr = row[0]
        start = row[1]
        end = row[2]
        df.loc[
            (df['Хромосома'] == chr)
            & (df['Позиция'] <= end)
            & (df['Позиция'] >= start),
        'Target_region'] = 'ДА'
    return df

if __name__ == "__main__":
    df = pd.read_csv(os.path.join(wd, 'xl_results', 'results.tsv'), sep = '\t')

    samples = list(df['Проба'].unique())
    samples = [i + '.xlsx' for i in samples]

    for s in samples:
        print('Adding OMIM info to', s)
        df_map = pd.read_excel(os.path.join(wd, 'xl_results', s), sheet_name=None)
        with pd.ExcelWriter(
            os.path.join(wd, 'xl_results', s),
            mode='w'
        ) as writer:
            for k in df_map:
                print('List', k)
                # if k != 'all':
                df_map[k] = add_omim(df_map[k])
                df_map[k] = get_region_info(df_map[k])
                if len(df_map[k].columns.tolist()) != len(order):
                    raise Exception('Change order list!!!')
                df_map[k] = df_map[k][order]
                vardict = get_vardict(df_map[k])
                df_map[k] = df_map[k].fillna(".")
                df_map[k] = df_map[k].style \
                    .applymap(gene_style, vard = vardict, subset=pd.IndexSlice[:, ['Ген']])\
                    .applymap(genotype_style, subset=pd.IndexSlice[:, ['Генотип']])\
                    .applymap(omim_style, subset=pd.IndexSlice[:, ['OMIM']])
                df_map[k].to_excel(writer, sheet_name = k, index = None)
