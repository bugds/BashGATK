import os
import sys
import pandas as pd
import acmg_classifier as ac

depth_limit = int(os.environ['depth_limit'])
print('Depth limit:', depth_limit)
af_limit = float(os.environ['af_limit'])
print('Allele frequency limit:', af_limit)
paf_limit = float(os.environ['paf_limit'])
print('Minor allele frequency limit:', paf_limit)
wd = os.environ['outputFolder']
freq_file = os.environ['freq_file']

rusDict = {
    'SAMPLE': 'Проба',
    '#CHROM': 'Хромосома',
    'POS': 'Позиция',
    'REF': 'Реф_аллель',
    'ALT': 'Альт_аллель',
    'ID': 'ID',
    'FORMAT_AF': 'Аллельная_частота',
    'FORMAT_VAF': 'Аллельная_частота',
    'FORMAT_DP': 'Глубина_прочтения',
    'FORMAT_GT': 'Генотип',
    'INFO_ANNO_Gene.refGene': 'Ген',
    'INFO_ANNO_Func.refGene': 'Последствие',
    'INFO_ANNO_ExonicFunc.refGene': 'Кодирующее_последствие',
    'INFO_ANNO_AAChange.refGene': 'Полная_запись',
    'INFO_ANNO_AF_popmax': 'Макс_попул_ч-та_1',
    'INFO_ANNO_AF': 'Макс_попул_ч-та_1_alt',
    'INFO_VEP_MAX_AF': 'Макс_попул_ч-та_2',
    'Макс_попул_ч-та': 'Макс_попул_ч-та',
    'In_silico_прогноз': 'In_silico_прогноз',
    'INFO_ANNO_SIFT_pred': 'PredSIFT',
    'INFO_ANNO_SIFT4G_pred': 'PredSIFT4G',
    'INFO_ANNO_Polyphen2_HDIV_pred': 'PredPP2HDIV',
    'INFO_ANNO_Polyphen2_HVAR_pred': 'PredPP2HVAR',
    'INFO_ANNO_LRT_pred': 'PredLRT',
    'INFO_ANNO_MutationTaster_pred': 'PredMT',
    'INFO_ANNO_MutationAssessor_pred': 'PredMA',
    'INFO_ANNO_FATHMM_pred': 'PredFATHMM',
    'INFO_ANNO_PROVEAN_pred': 'PredPROVEAN',
    'INFO_ANNO_MetaSVM_pred': 'PredMetaSVM',
    'INFO_ANNO_MetaLR_pred': 'PredMetaLR',
    'INFO_ANNO_M-CAP_pred': 'PredM-CAP',
    'INFO_ANNO_fathmm-MKL_coding_pred': 'PredFATHMM-MKL',
    'INFO_ANNO_fathmm-XF_coding_pred': 'PredFATHMM-XF',
    'INFO_ANNO_cosmic95_coding': 'COSMIC_кодир',
    'INFO_ANNO_cosmic95_noncoding': 'COSMIC_некодир',
    'INFO_ANNO_CLNSIG': 'Клин_знач',
    'INFO_ANNO_CLNDN': 'Клин_диаг',
    'INFO_ANNO_avsnp150': 'rsID'
}

preds = [
    'PredSIFT',
    'PredSIFT4G',
    'PredPP2HDIV',
    'PredPP2HVAR',
    'PredLRT',
    'PredMT',
    'PredMA',
    'PredFATHMM',
    'PredPROVEAN',
    'PredMetaSVM',
    'PredMetaLR',
    'PredM-CAP',
    'PredFATHMM-MKL',
    'PredFATHMM-XF'
]

def replace_x3b(string):
    if isinstance(string, str):
        string = string.replace('\\', '/')
        return string.replace('/x3b', ';')

def manage_x3d(string):
    if 'x3d' in string:
        string = string.replace('\\', '/')
        #string = string.split('/x3d')[1].split(';')[0]
    return string.replace('/x3d', '=')

def manage_preds(curr_preds):
    if ('D' in curr_preds) and not ('T' in curr_preds):
        if curr_preds.count('D') >= 3:
            return 'Поврежд'
    elif ('T' in curr_preds) and not ('D' in curr_preds):
        if curr_preds.count('T') >= 3:
            return 'Нейтрал'
    return '.'

def cosmsum(cosmlist):
    sums = list()
    for cosm in cosmlist:
        curr_sum = 0
        if '=' in cosm:
            cosmdigits = cosm.split('=')[-1].split(',')
            cosmdigits = [int(d.split('(')[0]) for d in cosmdigits]
            curr_sum = sum(cosmdigits)
        sums.append(curr_sum)
    return max(sums)

def check_BA1(maf):
    if maf > 0.03:
        return 1
    return 0

def check_BS1(maf):
    if maf > 0.005:
        return 1
    return 0

def check_BP4(insilico):
    if insilico == 'Нейтрал':
        return 1
    return 0

def check_BP6(clinvar):
    if ('benign' in clinvar.lower()) and (not ('pathogenic') in clinvar.lower()):
        return 1
    return 0

def checkB(row):
    if row['BA1'] == 1:
        return 1
    else:
        if sum([row['BS1'], row['BP4'], row['BP6']]) >= 2:
            return 1
    return "."

def add_id(df):
    df['ID'] = df['#CHROM'] \
        + '-' + df['POS'].astype(str) \
        + '-' + df['REF'] \
        + '-' + df['ALT']
    return df

def correct_GT(df):
    df['FORMAT_GT'] = '"' + df['FORMAT_GT'] + '"'
    return df

def add_predictions(df):
    for p in preds:
        df[p] = df[p].str.replace('H', 'D')
        df[p] = df[p].str.replace('B', 'T')
        df[p] = df[p].str.replace('N', 'T')
    print('Predictions reassessed')
    df['In_silico_прогноз'] = ''
    df['In_silico_прогноз'] = df[preds].astype(str).apply(''.join, axis=1)
    df['In_silico_прогноз'] = df['In_silico_прогноз'].map(manage_preds)
    df = df.drop(columns = preds)
    print('Single prognosis done')
    return df

def replace_hexs(df):
    for c in [i for i in rusDict.values() if not (i in ['Позиция', 'Аллельная_частота', 'Глубина_прочтения', 'Макс_попул_ч-та', *preds])]:
        df[c] = df[c].map(replace_x3b)
    print('Replaced HEX ;')
    df['COSMIC_кодир'] = df['COSMIC_кодир'].map(manage_x3d)
    df['COSMIC_некодир'] = df['COSMIC_некодир'].map(manage_x3d)
    df['Сумма_COSMIC'] = df[['COSMIC_кодир', 'COSMIC_некодир']].astype(str).apply(cosmsum, axis=1)
    print('Replaced HEX =')
    return df

def add_popfreq(df):
    df['Макс_попул_ч-та'] = ''
    df['Макс_попул_ч-та_1'] = df['Макс_попул_ч-та_1'].replace('.', -1)
    df['Макс_попул_ч-та_1_alt'] = df['Макс_попул_ч-та_1_alt'].replace('.', -1)
    df['Макс_попул_ч-та_2'] = df['Макс_попул_ч-та_2'].replace('.', -1)
    df['Макс_попул_ч-та'] = df[['Макс_попул_ч-та_1', 'Макс_попул_ч-та_2', 'Макс_попул_ч-та_1_alt']].astype(float).apply(max, axis=1)
    df = df.drop(columns = ['Макс_попул_ч-та_1', 'Макс_попул_ч-та_1_alt', 'Макс_попул_ч-та_2'])
    print('Single MAF done')
    return df

def add_trust(df):
    df['Доверие'] = ''
    df.loc[df['Глубина_прочтения'] < depth_limit, 'Доверие'] = 'Низк'
    df.loc[df['Аллельная_частота'] < af_limit ,'Доверие'] = 'Низк'
    df.loc[df['Макс_попул_ч-та'].astype(float) > paf_limit, 'Доверие'] = 'Низк'
    df.loc[df['Доверие'] == '', 'Доверие'] = 'Выс'
    return df

def add_inhouse_freq(df):
    df['temp_id'] = df['Хромосома'] + ':' + df['Позиция'].astype(str) + ':' + df['Реф_аллель'] + '>' + df['Альт_аллель'] + '_' + df['Коллер']
    freq = pd.read_csv(freq_file, sep = '\t')
    freq = pd.concat([freq, df[['temp_id', 'Проба']]])
    freq = freq.drop_duplicates()
    freq.to_csv(freq_file, index = False, sep = '\t')
    counts = freq['temp_id'].value_counts()
    df['Число_проб_с_вариантом'] = df['temp_id'].map(counts)
    sample_num = len(freq['Проба'].unique())
    print('Total number of samples:', sample_num)
    df['Доля_проб_с_вариантом'] = df['Число_проб_с_вариантом'] / sample_num
    print('Added in-house frequency')
    return df

def order_cols(df):
    rusDictValues = [v for v in df.columns if v in rusDict.values()]
    orderedValues = list()
    for v in rusDict.values():
        if v in rusDictValues:
            if not (v in orderedValues):
                orderedValues.append(v)
                if v == 'COSMIC_некодир':
                    orderedValues.append('Сумма_COSMIC')
    df['Добро'] = '.'
    df['Маска'] = ""
    orderedValues.append('Коллер')
    orderedValues.append('Доверие')
    orderedValues.append('Число_проб_с_вариантом')
    orderedValues.append('Доля_проб_с_вариантом')
    orderedValues.append('Добро')
    orderedValues.append('Маска')
    orderedValues.append('Missense_Z')
    orderedValues.append('pLI')
    orderedValues.append('PVS1')
    orderedValues.append('PS1')
    orderedValues.append('PS2')
    orderedValues.append('PS3')
    orderedValues.append('PS4')
    orderedValues.append('PM1')
    orderedValues.append('PM2')
    orderedValues.append('PM3')
    orderedValues.append('PM4')
    orderedValues.append('PM5')
    orderedValues.append('PM6')
    orderedValues.append('PP1')
    orderedValues.append('PP2')
    orderedValues.append('PP3')
    orderedValues.append('PP4')
    orderedValues.append('PP5')
    df = df[orderedValues]
    return df

def classify_benign(dfd):
    dfd['germ']['BA1'] = dfd['germ']['Макс_попул_ч-та'].astype(float).map(check_BA1)
    dfd['germ']['BS1'] = dfd['germ']['Макс_попул_ч-та'].astype(float).map(check_BS1)
    dfd['germ']['BP4'] = dfd['germ']['In_silico_прогноз'].map(check_BP4)
    dfd['germ']['BP6'] = dfd['germ']['Клин_знач'].map(check_BP6)
    dfd['germ']['Добро'] = dfd['germ'].apply(checkB, axis = 1)
    dfd['soma']['Добро'] = "."
    dfd['germ'] = dfd['germ'].drop(columns = ['BA1', 'BS1', 'BP4', 'BP6'])
    dfd['germ']['Макс_попул_ч-та'] = dfd['germ']['Макс_попул_ч-та'].replace(-1, '.')
    dfd['soma']['Макс_попул_ч-та'] = dfd['soma']['Макс_попул_ч-та'].replace(-1, '.')
    return dfd

dfd = dict()
if os.path.exists(wd + '/combined_passed_anno_germ.tsv'):
    dfd['germ'] = pd.read_csv(wd + '/combined_passed_anno_germ.tsv', sep = '\t')
    dfd['germ']['Коллер'] = 'DeepVariant'
if os.path.exists(wd + '/combined_passed_anno_soma.tsv'):
    dfd['soma'] = pd.read_csv(wd + '/combined_passed_anno_soma.tsv', sep = '\t')
    dfd['soma']['Коллер'] = 'Mutect2'

if len(dfd) == 0:
    print('No files to analyze')
else:
    for k in dfd:
        dfd[k] = add_id(dfd[k])
        dfd[k] = correct_GT(dfd[k])
        dfd[k] = dfd[k].rename(columns = rusDict)
        dfd[k] = add_predictions(dfd[k])
        dfd[k] = replace_hexs(dfd[k])
        dfd[k] = add_popfreq(dfd[k])
        dfd[k] = add_trust(dfd[k])
        dfd[k] = add_inhouse_freq(dfd[k])
        dfd[k] = ac.classify(dfd[k])
        dfd[k] = dfd[k].fillna(".")
        dfd[k] = order_cols(dfd[k])
    dfd = classify_benign(dfd)
    df = pd.concat(dfd.values())
    df = df.sort_values(by = ['Хромосома', 'Позиция'])
    df.to_csv(wd + '/rus2.tsv', sep = '\t', index = False)
