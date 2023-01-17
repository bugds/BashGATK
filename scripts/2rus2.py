import os
import sys
import pandas as pd

depth_limit = int(os.environ['depth_limit'])
print('Depth limit:', depth_limit)
af_limit = float(os.environ['af_limit'])
print('Allele frequency limit:', af_limit)
paf_limit = float(os.environ['paf_limit'])
print('Minor allele frequency limit:', paf_limit)
wd = os.environ['outputFolder']

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
        return 'Поврежд'
    elif ('T' in curr_preds) and not ('D' in curr_preds):
        return 'Нейтрал'
    else:
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
        dfd[k]['ID'] = dfd[k]['#CHROM'] \
            + '-' + dfd[k]['POS'].astype(str) \
            + '-' + dfd[k]['REF'] \
            + '-' + dfd[k]['ALT']
        dfd[k]['FORMAT_GT'] = '"' + dfd[k]['FORMAT_GT'] + '"'
        dfd[k] = dfd[k].rename(columns = rusDict)
        print('Columns renamed')
        for p in preds:
            dfd[k][p] = dfd[k][p].str.replace('H', 'D')
            dfd[k][p] = dfd[k][p].str.replace('B', 'T')
            dfd[k][p] = dfd[k][p].str.replace('N', 'T')
        print('Predictions reassessed')
        dfd[k]['In_silico_прогноз'] = ''
        dfd[k]['Макс_попул_ч-та'] = ''
        for c in [i for i in rusDict.values() if not (i in ['Позиция', 'Аллельная_частота', 'Глубина_прочтения'])]:
            dfd[k][c] = dfd[k][c].map(replace_x3b)
        print('Replaced HEX ;')
        dfd[k]['COSMIC_кодир'] = dfd[k]['COSMIC_кодир'].map(manage_x3d)
        dfd[k]['COSMIC_некодир'] = dfd[k]['COSMIC_некодир'].map(manage_x3d)
        dfd[k]['Сумма_COSMIC'] = dfd[k][['COSMIC_кодир', 'COSMIC_некодир']].astype(str).apply(cosmsum, axis=1)
        print('Replaced HEX =')
        dfd[k]['In_silico_прогноз'] = dfd[k][preds].astype(str).apply(''.join, axis=1)
        dfd[k]['In_silico_прогноз'] = dfd[k]['In_silico_прогноз'].map(manage_preds)
        print('Single prognosis done')
        dfd[k]['Макс_попул_ч-та_1'] = dfd[k]['Макс_попул_ч-та_1'].replace('.', -1)
        dfd[k]['Макс_попул_ч-та_1_alt'] = dfd[k]['Макс_попул_ч-та_1_alt'].replace('.', -1)
        dfd[k]['Макс_попул_ч-та_2'] = dfd[k]['Макс_попул_ч-та_2'].replace('.', -1)
        dfd[k]['Макс_попул_ч-та'] = dfd[k][['Макс_попул_ч-та_1', 'Макс_попул_ч-та_2', 'Макс_попул_ч-та_1_alt']].astype(float).apply(max, axis=1)
        print('Single MAF done')
        dfd[k] = dfd[k].drop(columns = preds)
        dfd[k] = dfd[k].drop(columns = ['Макс_попул_ч-та_1', 'Макс_попул_ч-та_1_alt', 'Макс_попул_ч-та_2'])
        print('Excessive columns dropped')
        dfd[k]['Доверие'] = ''
        dfd[k].loc[dfd[k]['Глубина_прочтения'] < depth_limit, 'Доверие'] = 'Низк'
        dfd[k].loc[dfd[k]['Аллельная_частота'] < af_limit ,'Доверие'] = 'Низк'
        dfd[k].loc[dfd[k]['Макс_попул_ч-та'].astype(float) > paf_limit, 'Доверие'] = 'Низк'
        dfd[k].loc[dfd[k]['Доверие'] == '', 'Доверие'] = 'Выс'
        dfd[k]['Макс_попул_ч-та'] = dfd[k]['Макс_попул_ч-та'].replace(-1, '.')
        print('Evaluation complete with:')
        print('Minimum depth', str(depth_limit))
        print('Minimum AF', str(af_limit))
        print('Maximum population AF', str(paf_limit))
        dfd[k]['temp_id'] = dfd[k]['Хромосома'] + ':' + dfd[k]['Позиция'].astype(str) + ':' + dfd[k]['Реф_аллель'] + '>' + dfd[k]['Альт_аллель'] + '_' + dfd[k]['Коллер']
        counts = dfd[k]['temp_id'].value_counts()
        dfd[k]['Число_проб_с_вариантом'] = dfd[k]['temp_id'].map(counts)
        sample_num = len(dfd[k]['Проба'].unique())
        dfd[k]['Доля_проб_с_вариантом'] = dfd[k]['Число_проб_с_вариантом'] / sample_num
        print('Added in-house frequency')
        rusDictValues = [v for v in dfd[k].columns if v in rusDict.values()]
        orderedValues = list()
        for v in rusDict.values():
            if v in rusDictValues:
                if not (v in orderedValues):
                    orderedValues.append(v)
                    if v == 'COSMIC_некодир':
                        orderedValues.append('Сумма_COSMIC')
        orderedValues.append('Коллер')
        orderedValues.append('Доверие')
        orderedValues.append('Число_проб_с_вариантом')
        orderedValues.append('Доля_проб_с_вариантом')
        dfd[k] = dfd[k][orderedValues]

    df = pd.concat(dfd.values())
    df = df.sort_values(by = ['Хромосома', 'Позиция'])

    df.to_csv(wd + '/rus2.tsv', sep = '\t', index = False)
