import os
import sys
import pandas as pd

depth_limit = 100
af_limit = 0.02
paf_limit = 0.01

wd = os.path.abspath(sys.argv[1])

rusDict = {
    'SAMPLE': 'Проба',
    '#CHROM': 'Хромосома',
    'POS': 'Позиция',
    'REF': 'Реф_аллель',
    'ALT': 'Альт_аллель',
    'FORMAT_AF': 'Аллельная_частота',
    'FORMAT_VAF': 'Аллельная_частота',
    'FORMAT_DP': 'Глубина_прочтения',
    'INFO_ANNO_Gene.refGene': 'Ген',
    'INFO_ANNO_Func.refGene': 'Последствие',
    'INFO_ANNO_ExonicFunc.refGene': 'Кодирующее_последствие',
    'INFO_ANNO_AAChange.refGene': 'Полная_запись',
    'INFO_ANNO_AF_popmax': 'Макс_попул_ч-та_1',
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
    'INFO_ANNO_CLNSIG': 'Клин_знач_1',
    'INFO_VEP_CLIN_SIG': 'Клин_знач_2',
    'Клин_знач': 'Клин_знач',
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

dfd = dict()
dfd['germ'] = pd.read_csv(wd + '/combined_passed_anno_germ.tsv', sep = '\t')
dfd['germ']['Коллер'] = 'DeepVariant'
dfd['soma'] = pd.read_csv(wd + '/combined_passed_anno_soma.tsv', sep = '\t')
dfd['soma']['Коллер'] = 'Mutect2'

for k in dfd:
    dfd[k] = dfd[k].rename(columns = rusDict)
    for p in preds:
        dfd[k][p] = dfd[k][p].str.replace('H', 'D')
        dfd[k][p] = dfd[k][p].str.replace('B', 'T')
        dfd[k][p] = dfd[k][p].str.replace('N', 'T')
    dfd[k]['In_silico_прогноз'] = ''
    dfd[k]['Макс_попул_ч-та'] = ''
    dfd[k]['Клин_знач'] = ''
    for i, r in dfd[k].iterrows():
        for c in rusDict.values():
            if isinstance(r[c], str):
                dfd[k].loc[i, c] = r[c].replace('\\x3b', ';')
        for cosm in ['COSMIC_кодир', 'COSMIC_некодир']:
            if '\\x3d' in r[cosm]:
                dfd[k].loc[i, cosm] = r[cosm].split('\\x3d')[1]
            dfd[k].loc[i, cosm] = dfd[k].loc[i, cosm].split('\\x3b')[0]
        result = '.'
        curr_preds = set([r[p] for p in preds])
        if ('D' in curr_preds) and not ('T' in curr_preds):
            result = 'Поврежд'
        elif ('T' in curr_preds) and not ('D' in curr_preds):
            result = 'Толерантн'
        dfd[k].loc[i, 'In_silico_прогноз'] = result
        if r['Макс_попул_ч-та_1'] == '.':
            dfd[k].loc[i, 'Макс_попул_ч-та'] = r['Макс_попул_ч-та_2']
        elif r['Макс_попул_ч-та_2'] == '.':
            dfd[k].loc[i, 'Макс_попул_ч-та'] = r['Макс_попул_ч-та_1']
        else:
            dfd[k].loc[i, 'Макс_попул_ч-та'] = max([float(r['Макс_попул_ч-та_1']), float(r['Макс_попул_ч-та_2'])])
        if r['Клин_знач_2'].lower() == r['Клин_знач_1'].lower():
            dfd[k].loc[i, 'Клин_знач'] = r['Клин_знач_1'].lower()
        elif r['Клин_знач_2'] == '.':
            dfd[k].loc[i, 'Клин_знач'] = '?' + r['Клин_знач_1'].lower()
        elif r['Клин_знач_1'] == '.':
            dfd[k].loc[i, 'Клин_знач'] = '?' + r['Клин_знач_2'].lower()
        else:
            dfd[k].loc[i, 'Клин_знач'] = '?' + r['Клин_знач_1'].lower()
    dfd[k] = dfd[k].drop(columns = preds)
    dfd[k] = dfd[k].drop(columns = ['Макс_попул_ч-та_1', 'Макс_попул_ч-та_2', 'Клин_знач_1', 'Клин_знач_2'])
    dfd[k]['Доверие'] = ''
    dfd[k]['Макс_попул_ч-та'] = dfd[k]['Макс_попул_ч-та'].replace('.', -1)
    dfd[k].loc[dfd[k]['Глубина_прочтения'] < depth_limit, 'Доверие'] = 'Низк'
    dfd[k].loc[dfd[k]['Аллельная_частота'] < af_limit ,'Доверие'] = 'Низк'
    dfd[k].loc[dfd[k]['Макс_попул_ч-та'].astype(float) > paf_limit, 'Доверие'] = 'Низк'
    dfd[k].loc[dfd[k]['Доверие'] == '', 'Доверие'] = 'Выс'
    dfd[k]['Макс_попул_ч-та'] = dfd[k]['Макс_попул_ч-та'].replace(-1, '.')
    rusDictValues = [v for v in dfd[k].columns if v in rusDict.values()]
    orderedValues = list()
    for v in rusDict.values():
        if v in rusDictValues:
            if not (v in orderedValues):
                orderedValues.append(v)
    orderedValues.append('Коллер')
    orderedValues.append('Доверие')
    dfd[k] = dfd[k][orderedValues]

df = pd.concat(dfd.values())
df = df.sort_values(by = ['Хромосома', 'Позиция'])

df.to_csv(wd + '/rus2.tsv', sep = '\t', index = False)