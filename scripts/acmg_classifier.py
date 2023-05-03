import os
from io import StringIO
import subprocess
import pandas as pd

clinvar = os.environ['clinvar']
bedtools = os.environ['bedtools']
constraint = os.environ['constraint']

class ACMG_criteria(dict):
    criteria = {
        'pv':
        [
            'PVS1'
        ],
        'ps':
        [
            'PS1',
            'PS2',
            'PS4',
            'PP1_strong'
        ],
        'pm':
        [
            'PM1',
            'PM2',
            'PM3',
            'PM4',
            'PM5',
            'PM6',
            'PP1_moderate'
        ],
        'pp':
        [
            'PP1',
            'PP2',
            'PP3',
            'PP4',
            'PP5'
        ],
        'ba':
        [
            'BA1'
        ],
        'bs':
        [
            'BS1',
            'BS2',
            'BS3',
            'BS4'
        ],
        'bp':
        [
            'BP1',
            'BP2',
            'BP3',
            'BP4',
            'BP5',
            'BP6',
            'BP7'
        ],
    }
    def __init__(self):
        for k in ACMG_criteria.criteria:
            for c in ACMG_criteria.criteria[k]:
                self[c] = False
    
    def combine(self):
        pv = sum(self['pv'])
        ps = sum(self['ps'])
        pm = sum(self['pm'])
        pp = sum(self['pp'])
        ba = sum(self['ba'])
        bs = sum(self['bs'])
        bp = sum(self['bp'])

        if (pv >= 1):
            if (ps >= 1):
                return 'Pathogenic'
            if (pm + pp >= 2):
                return 'Pathogenic'
            if (pm == 1):
                return 'Likely Pathogenic'
        if (ps >= 2):
            return 'Pathogenic'
        if (ps == 1):
            if (pm + (pp/2) >= 3):
                return 'Pathogenic'
            if (pm >= 1):
                return 'Likely Pathogenic'
            if (pp >= 2):
                return 'Likely Pathogenic'
        if (pm + (pp/2) >= 3):
            return 'Likely Pathogenic'
        
        if (ba >= 1):
            return 'Benign'
        if (bs >= 2):
            return 'Benign'
        if (bp >= 1):
            if (bs >= 1):
                return 'Likely Benign'
            if (bp >= 2):
                return 'Likely Benign'
        
        return 'VUS'

def process_bed(clinvar):
    cv = pd.read_csv(clinvar, sep = '\t')
    cv = cv[['#Chr', 'Start', 'End']]
    cv['#Chr'] = 'chr' + cv['#Chr'].astype(str)
    cv['Start'] = cv['Start'].astype(int)
    cv['End'] = cv['End'].astype(int)
    cv = cv.dropna()
    cv = cv.drop_duplicates()
    cv['Start'] = cv['Start'] - 2
    cv['End'] = cv['End'] + 3
    bed = clinvar + '.bed'
    cv.to_csv(bed, sep = '\t', index = False)
    return bed

def sort_and_merge(bedtools, bed):
    sorting = subprocess.Popen(
        [bedtools, "sort", "-i", bed],
        stdout = subprocess.PIPE
    )
    processed_bed = subprocess.check_output(
        [bedtools, "merge", "-i", "-"],
        stdin = sorting.stdout
    )
    sorting.wait()
    with open(bed, 'wb') as out:
        out.write(processed_bed)

def bed_intersect(bedtools, bed, temp):
    intersect_df = pd.DataFrame()

    intersection = subprocess.Popen(
        [bedtools, "intersect", "-a", bed, "-b", temp],
        stdout = subprocess.PIPE
    )

    temp_str = StringIO(intersection.communicate()[0].decode('utf-8'))

    intersect_df = pd.read_csv(temp_str, sep = '\t', header = None)

    return intersect_df

def add_columns(df):
    for k in ACMG_criteria.criteria:
        for c in ACMG_criteria.criteria[k]:
            df[c] = "."
    return df

def quest_pvs1_a(conseq):
    if not ('nonframeshift') in conseq:
        if 'frameshift' in conseq:
            return 1
    if 'splicing' in conseq:
        return 1
    if 'startloss' in conseq:
        return 1
    if 'stoploss' in conseq:
        return 1
    if 'stopgain' in conseq:
        return 1
    return 0

def quest_pvs1_b(conseq):
    if ('HIGH') in conseq:
        return 1
    return 0

def quest_pm1(domain):
    if domain != '.':
        return 1
    return 0

def quest_pm2(paf):
    if paf <= 0.005:
        return 1
    return 0

def quest_pm4(conseq):
    if 'inframe_insertion' in conseq:
        return 1
    if 'inframe_deletion' in conseq:
        return 1
    if 'stop_lost' in conseq:
        return 1
    return 0

def quest_pm5(conseq):
    if 'missense' in conseq:
        return 1
    return 0

def add_mis_z(row, mis_z_dict):
    if ('missense' in row['INFO_VEP_Consequence']):
        transcripts = row['INFO_VEP_Feature']
        z_scores = list()
        for t in transcripts.split(','):
            if t in mis_z_dict:
                z_scores.append(mis_z_dict[t])
        return ','.join([str(i) for i in sorted(list(set(z_scores)))])

def add_pLI(row, pLI_dict):
    transcripts = row['INFO_VEP_Feature']
    result = list()
    for t in transcripts.split(','):
        if t in pLI_dict:
            result.append(pLI_dict[t])
    return ','.join([str(i) for i in sorted(list(set(result)))])

def quest_pp2(mis_z):
    if not (mis_z is None):
        mis_z = [float(i) for i in mis_z.split(',') if i != ""]
        if len(mis_z) > 0:
            mis_z = max(mis_z)
            if float(mis_z) >= 3:
                return 1
    return 0

def quest_pp3(insilico):
    if insilico == 'Поврежд':
        return 1
    return 0

def quest_pp5(significance):
    if significance != '.':
        if not ('benign' in significance.lower()):
            return 1
    return 0

def classify(df):
    df = add_columns(df)
    df['PVS1'] = 0
    df['PVS1'] = df['Кодирующее_последствие'].map(quest_pvs1_a)
    df['PVS1'] = df['INFO_VEP_IMPACT'].map(quest_pvs1_b)
    print('PVS ready')
    bed = process_bed(clinvar)
    sort_and_merge(bedtools, bed)
    cv = pd.read_csv(bed, header = None, sep = '\t')
    cv[0] = cv[0].astype(str)
    coord = df[df['Кодирующее_последствие'] != '.'][['Хромосома', 'Позиция']]
    coord['End'] = coord['Позиция'] + 1
    coord.to_csv(clinvar + '.temp', index = False, header = False, sep = '\t')
    intersect_df = bed_intersect(bedtools, bed, clinvar + '.temp')
    intersect_df.columns = ['Хромосома', 'Позиция', 'Выбросить']
    ids = list(intersect_df.merge(df, on = ['Хромосома', 'Позиция'], how = 'left')['ID'].unique())
    df['PS1'] = 0
    df.loc[df['ID'].isin(ids), 'PS1'] = 1
    os.remove(clinvar + '.temp')
    os.remove(clinvar + '.bed')
    df['PS2'] = '?'
    df['PS3'] = '?'
    df['PS4'] = '?'
    print('PS ready')
    df['PM1'] = 0
    df['PM1'] = df['INFO_VEP_DOMAINS'].map(quest_pm1)
    df['PM2'] = 0
    df['PM2'] = df['Макс_попул_ч-та'].map(quest_pm2)
    df['PM3'] = '?'
    df['PM4'] = 0
    df['PM4'] = df['INFO_VEP_Consequence'].map(quest_pm4)
    df['PM5'] = 0
    df['PM5'] = df['INFO_VEP_Consequence'].map(quest_pm5)
    df['PM6'] = '?'
    print('PM ready')
    df['PP1'] = '?'
    cdf = pd.read_csv(constraint, sep = '\t')
    cdf = cdf[['transcript', 'mis_z', 'pLI']]
    mis_z_dict = dict(zip(cdf['transcript'], cdf['mis_z']))
    df['Missense_Z'] = df.apply(add_mis_z, args = (mis_z_dict, ), axis = 1)
    pLI_dict = dict(zip(cdf['transcript'], cdf['pLI']))
    df['pLI'] = df.apply(add_pLI, args = (pLI_dict, ), axis = 1)
    df['PP2'] = 0
    df['PP2'] = df['Missense_Z'].map(quest_pp2)
    df['PP3'] = 0
    df['PP3'] = df['In_silico_прогноз'].map(quest_pp3)
    df['PP4'] = '?'
    df['PP5'] = 0
    df['PP5'] = df['Клин_знач'].map(quest_pp5)
    print('PP ready')
    df['Баллы'] = df['PVS1']*4 \
    + df['PS1']*3 \
    + (df['PM1'] + df['PM2'] + df['PM4'] + df['PM5'])*2 \
    + df['PP2'] + df['PP3'] + df['PP5']
    return df
