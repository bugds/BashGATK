import os
import sys
import pandas as pd

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

def add_columns(df):
    for k in ACMG_criteria.criteria:
        for c in ACMG_criteria.criteria[k]:
            df[c] = False

def quest_pvs1(conseq):
    if not ('nonframeshift') in conseq:
        if 'frameshift' in conseq:
            return True
    if 'splicing' in conseq:
        return True
    if 'startloss' in conseq:
        return True
    if 'stoploss' in conseq:
        return True
    if 'stopgain' in conseq:
        return True

def quest_pm2(paf):
    if paf == '.':
        return True

def quest_pp3(insilico):
    if insilico == 'Поврежд':
        return True

def quest_ba1(paf):
    if paf != '.':
        if float(paf) > 0.05:
            return True

def quest_bp4(insilico):
    if insilico == 'Нейтрал':
        return True

def main():
    wd = os.path.abspath(sys.argv[1])
    df = pd.read_csv(wd + '/rus2.tsv', sep = '\t')
    #add_columns(df)
    df['PVS1'] = df['Кодирующее_последствие'].map(quest_pvs1)
    df['PM2'] = df['Макс_попул_ч-та'].map(quest_pm2)
    df['PP3'] = df['In_silico_прогноз'].map(quest_pp3)
    df['BA1'] = df['Макс_попул_ч-та'].map(quest_ba1)
    df['BP4'] = df['In_silico_прогноз'].map(quest_bp4)
    df.to_csv(wd + '/acmg.tsv', sep = '\t', index = False)
main()