import os
import glob
import subprocess
import pandas as pd

wd = os.environ['outputFolder']
ref = os.environ['refFasta']
anno = os.environ['annotation']
spai_scores = os.environ['knownScores']

if __name__ == "__main__":
    df = pd.read_csv(os.path.join(wd, 'xl_results', 'results.tsv'), sep = '\t')

    samples = list(df['Проба'].unique())
    samples = [i + '.xlsx' for i in samples]

    df = df[['Хромосома', 'Позиция', 'ID', 'Реф_аллель', 'Альт_аллель']]
    df = df.drop_duplicates()
    df.columns = ['#CHROM', 'POS', 'ID', 'REF', 'ALT']
    df['QUAL'] = '.'
    df['FILTER'] = 'PASS'
    df['INFO'] = '.'

    known_scores = pd.read_csv(spai_scores, sep = '\t')
    known_ids = list(known_scores['ID'])
    from_before = df[df['ID'].isin(known_ids)]
    df = df[~df['ID'].isin(known_ids)]

    df.to_csv(os.path.join(wd, 'xl_results', 'CE.vcf'), index = False, sep = '\t')

    if os.path.exists(os.path.join(wd, 'deepvariant')):
        for filename in glob.glob(os.path.join(wd, 'deepvariant', '*.vcf')):
            with open(filename, 'r') as inp:
                vcflines = [i for i in inp.readlines()if i.startswith('##')]
            break
    elif os.path.exists(os.path.join(wd, 'mutect2')):
        for filename in glob.glob(os.path.join(wd, 'mutect2', '*', '*.vcf')):
            with open(filename, 'r') as inp:
                vcflines = [i for i in inp.readlines()if i.startswith('##')]
            break
    elif os.path.exists(os.path.join(wd, 'haplotype_caller')):
        for filename in glob.glob(os.path.join(wd, 'haplotype_caller', '*.vcf')):
            with open(filename, 'r') as inp:
                vcflines = [i for i in inp.readlines()if i.startswith('##')]
            break
    else:
        raise Exception('No VCF files!!!')

    with open(os.path.join(wd, 'xl_results', 'CE.vcf'), 'r') as inp:
        lines = inp.readlines()

    with open(os.path.join(wd, 'xl_results', 'CE.vcf'), 'w') as out:
        out.write(''.join(vcflines))

    with open(os.path.join(wd, 'xl_results', 'CE.vcf'), 'a') as out:
        out.write(''.join(lines))

    subprocess.run(
        [
            "spliceai",
            "-I", os.path.join(wd, 'xl_results', 'CE.vcf'),
            "-O", os.path.join(wd, 'xl_results', 'CE_spai.txt'),
            "-R", ref,
            "-A", anno
        ]
    )

    lines_list = list()
    for filename in ['CE_spai.txt', spai_scores]:
        with open(os.path.join(wd, 'xl_results', filename), 'r') as inp:
            lines_list.append([i for i in inp.readlines() if (not (i.startswith('#')))])

    lines = lines_list[0] + lines_list[1]

    with open(spai_scores, 'w') as out:
        out.write('\t'.join(df.columns) + '\n')
        score_dict = dict()
        for i in lines:
            out.write(i)
            if 'SpliceAI=' in i:
                id = i.split('\t')[2]
                scores = i.split('SpliceAI=')[1].split('|')[2:6]
                scores = [float(i) for i in scores if i != '.']
                if scores:
                    score_dict[id] = str(max(scores))

    for s in samples:
        df_map = pd.read_excel(os.path.join(wd, 'xl_results', s), sheet_name=None)
        with pd.ExcelWriter(
            os.path.join(wd, 'xl_results', s),
            mode='w'
        ) as writer:
            for k in df_map:
                df_map[k]['SpliceAI'] = df_map[k]['ID'].map(score_dict)
                df_map[k]['SpliceAI'] = df_map[k]['SpliceAI'].fillna('.')
                df_map[k].to_excel(writer, sheet_name = k, index = None)
