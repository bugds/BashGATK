import os
import glob
import subprocess
import pandas as pd

reference = os.environ['reference']
maskGenes = os.environ['maskGenes']
bedtools = os.environ['bedtools']
wd = os.environ['outputFolder']

print('Make sure all your files have ".txt" extension')
print('Bed-files will be generated here:', maskGenes)

def get_ref_gene(annotation):
    return annotation.split('gene_name=')[1].split(';')[0]

def read_refGene(referencePath):
    return pd.read_csv(referencePath, header = None, sep = '\t')

def read_mask_files():
    mask_dict = dict()
    for filename in glob.glob(os.path.join(maskGenes, '*.txt')):
        path = os.path.join(maskGenes, filename)
        name = os.path.basename(path)
        with open(path, 'r') as inp:
            mask_dict[name] = inp.read().split('\n')
    return mask_dict

def get_masks(ref_df, mask_dict, maskGenes):
    ref_df['gene'] = ref_df[8].map(get_ref_gene)
    bed_names = list()
    for k in mask_dict:
        curr_df = pd.DataFrame()
        for g in mask_dict[k]:
            if g != '':
                mask_df = ref_df[ref_df['gene'] == g]
                if len(mask_df) == 0:
                    print('Please, replace', g, 'in', k)
                curr_df = pd.concat([curr_df, mask_df[[0, 3, 4]]])
        bed_name = os.path.join(maskGenes, k[:-3] + 'bed')
        bed_names.append(bed_name)
        curr_df.to_csv(bed_name, sep = '\t', header = None, index = False)
    return bed_names

def sort_and_merge(bedtools, bed_names):
    for bed in bed_names:
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

def filter_variants(wd, maskGenes):
    var_df = pd.read_csv(os.path.join(wd, 'rus2.tsv'), sep = '\t')
    var_df['Маска'] = ""
    for bed in sorted(glob.glob(os.path.join(maskGenes, '*.bed'))):
        path = os.path.join(maskGenes, bed)
        name = os.path.basename(path)[:-4] + ';'
        print("Creating mask", name[:-1])
        with open(bed, 'r') as inp:
            bed = [l.strip().split('\t') for l in inp.readlines()]
        for gene_string in bed:
            chr = gene_string[0]
            start = int(gene_string[1])
            end = int(gene_string[2])
            var_df.loc[
                (var_df['Хромосома'] == chr)
                & (var_df['Позиция'] <= end)
                & (var_df['Позиция'] >= start),
            'Маска'] += name
    var_df['Маска'] = var_df['Маска'].replace("", ".")
    return var_df

def save_sample_df(var_df):
    for sample in var_df['Проба'].unique():
        print("Saving", sample)
        try:
            with pd.ExcelWriter(os.path.join(wd, sample + '.xlsx')) as writer:
                sample_df = var_df[var_df['Проба'] == sample]
                sample_df[
                    (sample_df['Доверие'] == 'Выс') \
                    & (sample_df['Добро'] == '.') \
                    & (sample_df['Маска'] != '.') \
                    & (sample_df['Маска'] != 'CE;')
                ].to_excel(writer, sheet_name = 'Haematology', index = None)
                sample_df[
                    (sample_df['Доверие'] == 'Выс') \
                    & (sample_df['Добро'] == '.') \
                    & (sample_df['Маска'] != '.')
                ].to_excel(writer, sheet_name = 'Haematology+CE', index = None)
                sample_df.to_excel(writer, sheet_name = 'all', index = None)
            print("Excel saved")
        except:
            var_df[var_df['Проба'] == sample].to_csv(os.path.join(wd, sample + '.tsv'), sep = '\t', index = None)
            print("Excel saving failed, saved as tsv")

def main():
    ref_df = read_refGene(reference)
    mask_dict = read_mask_files()
    bed_names = get_masks(ref_df, mask_dict, maskGenes)
    sort_and_merge(bedtools, bed_names)
    var_df = filter_variants(wd, maskGenes)
    save_sample_df(var_df)

main()