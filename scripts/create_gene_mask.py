import os
import glob
import subprocess
import pandas as pd

reference = os.environ['reference']
maskGenes = os.environ['maskGenes']
bedtools = os.environ['bedtools']
generateBed = os.environ['generateBed']
wd = os.environ['outputFolder']

print('Make sure all your files have ".txt" extension')
print('Bed-files will be generated here:', maskGenes)

def get_ref_gene(annotation):
    return annotation.split('gene_name=')[1].split(';')[0]

def read_refGene(referencePath):
    return pd.read_csv(referencePath, header = None, sep = '\t')

def read_mask_file():
    mask_dict = dict()
    for filename in glob.glob(os.path.join(maskGenes, '*.txt')):
        path = os.path.join(maskGenes, filename)
        name = os.path.basename(path)
        with open(path, 'r') as inp:
            mask_dict[name] = inp.read().split('\n')
    return mask_dict

def get_masks(ref_df, mask_dict, maskGenes):
    bed_names = list()
    for k in mask_dict:
        curr_df = pd.DataFrame()
        for g in mask_dict[k]:
            if g != '':
                ref_df['gene'] = ref_df[8].map(get_ref_gene)
                mask_df = ref_df[ref_df['gene'] == g]
                if len(mask_df) == 0:
                    print('Please, replace', g, 'in', k)
                curr_df = curr_df.append(mask_df[[0, 3, 4]])
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
        print(processed_bed)
        with open(bed, 'wb') as out:
            out.write(processed_bed)

def filter_variants(wd, maskGenes):
    var_df = pd.read_csv(wd + '/rus2.tsv', sep = '\t')
    var_df_dict = dict()
    for bed in glob.glob(os.path.join(maskGenes, '*.bed')):
        path = os.path.join(maskGenes, bed)
        name = os.path.basename(path)[:-4]
        var_df_dict[name] = pd.DataFrame()
        with open(bed, 'r') as inp:
            bed = [l.strip().split('\t') for l in inp.readlines()]
        for gene_string in bed:
            chr = gene_string[0]
            start = int(gene_string[1])
            end = int(gene_string[2])
            var_df_dict[name] = var_df_dict[name].append(
                var_df[
                    (var_df['Хромосома'] == chr)
                    & (var_df['Позиция'] <= end)
                    & (var_df['Позиция'] >= start)
                ]
            )
    for k in var_df_dict:
        var_df_dict[k].to_csv(os.path.join(wd, k + '.tsv'), index = False, sep = '\t')

def main():
    if generateBed == 'yes':
        ref_df = read_refGene(reference)
        mask_dict = read_mask_file()
        bed_names = get_masks(ref_df, mask_dict, maskGenes)
        sort_and_merge(bedtools, bed_names)
    filter_variants(wd, maskGenes)

main()