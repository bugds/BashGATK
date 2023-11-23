import os

wd = os.environ['outputFolder']
ref_fasta_path = os.environ['refFasta']
R_script = os.environ['Rscript']
exome_table_file_path = os.environ['exomeTableFilePath']
gene_table_file_path = os.environ['geneTableFilePath']
dgv_merged = os.environ['DGVMerged']

bam_files = ['"' + bam + '"' for bam in list(os.listdir(wd + '/recalibrated')) if bam.endswith('.bam')]

tests = ','.join(bam_files)
tests = 'c(' + tests + ')'

with open(R_script, 'r') as inp:
    script = inp.read()

script = script.replace('curr.recalibrated.folder', '"' + wd + 'recalibrated"')
script = script.replace('exome_table_file_path', exome_table_file_path)
script = script.replace('gene_table_file_path', gene_table_file_path)
script = script.replace('DGV_merged', dgv_merged)
script = script.replace('bam_files', ','.join(bam_files))
script = script.replace('ref_fasta_path', ref_fasta_path)
script = script.replace('tests', tests)
script = script.replace('exome_calls_path', wd + 'CNV/')

os.makedirs(wd + '/CNV', exist_ok = True)

with open(wd + '/CNV/script.R', 'w') as out:
    out.write(script)