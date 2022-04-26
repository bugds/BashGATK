import os
import re
import sys
import pandas as pd
#import matplotlib.pyplot as plt

depth_limit = 100
af_limit = 0.02
paf_limit = 0.01

wd = os.path.abspath(sys.argv[1])

rusDict = {
    'SAMPLE': 'Проба',
    '#CHROM': 'Хромосома',
    'POS': 'Позиция',
    'ID': 'ID варианта',
    'REF': 'Реф_аллель',
    'ALT': 'Альт_аллель',
    'QUAL': 'Кач-во',
    'FILTER': 'Уровень кач-ва',
    'FORMAT_AD': 'Глуб_прочтения_аллелей',
    'FORMAT_GT': 'Генотип',
    'FORMAT_VF': 'Альт_аллел_ч-та',
    'FORMAT_AF': 'Аллельная_частота',
    'FORMAT_VAF': 'Аллельная_частота',
    'FORMAT_DP': 'Глубина_прочтения',
    'INFO_ANNO_GENE': 'Q:_Ген',
    'INFO_ANNO_GENEINFO': 'Q:_ID_гена',
    'INFO_ANNO_ANN': 'Q:_Аннотация',
    'INFO_ANNO_CDS': 'Q:_Поз_в_кодирующ_посл-ти',
    'INFO_ANNO_AA': 'Q:_Белк_поз',
    'INFO_ANNO_STRAND': 'Q:_Направление',
    'INFO_ANNO_CAF': 'Q:_Аллел_ч-ты',
    'INFO_ANNO_DP': 'Q:_Глубина_прочт',
    'INFO_ANNO_UMT': 'Q:_Глубина_(молекуляр_тэгов)',
    'INFO_ANNO_VMF': 'Q:_Доля_вариантн_молекул_тэгов',
    'INFO_ANNO_VMT': 'Q:_Число_вариантн_молекул_тэгов',
    'INFO_ANNO_TYPE': 'Q:_Тип_варианта',
    'INFO_ANNO_VC': 'Q:_Класс_варианта',
    'INFO_ANNO_LOF': 'Q:_Loss-of-function',
    'INFO_ANNO_NMD': 'Q:_Нонсенс-опосред_распад',
    'INFO_ANNO_NSF': 'Q:_Несиноним_фреймшифт',
    'INFO_ANNO_NSN': 'Q:_Несиноним_нонсенс',
    'INFO_ANNO_NSM': 'Q:_Несиноним_миссенс',
    'INFO_ANNO_ASS': 'Q:_Акц_сплайс',
    'INFO_ANNO_DSS': 'Q:_Донор_сплайс',
    'INFO_ANNO_R3': 'Q:_В_3’-регионе',
    'INFO_ANNO_R5': 'Q:_В_5’-регионе',
    'INFO_ANNO_U3': 'Q:_3’-UTR_вариант',
    'INFO_ANNO_U5': 'Q:_5’-UTR_вариант',
    'INFO_ANNO_INT': 'Q:_Интронн',
    'INFO_ANNO_SYN': 'Q:_Синонимич',
    'INFO_ANNO_RepRegion': 'Q:_Регион_с_повторами',
    'INFO_ANNO_KGPhase1': 'Q:_1000_Геномов_фаза_1',
    'INFO_ANNO_KGPhase3': 'Q:_1000_Геномов_фаза_3',
    'INFO_ANNO_G5': 'Q:_>5%_в_одной попул',
    'INFO_ANNO_G5A': 'Q:_>5%_во_всех_попул',
    'INFO_ANNO_COMMON': 'Q:_Полиморфизм(1%)',
    'INFO_ANNO_PM': 'Q:_Цитируемый_вар-т',
    'INFO_ANNO_OM': 'Q:_Есть_OMIM',
    'INFO_ANNO_PMC': 'Q:_Статья_в_PMC',
    'INFO_ANNO_MUT': 'Q:_Доказ-я_мутация',
    'INFO_ANNO_RS': 'Q:_dbSNP_ID',
    'INFO_ANNO_RSPOS': 'Q:_dbSNP_позиция',
    'INFO_ANNO_SAO': 'Q:_Происхождение_варианта',
    'INFO_ANNO_S3D': 'Q:_Есть_3D_структура',
    'INFO_ANNO_VLD': 'Q:_Встречался_ранее',
    'INFO_ANNO_HD': 'Q:_Ранее_обнар-н_с_высок_ур-м_прочтения',
    'INFO_ANNO_WGT': 'Q:_Вес_(число_встреч)',
    'INFO_ANNO_OTH': 'Q:_Есть_иной_вар-т_в_этих_позициях',
    'INFO_ANNO_SLO': 'Q:_SubmitterLinkOut',
    'INFO_ANNO_dbSNPBuildID': 'Q:_dbSNP_ID_версии',
    'INFO_ANNO_VP': 'Q:_dbSNP_bitField_код',
    'INFO_ANNO_WTD': 'Q:_dbSNP_отклонён',
    'INFO_ANNO_SSR': 'Q:_Причина_подозрит_вар-та',
    'INFO_ANNO_CLNHGVS': 'Q:_HGVS_формат',
    'INFO_ANNO_CLNACC': 'Q:_ClinVar_ID',
    'INFO_ANNO_CLNALLE': 'Q:_ClinVar_номер_аллеля',
    'INFO_ANNO_CLNDBN': 'Q:_ClinVar_заболевание',
    'INFO_ANNO_CLNDSDB': 'Q:_ClinVar_заболевание_в_базе',
    'INFO_ANNO_CLNDSDBID': 'Q:_ClinVar_ID_заболевания',
    'INFO_ANNO_CLNORIGIN': 'Q:_ClinVar_источник',
    'INFO_ANNO_CLNREVSTAT': 'Q:_ClinVar_статус',
    'INFO_ANNO_CLNSIG': 'Q:_ClinVar_значение',
    'INFO_ANNO_CLNSRC': 'Q:_ClinVar_организация',
    'INFO_ANNO_CLNSRCID': 'Q:_ClinVar_ID_организации',
    'INFO_ANNO_CDA': 'Q:_Клин_иссл',
    'INFO_ANNO_CNT': 'Q:_ClinVar_Встреч-ть',
    'INFO_ANNO_GNO': 'Q:_ClinVar_наличие_генотипа',
    'INFO_ANNO_ASP': 'Q:_Един_референс',
    'INFO_ANNO_CFL': 'Q:_Реф_конфликт',
    'INFO_ANNO_REF': 'Q:_Один_из_аллелей_совпад_с_реф',
    'INFO_ANNO_NOC': 'Q:_Реф_не_среди_аллелей',
    'INFO_ANNO_NOV': 'Q:_Непересек-ся_аллели',
    'INFO_ANNO_LSD': 'Q:_Из_локус_специф_БД',
    'INFO_ANNO_MTP': 'Q:_Доп_аннотац',
    'INFO_ANNO_TPA': 'Q:_Аннотац_с_фенотипом',
    'INFO_ANNO_RV': 'Q:_dbSNP_реверсн_ориентац',
    'INFO_ANNO_ReverseComplementedAlleles': 'Q:_Реверс_комплемент_между_реф_и_аннотац',
    'INFO_ANNO_SwappedAlleles': 'Q:_Реф_и_альт_аллели_были_спутаны',
    'INFO_ANNO_ANNOVAR_DATE': 'A:_Дата_выпуска_базы',
    'INFO_ANNO_Func.refGene': 'A:_Функц_послед_RG',
    'INFO_ANNO_Gene.refGene': 'A:_Ген_RG',
    'INFO_ANNO_GeneDetail.refGene': 'A:_Детально_ген_RG',
    'INFO_ANNO_ExonicFunc.refGene': 'A:_Экзон_функц_RG',
    'INFO_ANNO_AAChange.refGene': 'A:_Смена_аминок-т_RG',
    'INFO_ANNO_Func.ensGene': 'A:_Функц_послед_EG',
    'INFO_ANNO_Gene.ensGene': 'A:_Ген_EG',
    'INFO_ANNO_GeneDetail.ensGene': 'A:_Детально_ген_EG',
    'INFO_ANNO_ExonicFunc.ensGene': 'A:_Экзон_функц_EG',
    'INFO_ANNO_AAChange.ensGene': 'A:_Смена_аминок-т_EG',
    'INFO_ANNO_AF': 'A:_Попул_ч-та',
    'INFO_ANNO_AF_popmax': 'A:_Макс_попул_ч-та',
    'INFO_ANNO_AF_male': 'A:_Муж_попул_ч-та',
    'INFO_ANNO_AF_female': 'A:_Жен_попул_ч-та',
    'INFO_ANNO_AF_raw': 'A:_Попул_ч-та_до_оц_кач-ва',
    'INFO_ANNO_AF_afr': 'A:_Попул_ч-та_Афр',
    'INFO_ANNO_AF_sas': 'A:_Попул_ч-та_Юазия',
    'INFO_ANNO_AF_amr': 'A:_Попул_ч-та_Амер',
    'INFO_ANNO_AF_eas': 'A:_Попул_ч-та_Вазия',
    'INFO_ANNO_AF_nfe': 'A:_Попул_ч-та_Евр_не_Фин',
    'INFO_ANNO_AF_fin': 'A:_Попул_ч-та_Фин',
    'INFO_ANNO_AF_asj': 'A:_Попул_ч-та_Ашк',
    'INFO_ANNO_AF_oth': 'A:_Попул_ч-та_другие',
    'INFO_ANNO_non_topmed_AF_popmax': 'A:_Попул_ч-та_не_TOPMed',
    'INFO_ANNO_non_neuro_AF_popmax': 'A:_Попул_ч-та_не_неврол',
    'INFO_ANNO_non_cancer_AF_popmax': 'A:_Попул_ч-та_не_онко',
    'INFO_ANNO_controls_AF_popmax': 'A:_Попул_ч-та_контрольн_группа',
    'INFO_ANNO_AF_ami': 'A:_Попул_ч-та_Амиши',
    'INFO_ANNO_DamagePredCount': 'A:_Предсказ_патогенности',
    'INFO_ANNO_SIFT_pred': 'A:_Предсказ_SIFT',
    'INFO_ANNO_SIFT4G_pred': 'A:_Предсказ_SIFT4G',
    'INFO_ANNO_Polyphen2_HDIV_pred': 'A:_Предсказ_PolyPhen2_HDIV',
    'INFO_ANNO_Polyphen2_HVAR_pred': 'A:_Предсказ_PolyPhen2_HVAR',
    'INFO_ANNO_LRT_pred': 'A:_Предсказ_LRT',
    'INFO_ANNO_MutationTaster_pred': 'A:_Предсказ_MutationTaster',
    'INFO_ANNO_MutationAssessor_pred': 'A:_Предсказ_MutationAssessor',
    'INFO_ANNO_FATHMM_pred': 'A:_Предсказ_FATHMM',
    'INFO_ANNO_PROVEAN_pred': 'A:_Предсказ_PROVEAN',
    'INFO_ANNO_VEST4_score': 'A:_Предсказ_VEST4',
    'INFO_ANNO_MetaSVM_pred': 'A:_Предсказ_MetaSVM',
    'INFO_ANNO_MetaLR_pred': 'A:_Предсказ_MetaLR',
    'INFO_ANNO_M-CAP_pred': 'A:_Предсказ_M-CAP',
    'INFO_ANNO_REVEL_score': 'A:_Предсказ_REVEL',
    'INFO_ANNO_MutPred_score': 'A:_Предсказ_MutPred',
    'INFO_ANNO_MVP_score': 'A:_Предсказ_MVP',
    'INFO_ANNO_MPC_score': 'A:_Предсказ_MPC',
    'INFO_ANNO_PrimateAI_pred': 'A:_Предсказ_PrimateAI',
    'INFO_ANNO_DEOGEN2_pred': 'A:_Предсказ_DEOGEN2',
    'INFO_ANNO_BayesDel_addAF_pred': 'A:_Предсказ_BayesDel_addAF',
    'INFO_ANNO_BayesDel_noAF_pred': 'A:_Предсказ_BayesDel_noAF',
    'INFO_ANNO_ClinPred_pred': 'A:_Предсказ_ClinPred',
    'INFO_ANNO_LIST-S2_pred': 'A:_Предсказ_LIST-S2',
    'INFO_ANNO_CADD_raw': 'A:_Предсказ_CADD_raw',
    'INFO_ANNO_CADD_phred': 'A:_Предсказ_CADD_phred',
    'INFO_ANNO_DANN_score': 'A:_Предсказ_DANN',
    'INFO_ANNO_fathmm-MKL_coding_pred': 'A:_Предсказ_FATHMM-MKL',
    'INFO_ANNO_fathmm-XF_coding_pred': 'A:_Предсказ_FATHMM-XF',
    'INFO_ANNO_Eigen-raw_coding': 'A:_Предсказ_Eigen_raw',
    'INFO_ANNO_Eigen-phred_coding': 'A:_Предсказ_Eigen_phred',
    'INFO_ANNO_Eigen-PC-raw_coding': 'A:_Предсказ_Eigen_PC_raw',
    'INFO_ANNO_Eigen-PC-phred_coding': 'A:_Предсказ_Eigen_PC_phred',
    'INFO_ANNO_GenoCanyon_score': 'A:_Предсказ_GenoCanyon',
    'INFO_ANNO_integrated_fitCons_score': 'A:_Предсказ_fitCons_integrated',
    'INFO_ANNO_GM12878_fitCons_score': 'A:_Предсказ_fitCons_GM12878',
    'INFO_ANNO_H1-hESC_fitCons_score': 'A:_Предсказ_fitCons_H1-hESC',
    'INFO_ANNO_HUVEC_fitCons_score': 'A:_Предсказ_fitCons_HUVEC',
    'INFO_ANNO_LINSIGHT': 'A:_Предсказ_LINSIGHT',
    'INFO_ANNO_GERP++_NR': 'A:_Предсказ_GEPR++_NR',
    'INFO_ANNO_GERP++_RS': 'A:_Предсказ_GEPR++_RS',
    'INFO_ANNO_phyloP100way_vertebrate': 'A:_Предсказ_phyloP100way_vertebrate',
    'INFO_ANNO_phyloP30way_mammalian': 'A:_Предсказ_phyloP30way_mammalian',
    'INFO_ANNO_phyloP17way_primate': 'A:_Предсказ_phyloP17way_primate',
    'INFO_ANNO_phastCons100way_vertebrate': 'A:_Предсказ_phastCons100way_vertebrate',
    'INFO_ANNO_phastCons30way_mammalian': 'A:_Предсказ_phastCons30way_mammalian',
    'INFO_ANNO_phastCons17way_primate': 'A:_Предсказ_phastCons17way_primate',
    'INFO_ANNO_bStatistic': 'A:_Предсказ_bStatistic',
    'INFO_ANNO_Interpro_domain': 'A:_InterPro_домен',
    'INFO_ANNO_GTEx_V8_gene': 'A:_GTEx_ген',
    'INFO_ANNO_GTEx_V8_tissue': 'A:_GTEx_ткань_(экспрессия)',
    'INFO_ANNO_dbscSNV_ADA_SCORE': 'A:_dbscSNV_ADA_SCORE',
    'INFO_ANNO_dbscSNV_RF_SCORE': 'A:_dbscSNV_RF_SCORE',
    'INFO_ANNO_CLNALLELEID': 'A:_ClinVar_ID',
    'INFO_ANNO_CLNDN': 'A:_ClinVar_заболевание',
    'INFO_ANNO_CLNDISDB': 'A:_ClinVar_ID_заболевания',
    'INFO_ANNO_cosmic95_coding': 'A:_COSMIC92_кодир',
    'INFO_ANNO_cosmic95_noncoding': 'A:_COSMIC92_все',
    'INFO_ANNO_avsnp150': 'A:_avsnp150_ID',
    'INFO_VEP_Allele': 'V:_Считанный_аллель',
    'INFO_VEP_Consequence': 'V:_Последствия',
    'INFO_VEP_IMPACT': 'V:_Степень_влияния',
    'INFO_VEP_SYMBOL': 'V:_Наим_гена',
    'INFO_VEP_Gene': 'V:_ID_гена',
    'INFO_VEP_Feature_type': 'V:_Функция_продукта',
    'INFO_VEP_Feature': 'V:_ID_продукта',
    'INFO_VEP_BIOTYPE': 'V:_Тип_варианта',
    'INFO_VEP_EXON': 'V:_Экзон',
    'INFO_VEP_INTRON': 'V:_Интрон',
    'INFO_VEP_HGVSc': 'V:_HGVSc_формат',
    'INFO_VEP_HGVSp': 'V:_HGVSp_формат',
    'INFO_VEP_cDNA_position': 'V:_кДНК_поз',
    'INFO_VEP_CDS_position': 'V:_Поз_в_кодирующ_посл-ти',
    'INFO_VEP_Protein_position': 'V:_Белк_поз',
    'INFO_VEP_Amino_acids': 'V:_Аминок-ты',
    'INFO_VEP_Codons': 'V:_Кодоны',
    'INFO_VEP_Existing_variation': 'V:_Сущ_варианты',
    'INFO_VEP_DISTANCE': 'V:_Ск-ко_п.о._от_вар-та_до_транскрипта',
    'INFO_VEP_STRAND': 'V:_Направление',
    'INFO_VEP_FLAGS': 'V:_Флаги_кач-ва_транскрипта',
    'INFO_VEP_VARIANT_CLASS': 'V:_Класс_варианта',
    'INFO_VEP_SYMBOL_SOURCE': 'V:_Источник_наим_гена',
    'INFO_VEP_HGNC_ID': 'V:_HGNC_ID',
    'INFO_VEP_CANONICAL': 'V:_Канонич',
    'INFO_VEP_MANE_SELECT': 'V:_кДНК_ID',
    'INFO_VEP_MANE_PLUS_CLINICAL': 'V:_кДНК_ID+',
    'INFO_VEP_TSL': 'V:_Ур-нь_поддержки_транскр-та',
    'INFO_VEP_APPRIS': 'V:_Оц-ка_первичности_транскр-та',
    'INFO_VEP_CCDS': 'V:_CCDS_ID',
    'INFO_VEP_ENSP': 'V:_ENSP_ID',
    'INFO_VEP_SWISSPROT': 'V:_Swiss-Prot_ID',
    'INFO_VEP_TREMBL': 'V:_TrEMBL_ID',
    'INFO_VEP_UNIPARC': 'V:_UniParc_ID',
    'INFO_VEP_UNIPROT_ISOFORM': 'V:_UniProt_ID_изоформы',
    'INFO_VEP_GENE_PHENO': 'V:_Фенотип_гена',
    'INFO_VEP_SIFT': 'V:_Предсказ_SIFT',
    'INFO_VEP_PolyPhen': 'V:_Предсказ_PolyPhen',
    'INFO_VEP_DOMAINS': 'V:_Домен_белка',
    'INFO_VEP_miRNA': 'V:_miRNA',
    'INFO_VEP_HGVS_OFFSET': 'V:_Смещение_HGVS',
    'INFO_VEP_AF': 'V:_Попул_ч-та_1000_Геномов',
    'INFO_VEP_AFR_AF': 'V:_Попул_ч-та_Афр',
    'INFO_VEP_AMR_AF': 'V:_Попул_ч-та_Амер',
    'INFO_VEP_EAS_AF': 'V:_Попул_ч-та_Вазия',
    'INFO_VEP_EUR_AF': 'V:_Попул_ч-та_Евр',
    'INFO_VEP_SAS_AF': 'V:_Попул_ч-та_Юазия',
    'INFO_VEP_AA_AF': 'V:_Попул_ч-та_Афроамер',
    'INFO_VEP_EA_AF': 'V:_Попул_ч-та_Евроамер',
    'INFO_VEP_gnomAD_AF': 'V:_Попул_ч-та_gnomAD',
    'INFO_VEP_gnomAD_AFR_AF': 'V:_Попул_ч-та_gnomAD_Афр',
    'INFO_VEP_gnomAD_AMR_AF': 'V:_Попул_ч-та_gnomAD_Амер',
    'INFO_VEP_gnomAD_ASJ_AF': 'V:_Попул_ч-та_gnomAD_Ашк',
    'INFO_VEP_gnomAD_EAS_AF': 'V:_Попул_ч-та_gnomAD_Вазия',
    'INFO_VEP_gnomAD_FIN_AF': 'V:_Попул_ч-та_gnomAD_Фин',
    'INFO_VEP_gnomAD_NFE_AF': 'V:_Попул_ч-та_gnomAD_Евр_не_Фин',
    'INFO_VEP_gnomAD_OTH_AF': 'V:_Попул_ч-та_gnomAD_другие',
    'INFO_VEP_gnomAD_SAS_AF': 'V:_Попул_ч-та_gnomAD_Юазия',
    'INFO_VEP_MAX_AF': 'V:_Макс_попул_ч-та',
    'INFO_VEP_MAX_AF_POPS': 'V:_Популяц_с_макс_попул_ч-той',
    'INFO_VEP_CLIN_SIG': 'V:_ClinVar_значение',
    'INFO_VEP_SOMATIC': 'V:_Соматич_статус',
    'INFO_VEP_PHENO': 'V:_Фенотип_вар-та',
    'INFO_VEP_PUBMED': 'V:_Pubmed_ID',
    'INFO_VEP_MOTIF_NAME': 'V:_Профиль_связывания_с_ф-ми_транскрипц',
    'INFO_VEP_MOTIF_POS': 'V:_Позиц_вар-та_в_профиле_связ_с_ф-ми_транскр',
    'INFO_VEP_HIGH_INF_POS': 'V:_Вар-т_в_высокоинформативн_поз_профиля_связ_с_ф-ми_транскр',
    'INFO_VEP_MOTIF_SCORE_CHANGE': 'V:_Измен_профиля_связ_с_ф-ми_транскр',
    'INFO_VEP_TRANSCRIPTION_FACTORS': 'V:_Ф-ры_транскрипц'
}

for folder in ['_anno_soma', '_anno_germ', '']:
    for filename in ['/combined' + folder + '.csv', '/combined_passed_' + folder + '.csv']:
        if os.path.exists(filename):
            df = pd.read_csv(wd + filename, sep = '\t')
            short_dict = dict()
            for k in rusDict:
                if k in df.columns:
                    short_dict[k] = rusDict[k]
            df = df.rename(columns = short_dict)
            df = df[short_dict.values()]
            df.to_csv(wd + filename.replace('combined', 'rus_combined'), sep = '\t', index = False, encoding = 'utf-8')

    if os.path.exists(wd + '/combined_passed' + folder + '.csv'):
        report_df = pd.read_csv(wd + '/combined_passed' + folder + '.csv', sep = '\t')
        good_functions = ['exonic', 'splicing']
        # good_functions = '|'.join(good_functions)
        
        # report_df = report_df[report_df['INFO_ANNO_Func.refGene'].str.contains(good_functions)]
        report_df = report_df[report_df['INFO_ANNO_Func.refGene'].isin(good_functions)]
        if 'INFO_ANNO_UMT' in report_df.columns:
            depth = 'INFO_ANNO_UMT'
            af = 'INFO_ANNO_VMF'
        elif 'FORMAT_VAF' in report_df.columns:
            depth = 'FORMAT_DP'
            af = 'FORMAT_VAF'
            vaf_df = report_df[report_df['FORMAT_VAF'].str.contains(',')].reset_index(drop = True)
            alt_df = report_df[report_df['ALT'].str.contains(',')].reset_index(drop = True)
            if len(vaf_df) != len(alt_df):
                dif_df = pd.concat([vaf_df, alt_df]).drop_duplicates(keep = False)
                report_df = pd.concat([report_df, dif_df]).drop_duplicates(keep = False)
                for lst_col in ['FORMAT_VAF', 'ALT']:
                    report_df = report_df.assign(**{lst_col:report_df[lst_col].astype(str).str.split(',')})
                report_df = report_df.explode(['FORMAT_VAF', 'ALT'])
                report_df = pd.concat([report_df, dif_df])
            elif len(vaf_df) == 0:
                pass
            elif vaf_df != alt_df:
                dif_df = pd.concat([vaf_df, alt_df]).drop_duplicates(keep = False)
                report_df = pd.concat([report_df, dif_df]).drop_duplicates(keep = False)
                for lst_col in ['FORMAT_VAF', 'ALT']:
                    report_df = report_df.assign(**{lst_col:report_df[lst_col].astype(str).str.split(',')})
                report_df = report_df.explode(['FORMAT_VAF', 'ALT'])
                report_df = pd.concat([report_df, dif_df])
            else:
                for lst_col in ['FORMAT_VAF', 'ALT']:
                    report_df = report_df.assign(**{lst_col:report_df[lst_col].astype(str).str.split(',')})
                report_df = report_df.explode(['FORMAT_VAF', 'ALT'])
        elif 'FORMAT_AF' in report_df.columns:
            depth = 'FORMAT_DP'
            af = 'FORMAT_AF'
        report_df = report_df[report_df[depth] != '.']
        report_df = report_df[report_df[af] != '.']
        report_df = report_df[report_df[depth].astype(int) > depth_limit]
        report_df = report_df[report_df[af].astype(float) > af_limit]
        paf = 'INFO_ANNO_AF_popmax'
        numeric_df = report_df[pd.to_numeric(report_df[paf], errors = 'coerce').notnull()]
        non_numeric_df = pd.concat([numeric_df, report_df]).drop_duplicates(keep = False)
        numeric_df = numeric_df[numeric_df[paf].astype(float) < paf_limit]
        report_df = pd.concat([numeric_df, non_numeric_df])
        
        class ReportedMutation():
            def __init__(self, name, paf):
                self.name = name
                self.paf = paf
        
        report_dict = dict()
        for s in report_df['SAMPLE'].unique():
            sample_df = report_df[report_df['SAMPLE'] == s]
            report_dict[s] = dict()
            genes = sample_df['INFO_ANNO_Gene.refGene'].unique()
            for g in genes:
                genes_df = sample_df[sample_df['INFO_ANNO_Gene.refGene'] == g]
                report_dict[s][g] = list()
                for index, row in genes_df.iterrows():
                    report_dict[s][g].append(ReportedMutation(row['INFO_ANNO_AAChange.refGene'].split(',')[0].split(':')[-1], row[af]))
        
        #for s in report_dict:
        #    print(s)
        #    plt.figure(s)
        #    i = 0
        #    for g in report_dict[s]:
        #        i += 1
        #        ax = plt.subplot(len(report_dict[s]),1,i)
        #        fig = plt.gcf()
        #        fig.set_size_inches(10, 10)
        #        ax.bar([m.name for m in report_dict[s][g]], [m.paf for m in report_dict[s][g]])
        #        ax.set_yticklabels([0, 1])
        #        ax.set_ylim([-0.1, 1.1])
        #        ax.set_title(g)
        #    plt.tight_layout(pad = 1)
        #    plt.show()
        #    plt.close()
        #    break
        
        report_df = report_df.rename(columns = rusDict)
        rusDictValues = [v for v in report_df.columns if v in rusDict.values()]
        report_df = report_df[rusDictValues]
        report_df.to_csv(wd + '/report' + folder + '.csv', sep = '\t', index = False)
        
        report_df['the_id'] = report_df[rusDict['#CHROM']] \
            + ':' \
            + report_df[rusDict['POS']].astype(str) \
            + ':' \
            + report_df[rusDict['REF']] \
            + '/' \
            + report_df[rusDict['ALT']]
        
        new_columns = list(report_df[rusDict['SAMPLE']].unique())
        short_report_df = report_df[[rusDict['INFO_ANNO_Gene.refGene'], rusDict[paf], rusDict[af], rusDict['SAMPLE'], 'the_id']]
        #print(short_report_df.to_string())
        short_report_df = short_report_df.groupby('the_id').agg(list)
        
        for column in new_columns:
            short_report_df[column] = 0
        
        additional_new_columns = [
            rusDict['INFO_ANNO_Gene.refGene'],
            rusDict[paf]
        ]
        
        new_columns = additional_new_columns + new_columns
        
        for index, row in short_report_df.iterrows():
            i = 0
            for column in row[rusDict['SAMPLE']]:
                short_report_df.loc[index, column] = row[rusDict[af]][i]
                i += 1
            for column in additional_new_columns:
                if len(set(row[column])) == 1:
                    short_report_df.loc[index, column] = row[column][0]
                else:
                    raise Exception('BAD ADDITIONAL COLUMN')
        
        short_report_df = short_report_df[new_columns]
        
        short_report_df.to_csv(wd + '/short_report' + folder + '.csv', sep = '\t')
