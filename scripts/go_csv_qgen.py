import os
import sys

wd = os.environ['outputFolder']

global innerDelimeter
global outerDelimeter
global extensionFilter
global outputName

endingAnno = 'ALLELE_END'
innerDelimeter = ','
outerDelimeter = '\t'

def lineListCheck(someList, filename):
    if len(someList) == 1:
        if 'Format:' in next(iter(someList)):
            someLine = next(iter(someList)).split('Format: ')[1].replace('">', '')
            return someLine.split('|')
        else:
            someLine = next(iter(someList))
            return someLine.split('\t')
    else:
        raise Exception('VEP or header line not identifiable in file ' + filename)

def dataSizeCheck(varData, headerList):
    if not len(varData) == len(headerList):
        raise Exception('Header size different')

def getData(wd, filename):
    with open(wd + '/annotation/' + filename, 'r') as vcfFile:
        lines = vcfFile.readlines()

        preInfoLines = '\n'.join(lines).split('##reference=')[-1]
        infoLines = [l.replace('\n', '') for l in preInfoLines.split('\n') if l.startswith('##INFO=')]
        formatLines = [l.replace('\n', '') for l in lines if l.startswith('##FORMAT=')]
        
        formatList = [l.split('##FORMAT=<ID=')[1].split(',')[0] for l in formatLines]
        
        annoList = [l.split('##INFO=<ID=')[1].split(',')[0] for l in infoLines][:-2]
        
        vepList = [l for l in infoLines if l.startswith('##INFO=<ID=CSQ')]
        vepList = lineListCheck(vepList, filename)
        
        headerList = [l.replace('\n', '') for l in lines if l.startswith('#CHROM')]
        headerList = lineListCheck(headerList, filename)
        headerList[-1] = 'SAMPLE_DATA'
        
        dataLines = [l.replace('\n', '') for l in lines if not l.startswith('#')]

        return formatList, annoList, vepList, headerList, dataLines

def parseFormatColumns(line, headerList, formatList, filename):
    varData = line.split('\t')
    dataSizeCheck(varData, headerList)
    varDict = dict()
    for i in zip(headerList, varData):
        varDict[i[0]] = i[1]
    
    for f in formatList:
        if f in varDict['FORMAT'].split(':'):
            index = varDict['FORMAT'].split(':').index(f)
            varDict['FORMAT_' + f] = varDict['SAMPLE_DATA'].split(':')[index]
        else:
            varDict['FORMAT_' + f] = '.'
    varDict.pop('FORMAT')
    varDict.pop('SAMPLE_DATA')
    varDict['SAMPLE'] = filename.split('.')[0]
    return varDict

def parseAnnoColumns(varDict, annoList, delimeter=innerDelimeter):
    annoData = varDict['INFO'].split(';')[:-1]
    annoChunk = annoData[:annoData.index(endingAnno)]
    for c in range(0, annoData.count(endingAnno)):
        for key in annoList:
            valueFound = False
            for element in annoChunk:
                if key == element.split('=')[0]:
                    valueFound = True
                    if c > 0:
                        try:
                            varDict['INFO_ANNO_' + key] += delimeter + element.split('=')[1]
                        except IndexError:
                            varDict['INFO_ANNO_' + key] += delimeter + element
                    else:
                        try:
                            varDict['INFO_ANNO_' + key] = element.split('=')[1]
                        except IndexError:
                            varDict['INFO_ANNO_' + key] = element
                    # break
            if not valueFound:
                varDict['INFO_ANNO_' + key] = '.'
        annoData = annoData[annoData.index(endingAnno) + 1:]
        try:
            annoChunk = annoData[:annoData.index(endingAnno)]
        except:
            annoChunk = []
    return varDict

def parseVepColumns(varDict, vepList, delimeter=innerDelimeter):
    vepData = varDict['INFO'].split(';')[-1].replace('CSQ=', '').split(',')
    varDict.pop('INFO')
    for vepChunk in vepData:
        vepChunk = vepChunk.split('|')
        for i in range(0, len(vepChunk)):
            if vepChunk[i] == '':
                vepChunk[i] = '.'
        for i in zip(vepList, vepChunk):
            if ('INFO_VEP_' + i[0]) in varDict:
                varDict['INFO_VEP_' + i[0]] += delimeter + i[1]
            else:
                varDict['INFO_VEP_' + i[0]] = i[1]
    return varDict

def writeToCsv(varDict, headerDone, csvFile, delimeter=outerDelimeter):
    if not headerDone:
        csvFile.write('\t'.join(k for k in varDict.keys()) + '\n')
        headerDone = True

    for i in varDict:
        setOfElements = set(varDict[i].split(','))
        setOfElements.discard('.')
        if len(setOfElements) == 1:
            varDict[i] = next(iter(setOfElements))
        elif len(setOfElements) == 0:
            varDict[i] = '.'
    csvFile.write(delimeter.join(v for v in varDict.values()) + '\n')
    return headerDone

def createCsv(extensionFilter, outputName):
    headerDone = False
    
    with open(wd + outputName, 'w') as csvFile:
        for filename in os.listdir(wd + '/annotation'):
                if filename.endswith(extensionFilter):
                    formatList, annoList, vepList, headerList, dataLines = getData(wd, filename)
                    for line in dataLines:
                        varDict = parseFormatColumns(line, headerList, formatList, filename)
                        varDict = parseAnnoColumns(varDict, annoList)
                        varDict = parseVepColumns(varDict, vepList)
                        
                        headerDone = writeToCsv(varDict, headerDone, csvFile)

def add_inhouse_freq():
    import pandas as pd

    freq_file = os.environ['freq_file']
    
    with open(wd + '/combined.csv', 'r') as inpObj:
        df = pd.read_csv(inpObj, sep='\t')
    
    df['temp_id'] = df['#CHROM'] + ':' + df['POS'].astype(str) + ':' + df['REF'] + '>' + df['ALT'] + '_' + 'smCounter2'
    if os.path.exists(freq_file):
        freq = pd.read_csv(freq_file, sep = '\t')
        freq = pd.concat([freq, df[['temp_id', 'SAMPLE']]])
    else:
        freq = df[['temp_id', 'SAMPLE']]
    freq = freq.drop_duplicates()
    freq.to_csv(freq_file, index = False, sep = '\t')
    counts = freq['temp_id'].value_counts()
    df['InHouse_calls'] = df['temp_id'].map(counts)
    sample_num = len(freq['SAMPLE'].unique())
    print('Total number of samples:', sample_num)
    df['InHouse_freqs'] = df['InHouse_calls'] / sample_num
    print('Added in-house frequency')

    tableP = pd.read_csv(wd + '/combined_passed.csv', sep = '\t')
    tableP['temp_id'] = tableP['#CHROM'] + ':' + tableP['POS'].astype(str) + ':' + tableP['REF'] + '>' + tableP['ALT'] + '_' + 'smCounter2'
    tableP['InHouse_calls'] = tableP['temp_id'].map(counts)
    tableP['InHouse_freqs'] = tableP['InHouse_calls'] / sample_num
    tableP = tableP.drop(columns=['temp_id'])
    tableP['ID'] = tableP['#CHROM'] + '-' + tableP['POS'].astype(str) + '-' + tableP['REF'] + '-' + tableP['ALT']
    tableP['INFO_ANNO_gnomad40_genome_AF'] = tableP['INFO_ANNO_gnomad40_genome_AF'].replace('.', -1).astype(float)
    tableP['INFO_ANNO_gnomad40_exome_AF'] = tableP['INFO_ANNO_gnomad40_exome_AF'].replace('.', -1).astype(float)
    tableP['popAF'] = tableP[['INFO_ANNO_gnomad40_genome_AF', 'INFO_ANNO_gnomad40_exome_AF']].max(axis=1)
    tableP['INFO_ANNO_gnomad40_genome_AF'] = tableP['INFO_ANNO_gnomad40_genome_AF'].replace(-1, '.')
    tableP['INFO_ANNO_gnomad40_exome_AF'] = tableP['INFO_ANNO_gnomad40_exome_AF'].replace(-1, '.')
    tableP['popAF'][tableP['popAF'] > 0.03] = 'common'
    tableP['popAF'][tableP['popAF'] != 'common'] = 'uncommon'
    tableP.to_csv(wd + '/combined_passed.csv', sep = '\t', index = False)

if __name__ == "__main__":
    createCsv('vep.vcf.pass.vcf', '/combined_passed.csv')
    createCsv('vep.vcf', '/combined.csv')
    add_inhouse_freq()
