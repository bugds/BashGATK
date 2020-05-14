import os
import sys

wd = os.path.abspath(sys.argv[1])

global innerDelimeter
global outerDelimeter
global extensionFilter
global outputName

innerDelimeter = ','
outerDelimeter = '\t'
extensionFilter = 'vep.vcf.pass.vcf'
outputName = '/combined_passed.csv'

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

        infoLines = [l.replace('\n', '') for l in lines if l.startswith('##INFO=')]
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
    annoChunk = annoData[:annoData.index('ALLELE_END')]
    for c in range(0, annoData.count('ALLELE_END')):
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
                    break
            if not valueFound:
                varDict['INFO_ANNO_' + key] = '.'
        annoData = annoData[annoData.index('ALLELE_END') + 1:]
        try:
            annoChunk = annoData[:annoData.index('ALLELE_END')]
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

def main():
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
                    
main()
