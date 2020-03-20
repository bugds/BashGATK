import os
import sys

wd = os.path.abspath(sys.argv[1])

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
                
def main():
    headerDone = False
    
    with open(wd + '/combined.csv', 'w') as csvFile:
        for filename in os.listdir(wd + '/annotation'):
                if filename.endswith('vep.vcf'):
                    with open(wd + '/annotation/' + filename, 'r') as oneFile:
                        lines = oneFile.readlines()
                    
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
                    
                    for line in dataLines:
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
                                                varDict['INFO_ANNO_' + key] += ',' + element.split('=')[1]
                                            except IndexError:
                                                varDict['INFO_ANNO_' + key] += ',' + element
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
                        
                        vepData = varDict['INFO'].split(';')[-1].replace('CSQ=', '').split('|')
                        varDict.pop('INFO')
                        for i in range(0, len(vepData)):
                            if vepData[i] == '':
                                vepData[i] = '.'
                        vepChunk = vepData[:len(vepList)-1]
                        while len(vepChunk) > len(vepList)-2:
                            for i in zip(vepList, vepChunk):
                                if ('INFO_VEP_' + i[0]) in varDict:
                                    if i[1][0] == ',':
                                        varDict['INFO_VEP_' + i[0]] += i[1]
                                    else:
                                        varDict['INFO_VEP_' + i[0]] += ',' + i[1]
                                else:
                                    varDict['INFO_VEP_' + i[0]] = i[1]
                            if not headerDone:
                                csvFile.write('\t'.join(k for k in varDict.keys()) + '\n')
                                headerDone = True
                            vepData = vepData[len(vepList)-1:]
                            vepChunk = vepData[:len(vepList)-1]
                        for i in varDict:
                            setOfElements = set(varDict[i].split(','))
                            setOfElements.discard('.')
                            if len(setOfElements) == 1:
                                varDict[i] = next(iter(setOfElements))
                            elif len(setOfElements) == 0:
                                varDict[i] = '.'
                        csvFile.write('\t'.join(v for v in varDict.values()) + '\n')
                    
main()
