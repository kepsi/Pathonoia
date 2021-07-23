#!/usr/bin/env python
# coding: utf-8

#import csv
import pandas as pd

def evalKrakenAlign(outdir_prefix, unmappedreads, nucCutOff, logfile):
    krakenreportfile = outdir_prefix + "_krakenreport.txt"
    krakenalignfile = outdir_prefix + "_krakenalign.txt"
    kmermap = xmerMapFromKrakenAligned(krakenalignfile, unmappedreads, logfile)
    #savedictToCSVFile(kmermap, outdir_prefix+ "_kmerTaxMap.csv")
    nucPerTax = evaluateXmerMap(kmermap,krakenreportfile, nucCutOff)
    nucPerTax.to_csv(outdir_prefix+ "_nucleotiesPerTaxID.csv")
    return(nucPerTax)

def evaluateXmerMap(xmerMap, krakenReportFile, nucCutOff = 50):
    
    taxInfo = getTaxInfoFromKrakenReportFile(krakenReportFile) 
    taxMap = dict()
    for xmer, taxId in xmerMap.items():
        if(taxId in taxMap):
            taxMap[taxId] = taxMap[taxId] + len(xmer)
        else:
            taxMap[taxId] = len(xmer)
    
    df = pd.DataFrame.from_dict(taxMap, orient='index', columns=['nucleotides'])
    df.sort_values(by = "nucleotides", ascending=False, inplace=True)
    
    
    # NEW PART: integrate parent information
    df_all = df.join(taxInfo)
    df = df_all.loc[df_all.phylo_level.str.contains("S", na=False),]
    df = df.loc[df["nucleotides"]>nucCutOff,]

    for taxId in df.index:
        species_nucs = df.loc[taxId,"nucleotides"]
        parent = df.loc[taxId,"parent"]
        level = taxInfo.loc[parent,"phylo_level"]
        while("S" in level or "G" in level or "F" in level):
            if(parent in df_all.index):
                species_nucs = species_nucs + df_all.loc[parent,"nucleotides"]
            parent = taxInfo.loc[parent,"parent"]
            level = taxInfo.loc[parent,"phylo_level"]
        df.loc[taxId,"nucleotides"] = species_nucs
        
    return(df)

def xmerMapFromKrakenAligned(inputfile_krakenAlign, inputfile_unmappedreads, logfile):    
    if(testFilesCorrespondingReads(inputfile_krakenAlign, inputfile_unmappedreads, logfile)):
        kmer_map = dict()
        with open(inputfile_krakenAlign) as file_kraken:
            with open(inputfile_unmappedreads) as file_unmapped:
                kraken_line = file_kraken.readline().rstrip()                    
                #calculate K from first line
                readlength = int(kraken_line.split("\t")[3])
                aligntuple = kraken_line.split("\t")[4].split(" ")
                total_kmers = 0
                for tup in aligntuple:
                    total_kmers = total_kmers + int(tup.split(":")[1])
                K = readlength - total_kmers +1
                print("kmer-length in Kraken Report is: " + str(K), file = logfile)
                count_reads_to_check = 0
                while kraken_line:
                    file_unmapped.readline() #let go, this is just the read identifier and was tested
                    sequence = file_unmapped.readline()
                    line_parts = kraken_line.split("\t")
                    readlength = int(line_parts[3])
                    aligntuple = line_parts[4].split(" ")
                    species = line_parts[2]
                    position = 0
                    for tup in aligntuple:
                        taxAndNumber = tup.split(":")
                        tax = taxAndNumber[0]
                        number = int(taxAndNumber[1])
                        if (tax not in ("0", "28384", "1") and number > 4):
                            xmer = sequence[position:position + K+number-1]
                            if(xmer in kmer_map):
                                first_tax = kmer_map[xmer]
                                if(first_tax != tax):
                                    print("this kmer was assigned to a different taxID earlier:" + first_tax + " first tax vs new tax " + tax)
                            kmer_map[xmer] = tax
                        position = position + number

                    kraken_line = file_kraken.readline().rstrip()
        return(kmer_map)
        
#def savedictToCSVFile(dictToSave, outputfile):
#    w = csv.writer(open(outputfile, "w"))
#    w.writerow(["xmer", "taxId"])
#    for key, val in dictToSave.items():
#        w.writerow([key, val])
        
def getTaxInfoFromKrakenReportFile(reportFile):
    taxInfo = dict()
    phyloDict = dict()
    phyloDict[-2] = "NA"
    with open(reportFile) as reportFile:
        line = reportFile.readline()
        while line:
            p = line.split("\t")
            taxId = p[4]
            taxlevel = p[3]
            organism = p[5]
            species_name = organism.lstrip().rstrip()
            phyloDict[level(organism)] = taxId            
            parent = phyloDict[level(organism)-2]            
            taxInfo[taxId] = [species_name,taxlevel,parent]          
            line = reportFile.readline()
    taxInfo = pd.DataFrame.from_dict(taxInfo, orient='index', columns=['species_name','phylo_level','parent'])        
    return(taxInfo)

def level(organism):
    a = organism
    return len(a) - len(a.lstrip(' '))

def testFilesCorrespondingReads(inputfile_krakenAlign, inputfile_unmappedreads, logfile, numberLinesToTest = 500):
    lines_tested = 0
    with open(inputfile_krakenAlign) as file_kraken:
        with open(inputfile_unmappedreads) as file_unmapped:
            readname = file_kraken.readline().split()[1]
            unmapped_line = file_unmapped.readline()
            while(lines_tested<numberLinesToTest):
                if(readname not in unmapped_line):
                    print("ERROR: corresponding test failed for files: " + inputfile_krakenAlign + " and " + inputfile_unmappedreads, file = logfile)
                    return(False)
                lines_tested = lines_tested + 1                
                readname = file_kraken.readline().split()[1]
                unmapped_line = file_unmapped.readline()
                unmapped_line = file_unmapped.readline()
    return(True)

def mergeNucTaxDataframes(allNucPerTax):
    infocols = ['species_name', 'phylo_level', 'parent']
    
    sampleList = list(allNucPerTax.keys())
    listOfTaxInfo = []
    for sample, df in allNucPerTax.items():
        listOfTaxInfo.append(df[infocols])

    taxInfoDF = pd.concat(listOfTaxInfo)
    taxInfoDF = taxInfoDF.drop_duplicates()
    
    workingDF = allNucPerTax[sampleList[0]]    
    workingDF.drop(infocols, axis=1, inplace=True)
    workingDF.columns = [sampleList[0]]
    for i in range(1,len(sampleList)):
        toBeAddedDF = allNucPerTax[sampleList[i]]
        toBeAddedDF.drop(infocols, axis=1, inplace=True)
        toBeAddedDF.columns = [sampleList[i]]
        workingDF = pd.concat([workingDF, toBeAddedDF], axis=1, sort=False)

    workingDF["rowSum"] = workingDF.sum(axis = 1)
    workingDF = pd.concat([taxInfoDF,workingDF], axis = 1, sort=False)
    workingDF.sort_values(by = 'rowSum', ascending=0, inplace = True)
    workingDF.drop("rowSum", axis = 1, inplace=True)
    return(workingDF)
    

