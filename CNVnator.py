import copy
import util.count_util as count_util
import util.helper as helper
import argparse

def filter_1(initialFilename):
    """
    Summary:
    - Remove mutant with 'scaff', and only include those with 'Chr'
    - Check the line has exactly one 0/1 and 1/1, and parents contain or 0/0 or ./
    - Check if natorRD <= 0.7 and natorP1 <= 0.01
    - Check if the line contains 'DEL' in TYPE and 'PASS' in FILTER
    - Write a temporary result file for the next step
    """
    print("Start CNVnator 1:")

    inputFile = open(initialFilename, "r")
    print("Opening " + initialFilename)

    outputFileName = initialFilename.split(".vcf")[0] + "_Searched.vcf"
    print("Writing to " + str(outputFileName))
    outputFile = open(outputFileName, "w")

    parentColumns = []

    for line in inputFile:
        # Skip the description lines at the beginning
        if line.startswith("##"):
            continue

        # Split the line into columns
        parsingLine = line.strip().split("\t")

        # Header rows
        if (parsingLine[0].startswith("#")):
            columnHeader = parsingLine

            # Check which indices correspond to given columns, assuming they're not at the beginning.
            # If we're searching for more/other headers, they will need to be added here.
            # Old header names. Before, the indices of interest would be found here, but if they're named differently, it wouldn't find them.
            # Now I'm having it be hardcoded with this as a safety.
            # CHROM  POS ID  REF ALT QUAL    FILTER  INFO    FORMAT
            parentNames = ["953_1013647", "953_1013650", "953_1018153", "953_1018156"]

            for i in range(len(columnHeader)):
                currHeader = columnHeader[i]
                for parent in parentNames:
                    if currHeader == parent:
                        parentColumns.append(i) # UPDATE to the new parents
                        break

            columnHeaderORIG = copy.deepcopy(columnHeader)
            columnHeader.insert(9, "MutantID")
            columnHeader.insert(10, "Genotype")
            columnHeader.insert(11, "SVLEN")
            # columnHeader.insert(13, "RD")
            # columnHeader.insert(14, "P1")

            outputFile.write("\t".join(columnHeader) + "\n")

        # Data rows
        else:
            # Skip ChrUn or ChrSyn or not Chr chromosome
            if "Chr" not in parsingLine[0] or parsingLine[0].startswith("ChrUn") or parsingLine[0].startswith("ChrSy"):
                continue

            # Skip SCAFF chromosome
            if "scaff" in parsingLine[0].lower():
                continue

            # Check the condition of 0/1 and 1/1
            if not count_util.isZeroOneAccepted(parsingLine, parentColumns):
                continue

            # INFO is at column 7th. INFO includes multiple fields separated by ';'
            INFO_ROW = 7
            infoField = parsingLine[INFO_ROW].split(";")

            # Iterate through infoField
            for word in infoField:
                # Set SVLEN
                if word.startswith("SVLEN"):
                    svlen = int(word.split("=")[1])
                # Set SVTYPE
                elif word.startswith("SVTYPE"):
                    svtype = word.split("=")[1]
                # Set NatorRD
                elif word.startswith("natorRD"):
                    natorRD = float(word.split("=")[1])
                # Set NatorP1
                elif word.startswith("natorP1"):
                    natorP1 = float(word.split("=")[1])

            # Check conditions for NatorRD and NatorP1
            if (natorRD > 0.7 or natorP1 > 0.01):
                continue

            # Get DEL only. Get PASS only. Others = discard
            if (svtype.upper() != "DEL" or parsingLine[6] != "PASS"):
                continue

            # Get type hit (1/1 or 0/1) and what position/index it happens
            typeHit, indexHit = count_util.findIndexHit(parsingLine)

            parsingLine.insert(9, columnHeaderORIG[indexHit])
            parsingLine.insert(10, typeHit)
            parsingLine.insert(11, str(svlen))

            outputFile.write("\t".join(parsingLine))
            outputFile.write("\n")

    inputFile.close()
    outputFile.close()
    print("Successful. Closing file", initialFilename)
    print("Done with CNVnator 1\n")
    return outputFileName

def filter_2(inputFilename, parentFilename, depthfile):
    """
    Summary:
        - Remove lines that are included within other row deletions
        - Remove lines that overlap any portion with parent lines
        - Fetch genes from TABIX
        - Write a final mutant gene table for big deletions from CNVnator
    """
    print("Start CNVnator 2")
    
    parentFile = open(parentFilename)
    inputFile = open(inputFilename)
    print("Opening " + parentFilename + " and " + inputFilename)

    parentList = parentFile.readlines()
    outputFilename = inputFilename.split(".vcf")[0] + "_Mutant_Gene_Table.vcf"
    outputFile = open(outputFilename, "w")

    def currPos(line: str):
        """
        Return startPos, endPos, chrome, mutantId
        """
        col = line.split('\t')
        currStart = int(col[1])
        endPos = col[7].split(';')[0]
        endPos = int(endPos.split('=')[1])
        return currStart, endPos, col[0], col[9]

    CHROM_COL = 0
    POS_COL = 1
    INFO_COL = 7
    MUTANT_COL = 9
    GENOTYPE_COL = 10
    SVLEN_COL = 11

    outputFile.write("#mutant_id\tsynonym\tgeneration\tchrom_num\tstart_pos\tend_pos\tchrom_num2\tstart_pos2\tend_pos2\tmut_type\tmut_size\teffect\tgenotype\tgene_name\tgene_desc\n")
    for line in inputFile:
        line = line.strip()
        parsingLine = line.split("\t")

        # Skip header row
        if line.startswith("#"):
            continue

        # Data rows

        # Parse current start position, end position, chromosome, and mutantId
        currStart, currEnd, currChrom, currMutant = currPos(line)
        needRemove = False

        # Check if it's a part of other big deletion
        with open(inputFilename, "r") as compareInput:
            for compareLine in compareInput:
                # Skip column header row
                if compareLine.startswith("#"):
                    continue

                # Data row
                compareLine = compareLine.strip()
                compareStart, compareEnd, compareChrom, compareMutant = currPos(compareLine)
                if (line != compareLine     # Ignore if they are the same
                    and compareMutant == currMutant and compareChrom == currChrom
                    and compareStart <= currStart and compareEnd >= currEnd):
                    needRemove = True
                    break

        if needRemove:
            continue

        # Check if it overlaps with any parent row
        for parentLine in parentList:
            if "#" in parentLine:
                continue
            parentLine = parentLine.strip()
            col = parentLine.split('\t')
            parentChrom = col[0]
            parentStart = int(col[1])
            parentEnd = int(col[2])
            if currChrom == parentChrom:
                if (parentStart <= currStart and parentEnd >= currEnd)\
                or (parentStart >= currStart and parentStart < currEnd)\
                or (parentEnd > currStart and parentEnd <= currEnd):
                    needRemove = True
                    break

        if needRemove:
            continue

        chrom = parsingLine[CHROM_COL]
        pos = parsingLine[POS_COL]
        mutantId = parsingLine[MUTANT_COL]
        endPos = parsingLine[INFO_COL].split(";")[0].split("=")[1]
        mut_size = str(abs(int(parsingLine[SVLEN_COL])))
        genotype = helper.get_genotype(parsingLine[GENOTYPE_COL])

        # Define coordinate and get gene list if any
        coordinate = f"{chrom}:{pos}-{endPos}"
        gene_list = helper.tabix_genes(coordinate, depthfile)

        output_arr = [
            mutantId,    # mutant_id
            "",          # synonym
            "",          # generation
            chrom,       # chrom_num
            pos,         # start_pos
            endPos,      # end_pos
            "",          # chrom_num2
            "0",         # start_pos2
            "0",         # end_pos2
            "Deletion",  # mut_type
            mut_size,    # mut_size
            "",          # effect
            genotype,    # genotype
            "",          # gene_name - index 13
            "",          # gene_desc
            ""           # pheno_desc
        ]
        GENE_INDEX = 13
        if len(gene_list) == 0:
            outputFile.write("\t".join(output_arr) + "\n")
        else:
            for i in gene_list:
                output_arr[GENE_INDEX] = i
                outputFile.write("\t".join(output_arr) + "\n")

    inputFile.close()
    outputFile.close()
    parentFile.close()
    print("Done. Close", inputFilename, parentFilename)
    print("Successful produce", outputFilename)
    print("Done with CNVnator 2\n")

def main():
    parser = argparse.ArgumentParser(description = "Analyze deletion in CNVnator")
    parser.add_argument('-f', '--file', type=str, required=True)
    parser.add_argument('-p', '--parent', type=str, required=True)
    parser.add_argument('-d', '--depthfile', type=str, required=True)

    args = parser.parse_args()
    first_filename = filter_1(args.file)
    filter_2(first_filename, args.parent, args.depthfile)

if __name__ == '__main__':
    main()
