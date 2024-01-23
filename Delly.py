import copy
import util.count_util as count_util
import util.helper as helper
import argparse

def delly_1(input_filename):
    """
    Summary:
    - Remove mutant with 'scaff', and only include those with 'Chr'
    - Check the line has exactly one 0/1 and 1/1, and parents contain or 0/0 or ./
    - Check if type is either INV (inversion) or BND (translocation)
    - Check if the line contains 'PASS' in FILTER
    - Write a temporary result file for the next step
    """
    print("Start delly 1:")

    inputFile = open(input_filename, "r")
    print("Opening " + input_filename)

    outputFileName = input_filename.split(".vcf")[0] + "_Searched.vcf"
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
            # Now I'm' having it be hardcoded with this as a safety.
            # CHROM  POS ID  REF ALT QUAL    FILTER  INFO    FORMAT
            parentNames = ["Kit-Oryza_sativa_Kitaake_ubi-XA21_sample_A",
                           "Kit-Oryza_sativa_Kitaake_ubi-XA21_sample_B",
                           "Kit-FN_ubi-XA21_5-4",
                           "Kit-FN_ubi-XA21"]
            for i in range(0, len(columnHeader)):
                currHeader = columnHeader[i]
                for parent in parentNames:
                    if currHeader == parent:
                        parentColumns.append(i) # UPDATE to the new parents
                        break

            columnHeaderORIG = copy.deepcopy(columnHeader)
            columnHeader.insert(9, "MutantID")
            columnHeader.insert(10, "Genotype")

            # Write header to output file
            outputFile.write("\t".join(columnHeader) + "\n")

        # Data rows
        else:
            # Skip ChrUn or ChrSyn or not Chr chromosome
            if "Chr" not in parsingLine[0] or parsingLine[0].startswith("ChrUn") or parsingLine[0].startswith("ChrSy"):
                continue

            # Skip SCAFF chromosome
            if "scaff" in parsingLine[0].lower():
                continue

            # Get INV and TRA only. Others = discard
            if ("INV" not in line and "BND" not in line):
                continue

            # Check if the FILTER is PASS. If not, discard
            if (parsingLine[6] != "PASS"):
                continue

            # Check the condition of 0/1 and 1/1
            if not count_util.isZeroOneAccepted(parsingLine, parentColumns):
                continue

            # Get type hit (1/1 or 0/1) and what position/index it happens
            typeHit, indexHit = count_util.findIndexHit(parsingLine)

            parsingLine.insert(9, columnHeaderORIG[indexHit])
            parsingLine.insert(10, typeHit)

            outputFile.write("\t".join(parsingLine) + "\n")

    inputFile.close()
    outputFile.close()
    print("Successful. Closing file", input_filename)
    print("Done with Delly 1\n")
    return outputFileName

def delly_2(inputFileName, depthfile):
    print("Start delly 2:")
    print("Opening", inputFileName)
    delly_file = open(inputFileName)

    # Open 2 output files: INV and Translocation
    inversion_filename = inputFileName.split(".vcf")[0] + "_Inv_Mutant_Gene_Table.vcf"
    trans_filename = inputFileName.split(".vcf")[0] + "_Trans_Mutant_Gene_Table.vcf"
    inversion_file = open(inversion_filename, "w")
    trans_file = open(trans_filename, "w")

    inversion_file.write("#mutant_id\tsynonym\tgeneration\tchrom_num\tstart_pos\tend_pos\tchrom_num2\tstart_pos2\tend_pos2\tmut_type\tmut_size\teffect\tgenotype\tgene_name\tgene_desc\n")
    trans_file.write("#mutant_id\tsynonym\tgeneration\tchrom_num\tstart_pos\tend_pos\tchrom_num2\tstart_pos2\tend_pos2\tmut_type\tmut_size\teffect\tgenotype\tgene_name\tgene_desc\n")

    CHROM_COL = 0
    POS_COL = 1
    ID_COL = 2
    INFO_COL = 7
    MUTANT_COL = 9
    GENOTYPE_COL = 10

    for line in delly_file:
        # Skip header row
        if line.startswith("#"):
            continue

        # Split the line into columns
        parsingLine = line.strip().split("\t")

        # Define data in mutant gene table row
        chrom = parsingLine[CHROM_COL]
        pos = parsingLine[POS_COL]
        mutantId = parsingLine[MUTANT_COL]

        info_parse = parsingLine[INFO_COL].split(";")

        genotype = helper.get_genotype(parsingLine[GENOTYPE_COL])

        # Split the line into columns
        parsingLine = line.strip().split("\t")
        id = parsingLine[ID_COL]

        def combine_list(list1, list2):
            return list(set(list1 + list2))

        # Inversion type
        if id.upper().startswith("INV"):
            for info in info_parse:
                map = info.split("=")
                if len(map) == 2:
                    key, value = map
                    if key.upper() == "END":
                        endPos = value

            mut_size = str(int(endPos) - int(pos))

            # Define coordinate and get gene list if any
            coordinate1 = f"{chrom}:{pos}-{pos}"
            coordinate2 = f"{chrom}:{endPos}-{endPos}"
            gene_list1 = helper.tabix_genes(coordinate1, depthfile)
            gene_list2 = helper.tabix_genes(coordinate2, depthfile)

            final_gene_list = combine_list(gene_list1, gene_list2)

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
                "Inversion", # mut_type
                mut_size,    # mut_size
                "",          # effect
                genotype,    # genotype
                "",          # gene_name - index 13
                "",          # gene_desc
                ""           # pheno_desc
            ]
            GENE_INDEX = 13

            if len(final_gene_list) == 0:
                inversion_file.write("\t".join(output_arr) + "\n")
            else:
                for i in final_gene_list:
                    output_arr[GENE_INDEX] = i
                    inversion_file.write("\t".join(output_arr) + "\n")

        # Translocation type
        elif id.upper().startswith("BND"):
            # Get END position in INFO
            for info in info_parse:
                map = info.split("=")
                if len(map) == 2:
                    key, value = map
                    if key.upper() == "POS2":
                        pos2 = value
                    elif key.upper() == "CHR2":
                        chrom2 = value

            # Define coordinate and get gene list if any
            coordinate1 = f"{chrom}:{pos}-{pos}"
            coordinate2 = f"{chrom2}:{pos2}-{pos2}"
            gene_list1 = helper.tabix_genes(coordinate1, depthfile)
            gene_list2 = helper.tabix_genes(coordinate2, depthfile)

            final_gene_list = combine_list(gene_list1, gene_list2)

            output_arr = [
                mutantId,         # mutant_id
                "",               # synonym
                "",               # generation
                chrom,            # chrom_num
                pos,              # start_pos
                "",               # end_pos
                chrom2,           # chrom_num2
                pos2,             # start_pos2
                "",               # end_pos2
                "Translocation",  # mut_type
                "0",              # mut_size
                "",               # effect
                genotype,         # genotype
                "",               # gene_name - index 13
                "",               # gene_desc
                ""                # pheno_desc
            ]
            GENE_INDEX = 13

            if len(final_gene_list) == 0:
                trans_file.write("\t".join(output_arr) + "\n")
            else:
                for i in final_gene_list:
                    output_arr[GENE_INDEX] = i
                    trans_file.write("\t".join(output_arr) + "\n")

    delly_file.close()
    inversion_file.close()
    trans_file.close()
    print("Done. Close", inputFileName)
    print(f"Successful produce {inversion_filename} and {trans_filename}")
    print("Done with Delly 2\n")

def main():
    parser = argparse.ArgumentParser(description = "Analyze translocation and inversion in Delly")
    parser.add_argument('-f', '--file', type=str, required=True)
    parser.add_argument('-d', '--depthfile', type=str, required=True)

    args = parser.parse_args()
    first_filename = delly_1(args.file)
    delly_2(first_filename, args.depthfile)

if __name__ == '__main__':
    main()
