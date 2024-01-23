import copy
import util.count_util as count_util
import util.helper as helper
import os
import argparse

def get_samtool(samtool_filename):
    samtool_file = open(samtool_filename)
    samtool_set = set()
    for line in samtool_file:
        # Skip description and header lines
        if line.startswith("#"):
            continue
        parsingLine = line.strip().split("\t")
        chrom = parsingLine[0]
        pos = parsingLine[1]
        samtool_set.add(f"{chrom}\t{pos}")
    samtool_file.close()
    return samtool_set

def pindel_filter(inputFileName, samtool_filename):
    print("Start analyzing Pindel file")
    inputFile = open(inputFileName)

    outputFileName = inputFileName.split(".vcf")[0] + "_Searched.vcf"
    print("Writing to " + str(outputFileName))
    outputFile = open(outputFileName, "w")

    samtool_set = get_samtool(samtool_filename)
    parentColumns = []

    CHROM_COL = 0
    POS_COL = 1
    ID_COL = 2
    INFO_COL = 7

    for line in inputFile:
        # Skip description lines at the beginning
        if line.startswith("##"):
            continue

        # Split the line into columns
        parsingLine = line.strip().split("\t")

        # Header rows
        if (parsingLine[0].startswith("#")):
            columnHeader = parsingLine

            ## Check which indices correspond to given columns, assuming they're not at the beginning.
            ## If we're searching for more/other headers, they will need to be added here.
            ## Old header names. Before, the indices of interest would be found here, but if they're named differently, it wouldn't find them.
            ## Now I'm' having it be hardcoded with this as a safety.
            # CHROM  POS ID  REF ALT QUAL    FILTER  INFO    FORMAT
            parentNames = ["Oryza_sativa_Kitaake_ubi-XA21_sample_A",
                "Oryza_sativa_Kitaake_ubi-XA21_sample_B",
                "FN_ubi-XA21_5-4",
                "FN_ubi-XA21"]

            for i in range(len(columnHeader)):
                currHeader = columnHeader[i]
                for parent in parentNames:
                    if currHeader == parent:
                        parentColumns.append(i) # UPDATE to the new parents
                        break

            columnHeaderORIG = copy.deepcopy(columnHeader)
            outputFile.write("\t".join(columnHeader[:(INFO_COL + 1)]) + "\n") # Cut file upto INFO

        # Data rows
        else:
            chrom = parsingLine[CHROM_COL]

            # Skip ChrUn or ChrSyn or not Chr chromosome
            if "Chr" not in chrom or chrom.startswith("ChrUn") or chrom.startswith("ChrSy"):
                continue

            # Skip SCAFF chromosome
            if "scaff" in chrom.lower():
                continue

            # Check if the FILTER is PASS. If not, discard
            if (parsingLine[6] != "PASS"):
                continue

            # Check if it's NOT in samtool
            pos = parsingLine[POS_COL]
            if f"{chrom}\t{pos}" in samtool_set:
                continue

            # Check the condition of 0/1 and 1/1. Discard if it's False - meaning not meet requirement
            if (count_util.isZeroOneAccepted(parsingLine, parentColumns) == False):
                continue
            # If it contains any 0/2, 1/2, 2/2, then discard
            if (count_util.isNotContainZeroTwo(parsingLine) == False):
                continue

            # INFO is at column 7th. INFO includes multiple fields separated by ';'
            infoField = parsingLine[INFO_COL].split(";")
            for word in infoField:
                # Set SVLEN
                if word.startswith("SVLEN"):
                    svlen = int(word.split("=")[1])
                # Set SVTYPE
                elif word.startswith("SVTYPE"):
                    svtype = word.split("=")[1]

            # Get INS or DEL only. Others = discard
            if (svtype != "INS" and svtype != "DEL"):
                continue
            # Only get INS with SVLEN > 0. Others = discard
            if (svtype == "INS" and svlen == 0):
                continue
            if abs(svlen) >= 1000:
                if svtype != "DEL":
                    print(line)
                continue

            # Get type hit (1/1 or 0/1) and what position/index it happens
            typeHit, indexHit = count_util.findIndexHit(parsingLine)

            # Get values of AD and DP - 0/0:13,1:65 => AD = 13, DP = 65
            ad = float(parsingLine[indexHit].split(":")[1].split(",")[0])
            dp = float(parsingLine[indexHit].split(":")[2])

            # Check condition of AD, DP
            if not (ad >= 8 and dp >= 10):
                continue
            ad_dp_percent = ad / dp * 100
            if not (ad_dp_percent >= 30):
                continue

            parsingLine[ID_COL] = columnHeaderORIG[indexHit] # set ID to MUTANT ID
            parsingLine[INFO_COL] += f";GENOTYPE={helper.get_genotype(typeHit)}" # set genotype into INFO

            outputFile.write("\t".join(parsingLine[:(INFO_COL + 1)]) + "\n")

    print("Successful. Closing file " + inputFileName)
    inputFile.close()
    outputFile.close()
    return outputFileName

def run_snpEff(database, input_pindel_file):
    print("Perform snpEff")
    output_file = input_pindel_file.split(".vcf")[0] + "_SnpEff.vcf"
    command = f"java -Xmx8g -jar ~/snpEff/snpEff.jar {database} {input_pindel_file} > {output_file}"
    return_code = os.system(command)

    # Handle error when running snpEff
    if return_code != 0:
        print("SnpEff failed")
        exit()
    print("Done snpEff")
    return output_file

def produce_mutant_gene(input_filename):
    print("Start producing mutant gene table")
    input_file = open(input_filename)

    insertion_filename = input_filename.split(".vcf")[0] + "_Insertion_Mutant_Gene_Table.vcf"
    deletion_filename = input_filename.split(".vcf")[0] + "_Deletion_Mutant_Gene_Table.vcf"
    insertion_file = open(insertion_filename , "w")
    deletion_file = open(deletion_filename , "w")

    insertion_file.write("#mutant_id\tsynonym\tgeneration\tchrom_num\tstart_pos\tend_pos\tchrom_num2\tstart_pos2\tend_pos2\tmut_type\tmut_size\teffect\tgenotype\tgene_name\tgene_desc\n")
    deletion_file.write("#mutant_id\tsynonym\tgeneration\tchrom_num\tstart_pos\tend_pos\tchrom_num2\tstart_pos2\tend_pos2\tmut_type\tmut_size\teffect\tgenotype\tgene_name\tgene_desc\n")

    CHROM_COL = 0
    POS_COL = 1
    ID_COL = 2
    INFO_COL = 7

    for line in input_file:
        if line.startswith("#"):
            continue
        parsingLine = line.strip().split("\t")

        mutantId = parsingLine[ID_COL]
        chrom = parsingLine[CHROM_COL]
        pos = parsingLine[POS_COL]
        info = parsingLine[INFO_COL]

        parsingInfo = info.split(";")
        for i in parsingInfo:
            key, value = i.split("=")
            if key == "SVLEN":
                mut_size = abs(int(value))
            elif key == "SVTYPE":
                if value == "INS":
                    mut_type = "Insertion"
                elif value == "DEL":
                    mut_type = "Deletion"
                else:
                    print("Error", value, "type is not deletion or insertion")

            elif key == "GENOTYPE":
                genotype = value

        if mut_type == "Insertion":
            endPos = pos
            output_file = insertion_file
        elif mut_type == "Deletion":
            endPos = str(int(pos) + mut_size - 1)
            output_file = deletion_file
        mut_size = str(mut_size)
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
            mut_type,    # mut_type
            mut_size,    # mut_size
            "",          # effect
            genotype,    # genotype
            "",          # gene_name - index 13
            "",          # gene_desc
            ""           # pheno_desc
        ]

        GENE_INDEX = 13
        def extract_gene(info):
            parsingInfo = info.split(";")
            gene_set = set()
            for info_elem in parsingInfo:
                # Get gene in ANN
                if info_elem.startswith("ANN"):
                    value = info_elem.split("=")[1]
                    effects = value.split(",")
                    for effect in effects:
                        parsingEffect = effect.split("|")
                        if "HIGH" in parsingEffect[2] or "MODERATE" in parsingEffect[2]:
                            gene_set.add(parsingEffect[3].split(".")[0])
                    break
            return list(gene_set)

        gene_list = extract_gene(parsingLine[INFO_COL])
        if len(gene_list) == 0:
            output_file.write("\t".join(output_arr) + "\n")
        else:
            for i in gene_list:
                output_arr[GENE_INDEX] = i
                output_file.write("\t".join(output_arr) + "\n")

    print("Successful. Done producing mutant gene table")
    input_file.close()
    deletion_file.close()
    insertion_file.close()

def main():
    parser = argparse.ArgumentParser(description = "Analyze Deletion and Insertion in Pindel")
    parser.add_argument('-f', '--file', type=str, required=True)
    parser.add_argument('-d', '--database', type=str, required=True)
    parser.add_argument('-s', '--samtool', type=str, required=True)

    args = parser.parse_args()
    filter_file = pindel_filter(args.file, args.samtool)
    snpEff_file = run_snpEff(args.database, filter_file)
    print(snpEff_file)
    produce_mutant_gene(snpEff_file)

if __name__ == '__main__':
    main()
