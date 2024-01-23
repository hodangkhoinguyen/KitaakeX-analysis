import copy
import util.count_util as count_util
import util.helper as helper
import argparse

def samtool(inputFilename):
    print("Start analyzing Samtools file")

    inputFile = open(inputFilename)
    print("Opening " + inputFilename)

    insertion_filename = inputFilename.split(".vcf")[0] + "_Ins_Mutant_Gene_Table.vcf"
    deletion_filename = inputFilename.split(".vcf")[0] + "_Del_Mutant_Gene_Table.vcf"
    snp_filename = inputFilename.split(".vcf")[0] + "_Snp_Mutant_Gene_Table.vcf"

    insertion_file = open(insertion_filename, "w")
    deletion_file = open(deletion_filename, "w")
    snp_file = open(snp_filename, "w")

    insertion_file.write("#mutant_id\tsynonym\tgeneration\tchrom_num\tstart_pos\tend_pos\tchrom_num2\tstart_pos2\tend_pos2\tmut_type\tmut_size\teffect\tgenotype\tgene_name\tgene_desc\n")
    deletion_file.write("#mutant_id\tsynonym\tgeneration\tchrom_num\tstart_pos\tend_pos\tchrom_num2\tstart_pos2\tend_pos2\tmut_type\tmut_size\teffect\tgenotype\tgene_name\tgene_desc\n")
    snp_file.write("#mutant_id\tsynonym\tgeneration\tchrom_num\tstart_pos\tend_pos\tchrom_num2\tstart_pos2\tend_pos2\tmut_type\tmut_size\teffect\tgenotype\tgene_name\tgene_desc\n")

    parentColumns = []

    CHROM_COL = 0
    POS_COL = 1
    REF_COL = 3
    ALT_COL = 4
    QUAL_COL = 5
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
            #CHROM  POS ID  REF ALT QUAL    FILTER  INFO    FORMAT
            parentNames = ["Kit-Oryza_sativa_Kitaake_ubi-XA21_sample_A",
                           "Kit-Oryza_sativa_Kitaake_ubi-XA21_sample_B",
                           "Kit-FN_ubi-XA21_5-4",
                           "Kit-FN_ubi-XA21"]

            for i in range(len(columnHeader)):
                currHeader = columnHeader[i]
                for parent in parentNames:
                    if currHeader == parent:
                        parentColumns.append(i) # UPDATE to the new parents
                        break

            columnHeaderORIG = copy.deepcopy(columnHeader)
            # columnHeader.insert(9, "MutantID")
            # columnHeader.insert(10, "Genotype")
            
            ## Writing header to 2 output files. SNP first: doesn't need insertion/deletion and length columns
            # outputNO_INDEL.write("\t".join(columnHeader) + "\n")

            # columnHeader.insert(5, "I or D")
            # columnHeader.insert(6, "Length")
            # outputINDEL.write("\t".join(columnHeader) + "\n")

        # Data rows
        else:
            chrom = parsingLine[CHROM_COL]

            # Skip ChrUn or ChrSyn or not Chr chromosome
            if "Chr" not in chrom or chrom.startswith("ChrUn") or chrom.startswith("ChrSy"):
                continue

            # Skip SCAFF chromosome
            if "scaff" in chrom.lower():
                continue

            # Check the condition of 0/1 and 1/1
            if not count_util.isZeroOneAccepted(parsingLine, parentColumns):
                continue
            if not count_util.isNotContainZeroTwo(parsingLine):
                continue

            # Make sure quality is above 100, or else discard the row.
            try:
                quality = float(parsingLine[QUAL_COL])
            except ValueError:
                print("Quality of row isn't a number!")
                continue
            if quality < 100:
                continue

            # PASS only. Others = discard
            if parsingLine[6] != "PASS":
                continue

            # Get type hit (1/1 or 0/1) and what position/index it happens
            typeHit, indexHit = count_util.findIndexHit(parsingLine)

            mutantId = columnHeaderORIG[indexHit]
            genotype = helper.get_genotype(typeHit)

            pos = parsingLine[POS_COL]

            isIndel = "INDEL" in line
            # Indel: Insertion or Deletion
            if isIndel:
                ref = parsingLine[REF_COL]
                alt = parsingLine[ALT_COL].split(',')[0]
                length = len(alt) - len(ref) # Compare REF with ALT.

                if length > 0:
                    mut_type = "Insertion"
                    endPos = pos
                    outputFile = insertion_file
                elif length < 0:
                    mut_type = "Deletion"
                    endPos = str(int(pos) + abs(length) - 1)
                    outputFile = deletion_file
                else:
                    print("This is not INDEL because of no length difference")
                    continue
                mut_size = abs(length)

            # Single Base Substitution
            else:
                mut_type = "Single Base Substitution"
                mut_size = 0
                endPos = f"{parsingLine[REF_COL]}-{parsingLine[ALT_COL]}"
                outputFile = snp_file

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
                    # Get gene in EFF
                    if info_elem.startswith("EFF"):
                        value = info_elem.split("=")[1]
                        effects = value.split(",")
                        for effect in effects:
                            parsingEffect = effect.split("|")
                            if "HIGH" in parsingEffect[0] or "MODERATE" in parsingEffect[0]:
                                gene_set.add(parsingEffect[5])
                        break
                return list(gene_set)

            gene_list = extract_gene(parsingLine[INFO_COL])
            if len(gene_list) == 0:
                outputFile.write("\t".join(output_arr) + "\n")
            else:
                for i in gene_list:
                    output_arr[GENE_INDEX] = i
                    outputFile.write("\t".join(output_arr) + "\n")

    print("Successful. Closing file " + inputFilename)
    inputFile.close()
    insertion_file.close()
    deletion_file.close()
    snp_file.close()

def main():
    parser = argparse.ArgumentParser(description = "Analyze Deletion, Insertion, and SBS in Samtools")
    parser.add_argument('-f', '--file', type=str, required=True)

    args = parser.parse_args()
    samtool(args.file)

if __name__ == '__main__':
    main()
