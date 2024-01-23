import argparse

def get_gene_desc_dict(filename: str) -> dict:
    file = open(filename)
    gene_dict = dict()
    GENE_COL = 1
    GENE_DESC_COL = 9
    for line in file:
        # Skip the header lines
        if line.startswith("#"):
            continue

        parsingLine = line.split("\t")
        gene = parsingLine[GENE_COL]
        desc = parsingLine[GENE_DESC_COL].strip()
        if gene != "" and desc != "":
            gene_dict[gene] = desc
    file.close()
    return gene_dict

def get_gene_desc_mutant_table(filename: str) -> dict:
    file = open(filename)
    GENE_COL = 13
    GENE_DESC_COL = 14
    gene_dict = dict()
    for line in file:
        # Skip the header lines
        if line.startswith("#"):
            continue

        parsingLine = line.split("\t")
        gene = parsingLine[GENE_COL]
        desc = parsingLine[GENE_DESC_COL]
        if gene != "" and desc != "":
            gene_dict[gene] = desc
    file.close()
    return gene_dict

def add_gene_desc(filename: str, gene_dict: dict):
    GENE_COL = -3
    GENE_DESC_COL = -2
    file = open(filename)
    output = open(filename.split(".vcf")[0] + "_Gene_Desc.vcf", "w")
    for line in file:
        if line.startswith("#"):
            output.write(line)
            continue
        parsingLine = line.split("\t")
        gene = parsingLine[GENE_COL]
        desc_from_table = parsingLine[GENE_DESC_COL].strip()
        if gene != "":
            desc = gene_dict.get(gene)
            if desc == None:
                print(gene, "is not found")
            else:
                if desc_from_table != "" and desc_from_table != desc and desc_from_table not in desc:
                    print(gene, "has different descriptions")
                    print(desc_from_table)
                    print(desc)
                parsingLine[GENE_DESC_COL] = desc
        output.write("\t".join(parsingLine))

    file.close()
    output.close()

def main():
    parser = argparse.ArgumentParser(description = "Add gene description into mutant gene table")
    parser.add_argument('-f', '--file', type=str, required=True, help="Mutant gene table file")
    parser.add_argument('-g', '--gene', type=str, required=True, help="Gene file")

    args = parser.parse_args()
    gene_dict = get_gene_desc_dict(args.gene)
    add_gene_desc(args.file, gene_dict)

if __name__ == '__main__':
    main()
