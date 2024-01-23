import subprocess

def tabix_genes(coordinate: str, depthfile: str) -> list[str]:
    """
    Return a list of genes based on provided coordinate and depth file of genes
    coordinate (str): [chromosome]:[start]-[end]
    depthfile (str): a sorted gff.gz file
    """
    # Define the command of tabix
    command = f"tabix {depthfile} {coordinate}"

    # Perform shell command
    p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output, _ = p.communicate() # Genereate return code
    output = output.decode("UTF-8").strip()

    # Check if the shell command run well
    if p.returncode != 0:
        raise Exception(f"Run tabix failed\n{output}")

    # When there is no gene
    if output == "":
        return []

    # Parse the gene list containing unique genes
    line_list = output.split("\n")
    gene_set = set()
    for line in line_list:
        # line = Chr3    MSU_osa1r7      exon    1367718 1368572 .       +       .       ID=LOC_Os03g03210.1:exon_6;Parent=LOC_Os03g03210.1
        parsingLine = line.split("\t")
        # id = ID=LOC_Os03g03210.1:exon_6;Parent=LOC_Os03g03210.1
        id = parsingLine[-1]
        # gene = LOC_Os03g03210
        gene = id.split(";")[0].split("=")[1].split(".")[0]
        gene_set.add(gene)

    return list(gene_set)

def get_genotype(input: str) -> str:
    """
    Get genotype as HET (Heterozygous) or HOMO (Homozygous)

    Args:
        input (str): 0/1 or 1/1 (a.k.a 1/1.)
    """
    if input == "0/1":
        return "HET"
    if input == "1/1" or "1/1.":
        return "HOMO"
    raise Exception(f"{input} is not 0/1, 1/1, or 1/1")
