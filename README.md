# KitaakeX Analysis

## Purpose

This project is to perform rice analysis on KitaakeX from CNVnator, Delly, Pindel, and Samtools to output Single Base Substitution, Insertion, Deletion, Inversion, and Translocation.

For each of CNVnator, Delly, Pindel, and Samtools files, remember to update parent names accordingly (e.g. Kit-Oryza_sativa_Kitaake_ubi-XA21_sample_A in Delly)

## Requirements
- A [parent file](./parent-coverage-CNVnator-nipponbare.txt) for CNVnator.
- Install `tabix` following [this guidance](https://www.htslib.org/download).
- A `depth file` to fetch gene using `tabix`.
- Refer [the manual document](https://docs.google.com/document/d/1n3o6WMmIaKRRSO86NnCE1AoAldFIiRo4/edit?usp=sharing&ouid=106974761554491462875&rtpof=true&sd=true) for more information.


## How to run

### CNVnator
Produce big deletion (with length of > 1000) mutant gene table.
```
python3 CNVnator.py --file <CNVnator_input_file> --parent <parent-coverage-file> --depthfile <depth_file.gz>
```

### Delly
Produce translocation and inversion mutant gene tables.
```
python3 Delly.py --file <Delly_input_file> --depthfile <depth_file.gz>
```

### Samtools
Produce single base substitution, small deletion, and insertion mutant gene tables.
```
python3 Samtools.py --file <Samtools_input_file>
```

### Pindel
Produce small deletion and insertion mutant gene tables. The samtools input file is to remove duplicates between Pindel and Samtools files.
```
python3 Pindel.py --file <Pindel_input_file> --database <snpEff_database_name> --samtool <Samtool_input_file>
```

For example, I used `Nipponbare` as `snpEff database name`. So:
```
python3 Pindel.py --file <Pindel_input_file> --database Nipponbare --samtool <Samtool_input_file>
```

### Gene description
After running all scripts above, combine them together as a complete mutant gene table. To fetch gene description mapping to gene name, now run `gene_desc.py` with the help of [gene data file](./all.locus_brief_info.7.0) obtained from [RGAP](http://rice.uga.edu/).
```
python3 gene_desc.py --file <mutant_gene_table_file> --gene <gene_file>
```

Example:
```
python3 gene_desc.py --file Combined_Mutant_Gene_Table.tsv --gene all.locus_brief_info.7.0
```
