from cyvcf2 import VCF, Writer
from gtfparse import read_gtf
import argparse
import os

# Parse GTF
def read_gtf_in(gtf):
    return(read_gtf(gtf))

# Modify VCF to include GENES_IN, GENES_200KB & GENE_NEAREST in INFO field
def modify_vcf(vcf_reader, gtf_df, out_vcf):
    vcf_reader.add_info_to_header({'ID': 'GENES_IN', 'Description': 'overlapping gene', 'Type':'String', 'Number': '1'})
    vcf_reader.add_info_to_header({'ID': 'GENES_200KB', 'Description': 'Genes within 200KB', 'Type':'String', 'Number': '1'})
    vcf_reader.add_info_to_header({'ID': 'GENE_NEAREST', 'Description': 'Nearest gene', 'Type':'String', 'Number': '1'})

    w = Writer(out_vcf, vcf_reader)

    for variant in vcf_reader:
        # Overlapping genes & nearest gene
        overlap_df = gtf_df[(variant.POS>=gtf_df['start']) & (variant.POS<=gtf_df['end'])]
        if overlap_df.empty:
            lowest = 0
            variant.INFO["GENES_IN"] = "."
            for i, end_pos in enumerate(gtf_df['end']):
                tmp_low = variant.POS - end_pos
                if tmp_low < lowest:
                    lowest = tmp_low
                    gene_nearest_id = gtf_df.loc[i, "gene_id"]
            variant.INFO["GENE_NEAREST"] = gene_nearest_id      
        else:    
            variant.INFO["GENES_IN"] = ",".join(set(list(overlap_df.gene_id)))
            variant.INFO["GENE_NEAREST"] = "."

        # Genes within 200kb
        variant_start = variant.POS-200000
        variant_end = variant.POS+200000
        genes_200kb = gtf_df[(gtf_df['start']>=variant_start) & (gtf_df['end']<=variant_end)]
        if genes_200kb.empty:
            variant.INFO["GENES_200KB"] = "."
        else:   
            variant.INFO["GENES_200KB"] = ",".join(set(list(genes_200kb.gene_id)))
 
        w.write_record(variant)
    
    w.close() 
    vcf_reader.close()

if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(prog='fetch_overlap',
                                    usage='%(prog)s --variant_file <variants_file.vcf.gz> --reference_file <reference_file.gtf.gz> --out_file <out_file.vcf>',
                                    description='Fetch overlapping genes, ±200KB genes & nearest genes')
    parser.add_argument('--variant_file', type=str, required=True)
    parser.add_argument('--reference_file', type=str, required=True)
    parser.add_argument('--out_file', type=str, required=True)
    
    args = parser.parse_args()
    
    gtf_df = read_gtf_in(args.reference_file)
    vcf_reader = VCF(args.variant_file)
    modify_vcf(vcf_reader, gtf_df, args.out_file)

    if os.path.exists(args.out_file):
        print("VCF successfully written")