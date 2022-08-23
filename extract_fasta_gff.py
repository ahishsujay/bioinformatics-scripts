"""
Script that contains 'extract_gene_fasta_gff' function to extract a Gene ID's
sequence and information from a genomic FASTA and GFF as a FASTA and GFF file,
respectively.

Execution:
python3 extract_fasta_gff.py --fasta <path_to_reference_fasta_file> --gff <path_to_reference_gff_file> --geneid <GeneID>

Download test data:
wget http://ftp.ensembl.org/pub/release-107/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.1.fa.gz
wget http://ftp.ensembl.org/pub/release-107/gff3/homo_sapiens/Homo_sapiens.GRCh38.107.chromosome.1.gff3.gz
"""

import argparse
from BCBio import GFF
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def extract_gene_fasta_gff(ref_fasta_file=str, ref_gff_file=str, gene_id=str):
    """
    Function that extracts respective Gene ID's sequence and information from a
    FASTA and GFF file respectively.

        Parameters:
        - ref_fasta_file : Path of genomic FASTA file
        - ref_gff_file : Path of GFF file
        - gene_id : Gene ID to extract sequence and information

        Returns:
        - None
        - <gene_id>_extracted.fasta : Gene ID sequence
        - <gene_id>_extracted.gff3 : Gene ID information
    """
    #Retrieving reference FASTA file as a dict:
    with open(ref_fasta_file, encoding='UTF-8') as ref_fasta_handle:
        ref_fasta_seqs = SeqIO.to_dict(SeqIO.parse(ref_fasta_file, 'fasta'))

    gene_id_flag = False
    #Reading reference GFF file and extracting Gene ID sequence and info:
    with open(ref_gff_file, encoding='UTF-8') as ref_gff_handle:
        for rec in GFF.parse(ref_gff_handle, base_dict=ref_fasta_seqs,
                             target_lines=10000):
            for feature in rec.features:
                #Check if gene_id entered in feature.id:
                if gene_id in feature.id:
                    gene_id_flag = True #Update flag to True
                    gene_seq = Seq(rec.seq[
                                   feature.location.start:feature.location.end])
                    #Checking if on negative strand:
                    if feature.strand == -1:
                        gene_seq = gene_seq.reverse_complement()

                    if 'Name' in feature.qualifiers:
                        fasta_gff_header = feature.qualifiers['Name'][0]
                    elif 'gene_id' in feature.qualifiers:
                        fasta_gff_header = feature.qualifiers['gene_id'][0]
                    #Writing FASTA file:
                    with open(fasta_gff_header + '_extracted.fasta','w',
                              encoding='UTF-8') as gene_fasta_file:
                        SeqIO.write(SeqRecord(gene_seq, fasta_gff_header, '', ''),
                                              gene_fasta_file, 'fasta')
                    #Writing GFF3 file:
                    rec.features = [feature]
                    with open(fasta_gff_header + '_extracted.gff3', 'w',
                              encoding='UTF-8') as gene_gff_file:
                        GFF.write([rec], gene_gff_file)

    if not gene_id_flag:
        print('Entered Gene ID not found.')
        raise Exception('Entered Gene ID not found.')

def main():

    #Argparse code:
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta', help='Enter complete path of genomic FASTA \
                        file.', type=str, required=True)
    parser.add_argument('--gff', help='Enter complete path of GFF file.',
                        type=str, required=True)
    parser.add_argument('--geneid', help='Enter Gene ID.', type=str, required=True)
    args = vars(parser.parse_args())

    #Populating variables:
    reference_fasta_file = args['fasta']
    reference_gff_file = args['gff']
    gene_id = args['geneid']

    #Calling function:
    extract_gene_fasta_gff(reference_fasta_file, reference_gff_file, gene_id)

if __name__ == "__main__":
    main()
