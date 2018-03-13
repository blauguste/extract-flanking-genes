from Bio import SeqIO, SeqFeature, Seq
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
from Bio.Alphabet import IUPAC
import sys

def get_flanking_genes(ref_fn, feature_table, base_fn):
    goi = SeqIO.read(open(ref_fn, 'r'), 'genbank')
    
    print("Parsing %s: %s...", (goi.id, goi.description))

    excel_fn = base_fn + '_flanking_gene_analysis.xlsx'
    igr_fn = base_fn + '_igrs.fa'
    pseudo_fn = base_fn + '_flanking_pseudos.fa'
    prot_fn = base_fn + '_flanking_cds.fa'

    # Throw all the CDS genes in a list
    cdsl = [f for f in goi.features if f.type == 'CDS']

    df = pd.read_csv(feature_table, index_col=0)

    soi_flank = {}
    pseudos = []
    igr_seqs = []
    flanking_prot_seqs = []

    for i, r in df.iterrows():
        
        soi_start = r['start']
        soi_end = r['end']
        soi_name = i

        located_igr = False

        for index, cds in enumerate(cdsl):
            
            if index < (len(cdsl) - 1):
                next_ = cdsl[index + 1]
            prev = cdsl[index - 1]
            
            if soi_start >= (int(cds.location.end) - 10) and soi_end <= (int(next_.location.start) + 10):
                flank_cds = cdsl[index - 2: index + 4]
                soi_flank[soi_name] = {}
                located_igr = True

                # Save the entire IGR seq
                igr_id = goi.id + ':' + str(int(cds.location.end + 1)) + '-' + str(int(next_.location.start))
                igr_desc = i + '_' + 'IGR'
                igr_seq = goi.seq[cds.location.end:next_.location.start]
                igr_seqs.append(SeqRecord(igr_seq, id=igr_desc, description=igr_id))
                
                for p, c in enumerate(flank_cds):
                    
                    # Mark the gene position relative to the IGR of interest
                    if p < 3:
                        pos = p - 3
                    else:
                        pos = p - 2
                    pos = str(pos)
                    
                    # Add the gene name and gene product to the dictionary
                    if 'gene' in c.qualifiers:
                        name = c.qualifiers['gene'][0]
                    else:
                        name = c.qualifiers['locus_tag'][0]
                    if 'product' in c.qualifiers:
                        prod = c.qualifiers['product'][0]
                    else:
                        prod = 'unknown'

                    # Deal with pseudogenes

                    if 'pseudo' in c.qualifiers:
                        prod = 'PSEUDOGENE'
                        cds_seq = goi.seq[c.location.start:c.location.end]
                        pid = c.qualifiers['locus_tag'][0]
                        pseudos.append(SeqRecord(cds_seq, id=pid, description=name + ': ' + prod))
                    else:
                        cds_seq = Seq(c.qualifiers['translation'][0], IUPAC.protein)
                        pid = c.qualifiers['protein_id'][0]
                        flanking_prot_seqs.append(SeqRecord(cds_seq, id=pid, description=name + ': ' + prod))

                    # Save the protein seq
                    
                    soi_flank[soi_name][pos + '_id'] = pid
                    soi_flank[soi_name][pos + '_desc'] = name + ': ' + prod

                break

        if not located_igr:
            print("Hey! FYI, %s wasn't in an intergenic region!", soi_name)
    
    # Write the IGRs and flanking proteins to file
    
    with open(igr_fn, 'w') as igr_out:
        SeqIO.write(igr_seqs, igr_out, 'fasta')

    with open(prot_fn, 'w') as prot_out:
        SeqIO.write(flanking_prot_seqs, prot_out, 'fasta')

    with open(pseudo_fn, 'w') as pseudo_out:
        SeqIO.write(pseudos, pseudo_out, 'fasta')
    
    # Join the input feature table to the gathered info
    
    fg = pd.DataFrame.from_dict(soi_flank, orient='index')
    df = df.join(fg)
    print(df)

    # Write the info about flanking seqs to file

    df.to_excel(excel_fn)

if __name__ == '__main__':
    if len(sys.argv) == 4:
         get_flanking_genes(sys.argv[1], sys.argv[2], sys.argv[3])
    else:
         print("Usage: flanking_gene_extractor.py reference_genome.gb feature_table.csv basename_for_outfiles")
         sys.exit(0)

