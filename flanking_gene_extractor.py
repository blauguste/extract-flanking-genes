from Bio import SeqIO, SeqFeature, Seq
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Entrez
import pandas as pd
from Bio.Alphabet import IUPAC
import sys
import os

def get_gene_info(goi, feat):
    """
    Input a genbank feature to return its name, id, product, and sequence
    """
    
    if 'gene' in feat.qualifiers:
        name = feat.qualifiers['gene'][0]
    else:
        name = feat.qualifiers['locus_tag'][0]
    if 'product' in feat.qualifiers:
        prod = feat.qualifiers['product'][0]
    else:
        prod = 'unknown'

    if 'pseudo' in feat.qualifiers: # Deal with pseudogenes
        prod = 'PSEUDOGENE'
        cds_seq = goi.seq[feat.location.start:feat.location.end]
        pid = feat.qualifiers['locus_tag'][0]

    else:
        cds_seq = Seq(feat.qualifiers['translation'][0], IUPAC.protein)
        pid = feat.qualifiers['protein_id'][0]

    feat_dict = {'name': name, 'prod': prod, 'cds_seq': cds_seq, 'pid': pid}

    return feat_dict

def overlap(start1, end1, start2, end2):
    """
    Does the range (start1, end1) overlap with (start2, end2)?
    """
    return end1 >= start2 and end2 >= start1

def get_flanking_genes(feature_table, base_fn, feat_name_field, gb_acc_field, email):
    
    Entrez.email = email

    # Load the table containing the coordinates and accessions for sRNAs of interest 
    df = pd.read_pickle(feature_table)

    for a in df[gb_acc_field].unique():
        
        if os.path.isfile(a + '.gb'):
            foi = open(a + '.gb', 'r')
        else:
            foi = Entrez.efetch(db='nucleotide', id=a, rettype='gbwithparts', retmode='text')
        
        with foi as handle:
            print('parsing genbank %s' % a)
    
            goi = SeqIO.read(handle, 'genbank')

            print("Parsing %s: %s..." % (goi.id, goi.description))

            excel_fn = base_fn + '_flanking_gene_analysis.xlsx'
            igr_fn = base_fn + '_igrs.fa'
            pseudo_fn = base_fn + '_flanking_pseudos.fa'
            prot_fn = base_fn + '_flanking_cds.faa'

            # Throw all the CDS genes in a list
            cdsl = [f for f in goi.features if f.type == 'CDS']

            soi_flank = {}
            pseudos = []
            igr_seqs = []
            flanking_prot_seqs = []

            for i, r in df[df[gb_acc_field] == a].iterrows():
                
                print(r[feat_name_field])

                soi_start = r['start']
                soi_end = r['end']
                soi_name = r[feat_name_field] + '_' + goi.id

                for index, cds in enumerate(cdsl):
                    
                    if str(cds.location).startswith('join'): # Skip compound locations
                        continue

                    if index < (len(cdsl) - 1):
                        next_ = cdsl[index + 1]
                    
                    if soi_start >= (int(cds.location.end) - 10) and soi_end <= (int(next_.location.start) + 10):
                        print("IGR found.")
                        flank_cds = cdsl[index - 2: index + 4]
                        soi_flank[soi_name] = {}

                        # Save the entire IGR seq
                        igr_id = goi.id + ':' + str(int(cds.location.end + 1)) + '-' + str(int(next_.location.start))
                        igr_desc = soi_name + '_' + 'IGR'
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
                            
                            print(type(c))
                            cinfo = get_gene_info(goi, c)

                            if 'pseudo' in c.qualifiers:
                                pseudos.append(SeqRecord(cinfo['cds_seq'], id=cinfo['pid'], description=cinfo['name'] + ': ' + cinfo['prod']))
                            else:
                                flanking_prot_seqs.append(SeqRecord(cinfo['cds_seq'], id=cinfo['pid'], description=cinfo['name'] + ': ' + cinfo['prod']))

                            # Save the protein seq
                            
                            soi_flank[soi_name][pos + '_id'] = cinfo['pid']
                            soi_flank[soi_name][pos + '_desc'] = cinfo['name'] + ': ' + cinfo['prod']

                        break

                    elif overlap(soi_start, soi_end, int(cds.location.start), int(cds.location.end)):
                        
                        print("overlapping gene found.")
                        # Don't count the overlapping gene as part of the "flanking genes"
                        soi_flank[soi_name] = {}

                        cinfo = get_gene_info(goi, cds)
                        if 'pseudo' in cds.qualifiers:
                            pseudos.append(SeqRecord(cinfo['cds_seq'], id=cinfo['pid'], description=cinfo['name'] + ': ' + cinfo['prod']))
                        else:
                            flanking_prot_seqs.append(SeqRecord(cinfo['cds_seq'], id=cinfo['pid'], description=cinfo['name'] + ': ' + cinfo['prod']))

                        soi_flank[soi_name]['overlap_id'] = cinfo['pid']
                        soi_flank[soi_name]['overlap_desc'] = cinfo['name'] + ': ' + cinfo['prod']
                        
                        flank_cds = []
                        flank_cds.extend(cdsl[index - 3: index])
                        flank_cds.extend(cdsl[index + 1: index + 4])

                        # Save the entire IGR seq not including the overlapping gene
                        igr_id = goi.id + ':' + str(int(flank_cds[2].location.end + 1)) + '-' + str(int(flank_cds[3].location.start))
                        igr_desc = soi_name + '_' + 'IGR'
                        igr_seq = goi.seq[flank_cds[2].location.end:flank_cds[3].location.start]
                        igr_seqs.append(SeqRecord(igr_seq, id=igr_desc, description=igr_id))
                        
                        for p, c in enumerate(flank_cds):
                            
                            # Mark the gene position relative to the IGR of interest
                            if p < 3:
                                pos = p - 3
                            else:
                                pos = p - 2
                            pos = str(pos)
                            
                            # Add the gene name and gene product to the dictionary
                            cinfo = get_gene_info(goi, c)

                            if 'pseudo' in c.qualifiers:
                                pseudos.append(SeqRecord(cinfo['cds_seq'], id=cinfo['pid'], description=cinfo['name'] + ': ' + cinfo['prod']))
                            else:
                                flanking_prot_seqs.append(SeqRecord(cinfo['cds_seq'], id=cinfo['pid'], description=cinfo['name'] + ': ' + cinfo['prod']))

                            # Save the protein seq
                            
                            soi_flank[soi_name][pos + '_id'] = cinfo['pid']
                            soi_flank[soi_name][pos + '_desc'] = cinfo['name'] + ': ' + cinfo['prod']
                        
                        break

        
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

    # Write the info about flanking seqs to file

    df.to_excel(excel_fn)

if __name__ == '__main__':
    if len(sys.argv) == 6:
         get_flanking_genes(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
    else:
         print("Usage: flanking_gene_extractor.py feature_table.csv basename_for_outfiles feat_name_field field_containing_gb_accession email\n \
            Feature table should be a pickled dataframe\n \
            Start/end loci should correspond directly to the input reference genome!")
         sys.exit(0)

