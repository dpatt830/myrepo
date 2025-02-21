from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import Bio.SeqRecord

email = "dpatterson3@luc.edu"
accession = "NC_006273.2"


def extract_cds_to_fasta(accession, email):
    '''With valid email and accession number of genome, retrieve CDS seq,
        protein ID, and create fasta file'''

    # Entrez query initiation
    Entrez.email = email
    stream = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
    record = SeqIO.read(stream, format = "genbank")
    
    # List comp. to hold all features that are labeled as 'CDS'
    cds_content = [feature for feature in record.features if feature.type == "CDS"]
    
    # Initialize first cds seq
    first_cds = cds_content[0]

    # looking thru each cds content/feature extract the cds seq and protein ID header
    seq = (first_cds.extract(record.seq))
    protein_seq = seq.translate(to_stop=True)
    prot_id = first_cds.qualifiers.get("protein_id", ["unknown_protien"])[0]

    # creating fasta record with translated seq and protein identifier/header
    fasta_record = SeqRecord(protein_seq, id=prot_id)

    # Initializing CDS feature length and fasta file name
    num_CDS_feat = len(cds_content)
    fasta_name = f"{accession}-HCMV-transcriptome.fasta"

    # Using SeqIO to write a fasta file of the cds seq w/ identifier
    SeqIO.write(fasta_record, fasta_name, 'fasta')
    
    # Returning file name and num of cds features
    return(fasta_name, num_CDS_feat)
    

# Example usage:
print((extract_cds_to_fasta(accession, email)))