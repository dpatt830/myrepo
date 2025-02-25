from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import Bio.SeqRecord
import pandas as pd
import os
import subprocess
import argparse
import sys

#function to parse command line arguments
def check_arg(args=None):
    '''Command line arguments/flags'''

    parser = argparse.ArgumentParser(description='Python wrapper for transcriptome analysis.')
    parser.add_argument('-i', '--input',
		help = 'Path to input sequence data directory',
		required = 'True'
		)	
    parser.add_argument("-o", '--output',
        help= 'Path to output directory',
        required=True
        )
    parser.add_argument('-e', '--email',
		help = 'Email Address',
		required = 'True'
		)
    parser.add_argument('-a', '--accession',
		help = 'NCBI Genome Accession Number',
		required = 'True'
		)
    parser.add_argument('-l', '--logFileName',
		help = 'Desired Log File Name',
		required = 'True'
		)
    parser.add_argument('-s', '--subfamily',
		help = 'Subfamily',
		required = 'True'
		)
    return parser.parse_args(args)

#retrieve command line arguments and assign to variables
args = check_arg(sys.argv[1:])
inputPath = args.input
outputPath = args.output
email = args.email
accession = args.accession
logPath = args.logFileName
subfamily = args.subfamily

def outputDirectory(dir_path):
    '''Creating an output directory'''
    os.makedirs(dir_path)
    os.chdir(dir_path)

def buildLogFile(logFileName):
    '''Create out Log File/Output File'''
    log_file = open(f"logFileName", "w")

def extract_cds_to_fasta(accession, email):
    '''With valid email and accession number of genome, retrieve CDS seq,
        protein ID, and create fasta file'''

    # Entrez query initiation
    Entrez.email = email
    stream = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="text")
    record = SeqIO.read(stream, format = "genbank")
    
    # List comp. to hold all features that are labeled as 'CDS'
    cds_content = [feature for feature in record.features if feature.type == "CDS"]
    
    # initialize a lsit to hold all cds content
    fasta_file = []

    for cds in cds_content:
        # looking thru each cds content/feature extract the cds seq and protein ID header
        seq = (cds.extract(record.seq))
        prot_id = cds.qualifiers.get("protein_id", ["unknown_protien"])[0] # getting protein ID

        # cappending to fasta list
        fasta_file.append(SeqRecord(seq, id=prot_id))

    # Initializing CDS feature length and fasta file name
    num_CDS_feat = len(cds_content)
    fasta_name = f"{accession}-HCMV-transcriptome.fasta"

    # Using SeqIO to write a fasta file of the cds seq w/ identifier
    SeqIO.write(fasta_file, fasta_name, 'fasta')
    
    # Returning file name and num of cds features
    return(fasta_name, num_CDS_feat)
    

def buildIndex(transcript_file, num_cds, log_path):
    '''With our transcriptome fasta file and num of cds features of the reference genome'''

    index_name = 'index.idx'  # Initialize index file name

    kallisto_command = f'kallisto index -i {index_name} {transcript_file}' # Initialize kallisto command

    os.system(kallisto_command) # run command

    # Writing output to log file
    with open(log_path, "w") as log:
        log.write(f"The HCMV genome (NC_006273.2) has {num_cds} CDS."+"\n")
        log.write("\n")

    return index_name

def kallistoTPM(index, logPath, inputPath):
    '''Initializing Kallisto script'''

    # Create a kallisto output directory specific to each donor fastq pair
    kallistoDir = "kallisto"
    os.makedirs(kallistoDir)

    # iterating thru each donor data folder
    for donorFolder in os.listdir(f"../{inputPath}/"):

        # finding full fast file name
        for file_name in os.listdir(f'../{inputPath}/{donorFolder}/'):

            if len(file_name)>10:  # SRA fastq file name

                # accessing the SRA run accession ID
                SRArun = file_name[:-8] # everything before _#.fastq to get base name

        # getting forward and reverse fastq files
        fwd = f"../{inputPath}/{donorFolder}/{SRArun}_1.fastq"
        rev = f"../{inputPath}/{donorFolder}/{SRArun}_2.fastq"

        # Create a kallisto output directory specific to each donor fastq pair
        kallistoOutput = f"kallisto/{donorFolder}"
        os.makedirs(kallistoOutput)

        # creating kallisto command
        # change index file command to {}
        kallisto_command = f'time kallisto quant -i {index} -o {kallistoOutput} -b 30 -t 2 {fwd} {rev}'
        os.system(kallisto_command)
    
    # preparing the header for tpm output
    with open(logPath, 'a') as log:
        log.write("sample    condition    min_tpm    med_tpm    mean_tpm    max_tpm")

    return kallistoDir

def quantifyTPM(logPath, kallistoOutput):
    '''Writing tpm output to log file'''

    # Iterating each folder in the results folder of our kallisto folder
    for donor_folder in os.listdir(f'./{kallistoOutput}/'):
        
        # create a data frame of each abundance.tsv file & record the specific values within each donor
        df = pd.read_csv(f"./{kallistoOutput}/{donor_folder}/abundance.tsv", sep='\t')
        min_tpm = df["tpm"].min()
        med_tpm = df["tpm"].median()
        mean_tpm = df["tpm"].mean()
        max_tpm = df["tpm"].max()

        # Writing output to log file
        with open(logPath, "a") as log:
            log.write("\n")
            log.write(f"{donor_folder[:10]}    {donor_folder[11:]}    {min_tpm}    {med_tpm}    {mean_tpm}    {max_tpm}")

def metaData(kallistoOutput):
    '''With kallisto result folders, create a meta data txt file for Sleuth'''

    # create metadata txt file
    with open("metadata.txt", 'w') as meta:
        meta.write("sample condition path")
        meta.write("\n")

        # writing to the metadata txt file our sample, cond., and path
        for donor in os.listdir(f'./{kallistoOutput}/'):
            meta.write(f'{donor} {donor[-4:]} {kallistoOutput}/{donor}')
            meta.write("\n")

def sleuthRun(logPath):
    '''Running sleuth script and writing to logfile'''

    # run our sleuth.R script
    os.system("Rscript ../sleuth.R")

    # with the fdr results txt, make a df
    df = pd.read_csv("fdr05_results.txt",delimiter=" ")

    # keep only the columns we want to write to the log file
    columns_to_keep = ['target_id', "test_stat", "pval", "qval"]
    df_select = df[columns_to_keep]

    # write df results to log file
    with open(logPath, "a") as log:
        log.write("\n")
        log.write("\n")
        df_select.to_csv(log, header=True, index=False, sep="\t")

def HCMVgenome(accession, email):
    '''Fetching the full HCMV genome from NCBI and saving it as a fasta file'''

    Entrez.email = email

    stream = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
    record = SeqIO.read(stream, format="fasta")

    fasta_output = "HCMV_genome.fasta"
    # Save to a FASTA file
    SeqIO.write(record, fasta_output, "fasta")

    return fasta_output

def bowtie(fastaFile, inputPath):
    '''With the full HCMV genome, build index using Bowtie 2,
        and see how many seqs map to it of our fastq seqs'''

    # create new output directory for bowtie2 index and for mapping
    bowtie2_output = "bowtie2-Output"
    bowtie2_index = "bowtie2-Index"
    os.makedirs(bowtie2_output)
    os.makedirs(bowtie2_index)

    # initialize command for building index
    bowtieIndexCommand = f"bowtie2-build {fastaFile} ./{bowtie2_index}/HCMV"
    
    # run bowtie 2 command to build index with full genome
    os.system(bowtieIndexCommand)

    # iterating thru each donor data folder; similar code to Kallisto run
    for donorFolder in os.listdir(f"../{inputPath}/"):

        # finding full fast file name
        for file_name in os.listdir(f'../{inputPath}/{donorFolder}/'):

            if len(file_name)>10:  # SRA fastq file name

                # accessing the SRA run accession ID
                SRArun = file_name[:-8] # everything before _#.fastq to get base name

        # getting forward and reverse fastq files
        fwd = f"../{inputPath}/{donorFolder}/{SRArun}_1.fastq"
        rev = f"../{inputPath}/{donorFolder}/{SRArun}_2.fastq"

        # creating bowtie2 mapping command
        bowtieMapCommand = f'bowtie2 --quiet -x ./{bowtie2_index}/HCMV -1 {fwd} -2 {rev} \
                            -S ./{bowtie2_output}/{donorFolder}-map.sam --al-conc ./{bowtie2_output}/{donorFolder}-mapped_reads.fastq'
        os.system(bowtieMapCommand)
    
def bowtieLogFile(logPath, inputPath):
    '''With the SAM files output from our bowtie2 command, 
        write mapped reads to log file'''
    
    # initialize lists to hold both read pairs for before and after Bowtie
    before = []
    after = []

    # iterating our donor folders with the fastq folders
    for donorFolder in os.listdir(f"../{inputPath}/"):
        # finding full fast file name
        for file_name in os.listdir(f'../{inputPath}/{donorFolder}/'):

            if len(file_name)>10:  # SRA fastq file name

                # accessing the SRA run accession ID
                SRArun = file_name[:-8] # everything before _#.fastq to get base name
        
        # getting forward fastq file of before and after bowtie mapping
        fwd = f"../{inputPath}/{donorFolder}/{SRArun}_1.fastq"
        fwd_mapped = f"./bowtie2-Output/{donorFolder}-mapped_reads.1.fastq"
        
        # using subprocess to hold the output of the number of read pairs before for each donor
        reads_before = (int(subprocess.check_output(f"wc -l < {fwd}", shell=True).strip())) / 4

        # using samtools to count the num of read pairs in our sam files and appending to after list
        mapped_reads = (int(subprocess.check_output(f"wc -l < {fwd_mapped}", shell=True).strip())) / 4

        # appending results to before and after list
        before.append(reads_before)
        after.append(mapped_reads)
    
    # writing to log file our read pairs from the lists
    with open(logPath, 'a') as log:
        log.write("\n")
        log.write(f'Donor 1 (2dpi) had {int(before[0])} read pairs before Bowtie2 filtering and {int(after[0])} read pairs after.'+"\n")
        log.write(f'Donor 1 (6dpi) had {int(before[1])} read pairs before Bowtie2 filtering and {int(after[1])} read pairs after.'+"\n")
        log.write(f'Donor 3 (2dpi) had {int(before[2])} read pairs before Bowtie2 filtering and {int(after[2])} read pairs after.'+"\n")
        log.write(f'Donor 3 (6dpi) had {int(before[3])} read pairs before Bowtie2 filtering and {int(after[3])} read pairs after.'+"\n")
        log.write("\n")

def spades(logPath, inputPath):
    '''Running the spades program through the command line,
        using the mapped fastq files'''
    
    # initializing spades directory
    spades_dir = "./spades"

    # Initialize donor set
    donor_set = set()

    # iterating thru all the files in the bowtie2-Output directory
    for donor in os.listdir(f"../{inputPath}/"):
        
        # initialize a single Donor/Patient
        donorSingle = donor[:6]

        donor_set.add(donorSingle)

    # Iterating over set that will hold {"donor1", "donor3"}
    for patient in donor_set:
        # initialize forward and reverse mapped fastq file pairs
        fwd_2dpi = f"{patient}-2dpi-mapped_reads.1.fastq"
        rev_2dpi = f"{patient}-2dpi-mapped_reads.2.fastq"
        fwd_6dpi = f"{patient}-6dpi-mapped_reads.1.fastq"
        rev_6dpi = f"{patient}-6dpi-mapped_reads.2.fastq"

        # create new donor directory for spades outputs for each donor
        donor_output_dir = f"{spades_dir}/{patient}"
        os.makedirs(donor_output_dir)

        # create spades command
        spades_command = f"spades.py -k 77 -t 2 --only-assembler --pe-1 1 ./bowtie2-Output/{fwd_2dpi} --pe-2 1 ./bowtie2-Output/{rev_2dpi} --pe-1 2 ./bowtie2-Output/{fwd_6dpi} --pe-2 2 ./bowtie2-Output/{rev_6dpi} -o {donor_output_dir}"

        # run spades command
        os.system(spades_command)

        # write spades commands to log file
        with open(logPath, 'a') as log:
            log.write(spades_command+"\n")

def contigs():
    '''Get the longest contigs from each spades contig.fasta file'''

    # make blast directory
    blast_dir = "./blast/"

    # run command
    os.makedirs(blast_dir)

    # for each donor prefix name
    for donor in os.listdir("./spades/"):

        # creating command to grab the first seq in contigs.fasta file for each donor and create fasta file in blast directory
        long_contig_command = "awk '/^>/ {if (seq) exit; seq=1} {print}' ./spades/"+donor+"/contigs.fasta > ./blast/"+donor+"-long_contig.fasta"

        # run command
        os.system(long_contig_command)


def blast(subFamily):
    '''With a given sub family and our long contig fasta files,
        run blast commands'''
    
    # initializing commands to download and build the sub family database thru blast
    download_genome_command = f"datasets download virus genome taxon {subFamily}"
    unzip_command = "unzip ncbi_dataset.zip"
    make_db_command = f"makeblastdb -in ncbi_dataset/data/genomic.fna -out ./blast/{subFamily} -title {subFamily} -dbtype nucl"

    # run db commands
    os.system(download_genome_command)
    os.system(unzip_command)
    os.system(make_db_command)

    # iterating donor prefix names
    for donor in os.listdir("./spades/"):

        # make directories for each donor in blast directory
        os.makedirs(f"./blast/{donor}")

        # initialize blast command
        blast_command = f'blastn -query ./blast/{donor}-long_contig.fasta -db ./blast/{subFamily} \
                        -out ./blast/{donor}/{subFamily}_blastn_results.csv \
                        -outfmt "10 sacc pident length qstart qend sstart send bitscore evalue stitle"'
        
        os.system(blast_command)
    
def blastLogFile(logPath):
    '''Write results from BLAST to log file'''

    # open log file 
    with open(logPath, "a") as log:
        log.write("\n")

        # Initialize csv's of both blsatn csv results
        csv1 = f"./blast/donor1/Betaherpesvirinae_blastn_results.csv"
        csv3 = f"./blast/donor3/Betaherpesvirinae_blastn_results.csv"

	# begin to write BLAST hits
        log.write("Donor1:"+"\n")

        # read lines of BLAST hits for Donor 1
        with open(csv1, 'r') as f1:
            blast_results1 = f1.readlines()

            
            # Sometimes there are not 10 or more hits #
            # If there are write only the top ten hits
            if (len(blast_results1)) >= 10:
                for i in range(10):
                    line1 = blast_results1[i].split(",")
                    string1 = "    ".join(line1)
                    log.write(string1)

            # If there is less than ten, write them all
            else:
                for i in blast_results1:
                    line1 = i.split(',')
                    string1 = "    ".join(line1)
                    log.write(string1)
		
        # writing to log file for nice format
        log.write("\n")
        log.write("Donor3:"+"\n")

	    # read lines of BLAST hits for Donor 3                
        with open(csv3, 'r') as f3:
            blast_results3 = f3.readlines()

            
            # Sometimes there are not 10 or more hits #
            # If there are write only the top ten hits
            if (len(blast_results3)) >= 10:
                for k in range(10):
                    line3 = blast_results3[k].split(",")
                    string3 = "    ".join(line3)
                    log.write(string3)

            # If there is less than ten, write them all
            else:
                for k in blast_results3:
                    line3 = k.split(',')
                    string3 = "    ".join(line3)
                    log.write(string3)

# Initialize functions
outputDirectory(outputPath)
buildLogFile(logPath)
transcriptomeFasta, numCDS = extract_cds_to_fasta(accession, email)
indexName = buildIndex(transcriptomeFasta, numCDS, logPath)
kallitstoPath = kallistoTPM(indexName, logPath, inputPath)
quantifyTPM(logPath, kallitstoPath)
metaData(kallitstoPath)
sleuthRun(logPath)
genomeFasta = HCMVgenome(accession, email)
bowtie(genomeFasta, inputPath)
bowtieLogFile(logPath, inputPath)
spades(logPath, inputPath)
contigs()
blast(subfamily)
blastLogFile(logPath)
