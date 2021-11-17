import argparse
from Bio import SeqIO
from multiprocessing import Pool
import re
import zlib
import os.path
import os
import gzip
import glob
import wget
import hashlib
import sys
import calendar
import time
import shutil
from datetime import datetime
import sqlite3 as lite
import numpy as np
import tempfile
import subprocess
import warnings
from collections import defaultdict
import signal
import math


information = """
Usage:
Output:
* default mySQLiteDB.db in working directory.
* DB contains the following 3 tables:-
1. DNA sequence
2. Proteins/Genes
3. 16srRNA
* Provides RefSeq IDs
* Partial DNA sequence

Example1 (updates DB with DNA, Proteins and 16sRNA table):
- python genbank_groper_sqliteDB.py -g NC_000913.faa -s sqlilte.db
Example2 (provides RefSeq IDs - accepts multiple input organisms):
- python genbank_groper_sqliteDB.py -rs CP025541.3 CP025541.2 CP009125.1
Example3 (provides proteins/genes - accepts single input organisms with
            -r and -c parameters):
- python genbank_groper_sqliteDB.py -a CP025545.1 -r 1000 -c 500 -s sqlilte.db
Example4 (provides proteins/genes - accepts a tab separted file ):
- python genbank_groper_sqliteDB.py -a acc_file.txt -s sqlilte.db
Example5 (provides 16sRNA sequences - accepts multiple input organisms):
- python genbank_groper_sqliteDB.py -rRNA CP025541.3 CP025541.2 CP009125.1
            -s sqlilte.db
Example6 (provides partial DNA sequence(s) - accepts a space separted input
            string with @ separated sequence start, end and strands):
- python genbank_groper_sqliteDB.py -pdna CP011246.1@30@50@- CP011247.1@10@20@+
            -s sqlilte.db

Note: Please install the required dependencies
"""

# absolute path of this script
scriptPath = os.path.dirname(os.path.realpath(__file__))

# GLOBAL VARIABLES
genbank_summary_url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/" +\
                    "prokaryotes.txt"
ncbi_folder_name = scriptPath + "/myNCBI-LOOKUP"
ncbi_master_table = "NCBI-MASTER.table"

# disabling ctrl z option
signal.signal(signal.SIGTSTP, signal.SIG_IGN)

# SQLite connection object
con = None

cwd = os.getcwd()
organism = ""
id_16s = ""
sum_of_columns = []


# def organims_name(organism_name):
#     if "/" in organism_name:
#         organism = organism_name.strip().split("/")[-1]
#     else:
#         organism = organism_name.strip()

#     return organism


def barrbap_dict_creation(barrnap_file):
    flag = False
    with open(barrnap_file, 'r') as document:
        barrnap_dict = {}
        keys = ''
        for line in document:
            if line.startswith(">16S"):
                flag = True
                if not line:  # empty line?
                    continue
                # discarding the >16S_rRNA:: from line
                keys = (line.strip().split("::")[1])
            elif flag:
                values = line.strip()
                barrnap_dict[keys] = values
                flag = False
    return barrnap_dict


# def seq_check_db(con, organism, dna_file):
#     organism_in_table = []
#     flag = False
#     cursor = con.cursor()
#     if os.path.isfile("./mySQLiteDB.db"):
#         cursor.execute("SELECT acc From genome_dna WHERE acc=? AND [16sRNA] \
#                     IS NULL",
#                        (organism,))
#         organism_in_table = cursor.fetchall()
#         con.commit()

#     if organism_in_table:
#         flag = True
#         # return organism_in_table, flag
#         return flag
#     else:
#         cursor.execute("SELECT acc From genome_dna WHERE acc=?",
#                        (organism,))
#         organism_in_table = cursor.fetchall()
#         con.commit()  # committing changes to database
#         if organism_in_table:
#             print("16sRNA already exist for", organism)
#             sys.exit(1)
#         else:
#             flag = False
#             print(organism, "does not exist in DB")
#             # return organism_in_table, flag
#             return flag


def execute_barrbap(organism, dna_file, con):
    """determines the 16sRNA sequences using barrnap tool"""
    # barrnap output file name e.g. barrnap.NC_000913
    barrnap_out = cwd + "/barrnap." + organism
    # > /dev/null 2>&1 is to disable stdout from displaying on terminal
    barrnap_cmd = "barrnap " + str(dna_file) + " --quiet --outseq " +\
                    barrnap_out + " > /dev/null 2>&1"
    try:
        # print('Barrnap is RUNNING for ', organism, ' !!!!!')
        if os.stat(dna_file).st_size !=0:
            os.system(barrnap_cmd)
        if not os.path.isfile(barrnap_out):
            file_cmd = "touch " + barrnap_out
            try:
                os.system(file_cmd)
                # subprocess.run([file_cmd], check = True)
            except:
                # print('Ctrl C pressed, program safely exited !!!!!!!!?')
                os.system('rm ' + str(organism) + '*')
                # barrnap_out is not yet created, no removal required !!!
                # os.system('rm -r ' + barrnap_out)
                if con:
                    con.close()
                os._exit(0)
        return barrnap_out
    except:
        print('Ctrl C pressed, program safely exited !!!!!!!!! ###')
        if con:
            con.close()
        os._exit(0)


def hamming_distance(seq1, seq2):
    dist_counter = 0
    seq_len = 0
    diff_len = 0
    if seq2 == seq1:
        dist_counter = 0
    else:
        if len(seq1) == len(seq2):
            seq_len = len(seq1)
            for n in range(seq_len):
                if seq1[n] != seq2[n]:
                    dist_counter += 1
        else:
            if len(seq1) > len(seq2):
                seq_len = len(seq2)
                diff_len = len(seq1) - len(seq2)
            else:
                seq_len = len(seq1)
                diff_len = len(seq2) - len(seq1)
            for n in range(seq_len):
                if seq1[n] != seq2[n]:
                    dist_counter += 1
            dist_counter += diff_len
    return dist_counter


def matrix_fill(dim, arr, sequences):
    try:
        for i in range(dim):
            if i < dim:
                for j in range(dim):
                    if i + j < dim:
                        hamm_dist = hamming_distance(sequences[i],
                                                        sequences[i + j])
                        arr[i, i + j] = hamm_dist
                        arr[i + j, i] = hamm_dist
        return arr
    except:
        return '0'


def sumColumn(dim, matrix_nd, column):
    total = 0
    for row in range(dim):
        total += matrix_nd[row][column]
    return total


def insert_db(con, organism, dna_seq, consensus_seq, chr_ref):
    """Insert record into SqliteDB"""

    # print("insert_db function !!!!")
    
    refseq_id = ""
    if "NC_" in organism or "NZ_" in organism:
        refseq_id = organism
    else:
        refseq_id = find_refseq(ncbi_file_path, organism)
    cursor = con.cursor()
    org_id = chr_ref.strip().split(".")[0]
    cursor.execute("SELECT genome_dna.[16sRNA_Ref], [16sRNA].acc FROM genome_dna, [16sRNA] WHERE genome_dna.[16sRNA_Ref] like ?1 AND [16sRNA].acc like ?1", (org_id + "._%",))
    data = cursor.fetchall()
    con.commit()
    # print(data)
    # print(organism, refseq_id, dna_seq, chr_ref)
    # print(organism, consensus_seq)


    if len(data) == 0:
        cursor.execute("""INSERT into [16sRNA](acc, [16sRNA_Seq])
                    VALUES(?,?)""", (organism, consensus_seq))
        con.commit()

    cursor.execute("""INSERT into genome_dna(acc, refseq, dna, [16sRNA_Ref])
                VALUES(?,?,?,?)""", (organism, refseq_id, dna_seq, chr_ref))
    con.commit()
    # print("Exiting insert_db !!!!")


def consensus_sequence(dim, sequences, sequences_dict, keys):
    # an emtpy matrix created
    consensus_seq = ''
    sum_of_columns = []
    arr = np.zeros((dim, dim), dtype=int)
    arr = matrix_fill(dim, arr, sequences)
    # added to deal with the warning raises due to below comparison
    warnings.simplefilter(action='ignore', category=FutureWarning)
    array_data = np.all((arr == 0), axis=1)
    if np.all(array_data):
        consensus_seq = sequences[0]
    else:
        for column in range(dim):
            sum_of_columns.append(sumColumn(dim, arr, column))
        # consensus sequence index and sequence
        consensus_seq_indx = sum_of_columns.index(min(sum_of_columns))
        consensus_seq = sequences_dict[keys[consensus_seq_indx]]
    return consensus_seq


def md5(fname):
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


def extract_coordinates(record):
    start = 0
    end = 0
    ori = "+"
    range = re.findall("\[(.*?)]", str(record))[0]
    ori = re.findall("\((.*?)\)", str(record))[0]
    range_arr = range.split(":")

    if int(range_arr[0]) <= int(range_arr[1]):
        start = int(range_arr[0]) + 1
        end = int(range_arr[1])
    else:
        end = int(range_arr[0])
        start = int(range_arr[1])
    return start, end, ori


def init_sqlite_db(db_path):
    try:
        con = lite.connect(db_path)
        con.execute("PRAGMA journal_mode=WAL")
        with con:
            cur = con.cursor()
            # genbank generell table
            """
            cur.execute('''CREATE TABLE genbank_generell(acc TEXT PRIMARY KEY,
            #genes_total INTEGER, cds_total INTEGER, genome_size INTEGER,
            #circular_linear TEXT)''')
            """
            # protein table
            cur.execute('''CREATE TABLE proteins(id INTEGER
                    PRIMARY KEY AUTOINCREMENT, acc TEXT,
                    gen_locus TEXT, gene_name TEXT,
                    start_pos INTEGER, end_pos INTEGER,
                    orientation TEXT, sequence TEXT)''')
            # gene name table
            """
            cur.execute('''CREATE TABLE gene_names(id INTEGER
                    PRIMARY KEY AUTOINCREMENT, gen_locus TEXT,
                    gen_name TEXT)''')
            """
            # 16sRNA table
            cur.execute('''CREATE TABLE [16sRNA](
                    acc TEXT PRIMARY KEY,
                    [16sRNA_Seq] TEXT)''')
            # 16S table todo
            # genome DNA content todo
            cur.execute("""CREATE TABLE genome_dna(acc TEXT PRIMARY KEY,
                    refseq TEXT, dna TEXT, [16sRNA_Ref] TEXT)""")
            # create indexes for fast lookup
            cur.execute("CREATE INDEX acc ON proteins (acc)")
            cur.execute("CREATE INDEX prot_start_pos ON proteins (start_pos)")
            cur.execute("CREATE INDEX prot_end_pos ON proteins (end_pos)")
            # cur.execute("CREATE INDEX gen_locus ON gene_names (gen_locus)")
            con.commit()
    except lite.Error:
        sys.stderr.write("SQLite Error " + str(lite.Error) + "\n")
        sys.exit(1)
    finally:
        if con:
            con.close()
    return 0


def insert_genbank(genbank, con, user_acc, chr_ref):
    """
    1 - check if genbank is a file or folder and download file if the data
    is not in the db
    """
    # print("Insert_genbank function ##################")
    # print(genbank)
    state = ""
    dna_seq = ""
    consensus_seq = ""
    file_arr = []
    if os.path.isfile(genbank):
        file_arr = glob.glob(genbank)
    if os.path.isdir(genbank):
        genbank = str(genbank) + "*"
        file_arr = glob.glob(genbank)
    # print(file_arr)
    with con:
        cur = con.cursor()
        user_acc_flag = False
        try:
            for tmp_path in file_arr:
                for gb_record in SeqIO.parse(open(tmp_path, "r"), "genbank"):
                    # check if record exists. If, then take the next entry
                    cur.execute("SELECT acc FROM proteins WHERE acc=?",
                                (gb_record.id,))
                    data = cur.fetchall()
                    con.commit()
                    gb_record_desc = gb_record.description
                    # print(gb_record_desc)
                    # print(len(data))
                    if len(data) != 0:
                        state += str(gb_record.id) + " - exists!\n"
                        # print(state)
                        if user_acc in gb_record.id:
                            user_acc = str(gb_record.id)
                        continue

                    rows_genbank = list()
                    rows_proteins = list()
                    rows_genome_dna = list()
                    
                    gb_id = (gb_record.id).strip()
                    # print(gb_id, user_acc, chr_ref)

                    if args.genbank:
                        user_acc = gb_id
                        chr_ref = gb_id

                    # where we have version of gbk file and not in input ID.
                    if not '.' in user_acc:
                        gb_id = gb_id.strip().split('.')[0]

                    if gb_id == user_acc or gb_id == chr_ref:
                        user_acc_flag = True
                        # get data from genbank and rearrange the structure
                        for entry in gb_record.features:
                            if entry.type == "tRNA" or entry.type == "ncRNA" or entry.type == "rRNA":
                                location = entry.location
                                # check if location is unfinished -> [<0:687](-)
                                # or 7251..>7537
                                if len(re.findall("(<|>)", str(location))) == 0:
                                    # multiple hits for the same sequence
                                    if str(location).startswith("join"):
                                        locs = re.findall("{(.*?)}", str(location))[0]
                                        locs_arr = locs.split(",")
                                        # number of coordinates
                                        for coords in locs_arr:
                                            start, end, ori = extract_coordinates(coords)
                                            # print(start, end, ori)
                                            mito_frag = gb_record.seq[start:end]
                                            # RNA sequence starts with * to differentiate from protein
                                            ext_seq = '*' + mito_frag
                                            # DNA to RNA conversion
                                            rna = str(ext_seq).replace("T", "U")
                                            seq = zlib.compress(rna.
                                            encode('utf-8'))
                                            # print(zlib.decompress(seq).decode('utf-8') )
                                            if "locus_tag" in entry.qualifiers:
                                                gene_loc = entry.qualifiers["locus_tag"][0]
                                            else:
                                                gene_loc = "na"
                                            if "gene" in entry.qualifiers:
                                                gene = entry.qualifiers["gene"][0]
                                            else:
                                                gene = "na"

                                            rows_proteins.append([gb_record.id, gene_loc,
                                                                gene, start, end,
                                                                ori, seq])
                                    else:
                                        start, end, ori = extract_coordinates(location)
                                        # print(start, end, ori)
                                        mito_frag = gb_record.seq[start:end]
                                        ext_seq = '*' + mito_frag
                                        rna = str(ext_seq).replace("T", "U")
                                        seq = zlib.compress(rna.
                                            encode('utf-8'))
                                        # print(zlib.decompress(seq).decode('utf-8') )

                                        if "locus_tag" in entry.qualifiers:
                                            gene_loc = entry.qualifiers["locus_tag"][0]
                                        else:
                                            gene_loc = "na"
                                        if "gene" in entry.qualifiers:
                                            gene = entry.qualifiers["gene"][0]
                                        else:
                                            gene = "na"

                                        rows_proteins.append([gb_record.id, gene_loc, gene,
                                                            start, end, ori, seq])
                            if entry.type == "CDS" and ("translation" in entry.qualifiers):
                                location = entry.location
                                # check if location is unfinished -> [<0:687](-)
                                # or 7251..>7537
                                if len(re.findall("(<|>)", str(location))) == 0:
                                    # multiple hits for the same sequence
                                    if str(location).startswith("join"):
                                        locs = re.findall("{(.*?)}", str(location))[0]
                                        locs_arr = locs.split(",")
                                        # number of coordinates
                                        for coords in locs_arr:
                                            start, end, ori = extract_coordinates(coords)
                                            seq = zlib.compress(
                                                entry.qualifiers["translation"][0].
                                                encode('utf-8'))

                                            if "locus_tag" in entry.qualifiers:
                                                gene_loc = entry.qualifiers["locus_tag"][0]
                                            else:
                                                gene_loc = "na"
                                            if "gene" in entry.qualifiers:
                                                gene = entry.qualifiers["gene"][0]
                                            else:
                                                gene = "na"

                                            rows_proteins.append([gb_record.id, gene_loc,
                                                                gene, start, end,
                                                                ori, seq])
                                    else:
                                        start, end, ori = extract_coordinates(location)
                                        seq = zlib.compress(
                                            entry.qualifiers["translation"][0].
                                            encode('utf-8'))

                                        if "locus_tag" in entry.qualifiers:
                                            gene_loc = entry.qualifiers["locus_tag"][0]
                                        else:
                                            gene_loc = "na"
                                        if "gene" in entry.qualifiers:
                                            gene = entry.qualifiers["gene"][0]
                                        else:
                                            gene = "na"

                                        rows_proteins.append([gb_record.id, gene_loc, gene,
                                                            start, end, ori, seq])
                        try:
                            if len(str(gb_record.seq)) > 0:
                                # print(len(str(gb_record.seq)))
                                # creates a temporary file starting with name of organism
                                temp_dna = tempfile.NamedTemporaryFile(prefix=gb_record.id,
                                                                    dir='./')
                                # writes sequence along with its ID in DNA fasta format
                                # encode converts sequence to byte like object
                                # as required by temporary File
                                temp_dna.write((">" + gb_record.id + "\n" +
                                                str(gb_record.seq)).encode())
                                try:
                                    # print('before barrnap !!!! ')
                                    # executes barrnap and returns the output file
                                    barrnap_out = execute_barrbap(gb_record.id, temp_dna.name, con)
                                    # makes dictionary of 16sRNA seq from barrnap output file
                                    if os.stat(barrnap_out).st_size != 0:
                                        sequences_dict = barrbap_dict_creation(barrnap_out)
                                        # provides keys only in the dictionary
                                        keys = list(sequences_dict.keys())
                                        # extracts sequences from the dictionary
                                        sequences = list(sequences_dict.values())
                                        # earlier below was keys != 0:
                                        if keys != []:
                                            chr_ref = gb_record.id
                                        # array dimensions [dim x dim]
                                        dim = len(keys)
                                        # checks if the detected 16sRNA sequences are two or less
                                        if dim > 2:
                                            """ If 16sRNA sequences are more than two then search for
                                            more similar to each other"""
                                            # print("detected 16sRNA are more than 2")
                                            # determines consensus sequences
                                            consensus_seq = consensus_sequence(dim, sequences,
                                                                            sequences_dict,
                                                                            keys)
                                        elif dim == 2 or dim == 1:
                                            # print("\ndetected 16sRNA are 2 or less !!!!\n")
                                            consensus_seq = sequences[0]

                                    else:
                                        consensus_seq = "No 16sRNA sequence found"
                                        chr_ref = 'NA'
                                        # print("No 16sRNA sequence found")
                                    # dna_seq is being encoded to save space i.e. BLOB entry in DB
                                    dna_seq = zlib.compress(str(gb_record.seq).encode('utf-8'))
                                    insert_db(con, str(gb_record.id), dna_seq, consensus_seq, chr_ref)
                                    # print("SqliteDB is successfuly updated !!!!\n")
                                    # remove the redundant/temp file
                                    # remove_file()
                                except:
                                    print('Ctrl C pressed, program safely exited !!!!!!')
                                    os.system('rm ' + str(gb_record.id) + '*')
                                    if con:
                                        con.close()
                                    os._exit(0)

                            else:
                                dna_seq = "NA"
                                consensus_seq = "NA"
                                chr_ref = "NA"
                                insert_db(con, str(gb_record.id), dna_seq, consensus_seq, chr_ref)
                            rows_genome_dna.append([gb_record.id, dna_seq])
                            try:
                                cur.executemany("""INSERT INTO proteins (acc, gen_locus,
                                                gene_name, start_pos, end_pos, orientation,
                                                sequence) VALUES (?,?,?,?,?,?,?)""",
                                                rows_proteins)
                                con.commit()
                            except con.IntegrityError:
                                sys.stderr.write("ERROR: ID already exists in PRIMARY KEY\
                                                column " + str(con.IntegrityError) + "\n")
                                continue
                        except:
                            print("len(str(gb_record.seq)) in exception:" , len(str(gb_record.seq)))
                            print("something went wrong with extraction of data")
                            # ctrl c pressed
                            if con:
                                con.close()
                            os._exit(0)
        
        except:
            print("ERROR in extraction and insertion of data")
            os.system("rm -rf *.fai barrnap.* ./Temp_dir")
            # ctrl c pressed
            if con:
                con.close()
            os._exit(0)

    # remove the output files created by barrnap
    os.system("rm -f *.fai barrnap.*")
    user_acc = chr_ref
    # print("Exiting insert_genbank !!!")
    return user_acc

def get_data_from_ncbi(ftp_lnk, chr_ref, user_acc):
    # print("!!!!!!!!! get_data_from_ncbi_1 function !!!!")
    # print('user_acc: ', user_acc)
    ts = 'Temp_dir/' + user_acc
    tmp_download_path = os.getcwd() + "/" + str(ts) + "_TMP"
    # print(ftp_lnk)
    # os.makedirs(tmp_download_path)
    try:
        if ftp_lnk:
            try:
                # without creating dir, wget creates by iteself with -P parameter
                # os.makedirs(tmp_download_path)
                ftp_cmd = "wget -c --retry-connrefused -nv --show-progress --continue --read-timeout=20 --tries=40 --wait=10 --timeout=15 " + str(ftp_lnk.strip()) + " -P " + str(tmp_download_path)
                # print(ftp_cmd)
                os.system("wget -c --retry-connrefused -nv --show-progress --continue\
                        --read-timeout=20 --tries=40 --wait=10 --timeout=15 " +
                        str(ftp_lnk.strip()) + " -P " +
                        str(tmp_download_path))

                file_path = str(tmp_download_path) + "/" + \
                                str(ftp_lnk.split('/')[-1])
                # uncompress downloaded *.gz file
                inF = gzip.open(file_path, 'rb')
                uncompressed_gb = str(tmp_download_path) + "/" + str(user_acc) + ".gbk"
                outF = open(uncompressed_gb, 'wb')
                outF.write(inF.read())
                inF.close()
                outF.close()
            except:
                print('ERROR in download !!!!!!!!!!!!', tmp_download_path)
                if os.path.isdir(tmp_download_path):
                    # print('file not removed !!!!')
                    os.system('rm -r ' + tmp_download_path)    
                if con:
                    print('DB connection closed !!!!')
                    con.close()
                    os._exit()
                print('EXITING FTP !!!!!!!!!!!')
                #sys.exit(0)

        else:
            # print('FTP link does not exists !!!! ')
            # double quotes removed from link due to subprocess.run
            org_link = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=' + user_acc + '&rettype=gbwithparts&retmode=text'
            # output file name with location
            try:
                os.makedirs(tmp_download_path)
                os.system('ls ' + tmp_download_path)
            except:
                print(user_acc, ': directory exists !!!! ')

            download_path = str(tmp_download_path) + "/" + user_acc + ".gbk"
            cmd_wget = ['wget', '-c', '--retry-connrefused', '--waitretry=6', '--retry-on-http-error=429', '-nv', '--show-progress', '--continue', '--read-timeout=20', '--tries=40', '--wait=10', '--timeout=15', str(org_link),     '-O', str(download_path)]
            # cmd_wget = ['curl ']
            if os.path.exists(download_path):
                last_line = os.popen('tail -n 2 ' + download_path).read()
                # print('line: ', last_line)
                if not '//' in last_line:
                    # print(last_line)
                    os.system('rm -r ' + download_path)
                    # p = subprocess.run(cmd_wget, check=True)
                    subprocess.call(cmd_wget)
            else:
                try:
                    print()
                    p = subprocess.run(cmd_wget, check=True)
                    # subprocess.call(cmd_wget)
                except:
                    print('ERROR in eutils download !!!!!!!!!!!!', tmp_download_path)
                    p.send_signal(signal.SIGINT)
                    if os.path.isdir(tmp_download_path):
                        # print('file removed !!!!')
                        os.system('rm -r ' + tmp_download_path)    
                    if con:
                        # print('DB connection closed !!!!')
                        con.close()
                    os._exit(0)        
    except:
        # print('Termination get_data_from_ncbi_1 END !!!!')
        # if os.path.isdir(tmp_download_path):
            # print('file removed !!!!')
            # os.system('rm -r ' + tmp_download_path)
        # os.system('rm -r wget-log*')
        if con:
            con.close()
        os._exit(0)


def process_NCBI_lookup(in_file, directory, out_file):
    result = list()
    handle = open(in_file, "r", encoding='ascii', errors="replace")
    for line in handle:
        line = line.rstrip()
        tmp_arr = line.split("\t")
        # choosing only the entries with chromosomes. 
        if "chromosome" in line and not line.startswith("#") and (len(tmp_arr) == 23):
            try:
                # Replicons[8] + FTP Path [20]
                # extract identifiers; only ACC
                # splits the replicons into chromosomes and plasmids
                tmp_identifiers_arr = tmp_arr[8].split("; ")
                file_identifiers = ""
                chromosome = []
                plasmid = []
                for line in tmp_identifiers_arr:
                    # splitting the replicons to later write to a file
                    replicons = line.strip().split(":")
                    # if list item has chromosome
                    if "chromosome" in replicons[0].strip().split()[0]:
                        chromosome.append(replicons[1].strip())
                    # if list item has plasmid
                    if "plasmid" in replicons[0].strip().split()[0]:
                        plasmid.append(replicons[1].strip())
                # if no plasmid or chromosome entry, then put "-"
                if len(plasmid) == 0:
                    plasmid.append("-")
                if len(chromosome) == 0:
                    chromosome.append("-")
                # prefix of file to download
                missing_part = tmp_arr[20].split("/")[-1]
                file_identifiers = ','.join(chromosome) + "\t" + ','.join(plasmid) +  "\t" + str(tmp_arr[20]) + "/" + str(missing_part) + "_genomic.gbff.gz" + "\n"
                result.append(file_identifiers)
            except:
                pass

    # write output
    out_path = str(directory) + "/" + str(out_file)
    handle = open(out_path, "w")
    for entry in result:
        handle.write(entry)
    handle.close()
    # compute MD5SUM prokaryotes.txt
    md5value = md5(in_file)
    md5file = str(in_file) + ".md5"
    handle = open(md5file, "w")
    handle.write(md5value)
    handle.close()
    # clean prokayrotes.txt file
    os.remove(in_file)


def find_refseq(ncbi_file_path, input_param):
    """finds the refseq_id of an accession in master lookup table"""
    handle = open(ncbi_file_path, "r")
    refseq = ""
    for line in handle:
        line = line.rstrip()
        line_arr = line.split("\t")
        line_chr = re.split(',', line_arr[0])
        for i, val in enumerate(line_chr):
            if input_param in val:
                # taking 1st id as chr_ref bcz all of them has same 16sRNA seq
                if "/" in val:
                    refseq = re.split('/',val)[0]
                else: 
                    refseq = "NA"
        line_plas = re.split(',', line_arr[1])
        for i, val in enumerate(line_plas): 
            if input_param in val:
                if "/" in val:
                    refseq = (re.split('/',(val)))[0]
                else: 
                    refseq = "NA"
    if refseq == "":
        refseq = "NA"
    
    handle.close()
    return refseq


def seq_extraction(dseq, positions, complement):
    short_seq_all = dict()
    dna_seq = ''
    for row in dseq:
        dna_seq = zlib.decompress(row[0]).decode('utf-8')
    for name in positions:
        for (start, stop, strand) in positions[name]:
            short_seq = str(dna_seq)[start-1:stop]
            # name = str(name) + "_" + str(start) + "_" + str(stop)
            if strand == "+":
                name = str(name) + ":" + str(start) + "-" + str(stop)
                short_seq_all[name] = short_seq
            elif strand == "-":
                name = str(name) + ":c" + str(stop) + "-" + str(start)
                short_seq_neg = short_seq
                bases = list(short_seq_neg) 
                bases = reversed([complement.get(base,base) for base in bases])
                bases = ''.join(bases)
                short_seq_all[name] = bases
    return short_seq_all

def chrom_plas(input_param):
    # User starts a search - check if the ACC is in the database or not.
    # If not, download the genbank from ncbi
    handle = open(ncbi_file_path, "r")
    organism_col = ''
    ftp_lnk = ''
    chr_ref = ''
    organism_arr = ''
    refseq = ''
    for line in handle:
        line = line.rstrip()
        line_arr = line.split("\t")
        line_chr = re.split(',', line_arr[0])
        for i, val in enumerate(line_chr):
            # ensures refseq IDs are downloaded from NCBI directly, 
            # not from the ftp link
            if "/" in val:
                lookup_id = val.strip().split("/")[1]
            else:
                lookup_id = val.strip().split()[0]
            if input_param in lookup_id:
                # # we can use index to ensure the first organism is choosen as the chr_ref
                # array is here to make the list of all relevant chromosmes
                # if we don't find 16srRNA in any, we can have a look in another chromosome
                organism_arr = ','.join(line_chr)
                organism_col = lookup_id
                ftp_lnk = line_arr[2].strip()
                chr_ref = organism_col

        line_plas = re.split(',', line_arr[1])
        for i, val in enumerate(line_plas): 
            if "/" in val:
                lookup_id = val.strip().split("/")[1]
            else:
                lookup_id = val.strip().split()[0]

            if input_param in lookup_id:
                organism_arr = ','.join(line_plas)
                organism_col = line_plas[i]
                ftp_lnk = line_arr[2].strip()
                if "/" in organism_col:
                    chr_ref = (re.split('/',(line_chr[0])))[1]
                else: 
                    chr_ref = line_chr[0].strip()
    handle.close()
    return organism_col, ftp_lnk, organism_arr, chr_ref


def download_organisms(input_param):

    try:
        organism_col, ftp_lnk, organism_arr, chr_ref = chrom_plas(input_param)
        if ftp_lnk:
            chr_ref_dict[organism_col] = chr_ref
            get_data_from_ncbi(ftp_lnk, chr_ref, organism_col)
        else:
            # print('ELSE called in download !!!!')
            ftp_lnk = ''
            organism_col = input_param
            chr_ref = 'NA'
            chr_ref_dict[organism_col] = chr_ref
            get_data_from_ncbi(ftp_lnk, chr_ref, organism_col)
            # print('after exit !!!!')
            
    except:
        print(' In except state')
        if con:
            con.close()
            sys.exit()
        os._exit(0)
        


def data_check_in_db(args_rRNA):
    download_orgs = []
    db_path = args.sqlite

    try:
        con = lite.connect(db_path)
        con.execute("PRAGMA journal_mode=WAL")
        with con:
            cur = con.cursor()
            for line in args_rRNA:
                input_param = line.strip()
                if not '.' in input_param:
                    cur.execute("SELECT [16sRNA_Ref] FROM genome_dna WHERE acc LIKE ?", (input_param + "._%",))
                else:
                    cur.execute("SELECT [16sRNA_Ref] FROM genome_dna WHERE acc=?", (input_param,))
                chr_ref = cur.fetchall()
                con.commit()
                if chr_ref == []:
                    download_orgs.append(input_param)
                else:
                    pass
        con.close()

    except:
        print("\nctrl c pressed !!!!!\n")
        print("find_rRNA_gene safely exited !!!!!!")
        # os.system("rm -rf *.fai barrnap.* ./Temp_dir")
        if con:
            con.close()
        sys.exit(0)

    return download_orgs


def data_insertion_after_download(db_path):
    dir_list = os.listdir('./Temp_dir')
    con = lite.connect(db_path)
    con.execute("PRAGMA journal_mode=WAL")
    with con:
        cur = con.cursor()
        if dir_list:
            for line in dir_list:
                gen_bank = './Temp_dir/' + line + '/' + str(line.replace('_TMP', '')) + '.gbk'
                # changed chr_ref_dict[line.strip('_TMP')] to gen_bank since we dont need this info for this calculation.
                _ = insert_genbank(gen_bank, con, line.replace('_TMP', ''), gen_bank)
                cmd_rm = 'rm -r ./Temp_dir/' + line
                # print(cmd_rm)
                os.system(cmd_rm)
    con.close()


def find_rRNA_gene(args_rRNA):
    db_path = args.sqlite
    download_orgs = data_check_in_db(args_rRNA)

    no_of_orgs = args.norgs
    range_list = math.ceil(len(download_orgs) / no_of_orgs)
    orgs_list = []
    k = 0

    for i in range(range_list):
        orgs_list.append(download_orgs[k:k + no_of_orgs])
        k = k + no_of_orgs

    try:
        for lst in orgs_list:
            # print(lst)
            with Pool(args.cores) as p:
                p.map(download_organisms, lst)
            data_insertion_after_download(db_path)
    except:
        print('download interupted in pool !!!!')
        os._exit(0)

    for input_param in download_orgs:
        organism_col, ftp_lnk, organism_arr, chr_ref = chrom_plas(input_param)
        if ftp_lnk:
            # print(organism_col, ftp_lnk, organism_arr, chr_ref)
            chr_ref_dict[organism_col] = chr_ref
        else:
            chr_ref = 'NA'
            chr_ref_dict[input_param] = chr_ref

    con = lite.connect(db_path)
    con.execute("PRAGMA journal_mode=WAL")
    with con:
        cur = con.cursor()
        print()
        for line in args_rRNA:
            input_param = line.strip().split(".")[0]
            cur.execute("SELECT [16sRNA_Ref] FROM genome_dna WHERE acc LIKE ?", (input_param + "._%",))
            chr_ref = cur.fetchall()
            con.commit()
            if chr_ref != []:
                if 'NA' in chr_ref[0][0]:
                    print(">" + str(input_param))
                    print("No 16srRNA found !!!!!!!!!!!!!!")
                else:
                    input_param = line.strip()
                    cur.execute("SELECT * FROM [16sRNA] WHERE acc=?", (chr_ref[0][0],))
                    rows = cur.fetchall()
                    con.commit()
                    for i, line in enumerate(rows):
                        print(">" + str(input_param))
                        print(line[i+1])

    con.close()


def find_accession(id_container ):
    db_path = args.sqlite
    accesion_list = []
    # print(id_container)

    for i in id_container:
        accesion_list.append(i[1])

    download_orgs = data_check_in_db(accesion_list)

    no_of_orgs = args.norgs
    range_list = math.ceil(len(download_orgs) / no_of_orgs)
    orgs_list = []
    k = 0

    for i in range(range_list):
        orgs_list.append(download_orgs[k:k + no_of_orgs])
        k = k + no_of_orgs
    try:
        for lst in orgs_list:
            # print(lst)
            with Pool(args.cores) as p:
                p.map(download_organisms, lst)
            data_insertion_after_download(db_path)
    except:
        print('download interupted in pool !!!!')
        os._exit(0)
    
    final_results = ""
    try:
        con = lite.connect(db_path)
        con.execute("PRAGMA journal_mode=WAL")
        with con:
            cur =  con.cursor()
            for i in range(0, len(id_container)):
                # left_corner is start of the range (gene search region)
                left_corner = int(id_container[i][2]) - int(id_container[i][3])
                # print('left_corner: ', left_corner)
                if left_corner < 0:
                    left_corner = 0
                # right_corner is end of the range (gene search region)
                right_corner = int(id_container[i][2]) + int(id_container[i][3])
                # print('right_corner: ', right_corner)
                # input organism
                input_param = id_container[i][1].strip()
                # print('input_param: ', input_param)
                # check if the record is already in DB
                if not '.' in input_param:
                    cur.execute("SELECT acc FROM genome_dna WHERE acc LIKE ?", (input_param + "._%",))
                else:
                    cur.execute("SELECT acc FROM genome_dna WHERE acc=?", (input_param,))
                chr_ref = cur.fetchall()
                con.commit()
                # if record is found
                if chr_ref != []:
                    # chr_ref[0][0] is an input_organism
                    # extract the sequence within specified range
                    cur.execute("""SELECT * FROM proteins WHERE acc=?1 and
                                (start_pos<=?2 and end_pos>=?3)""",
                                (chr_ref[0][0], right_corner, left_corner))
                    rows = cur.fetchall()
                    con.commit()
                    # if no gene exists within the given range
                    if len(rows) == 0:
                        final_results += str(id_container[i][0]) + "\t" + str(id_container[i][1]) + "\t" \
                                    + "no annotation" + "\t" + "no annotation" \
                                    "\t" + "no annotation" + "\t" + "no annotation" \
                                    + "\t" + "no annotation" + "\t" + "no annotation" \
                                    + "\n"
                    # if there exist genes in organism DB
                    else:
                        for row in rows:
                            # DESCRIPTION => row[0]=unique_number row[1]=acc
                            # row[2]=locus_name row[3]=gene_name row[4]=start_pos
                            # row[5]=end_pos row[6]=ori row[7]=protein_seq
                            final_results += str(id_container[i][0]) + "\t" + str(row[1]) + "\t" + \
                                    str(row[3]) + "\t" + str(row[2]) + "\t" + \
                                    str(row[4]) + "\t" + str(row[5]) + "\t" + \
                                    str(row[6]) + "\t" + \
                                    str(zlib.decompress(row[7]).decode('utf-8') + \
                                    "\n")
                else:
                    print(input_param + ' not found in databse !!! ')
                    
    except KeyboardInterrupt:
        print("\nctrl c pressed !!!!!\n")
        print("program safely exited !!!!!!")
        if con:
            con.close()
        os._exit(0)
    except:
        print("ERROR IN update data function !!!!!!!!!!!!\n")    

    if args.output != "":
        # store data into file
        handle = open(args.output, "w")
        handle.write(final_results)
        handle.close()
    else:
        print()
        print(final_results)


def dna_extraction(args_extDNA, args_sqlite, ncbi_file_path):

    db_path = args.sqlite
    accesion_list = []

    for i in args_extDNA:
        accesion_list.append(i.split('@')[0])
    # print(accesion_list)
    download_orgs = data_check_in_db(accesion_list)

    no_of_orgs = args.norgs
    range_list = math.ceil(len(download_orgs) / no_of_orgs)
    orgs_list = []
    k = 0
    for i in range(range_list):
        orgs_list.append(download_orgs[k:k + no_of_orgs])
        k = k + no_of_orgs
    try:
        for lst in orgs_list:
            # print(args.cores)
            with Pool(args.cores) as p:
                p.map(download_organisms, lst)
            data_insertion_after_download(db_path)
    except:
        print('download interupted in pool !!!!')
        os._exit(0)
    
    # data_insertion_after_download(db_path)
    
    flag = False
    con = lite.connect(args.sqlite)
    con.execute("PRAGMA journal_mode=WAL")
    try:
        with con:
            cur = con.cursor()
            dna_flag = False
            # positions = defaultdict(list)
            short_seq_all = dict()
            complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
            for line in args_extDNA:
                positions = defaultdict(list)
                id, start, stop, strand = line.split("@")
                positions[id].append((int(start), int(stop), str(strand)))
                input_param = id.strip().split()[0]
                if not '.' in input_param:
                    cur.execute("SELECT dna FROM genome_dna WHERE acc LIKE ?", (input_param + "._%",))

                else:
                    cur.execute("SELECT dna FROM genome_dna WHERE acc=?", (input_param,))
                dseq = cur.fetchall()
                con.commit()
                if dseq != []:
                    data_extracted = seq_extraction(dseq, positions, complement)
                    short_seq_all.update(data_extracted)
                else:
                    print("no sequence found for ", id)

    except lite.IntegrityError:
        print("pdna could not be performed successfully !!!!!!!!!!!!\n")
    except:
        if con:
            con.close()
    print()
    for data in short_seq_all:
        print(">" + data )
        print(short_seq_all[data])
    print()
    # close db connection
    con.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--genbank", help="genbank file or path to \
                        genbank files", type=str, default="")
    parser.add_argument("-s", "--sqlite", help="Path to SQLite DB",
                        type=str, default="./mySQLiteDB.db")
    parser.add_argument("-a", "--accession", help="accession id\
                        e.g. CP021219.1", type=str, default="")
    parser.add_argument("-c", "--position", help="search position; \
                        not required if parameter -a is a file",
                        type=int, default=10)
    parser.add_argument("-r", "--range", help="considered range for extracting\
                        data around the search position; not required if \
                        parameter -a is a file", type=int, default=3000)
    parser.add_argument("-i", "--id", help="Set unique identifier. Used\
                        in batch-processing!", type=str, default="na")
    parser.add_argument("-o", "--output", help="Specify file name for\
                        the results. Otherwise, std is used",
                        type=str, default="")
    parser.add_argument("-u", "--update", help="auto updater checks every\
                        30 days for an updated NCBI lookup table \
                        (on-30 , every 30 days) -> values: on-30,\
                        off, man", type=str, default="on-60")
    parser.add_argument("-x", "--glassgo", help="Parameter takes a finite\
                        string and returns a string with dna sequences\
                        - replaces blastdbcmd!", type=str, default="")

    parser.add_argument("-rs", "--refseq", help="Parameter takes \
                        organism(s) separated by space and returns a refseq ID!",\
                        type=str, nargs="*")
    parser.add_argument("-rRNA", "--rRNA", help="Parameter takes\
                        organism(s) separated by space and returns 16srRNA(s)!",\
                        type=str, nargs="*")
    parser.add_argument("-pdna", "--extDNA", help="Parameter takes a\
                        space separted (if multiple organims)  input string with ID, start_seq, end_seq and strand (+/-); \
                        separated by @ (e.g. CP011246.1@30@50@+ CP011248.1@10@30@-) and returns a part of DNA accordingly!",\
                        type=str, nargs="*")
    parser.add_argument("-nr", "--norgs", help="No of organisms\
                        to be downaloded \
                        at once; default: 40", type=int, default=40)
    parser.add_argument("-cr", "--cores", help="No of cores\
                        to be used for downalod \
                        at once; default: 3", type=int, default=3)

    args = parser.parse_args()

    # chromosome reference records 
    chr_ref_dict = dict()
    # if sqlite-db-file not exists, then create new db
    if not os.path.exists(args.sqlite):
        init_sqlite_db(args.sqlite)

    # ncbi_file_path = ''
    ncbi_file_path = str(ncbi_folder_name) + "/" + str(ncbi_master_table)
    
    # check if there are new data available to store in the database
    gbk_file = args.genbank
    if gbk_file:
        if os.path.exists(gbk_file):
            print("genbank input provided !!!!")
            # start pickle of data
            con = lite.connect(args.sqlite)
            con.execute("PRAGMA journal_mode=WAL")
            # old call with False parameter
            # _ = insert_genbank(gbk_file, con, gbk_file.strip(), gbk_file.strip(), 'FALSE')
            _ = insert_genbank(gbk_file, con, gbk_file.strip(), gbk_file.strip())
            con.close()
        else:
            print('file does not exist !!!!!!')
    # check if ncbi lookup-file needs to be updated or not!
    update_mode = args.update.split("-")

    if update_mode[0] == "on" and os.path.exists(ncbi_folder_name):
        # get date of current version
        if os.path.isdir(ncbi_file_path):
            ncbi_file_date = os.path.getctime(ncbi_file_path)
            diff = datetime.now() - datetime.fromtimestamp(ncbi_file_date)
            days = re.findall("(.*) days", str(diff))
            if len(days) > 0:
                if int(days[0]) >= int(update_mode[1]):
                    # update -> delete content in folder /myNCBI-LOOKUP/ !!!!
                    to_delete_path = str(ncbi_folder_name) + "/"
                    shutil.rmtree(to_delete_path)
    elif update_mode[0] == "man" and os.path.exists(ncbi_folder_name):
        # update -> delete content in folder /myNCBI-LOOKUP/ !!!!
        to_delete_path = str(ncbi_folder_name) + "/"
        shutil.rmtree(to_delete_path)

    # check if ncbi lookup table exist,
    # otherwise download the table from ncbi's ftp server
    if os.path.isdir(ncbi_folder_name):
        ncbi_master_table_path = str(ncbi_folder_name) + "/" + str(ncbi_master_table)
        if os.path.isfile(ncbi_master_table_path):
            pass
        else:
            # get data and reprocess the file
            prok_file = str(ncbi_folder_name) + "/prokaryotes.txt"
            file_path = ''
            if not os.path.isfile(prok_file):
                file_path = wget.download(
                            genbank_summary_url, out=ncbi_folder_name
                            )
            else:
                file_path = prok_file
            process_NCBI_lookup(file_path, ncbi_folder_name, ncbi_master_table)
            print()
    else:
        os.mkdir(ncbi_folder_name)
        # get data and reprocess the file
        file_path = wget.download(genbank_summary_url, out=ncbi_folder_name)
        process_NCBI_lookup(file_path, ncbi_folder_name, ncbi_master_table)
        print()

    input_param = ""
    
    id_container = list()
    if args.accession:
        input_param = (args.accession).strip()
        # checks input is a string or file
        if os.path.isfile(args.accession):
            handle = open(args.accession)
            for line in handle:
                line = line.rstrip()
                line_arr = line.split("\t")
                line_arr[2] = int(line_arr[2])
                line_arr[3] = int(line_arr[3])
                id_container.append(line_arr)
        else:
            # when input is an organism, not a file
            id_container.append([args.id, input_param, int(args.position), int(args.range)])
        # determines the genes
        find_accession(id_container)
        # find_srRNA_gene(args.accession, args.rRNA, id_container)

    elif args.rRNA:
        # print('args.srRNA: ', args.rRNA)
        for line in args.rRNA:
            id_container.append([args.id, line, int(args.position), int(args.range)])
        # determines the rRNA sequences
        # print(args.rRNA)
        find_rRNA_gene(args.rRNA)
        # find_srRNA_gene(args.accession, args.srRNA, id_container)

    elif args.refseq:
        for line in args.refseq:
            refseq = find_refseq(ncbi_file_path, line.strip())
            print("refseq ID of ", line, " is ", refseq, "!!!!!!!!!!\n")

    elif args.extDNA:
        # extraction of part of a dna sequence within a given range
        dna_extraction(args.extDNA, args.sqlite, ncbi_file_path)
