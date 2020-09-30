import argparse
from Bio import SeqIO
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
import io
import warnings
from collections import defaultdict
import signal


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
- python genbank_groper_sqliteDB.py -rs CP025541.3 CP025541.2 CP009125.1 CP025541.1
Example3 (provides proteins/genes - accepts single input organisms with -r and -c parameters):
- python genbank_groper_sqliteDB.py -a CP025545.1 -s sqlilte.db
Example4 (provides proteins/genes - accepts a tab separted file with -r and -c parameters):
- python genbank_groper_sqliteDB.py -a acc_file.txt -s sqlilte.db
Example5 (provides 16sRNA sequences - accepts multiple input organisms):
- python genbank_groper_sqliteDB.py -rRNA CP025541.3 CP025541.2 CP009125.1 -s sqlilte.db
Example6 (provides partial DNA sequence(s) - accepts a tab separted file with tab separated sequence start, end and strands):
- python genbank_groper_sqliteDB.py -pdna file.bed -s sqlilte.db

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


def organims_name(organism_name):
    if "/" in organism_name:
        organism = organism_name.strip().split("/")[-1]
    else:
        organism = organism_name.strip()

    return organism


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


def seq_check_db(con, organism, dna_file):
    organism_in_table = []
    flag = False
    cursor = con.cursor()
    if os.path.isfile("./mySQLiteDB.db"):
        cursor.execute("SELECT acc From genome_dna WHERE acc=? AND [16sRNA] \
                    IS NULL",
                       (organism,))
        
        organism_in_table = cursor.fetchall()
        con.commit()

    if organism_in_table:
        flag = True
        # return organism_in_table, flag
        return flag
    else:
        cursor.execute("SELECT acc From genome_dna WHERE acc=?",
                       (organism,))
        organism_in_table = cursor.fetchall()
        con.commit()  # committing changes to database
        if organism_in_table:
            print("16sRNA already exist for", organism)
            sys.exit(1)
        else:
            flag = False
            print(organism, "does not exist in DB")
            # return organism_in_table, flag
            return flag


def execute_barrbap(organism, dna_file):
    """determines the 16sRNA sequences using barrnap tool"""
    # barrnap output file name e.g. barrnap.NC_000913
    barrnap_out = cwd + "/barrnap." + organism
    # > /dev/null 2>&1 is to disable stdout from displaying on terminal
    barrnap_cmd = "barrnap " + str(dna_file) + " --quiet --outseq " + barrnap_out + " > /dev/null 2>&1"
    os.system(barrnap_cmd)
    if not os.path.isfile(barrnap_out):
        file_cmd = "touch " + barrnap_out
        os.system(file_cmd)
    return barrnap_out


def hamming_distance(seq1, seq2):
    dist_counter = 0
    seq_len = 0
    diff_len = 0
    if seq2 == seq1:
        dist_counter = 0
    else:
        if len(seq1) ==  len(seq2):
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
    if len(data) == 0:
        cursor.execute("""INSERT into [16sRNA](acc, [16sRNA_Seq])
                    VALUES(?,?)""", (organism, consensus_seq))
        con.commit()                

    cursor.execute("""INSERT into genome_dna(acc, refseq, dna, [16sRNA_Ref])
                VALUES(?,?,?,?)""", (organism, refseq_id, dna_seq, chr_ref))
    con.commit()

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
        consensus_seq =  sequences[0]
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


def insert_genbank(genbank, con, user_acc, chr_ref, chr_flag):
    """
    1 - check if genbank is a file or folder and download file if the data
    is not in the db
    """
    state = ""
    dna_seq = ""
    consensus_seq = ""
    file_arr = []
    if os.path.isfile(genbank):
            file_arr = glob.glob(genbank)
    if os.path.isdir(genbank):
        genbank = str(genbank) + "*"
        file_arr = glob.glob(genbank)
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
                    if len(data) != 0:
                        state += str(gb_record.id) + " - exists!\n"
                        if user_acc in gb_record.id:
                            user_acc = str(gb_record.id)
                        continue

                    rows_genbank = list()
                    rows_proteins = list()
                    rows_genome_dna = list()
                    # todo not well implemented yet
                    rows_genbank.append((gb_record.id, 1, 2, 3, "foo"))
                    gb_id = (gb_record.id).strip()

                    if args.genbank:
                        user_acc = gb_id
                        chr_ref = gb_id

                    if gb_id == user_acc or gb_id == chr_ref:
                        user_acc_flag = True
                        # get data from genbank and rearrange the structure
                        for entry in gb_record.features:
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
                                # creates a temporary file starting with name of organism
                                temp_dna = tempfile.NamedTemporaryFile(prefix=gb_record.id,
                                                                    dir='./')
                                # writes sequence along with its ID in DNA fasta format
                                # encode converts sequence to byte like object
                                # as required by temporary File
                                temp_dna.write((">" + gb_record.id + "\n" +
                                                str(gb_record.seq)).encode())
                                # executes barrnap and returns the output file
                                barrnap_out = execute_barrbap(gb_record.id, temp_dna.name)
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
                                    # print("No 16sRNA sequence found")
                                # dna_seq is being encoded to save space i.e. BLOB entry in DB
                                dna_seq = zlib.compress(str(gb_record.seq).encode('utf-8'))

                                insert_db(con, str(gb_record.id), dna_seq, consensus_seq, chr_ref)
                                # print("SqliteDB is successfuly updated !!!!\n")
                                # remove the redundant/temp file
                                # remove_file()

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
        
        except:
            print("ERROR in extraction and insertion of data")
            os.system("rm -f *.fai barrnap.*")
    
    # remove the output files created by barrnap
    os.system("rm -f *.fai barrnap.*")
    user_acc = chr_ref
    return user_acc


def get_data_from_ncbi(ftp_lnk, chr_ref, user_acc, con, flag, chr_flag):
    ts = calendar.timegm(time.gmtime())
    tmp_download_path = os.getcwd() + "/" + str(ts) + "_TMP"
    os.makedirs(tmp_download_path)
    try:
        if flag:
            os.system("wget -c --retry-connrefused -nv --show-progress --continue\
                    --read-timeout=20 --tries=40 --wait=10 --timeout=15 " +
                    str(ftp_lnk.strip()) + " -P " +
                    str(tmp_download_path))

            file_path = str(tmp_download_path) + "/" + \
                            str(ftp_lnk.split('/')[-1])
            # uncompress downloaded *.gz file
            inF = gzip.open(file_path, 'rb')
            uncompressed_gb = str(file_path) + ".dec"
            outF = open(uncompressed_gb, 'wb')
            outF.write(inF.read())
            inF.close()
            outF.close()
            # integrate data into user's db
            user_acc = insert_genbank(uncompressed_gb, con, user_acc, chr_ref, chr_flag)
        else:
            org_link = '"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nucleotide&id=' + user_acc + '&rettype=gbwithparts&retmode=text"'
            # output file name with location
            download_path = str(tmp_download_path) + "/" + user_acc + ".gbk"

            os.system("wget -c --retry-connrefused -nv --show-progress --continue\
                    --read-timeout=20 --tries=40 --wait=10 --timeout=15 " +
                    str(org_link) + " -O " +
                    str(download_path))

            # integrate data into user's db
            chr_flag = "NA"
            chr_ref = "NA"
            organism_col = "NA"
            user_acc  = insert_genbank(download_path, con, user_acc, chr_ref, chr_flag)

        # clean tmp_download_path
        shutil.rmtree(tmp_download_path)

        return user_acc

    except:
        sys.stderr.write("***FTP/HTTP ERROR: " + user_acc +
                        " CANNOT BE DOWNLOADED! ***")
        shutil.rmtree(tmp_download_path)


def get_data(con, user_id, user_acc, user_center,
             user_range, organism_col, ftp_lnk, organism_arr, chr_ref):
    results = ""
    with con:
        cur = con.cursor()
        left_corner = int(user_center) - int(user_range)
        if left_corner < 0:
            left_corner = 0
        right_corner = int(user_center) + int(user_range)
        try:
            user_acc = user_acc.strip()
            chr_id = chr_ref.strip()
            cur.execute("SELECT acc FROM genome_dna WHERE acc=? OR acc=?", (user_acc, chr_id))
            rows = cur.fetchall()
            con.commit()
            print()
            user_flag = False
            chr_flag = False
            gb_id = 'NA'
            if len(rows) != 0:
                for line in range(len(rows)):
                    if user_acc in rows[line][0]:
                        user_flag = True
                    elif chr_id in rows[line][0]:
                        chr_flag = True
            if not chr_flag or not user_flag:
                if user_acc in organism_col:
                    flag = True
                    user_acc = get_data_from_ncbi(ftp_lnk, chr_ref, user_acc, con, flag, chr_flag)
                    rows = [1]
                    user_acc = user_acc.strip()
                    cur.execute("""SELECT * FROM proteins WHERE acc=?1 and
                                (start_pos<=?2 and end_pos>=?3)""",
                                (user_acc, right_corner, left_corner))
                    rows = cur.fetchall()
                    con.commit()
                    if len(rows) == 0:
                        results += str(user_id) + "\t" + str(user_acc) + "\t" \
                                    + "no annotation" + "\t" + "no annotation" \
                                    "\t" + "no annotation" + "\t" + "no annotation" \
                                    + "\t" + "no annotation" + "\t" + "no annotation" \
                                    + "\n"
                    else:
                        for row in rows:
                            # DESCRIPTION => row[0]=unique_number row[1]=acc
                            # row[2]=locus_name row[3]=gene_name row[4]=start_pos
                            # row[5]=end_pos row[6]=ori row[7]=protein_seq
                            results += str(user_id) + "\t" + str(row[1]) + "\t" + \
                                    str(row[3]) + "\t" + str(row[2]) + "\t" + \
                                    str(row[4]) + "\t" + str(row[5]) + "\t" + \
                                    str(row[6]) + "\t" + \
                                    str(zlib.decompress(row[7]).decode('utf-8') + \
                                    "\n")
                else:
                    flag = False
                    chr_ref =  user_acc
                    organism_col = user_acc
                    ftp_lnk = 'NA'
                    chr_flag = "NA"
                    user_acc = get_data_from_ncbi(ftp_lnk, chr_ref, user_acc, con, flag, chr_flag)
                    user_acc = user_acc.strip()
                    cur.execute("""SELECT * FROM proteins WHERE acc=?1 and
                                (start_pos<=?2 and end_pos>=?3)""",
                                (user_acc, right_corner, left_corner))
                    rows = cur.fetchall()
                    con.commit()
                    if len(rows) == 0:
                        results += str(user_id) + "\t" + str(user_acc) + "\t" \
                                    + "no annotation" + "\t" + "no annotation" \
                                    "\t" + "no annotation" + "\t" + "no annotation" \
                                    + "\t" + "no annotation" + "\t" + "no annotation" \
                                    + "\n"
                    else:
                        for row in rows:
                            # DESCRIPTION => row[0]=unique_number row[1]=acc
                            # row[2]=locus_name row[3]=gene_name row[4]=start_pos
                            # row[5]=end_pos row[6]=ori row[7]=protein_seq
                            results += str(user_id) + "\t" + str(row[1]) + "\t" + \
                                    str(row[3]) + "\t" + str(row[2]) + "\t" + \
                                    str(row[4]) + "\t" + str(row[5]) + "\t" + \
                                    str(row[6]) + "\t" + \
                                    str(zlib.decompress(row[7]).decode('utf-8') + \
                                    "\n")

        except:
            # entry not available
            # todo check if exception can happen and document what's the reason for
            print("Exception in get_data () :", user_id)
    return results


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


def update_db(con, id_container, acc_rrna):
    """update the record in database"""
    final_results = ""
    try:
        with con:
            cur =  con.cursor()
            for i in range(0, len(id_container)):
                # left_corner is start of the range (gene search region)
                left_corner = int(id_container[i][2]) - int(id_container[i][3])
                if left_corner < 0:
                    left_corner = 0
                # right_corner is end of the range (gene search region)
                right_corner = int(id_container[i][2]) + int(id_container[i][3])
                # input organism
                input_param = id_container[i][1].strip()
                # check if the record is already in DB
                cur.execute("SELECT acc FROM genome_dna WHERE acc=?", (input_param,))
                chr_ref = cur.fetchall()
                con.commit()
                # if input is -a paramter
                if acc_rrna:
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
                        # determines ftp link from Master Lookup file if exists
                        organism_col, ftp_lnk, organism_arr, chr_ref = chrom_plas(input_param)
                        # details obtained from master lookup shared to download organims and updaet the DB
                        final_results += get_data(con, id_container[i][0], id_container[i][1],
                                                id_container[i][2], id_container[i][3],
                                                organism_col, ftp_lnk, organism_arr, chr_ref)
                
                # if input is -rRNA paramter
                else:
                    # determines ftp link from Master Lookup file
                    organism_col, ftp_lnk, organism_arr, chr_ref = chrom_plas(input_param)
                    # details obtained from master lookup shared to download organims and updaet the DB
                    final_results += get_data(con, id_container[i][0], id_container[i][1],
                                            id_container[i][2], id_container[i][3],
                                            organism_col, ftp_lnk, organism_arr, chr_ref)
    
    except KeyboardInterrupt:
        print("\nctrl c pressed !!!!!\n")
        print("program safely exited !!!!!!")
        if con:
            con.close()
        os._exit(0)
    except:
        print("ERROR IN update data function !!!!!!!!!!!!\n")    
    return final_results

def find_srRNA_gene(args_accession, args_srRNA, id_container):
    """Return 16sRNA sequences if -rRNA paramter used 
    and with -a parameter returns genes"""
    
    db_path = args.sqlite
    final_results = ''
    if args_srRNA:
        con = lite.connect(db_path)
        con.execute("PRAGMA journal_mode=WAL")
        with con:
            cur = con.cursor()
            try:
                for line in args_srRNA:
                    input_param = line.strip()
                    cur.execute("SELECT [16sRNA_Ref] FROM genome_dna WHERE acc=?", (input_param,))
                    chr_ref = cur.fetchall()
                    con.commit()
                    if chr_ref == []:
                        input_param = line.strip().split(".")[0]
                        cur.execute("SELECT [16sRNA_Ref] FROM genome_dna WHERE acc LIKE ?", (input_param + "._%",))
                        chr_ref = cur.fetchall()
                        con.commit()
                        if chr_ref == []:
                            for data in id_container:
                                input_param = line.strip()
                                if input_param in data[1]:
                                    input_list = [data]
                                    _ = update_db(con, input_list, False)
                                    break
                print()
                for line in args_srRNA:
                    input_param = line.strip().split(".")[0]
                    cur.execute("SELECT [16sRNA_Ref] FROM genome_dna WHERE acc LIKE ?", (input_param + "._%",))
                    chr_ref = cur.fetchall()
                    con.commit()
                    if chr_ref != []:
                        if 'NA' in chr_ref[0][0]:
                            print(">" + str(input_param))
                            print("No 16srRNA found !!!!!!!!!!!!!!\n")
                        else:
                            input_param = line.strip()
                            cur.execute("SELECT * FROM [16sRNA] WHERE acc=?", (chr_ref[0][0],))
                            rows = cur.fetchall()
                            con.commit()
                            for i, line in enumerate(rows):
                                print(">" + str(input_param))
                                print(line[i+1])
            except:
                print("ERROR IN find_sRNA_gene function !!!!!!!!!!!\n")
        con.close()
    elif args_accession:
        con = lite.connect(db_path)
        con.execute("PRAGMA journal_mode=WAL")
        with con:
            final_results = update_db(con, id_container, True)
            if args.output != "":
                # store data into file
                handle = open(args.output, "w")
                handle.write(final_results)
                handle.close()
            else:
                print()
                print(final_results)

        con.close()


def input_check(input_param):
    """input check to ensure Id has correct length"""
    # with NZ length should be minimum 11
    if "NZ_" in input_param and len(input_param)>10:
        pass
    # with NC length should be minimum 9
    elif "NC_" in input_param and len(input_param)>8:
        pass
    # without NZ and NC should be minimum 9
    elif not "NC_" in input_param and not "NZ_" in input_param and len(input_param)>5:
        pass
    else:
        print("please check the organism ID")
        # system exit
        sys.exit()

def dna_extraction(args_extDNA, args_sqlite, ncbi_file_path):
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

            with open(args_extDNA, 'r') as fr:
                for line in fr:
                    positions = defaultdict(list)
                    id, start, stop, strand = line.split()
                    positions[id].append((int(start), int(stop), str(strand)))
                    input_param = line.strip().split()[0]
                    cur.execute("SELECT dna FROM genome_dna WHERE acc=?", (input_param,))
                    dseq = cur.fetchall()
                    con.commit()
                    if dseq != []:
                        data_extracted = seq_extraction(dseq, positions, complement)
                        short_seq_all.update(data_extracted)
                    else:
                        final_results = ""
                        organism_col, ftp_lnk, organism_arr, chr_ref = chrom_plas(input_param)
                        final_results += get_data(con, line[0], input_param,
                                                int(args.position), int(args.range),
                                                organism_col, ftp_lnk, organism_arr, chr_ref)
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
    print()
    for data in short_seq_all:
        print(">" + data )
        print(short_seq_all[data])

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
                        not needed if parameter -a is a file",
                        type=int, default=10)
    parser.add_argument("-r", "--range", help="considered range for extracting\
                        data around the search position; not needed if \
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
    parser.add_argument("-rRNA", "--srRNA", help="Parameter takes\
                        organism(s) separated by space and returns 16srRNA(s)!",\
                        type=str, nargs="*")
    parser.add_argument("-pdna", "--extDNA", help="Parameter takes a\
                        file with ID, start_seq, end_seq, strand in \
                        tab separated format (like bed format) and returns a part of DNA!",\
                        type=str, default="")

    args = parser.parse_args()

    # if sqlite-db-file not exists, then create new db
    if not os.path.exists(args.sqlite):
        init_sqlite_db(args.sqlite)

    # ncbi_file_path = ''
    ncbi_file_path = str(ncbi_folder_name) + "/" + str(ncbi_master_table)

    # check if there are new data available to store in the database
    org_name = args.genbank
    if os.path.exists(org_name):
        # start pickle of data
        con = lite.connect(args.sqlite)
        con.execute("PRAGMA journal_mode=WAL")
        _ = insert_genbank(org_name, con, org_name.strip(), org_name.strip(), 'FALSE')
        con.close()
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
    
    # input_param to avoid execution of -g *.gbk again to get_download()
    input_param = ""
    if args.genbank:
        input_param = (args.genbank).strip()
    if args.srRNA:
        input_param = args.srRNA
    if args.accession:
        input_param = (args.accession).strip()
        input_check(input_param)

    # check input is a string or file 
    id_container = list()
    if os.path.isfile(args.accession):
        handle = open(args.accession)
        for line in handle:
            line = line.rstrip()
            line_arr = line.split("\t")
            line_arr[2] = int(line_arr[2])
            line_arr[3] = int(line_arr[3])
            id_container.append(line_arr)
    elif args.srRNA:
        for line in args.srRNA:
            id_container.append([args.id, line, int(args.position), int(args.range)])
    else:
        # when input is an organism, not a file
        id_container.append([args.id, input_param, int(args.position), int(args.range)])

    if args.refseq:
        for line in args.refseq:
            refseq = find_refseq(ncbi_file_path, line.strip())
            print("refseq ID of " , line , " is " , refseq , "!!!!!!!!!!\n")
        sys.exit()
        
    if args.accession or args.srRNA:
        # determines the genes or rRNA sequences 
        find_srRNA_gene(args.accession, args.srRNA, id_container )

    if args.extDNA:
        # extraction of part of a dna sequence within a given range
        dna_extraction(args.extDNA, args.sqlite, ncbi_file_path)
        
