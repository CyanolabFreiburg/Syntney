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
# import time

information = """
Determines the 16sRNA for provided organism and stores
in mySQLiteDB.db (sqlite database).

Usage:
Input:
* Organism DNA file (Accession_No.faa)

Output:
* mySQLiteDB.db in working directory.

Example1:
- python 16sRNA_detection.py -g NC_000913.faa

Note: Please install barrnap before running this script.
Conda Installation:
- conda install -c bioconda -c conda-forge barrnap
"""

# absolute path of this script
scriptPath = os.path.dirname(os.path.realpath(__file__))

# GLOBAL VARIABLES
genbank_summary_url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/" +\
                    "prokaryotes.txt"
ncbi_folder_name = scriptPath + "/myNCBI-LOOKUP"
ncbi_master_table = "NCBI-MASTER.table"
# SQLite related
con = None


cwd = os.getcwd()
organism = ""
id_16s = ""
sum_of_columns = []


def organims_name(organism_name):
    if "/" in organism_name:
        organism = organism_name.strip().split("/")[-1].split(".")[0]
        # print(organism + " organism with / executed")
    else:
        organism = organism_name.strip().split(".")[0]
        # print(organism + " without / executed")

    return organism


def barrbap_dict_creation(barrnap_file):
    flag = False
    # print("In barrbap_dict\n")
    with open(barrnap_file, 'r') as document:
        barrnap_dict = {}
        keys = ''
        for line in document:
            if line.startswith(">16S"):
                # line = line.split()
                flag = True
                if not line:  # empty line?
                    continue
                # discarding the >16S_rRNA:: from line
                keys = (line.strip().split("::")[1])
                # print(keys)
            elif flag:
                values = line.strip()
                # print(values)
                barrnap_dict[keys] = values
                # chrm_details[key_name_strand_pos].append(chrm_pos_val)
                flag = False
    return barrnap_dict


def seq_check_db(con, organism, dna_file):
    organism_in_table = []
    flag = False
    # conn = lite.connect("./mySQLiteDB.db")
    cursor = con.cursor()
    if os.path.isfile("./mySQLiteDB.db"):
        cursor.execute("SELECT acc From genome_dna WHERE acc=? AND [16sRNA] \
                    IS NULL",
                       (organism,))
        
        organism_in_table = cursor.fetchall()
        con.commit()

    if organism_in_table:
        flag = True
        # con.close()  # closing database connection
        # return organism_in_table, flag
        return flag
    else:
        cursor.execute("SELECT acc From genome_dna WHERE acc=?",
                       (organism,))
        organism_in_table = cursor.fetchall()
        con.commit()  # committing changes to database
        # con.close()  # closing database connection
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
    barrnap_out = "barrnap." + organism
    # > /dev/null 2>&1 is to disable stdout from displaying on terminal
    barrnap_cmd = "barrnap " + str(dna_file) + " --quiet --outseq " + cwd + \
                  "/" + barrnap_out + " > /dev/null 2>&1"

    # barrnap_cmd = "barrnap " + str(dna_file) + " --quiet --outseq " + cwd + \
                #   "/" + barrnap_out
    os.system(barrnap_cmd)
    # print("\nbarrnap execution completed !!!")

    return barrnap_out


def hamming_distance(seq1, seq2):
    dist_counter = 0
    seq_len = 0
    diff_len = 0
    if seq2 == seq1:
        # print("BOTH SEQ ARE EQUAL !!!!!!!!!!!\n")
        dist_counter = 0
    else:
        # print("ELSE len of seq1", len(seq1), len(seq2))
        if len(seq1) ==  len(seq2):
            seq_len = len(seq1)
            for n in range(seq_len):
            # print(n)
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
            # print(n)
                if seq1[n] != seq2[n]:
                    dist_counter += 1
            dist_counter += diff_len
    return dist_counter


def matrix_fill(dim, arr, sequences):
    try:
        # print(dim)
        for i in range(dim):
            # print("i : ", i)
            if i < dim:
                for j in range(dim):
                    # print("i + j < dim: ", i + j < dim)
                    if i + j < dim:
                        # print(sequences[i], ":", sequences[i + j] )
                        hamm_dist = hamming_distance(sequences[i],
                                                    sequences[i + j])
                        # print("hammding distance: ", hamm_dist)
                        arr[i, i + j] = hamm_dist
                        arr[i + j, i] = hamm_dist
        # print(arr)
        return arr
    except:
        return '0'
    # return arr

def sumColumn(dim, matrix_nd, column):
    total = 0
    # print(dim, column)
    # print(matrix_nd)    
    for row in range(dim):
        total += matrix_nd[row][column]
    # print("total :" , total)
    return total


def insert_db(con, organism, dna_seq, consensus_seq, chr_ref):
    """Insert record into SqliteDB"""
    # print("Insert record into SqliteDB")
    refseq_id = ""
    if "NC_" in organism or "NZ_" in organism:
        refseq_id = organism
    else:
        refseq_id = find_refseq(ncbi_file_path, organism)
    # print("\nrefseq: ", refseq_id)
    cursor = con.cursor()
    cursor.execute("SELECT genome_dna.[16sRNA_Ref], [16sRNA].acc FROM genome_dna, [16sRNA] WHERE genome_dna.[16sRNA_Ref]=?1 AND [16sRNA].acc=?1", (chr_ref,))
    data = cursor.fetchall()
    con.commit()
    # print("data : ", (data))
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
    # print("array before fill : ", arr)
    # print("parameter to matrix-fill are dim, arr, sequences", dim, len(arr), len(sequences))
    arr = matrix_fill(dim, arr, sequences)
    # print("arr : ", arr)
    # added to deal with the warning raises due to below comparison
    warnings.simplefilter(action='ignore', category=FutureWarning)
    # print(type(arr[0][1]), arr[0][1])
    array_data = np.all((arr == 0), axis=1)
    # print("array_data", array_data)
    # print("sequences[0] : ", sequences[0])
    if np.all(array_data):
        # print(" its a zero matrix !!!!!!!!!!!\n")
        consensus_seq =  sequences[0]
    else:
        # print("else called for sum of column")
        # print(arr)
        for column in range(dim):
            sum_of_columns.append(sumColumn(dim, arr, column))

        # consensus sequence index and sequence
        # print("sum of columns :" , sum_of_columns)
        consensus_seq_indx = sum_of_columns.index(min(sum_of_columns))
        consensus_seq = sequences_dict[keys[consensus_seq_indx]]
        # print("sum of columns index :" , consensus_seq_indx)
    
    # print(consensus_seq)
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
    # print("args.sqlite in insert genbank", args.sqlite)
    # con = lite.connect(args.sqlite)
    if os.path.isfile(genbank):
            file_arr = glob.glob(genbank)
    if os.path.isdir(genbank):
        genbank = str(genbank) + "*"
        file_arr = glob.glob(genbank)
    with con:
        # print("Fayyaz INSERT_GENBANK function is called !!!!!!!!!!!\n")
        # print("genbank, con, user_acc, chr_ref, chr_flag : ", genbank, con, user_acc, chr_ref, chr_flag)
        # print("file_arr : ", file_arr)
        cur = con.cursor()
        # print(file_arr)
        user_acc_flag = False
        try:

            for tmp_path in file_arr:
                # print("temp_Path: " + tmp_path)
                for gb_record in SeqIO.parse(open(tmp_path, "r"), "genbank"):
                    # print("gb_record.id: " + gb_record.id)
                    # check if record exists. If, then take the next entry
                    cur.execute("SELECT acc FROM proteins WHERE acc=?",
                                (gb_record.id,))
                    data = cur.fetchall()
                    con.commit()
                    # print(data)
                    # print("gb_record.name:", gb_record.name)
                    # print("gb_record.desc:", gb_record.description)
                    gb_record_desc = gb_record.description
                    # print(type(gb_record_desc))
                    # if not "plasmids" in gb_record_desc:
                        # print("YES PLASMID FOUND !!!!!!!!!!!!!!\n")

                    # if not "plasmid" in gb_record_desc:
                        # print("NOT PLASMID FOUND !!!!!!!!!!!!!!\n")
                    # if not "plasmid" in gb_record.description: # or not "contigs" in gb_record.description or not "scaffold" in gb_record.description:
                    #     chr_ref = gb_record.id
                        # print("\nITS IS NOT A PLASMID !!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")
                    if len(data) != 0:
                        state += str(gb_record.id) + " - exists!\n"
                        # Fayyaz to make compatible both naming conventions referring to same gbk
                        # only name difference.
                        # print("chr_flag, organism_col : ", chr_flag, user_acc, gb_record.id)
                        if user_acc in gb_record.id:
                            # print("user acc and gb record", user_acc, gb_record.id)
                            user_acc = str(gb_record.id)
                        # print(state)
                        continue

                    rows_genbank = list()
                    rows_proteins = list()
                    rows_genome_dna = list()
                    # todo not well implemented yet
                    rows_genbank.append((gb_record.id, 1, 2, 3, "foo"))
                    gb_id = (gb_record.id).strip().split(".")[0]

                    # print("gb_id == user_acc, gb_id, user_acc: ",  gb_id == user_acc, gb_id, user_acc)
                    # print(gb_id == chr_ref.split(".")[0], gb_id, chr_ref.split(".")[0]  )

                    # if gb_id == user_acc or gb_id == chr_ref.split(".")[0]:
                    if gb_id == user_acc or gb_id == chr_ref.split(".")[0]:
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
                                # print("len(str(gb_record.seq)) :", len(str(gb_record.seq)))
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
                                # print("after barrnap results: does it have some data or not:", os.stat(barrnap_out).st_size != 0)
                                if os.stat(barrnap_out).st_size != 0:
                                    sequences_dict = barrbap_dict_creation(barrnap_out)
                                
                                    # remove the output files created by barrnap
                                    # os.system("rm -f " + str(gb_record.id) + "*.fai")
                                    # os.system("rm -f " + barrnap_out)

                                    # provides keys only in the dictionary
                                    keys = list(sequences_dict.keys())
                                    # extracts sequences from the dictionary
                                    sequences = list(sequences_dict.values())
                                    # print("keys are after barrnap execution : ", keys, len(keys), keys != [] )
                                    # earlier below was keys != 0:
                                    if keys != []:
                                        chr_ref = gb_record.id
                                    # array dimensions [dim x dim]
                                    dim = len(keys)
                                    # print("dimensions : ", dim)
                                    # checks if the detected 16sRNA sequences are two or less
                                    if dim > 2:
                                        """ If 16sRNA sequences are more than two then search for
                                        more similar to each other"""
                                        # print("detected 16sRNA are more than 2")
                                        # determines consensus sequences
                                        # print("sequence selection in progress !!!!")
                                        consensus_seq = consensus_sequence(dim, sequences,
                                                                        sequences_dict,
                                                                        keys)
                                        # print("sequence selected !!!!")
                                    elif dim == 2 or dim == 1:
                                        # print("\ndetected 16sRNA are 2 or less !!!!\n")
                                        consensus_seq = sequences[0]
                                    # else:
                                        # print("continue")
                                        # continue
                                
                                    
                                else:
                                    consensus_seq = "No 16sRNA sequence found"
                                    # print("No 16sRNA sequence found")
                                
                                # remove the output files created by barrnap
                                # os.system("rm -f " + str(gb_record.id) + "*.fai")
                                # os.system("rm -f " + barrnap_out)
                                # time.sleep(3)
                                
                                # dna_seq is being encoded to save space i.e. BLOB entry in DB
                                dna_seq = zlib.compress(str(gb_record.seq).encode('utf-8'))

                                # insert organism id, dna and 16sRNA sequence.
                                # print("input to insert_db: ", con)
                                # print("input to insert_db: ", gb_record.id)
                                # print("input to insert_db: ", consensus_seq)
                                # print("input to insert_db: ", chr_ref)
                                insert_db(con, str(gb_record.id), dna_seq, consensus_seq, chr_ref)
                                # print("SqliteDB is successfuly updated !!!!\n")
                                # remove the redundant/temp file
                                # remove_file()

                            else:
                                dna_seq = "NA"
                                consensus_seq = "NA"
                                chr_ref = "NA"
                                print("ELSE STATEMENT CALLED !!!!!!!!!\n")
                                insert_db(con, str(gb_record.id), dna_seq, consensus_seq, chr_ref)
                            # print("rows_genome_dna:", gb_re?cord.id, dna_seq )
                            rows_genome_dna.append([gb_record.id, dna_seq])
                            # print(rows_genome_dna)
                            try:
                                # cur.executemany("INSERT INTO genbank_generell
                                # (acc, genes_total, cds_total, genome_size, circular_linear)
                                # VALUES (?,?,?,?,?)", rows_genbank)
                                cur.executemany("""INSERT INTO proteins (acc, gen_locus,
                                                gene_name, start_pos, end_pos, orientation,
                                                sequence) VALUES (?,?,?,?,?,?,?)""",
                                                rows_proteins)
                                # cur.executemany("INSERT INTO genome_dna (acc, dna)
                                # VALUES (?,?)", rows_genome_dna)
                                con.commit()
                            except con.IntegrityError:
                                sys.stderr.write("ERROR: ID already exists in PRIMARY KEY\
                                                column " + str(con.IntegrityError) + "\n")
                                continue
                        except:
                            print("len(str(gb_record.seq)) in exception:" , len(str(gb_record.seq)))
                            print("something went wrong with extraction of data")
            # print(user_acc, chr_ref)
        
        except:
            print("ERROR in extraction and insertion of data")
    os.system("rm -f *.fai barrnap.*")
    return user_acc, chr_ref


def get_data_from_ncbi(ftp_lnk, chr_ref, user_acc, con, flag, chr_flag):
    ts = calendar.timegm(time.gmtime())
    tmp_download_path = os.getcwd() + "/" + str(ts) + "_TMP"
    os.makedirs(tmp_download_path)
    # print("tmp_download_path : ", tmp_download_path, "\n")
    # con = lite.connect(args.sqlite)
    # with con:
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
            # print("uncompressed_gb:", uncompressed_gb)
            outF = open(uncompressed_gb, 'wb')
            outF.write(inF.read())
            inF.close()
            outF.close()
            # integrate data into user's db
            user_acc, gb_id = insert_genbank(uncompressed_gb, con, user_acc, chr_ref, chr_flag)
        else:
            # print(flag)
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
            # print("ncbi_file_path in gbk: ", ncbi_file_path)
            user_acc, gb_id = insert_genbank(download_path, con, user_acc, chr_ref, chr_flag)

        # clean tmp_download_path
        # print("tmp delete !!!!!!!!!!\n")
        # print("tmp_download_path", tmp_download_path)
        shutil.rmtree(tmp_download_path)

        return user_acc, gb_id

    except:
        # error_list.append(user_acc)
        sys.stderr.write("***FTP/HTTP ERROR: " + user_acc +
                        " CANNOT BE DOWNLOADED! ***")
        # print(error_list)
        # fayyaz undo comments below
        shutil.rmtree(tmp_download_path)


def get_data(con, user_id, user_acc, user_center,
             user_range, organism_col, ftp_lnk, organism_arr, chr_ref):
    results = ""
    # con = lite.connect(args.sqlite)
    # print("con in get_data", con)
    with con:
        cur = con.cursor()
        left_corner = int(user_center) - int(user_range)
        if left_corner < 0:
            left_corner = 0
        right_corner = int(user_center) + int(user_range)

        try:
            # Fayyaz acc instead of * to save memory
            user_acc = user_acc.strip().split(".")[0]
            chr_id = chr_ref.strip().split(".")[0]
            # print("user_acc, chr_id: ",  user_acc, chr_id)

            # Fayyaz, please check if this needed anymore !!!!!!!!!!!
            cur.execute("SELECT acc FROM genome_dna WHERE acc LIKE ? OR acc LIKE ?", (user_acc + '._%', chr_id + '._%'))
            rows = cur.fetchall()
            con.commit()
            # print("rows :::", rows)
            print()
            user_flag = False
            chr_flag = False
            gb_id = 'NA'
            if len(rows) != 0:
                # print("lenght is not zero !!!!!!!!!! \n")
                for line in range(len(rows)):
                    if user_acc in rows[line][0]:
                        user_flag = True
                        # print("in line: ", line)
                    elif chr_id in rows[line][0]:
                        chr_flag = True
                        # print("else", line)
            if not chr_flag or not user_flag:
                if user_acc in organism_col:
                    flag = True
                    # print("before get_data_From_ncbi: ", ftp_lnk, chr_ref, user_acc, con, flag, chr_flag, "\n")
                    user_acc, gb_id = get_data_from_ncbi(ftp_lnk, chr_ref, user_acc, con, flag, chr_flag)
                    rows = [1]
                    user_acc = user_acc.strip().split(".")[0]

                    # cur.execute("""SELECT * FROM proteins WHERE acc LIKE ?1 and
                                # (start_pos<=?3 and end_pos>=?4) or acc LIKE ?2 and 
                                # (start_pos<=?3 and end_pos>=?4)""",
                                # (user_acc + '._%',  gb_id, right_corner, left_corner))
                    cur.execute("""SELECT * FROM proteins WHERE acc LIKE ?1 and
                                (start_pos<=?3 and end_pos>=?4) or acc LIKE ?1 and 
                                (start_pos<=?3 and end_pos>=?4)""",
                                (user_acc + '._%',  gb_id, right_corner, left_corner))
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
                # Fayyaz added to deal with organism not in the lookup file
                else:
                    flag = False
                    chr_ref =  user_acc
                    organism_col = user_acc
                    ftp_lnk = 'NA'
                    chr_flag = "NA"
                    # print(organism_col,ftp_lnk, chr_ref, user_acc, flag)
                    user_acc, gb_id = get_data_from_ncbi(ftp_lnk, chr_ref, user_acc, con, flag, chr_flag)

                    # Fayyaz added to show results when organism downloaded from ncbi directly.
                    user_acc = user_acc.strip().split(".")[0]
                    cur.execute("""SELECT * FROM proteins WHERE acc LIKE ?1 and
                                (start_pos<=?3 and end_pos>=?4) or acc LIKE ?1 and 
                                (start_pos<=?3 and end_pos>=?4)""",
                                (user_acc + '._%',  gb_id, right_corner, left_corner))
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
            # results += str(user_id)
            print("Exception in get_data () :", user_id)
            # pass
        # con.close()
    # print("exception in get_data :", results)
    return results


def process_NCBI_lookup(in_file, directory, out_file):
    # print("Fayyaz process ncbi lookup exe starts !!!!!!!!!!\n")
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
                # print("HEre is what you writing: " + file_identifiers)
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
    # os.remove(in_file)

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
        # print("name:", name)
        # print(records)
        for (start, stop, strand) in positions[name]:
            short_seq = str(dna_seq)[start-1:stop]
            if strand == "+":
                short_seq_all[name] = short_seq
            elif strand == "-":
                bases = list(short_seq) 
                bases = reversed([complement.get(base,base) for base in bases])
                bases = ''.join(bases)
                short_seq_all[name] = bases
    return short_seq_all

def chrom_plas(input_param):
    # User starts a search - check if the ACC is in the database or not.
    # If not, download the genbank from ncbi
    # print("opens the lookup table !!!!!!!!\n")
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
        # print()
        for i, val in enumerate(line_chr):
            # ensures refseq IDs are downloaded from NCBI directly, 
            # not from the ftp link
            if "/" in val:
                lookup_id = val.strip().split("/")[1]
            else:
                lookup_id = val.strip().split()[0]
            if input_param in lookup_id:
                # # we can use index to ensure the first organism is choosen as the chr_ref
                # print("ind_input_param: ", i)
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


def update_db(con, id_container):
    # print("id_container, args_accession: ", id_container, args_accession)
    final_results = ""
    # con = lite.connect(args.sqlite)
    # print("id_container:", id_container)
    # print("ncbi_file_path:", ncbi_file_path)
    with con:
        cur =  con.cursor()
        for i in range(0, len(id_container)):
            # print("len of len(id_container): ", len(id_container))
            left_corner = int(id_container[i][2]) - int(id_container[i][3])
            if left_corner < 0:
                left_corner = 0
            right_corner = int(id_container[i][2]) + int(id_container[i][3])
            input_param = id_container[i][1].strip().split(".")[0]
            try:
                cur.execute("SELECT acc FROM genome_dna WHERE acc LIKE ?", (input_param + "._%",))
                chr_ref = cur.fetchall()
                # con.commit()
                # input_param = chr_ref[0][0]
                # con.commit()
                if chr_ref != []:
                    # chr_ref[0][0] is an input_organism
                    cur.execute("""SELECT * FROM proteins WHERE acc=?1 and
                                (start_pos<=?2 and end_pos>=?3)""",
                                (chr_ref[0][0], right_corner, left_corner))
                    rows = cur.fetchall()
                    con.commit()
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
                    organism_col, ftp_lnk, organism_arr, chr_ref = chrom_plas(input_param)
                    # print("organism_col, ftp_lnk, organism_arr, chr_ref:", organism_col, ftp_lnk, organism_arr, chr_ref)
                    # print("\nget_data parameters :",  con, id_container[i][0], id_container[i][1],
                                            # id_container[i][2], id_container[i][3],
                                            # organism_col, ftp_lnk, organism_arr, chr_ref)
                    final_results += get_data(con, id_container[i][0], id_container[i][1],
                                            id_container[i][2], id_container[i][3],
                                            organism_col, ftp_lnk, organism_arr, chr_ref)
            except:
                print("ERROR IN update data function !!!!!!!!!!!!\n")    
    return final_results

def find_srRNA_gene(args_accession, args_srRNA, id_container):
    # srRNA_flag = False
    # print("FAYYAZ !!!!!!!!!!!!!!!!!\n")
    print()
    db_path = args.sqlite
    # db_path = os.getcwd() + "/" + str(args.sqlite).strip(".")
    # print("db_path", db_path)
    # print("args.sqlite: ", args.sqlite)
    final_results = ''
    if args_srRNA:
        con = lite.connect(db_path)
        with con:
            cur = con.cursor()
            try:
                for line in args_srRNA:
                    input_param = line.strip().split(".")[0]
                    print("input_param ", input_param, "!!!!!!!!!!!!!\n")
                    cur.execute("SELECT [16sRNA_Ref] FROM genome_dna WHERE acc LIKE ?", (input_param + "._%",))
                    chr_ref = cur.fetchall()
                    con.commit()
                    if chr_ref == []:
                        print("input_param not found in DB ", input_param, " !!!!!!!1\n")
                        for data in id_container:
                            if input_param in data[1]:
                                input_list = [data]
                                # print("input_list: ", input_list, "!!!!!!!!!\n")
                                final_results = update_db(con, input_list)
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
                            cur.execute("SELECT * FROM [16sRNA] WHERE acc LIKE ?", (chr_ref[0][0],))
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
        with con:
            final_results = update_db(con, id_container)

            if args.output != "":
                # store data into file
                handle = open(args.output, "w")
                handle.write(final_results)
                handle.close()
            else:
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
    # print(args_sqlite)
    con = lite.connect(args.sqlite)
    with con:

        cur = con.cursor()
        dna_flag = False

        positions = defaultdict(list)
        short_seq_all = dict()
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 

        with open(args_extDNA, 'r') as fr:
            for line in fr:
                id, start, stop, strand = line.split()
                # print(id, start, stop, strand)
                positions[id].append((int(start), int(stop), str(strand)))
                input_param = line.strip().split(".") [0]
                cur.execute("SELECT dna FROM genome_dna WHERE acc LIKE ?", (input_param + "._%",))
                dseq = cur.fetchall()
                con.commit()
                # print(dseq)
                if dseq != []:
                    # print("there is a sequence in DB")
                    short_seq_all = seq_extraction(dseq, positions, complement)
                else:
                    final_results = ""
                    organism_col, ftp_lnk, organism_arr, chr_ref = chrom_plas(input_param)
                    # print("return of chrom_plas: ", organism_col, ftp_lnk, organism_arr, chr_ref)
                    final_results += get_data(con, line[0], input_param,
                                            int(args.position), int(args.range),
                                            organism_col, ftp_lnk, organism_arr, chr_ref)
                    cur.execute("SELECT dna FROM genome_dna WHERE acc LIKE ?", (input_param + "._%",))
                    dseq = cur.fetchall()
                    con.commit()
                    if dseq != []:
                        short_seq_all = seq_extraction(dseq, positions, complement)

    print()
    for data in short_seq_all:
        print(">" + data)
        print(short_seq_all[data])

    # close db connection
    # con.close()



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
                        organism and returns a refseq ID!",\
                        type=str, default="")
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
        insert_genbank(org_name, con, org_name.strip().split(".")[0], org_name.strip().split(".")[0], 'FALSE')
        # con.close()
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
        input_param = (args.genbank).strip().split(".")[0]
    if args.srRNA:
        input_param = args.srRNA
    if args.accession:
        input_param = (args.accession).strip().split(".")[0]
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
        refseq = find_refseq(ncbi_file_path, args.refseq)
        print("refseq ID of " , input_param , " is " , refseq , "!!!!!!!!!!\n")
        sys.exit()
        
    if args.accession or args.srRNA:
        # print("args.accession, args.srRNA", args.accession, args.srRNA)
        find_srRNA_gene(args.accession, args.srRNA, id_container )

    if args.extDNA:
        # print(args.extDNA)
        dna_extraction(args.extDNA, args.sqlite, ncbi_file_path)
        
