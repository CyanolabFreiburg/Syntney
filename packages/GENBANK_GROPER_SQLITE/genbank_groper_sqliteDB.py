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
# SQLite related
import sqlite3 as lite


# absolute path of this script
scriptPath = os.path.dirname(os.path.realpath(__file__))

# GLOBAL VARIABLES
genbank_summary_url = "ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/prokaryotes.txt"
ncbi_folder_name = scriptPath + "/myNCBI-LOOKUP"
ncbi_master_table = "NCBI-MASTER.table"
# SQLite related
con = None


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
            #cur.execute("CREATE TABLE genbank_generell(acc TEXT PRIMARY KEY, genes_total INTEGER, cds_total INTEGER, genome_size INTEGER, circular_linear TEXT)")
            # protein table
            cur.execute("CREATE TABLE proteins(id INTEGER PRIMARY KEY AUTOINCREMENT, acc TEXT, gen_locus TEXT, gene_name TEXT, start_pos INTEGER, end_pos INTEGER, orientation TEXT, sequence TEXT)")
            # gene name table
            #cur.execute("CREATE TABLE gene_names(id INTEGER PRIMARY KEY AUTOINCREMENT, gen_locus TEXT, gen_name TEXT)")
            # 16S table todo
            # genome DNA content todo
            cur.execute("CREATE TABLE genome_dna(acc TEXT PRIMARY KEY, dna TEXT)")
            # create indexes for fast lookup
            cur.execute("CREATE INDEX acc ON proteins (acc)")
            cur.execute("CREATE INDEX prot_start_pos ON proteins (start_pos)")
            cur.execute("CREATE INDEX prot_end_pos ON proteins (end_pos)")
            #cur.execute("CREATE INDEX gen_locus ON gene_names (gen_locus)")
            con.commit()
    except lite.Error:
        sys.stderr.write("SQLite Error " + str(lite.Error) + "\n")
        sys.exit(1)

    finally:
        if con:
            con.close()
    return 0


def insert_genbank(genbank, con):
    state = ""
    # 1 - check if genbank is a file or folder and download file if the data is not in the db

    if os.path.isfile(genbank):
        file_arr = glob.glob(genbank)
    if os.path.isdir(genbank):
        genbank = str(genbank) + "*"
        file_arr = glob.glob(genbank)

    cur = con.cursor()
    for tmp_path in file_arr:
        for gb_record in SeqIO.parse(open(tmp_path, "r"), "genbank"):
            # check if record exists. If, then take the next entry
            cur.execute("SELECT acc FROM proteins WHERE acc=?", (gb_record.id,))
            data = cur.fetchall()

            if len(data) != 0:
                state += str(gb_record.id) + " - exists!\n"
                continue

            rows_genbank = list()
            rows_proteins = list()
            rows_genome_dna = list()

            # todo not well implemented yet
            rows_genbank.append((gb_record.id, 1, 2, 3, "foo"))
            # get data from genbank and rearrange the structure
            for entry in gb_record.features:
                if entry.type == "CDS" and ("translation" in entry.qualifiers):
                    location = entry.location
                    # check if location is unfinished -> [<0:687](-) or 7251..>7537
                    if len(re.findall("(<|>)", str(location))) == 0:
                        # multiple hits for the same sequence
                        if str(location).startswith("join"):
                            locs = re.findall("{(.*?)}", str(location))[0]
                            locs_arr = locs.split(",")
                            # number of coordinates
                            for coords in locs_arr:
                                start, end, ori = extract_coordinates(coords)
                                seq = zlib.compress(entry.qualifiers["translation"][0].encode('utf-8'))

                                if "locus_tag" in entry.qualifiers:
                                    gene_loc = entry.qualifiers["locus_tag"][0]
                                else:
                                    gene_loc = "na"
                                if "gene" in entry.qualifiers:
                                    gene = entry.qualifiers["gene"][0]
                                else:
                                    gene = "na"

                                rows_proteins.append([gb_record.id, gene_loc, gene, start, end, ori, seq])
                        else:
                            start, end, ori = extract_coordinates(location)
                            seq = zlib.compress(entry.qualifiers["translation"][0].encode('utf-8'))

                            if "locus_tag" in entry.qualifiers:
                                gene_loc = entry.qualifiers["locus_tag"][0]
                            else:
                                gene_loc = "na"
                            if "gene" in entry.qualifiers:
                                gene = entry.qualifiers["gene"][0]
                            else:
                                gene = "na"

                            rows_proteins.append([gb_record.id, gene_loc, gene, start, end, ori, seq])
            if len(str(gb_record.seq)) > 0:
                dna_seq = zlib.compress(str(gb_record.seq).encode('utf-8'))
            else:
                dna_seq = "na"
            rows_genome_dna.append([gb_record.id, dna_seq])
            try:
                #cur.executemany("INSERT INTO genbank_generell (acc, genes_total, cds_total, genome_size, circular_linear) VALUES (?,?,?,?,?)", rows_genbank)
                cur.executemany("INSERT INTO proteins (acc, gen_locus, gene_name, start_pos, end_pos, orientation, sequence) VALUES (?,?,?,?,?,?,?)", rows_proteins)
                cur.executemany("INSERT INTO genome_dna (acc, dna) VALUES (?,?)", rows_genome_dna)
                con.commit()
            except con.IntegrityError:
                sys.stderr.write("ERROR: ID already exists in PRIMARY KEY column " + str(con.IntegrityError) + "\n")
                continue
    return 0


def get_data_from_ncbi(ncbi_master_table, user_acc, con):
    ts = calendar.timegm(time.gmtime())
    tmp_download_path = "./" + str(ts) + "_TMP"
    os.makedirs(tmp_download_path)

    # external wget function is used
    try:
        os.system("wget -c --retry-connrefused -nv --show-progress --continue --read-timeout=20 --tries=40 --wait=10 --timeout=15 " +
                  str(ncbi_master_table[user_acc]) + " -P " + str(tmp_download_path))

        file_path = str(tmp_download_path) + "/" + str(ncbi_master_table[user_acc].split("/")[-1])

        # uncompress downloaded *.gz file
        inF = gzip.open(file_path, 'rb')
        uncompressed_gb = str(file_path) + ".dec"
        outF = open(uncompressed_gb, 'wb')
        outF.write(inF.read())
        inF.close()
        outF.close()

        # integrate data into user's db
        insert_genbank(uncompressed_gb, con)

        # clean tmp_download_path
        shutil.rmtree(tmp_download_path)
    except:
        sys.stderr.write("*** FTP ERROR: " + str(ncbi_master_table[user_acc]) + " CANNOT BE DOWNLOADED! ***")
        shutil.rmtree(tmp_download_path)


def get_data(conn, user_id, user_acc, user_center, user_range, ncbi_master_table):
    results = ""
    cur = conn.cursor()
    left_corner = int(user_center) - int(user_range)
    if left_corner < 0:
        left_corner = 0
    right_corner = int(user_center) + int(user_range)

    try:
        cur.execute("SELECT * FROM genome_dna WHERE acc=?", (user_acc,))
        rows = cur.fetchall()
        if len(rows) == 0:
            # try to get data from online repo
            if user_acc in ncbi_master_table:
                get_data_from_ncbi(ncbi_master_table, user_acc, conn)
                rows = [1]

        if len(rows) == 0:
            # DESCRIPTION => row[0]=unique_number row[1]=acc row[2]=locus_name row[3]=gene_name row[4]=start_pos row[5]=end_pos row[6]=ori row[7]=protein_seq
            results += str(user_id) + "\t" + str(user_acc) + "\t" + "missing entry in LUT" + "\t" + "missing entry in LUT" + \
                       "\t" + "missing entry in LUT" + "\t" + "missing entry in LUT" + "\t" + "missing entry in LUT" + "\t" + "missing entry in LUT" + "\n"
        else:
            cur.execute("SELECT * FROM proteins WHERE acc=? and (start_pos<=? and end_pos>=?)", (user_acc, right_corner, left_corner,))
            rows = cur.fetchall()
            if len(rows) == 0:
                results += str(user_id) + "\t" + str(user_acc) + "\t" + "no annotation" + "\t" + "no annotation" + "\t" + "no annotation" + "\t" + "no annotation" + "\t" + \
                              "no annotation" + "\t" + "no annotation" + "\n"
            else:
                for row in rows:
                    # DESCRIPTION => row[0]=unique_number row[1]=acc row[2]=locus_name row[3]=gene_name row[4]=start_pos row[5]=end_pos row[6]=ori row[7]=protein_seq
                    results += str(user_id) + "\t" + str(row[1]) + "\t" + str(row[3]) + "\t" + str(row[2]) + "\t" + \
                               str(row[4]) + "\t" + str(row[5]) + "\t" + str(row[6]) + "\t" + str(zlib.decompress(row[7]).decode('utf-8') + "\n")
    except:
        # entry not available
        # todo check if exception can happen and document what's the reason for
        results += str(user_id)
    cur.close()
    return results


def process_NCBI_lookup(in_file, directory, out_file):
    master_lookup = dict()
    result = list()
    handle = open(in_file, "r", encoding='ascii', errors="replace")
    for line in handle:
        line = line.rstrip()
        tmp_arr = line.split("\t")
        if not line.startswith("#") and (len(tmp_arr) == 23):
            try:
                # Replicons[8] + FTP Path [20]
                # extract identifiers; only ACC
                tmp_identifiers_arr = tmp_arr[8].split("; ")
                file_identifiers = ""
                for identifier in tmp_identifiers_arr:
                    tmp_acc = re.findall(".*/(.*)", str(identifier))

                    if len(tmp_acc) > 0:
                        if file_identifiers == "":
                            file_identifiers = tmp_acc[0]
                        else:
                            file_identifiers += "," + str(tmp_acc[0])
                # get real ftp path
                if file_identifiers != "":
                    missing_part = tmp_arr[20].split("/")[-1]
                    file_identifiers += "\t" + str(tmp_arr[20]) + "/" + str(missing_part) + "_genomic.gbff.gz" + "\n"
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
    # clean original file
    os.remove(in_file)
    return master_lookup


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--genbank", help="genbank file or path to genbank files", type=str, default="")
    parser.add_argument("-s", "--sqlite", help="Path to SQLite DB", type=str, default="./mySQLiteDB.db")
    parser.add_argument("-a", "--accession", help="accession id e.g. CP021219.1", type=str, default="")
    parser.add_argument("-c", "--position", help="search position; not needed if parameter -a is a file", type=int, default=10)
    parser.add_argument("-r", "--range", help="considered range for extracting data around the search position; "
                                              "not needed if parameter -a is a file", type=int, default=3000)
    parser.add_argument("-i", "--id", help="Set unique identifier. Used in batch-processing!", type=str, default="na")
    parser.add_argument("-o", "--output", help="Specify file name for the results. Otherwise, std is used", type=str, default="")
    parser.add_argument("-u", "--update", help="auto updater checks every 30 days for an updated NCBI "
                                               "lookup table (on-30 , every 30 days) -> values: on-30, off, man", type=str, default="on-60")
    parser.add_argument("-x", "--glassgo", help="Parameter takes a finite string and returns a string with dna sequences - replaces blastdbcmd!", type=str, default="")


    args = parser.parse_args()

    # if sqlite-db-file not exists, then create new db
    if not os.path.exists(args.sqlite):
        init_sqlite_db(args.sqlite)

    # check if there are new data available to store in the database
    if os.path.exists(args.genbank):
        # start pickle of data
        con = lite.connect(args.sqlite)
        insert_genbank(args.genbank, con)

    # check if ncbi lookup-file needs to be updated or not!
    update_mode = args.update.split("-")
    if update_mode[0] == "on" and os.path.exists(ncbi_folder_name):
        # get date of current version
        ncbi_file_path = str(ncbi_folder_name) + "/" + str(ncbi_master_table)
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

    # check if ncbi lookup table exist, otherwise download the table from ncbi's ftp server
    if os.path.isdir(ncbi_folder_name):
        ncbi_master_table_path = str(ncbi_folder_name) + "/" + str(ncbi_master_table)
        if os.path.isfile(ncbi_master_table_path):
            pass
        else:
            # get data and reprocess the file
            file_path = wget.download(genbank_summary_url, out=ncbi_folder_name)
            process_NCBI_lookup(file_path, ncbi_folder_name, ncbi_master_table)
            print()
    else:
        os.mkdir(ncbi_folder_name)
        # get data and reprocess the file
        file_path = wget.download(genbank_summary_url, out=ncbi_folder_name)
        process_NCBI_lookup(file_path, ncbi_folder_name, ncbi_master_table)
        print()

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
    else:
        tmp_arr = list()
        tmp_arr.append(args.id)
        tmp_arr.append(args.accession)
        tmp_arr.append(int(args.position))
        tmp_arr.append(int(args.range))
        id_container.append(tmp_arr)

    # User starts a search - check if the ACC is in the database or not. If not, download the genbank from ncbi
    # build lookup table
    ncbi_lookup_dict = dict()
    ncbi_path = str(ncbi_folder_name) + "/" + str(ncbi_master_table)
    handle = open(ncbi_path, "r")
    for line in handle:
        line = line.rstrip()
        line_arr = line.split("\t")
        id_arr = line_arr[0].split(",")
        for tmp_id in id_arr:
            ncbi_lookup_dict[tmp_id] = line_arr[1]
    handle.close()

    # get final results
    final_results = ""
    con = lite.connect(args.sqlite)
    for i in range(0, len(id_container)):
        final_results += get_data(con, id_container[i][0], id_container[i][1], id_container[i][2], id_container[i][3], ncbi_lookup_dict)
    # close db connection
    con.close()

    # store or print final results
    if args.output != "":
        # store data into file
        handle = open(args.output, "w")
        handle.write(final_results)
        handle.close()
    else:
        print(final_results)