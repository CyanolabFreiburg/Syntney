import os
from Bio import SeqIO
import numpy as np
#from colour import Color
import argparse
import subprocess
from ete3 import *

# produces synteny file and cluster file with an R script. Therefore uses the input network fasta file and test fasta
# file. Synteny window is used for the extraction window and is per default set to 5000
# returns a dictionary of network_ids and test_ids.

# input:
# network_file:     path to file used for network construction
# test_file:        path to fasta file with questionable sequences
# wdir: working     directory where temporary files get stored

# output:
# network_ids, test_ids:    {sequence_id: [seq_desciption, seq_sequence]}
# clusterfile:              clusterfile that was produced by R script
# syntenyfile:              syntenyfile that was produced by R script
# r_script_path:            path to the synteny clustering R script
# synteny_window:           up and downstream number of bp of sequence that is searched for protein coding sequences
def run_r_script(network_file, test_file, wdir, r_script_path, sql_db_path, sql_script_path, synteny_window=str(5000)):
    print("pla")
    seqdict = dict()
    network_ids = dict()
    for seq_record in SeqIO.parse(network_file, "fasta"):
        seqdict.update({seq_record.id: [seq_record.description, seq_record.seq]})
        seq_id = seq_record.description
        seq_id = seq_id.split("-")
        if len(seq_id) == 1:
            pass
        else:
            seq_id = seq_id[0]
            seq_id = seq_id.split(":")
            if len(seq_id) > 1:
                if seq_id[1].startswith("c"):
                    seq_id[1] = seq_id[1][1:]
                seq_id = seq_id[0] + "_" + seq_id[1]
                network_ids.update({seq_id: [seq_record.description, seq_record.seq]})
    test_ids = dict()
    if test_file is not None:
        for seq_record in SeqIO.parse(test_file, "fasta"):

            seq_id = seq_record.description
            seq_id = seq_id.split("-")
            if len(seq_id) == 1:
                pass
            else:
                seq_id = seq_id[0]
                seq_id = seq_id.split(":")
                if len(seq_id) > 1:
                    if seq_id[1].startswith("c"):
                        seq_id[1] = seq_id[1][1:]
                    seq_id = seq_id[0] + "_" + seq_id[1]
                    if seq_id not in network_ids:
                        test_ids.update({seq_id: [seq_record.description, seq_record.seq]})
                        seqdict.update({seq_record.id: [seq_record.description, seq_record.seq]})
    else:
        test_ids = None
    f = open(wdir + "syntenyfile.fasta", "w")
    f.write(">pseudoquery\nGTA\n")
    for element in seqdict:
        desc, seq = seqdict[element]
        f.write(">" + desc + "\n" + str(seq) + "\n")
    f.close()
    
    subprocess.run(["R", "--slave", "-f " + r_script_path, "--args", "filename=" + wdir + "syntenyfile.fasta", "duplicates_allowed=TRUE", "synteny_window=" + synteny_window, "name=" + wdir + "syntenyfile", "coprarna_compatible=FALSE", "script_path=" + sql_script_path, "db_path=" + sql_db_path])
    r_script_cluster_table = wdir + "syntenyfile_cluster_table.txt"
    r_script_synteny_table = wdir + "syntenyfile_synteny_table.txt"
    os.system("rm " + wdir + "syntenyfile.fasta")

    return r_script_cluster_table, r_script_synteny_table, network_ids, test_ids


# produces a dictionary from the identifiers of sRNAs (ids). Identifiers must be like "Accessionnumber" + underscore +
# "starting position of hit" (e.g. CP001291.1_4248628). The synteny dictionary (synteny_dict) contains sRNA surrounding
# proteins. The returned dict has the following topology:
# input:
# ids                       dict of ids for which a synteny_dict is created
# r_script_synteny_table    synteny master table produced by run_r_script()
#
#
# output:
# {CP001291.1_4248628: [{Protein1:  "position to sRNA",
#                       Protein2: 2,                        (upstream proteins)
#                       Protein3: 1},
#                      {...}}                               (downstream proteins)
# infile specifies a synteny_table file from the synteny clustering R script.
def get_synteny_dict(ids, r_script_synteny_table):
    synteny_dict = dict()
    with open(r_script_synteny_table) as f:
        next(f)
        for line in f:
            line = line.rstrip()
            handle = line.split("\t")
            seq_id = handle[0]
            if seq_id not in synteny_dict:
                if seq_id in ids:   # checks whether the id from master table should be added to the synteny_dict
                    synteny_dict.update({seq_id: []})
                    proteins = handle[4].split(",") # all surrounding proteins
                    positions = handle[5].split(",")  # positions of all surrounding proteins positions[x] is position of proteins[x]
                    downstream_dict = {}    # dict of proteins downstream the sRNA ({protein: position, ... })
                    upstream_dict = {}      # proteins upstream the sRNA
                    switch = 0              # needed to switch at position 1 from upstream to downstream
                    for x in range(len(proteins)): # adds proteins to down and upstream dict
                        if int(positions[x]) < 10:
                            if switch == 0:
                                if positions[x] == "1":
                                    switch = 1
                            if switch == 2:
                                downstream_dict.update({proteins[x]: positions[x]})

                            if switch < 2:

                                upstream_dict.update({proteins[x]: positions[x]})
                                if switch == 1:
                                    switch = 2
                    synteny_dict[seq_id].append(upstream_dict)
                    synteny_dict[seq_id].append(downstream_dict)

    return synteny_dict


# Returns a cluster_dict where proteins from the synteny_table point on their clusters.
# {Protein1: Cluster_1,
#  Protein2: Cluster_2}
# The infile is the cluster_table produced by the synteny R script.
def get_clusters(r_script_cluster_table):
    f = open(r_script_cluster_table)
    cluster_dict = dict()
    for line in f:
        line = line.rstrip()
        if line.startswith("cluster"):
            handle = line.split("\t")
            name = handle[0]
            cluster = handle[1].split(",")
            for element in cluster:
                cluster_dict.update({element: name})

    return cluster_dict


# Uses a synteny_dict and a cluster_dict to APPEND clusters matching the proteins of entries in a synteny_dict.

# {CP001291.1_4248628: [{Protein1: 3, ...}             upstream Proteins
#                       {Protein4: 1, ...}             downstream Proteins
#                       [Cluster_1, Cluster_2, ...]    upstream cluster - positions equal position in array - APPENDED
#                       [CLuster_5, ...]               downstream cluster - positions equal position in array - APPENDED
def add_cluster_to_synteny_table(synteny_dict, cluster_dict, number_of_clusters):
    count = 0
    for entry in synteny_dict:
        up_proteins = synteny_dict[entry][0]# the upstream proteins of the considered sRNA
        down_proteins = synteny_dict[entry][1]  # the downstream proteins of the considered sRNA
        up_cluster = []         # creates a list of upstream clusters
        down_cluster = ["sRNA"] # adds sRNA to downstream clusters
        for protein in up_proteins:
            try:
                cluster = cluster_dict[protein]
                up_cluster.append(cluster)
                count += 1
            except KeyError:
                print("Cluster not found") # proteins without annotated aminoacid sequence do not have clusters

        for protein in down_proteins:
            try:
                cluster = cluster_dict[protein]
                down_cluster.append(cluster)
                count += 1
            except KeyError:                # proteins without annotated aminoacid sequence do not have clusters
                print("Cluster not found")

        up_cluster.append("sRNA")
        down_cluster = down_cluster[0:number_of_clusters]   # the number of clusters used in this approach
        up_cluster = list(reversed(up_cluster))[0:number_of_clusters]   # the number of clusters used in this approach
        synteny_dict[entry].append(up_cluster)
        synteny_dict[entry].append(down_cluster)

    return synteny_dict


# produces a network from a synteny_dict.
# output:
# network:      {cluster: {connected_cluster1: [number of appearance of this node, [list of Accessions with this node]],
#                   connected_cluster2: [...]}
def build_network(synteny_dict):
    network = dict()
    for entry in synteny_dict:
        upcluster = synteny_dict[entry][2]  # upstream cluster of a sRNA in the synteny_dict
        prev_cluster = 0    # previous cluster is 0 for first iteration
        for cluster in upcluster: # assigns node connections to the network
            if prev_cluster != 0:
                if prev_cluster == cluster:  # prevents loops in the network
                    pass
                else: # assigns the connection between cluster and previous cluster to the network
                    if cluster not in network:
                        network.update({cluster: [{prev_cluster: [1, [entry]]}, [entry]]})

                    else:
                        network[cluster][1].append(entry)
                        if prev_cluster not in network[cluster][0]:
                            network[cluster][0].update({prev_cluster: [1, [entry]]})
                        else:
                            network[cluster][0][prev_cluster][0] += 1
                            network[cluster][0][prev_cluster][1].append(entry)
            prev_cluster = cluster  # previous cluster is the cluster of the earlier iteration

        prev_cluster = 0
        downcluster = synteny_dict[entry][3]  # downstream cluster of a sRNA in the synteny_dict
        for cluster in downcluster:
            if prev_cluster != 0:
                if prev_cluster == cluster:
                    pass
                else:
                    if cluster not in network:
                        network.update({cluster: [{prev_cluster: [1, [entry]]}, [entry]]})

                    else:
                        network[cluster][1].append(entry)

                        if prev_cluster not in network[cluster][0]:
                            network[cluster][0].update({prev_cluster: [1, [entry]]})
                        else:
                            network[cluster][0][prev_cluster][0] += 1
                            network[cluster][0][prev_cluster][1].append(entry)
            prev_cluster = cluster

    return network


# builds and returns a ete3 tree from the sRNA sequences from a "fasta" infile (should be the trustable GLASSgo file).
# as the tree is built with numbers instead of the identifier ("accession id"_"starting nucleotide"), also a tree_iddict
# is returned where the id points on the corresponding number.
def tree_construction(wdir, network_file):
    count = 0
    tree_iddict = dict()
    forbidden = set()

    # produces a FASTA with numbers instead of original headers and a tree_iddict that is used to get the number
    # from an identifier
    f = open(wdir + "treefasta.tmp", "w")
    for seq_record in SeqIO.parse(network_file, "fasta"):
        if seq_record.description not in forbidden:
            seq_id = seq_record.id.split(":")
            if len(seq_id) > 1:
                pos = seq_id[1].split("-")[0]
                if pos.startswith("c"):
                    pos = pos[1::]
                seq_id = seq_id[0] + "_" + pos
                f.write(">" + str(count) + "\n" + str(seq_record.seq) + "\n")
                tree_iddict.update({seq_id: str(count)})
                count += 1
                forbidden.add(seq_record.description)
    f.close()
    # produces a distance matrix from the numbered FASTA via clustalo
    os.system("clustalo --infile=" + wdir + "treefasta.tmp -o align.out --distmat-out=" + wdir + "distmat.txt --full --force -v")
    os.system("rm align.out")
    # uses quicktree to built a tree from the distance matrix and removes the distance matrix
    os.system("quicktree -in m " + wdir + "distmat.txt >" + wdir + "tree.tmp")
    os.system("rm " + wdir + "distmat.txt")
    # produces a ete3 object from the tree and removes the treefile and the tree FASTA
    f = open(wdir + "tree.tmp")
    tree = str()
    for line in f:
        line = line.rstrip()
        tree = tree + line
    f.close()
    tree = Tree(tree, format=1)
    os.system("rm " + wdir + "tree.tmp")
    os.system("rm " + wdir + "treefasta.tmp")
    return tree_iddict, tree


# returns the whole length of an ete3 tree
def whole_tree_length(tree):
    leng = 0
    for node in tree.iter_descendants():
        leng += node.dist
    return leng


# calculates the sum of branches of a list with accessions. As the tree is built from identifiers also an tree iddict needs to be
# passed to convert the accession numbers into corresponding numbers. Returns the sum of branches containing all the edges to
# the lowest common ancestor of the passed accessionlist as well as the edge poitning to the parent node of the lca.
# identifier = "Accession_number1"_"startingnucleotide"
# input:
# tree          =   ete3.Tree
# accessionlist =   [identifier, ...]
# tree_iddict   =   {identifier: number that was used for the identifier in tree_construction}
# output:
# sob           =   float(sum of branches)
def sum_of_branches(tree, accessions_list, tree_iddict):
    acc_ids = []
    for entry in accessions_list: # writes identifier to numbers that were used in treeconstruction
        acc_ids.append(tree_iddict[entry])
    accessions = tuple(acc_ids)

    if len(accessions) > 1:
        n1 = tree.get_common_ancestor(accessions) # nl is the n(ode) of the l(ca)
        sob = n1.dist # adds distance to parent of lca

        lookup = set()

        for element in acc_ids: # sums up the branches to the leafs of the passes identifiers
            leaf = tree.get_leaves_by_name(element)[0]

            while leaf != n1 and leaf not in lookup:
                sob = sob + leaf.dist
                lookup.add(leaf)
                leaf = leaf.up
    else: # if only one identifier is passed, only the branch to its parent node is returned as sob

        node = tree.search_nodes(name=acc_ids[0])[0]
        parent = node.up
        dist = tree.get_distance(node, parent)
        sob = dist
    return sob


# needs a network as input and calculates the sum of outgoing connection weights for each node.
# This number is then added to each node in the network:
# input:
# {cluster: {connected_cluster1: [normalized weight of this connection, [list of Accessions with this node]],
#            connected_cluster2: [...]}
# output:
# {Cluster1: [outgoing connections, {connected_cluster1: ...}]}
# should be used after normalization of conections
def add_outgoing_connection_weights(network):
    for cluster in network:
        outgoing = 0
        for connected_cluster in network[cluster][0]:
            outgoing += network[cluster][0][connected_cluster][0]
        network.update({cluster: [outgoing] + network[cluster]})
    return network


# Normalizes the number of connections in a network with a sum of branches approach on a tree built from the
# sRNA sequences. infile is passed to the tree construction function and is therefore the fasta file of network sRNAs.
# The returned network has normalized outgoing connections and normalized number of connections.
# input:
# wdir: place where temporary files are stored
# network_file: FASTA file that was used for network construction
# network:  {cluster: {connected_cluster1: [number of appearance of this connection, [list of Accessions with this connection]],
#            connected_cluster2: [...]}, [list of Accessions with "cluster"]}
# output:
# network   {cluster: [normalized sum of outgoing connections, {connected_cluster1: [normalized weight of this connection,
#               [list of accessions with this connection], sum of branches weight], ...}[list of Accessions with "cluster"]}
def normalize_connections(wdir, network_file, network):
    tree_iddict, tree = tree_construction(wdir, network_file)
    treelength = whole_tree_length(tree) # whole tree length
    values = []
    zeroweights = []  # stores connections with a weight of 0
    for cluster in network:
        for connectedcluster in network[cluster][0]:
            accessions = network[cluster][0][connectedcluster][1]
            sob = sum_of_branches(tree, accessions, tree_iddict)
            if treelength == 0:
                value = 1
            else:
                value = sob/treelength

            if value != 0:
                values.append(value)
                network[cluster][0][connectedcluster][0] = value
            if value == 0:
                zeroweights.append([cluster, connectedcluster])
    if len(values) > 0:
        minimum = min(values)
    else:
        minimum = 1
    for entry in zeroweights:

        network[entry[0]][0][entry[1]][0] = minimum
    add_outgoing_connection_weights(network) # sums up the weights of outgoing connections and assigns them to the network
    for cluster in network: # splits up connection weights to an percentage value of their importance
        for connected_cluster in network[cluster][1]:
            network[cluster][1][connected_cluster].append(network[cluster][1][connected_cluster][0]) # appends the sob
            # weight that can be used for SV calculation instead of the weights used for PageRank
            network[cluster][1][connected_cluster][0] = \
                network[cluster][1][connected_cluster][0] / network[cluster][0]
    return network, tree, tree_iddict


# adds a teleport probability for PageRank usage to the network:
# input:
# tree:         ete3 tree from FASTA used for Network construction
# tree_iddict:  iddict for identifier number from headers of the Network FASTA file
# network:  {cluster: {connected_cluster1: [number of appearance of this connection, [list of Accessions with this connection]],
#            connected_cluster2: [...]}, [list of Accessions with "cluster"]]}
# output:
# network   {cluster: [normalized sum of outgoing connections, {connected_cluster1: [normalized weight of this connection,
#               [list of accessions with this connection]], ...},[list of Accessions with "cluster"], teleport prob.]}
def normalize_nodes(tree, tree_iddict, network):
    treelength = whole_tree_length(tree) # whole tree length
    values = []
    zeroweights = []  # stores connections with a weight of 0
    sum_of_clustervalues = 0
    for cluster in network:
        accessions = network[cluster][2]
        sob = sum_of_branches(tree, accessions, tree_iddict)
        if treelength == 0:
            value = 1
        else:
            value = sob / treelength
        if value != 0:
            values.append(value)
            sum_of_clustervalues += value
            network[cluster].append(value)
        if value == 0:
            zeroweights.append(cluster)
    if len(values) > 0:
        minimum = min(values)
    else:
        minimum = 1
    for entry in zeroweights:
        sum_of_clustervalues += minimum
        network[entry].append(minimum)
    for cluster in network:
        network[cluster][-1] /= sum_of_clustervalues
    return network


# edited from https://gist.github.com/joninvski/701720
# Step 1: For each node prepare the destination and predecessor
def initialize(graph, source):
    d = {}  # Stands for destination
    p = {}  # Stands for predecessor
    for node in graph:
        d[node] = 0  # We start admiting that the rest of nodes are not reachable and therefore have a value of zero
        p[node] = None
    d[source] = 1  # source has a distance of 1
    return d, p


# edited from https://gist.github.com/joninvski/701720
def relax(node, neighbour, graph, d, p, i):
    # If the distance between the node and the neighbour is higher than the one I have now
    i = i+1
    if d[neighbour] < (d[node] * graph[node][neighbour]) / i:
        # Record this higher distance
        d[neighbour] = (d[node] * graph[node][neighbour]) / i
        p[neighbour] = node


# edited from https://gist.github.com/joninvski/701720
# edited bellman ford alorithm multiplying edges with a weight between 0 and 1. Therefore the best path has a length
# of 1 and worse paths have a low weight.
# input:
# lite_network      network without outgoing edge weights
def bellman_ford(lite_network, source):
    d, p = initialize(lite_network, source)
    for i in range(len(lite_network)-1):  # Run this until is converges
        for u in lite_network:
            for v in lite_network[u]:  # For each neighbour of u
                relax(u, v, lite_network, d, p, i)  # Lets relax it

    return d, p


# uses the network to create a dictionary with the best paths:
# input:
# network           {cluster: [number outgoing connections, {connected_cluster1: [normalized number of connections,
#                       [list of Accessions with this node]], connected_cluster2: [...]], ...}
# sob_weights:      will use the sob weights for later SV calculation instead of the normalized values that were used
#                   for PageRank calculation if its set to True
# output:
# best_paths        {cluster: {connected_cluster: [distance between 0 and 1, first cluster on way to connected_cluster,
#                       number of steps]}}
def get_distances(network, sob_weights=False):
    # lite network is a data structure of a network without the sum of outgoing weights
    lite_network = dict()
    for cluster in network:
        if cluster not in lite_network:
            lite_network.update({cluster: dict()})
        for prev_cluster in network[cluster][1]:
            if sob_weights is True:
                value = network[cluster][1][prev_cluster][-1]
            else:
                value = network[cluster][1][prev_cluster][0]
            lite_network[cluster].update({prev_cluster: value})
            if prev_cluster not in lite_network:
                lite_network.update({prev_cluster: dict()})

    best_paths = dict()
    for entry in lite_network:
        distance, predecessor = bellman_ford(lite_network, source=entry)
        for dist in distance:
            if distance[dist] != 0:
                if entry != dist:
                    pred = predecessor[dist]
                    prevpred = dist
                    step = 1
                    while pred != entry:
                        prevpred = pred
                        pred = predecessor[pred]
                        step += 1
                    try:
                        # prevpred is the first cluster on the way to cluster(dist)
                        best_paths[entry].update({dist: [(distance[dist]), prevpred, step]})
                    except KeyError:
                        best_paths.update({entry: {dist: [(distance[dist]), prevpred, step]}})
    return best_paths


# uses a more complex approach to calculate the PageRanl of each cluster in a connectiondict (Network)
# the approach minimizes the used memory by not calculating the whole matrix. Therefore this approach is
# also able to handle big Networks.
# Changes the number of outgoing connection weights in the network to the pagerank_value.
# for a detailed description of the function look up in my master thesis chapter data structure and pagerank in the
# methods
#/media/cyano_share/documents/Bachelor- & Masterarbeiten/Master_Thesis_Dominik_Rabsch.pdf
# input:
# network:      {cluster: [number_outgoing_connections, {connected_cluster1: [normalized number of connections,
#                   [list of Accessions with this node]], connected_cluster2: [...]], ...}
# output:
# network:      {cluster: [pagerank_value, {connected_cluster1: [normalized number of connections,
#                   [list of Accessions with this node]], connected_cluster2: [...]], ...}
def pagerank(network, eps=1.0e-14, teleport=False):
    header = ["sRNA"]
    header = header + list(network)
    n = len(header)
    iddict = {}
    reverse_iddict = {}
    count = 0
    lines = []
    pagerank_vector = []
    teleport_vector = []
    for element in header:
        iddict.update({element: count})
        reverse_iddict.update({count: element})
        lines.append([])
        if teleport == True:
            teleport_vector.append(0)
        pagerank_vector.append(1/n)
        count += 1

    pagerank_vector = np.array(pagerank_vector, dtype="float64")
    i_table = []
    weights = []
    count = 0
    for cluster in network:
        if teleport == True:
            teleport_vector[iddict[cluster]] = network[cluster][-1]
        for connected_cluster in network[cluster][1]:
            i = iddict[cluster]
            i_table.append(i)
            value = network[cluster][1][connected_cluster][0]
            weights.append(value)

            lines[iddict[connected_cluster]].append(count)

            count += 1

    i_table = np.array(i_table)

    check = False
    count = 0
    while check is not True:
        check = True
        old_pagerank_vector = np.copy(pagerank_vector)
        for x in range(len(lines)):
            value = 0
            for index in lines[x]:
                i = i_table[index]
                weight = weights[index]
                value = value + weight * old_pagerank_vector[i]
            if teleport is True: #  teleport based on cluster occurences
                value += teleport_vector[x] * old_pagerank_vector[0]
            else:  # random teleport because sRNA column sums up to 0 (spider trap)
                value += (1/n) * old_pagerank_vector[0]
            diff = np.absolute(old_pagerank_vector[x] - value)
            pagerank_vector[x] = value
            if diff > eps:
                check = False
        count += 1
    pagerankdict = dict()

    for x in range(len(pagerank_vector)):

        pagerankdict.update({reverse_iddict[x]: pagerank_vector[x]})
    for entry in network:
        network[entry][0] = pagerankdict[entry]
    return network


# uses a synteny_dict dictionary as well as a  best_paths dict and a normalized Network that was used for pagerank
# calculation to calculate the synteny value of each sequence in the sequencedict. Afterwards the synteny value is
# stored in the sequencedict at position sequencedict[sequence_id][4]
# input:
# synteny_dict: {seq_id: [{upstream_Proteins},
#                      {downstream_ proteins},
#                      [upstream_Cluster],
#                      [downstream_Cluster]]}
# best_paths: cluster: {connected_cluster: [best_path , first cluster on way to connected_cluster, number of steps]}}
# network:
# {cluster: [pagerank_value, {connected_cluster1: [normalized number of connections, [list of Accessions with this node]],
#             connected_cluster2: [...]], ...}
# output:
# synteny_dict: {seq_id: [{upstream_Proteins},
#                      {downstream_ proteins},
#                      [upstream_Cluster],
#                      [downstream_Cluster],
#                       synteny_value]} "appends the synteny value here"
def calculate_synteny_value(synteny_dict, best_paths, network):

    for entry in synteny_dict:
        uppath = synteny_dict[entry][2] # upstream cluster of a considered entry in the synteny dict
        count = 0
        prevlist = ["sRNA"] # adds the sRNA to the already visited list and makes it possible to start from the sRNA
        synvalue = 0        # starting SV
        for z in range(len(uppath)):
            cluster = uppath[z]
            if count == 0: # does not calculate a value for the first cluster as it is the sRNA
                pass
            else:
                if cluster in prevlist:  # if the considered cluster is the same like a already visited cluster
                    synvalue += network[cluster][0] # add the pageRank of the cluster to the SV
                else:
                    if cluster in network:  # checks if the considered cluster is in the network
                        prevlist.append(cluster) # appends the cluster to list of already visited clusters
                        tmp = []
                        for cl in prevlist: # for cluster in already visited clusters
                            if cl in best_paths[cluster]: # checks if cluster is reachable from the already visited one
                                p = best_paths[cluster][cl][0] / best_paths[cluster][cl][2]
                                tmp.append(p) # appends the value of stepping to the cluster to a temporary list
                            elif cluster in best_paths[cl]: # checks if the cluster is reachable by stepping backwards
                                tar = 0 # sets up a value for the best path
                                div = 1 # is the number of edges that is used by stepping backwards
                                for clus in network[cl][1]: # for all clus cluster that are connected to already visited cluster
                                    if clus != "sRNA":
                                        if cluster in best_paths[clus]: # if target cluster is reachable from that clus by stepping backwards
                                            x = best_paths[cl][clus][0] / (best_paths[cl][clus][2] + 1) # x is the path value that is needed to reach this clus
                                            if x > tar: # tar stores the best of these paths to the target cluster
                                                tar = x
                                                div = (best_paths[cl][clus][2] + 1) # number of used edges
                                if z+1 != len(uppath): # need to add the best oputgoing edge of the target cluster
                                    if uppath[z + 1] in best_paths[cluster]: # if it is possible to reach that
                                        add = best_paths[cluster][best_paths[cluster][uppath[z+1]][1]][0] / div # add this path to the current path
                                    else:
                                        add = best_paths[cluster][best_paths[cluster]["sRNA"][1]][0] / div # add the edge of the cluster on the way to the sRNA
                                else:
                                    add = best_paths[cluster][best_paths[cluster]["sRNA"][1]][0] / div # if there is no next cluster in the considered synteny _ try to step towards the sRNA
                                if tar == 0:
                                    p = add
                                else:
                                    p = add * tar
                                tmp.append(p)
                        synvalue += max(tmp) * network[cluster][0]
                    else:
                        pass
            count += 1

        downpath = synteny_dict[entry][3]
        count = 0
        prevlist = ["sRNA"]
        for z in range(len(downpath)):
            cluster = downpath[z]
            if count == 0:
                pass
            else:
                if cluster in prevlist:
                    synvalue += network[cluster][0]
                else:
                    if cluster in network:
                        prevlist.append(cluster)
                        tmp = []
                        for cl in prevlist:
                            if cl in best_paths[cluster]:
                                p = best_paths[cluster][cl][0] / best_paths[cluster][cl][2]
                                tmp.append(p)
                            elif cluster in best_paths[cl]:
                                tar = 0
                                div = 1
                                for clus in network[cl][1]:
                                    if clus != "sRNA":
                                        if cluster in best_paths[clus]:
                                            x = best_paths[cl][clus][0] / (best_paths[cl][clus][2] + 1)
                                            if x > tar:
                                                tar = x
                                                div = (best_paths[cl][clus][2] + 1)
                                if z+1 != len(downpath):
                                    if downpath[z + 1] in best_paths[cluster]:
                                        add = best_paths[cluster][best_paths[cluster][downpath[z+1]][1]][0] / div
                                    else:
                                        add = best_paths[cluster][best_paths[cluster]["sRNA"][1]][0] / div
                                else:
                                    add = best_paths[cluster][best_paths[cluster]["sRNA"][1]][0] / div
                                if tar == 0:
                                    p = add
                                else:
                                    p = add * tar
                                tmp.append(p)
                        synvalue += max(tmp) * network[cluster][0]
                    else:
                        count -= 1
            count += 1
        synteny_dict[entry].append(synvalue)
    return synteny_dict


# must be used after the pagerank function. uses the connectiondict to create a svg outfile of the connectiondict
# Network with graphviz.
# input:
# {cluster: [pagerank, {connected_cluster1: [normalized number of connections, [list of Accessions with this node]],
#            connected_cluster2: [...]], ...}
def visualize_network(connectiondict, outfile):
    node_weights = []
    weightdict = dict()
    for cluster in connectiondict:
        weightdict.update({cluster[8:]: connectiondict[cluster][0]})
        node_weights.append(connectiondict[cluster][0])

    nodew = list(set(node_weights))
    red = Color("red")
    colors = list(red.range_to(Color("yellow"), len(nodew)))
    for x in range(len(colors)):
        colors[x] = str(Color(colors[x]).hex_l)

    nodew.sort(reverse=True)
    colordict = dict()
    for x in range(len(nodew)):
        colordict.update({nodew[x]: colors[x]})

    graph = Digraph()
    graph.node("sRNA", style="filled", fillcolor="green")
    for cluster in connectiondict:
        graph.node(cluster[8:], style="filled", fillcolor=colordict[connectiondict[cluster][0]])
    for cluster in connectiondict:
        for connected_cluster in connectiondict[cluster][1]:
            w = 3 * connectiondict[cluster][1][connected_cluster][0]
            if connected_cluster != "sRNA":
                connected_cluster = connected_cluster[8:]
            graph.edge(cluster[8:], connected_cluster, penwidth=str(w))
    try:
        graph.render(outfile, format="svg")
    except subprocess.CalledProcessError:
        pass


# must be used after the pagerank function. uses the network to create a cytoscape compatible comma seperated
# outfile.
# Network with graphviz.
# input:
# {cluster: [pagerank, {connected_cluster1: [normalized number of connections, [list of Accessions with this node]],
#            connected_cluster2: [...]], ...}, [list of accessions with "cluster"], teleport prob. to cluster], ...}
def visualize_cytoscape_network(network, outfile):
    f = open(outfile, "w")
    f.write("cluster,connected cluster,PageRank, connection weight\n")
    for cluster in network:
        pagerank = network[cluster][0]
        for connected_cluster in network[cluster][1]:
            weight = network[cluster][1][connected_cluster][0]
            f.write(cluster + "," + connected_cluster + "," + str(pagerank) + "," + str(weight) + "\n")


# produces an output file containing the sequence identifiers with their up and downstream cluster numbers.
# this file is used to observe the clusters of each protein in the corresponding network
# input:
# synteny_table: {seq_id: [{upstream_Proteins},
#                       {downstream_ proteins},
#                       [upstream_Cluster],
#                       [downstream_Cluster]]}
def output_cluster_synteny_file(syteny_table, outfile):
    f = open(outfile, "w")
    f.write("identiefier\tupstream cluster\tdownstream cluster\n")
    for entry in syteny_table:
        upstream = syteny_table[entry][2]
        downstream = syteny_table[entry][3]
        f.write(entry + "\t")
        for cluster in upstream:
            f.write(cluster + ",")
        f.write("\t")
        for cluster in downstream:
            f.write(cluster + ",")
        f.write("\n")
    f.close()


# normalizes the synteny value of the sequences used for network contruction to the max value of these.
# Also normalizes the synteny value of the tested sequences to the max value of the sequences used for network
# construiction if the test FASTA input file was added.
# input:
# *synteny_table: {seq_id: [{upstream_Proteins},
#                      {downstream_ proteins},
#                      [upstream_Cluster],
#                      [downstream_Cluster],
#                       synteny_value]}
# output:
# *synteny_table: {seq_id: [{upstream_Proteins},
#                      {downstream_ proteins},
#                      [upstream_Cluster],
#                      [downstream_Cluster],
#                      normalized synteny_value]}
def normalize_synteny_value(network_synteny_table, test_synteny_table):
    network_values = []
    for entry in network_synteny_table:
        network_values.append(network_synteny_table[entry][4])
    network_max = max(network_values)

    for entry in network_synteny_table:
        network_synteny_table[entry][4] /= network_max
    if test_synteny_table is not None: # tests if test sequences are added at the command line call
        for entry in test_synteny_table:
            test_synteny_table[entry][4] /= network_max


# uses a sequences dictionary where the synteny value is already calculated and a matching iddict to create an outfile.
# iddcit: {sequence_id: [seq_description, seq_sequence]}
# synteny_table: {seq_id: [{upstream_Proteins},
#                      {downstream_ proteins},
#                      [upstream_Cluster],
#                      [downstream_Cluster],
#                      normalized synteny_value]}
# outfile: path to outfile
def write_outfile_from_synteny_table(synteny_table, iddict, outfile):
    f = open(outfile, "w")
    for entry in synteny_table:
        desc, seq = iddict[entry]
        synvalue = synteny_table[entry][4]
        f.write(">" + desc + " synteny:" + str(synvalue) + "\n" + str(seq) + "\n")
    f.close()



def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--network_file", help="fasta file containing sequences used for network construction",
                        type=str)
    parser.add_argument("-t", "--test_file", help="optional fasta file containing sequences that are checked for network match",
                        type=str, default=None)
    parser.add_argument("-c", "--cluster_script", help="path to synteny clustering R script",
                        type=str,
                        default="/home/steffen/projects/Dominik_Master_Thesis/Synteny_Network/packages/Rscript/Synteny_Cluster_Script_sqlite.r")
    parser.add_argument("-p", "--w_dir", help="working directory where temporary files are stored default is the "
                                              "current directory", type=str, default="")
    parser.add_argument("-n", "--network", help="if set to svg, 'outfile'_Network.svg is produced as an output."
                                                "If set to cys, cytoscape compatible comma seperated Network 'outfile'_Network.txt is "
                                                "produced", type=str, default="False")
    parser.add_argument("-o", "--outfiles", help="path and name of outfiles. Will produce 'outfile'_network.fasta and "
                                                 "'outfile'_questionable.fasta ", type=str, default="")
    parser.add_argument("-w", "--synteny_window", help="synteny window used for extraction", type=int, default=5000)
    parser.add_argument("--protein_number", help="number of proteins up and downstream that should be used. default is 4",
                        type=int, default=4)
    parser.add_argument("--node_normalization",
                        help="If True uses a teleport at the sRNA based on a normalized number of cluster occurrences. "
                             "Default is False",type=bool, default=False)
    parser.add_argument("--use_sob_weights", help="If True uses sum of branch weights for Synteny Value calculation. "
                                                  "Default is False", type=bool, default=False)
    parser.add_argument("-d", "--sqlite_db", help="Path to SQLite DB", type=str, default="../Syntney_DB/mySQLiteDB_new.db")
    parser.add_argument("-s", "--sqlite_script", help="", type=str, default="./packages/GENBANK_GROPER_SQLITE/genbank_groper_sqliteDB.py")
    args = parser.parse_args()

    wdir = args.w_dir
    r_script_cluster_table, r_script_synteny_table, network_ids, test_ids = run_r_script(args.network_file, args.test_file,
                                                                                  wdir, args.cluster_script, args.sqlite_db, args.sqlite_script,
                                                                                  synteny_window=
                                                                                  str(args.synteny_window))
    number_of_clusters = args.protein_number + 1  # needs to be done as sRNA is also considered as a cluster
    network_synteny_table = get_synteny_dict(network_ids, r_script_synteny_table)
    cluster_dict = get_clusters(r_script_cluster_table)
    network_synteny_table = add_cluster_to_synteny_table(network_synteny_table, cluster_dict, number_of_clusters)
    network = build_network(network_synteny_table)
    network, tree, tree_iddict = normalize_connections(wdir, args.network_file, network)
    if args.node_normalization is True:
        normalize_nodes(tree, tree_iddict, network)
    best_paths = get_distances(network, sob_weights=args.use_sob_weights)
    network = pagerank(network, teleport=args.node_normalization)
    network_synteny_table = calculate_synteny_value(network_synteny_table, best_paths, network)
    if test_ids is not None:
        test_synteny_table = get_synteny_dict(test_ids, r_script_synteny_table)
        test_synteny_table = add_cluster_to_synteny_table(test_synteny_table, cluster_dict, number_of_clusters)
        test_synteny_table = calculate_synteny_value(test_synteny_table, best_paths, network)
    else:
        test_synteny_table = None
    normalize_synteny_value(network_synteny_table, test_synteny_table)
    write_outfile_from_synteny_table(network_synteny_table, network_ids, args.outfiles + "_network.fasta")
    if test_synteny_table is not None:
        write_outfile_from_synteny_table(test_synteny_table, test_ids, args.outfiles + "_questionable.fasta")
    if args.network == "svg":
        visualize_network(network, outfile=args.outfiles + "_Network.svg")
        output_cluster_synteny_file(test_synteny_table, outfile=args.outfiles + "cluster.txt")
    elif args.network == "cys":
        visualize_cytoscape_network(network, outfile=args.outfiles + "_Network.txt")
        if test_synteny_table is not None:
            output_cluster_synteny_file(test_synteny_table, outfile=args.outfiles + "test_cluster.txt")
        output_cluster_synteny_file(network_synteny_table, outfile=args.outfiles + "network_cluster.txt")
    else:
        pass


if __name__ == "__main__":
    main()
