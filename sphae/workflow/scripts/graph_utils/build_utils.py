import re

from Bio import SeqIO
from igraph import *
from collections import defaultdict

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2019, GraphBin Project"
__credits__ = "Vijini Mallawaarachchi"
__license__ = "BSD-3"
__type__ = "Support Script"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"

# -------------------------------------------------------------------
# Bidirectional map
# -------------------------------------------------------------------


class BidirectionalError(Exception):
    """Must set a unique value in a BijectiveMap."""

    def __init__(self, value):
        self.value = value
        msg = 'The value "{}" is already in the mapping.'
        super().__init__(msg.format(value))


class BidirectionalMap(dict):
    """Invertible map."""

    def __init__(self, inverse=None):
        if inverse is None:
            inverse = self.__class__(inverse=self)
        self.inverse = inverse

    def __setitem__(self, key, value):
        if value in self.inverse:
            raise BidirectionalError(value)

        self.inverse._set_item(value, key)
        self._set_item(key, value)

    def __delitem__(self, key):
        self.inverse._del_item(self[key])
        self._del_item(key)

    def _del_item(self, key):
        super().__delitem__(key)

    def _set_item(self, key, value):
        super().__setitem__(key, value)


# -------------------------------------------------------------------
# Get lengths and contig seqs
# -------------------------------------------------------------------


def get_edge_lengths(edge_file):

    contig_lengths = {}
    graph_contigs = {}

    for index, record in enumerate(SeqIO.parse(edge_file, "fasta")):
        contig_lengths[record.id] = len(record.seq)
        graph_contigs[record.id] = record.seq

    return contig_lengths, graph_contigs


# -------------------------------------------------------------------
# Unicycler graph functions
# -------------------------------------------------------------------


def get_links_unicycler(assembly_graph_file):

    node_count = 0

    graph_contigs = {}

    edge_depths = {}

    links = []

    my_map = BidirectionalMap()

    # Get links from .gfa file
    with open(assembly_graph_file) as file:
        for line in file.readlines():

            # Identify lines with link information
            if line.startswith("L"):
                link = []

                strings = line.split("\t")

                link1 = strings[1]
                link2 = strings[3]

                link1_orientation = strings[2]
                link2_orientation = strings[4]

                if link1_orientation == "+" and link2_orientation == "+":
                    link.append(link1)
                    link.append(link2)

                elif link1_orientation == "-" and link2_orientation == "-":
                    link.append(link2)
                    link.append(link1)

                elif link1_orientation == "+" and link2_orientation == "-":
                    link.append(link1)
                    link.append(link2)

                elif link1_orientation == "-" and link2_orientation == "+":
                    link.append(link1)
                    link.append(link2)

                links.append(link)

            elif line.startswith("S"):
                strings = line.strip().split()

                my_map[node_count] = strings[1]

                graph_contigs[strings[1]] = strings[2]

                depth = float(strings[3].split(":")[-1])
                edge_depths[strings[1]] = depth

                node_count += 1

            line = file.readline()

    return node_count, links, my_map


def get_graph_edges_unicycler(links, contig_names_rev):

    self_looped_nodes = []

    edge_list = []

    # Iterate links
    for link in links:
        # Remove self loops
        if link[0] != link[1]:
            # Add edge to list of edges
            edge_list.append((contig_names_rev[link[0]], contig_names_rev[link[1]]))
        else:
            self_looped_nodes.append(contig_names_rev[link[0]])

    return edge_list, self_looped_nodes


# -------------------------------------------------------------------
# SPAdes graph functions
# -------------------------------------------------------------------


def get_segment_paths_spades(contig_paths):

    paths = {}
    segment_contigs = {}
    node_count = 0

    my_map = BidirectionalMap()

    contig_names = BidirectionalMap()

    current_contig_num = ""

    with open(contig_paths) as file:
        name = file.readline().strip()
        path = file.readline().strip()

        while name != "" and path != "":

            while ";" in path:
                path = f"{path[:-2]},{file.readline()}"

            start = "NODE_"
            end = "_length_"
            contig_num = str(int(re.search("%s(.*)%s" % (start, end), name).group(1)))

            segments = path.rstrip().split(",")

            if current_contig_num != contig_num:
                my_map[node_count] = int(contig_num)
                contig_names[node_count] = name.strip()
                current_contig_num = contig_num
                node_count += 1

            if contig_num not in paths:
                paths[contig_num] = segments

            for segment in segments:

                if segment not in segment_contigs:
                    segment_contigs[segment] = set([contig_num])
                else:
                    segment_contigs[segment].add(contig_num)

            name = file.readline().strip()
            path = file.readline().strip()

    return paths, segment_contigs, node_count, my_map, contig_names


def get_graph_edges_spades(
    assembly_graph_file, contigs_map, contigs_map_rev, paths, segment_contigs
):

    links = []
    links_map = defaultdict(set)

    # Get links from assembly_graph_with_scaffolds.gfa
    with open(assembly_graph_file) as file:
        line = file.readline()

        while line != "":

            # Identify lines with link information
            if "L" in line:
                strings = line.split("\t")
                f1, f2 = f"{strings[1]}{strings[2]}", f"{strings[3]}{strings[4]}"
                links_map[f1].add(f2)
                links_map[f2].add(f1)
                links.append(f"{strings[1]}{strings[2]} {strings[3]}{strings[4]}")
            line = file.readline()

    # Create list of edges
    edge_list = []

    self_looped_nodes = []

    for i in range(len(paths)):
        segments = paths[str(contigs_map[i])]

        new_links = []

        for segment in segments:

            my_segment = segment

            my_segment_rev = ""

            if my_segment.endswith("+"):
                my_segment_rev = f"{my_segment[:-1]}-"
            else:
                my_segment_rev = f"{my_segment[:-1]}+"

            if segment in links_map:
                new_links.extend(list(links_map[segment]))

            if my_segment_rev in links_map:
                new_links.extend(list(links_map[my_segment_rev]))

        if my_segment in segment_contigs:
            for contig in segment_contigs[my_segment]:
                if i != contigs_map_rev[int(contig)]:
                    # Add edge to list of edges
                    edge_list.append((i, contigs_map_rev[int(contig)]))

        if my_segment_rev in segment_contigs:
            for contig in segment_contigs[my_segment_rev]:
                if i != contigs_map_rev[int(contig)]:
                    # Add edge to list of edges
                    edge_list.append((i, contigs_map_rev[int(contig)]))

        for new_link in new_links:
            if new_link in segment_contigs:
                for contig in segment_contigs[new_link]:
                    if i != contigs_map_rev[int(contig)]:
                        # Add edge to list of edges
                        edge_list.append((i, contigs_map_rev[int(contig)]))
                    else:
                        self_looped_nodes.append(i)

    return edge_list, self_looped_nodes


# -------------------------------------------------------------------
# Flye graph functions
# -------------------------------------------------------------------


def get_flye_contig_map(contigs_file):

    contig_names = BidirectionalMap()

    contig_num = 0

    for index, record in enumerate(SeqIO.parse(contigs_file, "fasta")):
        contig_names[contig_num] = record.id
        contig_num += 1

    return contig_names


def get_links_flye(contig_paths, contig_names_rev):

    paths = {}
    segment_contigs = {}

    my_map = BidirectionalMap()

    with open(contig_paths) as file:

        for line in file.readlines():

            if not line.startswith("#"):

                strings = line.strip().split()

                contig_name = strings[0]

                path = strings[-1]
                path = path.replace("*", "")

                if path.startswith(","):
                    path = path[1:]

                if path.endswith(","):
                    path = path[:-1]

                segments = path.rstrip().split(",")

                contig_num = contig_names_rev[contig_name]

                if contig_num not in paths:
                    paths[contig_num] = segments

                for segment in segments:

                    if segment not in segment_contigs:
                        segment_contigs[segment] = set([contig_num])
                    else:
                        segment_contigs[segment].add(contig_num)

    return paths, segment_contigs, len(contig_names_rev), my_map


def get_graph_edges_flye(
    assembly_graph_file, contigs_map, contigs_map_rev, paths, segment_contigs
):

    links_map = defaultdict(set)

    # Get links from assembly_graph_with_scaffolds.gfa
    with open(assembly_graph_file) as file:
        line = file.readline()

        while line != "":

            # Identify lines with link information
            if "L" in line:
                strings = line.split("\t")

                f1, f2 = "", ""

                if strings[2] == "+":
                    f1 = strings[1][5:]
                if strings[2] == "-":
                    f1 = f"-{strings[1][5:]}"
                if strings[4] == "+":
                    f2 = strings[3][5:]
                if strings[4] == "-":
                    f2 = f"-{strings[3][5:]}"

                links_map[f1].add(f2)
                links_map[f2].add(f1)

            line = file.readline()

    # Create list of edges
    edge_list = []

    self_looped_nodes = []

    for i in paths:
        segments = paths[i]

        new_links = []

        for segment in segments:

            my_segment = segment
            my_segment_num = ""

            my_segment_rev = ""

            if my_segment.startswith("-"):
                my_segment_rev = my_segment[1:]
                my_segment_num = my_segment[1:]
            else:
                my_segment_rev = f"-{my_segment}"
                my_segment_num = my_segment

            if my_segment in links_map:
                new_links.extend(list(links_map[my_segment]))

            if my_segment_rev in links_map:
                new_links.extend(list(links_map[my_segment_rev]))

            if my_segment in segment_contigs:
                for contig in segment_contigs[my_segment]:
                    if i != contig:
                        # Add edge to list of edges
                        edge_list.append((i, contig))

            if my_segment_rev in segment_contigs:
                for contig in segment_contigs[my_segment_rev]:
                    if i != contig:
                        # Add edge to list of edges
                        edge_list.append((i, contig))

            if my_segment_num in segment_contigs:
                for contig in segment_contigs[my_segment_num]:
                    if i != contig:
                        # Add edge to list of edges
                        edge_list.append((i, contig))

        for new_link in new_links:

            if new_link in segment_contigs:
                for contig in segment_contigs[new_link]:
                    if i != contig:
                        # Add edge to list of edges
                        edge_list.append((i, contig))
                    else:
                        self_looped_nodes.append(i)

            if new_link.startswith("-"):
                if new_link[1:] in segment_contigs:
                    for contig in segment_contigs[new_link[1:]]:
                        if i != contig:
                            # Add edge to list of edges
                            edge_list.append((i, contig))
                        else:
                            self_looped_nodes.append(i)

    return edge_list, self_looped_nodes


# -------------------------------------------------------------------
# MEGAHIT graph functions
# -------------------------------------------------------------------


def get_links_megahit(assembly_graph_file):

    node_count = 0

    graph_contigs = {}

    links = []

    my_map = BidirectionalMap()

    line_number = 0

    # Get links from .gfa file
    with open(assembly_graph_file) as file:

        line = file.readline()

        while line != "":

            if line_number % 4 == 0 and line.startswith(">"):

                # Identify contigs with link information
                if ":" in line:

                    strings = line.strip().split(":")

                    if strings[0].endswith("'"):
                        link1 = strings[0][1:-1]
                    else:
                        link1 = strings[0][1:]

                    if "," in strings[1]:
                        second_links = strings[1][:-1].split(",")

                        for second_link in second_links:

                            link = []

                            if second_link.endswith("'"):
                                link2 = second_link[:-1]
                            else:
                                link2 = second_link

                            link.append(link1)
                            link.append(link2)
                            links.append(link)

                # Identify main contigs
                if ":" in line:
                    contig_name = line.strip().split(":")[0][1:]
                else:
                    contig_name = line.strip()[1:-1]

                line = file.readline()
                line_number += 1
                contig_seq = line.strip()

                my_map[node_count] = contig_name
                graph_contigs[contig_name] = contig_seq

                node_count += 1

            elif line_number % 2 == 0 and ":" in line:

                # Identify contigs with link information

                strings = line.strip().split(":")

                if strings[0].endswith("'"):
                    link1 = strings[0][1:-1]
                else:
                    link1 = strings[0][1:]

                if "," in strings[1]:
                    second_links = strings[1][:-1].split(",")

                    for second_link in second_links:

                        link = []

                        if second_link.endswith("'"):
                            link2 = second_link[:-1]
                        else:
                            link2 = second_link

                        link.append(link1)
                        link.append(link2)
                        links.append(link)

            line_number += 1
            line = file.readline()

    return node_count, graph_contigs, links, my_map


def get_graph_edges_megahit(links, contig_names_rev):

    edge_list = []

    self_looped_nodes = []

    # Iterate links
    for link in links:
        # Remove self loops
        if link[0] != link[1]:
            # Add edge to list of edges
            edge_list.append((contig_names_rev[link[0]], contig_names_rev[link[1]]))
        else:
            self_looped_nodes.append(contig_names_rev[link[0]])

    return edge_list, self_looped_nodes


# -------------------------------------------------------------------
# Main graph building function
# -------------------------------------------------------------------


def build_assembly_graph(
    assembler,
    assembly_graph_file,
    contigs_file=None,
    contig_paths_file=None,
    is_directed=False,
):
    graph_contigs = []
    graph_to_contig_map = BidirectionalMap()

    if assembler == "unicycler":

        # Get unicycler links
        node_count, links, contig_names = get_links_unicycler(assembly_graph_file)

        # Get reverse mapping of contig identifiers
        contig_names_rev = contig_names.inverse

        # Get list of edges
        edge_list, self_looped_nodes = get_graph_edges_unicycler(
            links=links, contig_names_rev=contig_names_rev
        )

    elif assembler == "spades":

        # Get paths, segments, links and contigs of the assembly graph
        (
            paths,
            segment_contigs,
            node_count,
            contigs_map,
            contig_names,
        ) = get_segment_paths_spades(contig_paths_file)

        # Get reverse mapping of contig map
        contigs_map_rev = contigs_map.inverse

        # Get reverse mapping of contig identifiers
        contig_names_rev = contig_names.inverse

        # Get list of edges
        edge_list, self_looped_nodes = get_graph_edges_spades(
            assembly_graph_file=assembly_graph_file,
            contigs_map=contigs_map,
            contigs_map_rev=contigs_map_rev,
            paths=paths,
            segment_contigs=segment_contigs,
        )

    elif assembler == "flye":

        # Get contigs map
        contig_names = get_flye_contig_map(contigs_file)

        # Get reverse mapping of contig identifiers
        contig_names_rev = contig_names.inverse

        # Get links and contigs of the assembly graph
        (
            paths,
            segment_contigs,
            node_count,
            contigs_map,
        ) = get_links_flye(contig_paths_file, contig_names_rev)

        # Get reverse mapping of contig map
        contigs_map_rev = contigs_map.inverse

        # Get list of edges
        edge_list, self_looped_nodes = get_graph_edges_flye(
            assembly_graph_file=assembly_graph_file,
            contigs_map=contigs_map,
            contigs_map_rev=contigs_map_rev,
            paths=paths,
            segment_contigs=segment_contigs,
        )

    elif assembler == "megahit":

        original_contigs = {}
        contig_descriptions = {}

        # Get mapping of original contig identifiers with descriptions
        for index, record in enumerate(SeqIO.parse(contigs_file, "fasta")):
            original_contigs[record.id] = str(record.seq)
            contig_descriptions[record.id] = record.description

        # Get links and contigs of the assembly graph
        (
            node_count,
            graph_contigs,
            links,
            contig_names,
        ) = get_links_megahit(assembly_graph_file)

        # Get reverse mapping of contig identifiers
        contig_names_rev = contig_names.inverse

        # Map original contig identifiers to contig identifiers of MEGAHIT assembly graph
        graph_to_contig_map = BidirectionalMap()

        for (n, m), (n2, m2) in zip(graph_contigs.items(), original_contigs.items()):
            if m == m2:
                graph_to_contig_map[n] = n2

        # Get list of edges
        edge_list, self_looped_nodes = get_graph_edges_megahit(
            links=links, contig_names_rev=contig_names_rev
        )

    # Create graph
    assembly_graph = Graph(directed=is_directed)

    # Add vertices
    assembly_graph.add_vertices(node_count)

    # Name vertices with contig identifiers
    for i in range(node_count):
        assembly_graph.vs[i]["id"] = i
        assembly_graph.vs[i]["name"] = contig_names[i]
        assembly_graph.vs[i]["label"] = f"{contig_names[i]}\nID:{str(i)}"

    # Add edges to the graph
    assembly_graph.add_edges(edge_list)

    # Simplify the graph
    assembly_graph.simplify(multiple=True, loops=False, combine_edges=None)

    return (
        assembly_graph,
        contig_names,
        contig_names_rev,
        graph_to_contig_map,
        self_looped_nodes,
    )
