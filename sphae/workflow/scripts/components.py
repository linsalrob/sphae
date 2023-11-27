#!/usr/bin/python3

"""components.py: Obtain the sequences corresponding to edges in components of the SPAdes, Flye and Unicycler assembly graphs.

The assembly graph file  should be provided as inputs.

"""

# import click
import os
from graph_utils import build_utils

__author__ = "Vijini Mallawaarachchi"
__copyright__ = "Copyright 2019, GraphBin Project"
__credits__ = "Vijini Mallawaarachchi"
__license__ = "BSD-3"
__type__ = "Support Script"
__maintainer__ = "Vijini Mallawaarachchi"
__email__ = "viji.mallawaarachchi@gmail.com"

# Sample command
# -------------------------------------------------------------------
# python components.py  -a assembler_name
#                       -g /path/to/folder_with_binning_result
#                       -c /path/to/contigs
#                       -p /path/to/contig_paths
#                       -o /path/to/output_folder
# -------------------------------------------------------------------


# @click.command()
# @click.option(
#     "--assembler",
#     "-a",
#     required=True,
#     help="assembler name (flye, spades or unicycler)",
#     type=click.Choice(["flye", "spades", "megahit", "unicycler"], case_sensitive=False),
# )
# @click.option(
#     "--graph",
#     "-g",
#     required=True,
#     help="path to the assembly graph file",
#     type=click.Path(exists=True),
# )
# @click.option(
#     "--contigs",
#     "-c",
#     required=True,
#     help="path to the contigs file",
#     type=click.Path(exists=True),
# )
# @click.option(
#     "--paths",
#     "-p",
#     required=False,
#     help="path to the contig paths file",
#     type=click.Path(exists=True),
# )
# @click.option(
#     "--output",
#     "-o",
#     required=True,
#     help="path to the output folder",
#     type=click.Path(exists=True),
# )
def main(assembler, graph, contigs, path, output):

    # Get contig lengths
    # -------------------------------------------------------------------

    contig_lengths, graph_contigs = build_utils.get_edge_lengths(contigs)

    # Build the assembly graph
    # -------------------------------------------------------------------

    (
        assembly_graph,
        contig_names,
        contig_names_rev,
        graph_to_contig_map,
        self_looped_nodes,
    ) = build_utils.build_assembly_graph(
        assembler, graph, contigs, path, is_directed=False
    )

    # Get reverse contig mapping for megahit

    if assembler == "megahit":
        graph_to_contig_map_rev = graph_to_contig_map.inverse

    # Get circlar contigs
    # -------------------------------------------------------------------

    circular = []

    # print(self_looped_nodes)

    for contig_num in self_looped_nodes:

        if len(assembly_graph.neighbors(contig_num)) == 0:
            circular.append(contig_num)

    # Get non-isolated contigs with no neighbours
    # -------------------------------------------------------------------

    non_isolated = [v.index for v in assembly_graph.vs if v.degree() != 0]

    # Write non isolated and circular contigs and their details
    # -------------------------------------------------------------------

    to_write = non_isolated + circular

    if len(to_write) != 0:

        # Save details of important contigs
        with open(f"{output}/graph_seq_details_{assembler}.tsv", "w") as myfile:

            myfile.write(f"ContigID\tLength\tCircular\tConnections\n")
            for contig in [v.index for v in assembly_graph.vs]:
                is_circular = contig in circular

                if assembler == "megahit":
                    contig_identifier = graph_to_contig_map[contig_names[contig]]
                else:
                    contig_identifier = contig_names[contig]

                length = contig_lengths[contig_identifier]
                n_connections = len(assembly_graph.neighbors(contig))
                myfile.write(
                    f"{contig_identifier}\t{length}\t{is_circular}\t{n_connections}\n"
                )

        # Save important contigs and the rest in two separate FASTA files
        with open(f"{output}/important_seqs_{assembler}.fasta", "w") as myfile1:

            with open(f"{output}/other_seqs_{assembler}.fasta", "w") as myfile2:

                for contig in graph_contigs:

                    if assembler == "megahit":
                        contig_id = contig_names_rev[graph_to_contig_map_rev[contig]]
                    else:
                        contig_id = contig_names_rev[contig]

                    contig_seq = graph_contigs[contig]

                    if contig_id in to_write:
                        myfile1.write(f">{contig}\n{contig_seq}\n")
                    else:
                        myfile2.write(f">{contig}\n{contig_seq}\n")

    else:
        print("No circular contigs or components were found")

        with open(f"{output}/graph_seq_details_{assembler}.tsv", "w") as myfile:
            myfile.write(f"ContigID\tLength\tCircular\tConnections\n")
            for contig in [v.index for v in assembly_graph.vs]:

                if assembler == "megahit":
                    contig_identifier = graph_to_contig_map[contig_names[contig]]
                else:
                    contig_identifier = contig_names[contig]

                length = contig_lengths[contig_identifier]
                n_connections = len(assembly_graph.neighbors(contig))
                is_circular = False
                myfile.write(
                    f"{contig_identifier}\t{length}\t{is_circular}\t{n_connections}\n"
                )

        print("Exiting components.py... Bye...!")

    print("Thanks for using components.py... Bye...!")


main(
    snakemake.params.assembler,
    snakemake.input.graph,
    snakemake.input.contigs,
    snakemake.input.path,
    snakemake.params.o,
)