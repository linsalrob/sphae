import argparse
import os
import string
import sys

__author__ = 'Rob Edwards'

"""

Newick format:
((15:0.04110,((16:0.03869,17:0.11891):0.00888,14:0.10566):0.00609):0.09345,
(((6:0.00000,7:0.00000):0.01537,13:0.00839):0.00229,((1:0.00011,
4:0.00065):0.01946,(11:-0.01252,(10:0.00692,(8:0.00463,(((2:0.00178,
(3:0.00192,5:0.00082):0.00081):0.00254,9:0.00327):0.00157,
12:-0.00216):0.00082):0.00088):0.01351):0.05366):0.01353):0.01429,0:0.00383);

This is a parser that I wrote myself. (Mainly so that other people don't have to install biopython or something
similar). It works for vanilla newick trees, but not for more complex trees (e.g. that have trifurcating branches).
Use at your own risk!

"""


class Node(object):
    """A node object"""
    def __init__(self, id):
        self.id = id
        self.left = None
        self.right = None
        self.parent = None
        self.distance = ""
        self.name = ""
        self.side = None


class Newick_Tree(object):

    def __init__(self):
        pass

    def count_nodes(self, root):
        """ Count the number of nodes in the tree
        :param root: The root node
        :type root: Node
        :return: The number of nodes
        :rtype: int
        """

        def count(node):
            c = 1
            if node.left:
                c += count(node.left)
            if node.right:
                c += count(node.right)
            return c

        return count(root)

    def parse(self, tree, verbose=False):
        """
        Parse the string given by tree
        :param tree: the string to parse
        :param verbose: whether to make lots of output
        :return:
        """
        def process_tree(treestr, pos, node):
            # sys.stderr.write("At pos {} tree has depth {}\n".format(pos, Tree().count_nodes(root)))

            if treestr[pos] == '(':
                pos += 1
                # sys.stderr.write("LEFT: When adding node {} tree has depth {}\n".format(pos, Tree().count_nodes(root)))
                newnode = Node(pos)
                newnode.parent = node
                node.left = newnode
                newnode.side = "Left"
                # sys.stderr.write("ADDED NODE {} to the LEFT\n".format(newnode.id))
                return process_tree(treestr, pos, newnode)
            elif treestr[pos] == ',':
                pos += 1
                # sys.stderr.write("RIGHT: When adding node {} tree has depth {}\n".format(pos, Tree().count_nodes(root)))
                newnode = Node(pos)
                if node.parent.right:
                    newnode = node.parent
                else:
                    newnode.parent = node.parent
                    node.parent.right = newnode
                    newnode.side = 'Right'
                # sys.stderr.write("ADDED NODE {} to the RIGHT\n".format(newnode.id))
                return process_tree(treestr, pos, newnode)
            elif treestr[pos] == ')':
                pos += 1
                if pos >= len(treestr):
                    return
                while treestr[pos] in string.ascii_letters or treestr[pos] == '_' or treestr[pos] in string.digits or treestr[pos] in '-':
                    node.name += treestr[pos]
                    pos += 1
                if verbose:
                    sys.stderr.write("At pos {} set node {} ({}) to {}\n".format(pos, node.id, node.side, node.name))
                return process_tree(treestr, pos, node.parent)
            elif treestr[pos] == ':':
                pos += 1
                try:
                    while treestr[pos] in string.digits or treestr[pos] == '.' or treestr[pos] in ['-', 'e']:
                        node.distance += treestr[pos]
                        pos += 1
                except TypeError:
                    raise TypeError("TypeError: CANNOT ADD {} at POS {} to {} in node {}\n".format(treestr[pos], pos, node.distance, node.id))
                if verbose:
                    sys.stderr.write("Set NODE {} dist to {}\n".format(node.id, node.distance))
                node.distance = float(node.distance)
                return process_tree(treestr, pos, node)
            else:
                while treestr[pos] in string.ascii_letters or treestr[pos] == '_' or treestr[pos] in string.digits or \
                                treestr[pos] in '-':
                    node.name += treestr[pos]
                    pos += 1
                # sys.stderr.write("When adding node {} tree has depth {}\n".format(node.name, Tree().count_nodes(root)))
                if verbose:
                    sys.stderr.write("At pos {} Set node {} ({}) to {}\n".format(pos, node.id, node.side, node.name))
                return process_tree(treestr, pos, node)

        parent = Node("root")
        root = parent
        pos = 0
        tree = tree.rstrip(';')
        pos = process_tree(tree, pos, parent)
        if verbose:
            sys.stderr.write("TREE HAS DEPTH {}\n".format(Tree().count_nodes(root)))
        return parent

    def print_tree(self, root):
        """
        Print out a tree in newick format.

        :param root: The root node of the tree
        :type root: Node
        :return:
        :rtype:
        """

        def process_child(node):
            toreturn = ''
            if node.left or node.right:
                toreturn = '('
            if node.left and node.right:
                toreturn += process_child(node.left) + "," + process_child(node.right)
            elif node.left:
                toreturn += process_child(node.left)
            elif node.right:
                toreturn += process_child(node.right)
            if node.left and node.right and node.name:
                # the root node??
                toreturn += ','
            elif node.left or node.right:
                toreturn += ")"
            if node.name:
                toreturn += node.name
            toreturn += ":{}".format(node.distance)
            return toreturn

        print(process_child(root) + ");")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Parse a tree')
    parser.add_argument('-t', help='tree file', required=True)
    parser.add_argument('-v', help='verbose output during parsing', action='store_true')
    args = parser.parse_args()

    tre = []
    with open(args.t, 'r') as f:
        for l in f:
            tre.append(l.strip())

    root = Tree().parse(''.join(tre), args.v)
    print("PARSED\n\nThere are {} nodes\n\n".format(Tree().count_nodes(root)))
    print("Root has left {} and right {}\n".format(root.left, root.right))
    Tree().print_tree(root)

