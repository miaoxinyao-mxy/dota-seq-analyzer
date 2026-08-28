#!/usr/bin/env python3

class Node:
    """Defines the Node units of the taxonomy tree, each corresponding to one observed value of one taxonomic rank"""
    
    def __init__(self, tax_name: str, tax_lvl_num: int, tax_lvl_char: str, parent=None, tax_id=None):
        self.tax_name = tax_name
        self.tax_lvl_num = tax_lvl_num
        self.tax_lvl_char = tax_lvl_char
        self.parent = parent
        self.children = []
        # 2026-08-10: Keep Kraken2 taxids on taxonomy nodes for standard output parsing.
        # Reason: Kraken2 output supplies taxids, while the report supplies names and ranks.
        self.tax_id = tax_id

    def add_child(self, child):
        """Add parameter Node as a child to the current Node"""
        if isinstance(child, Node):
            self.children.append(child)
        else:
            print("Error - child is not an instance of Node")

    def __str__(self):
        return f"({self.tax_name}, {self.tax_lvl_char})"

def parse_report_line(r_line: str):
    """Parse the taxonomic name, level, & tax_ID from a line of the kraken report file"""

    # match taxonomic level to level number
    # note kingdom rank is not included, since only analyzing bacteria
    match_lvls = {'U':-1, 'R':0, 'D':1, 'R1':1, 'P':2, 'C':3, 'O':4, 'F':5, 'G':6, 'S':7}    

    r_line_data = r_line.split("\t")
    tax_lvl_char = r_line_data[3].strip()

    # 2026-08-10: Accept Kraken2's D (domain) rank as the project's R1 rank.
    # Reason: Kraken2 emits D, while downstream taxonomy objects use R1.
    if tax_lvl_char == "D":
        tax_lvl_char = "R1"
    tax_lvl_num = match_lvls[tax_lvl_char] # determine tax level number based on its letter character
    
    # deal with "unclassified" and "root" lines, which are edge cases
    if (tax_lvl_num == -1) or (tax_lvl_num == 0):
        tax_name = r_line_data[5].strip()
    else:
        tax_name = r_line_data[5].strip().split("__")[1]

    return tax_name, tax_lvl_num, tax_lvl_char, r_line_data[4].strip()


def create_taxonomy_tree(report_filename: str):
    """
    Create a taxonomy tree, based on the report file outputted by the Kraken2 analysis.
    Purpose: this will help later on in determining the full taxonomic path of a 16s read. 
     Specifically, given the Kraken2 classification for that read - which would only say the value of the most specific taxonomic rank 
     (e.g. "Genus - Staphylococcus") - we could traverse from the given node up the tree to obtain the full taxonomic path.
    Structure of tree:
      Each observed value of each taxonomic rank is represented by a Node object.
      In constructing the tree, two nodes will be connected if one is a direct descendant of the other 
       (e.g. one node represents a species that is a descendant of the genus represented by the other node).
      Nodes on the same "level" of the tree have the same taxonomic rank.
      The taxonomic ranks considered here are: domain (D), phylum (P), class (C), order (O), family (F), genus (G), species (S).
    """

    # create lists to store all nodes of a given taxonomy level (e.g. all species-level nodes, genus-level nodes, etc.)
    # to be used later in downstream processing
    nodes_lists = {'R1': [], 'P': [], 'C': [], 'O': [], 'F': [], 'G': [], 'S': []}
    # 2026-08-10: Add a taxid lookup for standard Kraken2 classification output.
    # Reason: create_ID_packets receives taxids, not taxonomy names.
    nodes_lists["_taxid"] = {}

    # preliminary steps
    r_file = open(report_filename, 'r')
    prev_node = None

    # iterate through each line of the report file, each of will be represented by a Node object
    for line in r_file:
        # parsing
        tax_name, tax_lvl_num, tax_lvl_char, tax_id = parse_report_line(line)

        # deal with "unclassified" and "root" lines, which are edge cases
        if tax_lvl_num == -1:
            continue
        elif tax_lvl_num == 0:
            root_node = Node(tax_name, tax_lvl_num, tax_lvl_char, tax_id=tax_id)
            prev_node = root_node
            continue

        # 2026-08-10: Attach a taxon to the nearest shallower rank, allowing skipped ranks.
        # Reason: Kraken2 reports can omit ranks, so requiring consecutive levels can dereference None.
        while prev_node is not None and prev_node.tax_lvl_num >= tax_lvl_num:
            prev_node = prev_node.parent
        if prev_node is None:
            raise ValueError("Kraken report has no valid parent")

        # create new node, and connect the parent & child node
        current_node = Node(tax_name, tax_lvl_num, tax_lvl_char, prev_node, tax_id)
        prev_node.add_child(current_node)

        # add node to correct node list
        nodes_lists[tax_lvl_char].append(current_node)
        nodes_lists["_taxid"][tax_id] = current_node

        # preparing for iteration over next line
        prev_node = current_node

    r_file.close()

    return nodes_lists
