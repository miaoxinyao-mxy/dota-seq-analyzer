
class Node:
    """Defines the Node units of the taxonomy tree"""
    def __init__(self, tax_name: str, tax_lvl_num: int, tax_lvl_char: str, parent=None):
        self.tax_name = tax_name
        self.tax_lvl_num = tax_lvl_num
        self.tax_lvl_char = tax_lvl_char
        self.parent = parent
        self.children = []

    def add_child(self, child):
        if isinstance(child, Node):
            self.children.append(child)
        else:
            print("Error - child is not an instance of Node")

    def __str__(self):
        return f"({self.tax_name}, {self.tax_lvl_char})"

def parse_report_line(r_line: str):

    """Parse the taxonomic name & level from a line of the kraken report file"""

    # match taxonomic level to level number
    # note kingdom rank is not included, since only analyzing bacteria
    match_lvls = {'U':-1, 'R':0, 'R1':1, 'P':2, 'C':3, 'O':4, 'F':5, 'G':6, 'S':7}    

    r_line_data = r_line.split("\t")
    tax_lvl_char = r_line_data[3].strip()
    tax_lvl_num = match_lvls[tax_lvl_char]
    
    # deal with "unclassified" and "root" lines, which are edge cases
    if (tax_lvl_num == -1) or (tax_lvl_num == 0):
        tax_name = r_line_data[5].strip()
    else:
        tax_name = r_line_data[5].strip().split("__")[1]

    return tax_name, tax_lvl_num, tax_lvl_char


def create_taxonomy_tree(report_filename: str):

    """Create a taxonomy tree, consisting of interconnected Node units"""

    # create lists to store all nodes of a given taxonomy level (e.g. all species-level nodes, genus-level nodes, etc.)
    # to be used later in downstream processing
    nodes_lists = {'R1': [], 'P': [], 'C': [], 'O': [], 'F': [], 'G': [], 'S': []}

    r_file = open(report_filename, 'r')
    prev_node = None

    for line in r_file:
        tax_name, tax_lvl_num, tax_lvl_char = parse_report_line(line)

        # deal with "unclassified" and "root" lines, which are edge cases
        if tax_lvl_num == -1:
            continue
        elif tax_lvl_num == 0:
            root_node = Node(tax_name, tax_lvl_num, tax_lvl_char)
            prev_node = root_node
            continue

        # move node to correct parent
        while tax_lvl_num != prev_node.tax_lvl_num + 1:
            prev_node = prev_node.parent

        # create new node, and connect the parent & child node
        current_node = Node(tax_name, tax_lvl_num, tax_lvl_char, prev_node)
        prev_node.add_child(current_node)

        # add node to correct node list
        nodes_lists[tax_lvl_char].append(current_node)

        # preparing for iteration over next line
        prev_node = current_node

    r_file.close()

    return nodes_lists