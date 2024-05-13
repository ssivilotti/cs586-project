import networkx as nx
from ete3 import Tree
from cassiopeia.data import CassiopeiaTree
from cassiopeia.tools import small_parsimony as sp
import pandas as pd


# for each tree:
# read in tree
# exclude if leaves all in same tumor
# infer metastasis labels, without adding new nodes
# new tree has node labels of loc.name
# use inferred labels to resolve polytomies by adding new nodes
# write resolved tree to new file in newick format (same as input file)
# compute metastasis score for initial tree and resolved tree
    # number of edges where metastsis occurs / total number of edges (use original number of edges for resolved tree)
def compute_met_score(tree: CassiopeiaTree, resolved_tree: CassiopeiaTree = None) -> float:
    # compute metastasis score
    met_edges = 0
    total_edges = len(tree.edges)

    for node in tree.nodes:
        # check if there is a metastasis event between node and children
        if tree.is_leaf(node):
            continue
        
        try:
            node_location = tree.get_attribute(node, 'met_location')
        except:
            node_location = node.split('.')[0]
        
        children = tree.children(node)
        if resolved_tree is not None:
            children = resolved_tree.children(node)

        for child in children:
            try:
                child_location = tree.get_attribute(child, 'met_location')
            except:
                child_location = child.split('.')[0]
            
            if node_location != child_location:
                met_edges += 1

    return met_edges / total_edges

def resolve_polytomies(tree: CassiopeiaTree, editable_tree: Tree) -> CassiopeiaTree:
    for node in tree.nodes:
        # check if there is a metastasis event between node and children
        if tree.is_leaf(node) or len(tree.children(node)) <= 2:
            continue
        
        node_location = tree.get_attribute(node, 'met_location')

        met_nodes = {key: [] for key in ["LL", "RE", "RW", "M1", "M2", "Liv"] if key != node_location}
        
        for child in tree.children(node):
            child_location = tree.get_attribute(child, 'met_location')
            if node_location != child_location:
                met_nodes[child_location].append(child)
        tree_node:Tree = editable_tree.search_nodes(name=node)[0]
        for key in met_nodes:
            if len(met_nodes[key]) > 1:
                poly_resolve_node = tree_node.add_child(name=f'{key}.{node}')
                for child in met_nodes[key]:
                    child_node = tree_node.search_nodes(name=child)[0].detach()
                    poly_resolve_node.add_child(child_node)

    resolved_tree = CassiopeiaTree(tree=editable_tree, cell_meta=tree.cell_meta)
    return resolved_tree


            

origin_tumor = 'LL'

scores = {}

trees = {}

for clonal_population in range(1,101):
    # read in tree from file
    try:
        with open(f'GSE161363_RAW/trees/m5k_lg{clonal_population}_tree_hybrid_priors.alleleThresh.processed.txt', 'r') as f:
            newick_str = f.read().strip()
    except FileNotFoundError:
        print(f'File not found: GSE161363_RAW/trees/m5k_lg{clonal_population}_tree_hybrid_priors.alleleThresh.processed.txt')
        continue

    try:

        t = Tree(newick_str, format=1)

        # convert to CassiopeiaTree

        # handle duplicate node names
        count = 0
        all_nodes_already_in_tree = set()
        # iterate through trees and rename edges
        for node in t.traverse():
            if node.name in all_nodes_already_in_tree:
                node.name = f"{node.name}_{count}"
                count += 1
            all_nodes_already_in_tree.add(node.name)

        leaves = [node.name for node in t.get_leaves()]

        cell_met_labels = [leaf.split('.')[0] for leaf in leaves]

        cell_meta_data = pd.DataFrame({'met_labels': cell_met_labels}, index=leaves).astype("category")

        cass_tree = CassiopeiaTree(tree=t, cell_meta=cell_meta_data)

        # compute metastasis labels

        fitch_tree = sp.fitch_hartigan(cass_tree, meta_item='met_labels', label_key='met_location', copy=True)

        resolved_tree = resolve_polytomies(fitch_tree, t)

        # compute scores

        original_met_score = compute_met_score(fitch_tree)

        resolved_met_score = compute_met_score(fitch_tree, resolved_tree)

        scores[clonal_population] = (original_met_score, resolved_met_score)

        trees[clonal_population] = (fitch_tree, resolved_tree)

        # print to file

        with open(f'resolvedtrees/m5k_lg{clonal_population}_tree_hybrid_priors.alleleThresh.processed_resolved.txt', 'w+') as f:
            f.write(resolved_tree.get_newick())
    except Exception as e:
        with open(f'resolvedtrees/error.txt', 'a+') as f:
            f.write(str(e)+'\n')
    
with open(f'resolvedtrees/score.txt', 'w+') as f:
    f.write(str(scores))


from ete3 import faces, AttrFace, TreeStyle, NodeStyle

met_positions = ["LL", "RE", "RW", "M1", "M2", "Liv"]
colors = ["Pink", "Salmon", "Crimson", "LawnGreen", "LimeGreen", "DodgerBlue"]

postion_to_style = {}

for i, pos in enumerate(met_positions):
    postion_to_style[pos] = NodeStyle()
    postion_to_style[pos]["fgcolor"] = colors[i]
    # postion_to_style[pos]["size"] = 50
    postion_to_style[pos]["vt_line_color"] = colors[i]
    postion_to_style[pos]["hz_line_color"] = colors[i]
    postion_to_style[pos]["hz_line_width"] = 15
    postion_to_style[pos]["vt_line_width"] = 15

ete_styled_trees = {}

def cass_to_ete3(tree: CassiopeiaTree, cass_t: CassiopeiaTree) -> Tree:
    t = Tree(f'{tree.root};', format=1)
    n = Tree()
    for parent, child in tree.breadth_first_traverse_edges():
        if n.name != parent:
            n = t.search_nodes(name=parent)[0]
            try:
                node_location = cass_t.get_attribute(parent, 'met_location')
            except:
                node_location = parent.split('.')[0]
            print(node_location)
            n.set_style(postion_to_style[node_location])
        n.add_child(name=child)
        if tree.is_leaf(child):
            try:
                node_location = cass_t.get_attribute(child, 'met_location')
            except:
                node_location = child.split('.')[0]
            n.search_nodes(name=child)[0].set_style(postion_to_style[node_location])
    return t

for key, val in trees.items():
    cass_t:CassiopeiaTree = val[0]
    cass_t_resolved:CassiopeiaTree = val[1]

    t = cass_to_ete3(cass_t, cass_t)
    t_resolved = cass_to_ete3(cass_t_resolved, cass_t)

    ete_styled_trees[key] = (t, t_resolved)

    ts = TreeStyle()
    ts.optimal_scale_level = 'full'
    ts.guiding_lines_type = 1
    ts.show_leaf_name = False
    ts.mode = "c"
    ts.root_opening_factor = 1

    val[0].render(f"resolvedtrees/cp{key}.png", w=1700, tree_style=ts)

    def layout(node):
        if node.is_leaf():
            N = AttrFace("name", fsize=5)
            faces.add_face_to_node(N, node, 0, position="aligned")

    ts = TreeStyle()
    ts.optimal_scale_level = 'full'
    ts.guiding_lines_type = 1
    ts.show_leaf_name = False
    ts.mode = "c"
    ts.root_opening_factor = 1
    
    val[1].render(f"resolvedtrees/cp{key}_resolved.png", w=1700, tree_style=ts)




def compute_sc_met_score(tree: CassiopeiaTree, cass_t: CassiopeiaTree) -> float:
    # compute metastasis score
    
    saved_node_info = {}

    # compute TreeMetRate for each subclade
    for node in tree.depth_first_traverse_nodes():
        met_edges = 0
        total_edges = 0
        # check if there is a metastasis event between node and children
        if tree.is_leaf(node):
            saved_node_info[node] = (0, 0, 0)
            continue
        
        try:
            node_location = cass_t.get_attribute(node, 'met_location')
        except:
            node_location = node.split('.')[0]
            saved_node_info[node] = None
            continue
        
        children = tree.children(node)

        for child in children:
            try:
                child_location = cass_t.get_attribute(child, 'met_location')
                met_edges += saved_node_info[child][1]
                total_edges += saved_node_info[child][2] + 1
            except:
                child_location = child.split('.')[0]
                for c2 in tree.children(child):
                    met_edges += saved_node_info[c2][1]
                    total_edges += saved_node_info[c2][2] + 1
            
            if node_location != child_location:
                met_edges += 1
        saved_node_info[node] = (met_edges / total_edges, met_edges, total_edges)
    
    scMetRates = {}
    for node in tree.leaves:
        ancestors = tree.get_all_ancestors(node)
        all_tmrs = [saved_node_info[ancestor][0] for ancestor in ancestors if saved_node_info[ancestor] != None]
        scMetRates[node] = sum(all_tmrs) / len(all_tmrs)

    return scMetRates

scMetRates = {}
for key, val in trees.items():
    # if key != 7:
    #     continue
    cass_t:CassiopeiaTree = val[0]
    cass_t_resolved:CassiopeiaTree = val[1]
    scMetRate = compute_sc_met_score(cass_t, cass_t)
    scMetRate_resolved = compute_sc_met_score(cass_t_resolved, cass_t)
    scMetRates[key] = (scMetRate, scMetRate_resolved)

import pickle

with open('resolvedtrees/resolved_trees_fixed.pkl', 'wb') as f:
    pickle.dump(trees, f)

with open('resolvedtrees/scores_fixed.pkl', 'wb') as f:
    pickle.dump(scores, f)

with open('resolvedtrees/ete_styled_trees_fixed.pkl', 'wb') as f:
    pickle.dump(ete_styled_trees, f)

with open('resolvedtrees/scMetRates.pkl', 'wb') as f:
    pickle.dump(scMetRates, f)