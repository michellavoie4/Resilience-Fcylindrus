# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
#
# Python 3 code with UTF-8 as default encoding
#
# Functions to be used with the main script called 'main_network_analysis'
# Function 'load_graph()' creates a bipartite network from a metabolic network model generated in the cobra toolbox
#
#
# authors: Antoine Allard (www.antoineallard.info) and Michel Lavoie (https://github.com/michellavoie4/)
# 
#
#
# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=


# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
# Importing modules

import networkx as nx
import numpy as np
import collections

# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
# This function generates a network from a list of steady-state fluxes calculated with the cobra toolbox
# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
def load_graph(name_source_file):

    # =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
    # Read the list of links (edges) between metabolites/reactions ignoring edges having production or consumption rate = 0.

    # Create an empty  directed graph.
    g = nx.MultiDiGraph()

    # Initializing counters for network validation below.
    nb_ignored_edges = 0
    nb_lines = 0

    # Read the text file with all edges generated with MATLAB.
    with open(name_source_file, "r", encoding="utf-8") as f:

        # Loop on each line
        for line in f:

            # Remove the symbol (\n) and separate the line in chunks considering comma as separators.
            sline = line.strip().split(",")

            # Ignore the line if it is a comment
            if sline[0][0] != "#":

                # Count the number of lines in the source file for validation below.
                nb_lines += 1

                # Add a link to the graph if production/consumption flux is not nul.
                if sline[2] != "0":

                    # Add the first node to the graph (nothing is done if this node has already been added).
                    node_name_1 = sline[0]
                    type_de_noeud_1 = node_name_1.split("_")[0]
                    g.add_node(node_name_1, type=type_de_noeud_1)

                    # Add a second node in the network (nothing is done if this node has already been added).
                    node_name_2 = sline[1]
                    type_de_noeud_2 = node_name_2.split("_")[0]
                    g.add_node(node_name_2, type=type_de_noeud_2)

                    # Flux.
                    flux = float(sline[2])

                    # Add the edge to the graph.
                    g.add_edge(node_name_1, node_name_2, weight=flux)

                # Counting ignored edges with a nul flux for validation.
                else:
                    nb_ignored_edges += 1


    # =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
    # Validating that the function works properly

    # Check the potential presence of double links
    for n1 in g.nodes():

        neigh1 = list(g.neighbors(n1))

        decupl = [x for x in neigh1 if neigh1.count(x) > 1]

        for n2 in decupl:

            print("Double links: " + str((n1, n2)))

    # Make sure that the number of ignored links and the number of graph edges = total number of non-zero links found int he metabolic network
    if g.number_of_edges() + nb_ignored_edges - nb_lines != 0:
        print("There is a problem with the number of links added to the network.")


    # =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
    # Return a network object as a simple, ponderated, and directed bipartite graph
    return nx.DiGraph(g)

# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=






# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
# Unit tests
# Can be run if the script is executed as a main script

if __name__ == '__main__':

    # read the txt file with metabolite/reaction links and flows and generate a network.
    g = load_graph("Table_bipartite.txt")


    # Create a list of all network nodes
    liste_noeuds = g.nodes()

    # Total number of nodes (reactions + metabolites).
    nb_noeuds = g.number_of_nodes()

    # Number of links
    nb_liens = g.number_of_edges()

    # List of active reactions (reactions for which fluxes are not equal to 0).
    liste_react_actives = [n for n in g.nodes() if g.node[n]['type'] == 'R']

    # Number of active reactions
    nb_react_actives = len(liste_react_actives)

    # List of active metabolites (The flux of at least one reaction in which the metabolite appears is not equal to 0)
    liste_metabo_actifs = [n for n in g.nodes() if g.node[n]['type'] == 'M']

    # Number of active metabolites
    nb_metabo_actifs = len(liste_metabo_actifs)

    # For instance, find the number of reactants of the reaction : 'R_PIPLC_HDE_PALM_c'
    g.in_degree('R_PIPLC_HDE_PALM_c')

    # Find the number of products of the reaction : 'R_PIPLC_HDE_PALM_c'
    g.out_degree('R_PIPLC_HDE_PALM_c')

    # The number of reactions generating the metabolite 'M_nh4_c'
    g.in_degree('M_nh4_c')

    # Number of reactions using 'M_nh4_c'
    g.out_degree('M_nh4_c')

    # A list of all the reactions generating each of the metabolites
    in_degree_metabo = []
    for n in liste_metabo_actifs:
        in_degree_metabo.append( g.out_degree(n) )
    #  "in_degree_metabo" is a vector of in-degree of each active metabolite. We can extract different statistics from this vector.
    in_degree_metabo_dist = collections.Counter( sorted(in_degree_metabo, reverse=True) )
    # see: https://networkx.github.io/documentation/stable/auto_examples/drawing/plot_degree_histogram.html
    # Number of metabolites generated by 3 reactions
    in_degree_metabo_dist[3]
    # Mean number of in-degree
    np.mean(in_degree_metabo)
    # Standard deviation of in-degree
    np.std(in_degree_metabo)
    # Median of in-degree
    np.median(in_degree_metabo)
    # 90th percentile of in-degree
    np.percentile(in_degree_metabo, 90)
    # In-degree variance
    np.var(in_degree_metabo)

    # Example: We can shorten the last example using 'list comprehension'
    in_degree_metabo_dist = collections.Counter( sorted([g.out_degree(n) for n in liste_metabo_actifs], reverse=True) )

    # Calculate "betweenness centrality" of nodes. This calculates the number of 'Shortest paths' through each node.
    # see: https://networkx.github.io/documentation/stable/reference/algorithms/generated/networkx.algorithms.centrality.betweenness_centrality.html#networkx.algorithms.centrality.betweenness_centrality
    bc = nx.algorithms.centrality.betweenness_centrality(g, normalized=False)
    # The number of shortest paths through the reaction "R_BIOMASS_CARB_c" is
    # bc['R_BIOMASS_CARB_c']

    # Closeness centrality
    # see: https://networkx.github.io/documentation/stable/reference/algorithms/generated/networkx.algorithms.centrality.closeness_centrality.html#networkx.algorithms.centrality.closeness_centrality
    cc = nx.algorithms.centrality.closeness_centrality(g)
    # "closeness centrality" of the node "R_DMPSS_c" is
    # bc['R_DMPSS_c']

    # Testing if the network is "weakly connected", which is expected.
    # This could be studied in the network with a deleted reaction.
    # see: https://networkx.github.io/documentation/stable/reference/algorithms/generated/networkx.algorithms.components.is_weakly_connected.html#networkx.algorithms.components.is_weakly_connected
    nx.algorithms.components.is_weakly_connected(g)
    
    nx.algorithms.components.number_weakly_connected_components(g)
    
    # Testing if the network is "strongly connected"
    # see: https://networkx.github.io/documentation/stable/reference/algorithms/generated/networkx.algorithms.components.is_strongly_connected.html#networkx.algorithms.components.is_strongly_connected
    nx.algorithms.components.is_strongly_connected(g)
    # Number of components that are "strongly connected"
    # see: https://networkx.github.io/documentation/stable/reference/algorithms/generated/networkx.algorithms.components.number_strongly_connected_components.html#networkx.algorithms.components.number_strongly_connected_components
    nx.algorithms.components.number_strongly_connected_components(g)
    
# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
