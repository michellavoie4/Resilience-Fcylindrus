# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
#
# Python 3 code with UTF-8 as default encoding
#
# Functions to be used with the main script called 'main_analyse_reseau'
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


    # Liste de tous les noeuds.
    liste_noeuds = g.nodes()

    # Nombre total de noeuds (reactions + metabolites).
    nb_noeuds = g.number_of_nodes()

    # Nombre de liens.
    nb_liens = g.number_of_edges()

    # Liste des reactions actives (reactions pour lesquelles le flux n'est pas nul).
    liste_react_actives = [n for n in g.nodes() if g.node[n]['type'] == 'R']

    # Nombre de reactions actives.
    nb_react_actives = len(liste_react_actives)

    # Liste des reactions actives (le flux d'au moins une des reactions dans lesquelles le
    #  metabolite apparait n'est pas nul).
    liste_metabo_actifs = [n for n in g.nodes() if g.node[n]['type'] == 'M']

    # Nombre de reactions actives.
    nb_metabo_actifs = len(liste_metabo_actifs)

    # Pour iterer sur tous les noeuds (n correspond au nom de chaque noeud et sera utilise pour
    #   aller chercher l'information voulue).
    for n in liste_noeuds:
        1# [remplacer le 1 par le code]

    # Pour iterer sur toutes les reactions actives.
    for n in liste_react_actives:
        1# [remplacer le 1 par le code]

    # Pour iterer sur toutes les reactions actives.
    for n in liste_metabo_actifs:
        1# [remplacer le 1 par le code]

    # Le nombre de reactifs de la reaction 'R_PIPLC_HDE_PALM_c'
    g.in_degree('R_PIPLC_HDE_PALM_c')

    # Le nombre de produits de la reaction 'R_PIPLC_HDE_PALM_c'
    g.out_degree('R_PIPLC_HDE_PALM_c')

    # Le nombre de reactions pouvant generer le metabolite 'M_nh4_c'
    g.in_degree('M_nh4_c')

    # Le nombre de reactions auxquel le metabolite 'M_nh4_c' participe
    g.out_degree('M_nh4_c')

    # Exemple: lister le nombre de reactions pouvant generer chacun des metabolites.
    in_degree_metabo = []
    for n in liste_metabo_actifs:
        in_degree_metabo.append( g.out_degree(n) )
    # L'objet "in_degree_metabo" est un vecteur contenant la valeur des in-degrees de
    #   chacun des metabolites actifs duquel il est possible d'extraire differentes
    #   quantites statistiques.
    in_degree_metabo_dist = collections.Counter( sorted(in_degree_metabo, reverse=True) )
    # voir: https://networkx.github.io/documentation/stable/auto_examples/drawing/plot_degree_histogram.html
    # Le nombre de metabolites etant generes par 3 reactions est donne par
    in_degree_metabo_dist[3]
    # Nombre moyen des in-degrees
    np.mean(in_degree_metabo)
    # Ecart-type des in-degrees
    np.std(in_degree_metabo)
    # Valeur mediane des in-degrees
    np.median(in_degree_metabo)
    # 90e percentile des in-degrees
    np.percentile(in_degree_metabo, 90)
    # Variance des in-degrees
    np.var(in_degree_metabo)

    # Exemple: il est possible de condenser l'exemple precedent en une seule ligne via
    #   le concept de "list comprehensions"
    in_degree_metabo_dist = collections.Counter( sorted([g.out_degree(n) for n in liste_metabo_actifs], reverse=True) )

    # Calcul la "betweenness centrality" des noeuds. Autrement dit, compte de
    #  nombre de "shortest paths" passant par chacun des noeuds.
    # voir: https://networkx.github.io/documentation/stable/reference/algorithms/generated/networkx.algorithms.centrality.betweenness_centrality.html#networkx.algorithms.centrality.betweenness_centrality
    bc = nx.algorithms.centrality.betweenness_centrality(g, normalized=False)
    # Le nombre de chemins passant par la reaction "R_BIOMASS_CARB_c" est
    # bc['R_BIOMASS_CARB_c']

    # Closeness centrality
    # voir: https://networkx.github.io/documentation/stable/reference/algorithms/generated/networkx.algorithms.centrality.closeness_centrality.html#networkx.algorithms.centrality.closeness_centrality
    cc = nx.algorithms.centrality.closeness_centrality(g)
    # La "closeness centrality" du noeud "R_DMPSS_c" est
    # bc['R_DMPSS_c']

    # Verification que le reseau est "weakly connected" ce qui devrait etre le cas.
    # De fait, ceci est une propriete que nous pourrions mesurer sur les reseaux pour
    #  lesquels une reaction a ete supprimee avant l'optimisation du modele.
    # voir: https://networkx.github.io/documentation/stable/reference/algorithms/generated/networkx.algorithms.components.is_weakly_connected.html#networkx.algorithms.components.is_weakly_connected
    nx.algorithms.components.is_weakly_connected(g)
    
    nx.algorithms.components.number_weakly_connected_components(g)
    
    # Verification si le reseau est "strongly connected"
    # voir: https://networkx.github.io/documentation/stable/reference/algorithms/generated/networkx.algorithms.components.is_strongly_connected.html#networkx.algorithms.components.is_strongly_connected
    nx.algorithms.components.is_strongly_connected(g)
    # Nombre de composantes "strongly connected"
    # voir: https://networkx.github.io/documentation/stable/reference/algorithms/generated/networkx.algorithms.components.number_strongly_connected_components.html#networkx.algorithms.components.number_strongly_connected_components
    nx.algorithms.components.number_strongly_connected_components(g)
    
# =~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=~=
  
    








