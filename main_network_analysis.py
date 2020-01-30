#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 23:10:35 2019

@author: mlavoie

Written in Python 3
"""

# Script creating a bipartite network of F. cylindrus metabolism and calculating different properties

def clear_all():
    """Clears all the variables from the workspace of the spyder application."""
    gl = globals().copy()
    for var in gl:
        if var[0] == '_': continue
        if 'func' in str(globals()[var]): continue
        if 'module' in str(globals()[var]): continue

        del globals()[var]
if __name__ == "__main__":
    clear_all()
    # insert here your code
    
import os
os.chdir('/Users/mlavoie/Documents/MATLAB/Paper_resilience')
cwd = os.getcwd()
print(cwd)
# Need to add 'network_metrics.py" and "Table_bipartite.txt" in the working directory
# "Table_bipartite.txt" is produced by the file 'bipartite.m' in MATLAB

########
# run network_metrics.py
# This runs a function and executes all the script as it was a main script
# exec(open("network_metrics.py").read())
# g
######

# Rather than execute directly different main script, I prefer to import the file as a module
import network_metrics

g = network_metrics.load_graph("Table_bipartite.txt")

# g['M_biomass_c'] # gives its neighbour and the associated weight
# g['M_tag1819Z160160_c']

# g.nodes() # Print all nodes of the network (metabolite and reaction names)
len(g.nodes()) # total number of nodes in the network


# If we want to perform the calculation without biomass components
# g = load_graph("Table_bipartite_bio.txt")

g.in_degree('R_PIPLC_HDE_PALM_c') # 2
g.out_degree('R_PIPLC_HDE_PALM_c') # 3
g.degree # gives total degrre for each node (reaction or metabolites) 'R_PIPLC_HDE_PALM_c': 5
len(g.degree)

############################################################################
# Draw the network
# nx.draw(g, with_labels=False, node_size=1, width=1, linewidths=0.3, alpha = 0.6, arrows = False) 

# width = Line width of edges (default = 1), font_weight='bold')
# linewidths = Line width of symbol border (default =1.0)
# alpha = The node transparency (default = 1)
# arrows = Truw is default = draw arrow head

# Produce a bipartite graph in .gml
# K_3_5 = nx.complete_bipartite_graph(3, 5)
# nx.write_gml(K_3_5, '/Users/mlavoie/Documents/MATLAB/Paper_resilience/bipartite.gml')

# https://networkx.github.io/documentation/networkx-1.10/reference/generated/networkx.algorithms.bipartite.generators.complete_bipartite_graph.html
# After that, in cytoscape, click 'file', 'import', 'network from file', and import the .gml graph

############################################################################


# x is a graph object
# y is a string : 'R' means only reaction node type. 'M' only metabolite node type and 'all' means all nodes

def stat_array(x, y):
    import numpy as np
    import statistics as stat
    g = x
        
    node_nb = 0
    deg_array = np.array([])    # Initialize an empty array, which will store number of degree per node
    indeg_array = np.array([]) 
    outdeg_array = np.array([]) 

    for n in g.nodes():
        if y == 'all': 
            deg = g.degree(n) # Calculate total degree per node
            deg_array = np.append(deg_array, deg) # create the array storing degree
            indeg = g.in_degree(n)
            indeg_array = np.append(indeg_array, indeg)
            outdeg = g.out_degree(n)
            outdeg_array = np.append(outdeg_array, outdeg)
            node_nb +=1
        elif g.node[n]['type'] == y:
            deg = g.degree(n) # Calculate total degree per node
            deg_array = np.append(deg_array, deg) # create the array storing degree
            indeg = g.in_degree(n)
            indeg_array = np.append(indeg_array, indeg)
            outdeg = g.out_degree(n)
            outdeg_array = np.append(outdeg_array, outdeg)
            node_nb +=1
        else:
            pass
        
    # print the number of nodes
    print("the number of node of type", y, "is", node_nb, "\n")
    
    # Create a matrix of degree calculations (3 rows. 1: total degree, 2: in degree, 3: out degree)
    mat_deg = np.vstack((deg_array, indeg_array, outdeg_array)) 
    mat_names = ['total degree', 'in degree', 'out degree']
    
    # Loop for each row (0,1,2) corresponding to total degree, in degree or out degree
    for i in range(0,3):
        print("--------------------------------")
        print("Statistics for :", mat_names[i])
        print("Mean :", np.mean(mat_deg[i]) ) # Compute the total mean of degree
        print("SD :", stat.stdev(mat_deg[i]) )
        print("Median :", np.median(mat_deg[i]) )
        print("IQR :", np.percentile(mat_deg[i], 75) - np.percentile(mat_deg[i], 25)) # Interquartile range, diference between the 75h and the 35th percentile
        print("Variance :", np.var(mat_deg[i]))
        print("--------------------------------")
    
    # Creating a dataframe with important statistics
    import pandas as pd   
    mean_deg_arr = []
    std_deg_arr = []
    med_deg_arr = []
    IQR_deg_arr = []
    for i in range(0,3):
        mean_deg = np.mean(mat_deg[i])
        std_deg = np.std(mat_deg[i])
        med_deg = np.median(mat_deg[i])
        IQR_deg = np.percentile(mat_deg[i], 75) - np.percentile(mat_deg[i], 25)
        mean_deg_arr = np.append(mean_deg_arr, mean_deg)
        std_deg_arr = np.append(std_deg_arr, std_deg)
        med_deg_arr = np.append(med_deg_arr, med_deg)
        IQR_deg_arr = np.append(IQR_deg_arr, IQR_deg)
    
    stat_deg = np.vstack((mean_deg_arr, std_deg_arr, med_deg_arr, IQR_deg_arr))
    df = pd.DataFrame(stat_deg.transpose())
    df.columns = ['mean_deg_arr', 'std_deg_arr', 'med_deg_arr', 'IQR_deg_arr']
    df.index = mat_names
    
    # Plot histogram of frequency or degree distribution
    import matplotlib.pyplot as plt

    plt.rcParams.update({'figure.figsize':(7,5), 'figure.dpi':100})

    plt.hist(deg_array, bins=50)
    plt.gca().set(title='Frequency Histogram of total degree', ylabel='Frequency');
    plt.show()

    plt.hist(indeg_array, bins=50)
    plt.gca().set(title='Frequency Histogram of in_degree', ylabel='Frequency');
    plt.show()

    plt.hist(outdeg_array, bins=50)
    plt.gca().set(title='Frequency Histogram of out_degree', ylabel='Frequency');
    plt.show()
    
    from scipy import stats
    import numpy as np
    res = stats.relfreq(deg_array, numbins=50)
    res.frequency
    np.sum(res.frequency)  
    x = res.lowerlimit + np.linspace(0, res.binsize*res.frequency.size, res.frequency.size)

    fig = plt.figure(figsize=(5, 4))
    ax = fig.add_subplot(1, 1, 1)
    ax.bar(x, res.frequency, width=res.binsize)
    ax.set_title('Relative frequency histogram of total degree')
    #ax.set_xlim([x.min(), x.max()])
    # ax.set_xlim([0, 5.0])
    plt.show()
    
    res = stats.relfreq(indeg_array, numbins=50)
    res.frequency
    np.sum(res.frequency)  
    x = res.lowerlimit + np.linspace(0, res.binsize*res.frequency.size, res.frequency.size)

    fig = plt.figure(figsize=(5, 4))
    ax = fig.add_subplot(1, 1, 1)
    ax.bar(x, res.frequency, width=res.binsize)
    ax.set_title('Relative frequency histogram of total degree')
    #ax.set_xlim([x.min(), x.max()])
    # ax.set_xlim([0, 5.0])
    plt.show()
    
    res = stats.relfreq(outdeg_array, numbins=50)
    res.frequency
    np.sum(res.frequency)  
    x = res.lowerlimit + np.linspace(0, res.binsize*res.frequency.size, res.frequency.size)

    fig = plt.figure(figsize=(5, 4))
    ax = fig.add_subplot(1, 1, 1)
    ax.bar(x, res.frequency, width=res.binsize)
    ax.set_title('Relative frequency histogram of total degree')
    #ax.set_xlim([x.min(), x.max()])
    # ax.set_xlim([0, 5.0])
    plt.show()
    
    return(df)

##############################################################################
# Computing total degree, in degree, out degree statistics for different subsets of the network
df_metab = stat_array(g, 'M') # For metabolite nodes only
stat_array(g, 'R') # For reaction nodes only
df_all = stat_array(g, 'all') # For both node types (i.e., all nodes)
#############################################################################

##############################################################################
# Calculate the degree numbers (and statistics about that) for a subset of METABOLITES
# INPUT: node_List_csv is a .csv file containing a subset of sensitive or robust metabolites
# OUTPUT: mat_deg : a matrix of all total degree (row 1), in degree (row 2) and out degree (row 3)
#         subset_mets: a numpy array showing all metabolite in the subsets
#         df : a dataframe 3 rows (total degree, in degree, out degree) and 4 columns (mean, std, median, IQR) summarizing the results
def degree_subset(node_List_csv):
    # Read a .csv file containing a subset of robust reactions
    import pandas as pd
    import numpy as np
    import statistics as stat
    df = pd.read_csv(node_List_csv) # One column only , no ',' adn no header    # df=pd.read_csv('Rxns_Robust.csv', sep=',',header=None)
    df.values # vector of robust reactions

    arr_nodes = np.array(g.nodes()) # convert g.nodes to a numpy array

    # Convert the array 'df.values' to 'robust_mets'
    # Necassary because the functions g.degree does not work on 'df.values'
    node_nb = 0
    subset_mets = np.array([])
    for n in arr_nodes:
        if np.any(df.values == n): # test wether or not any values in 'df.values' is equal to 'n'
            #deg = g.degree(n) # Calculate total degree per node
            subset_mets = np.append(subset_mets, n) # create the array 
            node_nb += 1
        
    deg_array = np.array([])
    indeg_array = np.array([])
    outdeg_array = np.array([])
    for n in subset_mets: # does not work with 'df.values:'
         deg_array = np.append(deg_array, g.degree(n))
         indeg_array = np.append(indeg_array, g.in_degree(n))
         outdeg_array = np.append(outdeg_array, g.out_degree(n))   
      
        # print the number of nodes
    print("the number of node is", node_nb, "\n")
    
    # Create a matrix of degree calculations (3 rows. 1: total degree, 2: in degree, 3: out degree)
    mat_deg = np.vstack((deg_array, indeg_array, outdeg_array)) 
    mat_names = ['total degree', 'in degree', 'out degree']
    
    # Loop for each row (0,1,2) corresponding to total degree, in degree or out degree
    for i in range(0,3):
        print("--------------------------------")
        print("Statistics for :", mat_names[i])
        print("Mean :", np.mean(mat_deg[i]) ) # Compute the total mean of degree
        print("SD :", stat.stdev(mat_deg[i]) )
        print("Median :", np.median(mat_deg[i]) )
        print("IQR :", np.percentile(mat_deg[i], 75) - np.percentile(mat_deg[i], 25)) # Interquartile range, diference between the 75h and the 35th percentile
        print("Variance :", np.var(mat_deg[i]))
        
        ind = np.where(mat_deg[i] > 20)
        print("Metabolites with highest > 20", mat_names[i], " : ", subset_mets[ind])
        print(mat_names[i], " : ", mat_deg[i][ind])
        
        ind = np.where(mat_deg[i] < 2.5)
        print("Metabolites with lowest < 2.5", mat_names[i], " : ", subset_mets[ind])
        print(mat_names[i], " : ", mat_deg[i][ind])
        print("Number of metabolites with lowest", mat_names[i], " : ", len(subset_mets[ind]))
        print("Number of metabolites with highest", mat_names[i], " : ", len(subset_mets[ind]))
        print("--------------------------------")
    
    # Creating a dataframe with important statistics
    import pandas as pd   
    mean_deg_arr = []
    std_deg_arr = []
    med_deg_arr = []
    IQR_deg_arr = []
    for i in range(0,3):
        mean_deg = np.mean(mat_deg[i])
        std_deg = np.std(mat_deg[i])
        med_deg = np.median(mat_deg[i])
        IQR_deg = np.percentile(mat_deg[i], 75) - np.percentile(mat_deg[i], 25)
        mean_deg_arr = np.append(mean_deg_arr, mean_deg)
        std_deg_arr = np.append(std_deg_arr, std_deg)
        med_deg_arr = np.append(med_deg_arr, med_deg)
        IQR_deg_arr = np.append(IQR_deg_arr, IQR_deg)
    
    stat_deg = np.vstack((mean_deg_arr, std_deg_arr, med_deg_arr, IQR_deg_arr))
    df = pd.DataFrame(stat_deg.transpose())
    df.columns = ['mean_deg_arr', 'std_deg_arr', 'med_deg_arr', 'IQR_deg_arr']
    df.index = mat_names
        
    return(mat_deg, subset_mets, df)



# Calculate the degree numbers for a subset of robust METABOLITES and return a matrix of degree numbers
# This also returns a vector of metabolites.
mat_rob, robust_mets, df_rob = degree_subset('Mets_Robust.csv')


# Calculate the degree numbers for a subset of sensitives and robust METABOLITES
mat_sens, sens_mets, df_sens = degree_subset('Mets_Less_Robust.csv')


##############################################################################
# function plotting the relative frequency of total, in or out degree number  for sensitive or robust metabolites
def histoMetDegree(mat_rob, mat_sens):
    import matplotlib.pyplot as plt
    from scipy import stats
    import numpy as np
        
    # Plot relative frequency of total degree: Robust or sensitive metabolites
    for i in range(0,3):
        plt.subplot(3,1,i+1)
        numbins=int(max(mat_sens[i])) # Same than numbins=int(max(indeg_array_rob))
        res = stats.relfreq(mat_sens[i], numbins=numbins, defaultreallimits=(0,numbins))
        res2 = stats.relfreq(mat_rob[i], numbins=numbins, defaultreallimits=(0,numbins))
        # res.frequency
        # np.sum(res.frequency)  
        x = res.lowerlimit + np.linspace(0, res.binsize*res.frequency.size, res.frequency.size)
        x2 = res2.lowerlimit + np.linspace(0, res2.binsize*res2.frequency.size, res2.frequency.size)
    
        #plt.rcParams.update({'figure.figsize':(7,5), 'figure.dpi':100})
        p1 = plt.bar(x, res.frequency, width=res.binsize, alpha = 0.4) # bar width should be decreased otherwise we cannot see properly both bars for each frequency
        p2 = plt.bar(x2, res2.frequency, width=res2.binsize, alpha=0.4)
        ylab = ['Total degree', 'In degree', 'Out degree']
        #plt.title(titles[i])
        plt.ylabel('Relative frequency') #, xlabel='Total degree');
        plt.xlabel(ylab[i])
        plt.legend((p1[0], p2[0]), ('Sensitive metabolites', 'Robust metabolites'))
        plt.show() # Do not add plt.show() before plt.savefig(), otherwise it won,t save
        #plt.savefig('total_degree_rel1.png')

##############################################################################

# Plot the relative frequency of total, in, or out degree for the sensitive metabolites or the robust metabolites (Figure 4 , A, B, and C)
histoMetDegree(mat_rob, mat_sens) 



###################################################################
# Betweeness centrality of METABOLITES

# A function calculating betweeness centrality (bc) or closeness centrality (cc) of Metabolites
# g is a graph object
# central is a string : 'bc' or 'cc' for betweeness or closeness centrality calculations
def centralMetabol(g, central):
    import numpy as np
    import matplotlib.pyplot as plt
    import networkx as nx
    
    mat_rob, robust_mets, df_rob = degree_subset('Mets_Robust.csv')
    mat_sens, sens_mets, df_sens = degree_subset('Mets_Less_Robust.csv')
    
    if central == 'bc':
        bc = nx.algorithms.centrality.betweenness_centrality(g, normalized=False)
        cent = bc
    elif central == 'cc':
        cc = nx.algorithms.centrality.closeness_centrality(g)
        cent = cc
    
    c_array = np.array([]) # array of bc for all nodes
    c_array_M = np.array([]) # array of bc for metabolites only
    c_array_rob = np.array([]) # array of cc for robust reactions only
    c_array_sens = np.array([]) # array of cc for sensitive reactions only
    for n, val in cent.items():  # iterate over item names and '.values' in the dictionary
        c_array = np.append(c_array, val) # Create bc_array
        if n.startswith('M_') == True: 
            c_array_M = np.append(c_array_M, val)        
        if np.any(robust_mets == n): 
            c_array_rob = np.append(c_array_rob, val)
        if np.any(sens_mets == n): 
            c_array_sens = np.append(c_array_sens, val)  
    
    # Create a matrix of degree calculations (3 rows. 1: total degree, 2: in degree, 3: out degree)
    # Since arrays of different length cannot be stacked, array must be padded with nan
    maxL = len(cent.items())
    
    c_array_M = np.pad(c_array_M, (0, maxL - len(c_array_M)), 'constant', constant_values=np.nan)
    assert(len(c_array_M == maxL))
    c_array_rob = np.pad(c_array_rob, (0, maxL - len(c_array_rob)), 'constant', constant_values=np.nan)
    assert(len(c_array_rob == maxL))
    c_array_sens = np.pad(c_array_sens, (0, maxL - len(c_array_sens)), 'constant', constant_values=np.nan)
    assert(len(c_array_sens == maxL))
    
    mat_cent = np.vstack((c_array_M, c_array_rob, c_array_sens))
    mat_names = ['all metabolites', 'robust metabolites', 'sensitive metabolites']
    
    # Loop for each row (0,1,2) corresponding to total degree, in degree or out degree
    for i in range(0,3):
        print("--------------------------------")
        print("Statistics for :", mat_names[i])
        print("Mean :", np.nanmean(mat_cent[i]) ) # Compute the total mean of degree
        print("SD :", np.nanstd(mat_cent[i]) )
        print("Median :", np.nanmedian(mat_cent[i]) )
        print("IQR :", np.nanpercentile(mat_cent[i], 75) - np.nanpercentile(mat_cent[i], 25)) # Interquartile range, diference between the 75h and the 35th percentile
        print("Variance :", np.nanvar(mat_cent[i]))
        
        '''
        ind = np.where(mat_deg[i] > 20)
        print("Metabolites with highest > 20", mat_names[i], " : ", subset_mets[ind])
        print(mat_names[i], " : ", mat_deg[i][ind])
        
        ind = np.where(mat_deg[i] < 2.5)
        print("Metabolites with lowest < 2.5", mat_names[i], " : ", subset_mets[ind])
        print(mat_names[i], " : ", mat_deg[i][ind])
        print("Number of metabolites with lowest", mat_names[i], " : ", len(subset_mets[ind]))
        print("Number of metabolites with highest", mat_names[i], " : ", len(subset_mets[ind]))
        '''
        print("--------------------------------")
    
    # Creating a dataframe with important statistics
    import pandas as pd   
    mean_cent_arr = []
    std_cent_arr = []
    med_cent_arr = []
    IQR_cent_arr = []
    for i in range(0,3):
        mean_c = np.nanmean(mat_cent[i])
        std_c = np.nanstd(mat_cent[i])
        med_c = np.nanmedian(mat_cent[i])
        IQR_c = np.nanpercentile(mat_cent[i], 75) - np.nanpercentile(mat_cent[i], 25)
        mean_cent_arr = np.append(mean_cent_arr, mean_c)
        std_cent_arr = np.append(std_cent_arr, std_c)
        med_cent_arr = np.append(med_cent_arr, med_c)
        IQR_cent_arr = np.append(IQR_cent_arr, IQR_c)
    
    stat_cent = np.vstack((mean_cent_arr, std_cent_arr, med_cent_arr, IQR_cent_arr))
    df = pd.DataFrame(stat_cent.transpose())
    df.columns = ['mean_cent_arr', 'std_cent_arr', 'med_cent_arr', 'IQR_cent_arr']
    df.index = mat_names
    
    return(mat_cent, df)

# calculate statistics of betweeness centrality
mat_bc, df_bc = centralMetabol(g, 'bc')

# calculate statistics of closeness centrality
mat_cc, df_cc = centralMetabol(g, 'cc')

# remove nan in the 2D array of bc or cc
#import numpy as np
#[row[~np.isnan(row)] for row in mat_cc]
#[row[~np.isnan(row)] for row in mat_bc]
# DOES NOT WORK: mat_bc = mat_bc[~np.isnan(mat_bc)]

#################################################################
# Creating a dataframe summarizing the main degree statistics, bc and cc
# This is some of the results shown in the Supporting information section S2 (Table SI.7)

import pandas as pd
import numpy as np
df_common = pd.concat([df_all, df_metab, df_rob, df_sens], axis=1)

'''
# Append a column with some notes describing the table
df_note = pd.DataFrame(['df_all', 'df_metab', 'df_rob', 'df_sens'])


x = []
for i in range(4):
    x.append('df_all')
for i in range(4):
    x.append('df_metab')
for i in range(4):
    x.append('df_rob')
for i in range(4):
    x.append('df_sens')

df2 = pd.DataFrame(np.array(x))
df2 = df2.T
df2.columns = df_common.columns # This is important, otherwise append does not work
df_common.append(df2) #, ignore_index=True)
'''

# Adding a note to the dataframe
df_common = df_common.assign(notes=np.array(['this table shows 4 successive dataframe, 1) for all nodes, 2) active metabolites, 3) robust metabolites, 4) for sensitive metabolites' ,np.nan,np.nan]))
    
# Export dataframe to csv
wd = os.getcwd()
export_csv = df_common.to_csv (wd+'/degree_stat.csv', index = True, header=True)


# Export also betweeness centrality and closeness centrality results to csv
wd = os.getcwd()
export_csv = df_bc.to_csv (wd+'/bc_stat.csv', index = True, header=True)
export_csv = df_cc.to_csv (wd+'/cc_stat.csv', index = True, header=True)
#############################################################################


#############################################################################
# Plot relative frequency of betweenness centrality: Robust or sensitive metabolites'
# Here are the results of Figure SI.1 in the SI section of the paper.

import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
res = stats.relfreq(mat_bc[1], numbins=100, defaultreallimits=(0,100000))
res2 = stats.relfreq(mat_bc[2], numbins=100, defaultreallimits=(0,100000))
# res.frequency
# np.sum(res.frequency)  
x = res.lowerlimit + np.linspace(0, res.binsize*res.frequency.size, res.frequency.size)
x2 = res2.lowerlimit + np.linspace(0, res2.binsize*res2.frequency.size, res2.frequency.size)
# Here we must also adjust "res2.lowerlimit" so that the bars do not overlay.

plt.rcParams.update({'figure.figsize':(7,5), 'figure.dpi':100})
p1 = plt.bar(x, res.frequency, width=res.binsize, alpha=0.4)
p2 = plt.bar(x2, res2.frequency, width=res2.binsize, alpha=0.4)
plt.title('Relative frequency histogram of betweenness centrality')
plt.ylabel('Relative frequency')
plt.xlabel('Betweenness centrality')
#plt.xlim(0,50000)
plt.legend((p1[0], p2[0]), ('Sensitive metabolites', 'Robust metabolites'))
plt.show() # Do not add plt.show() before plt.savefig(), otherwise it won,t save
#plt.savefig('Betweenness_central_rel1.png')

# Plot relative frequency of closeness centrality: Robust or sensitive metabolites'
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
res = stats.relfreq(mat_cc[1], numbins=50, defaultreallimits=(0,0.25))
res2 = stats.relfreq(mat_cc[2], numbins=50, defaultreallimits=(0,0.25))
# res.frequency
# np.sum(res.frequency)  
x = res.lowerlimit + np.linspace(0, res.binsize*res.frequency.size, res.frequency.size)
x2 = res2.lowerlimit + np.linspace(0, res2.binsize*res2.frequency.size, res2.frequency.size)
# Here we must also adjust "res2.lowerlimit" so that the bars do not overlay.

plt.rcParams.update({'figure.figsize':(7,5), 'figure.dpi':100})
p1 = plt.bar(x, res.frequency, width=res.binsize, alpha=0.4)
p2 = plt.bar(x2, res2.frequency, width=res2.binsize, alpha=0.4)
plt.title('Relative frequency histogram of closeness centrality')
plt.ylabel('Relative frequency')
plt.xlabel('Closeness centrality')
#plt.xlim(0,50000)
plt.legend((p1[0], p2[0]), ('Sensitive metabolites', 'Robust metabolites'))
plt.show() # Do not add plt.show() before plt.savefig(), otherwise it won,t save
#plt.savefig('Closeness_central_rel1.png')
 






'''
############################################################################
# Calculate the degree distribution of all active nodes in the network
import numpy as np
import statistics as stat
deg_array = np.array([])    # Initialize an empty array, which will store number of degree per node
indeg_array = np.array([]) 
outdeg_array = np.array([]) 
for n in g.nodes():
    deg = g.degree(n) # Calculate total degree per node
    deg_array = np.append(deg_array, deg) # create the array storing degree
    indeg = g.in_degree(n)
    indeg_array = np.append(indeg_array, indeg)
    outdeg = g.out_degree(n)
    outdeg_array = np.append(outdeg_array, outdeg)
'''

'''
print(np.mean(deg_array)) # Compute the total mean of degree
print(stat.stdev(deg_array))
print(np.median(deg_array))
print(np.percentile(deg_array, 75) - np.percentile(deg_array, 25)) # Interquartile range, diference between the 75h and the 35th percentile
print(np.var(deg_array))

print(np.mean(indeg_array)) # Compute the total mean of in_degree
print(stat.stdev(indeg_array))
print(np.median(indeg_array))
print(np.percentile(indeg_array, 75) - np.percentile(indeg_array, 25))

print(np.mean(outdeg_array)) # Compute the total mean of out_degree
print(stat.stdev(outdeg_array))
print(np.median(outdeg_array))
print(np.percentile(outdeg_array, 75) - np.percentile(outdeg_array, 25))
'''

'''
# Create a matrix of degree calculations (3 rows. 1: total degree, 2: in degree, 3: out degree)

mat_deg = np.vstack((deg_array, indeg_array, outdeg_array)) 
# inv = mat_deg.transpose() # if we want 3 columns instead of 3 rows.

mat_deg.shape
# Select elements of the multidimentional array
mat_deg[[1], [1]]
# Select the second row
mat_deg[[1]]
len(mat_deg[1])
len(mat_deg[0])    
'''


'''
# Plot Histogram of frequency or degree distribution
import matplotlib.pyplot as plt
#%matplotlib inline
plt.rcParams.update({'figure.figsize':(7,5), 'figure.dpi':100})

plt.hist(deg_array, bins=50)
plt.gca().set(title='Frequency Histogram of total degree', ylabel='Frequency');
plt.show()


plt.hist(indeg_array, bins=50)
plt.gca().set(title='Frequency Histogram of in_degree', ylabel='Frequency');
plt.show()

plt.hist(outdeg_array, bins=50)
plt.gca().set(title='Frequency Histogram of out_degree', ylabel='Frequency');
plt.show()

# Plot relative frequency of degree distribution
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
res = stats.relfreq(deg_array, numbins=50)
res.frequency
np.sum(res.frequency)  
x = res.lowerlimit + np.linspace(0, res.binsize*res.frequency.size, res.frequency.size)

fig = plt.figure(figsize=(5, 4))
ax = fig.add_subplot(1, 1, 1)
ax.bar(x, res.frequency, width=res.binsize)
ax.set_title('Relative frequency histogram of total degree')
#ax.set_xlim([x.min(), x.max()])
# ax.set_xlim([0, 5.0])
plt.show()


#######################################################

#########################################################
# Calculate the degree distribution of reaction nodes only (excluding metabolites)
import numpy as np
import statistics as stat
nb_lignes = 0
deg_array = np.array([])    # Initial an empty array, which will store number of degree per node
indeg_array = np.array([]) 
outdeg_array = np.array([])
for n in g.nodes():  
    if g.node[n]['type'] == 'R': # OR if n.startswith('R_') == True:
        #print(n)
        deg = g.degree(n) # Calculate total degree per node
        deg_array = np.append(deg_array, deg) # create the array storing degree
        indeg = g.in_degree(n)
        indeg_array = np.append(indeg_array, indeg)
        outdeg = g.out_degree(n)
        outdeg_array = np.append(outdeg_array, outdeg)
        nb_lignes += 1 # Counting the number of reactions
print(nb_lignes)     # total number of reactions

print(np.mean(deg_array)) # Compute the total mean of degree
print(stat.stdev(deg_array))
print(np.median(deg_array))
print(np.percentile(deg_array, 75) - np.percentile(deg_array, 25))
print(np.percentile(deg_array, 90) )
print(np.percentile(deg_array, 99) )

print(np.mean(indeg_array)) # Compute the total mean of in_degree
print(stat.stdev(indeg_array))
print(np.median(indeg_array))
print(np.percentile(indeg_array, 75) - np.percentile(indeg_array, 25))
print(np.percentile(indeg_array, 90))

print(np.mean(outdeg_array)) # Compute the total mean of out_degree
print(stat.stdev(outdeg_array))
print(np.median(outdeg_array))
print(np.percentile(outdeg_array, 75) - np.percentile(outdeg_array, 25))
print(np.percentile(outdeg_array, 90))

# Plot Histogram of frequency or degree distribution
import matplotlib.pyplot as plt
#%matplotlib inline
plt.rcParams.update({'figure.figsize':(7,5), 'figure.dpi':100})

plt.hist(deg_array, bins=50)
plt.gca().set(title='Frequency Histogram of total degree', ylabel='Frequency');

plt.hist(indeg_array, bins=50)
plt.gca().set(title='Frequency Histogram of in_degree', ylabel='Frequency');

plt.hist(outdeg_array, bins=50)
plt.gca().set(title='Frequency Histogram of out_degree', ylabel='Frequency');

# The removing technique does not work
# nb_lignes = 0
# met_list = np.array([])
# for n in g.nodes(): 
#    if ('M_' in n) == True:
#        nb_lignes += 1
#        met = n # print(n)
#        met_list = np.append(met_list, met)
        
# Remove nodes corresponding to metabolites   ## DOES NOT WORK !!     
# g.remove_node(met_list)
# len(g)  
#    
##########################################################################       


##########################################################################
# Calculate the degree distribution of metabolite nodes only (excluding reactions)
import numpy as np
import statistics as stat
nb_lignes = 0
deg_array = np.array([])    # Initial an empty array, which will store number of degree per node
indeg_array = np.array([]) 
outdeg_array = np.array([])
for n in g.nodes():  
    if n.startswith('M_') == True: # if ('M_' in n) == True:
        print(n)
        deg = g.degree(n) # Calculate total degree per node
        deg_array = np.append(deg_array, deg) # create the array storing degree
        indeg = g.in_degree(n)
        indeg_array = np.append(indeg_array, indeg)
        outdeg = g.out_degree(n)
        outdeg_array = np.append(outdeg_array, outdeg)
        nb_lignes += 1 # Counting the number of reactions
print(nb_lignes)     # total number of reactions

print(np.mean(deg_array)) # Compute the total mean of degree
print(stat.stdev(deg_array))
print(np.median(deg_array))
print(np.percentile(deg_array, 75) - np.percentile(deg_array, 25))
print(np.percentile(deg_array, 90))

print(np.mean(indeg_array)) # Compute the total mean of in_degree
print(stat.stdev(indeg_array))
print(np.median(indeg_array))
print(np.percentile(indeg_array, 75) - np.percentile(indeg_array, 25))
print(np.percentile(indeg_array, 90))

print(np.mean(outdeg_array)) # Compute the total mean of in_degree
print(stat.stdev(outdeg_array))
print(np.median(outdeg_array))
print(np.percentile(outdeg_array, 75) - np.percentile(outdeg_array, 25))
print(np.percentile(outdeg_array, 90))

# Plot Histogram of frequency or degree distribution
import matplotlib.pyplot as plt
# %matplotlib inline
plt.rcParams.update({'figure.figsize':(7,5), 'figure.dpi':100})

plt.hist(deg_array, bins=50)
plt.gca().set(title='Frequency Histogram of total degree', ylabel='Frequency');

plt.hist(indeg_array, bins=50)
plt.gca().set(title='Frequency Histogram of in_degree', ylabel='Frequency');

plt.hist(outdeg_array, bins=50)
plt.gca().set(title='Frequency Histogram of out_degree', ylabel='Frequency');
############################################################################
'''



'''
# Select metabolites with highest total degrees in robust_mets or sens_mets
robust_mets
ind = np.where(deg_array_rob > 20)
len(ind[0]) # number of elements in the first dimension of the tuple
print(robust_mets[ind])
print(deg_array_rob[ind])

sens_mets
ind = np.where(deg_array_sens > 20)
len(ind[0]) # number of elements in the first dimension of the tuple
print(sens_mets[ind])
print(deg_array_sens[ind])
# We see no differences between both groups

# Select metabolites with lower total degrees in robust_mets or sens_mets
robust_mets
indice = np.where(deg_array_rob < 2.5)
len(indice[0])
print(robust_mets[indice])
print(deg_array_rob[indice])

sens_mets
indice = np.where(deg_array_sens < 2.5)
len(indice[0])
# print(indice.size)
print(sens_mets[indice])
print(deg_array_sens[indice])
# We see that in the robust Rxns set, there are much less (2 times less?) metabolites with small degree number (2) than the sensitive set of rxns or metabolites

# Select metabolites with highest in degrees in robust_mets or sens_mets
robust_mets
ind = np.where(indeg_array_rob > 20)
len(ind[0])
print(robust_mets[ind])
print(indeg_array_rob[ind])

sens_mets
ind = np.where(indeg_array_sens > 20)
len(ind[0])
print(sens_mets[ind])
print(indeg_array_sens[ind])
# We see no differences

# Select metabolites with lower IN degrees in robust_mets or sens_mets
robust_mets
indice = np.where(indeg_array_rob < 2.5)
print(robust_mets[indice])
print(indeg_array_rob[indice])

sens_mets
indice = np.where(indeg_array_sens < 2.5)
# print(indice.size)
print(sens_mets[indice])
print(indeg_array_sens[indice])
# We see much less in degree equals 1 for the robust set of metabolites.

# Select metabolites with highest OUT degrees in robust_mets or sens_mets
robust_mets
ind = np.where(outdeg_array_rob > 20)
print(robust_mets[ind])
print(outdeg_array_rob[ind])

sens_mets
ind = np.where(outdeg_array_sens > 20)
print(sens_mets[ind])
print(outdeg_array_sens[ind])
# We see no differences

# Select metabolites with lower OUT degrees in robust_mets or sens_mets
robust_mets
indice = np.where(outdeg_array_rob < 2.5)
print(robust_mets[indice])
print(outdeg_array_rob[indice])

sens_mets
indice = np.where(outdeg_array_sens < 2.5)
# print(indice.size)
print(sens_mets[indice])
print(outdeg_array_sens[indice])
# We see much less OUT degree equals 1 for the robust set of metabolites.

############################################################################
'''











'''
###################################################################
# Calculate the degree numbers for a subset of robust METABOLITES

# Read a .csv file containing a subset of robust reactions
import pandas as pd
df = pd.read_csv('Mets_Robust.csv') # One column only , no ',' adn no header    # df=pd.read_csv('Rxns_Robust.csv', sep=',',header=None)
df.values # vector of robust reactions

arr_nodes = np.array(g.nodes()) # convert g.nodes to a numpy array

# Convert the array 'df.values' to 'robust_mets'
# Necassary because the functions g.degree does not work on 'df.values'
nb_lignes = 0
robust_mets = np.array([])
for n in arr_nodes:
    if np.any(df.values == n): # test wether or not any values in 'df.values' is equal to 'n'
        #deg = g.degree(n) # Calculate total degree per node
        robust_mets = np.append(robust_mets, n) # create the array 
        nb_lignes += 1
        
print(nb_lignes)   
print(len(robust_mets))

g.degree(robust_mets) # this works with 'robust_rxn' but not with 'df.values' !!
g.in_degree(robust_mets) # this works with 'robust_rxn' but not with 'df.values' !!
g.out_degree(robust_mets) # this works with 'robust_rxn' but not with 'df.values' !!

# Calculate total, in or out degree for robust Metabolites.
deg_array = np.array([])
indeg_array = np.array([])
outdeg_array = np.array([])
for n in robust_mets: # does not work with 'df.values:'
    deg_array = np.append(deg_array, g.degree(n))
    indeg_array = np.append(indeg_array, g.in_degree(n))
    outdeg_array = np.append(outdeg_array, g.out_degree(n))

print(np.mean(deg_array)) # Compute the total mean of degree
print(stat.stdev(deg_array))
print(np.median(deg_array))
print(np.percentile(deg_array, 75) - np.percentile(deg_array, 25))
print(np.percentile(deg_array, 90) )
print(np.percentile(deg_array, 99) )

print(np.mean(indeg_array)) # Compute the total mean of in_degree
print(stat.stdev(indeg_array))
print(np.median(indeg_array))
print(np.percentile(indeg_array, 75) - np.percentile(indeg_array, 25))
print(np.percentile(indeg_array, 90)) 
print(np.percentile(indeg_array, 99) )

print(np.mean(outdeg_array)) # Compute the total mean of in_degree
print(stat.stdev(outdeg_array))
print(np.median(outdeg_array))
print(np.percentile(outdeg_array, 75) - np.percentile(outdeg_array, 25))
print(np.percentile(outdeg_array, 90) )
print(np.percentile(outdeg_array, 99) )
    
# Plot histogram of frequency or degree distribution
import matplotlib.pyplot as plt
#%matplotlib inline
plt.rcParams.update({'figure.figsize':(7,5), 'figure.dpi':100})

plt.hist(deg_array, bins=50)
plt.gca().set(title='Frequency Histogram of total degree', ylabel='Frequency');

plt.hist(indeg_array, bins=50)
plt.gca().set(title='Frequency Histogram of in_degree', ylabel='Frequency');

plt.hist(outdeg_array, bins=50)
plt.gca().set(title='Frequency Histogram of out_degree', ylabel='Frequency');

deg_array_rob = deg_array
indeg_array_rob = indeg_array
outdeg_array_rob = outdeg_array
############################################################################
'''

'''
###################################################################
# Calculate the degree numbers for a subset of sensitive METABOLITES

# Read a .csv file containing a subset of sensitive reactions
import pandas as pd
df = pd.read_csv('Mets_Less_Robust.csv') # One column only , no ',' adn no header    # df=pd.read_csv('Rxns_Robust.csv', sep=',',header=None)
df.values # vector of robust reactions

arr_nodes = np.array(g.nodes()) # convert g.nodes to a numpy array

# Convert the array 'df.values' to 'robust.rxn'
# Necassary because the functions g.degree does not work on 'df.values'
nb_lignes = 0
sens_mets = np.array([])
for n in arr_nodes:
    if np.any(df.values == n): # test wether or not any values in 'df.values' is equal to 'n'
        #deg = g.degree(n) # Calculate total degree per node
        sens_mets = np.append(sens_mets, n) # create the array 
        nb_lignes += 1
        
print(nb_lignes)   
print(len(sens_mets))

g.degree(sens_mets) # this works with 'robust_rxn' but not with 'df.values' !!
g.in_degree(sens_mets) # this works with 'robust_rxn' but not with 'df.values' !!
g.out_degree(sens_mets) # this works with 'robust_rxn' but not with 'df.values' !!

# Calculate total, in or out degree for sensitive METABOLITES.
deg_array = np.array([])
indeg_array = np.array([])
outdeg_array = np.array([])
for n in sens_mets: # does not work with 'df.values:'
    deg_array = np.append(deg_array, g.degree(n))
    indeg_array = np.append(indeg_array, g.in_degree(n))
    outdeg_array = np.append(outdeg_array, g.out_degree(n))

print(np.mean(deg_array)) # Compute the total mean of degree
print(stat.stdev(deg_array))
print(np.median(deg_array))
print(np.percentile(deg_array, 75) - np.percentile(deg_array, 25))
print(np.percentile(deg_array, 90) )
print(np.percentile(deg_array, 99) )

print(np.mean(indeg_array)) # Compute the total mean of in_degree
print(stat.stdev(indeg_array))
print(np.median(indeg_array))
print(np.percentile(indeg_array, 75) - np.percentile(indeg_array, 25))
print(np.percentile(indeg_array, 90) )
print(np.percentile(indeg_array, 99) )

print(np.mean(outdeg_array)) # Compute the total mean of in_degree
print(stat.stdev(outdeg_array))
print(np.median(outdeg_array))
print(np.percentile(outdeg_array, 75) - np.percentile(outdeg_array, 25))
print(np.percentile(outdeg_array, 90) )
print(np.percentile(outdeg_array, 99) )
    
# Plot histogram of frequency or degree distribution
import matplotlib.pyplot as plt
#%matplotlib inline
plt.rcParams.update({'figure.figsize':(7,5), 'figure.dpi':100})

plt.hist(deg_array, bins=50)
plt.gca().set(title='Frequency Histogram of total degree', ylabel='Frequency');

plt.hist(indeg_array, bins=50)
plt.gca().set(title='Frequency Histogram of in_degree', ylabel='Frequency');

plt.hist(outdeg_array, bins=50)
plt.gca().set(title='Frequency Histogram of out_degree', ylabel='Frequency');

deg_array_sens = deg_array
indeg_array_sens = indeg_array
outdeg_array_sens = outdeg_array

# title='Frequency Histogram of total degree: Robust or sensitive metabolites'
plt.hist(deg_array_sens, bins=122, range=[0, 125], alpha=0.4) # max(deg_array_sens) = 122 max(deg_array_rob) = 122 
plt.hist(deg_array_rob, bins=122, range=[0, 125], alpha=0.4)
plt.gca().set(ylabel='Frequency', xlabel='Total degree');
plt.gca().legend(('sensitive metabolites','robust metabolites'))
plt.savefig('total_degree.png')

# title='Frequency Histogram of IN degree: Robust or sensitive metabolites'
plt.hist(indeg_array_sens, bins=58, range=[0, 60], alpha=0.4) # max(indeg_array_sens) = 58 max(indeg_array_rob) = 58
plt.hist(indeg_array_rob, bins=58, range=[0,60], alpha=0.4)
plt.gca().set(ylabel='Frequency', xlabel='In-degree');
plt.gca().legend(('sensitive metabolites','robust metabolites'))
plt.savefig('In_degree.png')

# title='Frequency Histogram of OUT degree: Robust or sensitive metabolites'
plt.hist(outdeg_array_sens, bins=80, range=[0, 80], alpha=0.4) # max(outdeg_array_sens) = 80 max(outdeg_array_rob) = 80
plt.hist(outdeg_array_rob, bins=80, range=[0, 80], alpha=0.4)
plt.gca().set(ylabel='Frequency', xlabel='Out-degree');
plt.gca().legend(('sensitive metabolites','robust metabolites'))
plt.savefig('Out_degree.png')

# Plot relative frequency of total degree: Robust or sensitive metabolites'
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
numbins=int(max(deg_array_sens)) # Same than numbins=int(max(indeg_array_rob))
res = stats.relfreq(deg_array_sens, numbins=numbins, defaultreallimits=(0,numbins))
res2 = stats.relfreq(deg_array_rob, numbins=numbins, defaultreallimits=(0,numbins))
# res.frequency
# np.sum(res.frequency)  
x = res.lowerlimit + np.linspace(0, res.binsize*res.frequency.size, res.frequency.size)
x2 = res2.lowerlimit + np.linspace(0, res2.binsize*res2.frequency.size, res2.frequency.size)

plt.rcParams.update({'figure.figsize':(7,5), 'figure.dpi':100})
p1 = plt.bar(x, res.frequency, width=res.binsize, alpha = 0.4) # bar width should be decreased otherwise we cannot see properly both bars for each frequency
p2 = plt.bar(x2, res2.frequency, width=res2.binsize, alpha=0.4)
plt.title('Relative frequency histogram of total degree')
plt.ylabel('Relative frequency') #, xlabel='Total degree');
plt.xlabel('Total degree')
plt.legend((p1[0], p2[0]), ('Sensitive metabolites', 'Robust metabolites'))
# plt.show() # Do not add plt.show() before plt.savefig(), otherwise it won,t save
plt.savefig('total_degree_rel1.png')

# Plot relative frequency of In-degree: Robust or sensitive metabolites'
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
numbins=int(max(indeg_array_sens)) # Same than numbins=int(max(indeg_array_rob))
res = stats.relfreq(indeg_array_sens, numbins=numbins, defaultreallimits=(0,numbins))
res2 = stats.relfreq(indeg_array_rob, numbins=numbins, defaultreallimits=(0,numbins))
# res.frequency
# np.sum(res.frequency)  
x = res.lowerlimit + np.linspace(0, res.binsize*res.frequency.size, res.frequency.size)
x2 = res2.lowerlimit + np.linspace(0, res2.binsize*res2.frequency.size, res2.frequency.size)

plt.rcParams.update({'figure.figsize':(7,5), 'figure.dpi':100})
p1 = plt.bar(x, res.frequency, width=res.binsize, alpha=0.4)
p2 = plt.bar(x2, res2.frequency, width=res2.binsize, alpha=0.4)
plt.title('Relative frequency histogram of in-degree')
plt.ylabel('Relative frequency')
plt.xlabel('In-degree')
plt.legend((p1[0], p2[0]), ('Sensitive metabolites', 'Robust metabolites'))
# plt.show() # Do not add plt.show() before plt.savefig(), otherwise it won,t save
plt.savefig('In_degree_rel1.png')

# Plot relative frequency of OUT-degree: Robust or sensitive metabolites'
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
numbins=int(max(outdeg_array_sens)) # Same than numbins=int(max(outdeg_array_rob))
res = stats.relfreq(outdeg_array_sens, numbins=numbins, defaultreallimits=(0,numbins))
res2 = stats.relfreq(outdeg_array_rob, numbins=numbins, defaultreallimits=(0,numbins))
# res.frequency
# np.sum(res.frequency)  
x = res.lowerlimit + np.linspace(0, res.binsize*res.frequency.size, res.frequency.size)
x2 = res2.lowerlimit + np.linspace(0, res2.binsize*res2.frequency.size, res2.frequency.size)
# Here we must also adjust "res2.lowerlimit" so that the bars do not overlay.

plt.rcParams.update({'figure.figsize':(7,5), 'figure.dpi':100})
p1 = plt.bar(x, res.frequency, width=res.binsize, alpha=0.4)
p2 = plt.bar(x2, res2.frequency, width=res2.binsize, alpha=0.4)
plt.title('Relative frequency histogram of out-degree')
plt.ylabel('Relative frequency')
plt.xlabel('Out-degree')
#plt.xlim(0,90)
plt.legend((p1[0], p2[0]), ('Sensitive metabolites', 'Robust metabolites'))
# plt.show() # Do not add plt.show() before plt.savefig(), otherwise it won,t save
plt.savefig('Out_degree_rel1.png')

'''



'''
###################################################################
# Calculate the degree numbers for a subset of robust reactions

# Read a .csv file containing a subset of robust reactions
import pandas as pd
df = pd.read_csv('Rxns_Robust.csv') # One column only , no ',' adn no header    # df=pd.read_csv('Rxns_Robust.csv', sep=',',header=None)
df.values # vector of robust reactions

arr_nodes = np.array(g.nodes()) # convert g.nodes to a numpy array

# Convert the array 'df.values' to 'robust.rxn'
# Necassary because the functions g.degree does not work on 'df.values'
nb_lignes = 0
robust_rxn = np.array([])
for n in arr_nodes:
    if np.any(df.values == n): # test wether or not any values in 'df.values' is equal to 'n'
        #deg = g.degree(n) # Calculate total degree per node
        robust_rxn = np.append(robust_rxn, n) # create the array 
        nb_lignes += 1
        
print(nb_lignes)   
print(len(robust_rxn))

g.degree(robust_rxn) # this works with 'robust_rxn' but not with 'df.values' !!
g.in_degree(robust_rxn) # this works with 'robust_rxn' but not with 'df.values' !!
g.out_degree(robust_rxn) # this works with 'robust_rxn' but not with 'df.values' !!

# Calculate total, in or out degree for robust reactions.
deg_array = np.array([])
indeg_array = np.array([])
outdeg_array = np.array([])
for n in robust_rxn: # does not work with 'df.values:'
    deg_array = np.append(deg_array, g.degree(n))
    indeg_array = np.append(indeg_array, g.in_degree(n))
    outdeg_array = np.append(outdeg_array, g.out_degree(n))

print(np.mean(deg_array)) # Compute the total mean of degree
print(stat.stdev(deg_array))
print(np.median(deg_array))
print(np.percentile(deg_array, 75) - np.percentile(deg_array, 25))
print(np.percentile(deg_array, 90) )
print(np.percentile(deg_array, 99) )

print(np.mean(indeg_array)) # Compute the total mean of in_degree
print(stat.stdev(indeg_array))
print(np.median(indeg_array))
print(np.percentile(indeg_array, 75) - np.percentile(indeg_array, 25))
print(np.percentile(indeg_array, 90)) 

print(np.mean(outdeg_array)) # Compute the total mean of in_degree
print(stat.stdev(outdeg_array))
print(np.median(outdeg_array))
print(np.percentile(outdeg_array, 75) - np.percentile(outdeg_array, 25))
print(np.percentile(outdeg_array, 90) )
    
# Plot histogram of frequency or degree distribution
import matplotlib.pyplot as plt
#%matplotlib inline
plt.rcParams.update({'figure.figsize':(7,5), 'figure.dpi':100})

plt.hist(deg_array, bins=50)
plt.gca().set(title='Frequency Histogram of total degree', ylabel='Frequency');

plt.hist(indeg_array, bins=50)
plt.gca().set(title='Frequency Histogram of in_degree', ylabel='Frequency');

plt.hist(outdeg_array, bins=50)
plt.gca().set(title='Frequency Histogram of out_degree', ylabel='Frequency');
############################################################################


###################################################################
# Calculate the degree numbers for a subset of sensitive reactions

# Read a .csv file containing a subset of sensitive reactions
import pandas as pd
df = pd.read_csv('Rxns_Sens.csv') # One column only , no ',' adn no header    # df=pd.read_csv('Rxns_Robust.csv', sep=',',header=None)
df.values # vector of robust reactions

arr_nodes = np.array(g.nodes()) # convert g.nodes to a numpy array

# Convert the array 'df.values' to 'robust.rxn'
# Necassary because the functions g.degree does not work on 'df.values'
nb_lignes = 0
sens_rxn = np.array([])
for n in arr_nodes:
    if np.any(df.values == n): # test wether or not any values in 'df.values' is equal to 'n'
        #deg = g.degree(n) # Calculate total degree per node
        sens_rxn = np.append(sens_rxn, n) # create the array 
        nb_lignes += 1
        
print(nb_lignes)   
print(len(sens_rxn))

g.degree(sens_rxn) # this works with 'robust_rxn' but not with 'df.values' !!
g.in_degree(sens_rxn) # this works with 'robust_rxn' but not with 'df.values' !!
g.out_degree(sens_rxn) # this works with 'robust_rxn' but not with 'df.values' !!

# Calculate total, in or out degree for robust reactions.
deg_array = np.array([])
indeg_array = np.array([])
outdeg_array = np.array([])
for n in sens_rxn: # does not work with 'df.values:'
    deg_array = np.append(deg_array, g.degree(n))
    indeg_array = np.append(indeg_array, g.in_degree(n))
    outdeg_array = np.append(outdeg_array, g.out_degree(n))

print(np.mean(deg_array)) # Compute the total mean of degree
print(stat.stdev(deg_array))
print(np.median(deg_array))
print(np.percentile(deg_array, 75) - np.percentile(deg_array, 25))
print(np.percentile(deg_array, 90) )
print(np.percentile(deg_array, 99) )

print(np.mean(indeg_array)) # Compute the total mean of in_degree
print(stat.stdev(indeg_array))
print(np.median(indeg_array))
print(np.percentile(indeg_array, 75) - np.percentile(indeg_array, 25))
print(np.percentile(indeg_array, 90) )
print(np.percentile(indeg_array, 99) )

print(np.mean(outdeg_array)) # Compute the total mean of in_degree
print(stat.stdev(outdeg_array))
print(np.median(outdeg_array))
print(np.percentile(outdeg_array, 75) - np.percentile(outdeg_array, 25))
print(np.percentile(outdeg_array, 90) )
print(np.percentile(outdeg_array, 99) )
    
# Plot histogram of frequency or degree distribution
import matplotlib.pyplot as plt
#%matplotlib inline
plt.rcParams.update({'figure.figsize':(7,5), 'figure.dpi':100})

plt.hist(deg_array, bins=50)
plt.gca().set(title='Frequency Histogram of total degree', ylabel='Frequency');

plt.hist(indeg_array, bins=50)
plt.gca().set(title='Frequency Histogram of in_degree', ylabel='Frequency');

plt.hist(outdeg_array, bins=50)
plt.gca().set(title='Frequency Histogram of out_degree', ylabel='Frequency');
############################################################################
'''

############################################################################
# Calculate centrality
############################################################################

'''
# Betweeness centrality or bc for REACTIONS
bc = nx.algorithms.centrality.betweenness_centrality(g, normalized=False)
bc_array = np.array([]) # array of bc for all nodes
bc_array_R = np.array([]) # array of bc for reaction only
bc_array_rob = np.array([]) # array of cc for robust reactions only
bc_array_sens = np.array([]) # array of cc for sensitive reactions only
for n, val in bc.items():  # iterate over item names and '.values' in the dictionary
    bc_array = np.append(bc_array, val) # Create bc_array
    if n.startswith('R_') == True: 
        bc_array_R = np.append(bc_array_R, val)        
    if np.any(robust_rxn == n): 
        bc_array_rob = np.append(bc_array_rob, val)
    if np.any(sens_rxn == n): 
        bc_array_sens = np.append(bc_array_sens, val)  
    
print(np.mean(bc_array)) # mean betweeness centrality for all nodes
print(np.std(bc_array))
print(np.median(bc_array))
print(np.percentile(bc_array, 75) - np.percentile(bc_array, 25))

print(np.mean(bc_array_R)) # mean betweeness centrality for reactions
print(np.std(bc_array_R))
print(np.median(bc_array_R))
print(np.percentile(bc_array_R, 75) - np.percentile(bc_array_R, 25))

print(np.mean(bc_array_rob)) # mean betweeness centrality for robust rxns only
print(np.std(bc_array_rob))
print(np.median(bc_array_rob))
print(np.percentile(bc_array_rob, 75) - np.percentile(bc_array_rob, 25))

print(np.mean(bc_array_sens)) # mean betweeness centrality for sensitive rxns only
print(np.std(bc_array_sens))
print(np.median(bc_array_sens))
print(np.percentile(bc_array_sens, 75) - np.percentile(bc_array_sens, 25))

# Closeness centrality or cc for REACTIONS
cc = nx.algorithms.centrality.closeness_centrality(g)
cc_array = np.array([]) # array of bc for all nodes
cc_array_R = np.array([]) # array of bc for reaction only
cc_array_rob = np.array([]) # array of cc for robust reactions only
cc_array_sens = np.array([]) # array of cc for sensitive reactions only
for n, val in cc.items():  # iterate over item names and '.values' in the dictionary
    cc_array = np.append(cc_array, val) # Create bc_array
    if n.startswith('R_') == True: 
        cc_array_R = np.append(cc_array_R, val)        
    if np.any(robust_rxn == n): 
        cc_array_rob = np.append(cc_array_rob, val)
    if np.any(sens_rxn == n): 
        cc_array_sens = np.append(cc_array_sens, val)  
    
print(np.mean(cc_array)) # mean closeness centrality for all nodes
print(np.std(cc_array))
print(np.median(cc_array))
print(np.percentile(cc_array, 75) - np.percentile(cc_array, 25))

print(np.mean(cc_array_R)) # mean closeness centrality for reactions
print(np.std(cc_array_R))
print(np.median(cc_array_R))
print(np.percentile(cc_array_R, 75) - np.percentile(cc_array_R, 25))

print(np.mean(cc_array_rob)) # mean closeness centrality for robust rxns only
print(np.std(cc_array_rob))
print(np.median(cc_array_rob))
print(np.percentile(cc_array_rob, 75) - np.percentile(cc_array_rob, 25))

print(np.mean(cc_array_sens)) # mean closeness centrality for sensitive rxns only
print(np.std(cc_array_sens))
print(np.median(cc_array_sens))
print(np.percentile(cc_array_sens, 75) - np.percentile(cc_array_sens, 25))
'''

'''
# Betweeness centrality or bc for METABOLITES
bc = nx.algorithms.centrality.betweenness_centrality(g, normalized=False)
bc_array = np.array([]) # array of bc for all nodes
bc_array_R = np.array([]) # array of bc for reaction only
bc_array_rob = np.array([]) # array of cc for robust reactions only
bc_array_sens = np.array([]) # array of cc for sensitive reactions only
for n, val in bc.items():  # iterate over item names and '.values' in the dictionary
    bc_array = np.append(bc_array, val) # Create bc_array
    if n.startswith('M_') == True: 
        bc_array_R = np.append(bc_array_R, val)        
    if np.any(robust_mets == n): 
        bc_array_rob = np.append(bc_array_rob, val)
    if np.any(sens_mets == n): 
        bc_array_sens = np.append(bc_array_sens, val)  
    
print(np.mean(bc_array)) # mean betweeness centrality for all nodes
print(np.std(bc_array))
print(np.median(bc_array))
print(np.percentile(bc_array, 75) - np.percentile(bc_array, 25))

print(np.mean(bc_array_R)) # mean betweeness centrality for metabolites
print(np.std(bc_array_R))
print(np.median(bc_array_R))
print(np.percentile(bc_array_R, 75) - np.percentile(bc_array_R, 25))

print(np.mean(bc_array_rob)) # mean betweeness centrality for robust METABO only
print(np.std(bc_array_rob))
print(np.median(bc_array_rob))
print(np.percentile(bc_array_rob, 75) - np.percentile(bc_array_rob, 25))

print(np.mean(bc_array_sens)) # mean betweeness centrality for sensitive METABO only
print(np.std(bc_array_sens))
print(np.median(bc_array_sens))
print(np.percentile(bc_array_sens, 75) - np.percentile(bc_array_sens, 25))

# Plot histogram of frequency of betweeneess centrality
import matplotlib.pyplot as plt
#%matplotlib inline
plt.rcParams.update({'figure.figsize':(7,5), 'figure.dpi':100})

# title='Frequency Histogram of betweenness centrality: Robust or sensitive metabolites'
plt.hist(bc_array_sens, bins=50, range=[0, 10000], alpha = 0.4)
plt.hist(bc_array_rob, bins=50, range=[0, 10000], alpha=0.4)
plt.gca().set(ylabel='Frequency', xlabel='Betweenness centrality');
plt.gca().legend(('sensitive metabolites','robust metabolites'))
plt.savefig('Betweenness_centrality_start.png')

# title='Frequency Histogram of betweenness centrality: Robust or sensitive metabolites'
plt.hist(bc_array_sens, bins=50, range=[0, 700000])
plt.hist(bc_array_rob, bins=50, range=[0, 700000])
plt.gca().set(ylabel='Frequency', xlabel='Betweenness centrality');
plt.gca().legend(('sensitive metabolites','robust metabolites'))
plt.savefig('Betweenness_centrality_all.png')

plt.hist(bc_array_sens, bins=50, range=[10000, 700000])
plt.hist(bc_array_rob, bins=50, range=[10000, 700000])
plt.gca().set(ylabel='Frequency', xlabel='Betweenness centrality');
plt.gca().legend(('sensitive metabolites','robust metabolites'))
axes = plt.gca()
axes.set_xlim([9900, 700000])
plt.savefig('Betweenness_centrality_end.png')

# Plot relative frequency of betweenness centrality: Robust or sensitive metabolites'
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
res = stats.relfreq(bc_array_sens, numbins=100, defaultreallimits=(0,100000))
res2 = stats.relfreq(bc_array_rob, numbins=100, defaultreallimits=(0,100000))
# res.frequency
# np.sum(res.frequency)  
x = res.lowerlimit + np.linspace(0, res.binsize*res.frequency.size, res.frequency.size)
x2 = res2.lowerlimit + np.linspace(0, res2.binsize*res2.frequency.size, res2.frequency.size)
# Here we must also adjust "res2.lowerlimit" so that the bars do not overlay.

plt.rcParams.update({'figure.figsize':(7,5), 'figure.dpi':100})
p1 = plt.bar(x, res.frequency, width=res.binsize, alpha=0.4)
p2 = plt.bar(x2, res2.frequency, width=res2.binsize, alpha=0.4)
plt.title('Relative frequency histogram of betweenness centrality')
plt.ylabel('Relative frequency')
plt.xlabel('Betweenness centrality')
#plt.xlim(0,50000)
plt.legend((p1[0], p2[0]), ('Sensitive metabolites', 'Robust metabolites'))
# plt.show() # Do not add plt.show() before plt.savefig(), otherwise it won,t save
plt.savefig('Betweenness_central_rel1.png')

# We see that almost no changes in betweeneess centrality occur at values > 10 000

len(bc_array_rob)
len(bc_array_sens)
len(bc_array_rob) - len(bc_array_sens)
# We see that the 230 more metabolites in the sensitive sets are metabolites with very low betweenees centrality (< 2)


#from scipy import stats
#stats.cumfreq(bc_array_rob, numbins=10)

# Closeness centrality or cc for METABOLITES
cc = nx.algorithms.centrality.closeness_centrality(g)
cc_array = np.array([]) # array of bc for all nodes
cc_array_R = np.array([]) # array of bc for reaction only
cc_array_rob = np.array([]) # array of cc for robust reactions only
cc_array_sens = np.array([]) # array of cc for sensitive reactions only
for n, val in cc.items():  # iterate over item names and '.values' in the dictionary
    cc_array = np.append(cc_array, val) # Create bc_array
    if n.startswith('M_') == True: 
        cc_array_R = np.append(cc_array_R, val)        
    if np.any(robust_mets == n): 
        cc_array_rob = np.append(cc_array_rob, val)
    if np.any(sens_mets == n): 
        cc_array_sens = np.append(cc_array_sens, val)  
    
print(np.mean(cc_array)) # mean closeness centrality for all nodes
print(np.std(cc_array))
print(np.median(cc_array))
print(np.percentile(cc_array, 75) - np.percentile(cc_array, 25))

print(np.mean(cc_array_R)) # mean closeness centrality for METABOLITES
print(np.std(cc_array_R))
print(np.median(cc_array_R))
print(np.percentile(cc_array_R, 75) - np.percentile(cc_array_R, 25))

print(np.mean(cc_array_rob)) # mean closeness centrality for robust METABOLITES only
print(np.std(cc_array_rob))
print(np.median(cc_array_rob))
print(np.percentile(cc_array_rob, 75) - np.percentile(cc_array_rob, 25))

print(np.mean(cc_array_sens)) # mean closeness centrality for sensitive METABOLITES only
print(np.std(cc_array_sens))
print(np.median(cc_array_sens))
print(np.percentile(cc_array_sens, 75) - np.percentile(cc_array_sens, 25))

# Plot histogram of frequency of closeness centrality
import matplotlib.pyplot as plt
#%matplotlib inline
plt.rcParams.update({'figure.figsize':(7,5), 'figure.dpi':100})

# title='Frequency Histogram of closeness centrality: Robust or sensitive metabolites'
plt.hist(cc_array_sens, bins=50, range=[0, 0.25], alpha=0.4) #, rwidth=0.5)
plt.hist(cc_array_rob, bins=50, range=[0, 0.25], alpha=0.4) #, rwidth=0.5)
plt.gca().set(ylabel='Frequency', xlabel='Closeness centrality');
plt.gca().legend(('sensitive metabolites','robust metabolites'))
plt.savefig('Closeness_centrality.png')

# Plot relative frequency of closeness centrality: Robust or sensitive metabolites'
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
res = stats.relfreq(cc_array_sens, numbins=50, defaultreallimits=(0,0.25))
res2 = stats.relfreq(cc_array_rob, numbins=50, defaultreallimits=(0,0.25))
# res.frequency
# np.sum(res.frequency)  
x = res.lowerlimit + np.linspace(0, res.binsize*res.frequency.size, res.frequency.size)
x2 = res2.lowerlimit + np.linspace(0, res2.binsize*res2.frequency.size, res2.frequency.size)
# Here we must also adjust "res2.lowerlimit" so that the bars do not overlay.

plt.rcParams.update({'figure.figsize':(7,5), 'figure.dpi':100})
p1 = plt.bar(x, res.frequency, width=res.binsize, alpha=0.4)
p2 = plt.bar(x2, res2.frequency, width=res2.binsize, alpha=0.4)
plt.title('Relative frequency histogram of closeness centrality')
plt.ylabel('Relative frequency')
plt.xlabel('Closeness centrality')
#plt.xlim(0,50000)
plt.legend((p1[0], p2[0]), ('Sensitive metabolites', 'Robust metabolites'))
# plt.show() # Do not add plt.show() before plt.savefig(), otherwise it won,t save
plt.savefig('Closeness_central_rel1.png')
 '''

# We see that robust metabolite sets have more much less (2? 3? times less?) metabolites of very low closeness centrality, but those metabolites represent 
# a small proportion of the total frequency of closeness centrality.

#from scipy import stats
#stats.describe(cc_array)
#############################################################################
