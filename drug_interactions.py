# We want to create a network visualization graph of all the known drug-drug interactions
# Based on tutorial here: https://towardsdatascience.com/network-analysis-and-visualization-of-drug-drug-interactions-1e0b41d0d3df
# Drug drug interaction dataset from here: http://snap.stanford.edu/biodata/datasets/10001/10001-ChCh-Miner.html
# Drug name lookup table from here: https://go.drugbank.com/releases/latest#open-data
# See also: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0247018

# First get the interaction data into a dataframe
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
df = pd.read_csv('ChCh-Miner_durgbank-chem-chem.tsv', sep='\t', names=["drug1", "drug2"])
# Each row is an interaction between two drugs

# THe names of the drugs are Drugbank IDs which are not very helpful at the moment
# LEts use a lookup table to convert these IDs to common names
# We will first change the lookup table from a csv file into a dictionary and use the dictionary to replace the drug IDs
df_lookup = pd.read_csv("drug_lookup.csv")
lookup_dict = dict(zip(df_lookup["DrugBank ID"], df_lookup["Common name"]))
df["drug1"] = df["drug1"].map(lookup_dict).fillna(df["drug1"])  # Replace all drug IDs using map and return NA for anything that is not stated in the dictionary 
df["drug2"] = df["drug2"].map(lookup_dict).fillna(df["drug2"])

# Manual checking shows that most of the conversion worked with a few only about 0.3% of drugs that were not present
print("Percentage of drug1 column that was not converted to common name: {}". format(len(df[df['drug1'].str.contains('DB\d*')==True])/len(df)))

# A downside of this dataset is that the severity of drug interactions is not provided. 
# Therefore, a custom column (named 'weight') was filled with 1’s was added to indicate ‘equal’ severity for all interactions.
df["weight"] = 1

#%%
# Basic analyses of drugs using networkx

# Generate a networkx graph from my pandas dataframe
G = nx.from_pandas_edgelist(df, 'drug1', 'drug2')
# Give the graph a name
G.name = 'Drug Interactions Network'
# Summary statistic of graph
print(nx.info(G))

# Lets make a function that finds the node(s) with the highest degree and betweeness centrality

def find_nodes_with_highest_deg_cent(G):
    """find the node that has the highest degree centrality"""
    # Compute the degree centrality of G: deg_cent
    deg_cent = nx.degree_centrality(G)
    # Compute the maximum degree centrality: max_dc
    max_dc = max(list(deg_cent.values()))
    nodes = set()
    # Iterate over the degree centrality dictionary
    for k, v in deg_cent.items():
        # Check if the current value has the maximum degree centrality
        if v == max_dc:
            # Add the current node to the set of nodes
            nodes.add(k)
    return nodes

# Find the node(s) that has the highest degree centrality: top_dc
top_dc = find_nodes_with_highest_deg_cent(G)
print("Node with the highest degree centrality:")
print(top_dc)


# Lets also look at betweenenss centrality
def find_node_with_highest_bet_cent(G):
    """find the node that has the highest betwenness centrality"""
    # Compute betweenness centrality: bet_cent
    bet_cent = nx.betweenness_centrality(G)
    # Compute maximum betweenness centrality: max_bc
    max_bc = max(list(bet_cent.values()))
    nodes = set()
    # Iterate over the betweenness centrality dictionary
    for k, v in bet_cent.items():
        # Check if the current value has the maximum betweenness centrality
        if v == max_bc:
            # Add the current node to the set of nodes
            nodes.add(k)
    return nodes

# Use that function to find the node(s) that has the highest betweenness centrality in the network: top_bc
top_bc = find_node_with_highest_bet_cent(G)
print("Node with the highest betweenness centrality:")
print(top_bc)

# Phenytoin has the highest degree centrality while warfarin has the highest betweenness centrality



#%%
# Network visualization
# Subset data with Kanamycin, Oxytetracycline and Streptomycin

