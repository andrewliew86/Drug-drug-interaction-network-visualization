# We want to create a network visualization graph of all the known drug-drug interactions
# Based on tutorial here: https://towardsdatascience.com/network-analysis-and-visualization-of-drug-drug-interactions-1e0b41d0d3df
# Drug drug interaction dataset from here: http://snap.stanford.edu/biodata/datasets/10001/10001-ChCh-Miner.html
# Drug name lookup table from here: https://go.drugbank.com/releases/latest#open-data
# See also: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0247018

# First get the interaction data into a dataframe
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from pyvis.network import Network

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
# Basic analyses (centrality) of drugs using networkx

# Generate a networkx graph from my pandas dataframe
G = nx.from_pandas_edgelist(df, 'drug1', 'drug2')
# Give the graph a name
G.name = 'Drug Interactions Network'
# Summary statistic of graph
print(nx.info(G))

# Lets also create a dataframe  that has the different centrality measures as different columns so we can have a look.
# Note that the Betweeness centrality might take a while to calculate
centrality_df = pd.DataFrame({'Degree':pd.Series(nx.degree_centrality(G)),
                              'Betweeness':pd.Series(nx.betweenness_centrality(G)), 
                              'Eigenvector':pd.Series(nx.eigenvector_centrality(G)), 
                              'Pagerank': pd.Series(nx.pagerank(G))})

# Phenytoin has the highest degree and eigenvector centrality 
# warfarin has the highest betweenness centrality and Pagerank score

#%%
# Network visualization 
# Lets focus on some antibiotics: Kanamycin, Oxytetracycline and Streptomycin
# I use only  small subset because plotting too many drugs will be a very bad visualization - too dense a network to make any sense!
subset = ["Kanamycin", "Oxytetracycline", "Streptomycin", "Clarithromycin"]
df_anti_small = df.loc[df['drug1'].isin(subset) | df['drug2'].isin(subset)]
df_anti_small = df_anti_small.reset_index(drop=True)



def generate_network_viz(df, source_col, target_col, weights):
    """ Generate an interative network graph"""
    
    # Generate a networkx graph
    G = nx.from_pandas_edgelist(df, source_col, target_col, weights)

    # Initiate PyVis network object
    drug_net = Network(
                       height='700px', 
                       width='100%',
                       bgcolor="white", 
                       font_color="black", 
                       notebook=True
                      )
    
    # Take Networkx graph and translate it to a PyVis graph format
    drug_net.from_nx(G)
    
    # Set repulsion so the network is spread out
    drug_net.repulsion(
                            node_distance=420, 
                            central_gravity=0.15, 
                            spring_length=100, 
                            spring_strength=0.15, 
                            damping=0.96
                           )
    return drug_net

# Generate a networkx graph based on the antibiotic subset
db_subset_net = generate_network_viz(df_anti_small, 'drug1', 'drug2', 'weight')

# Display interactive graph
db_subset_net.show('drug_interactions_network_subset_repulsion.html')

# There are actually very few overlaps with the different drugs. 
# other interestin observations are that Magnesium sulphate interaction occurs with both tetracycline and clarithromycin 
