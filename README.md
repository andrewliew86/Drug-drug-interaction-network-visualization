# Visualizing drug-drug interactions 

Background: Drug-drug interactions can often cause unwanted adverse events in patients. A library of drug-drug interactions would therefore be a useful resource for any clinician or scientist or pharmacist when choosing a drug for their patient.

Results: I created a network visualization of the drug-drug interactions using publically available data from the Stanford Biomedical Network Dataset Collection. The visualization is interactive and can be further explored by the user. I also report on network characteristics (degree centrality, eigenvector centrality, betweenness centrality, page rank)

Python libraries/tools: Standard Python libraries and Networkx for network analysis. Pyvis was used for visualization.
Code was modified from :https://towardsdatascience.com/network-analysis-and-visualization-of-drug-drug-interactions-1e0b41d0d3df and https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0247018

Figure: Network map of drug-drug interactions for 4 different antibiotics: Kanamycin, Oxytetracycline, Streptomycin, Clarithromycin
![Network map of antibiotics](https://github.com/andrewliew86/network_viz/blob/master/network_antibiotics.PNG)


