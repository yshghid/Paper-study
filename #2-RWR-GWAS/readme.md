Korean version is provided [here][2].

# Outline

1) Load GWAS catalog, perform SNP statistical tests
2) Set seed genes, target genes
3) Load graph
4) Set initial weights
5) Perform RWR algorithm

# Pipeline

## Load GWAS catalog, perform SNP statistical tests

Information on disease-associated SNPs can be downloaded in the format of a GWAS catalog. 
The catalog provides details such as the chromosomal position of each SNP, associated genes, statistical significance (p-value), and accession ID. 
In this data, there are 1283 SNPs and 20 genes (g1-g20).

```python
# Load gwas catalog
gwas_catalog = load_catalog()
```
```
gwas_catalog
>>
	SNP	pValue	mappedGenes
0	snp1	0.018565	g16, g20, g13, g6
1	snp2	0.006678	g2, g7, g16, g18
2	snp3	0.043292	g5, g11, g18, g19
3	snp4	0.031542	g10, g15, g20, g1
4	snp5	0.017214	g5, g16
...	...	...	...
1278	snp1279	0.045134	g7, g13, g19, g2
1279	snp1280	0.021664	g13
1280	snp1281	0.014565	g9
1281	snp1282	0.030025	g19, g1, g12, g8
1282	snp1283	0.045706	g10
1283 rows × 3 columns
```

There are various tools available for converting information on disease-associated SNPs into gene scores. 
These tools use slightly different statistical testing methods, and in this paper, PEGASUS was used. 
PEGASUS takes the GWAS catalog as input and generates a statistical significance (Pvalue) and score (Score) for each gene using the p-values of the disease-associated SNPs.

```
# Load pegasus output
pegasus_data = load_pegasus_results()
```
```
pegasus_data
>>
	Gene	Gene NCBI ID	Chr	nSNPs	Start	Stop	Score	Pvalue
0	g1	1	Chr7	73	5	10	1.147386	0.065060
1	g2	2	Chr20	69	15	20	2.949693	0.144366
2	g3	3	Chr15	30	25	30	1.400289	0.429186
3	g4	4	Chr11	42	35	40	2.579821	0.023038
4	g5	5	Chr8	85	45	50	2.040923	0.027445
5	g6	6	Chr21	67	55	60	1.351498	0.029751
6	g7	7	Chr7	31	65	70	0.039795	0.017000
7	g8	8	Chr19	98	75	80	2.826605	0.028952
8	g9	9	Chr11	58	85	90	1.689865	0.042382
9	g10	10	Chr11	100	95	100	1.156250	0.305897
10	g11	11	Chr21	68	105	110	0.047899	0.030745
11	g12	12	Chr4	51	115	120	0.692681	0.014524
12	g13	13	Chr8	69	125	130	0.723076	0.026521
13	g14	14	Chr3	89	135	140	2.049791	0.048298
14	g15	15	Chr22	24	145	150	1.829990	0.046406
15	g16	16	Chr21	71	155	160	2.499585	0.036636
16	g17	17	Chr2	71	165	170	0.520094	0.020362
17	g18	18	Chr12	56	175	180	1.173182	0.241895
18	g19	19	Chr6	71	185	190	0.546708	0.048097
19	g20	20	Chr2	60	195	200	2.266084	0.037619
```

## Set seed genes, target genes

In the PEGASUS data, genes with a p-value < 0.05 are selected as seed genes. 
Genes associated with SNPs mentioned in the GWAS catalog as related to the disease are selected as target genes. 
Since the p-values of 15 out of 20 genes in the PEGASUS data are <0.05, there are 15 seed genes and 20 target genes.

```python
def load_seeds_and_targets():
    # Load seeds
    pegasus_data = load_pegasus_results()
    gene_seeds_ncbi = pegasus_data.loc[pegasus_data['Pvalue'] <= 0.05, 'Gene NCBI ID'].tolist()

    # Load targets
    gwas_catalog = load_catalog()
    ncbi_targets = set() 
    for i, row in gwas_catalog.iterrows():
        gns = row["mappedGenes"].split(", ")
        for gn in gns:
            if gn in pegasus_data['Gene'].values:
                ncbi_targets.add(pegasus_data.loc[pegasus_data['Gene'] == gn, 'Gene NCBI ID'].iloc[0])
        ncbi_targets_sorted = sorted(list(ncbi_targets), key=lambda x: int(x[0:]))
    
    return gene_seeds_ncbi, ncbi_targets_sorted

seeds, targets = load_seeds_and_targets()
```
```
seeds
>>
['4',
 '5',
 '6',
 '7',
 '8',
 '9',
 '11',
 '12',
 '13',
 '14',
 '15',
 '16',
 '17',
 '19',
 '20']
```
```
targets
>>
['1',
 '2',
 '3',
 '4',
 '5',
 '6',
 '7',
 '8',
 '9',
 '10',
 '11',
 '12',
 '13',
 '14',
 '15',
 '16',
 '17',
 '18',
 '19',
 '20']
```

## Load graph

PPI network consists of information on source nodes and target nodes. 
The nodes are genes, typically represented by NCBI IDs. 
This network has about 100 edges for 20 nodes.

```
# Generate graph
def load_graph_nx():
    # Generate 100 rows of node pairs
    num_rows = 100
    num_nodes = 20
    min_frequency = 10
    node_pairs = generate_node_pairs(num_rows, num_nodes, min_frequency)

    # Convert to DataFrame
    node_data = pd.DataFrame(node_pairs, columns=["node1", "node2"])
    node_data_without_duplicates = node_data.drop_duplicates()

    graph = from_pandas_edgelist(node_data_without_duplicates, source="node1", target="node2")
    
    return graph, node_data_without_duplicates

graph, node = load_graph_nx()
node
```
```
	node1	node2
0	8	12
1	12	15
2	6	20
3	3	18
4	9	11
...	...	...
94	8	11
95	16	17
96	2	18
98	12	20
99	14	15
80 rows × 2 columns
```

## Set initial weights

Initial weights for the seeds are set in the loaded network using the Scores obtained from the PEGASUS data. 
The seeds were the 15 genes obtained from PEGASUS data, and the initial weight is set to zero for genes that are not seed genes.

```
def init_rwr_scores_nx(graph, data):   
    pegasus_genes = set(data['Gene NCBI ID'])
    pegasus_scores = dict(zip(data['Gene NCBI ID'], data['Score']))
    
    pagerank_seeds = {}
    for node in graph.nodes:
        if node in pegasus_genes:
            pagerank_seeds[node] = pegasus_scores[node]
        else:
            pagerank_seeds[node] = 0

    return pagerank_seeds

pagerank_seeds = init_rwr_scores_nx(graph, pegasus_data)
```
```
pagerank_seeds
>>
{'8': 2.8266052670545583,
 '12': 0.692681476866447,
 '15': 1.8299899733478626,
 '6': 1.351497755908629,
 '20': 2.2660842309529574,
 '3': 0,
 '18': 0,
 '9': 1.6898646535366177,
 '11': 0.04789875666064258,
 '16': 2.499584735208493,
 '17': 0.5200939605233162,
 '13': 0.7230763980780351,
 '19': 0.546708263364187,
 '7': 0.039794883479599585,
 '10': 0,
 '4': 2.5798212202089617,
 '14': 2.0497905564763745,
 '1': 0,
 '2': 0,
 '5': 2.040922615763339}
 ```

## Perform RWR algorithm

The Random Walk with Restart (RWR) algorithm used the PageRank function. 
PageRank is an algorithm developed by Google to calculate the importance of web pages, where each node (e.g., web page) receives weights from other connected nodes. 
In PageRank, initial weights are uniformly distributed across all nodes and the walk starts from all nodes, but in RWR, specific weights are assigned to seed nodes and the walk starts from these seed nodes.

In PageRank, α represents the probability of moving to a random node. A higher value allows more transitions between nodes. 
In RWR, α represents the restart probability, indicating the likelihood of returning to the seed nodes and emphasizing relationships stemming from the seed nodes. 
The paper mentions that with α = 0.1, the walk essentially propagated only to the one-hop neighborhood of the disease genes, i.e., the seed nodes, and with α = 0.3, the influence propagated to two- or three-hop neighborhoods, and it was performed with α = 0.3.

```
def perform_rwr_nx(alpha, graph, seeds):   
    rwr_scores = pagerank(graph, alpha=alpha, personalization=seeds)
    return rwr_scores

alpha = 0.3
rwr_scores = perform_rwr_nx(alpha, graph, pagerank_seeds)
```
```
rwr_scores
>>
{'8': 0.10426855672164884,
 '12': 0.04425145356045913,
 '15': 0.07432322860869753,
 '6': 0.06815092036489143,
 '20': 0.09004755193034297,
 '3': 0.013645421197650718,
 '18': 0.023150821508225422,
 '9': 0.07581465423125429,
 '11': 0.016890247939529115,
 '16': 0.09547155745649892,
 '17': 0.02786246477943842,
 '13': 0.0382453415081826,
 '19': 0.028869143467498616,
 '7': 0.011678794362311333,
 '10': 0.006413068843695989,
 '4': 0.09938439277021324,
 '14': 0.07466026734035734,
 '1': 0.0171008992699609,
 '2': 0.01208603866611786,
 '5': 0.07768517547302517}
```

The results are summarized as follows.

```
rwr_results = process_rwr_results_nx(rwr_scores, graph, pegasus_data, pagerank_seeds)
rwr_results
```
```
Gene NCBI ID	Gene	Initial Score	Final Score
0	8	g8	2.826605	0.104269
15	4	g4	2.579821	0.099384
9	16	g16	2.499585	0.095472
4	20	g20	2.266084	0.090048
19	5	g5	2.040923	0.077685
7	9	g9	1.689865	0.075815
16	14	g14	2.049791	0.074660
2	15	g15	1.829990	0.074323
3	6	g6	1.351498	0.068151
1	12	g12	0.692681	0.044251
11	13	g13	0.723076	0.038245
12	19	g19	0.546708	0.028869
10	17	g17	0.520094	0.027862
6	18	g18	0.000000	0.023151
17	1	g1	0.000000	0.017101
8	11	g11	0.047899	0.016890
5	3	g3	0.000000	0.013645
18	2	g2	0.000000	0.012086
13	7	g7	0.039795	0.011679
14	10	g10	0.000000	0.006413
```

Link: [Network propagation for GWAS analysis: a practical guide to leveraging molecular networks for disease gene discovery][1]

[1]: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10858647/
[2]: https://github.com/yshghid/Paper-study/blob/main/%232-RWR-GWAS/readme-kor.md#outline
