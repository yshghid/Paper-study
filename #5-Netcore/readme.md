Korean version is provided [here][2].

# Outline

1. Load PPI network
2. Normalize PPI network using node coreness
3. Run RWR, identify significantly re-ranked nodes
4. Select of parameter α
5. Comparison of four normalization methods through 5-fold cross-validation
6. Identify gene modules

# Pipeline

## Load PPI network
Load the PPI network for analysis. The PPI network consists of 20 nodes (g1-g20) and 44 edges.

```python
network = generate_network()
edges = generate_edges(network)
```
```
network.nodes()
>>
NodeView(('g1', 'g2', 'g3', 'g4', 'g5', 'g6', 'g7', 'g8', 'g9', 'g10', 'g11', 'g12', 'g13', 'g14', 'g15', 'g16', 'g17', 'g18', 'g19', 'g20'))
```
```
len(edges)
>>
44
```

The core number of each node is calculated. 
A k-core is the largest subgraph where all nodes have at least k neighbors, meaning each node's degree is at least k. 
The core number is the value of k in the largest k-core to which the node belongs.

While degree simply refers to the number of connected neighbors, the core number is a more complex indicator that shows how close a node is to the central part of the graph. 
The core number reflects how close a node is to the center of the graph, taking into consideration the overall structure of the graph and the position of each node.

```python
def generate_nodes(G):
    node_core_numbers = nx.core_number(G)
    return node_core_numbers

nodes = generate_nodes(network)
```
```
nodes
>>
{'g1': 4,
 'g2': 4,
 'g3': 4,
 'g4': 4,
 'g5': 4,
 'g6': 2,
 'g7': 4,
 'g8': 4,
 'g9': 2,
 'g10': 1,
 'g11': 3,
 'g12': 2,
 'g13': 3,
 'g14': 1,
 'g15': 3,
 'g16': 1,
 'g17': 2,
 'g18': 4,
 'g19': 3,
 'g20': 2}
```

### Normalize PPI network using node coreness

Generate the adjacency matrix for the PPI network and normalize the PPI network. 
There are two methods of normalization: one using the conventional degree and another using three types of coreness as proposed in the paper.

```python
adj_matrix = generate_adj_matrix(network)

# Function to normalize based on degree
def normalize_by_degree(adj_matrix):
    degree = adj_matrix.sum(axis=1)
    degree[degree == 0] = 1  # Avoid division by zero
    W_degree = adj_matrix.values @ np.diag(1 / degree)
    W_degree_df = pd.DataFrame(W_degree, index=adj_matrix.index, columns=adj_matrix.columns)
    return W_degree_df.div(W_degree_df.sum(axis=0), axis=1)

# Function to normalize based on core
def normalize_by_core(adj_matrix, cores):
    core_values = np.array([cores[node] for node in adj_matrix.columns])
    W_core = np.zeros(adj_matrix.shape)
    for i in range(len(adj_matrix.columns)):
        neighbors = adj_matrix.iloc[:, i] != 0
        core_sum = core_values[neighbors].sum()
        if core_sum > 0:
            W_core[:, i] = core_values / core_sum * neighbors
    W_core_df = pd.DataFrame(W_core, index=adj_matrix.index, columns=adj_matrix.columns)
    return W_core_df.div(W_core_df.sum(axis=0), axis=1)

# Function to normalize based on degree-core difference
def normalize_by_diff(adj_matrix, cores):
    degree = adj_matrix.sum(axis=1)
    W_diff = np.zeros(adj_matrix.shape)
    for i, node in enumerate(adj_matrix.columns):
        neighbors = adj_matrix.iloc[:, i] != 0
        for j, neighbor in enumerate(adj_matrix.index):
            if neighbors[j]:
                diff = degree[node] - cores[node]
                if diff > 0:
                    W_diff[j, i] = adj_matrix.iloc[j, i] / diff
                else:
                    W_diff[j, i] = adj_matrix.iloc[j, i]
        column_sum = W_diff[:, i].sum()
        if column_sum > 0:
            W_diff[:, i] /= column_sum

    W_diff_df = pd.DataFrame(W_diff, index=adj_matrix.index, columns=adj_matrix.columns)
    return W_diff_df.div(W_diff_df.sum(axis=0), axis=1)

# Function to normalize based on degree/core ratio
def normalize_by_ratio(adj_matrix, cores):
    degree = adj_matrix.sum(axis=1)
    W_ratio = np.zeros(adj_matrix.shape)
    for i, node in enumerate(adj_matrix.columns):
        W_ratio[:, i] = (np.array(list(cores.values())) / degree[node]) * adj_matrix.iloc[:, i]
    W_ratio /= W_ratio.sum(axis=0)
    W_ratio_df = pd.DataFrame(W_ratio, index=adj_matrix.index, columns=adj_matrix.columns)
    return W_ratio_df.div(W_ratio_df.sum(axis=0), axis=1)

# Function to normalize
def normalize_adj_matrix(adj_matrix, cores, norm_method):
    if norm_method == 'degree':
        norm_adj_matrix = normalize_by_degree(adj_matrix)
    elif norm_method == 'core':
        norm_adj_matrix = normalize_by_core(adj_matrix, cores)
    elif norm_method == 'diff':
        norm_adj_matrix = normalize_by_diff(adj_matrix, cores)
    elif norm_method == 'ratio':
        norm_adj_matrix = normalize_by_ratio(adj_matrix, cores)
    else:
        raise ValueError("Invalid normalization method.")
    
    return norm_adj_matrix
```

For example, the results of normalization using the core number would be as follows.

```python
norm_method = 'core'
norm_adj_matrix = normalize_adj_matrix(adj_matrix, nodes, norm_method)
norm_adj_matrix
```
```
	g1	g2	g3	g4	g5	g6	g7	g8	g9	g10	g11	g12	g13	g14	g15	g16	g17	g18	g19	g20
g1	0.000000	0.000000	0.121212	0.114286	0.129032	0.5	0.25	0.25	0.5	0.0	0.000000	0.0	0.000000	1.0	0.000000	0.0	0.5	0.25	0.000000	0.0
g2	0.000000	0.000000	0.121212	0.114286	0.129032	0.0	0.25	0.00	0.0	1.0	0.333333	0.5	0.333333	0.0	0.333333	0.0	0.0	0.00	0.333333	0.5
g3	0.129032	0.121212	0.000000	0.114286	0.000000	0.5	0.25	0.25	0.0	0.0	0.000000	0.5	0.000000	0.0	0.333333	0.0	0.0	0.25	0.000000	0.5
g4	0.129032	0.121212	0.121212	0.000000	0.000000	0.0	0.25	0.25	0.5	0.0	0.333333	0.0	0.333333	0.0	0.000000	0.0	0.0	0.25	0.333333	0.0
g5	0.129032	0.121212	0.000000	0.000000	0.000000	0.0	0.00	0.25	0.0	0.0	0.333333	0.0	0.333333	0.0	0.333333	1.0	0.5	0.25	0.333333	0.0
g6	0.064516	0.000000	0.060606	0.000000	0.000000	0.0	0.00	0.00	0.0	0.0	0.000000	0.0	0.000000	0.0	0.000000	0.0	0.0	0.00	0.000000	0.0
g7	0.129032	0.121212	0.121212	0.114286	0.000000	0.0	0.00	0.00	0.0	0.0	0.000000	0.0	0.000000	0.0	0.000000	0.0	0.0	0.00	0.000000	0.0
g8	0.129032	0.000000	0.121212	0.114286	0.129032	0.0	0.00	0.00	0.0	0.0	0.000000	0.0	0.000000	0.0	0.000000	0.0	0.0	0.00	0.000000	0.0
g9	0.064516	0.000000	0.000000	0.057143	0.000000	0.0	0.00	0.00	0.0	0.0	0.000000	0.0	0.000000	0.0	0.000000	0.0	0.0	0.00	0.000000	0.0
g10	0.000000	0.030303	0.000000	0.000000	0.000000	0.0	0.00	0.00	0.0	0.0	0.000000	0.0	0.000000	0.0	0.000000	0.0	0.0	0.00	0.000000	0.0
g11	0.000000	0.090909	0.000000	0.085714	0.096774	0.0	0.00	0.00	0.0	0.0	0.000000	0.0	0.000000	0.0	0.000000	0.0	0.0	0.00	0.000000	0.0
g12	0.000000	0.060606	0.060606	0.000000	0.000000	0.0	0.00	0.00	0.0	0.0	0.000000	0.0	0.000000	0.0	0.000000	0.0	0.0	0.00	0.000000	0.0
g13	0.000000	0.090909	0.000000	0.085714	0.096774	0.0	0.00	0.00	0.0	0.0	0.000000	0.0	0.000000	0.0	0.000000	0.0	0.0	0.00	0.000000	0.0
g14	0.032258	0.000000	0.000000	0.000000	0.000000	0.0	0.00	0.00	0.0	0.0	0.000000	0.0	0.000000	0.0	0.000000	0.0	0.0	0.00	0.000000	0.0
g15	0.000000	0.090909	0.090909	0.000000	0.096774	0.0	0.00	0.00	0.0	0.0	0.000000	0.0	0.000000	0.0	0.000000	0.0	0.0	0.00	0.000000	0.0
g16	0.000000	0.000000	0.000000	0.000000	0.032258	0.0	0.00	0.00	0.0	0.0	0.000000	0.0	0.000000	0.0	0.000000	0.0	0.0	0.00	0.000000	0.0
g17	0.064516	0.000000	0.000000	0.000000	0.064516	0.0	0.00	0.00	0.0	0.0	0.000000	0.0	0.000000	0.0	0.000000	0.0	0.0	0.00	0.000000	0.0
g18	0.129032	0.000000	0.121212	0.114286	0.129032	0.0	0.00	0.00	0.0	0.0	0.000000	0.0	0.000000	0.0	0.000000	0.0	0.0	0.00	0.000000	0.0
g19	0.000000	0.090909	0.000000	0.085714	0.096774	0.0	0.00	0.00	0.0	0.0	0.000000	0.0	0.000000	0.0	0.000000	0.0	0.0	0.00	0.000000	0.0
g20	0.000000	0.060606	0.060606	0.000000	0.000000	0.0	0.00	0.00	0.0	0.0	0.000000	0.0	0.000000	0.0	0.000000	0.0	0.0	0.00	0.000000	0.0
```

### Run RWR, identify significantly re-ranked nodes

Run the Random Walk with Restart (RWR) algorithm on the normalized PPI network. 
Select disease-associated genes as seed nodes, setting the initial weights of the seed nodes to 1 and all other nodes to 0. 
By comparing the initial weights (p0) with the final propagation weights (pk), nodes can be re-ranked.

To assign significance levels to the re-ranking results, compare the propagation weights obtained with those from randomized networks. 
The network randomization employs the double-edge swap algorithm to generate n=100 randomized networks. 
These are then compared with pk to assign significance levels to each node.

```python
# Disease related genes (seed genes)
seeds = ['g2', 'g4','g6','g8','g10']

# Function to run RWR
def rwr(norm_adj_matrix, seeds, alpha, max_iter=100, tol=1e-6):
    p0 = np.zeros(norm_adj_matrix.shape[0])
    seed_indices = [norm_adj_matrix.index.get_loc(seed) for seed in seeds]
    p0[seed_indices] = 1
    W = norm_adj_matrix.values
    n = W.shape[0]
    pk = p0.copy()
    for _ in range(max_iter):
        pk_new = alpha * p0 + (1 - alpha) * W @ pk
        if np.linalg.norm(pk_new - pk, 1) < tol:
            break
        pk = pk_new
    return pk

# Function to run random RWR
def random_rwr(network, seeds, alpha, n_random_networks=100):
    pks_random = []
    for _ in range(n_random_networks):
        G_random = network.copy()
        nx.double_edge_swap(G_random, nswap=len(G_random.edges()), max_tries=len(G_random.edges()) * 10)

        core_random = generate_nodes(G_random)
        adj_matrix_random = generate_adj_matrix(G_random)
        norm_adj_random = normalize_adj_matrix(adj_matrix_random, core_random, norm_method)
        
        pk_random = rwr(norm_adj_random, seeds, alpha, max_iter=100, tol=1e-6)
        pks_random.append(pk_random)

    pks_random = np.array(pks_random)
    return pks_random

# Function to calculate P-values for each node
def calculate_p_values(network, pk_original, pks_random):
    nodes = generate_nodes(network)
    p_values = []
    for i, node in enumerate(nodes):
        original_weight = pk_original[i]
        random_weights = pks_random[:, i]
        p_value = np.sum(random_weights >= original_weight) / pks_random.shape[0]
        p_values.append(p_value)
    return p_values
```

For example, perform the RWR with α = 0.1 and identify significantly re-ranked nodes as follows. 
While the paper applies a p-value threshold of <0.0099, the code uses a threshold of <0.3.

```python
# Run RWR algorithm for α = 0.1
alpha = 0.1

pk = rwr(norm_adj_matrix, seeds, alpha=alpha)
pks_random = random_rwr(network, seeds, alpha=alpha)
p_values = calculate_p_values(network, pk, pks_random)
    
print(f"\nRWR Result (α = {alpha}):")
for node, score in zip(norm_adj_matrix.index, pk):
        print(f"{node}: {score:.4f}")
    
# Identify significantly re-ranked nodes
significance_level = 0.3
significant_nodes = [node for node, p_value in zip(norm_adj_matrix.index, p_values) if p_value < significance_level]
        
print(f"\nSignificantly re-ranked nodes (α = {alpha}):")
for node in significant_nodes:
    print(node)
```
```
RWR Result (α = 0.1):
g1: 0.5260
g2: 0.6989
g3: 0.5806
g4: 0.6608
g5: 0.4917
g6: 0.1622
g7: 0.2686
g8: 0.3495
g9: 0.0645
g10: 0.1191
g11: 0.1510
g12: 0.0698
g13: 0.1510
g14: 0.0153
g15: 0.1475
g16: 0.0143
g17: 0.0591
g18: 0.2495
g19: 0.1510
g20: 0.0698

Significantly re-ranked nodes (α = 0.1):
g7
```

### Select of parameter α

The algorithm was performed with α values of 0.3, 0.5, and 0.8. 

In the Random Walk with Restart (RWR), α represents the restart parameter; the closer the value is to 1, the higher the probability that the walk will return to the seed node. 
Thus, the weight of the seed node remains close to 1, and the influence of the seed node is primarily propagated to nodes that are near the seed.

```python
# Run RWR algorithm for α = 0.3, 0.5, 0.8
alphas = [0.3, 0.5, 0.8]
for alpha in alphas:
    pk = rwr(norm_adj_matrix, seeds, alpha=alpha)
    pks_random = random_rwr(network, seeds, alpha=alpha)
    p_values = calculate_p_values(network, pk, pks_random)
    
    print(f"\nRWR Result (α = {alpha}):")
    #for node, score in zip(norm_adj_matrix.index, pk):
    #    print(f"{node}: {score:.4f}")
    
    # Identify significantly re-ranked nodes
    significance_level = 0.3
    significant_nodes = [node for node, p_value in zip(norm_adj_matrix.index, p_values) if p_value < significance_level]
        
    print(f"\nSignificantly re-ranked nodes (α = {alpha}):")
    for node in significant_nodes:
        print(node)
```
```
RWR Result (α = 0.3):

Significantly re-ranked nodes (α = 0.3):
g7
g11
g12
g13
g19
g20

RWR Result (α = 0.5):

Significantly re-ranked nodes (α = 0.5):
g2
g3
g7
g11
g12
g13
g19
g20

RWR Result (α = 0.8):

Significantly re-ranked nodes (α = 0.8):
g2
g7
g11
g12
g19
g20
```

Adjusting α allows for balancing between finding novel disease-associated genes and including potential false predictions. 
In this paper, an α value of 0.8 was chosen to minimize false positives, optimizing the likelihood that the influence remains concentrated around the seed nodes, thereby reducing the chance of erroneously identifying irrelevant nodes as significant.

### Comparison of four normalization methods through 5-fold cross-validation

To compare the performance of four different normalization methods in identifying disease-associated genes, the disease genes were split in a 1:4 ratio to create training sets. 
Given that there were five disease genes, four genes were used as the training set in each of the four network propagations. 
The α parameter was set to 0.8 for these experiments, which aimed to maximize the specificity of the propagation and minimize false positives by focusing the influence primarily around the training set genes.

```python
# Function to perform 5-fold cross-validation and calculate ROC curves for 4 normalization methods
def cross_validation_and_roc(network, disease_genes, alpha, n_splits=5):
    
    edges = generate_edges(network)
    nodes = generate_nodes(network)
    adj_matrix = generate_adj_matrix(network)

    norm_methods = ['degree', 'core', 'diff', 'ratio']
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)
    
    all_fpr = []
    all_tpr = []
    all_auroc = []

    for norm_method in norm_methods:
        norm_adj_matrix = normalize_adj_matrix(adj_matrix, nodes, norm_method=norm_method)
        
        tprs = []
        fprs = []
        aurocs = []
        
        for train_index, test_index in kf.split(disease_genes):
            train_genes = [disease_genes[i] for i in train_index]
            test_genes = [disease_genes[i] for i in test_index]
            
            p0 = np.zeros(len(norm_adj_matrix))
            p0[[norm_adj_matrix.index.get_loc(gene) for gene in train_genes]] = 1.0 / len(train_genes)
            
            pk = rwr(norm_adj_matrix, seeds, alpha=alpha, max_iter=100, tol=1e-6)
            
            pks_random = random_rwr(network, seeds, alpha=alpha, n_random_networks=100)
            
            p_values = calculate_p_values(network, pk, pks_random)
            
            y_true = np.isin(norm_adj_matrix.index, test_genes).astype(int)
            y_scores = -np.array(p_values)
            
            fpr, tpr, _ = roc_curve(y_true, y_scores)
            interp_fpr = np.linspace(0, 1, 100)
            interp_tpr = interp1d(fpr, tpr, kind='linear')(interp_fpr)
            
            tprs.append(interp_tpr)
            fprs.append(interp_fpr)
            aurocs.append(auc(fpr, tpr))

        mean_tpr = np.mean(tprs, axis=0)
        mean_fpr = np.mean(fprs, axis=0)
        mean_auc = np.mean(aurocs)
        
        all_tpr.append(mean_tpr)
        all_fpr.append(mean_fpr)
        all_auroc.append(mean_auc)

        plt.plot(mean_fpr, mean_tpr, label=f'{norm_method} (AUC = {mean_auc:.2f})')

    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('ROC Curve for Different Normalization Methods')
    plt.legend()
    plt.show()
    
    return all_fpr, all_tpr, all_auroc

# Disease-related genes (seed genes)
seeds = ['g2', 'g4','g6','g8','g10']

# Perform cross-validation and plot ROC curves
all_fpr, all_tpr, all_auroc = cross_validation_and_roc(network, seeds, alpha=0.8)
```
![image](https://github.com/yshghid/Paper-study/assets/153489198/12e7ea13-8aea-4a3f-842c-a43cd9a81c55)

```python
# Display results
for norm_method, auroc in zip(['degree', 'core', 'diff', 'ratio'], all_auroc):
    print(f"{norm_method} normalization AUROC: {auroc:.4f}")
```
```
degree normalization AUROC: 0.3632
core normalization AUROC: 0.4000
diff normalization AUROC: 0.4105
ratio normalization AUROC: 0.3684
```

The analysis results showed that the 'diff' normalization, which applies the difference between the node's degree and core number, exhibited the highest AUROC. 
However, in this paper, it was the 'core' normalization that achieved the highest AUROC, suggesting that considering the core number as a primary factor in the normalization process might be particularly effective in identifying disease-associated genes.

### Identify gene modules

Initially, a seed subnetwork is created, consisting of nodes that have seeds and edges. 
After performing network propagation, if there are nodes that meet specific criteria (p-value < 0.01, w > w_min) and have edges connecting to the initial seed subnetwork, these nodes are added. 
This process results in the formation of various modules within the final seed subnetwork.

```python
# Extension of seed modules using P-values and propagation weights
def identify_modules(network, seeds, norm_method, alpha, p_threshold, wmin_percentile):

    edges = generate_edges(network)
    nodes = generate_nodes(network)
    adj_matrix = generate_adj_matrix(network)
    norm_adj_matrix = normalize_adj_matrix(adj_matrix, nodes, norm_method)

    pk = rwr(norm_adj_matrix, seeds, alpha=alpha)
    pks_random = random_rwr(network, seeds, alpha=alpha, n_random_networks=100)
    p_values = calculate_p_values(network, pk, pks_random)
    
    # Step i: Extract seed-induced sub-network
    seed_subnetwork = network.subgraph(seeds).copy()
    
    # Step ii: Extend seed-induced sub-network
    extended_subnetwork = seed_subnetwork.copy()
    
    # Get propagation weights and significant nodes
    significant_nodes = [node for node, p_val in zip(network.nodes, p_values) if p_val < p_threshold]
    
    # Calculate wmin based on the 75th percentile of propagation weights of significant nodes
    if significant_nodes:
        propagation_weights = np.array([pk[list(network.nodes).index(node)] for node in significant_nodes])
        wmin = np.percentile(propagation_weights, wmin_percentile)
    else:
        wmin = 0
    
    # Add nodes to the extended sub-network
    for node in significant_nodes:
        if node not in seeds and pk[list(network.nodes).index(node)] > wmin:
            neighbors = set(network.neighbors(node))
            if neighbors & set(seeds):
                extended_subnetwork.add_node(node)
                for neighbor in neighbors:
                    if neighbor in seeds or neighbor in extended_subnetwork.nodes:
                        extended_subnetwork.add_edge(node, neighbor, weight=network.edges[node, neighbor]['weight'])
    
    # Step iii: Identify modules as connected components
    modules = [component for component in nx.connected_components(extended_subnetwork)]
    
    return modules

# Set disease genes as seed genes
seeds = ['g2', 'g4','g6','g8','g10']

# Set normalization method
norm_method = 'core'

# Identify modules
modules = identify_modules(network, seeds, norm_method, alpha = 0.8, p_threshold=0.3, wmin_percentile=75)

for i, module in enumerate(modules):
    print(f"Module {i+1}: {module}")
```
```
Module 1: {'g8', 'g10', 'g3', 'g6', 'g4', 'g2'}
```

In the code, the criteria were set to a p-value < 0.3 and w_min = 0.75, resulting in the formation of a module containing six genes. 
The module includes disease-associated genes g2, g4, g6, g8, and g10, with g3 being added as well.

Initially, the seed subnetwork did not include any connections for g6 with other disease genes. 
However, after network propagation, g3, which was identified as having a significantly changed weight, was added to the network. 
This addition allowed g6 to form a single module with the other disease genes, effectively integrating it into a connected group within the network.

```
Module 1: {'g4', 'g10', 'g8', 'g2'}
Module 2: {'g6'}
```


Link: [NetCore: a network propagation approach using node coreness][1]





[1]: https://academic.oup.com/nar/article/48/17/e98/5879427?login=false
[2]: https://github.com/yshghid/Paper-study/blob/main/%235-Netcore/readme-kor.md
