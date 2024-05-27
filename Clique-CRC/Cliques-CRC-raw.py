import numpy as np
import pandas as pd
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests
import networkx as nx

# Set random seed for reproducibility
np.random.seed(42)

# Generate synthetic data with baseline and noise
def generate_data(num_genes=100, num_samples=10, patients=5, controls=5):
    assert patients + controls == num_samples
    # Baseline expression levels for patients and controls
    baseline_patient = np.random.normal(loc=50, scale=10, size=num_genes)
    baseline_control = np.random.normal(loc=30, scale=10, size=num_genes)
    
    # Generate data for patients and controls with noise
    patient_data = np.array([np.random.negative_binomial(max(1, int(base)), 0.5, size=patients) for base in baseline_patient])
    control_data = np.array([np.random.negative_binomial(max(1, int(base)), 0.5, size=controls) for base in baseline_control])
    
    data = np.concatenate([patient_data, control_data], axis=1)
    columns = [f'p{i+1}' for i in range(patients)] + [f'c{i+1}' for i in range(controls)]
    df = pd.DataFrame(data, index=[f'g{i+1}' for i in range(num_genes)], columns=columns)
    return df

# Create datasets
datasets = {
    'USA': generate_data(),
    'CHN': generate_data(),
    'SA': generate_data(),
    'GER': generate_data()
}

# Normalize the data using log transformation (simulating RMA normalization)
def normalize_data(df):
    return np.log1p(df)

normalized_datasets = {k: normalize_data(v) for k, v in datasets.items()}

# Perform two-sample t-test
def perform_t_test(df):
    patients = df.columns[:5]
    controls = df.columns[5:]
    t_results = []
    for gene in df.index:
        t_stat, p_val = stats.ttest_ind(df.loc[gene, patients], df.loc[gene, controls])
        t_results.append((gene, t_stat, p_val))
    t_results_df = pd.DataFrame(t_results, columns=['Gene', 't_stat', 'p_val'])
    return t_results_df

t_test_results = {k: perform_t_test(v) for k, v in normalized_datasets.items()}

# Apply multiple testing correction (FDR)
def apply_fdr(df):
    p_vals = df['p_val'].values
    _, q_vals, _, _ = multipletests(p_vals, alpha=0.1, method='fdr_bh')
    df['q_val'] = q_vals
    return df

corrected_results = {k: apply_fdr(v) for k, v in t_test_results.items()}

# Filter genes based on p-value and q-value thresholds
def filter_genes(df, p_val_threshold=0.05, q_val_threshold=0.1):
    return df[(df['p_val'] < p_val_threshold) & (df['q_val'] < q_val_threshold)]

filtered_genes = {k: filter_genes(v) for k, v in corrected_results.items()}

# Create graphs for each dataset
#def create_graph(filtered_genes, normalized_data):
#    graph = nx.Graph()
#    genes = filtered_genes['Gene']
#    for gene in genes:
#        graph.add_node(gene)
#    for i in range(len(genes)):
#        for j in range(i + 1, len(genes)):
#            gene1 = genes.iloc[i]
#            gene2 = genes.iloc[j]
            # Adding edge with weight based on correlation of expression levels
#            weight = np.corrcoef(normalized_data.loc[gene1], normalized_data.loc[gene2])[0, 1]
#            graph.add_edge(gene1, gene2, weight=weight)
#    return graph

def create_graph(filtered_genes, normalized_data, threshold=0.7):
    graph = nx.Graph()
    genes = filtered_genes['Gene']
    for gene in genes:
        graph.add_node(gene)
    for i in range(len(genes)):
        for j in range(i + 1, len(genes)):
            gene1 = genes.iloc[i]
            gene2 = genes.iloc[j]
            # Adding edge with weight based on correlation of expression levels if above threshold
            weight = np.corrcoef(normalized_data.loc[gene1], normalized_data.loc[gene2])[0, 1]
            if abs(weight) > threshold:
                graph.add_edge(gene1, gene2, weight=weight)
    return graph

graph_dic = {k: create_graph(filtered_genes[k], normalized_datasets[k]) for k, v in filtered_genes.items()}

# Node strength analysis
def node_strength_analysis(graph):
    # Calculate Degree
    degree = nx.degree_centrality(graph)

    # Calculate Eccentricity
    eccentricity = nx.eccentricity(graph)
    #try:
    #    eccentricity = nx.eccentricity(graph)
    #except nx.NetworkXError:
    #    eccentricity = {node: 0 for node in graph.nodes()}

    # Calculate Closeness
    closeness = nx.closeness_centrality(graph)

    # Calculate Betweenness
    betweenness = nx.betweenness_centrality(graph)

    # Normalize values
    degree_max = max(degree.values())
    eccentricity_max = max(eccentricity.values())
    closeness_max = max(closeness.values())
    betweenness_max = max(betweenness.values())
    #degree_max = max(degree.values()) if degree else 1
    #eccentricity_max = max(eccentricity.values()) if eccentricity else 1
    #closeness_max = max(closeness.values()) if closeness else 1
    #betweenness_max = max(betweenness.values()) if betweenness else 1

    #degree_normalized = {k: v / degree_max for k, v in degree.items()}
    #eccentricity_normalized = {k: v / eccentricity_max for k, v in eccentricity.items()}
    #closeness_normalized = {k: v / closeness_max for k, v in closeness.items()}
    #betweenness_normalized = {k: v / betweenness_max for k, v in betweenness.items()}
    degree_normalized = {k: (v / degree_max if degree_max != 0 else 0) for k, v in degree.items()}
    eccentricity_normalized = {k: (v / eccentricity_max if eccentricity_max != 0 else 0) for k, v in eccentricity.items()}
    closeness_normalized = {k: (v / closeness_max if closeness_max != 0 else 0) for k, v in closeness.items()}
    betweenness_normalized = {k: (v / betweenness_max if betweenness_max != 0 else 0) for k, v in betweenness.items()}


    # Calculate Node Strength
    node_strength = {
        node: (degree_normalized[node] + eccentricity_normalized[node] + closeness_normalized[node] + betweenness_normalized[node]) / 4
        for node in graph.nodes()
    }

    # Convert results to DataFrame for better readability
    node_strength_df = pd.DataFrame({
        'Node': list(node_strength.keys()),
        'Degree': [degree[node] for node in node_strength.keys()],
        'Eccentricity': [eccentricity[node] for node in node_strength.keys()],
        'Closeness': [closeness[node] for node in node_strength.keys()],
        'Betweenness': [betweenness[node] for node in node_strength.keys()],
        'NodeStrength': list(node_strength.values())
    })

    return node_strength_df

# Perform node strength analysis for each graph
node_strength_results = {k: node_strength_analysis(v) for k, v in graph_dic.items()}
