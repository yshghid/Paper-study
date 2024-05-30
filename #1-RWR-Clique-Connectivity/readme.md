### 데이터셋 생성
   
4개 국가(USA, GER, CHN, SA)의 대장암(CRC) 데이터셋을 가져온다. 4개 데이터셋은 모두 행은 100개의 유전자(g1-g100), 열은 10개의 샘플로 구성되었으며 CRC 환자(p1-p5)와 대조군(c1-c5)으로 구성해줬다.

```python
# Create datasets
datasets = {
    'USA': generate_data(),
    'CHN': generate_data(),
    'SA': generate_data(),
    'GER': generate_data()
}
datasets['USA']
```
```
	  p1	p2	p3	p4	p5	c1	c2	c3	c4	c5
g1	61	65	61	101	49	21	16	14	12	22
g2	25	45	58	51	44	21	26	34	17	17
g3	38	71	42	57	66	28	28	39	17	19
g4	53	55	62	58	54	24	18	16	27	14
g5	43	64	59	32	45	23	30	18	17	40
...	...	...	...	...	...	...	...	...	...	...
g96	47	31	45	31	31	29	19	44	23	27
g97	71	45	51	56	57	19	27	14	13	27
g98	40	46	48	59	54	30	21	35	40	26
g99	55	59	48	54	42	29	20	23	36	40
g100	53	37	53	39	64	25	12	20	14	18
100 rows × 10 columns
```


### 데이터 정규화, t-검정

R 패키지 RMA 알고리즘과 동일하게 정규화가 수행되었다. CRC 환자, 대조군 사이 차등 발현 유전자를 t-test로 식별했다.

```python
# Normalize the data using log transformation (simulating RMA normalization)
normalized_datasets = {k: normalize_data(v) for k, v in datasets.items()}

# Perform two-sample t-test
t_test_results = {k: perform_t_test(v) for k, v in normalized_datasets.items()}

# Apply multiple testing correction (FDR)
corrected_results = {k: apply_fdr(v) for k, v in t_test_results.items()}

corrected_results['USA']
```
```
  Gene	t_stat	p_val	q_val
0	g1	8.293142	0.000034	0.000425
1	g2	3.379468	0.009649	0.016404
2	g3	3.894228	0.004582	0.010182
3	g4	8.650333	0.000025	0.000413
4	g5	3.287387	0.011065	0.018140
...	...	...	...	...
95	g96	1.676895	0.132086	0.159140
96	g97	6.243460	0.000247	0.001025
97	g98	3.841077	0.004940	0.010611
98	g99	4.073805	0.003565	0.008897
99	g100	6.243827	0.000247	0.001025
100 rows × 4 columns
```

p-value < 0.05 and q-value(FDR) < 0.1 인 유전자만 네트워크 노드로 사용된다. 네트워크의 에지는 각 노드(유전자)의 발현량의 상관계수에 따라 가중치가 매겨져있다. 

```python
# Filter genes based on p-value and q-value thresholds
filtered_genes = {k: filter_genes(v) for k, v in corrected_results.items()}

# Create graphs for each dataset
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
	     # 데이터 생성 과정에서 그룹 내 노이즈를 너무 작게 줬는지 모든 유전자 사이에 에지가 존재해버려서, 두 노드 사이의 상관계수 절대값이 0.7을 초과하는 경우에만 에지를 추가했다
            if abs(weight) > threshold:
                graph.add_edge(gene1, gene2, weight=weight)
    return graph

graph_dic = {k: create_graph(filtered_genes[k], normalized_datasets[k]) for k, v in filtered_genes.items()}

for key, value in filtered_genes.items():
    print(f"{key}-{value.shape}")
```
```
# 국가별 네트워크 노드 개수
USA-(73, 4) 
CHN-(72, 4)
SA-(69, 4)
GER-(71, 4)
```


### 노드 강도 계산

3개의 위상 매개변수(편심도, 근접성, 매개중심성)를 계산해서 각 그래프의 노드 강도(node strength)를 계산한다.

```python
# Node strength analysis
def node_strength_analysis(graph):
    # Calculate Degree
    degree = nx.degree_centrality(graph)

    # Calculate Eccentricity
    eccentricity = nx.eccentricity(graph)

    # Calculate Closeness
    closeness = nx.closeness_centrality(graph)

    # Calculate Betweenness
    betweenness = nx.betweenness_centrality(graph)

    # Normalize values
    degree_max = max(degree.values())
    eccentricity_max = max(eccentricity.values())
    closeness_max = max(closeness.values())
    betweenness_max = max(betweenness.values())
    
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

node_strength_results['USA']
```
```
Node	Degree	Eccentricity	Closeness	Betweenness	NodeStrength
0	g1	0.888889	2	0.900000	0.029281	0.902841
1	g2	0.444444	2	0.642857	0.004125	0.497205
2	g3	0.500000	2	0.666667	0.004598	0.522846
3	g4	0.833333	2	0.857143	0.016896	0.770340
4	g5	0.375000	2	0.615385	0.001097	0.444975
...	...	...	...	...	...	...
68	g95	0.347222	2	0.605042	0.001026	0.433993
69	g97	0.638889	2	0.734694	0.005138	0.583759
70	g98	0.458333	2	0.648649	0.002088	0.485173
71	g99	0.569444	2	0.699029	0.004900	0.553126
72	g100	0.708333	2	0.774194	0.009878	0.653862
73 rows × 6 columns
```


### 에지 강도 계산
   
3가지 생물학적 특성 (피어슨 상관계수, Gene Ontology 거리, 경로 유사성 점수)에 기반하여 에지 강도를 계산한다. 

```python
# Calculate PCC
def calculate_pcc(gene1_data, gene2_data):
    return np.corrcoef(gene1_data, gene2_data)[0, 1]

# Example for one graph (USA)
usa_graph = graph_dic['USA']
normalized_data = normalized_datasets['USA']

# Initialize edge strength dictionary
edge_strengths = {}

for edge in usa_graph.edges():
    gene1, gene2 = edge
    pcc = calculate_pcc(normalized_data.loc[gene1], normalized_data.loc[gene2])
    edge_strengths[(gene1, gene2)] = {'PCC': pcc}

print("\n".join([f"{k}: {v}" for k, v in list(edge_strengths.items())[:10]]))
```
```
('g1', 'g2'): {'PCC': 0.7282491388045326, 'GO': 1.0, 'Pathway': 0.20689655172413793, 'PCC_Norm': 0.8681905318695081, 'GO_Norm': 1.0, 'Pathway_Norm': 0.5172413793103448, 'EdgeStrength': 0.7951439703932843}
('g1', 'g3'): {'PCC': 0.7653527213164605, 'GO': 1.0, 'Pathway': 0.13333333333333333, 'PCC_Norm': 0.8873800558719053, 'GO_Norm': 1.0, 'Pathway_Norm': 0.3333333333333333, 'EdgeStrength': 0.7402377964017463}
('g1', 'g4'): {'PCC': 0.8789570085087596, 'GO': 1.0, 'Pathway': 0.03125, 'PCC_Norm': 0.9461348188425974, 'GO_Norm': 1.0, 'Pathway_Norm': 0.078125, 'EdgeStrength': 0.6747532729475325}
('g1', 'g5'): {'PCC': 0.7507769564392135, 'GO': 1.0, 'Pathway': 0.0625, 'PCC_Norm': 0.8798416466124945, 'GO_Norm': 1.0, 'Pathway_Norm': 0.15625, 'EdgeStrength': 0.6786972155374982}
('g1', 'g7'): {'PCC': 0.89352960600505, 'GO': 0.7142857142857143, 'Pathway': 0.19230769230769232, 'PCC_Norm': 0.9536715899708958, 'GO_Norm': 0.4285714285714286, 'Pathway_Norm': 0.4807692307692308, 'EdgeStrength': 0.6210040831038518}
('g1', 'g10'): {'PCC': 0.7938443708009103, 'GO': 0.7714285714285715, 'Pathway': 0.19230769230769232, 'PCC_Norm': 0.9021155922981782, 'GO_Norm': 0.5428571428571429, 'Pathway_Norm': 0.4807692307692308, 'EdgeStrength': 0.6419139886415173}
('g1', 'g11'): {'PCC': 0.9128834471868577, 'GO': 0.9333333333333333, 'Pathway': 0.13636363636363635, 'PCC_Norm': 0.9636811624906784, 'GO_Norm': 0.8666666666666667, 'Pathway_Norm': 0.3409090909090909, 'EdgeStrength': 0.723752306688812}
('g1', 'g12'): {'PCC': 0.8791171328799973, 'GO': 1.0, 'Pathway': 0.13636363636363635, 'PCC_Norm': 0.9462176332302443, 'GO_Norm': 1.0, 'Pathway_Norm': 0.3409090909090909, 'EdgeStrength': 0.7623755747131118}
('g1', 'g13'): {'PCC': 0.8601812558178221, 'GO': 1.0, 'Pathway': 0.23809523809523808, 'PCC_Norm': 0.9364242266966776, 'GO_Norm': 1.0, 'Pathway_Norm': 0.5952380952380951, 'EdgeStrength': 0.8438874406449243}
('g1', 'g18'): {'PCC': 0.9459181554886529, 'GO': 1.0, 'Pathway': 0.18181818181818182, 'PCC_Norm': 0.9807663139761367, 'GO_Norm': 1.0, 'Pathway_Norm': 0.45454545454545453, 'EdgeStrength': 0.811770589507197}
```


### Gene Ontology 거리, pathway 유사성 점수 계산

각 유전자에 대한 GO term과 pathway 정보는 임의로 생성했다. Gene Ontology 거리 계산, Pathway 유사성 점수 계산 로직 또한 실제 로직이 아닌 (오직 결과만을 위한) 임시 로직을 사용했다.

```python
# Calculate Gene Ontology Distance 
def calculate_go_distance(go_annotations1, go_annotations2): #임시 로직
    go_union = go_annotations1.union(go_annotations2)
    go_intersection = go_annotations1.intersection(go_annotations2)
    go_symmetric_difference = go_annotations1.symmetric_difference(go_annotations2)
    if len(go_union) == 0:
        return 1
    return len(go_symmetric_difference) / (len(go_union) + len(go_intersection))

with open('GO-term.dic', 'r') as file:
    gene_go_annotations = json.load(file)
gene_go_annotations = {k: set(v) for k, v in gene_go_annotations.items()}
    
for edge in edge_strengths.keys():
    gene1, gene2 = edge
    go_distance = calculate_go_distance(gene_go_annotations[gene1], gene_go_annotations[gene2])
    edge_strengths[edge]['GO'] = go_distance

# Calculate Pathway Similarity 
def calculate_pathway_similarity(pathways1, pathways2): #임시 로직
    common_pathways = pathways1.intersection(pathways2)
    unique_pathways = pathways1.union(pathways2)
    if len(unique_pathways) == 0:
        return 0
    return len(common_pathways) / len(unique_pathways)

# Calculate Pathway Distance
with open('gene-pathways.dic', 'r') as file:
    gene_pathways = json.load(file)
gene_pathways = {k: set(v) for k, v in gene_pathways.items()}

for edge in edge_strengths.keys():
    gene1, gene2 = edge
    pathway_similarity = calculate_pathway_similarity(gene_pathways[gene1], gene_pathways[gene2])
    edge_strengths[edge]['Pathway'] = pathway_similarity

# Normalize and Combine Features to Calculate Edge Strength
def normalize(values):
    min_val = min(values)
    max_val = max(values)
    return [(v - min_val) / (max_val - min_val) for v in values]

# Normalize the features
pcc_values = [v['PCC'] for v in edge_strengths.values()]
go_values = [v['GO'] for v in edge_strengths.values()]
pathway_values = [v['Pathway'] for v in edge_strengths.values()]

pcc_normalized = normalize(pcc_values)
go_normalized = normalize(go_values)
pathway_normalized = normalize(pathway_values)

for i, edge in enumerate(edge_strengths.keys()):
    edge_strengths[edge]['PCC_Norm'] = pcc_normalized[i]
    edge_strengths[edge]['GO_Norm'] = go_normalized[i]
    edge_strengths[edge]['Pathway_Norm'] = pathway_normalized[i]
    edge_strengths[edge]['EdgeStrength'] = (pcc_normalized[i] + go_normalized[i] + pathway_normalized[i]) / 3
```


### 최대 클리크 식별, 클리크 연결성 점수 계산, CCP(클리크 연결성 프로파일) 식별

최대 크기의 클리크(5개)가 시드로 선택되고, 클리크 연결성 점수가 계산된다.

```python
# Identify maximum cliques
#cliques = [clique for size in range(3, 8) for clique in nx.find_cliques(usa_graph) if len(clique) == size]
cliques = sorted(nx.find_cliques(usa_graph), key=len, reverse=True)[:5]

# Calculate clique strength
#def calculate_clique_strength(clique, node_strength, edge_strengths):
#    node_strength_sum = sum(node_strength[node] for node in clique)
#    edge_strength_sum = sum(edge_strengths[edge]['EdgeStrength'] for edge in combinations(clique, 2) if edge in edge_strengths)
#    return node_strength_sum + edge_strength_sum

# Calculate clique strength
def calculate_clique_strength(clique, node_strength, edge_strengths):
    node_strength_sum = sum(node_strength[node] for node in clique)
    edge_strength_sum = 0
    for i in range(len(clique)):
        for j in range(i + 1, len(clique)):
            edge = (clique[i], clique[j])
            if edge in edge_strengths:
                edge_strength_sum += edge_strengths[edge]['EdgeStrength']
            else:
                # Also consider reverse order of edge
                edge = (clique[j], clique[i])
                if edge in edge_strengths:
                    edge_strength_sum += edge_strengths[edge]['EdgeStrength']
    return node_strength_sum + edge_strength_sum

# Node strength (previously calculated)
node_strength = node_strength_results['USA'].set_index('Node')['NodeStrength'].to_dict()

# Calculate strength for each clique
clique_strengths = {tuple(clique): calculate_clique_strength(clique, node_strength, edge_strengths) for clique in cliques}

[print(len(cliques[_])) for _ in range(5)]
print("\n".join([f"{k}: {v}" for k, v in list(edge_strengths.items())[:10]]))
```
```
30
30
29
28
27

('g1', 'g2'): {'PCC': 0.7282491388045326, 'GO': 1.0, 'Pathway': 0.20689655172413793, 'PCC_Norm': 0.8681905318695081, 'GO_Norm': 1.0, 'Pathway_Norm': 0.5172413793103448, 'EdgeStrength': 0.7951439703932843}
('g1', 'g3'): {'PCC': 0.7653527213164605, 'GO': 1.0, 'Pathway': 0.13333333333333333, 'PCC_Norm': 0.8873800558719053, 'GO_Norm': 1.0, 'Pathway_Norm': 0.3333333333333333, 'EdgeStrength': 0.7402377964017463}
('g1', 'g4'): {'PCC': 0.8789570085087596, 'GO': 1.0, 'Pathway': 0.03125, 'PCC_Norm': 0.9461348188425974, 'GO_Norm': 1.0, 'Pathway_Norm': 0.078125, 'EdgeStrength': 0.6747532729475325}
('g1', 'g5'): {'PCC': 0.7507769564392135, 'GO': 1.0, 'Pathway': 0.0625, 'PCC_Norm': 0.8798416466124945, 'GO_Norm': 1.0, 'Pathway_Norm': 0.15625, 'EdgeStrength': 0.6786972155374982}
('g1', 'g7'): {'PCC': 0.89352960600505, 'GO': 0.7142857142857143, 'Pathway': 0.19230769230769232, 'PCC_Norm': 0.9536715899708958, 'GO_Norm': 0.4285714285714286, 'Pathway_Norm': 0.4807692307692308, 'EdgeStrength': 0.6210040831038518}
('g1', 'g10'): {'PCC': 0.7938443708009103, 'GO': 0.7714285714285715, 'Pathway': 0.19230769230769232, 'PCC_Norm': 0.9021155922981782, 'GO_Norm': 0.5428571428571429, 'Pathway_Norm': 0.4807692307692308, 'EdgeStrength': 0.6419139886415173}
('g1', 'g11'): {'PCC': 0.9128834471868577, 'GO': 0.9333333333333333, 'Pathway': 0.13636363636363635, 'PCC_Norm': 0.9636811624906784, 'GO_Norm': 0.8666666666666667, 'Pathway_Norm': 0.3409090909090909, 'EdgeStrength': 0.723752306688812}
('g1', 'g12'): {'PCC': 0.8791171328799973, 'GO': 1.0, 'Pathway': 0.13636363636363635, 'PCC_Norm': 0.9462176332302443, 'GO_Norm': 1.0, 'Pathway_Norm': 0.3409090909090909, 'EdgeStrength': 0.7623755747131118}
('g1', 'g13'): {'PCC': 0.8601812558178221, 'GO': 1.0, 'Pathway': 0.23809523809523808, 'PCC_Norm': 0.9364242266966776, 'GO_Norm': 1.0, 'Pathway_Norm': 0.5952380952380951, 'EdgeStrength': 0.8438874406449243}
('g1', 'g18'): {'PCC': 0.9459181554886529, 'GO': 1.0, 'Pathway': 0.18181818181818182, 'PCC_Norm': 0.9807663139761367, 'GO_Norm': 1.0, 'Pathway_Norm': 0.45454545454545453, 'EdgeStrength': 0.811770589507197}
```

클리크 연결성 점수를 사용해서, 최대 공통 노드와 최고 클리크 강도를 기준으로 CCP가 식별된다.

```python
# Clique Connectivity Profile Algorithm (CCP)
def calculate_clique_connectivity(clique1, clique2, clique_strengths):
    common_nodes = set(clique1).intersection(clique2)
    if len(common_nodes) > 0:
        return (clique_strengths[tuple(clique1)] + clique_strengths[tuple(clique2)]) / 2
    return 0

ccp = []
for i, clique1 in enumerate(cliques):
    for j, clique2 in enumerate(cliques):
        if i != j:
            connectivity_score = calculate_clique_connectivity(clique1, clique2, clique_strengths)
            if connectivity_score > 0:
                ccp.append((clique1, clique2, connectivity_score))

ccp_df = pd.DataFrame(ccp, columns=['Clique1', 'Clique2', 'ConnectivityScore'])

ccp_df
print(len(ccp[0][0]),len(ccp[1][0])) # 0행의 클리크 크기
```
```
	Clique1							Clique2						ConnectivityScore
0	[g40, g90, g4, g94, g72, g32, g13, g69, g71, g...	[g40, g90, g4, g94, g72, g32, g13, g74, g18, g...	312.641903
1	[g40, g90, g4, g94, g72, g32, g13, g69, g71, g...	[g40, g90, g47, g97, g56, g33, g88, g1, g35, g...	302.714158
2	[g40, g90, g4, g94, g72, g32, g13, g69, g71, g...	[g40, g90, g4, g94, g72, g32, g13, g75, g88, g...	289.605115
3	[g40, g90, g4, g94, g72, g32, g13, g69, g71, g...	[g40, g90, g47, g97, g56, g33, g88, g1, g35, g...	280.067959
4	[g40, g90, g4, g94, g72, g32, g13, g74, g18, g...	[g40, g90, g4, g94, g72, g32, g13, g69, g71, g...	312.641903
5	[g40, g90, g4, g94, g72, g32, g13, g74, g18, g...	[g40, g90, g47, g97, g56, g33, g88, g1, g35, g...	303.263337
6	[g40, g90, g4, g94, g72, g32, g13, g74, g18, g...	[g40, g90, g4, g94, g72, g32, g13, g75, g88, g...	290.154293
7	[g40, g90, g4, g94, g72, g32, g13, g74, g18, g...	[g40, g90, g47, g97, g56, g33, g88, g1, g35, g...	280.617137
8	[g40, g90, g47, g97, g56, g33, g88, g1, g35, g...	[g40, g90, g4, g94, g72, g32, g13, g69, g71, g...	302.714158
9	[g40, g90, g47, g97, g56, g33, g88, g1, g35, g...	[g40, g90, g4, g94, g72, g32, g13, g74, g18, g...	303.263337
10	[g40, g90, g47, g97, g56, g33, g88, g1, g35, g...	[g40, g90, g4, g94, g72, g32, g13, g75, g88, g...	280.226548
11	[g40, g90, g47, g97, g56, g33, g88, g1, g35, g...	[g40, g90, g47, g97, g56, g33, g88, g1, g35, g...	270.689393
12	[g40, g90, g4, g94, g72, g32, g13, g75, g88, g...	[g40, g90, g4, g94, g72, g32, g13, g69, g71, g...	289.605115
13	[g40, g90, g4, g94, g72, g32, g13, g75, g88, g...	[g40, g90, g4, g94, g72, g32, g13, g74, g18, g...	290.154293
14	[g40, g90, g4, g94, g72, g32, g13, g75, g88, g...	[g40, g90, g47, g97, g56, g33, g88, g1, g35, g...	280.226548
15	[g40, g90, g4, g94, g72, g32, g13, g75, g88, g...	[g40, g90, g47, g97, g56, g33, g88, g1, g35, g...	257.580349
16	[g40, g90, g47, g97, g56, g33, g88, g1, g35, g...	[g40, g90, g4, g94, g72, g32, g13, g69, g71, g...	280.067959
17	[g40, g90, g47, g97, g56, g33, g88, g1, g35, g...	[g40, g90, g4, g94, g72, g32, g13, g74, g18, g...	280.617137
18	[g40, g90, g47, g97, g56, g33, g88, g1, g35, g...	[g40, g90, g47, g97, g56, g33, g88, g1, g35, g...	270.689393
19	[g40, g90, g47, g97, g56, g33, g88, g1, g35, g...	[g40, g90, g4, g94, g72, g32, g13, g75, g88, g...	257.580349

30 30
```

위 결과는 US 데이터로 식별한 최대 클리크들과 그들의 연결성 점수이다. 이를 국가마다 수행하고, 공통된 최대 클리크(top common scoring cliques)를 찾아서 그 클리크의 연결성 프로필(connectivity profile for one of the top common scoring cliques)을 국가마다 비교해보면, 즉 동일한 시드에서 시작하는 각 국가의 CCP를 비교하면 공통된 최대 클리크에 대한 클리크 연결성 인사이트를 얻을 수 있다. 


원문: [Cliques for the identification of gene signatures for colorectal cancer across population][1]

[1]: https://bmcsystbiol.biomedcentral.com/articles/10.1186/1752-0509-6-S3-S17
