# Outline

ppi 네트워크 로드
코어성을 활용한 ppi 네트워크 정규화
RWR 실행, 유의미하게 재순위화된 노드 식별
파라미터 α 선택
5-fold 교차 검증을 통한 4가지 정규화 방법 비교
유전자 모듈 식별

# Code

### ppi 네트워크 로드

분석에 사용할 ppi 네트워크를 가져온다. 해당 ppi 네트워크는 20개의 노드(g1-g20)와 44개의 에지를 갖는다.

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

각 노드의 코어 숫자가 계산된다. k-코어는 모든 노드가 최소 k개의 이웃을 갖는, 즉 모든 노드의 차수가 k 이상인 최대 부분 그래프이다. 코어 숫자는 그 노드가 속하는 가장 큰 k-코어에서의 k 값이다. 

차수는 단순히 연결된 이웃의 수를 의미하지만, 코어 숫자는 그 노드가 그래프의 중심부에 얼마나 가까운지를 나타내는 더 복잡한 지표이다. 코어 숫자는 그래프 전체의 구조와 각 노드의 위치를 고려하여 노드가 그래프의 중심부에 얼마나 가까운지를 나타낸다.

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

### 코어성을 활용한 ppi 네트워크 정규화

ppi 네트워크의 인접 행렬을 생성하고 ppi 네트워크를 정규화한다. 정규화 방식으로는 일반적인 차수를 사용한 방식과, 논문에서 제시하는 3가지의 코어성을 사용하는 방식이 있다.

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

예시로, 코어 숫자를 사용한 정규화 결과를 나타내면 다음과 같다.

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

### RWR 실행, 유의미하게 재순위화된 노드 식별

정규화한 ppi 네트워크에서 RWR 알고리즘을 실행한다. 질병 연관 유전자를 시드 노드로 선택하고, 노드의 초기 가중치를 시드 노드는 1, 나머지 노드는 0으로 설정한다. 초기 가중치(p0)와 최종 전파 가중치(pk)를 비교하여 재순위화된 노드를 식별할 수 있다.

재순위화 결과에 유의성 수준을 할당하기 위해 무작위화된 네트워크를 사용하여 얻은 전파 가중치와 비교한다. 네트워크 무작위화에는 double-edge swap 알고리즘을 사용하여 n=100개의 무작위 네트워크를 생성하고 이를 pk와 비교하여 각 노드에 유의성 수준을 할당하였다. 

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

예시로, α = 0.1 에 대해 RWR를 수행하고 유의미하게 재순위화된 노드를 식별하면 다음과 같다. 논문에서는 p-value<0.0099 를 적용하였으나, 코드에서는 <0.3을 적용하였다.

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

### 파라미터 α 선택

α = 0.3, 0.5, 0.8 에 대해 해당 알고리즘을 수행하였다. RWR 에서 α는 재시작 파라미터로, 값이 1에 가까울수록 해당 워크가 시드 노드로 돌아갈 확률이 높다. 즉, 시드 노드의 가중치는 1에 가깝게 유지되며, 시드 노드의 영향력이 시드에 가까운 노드에만 전파되게 된다. 

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

α를 조절해서 새로운 질병 관련 유전자를 찾는 것(finding novel disease-associated genes)과 잠재적인 거짓 예측을 포함시키는 것(including potential false predictions) 사이의 균형을 조절할 수 있다. 본 논문에서는 α = 0.8을 적용하여 false positive 를 줄이는 것을 선택하였다.

### 5-fold 교차 검증을 통한 4가지 정규화 방법 비교

4가지 정규화 방법의 질병 연관 유전자를 식별하는 성능을 비교하기 위해, 질병 유전자를 1:4 로 나누어 훈련 세트를 생성하고, 훈련 세트가 시드 유전자로 설정된 네 번의 네트워크 전파를 실행하였다. 질병 유전자를 5개로 정했기 때문에 4번의 네트워크 전파 각각에서 4개 유전자가 훈련 세트가 되었다. α는 0.8 로 설정되었다. 

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

분석 결과 노드의 차수와 코어 숫자 간의 차이를 적용한 'diff' 정규화가 가장 높은 AUROC를 보인다. 본 논문에서는 'core' 정규화가 가장 높은 AUROC를 보였다.

### 유전자 모듈 식별

먼저 시드와 에지를 갖는 노드들로 구성된 초기 시드 서브 네트워크를 생성한다. 네트워크 전파 수행 후 특정 조건(pvalue<0.01, w>wmin)을 만족하는 노드 중 시드 서브 네트워크와 에지를 갖는 노드가 있다면 추가하여, 최종 시드 서브 네트워크 내에 각 모듈이 생성된다.

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

코드에서는 p-value<0.3, wmin=0.75로 설정하여 6개의 유전자가 포함된 모듈이 생성되었다. 모듈에는 질병 유전자 g2, g4, g6, g8, g10이 포함되어 있으며, g3이 추가되었다. 

참고로, 초기 시드 서브 네트워크는 아래와 같았다. g6은 다른 질병 유전자와 노드를 갖지 않았으나, 네트워크 전파 후 유의미하게 가중치가 변화한 것으로 식별된 g3이 추가됨으로써 나머지 질병 유전자와 하나의 모듈을 형성하게 되었다.

```
Module 1: {'g4', 'g10', 'g8', 'g2'}
Module 2: {'g6'}
```


원문: [NetCore: a network propagation approach using node coreness][1]

[1]: https://academic.oup.com/nar/article/48/17/e98/5879427?login=false
