# Outline

PPI 네트워크 로드, 코어 숫자 계산

# Code

## ppi 네트워크 로드, 인접 행렬 생성

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

## 코어성을 활용한 ppi 네트워크 정규화

ppi 네트워크의 인접 행렬을 생성하고 ppi 네트워크를 정규화한다. 정규화 방식으로는 일반적인 차수를 사용한 방식과, 논문에서 제시하는 3가지의 코어성을 사용하는 방식이 있다.

```
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

### RWR 알고리즘 수행, 유의미하게 재순위화된 노드 식별

네트워크 전파 기법으로는 랜덤 워크 재시작 알고리즘이 사용되었다. 최종 전파 가중치와 네트워크 무작위화에서 유도된 P-값을 사용하여 유의미하게 재순위화된 노드를 식별한다. 


NetCore의 재순위화 결과에 유의성 수준을 할당하기 위해 우리는 무작위 차수-보존 네트워크(RDPN)을 사용하여 정규화를 적용했습니다. 이 방법은 입력 네트워크의 무작위화를 기반으로 하여 각 노드의 전파 가중치가 무작위 차수-보존 네트워크를 사용하여 얻은 전파 가중치와 비교되도록 합니다. 이러한 네트워크를 생성하기 위해 우리는 Networkx의 double-edge swap 알고리즘을 사용했습니다. 알고리즘은 각 단계에서 무작위로 두 개의 엣지 (u, v)와 (x, y)를 선택하여 제거한 다음, 새로운 엣지 (u, y)와 (x, v)를 생성하며 최대 n개의 무작위 교환을 실행할 수 있도록 합니다. 스왑은 네트워크가 연결된 상태를 유지할 때만 유지됩니다. 이렇게 n개의 무작위 네트워크가 생성된 후, 이러한 무작위 네트워크로 얻어진 전파 가중치를 사용하여 유의성 수준을 계산합니다. 따라서, 노드 v의 전파 가중치 w(v)에 대한 P-값 pv는 n개의 무작위 네트워크와 해당 전파 가중치 w1(v), ..., wn(v)로 정의됩니다. 분석에서 우리는 n = 100개의 무작위 네트워크를 생성하여 달성할 수 있는 최소 P-값이 P = 0.0099가 되도록 했습니다. 이 임계값은 기본값으로 선택되지만 다른 설정에서는 사용자가 변경할 수 있습니다.

