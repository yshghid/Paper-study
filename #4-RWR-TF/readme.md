## 데이터 로드

시계열 발현량 데이터, 전사 인자 네트워크, 관심 유전자 리스트를 가져온다. 발현량 데이터의 row는 g1-g10, tf1-tf7로 17개 유전자이고, column은 10개 샘플이다. 전사 인자 네트워크는 node1이 전사인자, node2가 유전자로 구성된 데이터프레임이다. 관심 유전자는 g2, g4, g6, g8로 가정하였다. 발현량 데이터는 시계열로서 7개 시점(t1-t7)을 가져서, row의 총 개수는 17*7=189 이다. 

```python
expression = load_expression()
tf_network = load_tf_network()
target_genes = load_target_genes()
```
```
expression
>>
	p1	p2	p3	p4	p5	p6	p7	p8	p9	p10	time
g1	7	7	9	5	4	5	4	10	8	7	t1
g2	14	5	14	14	14	10	10	16	14	16	t1
g3	8	15	10	15	7	7	11	9	6	7	t1
g4	14	16	10	9	14	13	12	15	15	15	t1
g5	7	5	14	10	4	7	8	8	8	5	t1
...	...	...	...	...	...	...	...	...	...	...	...
tf3	8	8	2	17	10	7	7	10	6	7	t7
tf4	7	8	18	7	6	4	9	8	8	6	t7
tf5	8	18	4	24	4	10	5	8	8	10	t7
tf6	7	14	5	6	3	7	4	6	9	9	t7
tf7	21	6	18	18	22	19	20	20	19	22	t7
189 rows × 11 columns
```
```
tf_network
>>
	TF	Target
0	tf1	g2
1	tf1	g5
2	tf1	g16
3	tf2	g1
4	tf2	g4
5	tf2	g8
6	tf2	g11
7	tf2	g12
8	tf2	g13
9	tf3	g10
10	tf4	g1
11	tf4	g2
12	tf4	g6
13	tf4	g7
14	tf4	g14
15	tf4	g18
16	tf4	g19
17	tf4	g20
18	tf5	g9
19	tf5	g15
20	tf6	g2
21	tf6	g17
22	tf7	g3
23	tf7	g16
24	tf7	g18
```
```
target_genes
>>
['g2', 'g4', 'g6', 'g8']
```

## 시점별 차등 발현 유전자 식별

시계열 발현량 데이터로 시간 흐름에 따른 차등 발현 유전자를 식별한다. t2-t7 6개 시점에 대해, t1과의 발현량이 유의미하게 변화한 유전자를 식별하였다.

```
deg_sets = identify_degs(expression)
deg_sets
```
```
{'t2': ['g1', 'g4', 'g13', 'g16', 'g17', 'tf2', 'tf3', 'tf6'],
 't3': ['g1',
  'g2',
  'g4',
  'g6',
  'g7',
  'g9',
  'g11',
  'g17',
  'g18',
  'g19',
  'tf2',
  'tf3',
  'tf4',
  'tf5',
  'tf6'],
 't4': ['g1',
  'g2',
  'g4',
  'g6',
  'g7',
  'g8',
  'g10',
  'g11',
  'g13',
  'g16',
  'g17',
  'g18',
  'g19',
  'g20',
  'tf2',
  'tf6'],
 't5': ['g1',
  'g2',
  'g3',
  'g4',
  'g10',
  'g11',
  'g12',
  'g13',
  'g16',
  'g18',
  'g19',
  'g20',
  'tf2'],
 't6': ['g4', 'g5', 'g8', 'g10', 'g12', 'g14', 'g16', 'g18', 'g19', 'tf1'],
 't7': ['g1',
  'g2',
  'g3',
  'g4',
  'g5',
  'g10',
  'g16',
  'g18',
  'g19',
  'tf3',
  'tf6',
  'tf7']}
```

## 시점별 전사 인자 네트워크 생성

시점별 차등 발현 유전자(DEG)와 타겟 유전자(TG)를 사용해서 시점별 초기 네트워크를 생성한다. 그래프의 노드는 DEG와 TG에서 겹치는 유전자와 TF 유전자이다. 

```python
# Construct time-specific networks
def construct_time_specific_networks(tf_network, deg_sets, target_genes):
    time_specific_networks = {}
    
    for time_point, deg_set in deg_sets.items():
        # Intersected gene set
        vj = (set(deg_set) & set(target_genes)) | set(tf_network['TF'])
        
        # Construct network
        edges = []
        for _, row in tf_network.iterrows():
            if row['TF'] in vj and row['Target'] in vj:
                edges.append((row['TF'], row['Target']))
        
        time_specific_networks[time_point] = {
            'nodes': list(vj),
            'edges': edges
        }
    
    return time_specific_networks

# Main pipeline function
def propa_net_pipeline():
    expression = load_expression()
    tf_network = load_tf_network()
    target_genes = load_target_genes()
    deg_sets = identify_degs(expression)
    time_specific_networks = construct_time_specific_networks(tf_network, deg_sets, target_genes)
    
    # Populate the dictionaries
    time_specific_nodes = {}
    time_specific_graphs = {}
    
    for key, value in time_specific_networks.items():
        time_specific_nodes[key] = value['nodes']
    
        # Create DataFrame for edges
        df_edges = pd.DataFrame(value['edges'], columns=['node1', 'node2'])
        time_specific_graphs[key] = df_edges
        
    return time_specific_nodes, time_specific_graphs

time_specific_nodes, time_specific_graphs = propa_net_pipeline()
```
```
time_specific_nodes
>>
{'t2':   node1 node2
 0   tf1    g2
 1   tf2    g8
 2   tf4    g2
 3   tf6    g2,
 't3':   node1 node2
 0   tf1    g2
 1   tf2    g8
 2   tf4    g2
 3   tf6    g2,
 't4':   node1 node2
 0   tf1    g2
 1   tf4    g2
 2   tf4    g6
 3   tf6    g2,
 't5':   node1 node2
 0   tf2    g8,
 't6':   node1 node2
 0   tf4    g6,
 't7':   node1 node2
 0   tf1    g2
 1   tf4    g2
 2   tf6    g2}
```
```
time_specific_graphs
>>
{'t2':   node1 node2
 0   tf1    g2
 1   tf2    g4
 2   tf2    g8
 3   tf4    g2
 4   tf4    g6
 5   tf6    g2,
 't3':   node1 node2
 0   tf1    g2
 1   tf4    g2
 2   tf4    g6
 3   tf6    g2,
 't4':   node1 node2
 0   tf1    g2
 1   tf2    g8
 2   tf4    g2
 3   tf4    g6
 4   tf6    g2,
 't5':   node1 node2
 0   tf2    g4
 1   tf2    g8,
 't6':   node1 node2
 0   tf1    g2
 1   tf2    g4
 2   tf2    g8
 3   tf4    g2
 4   tf4    g6
 5   tf6    g2,
 't7':   node1 node2
 0   tf1    g2
 1   tf4    g2
 2   tf4    g6
 3   tf6    g2}
```

## 전사 인자 영향력 순위 측정

시점(j)별 전사 인자 초기 네트워크를 형성했다. 이때 노드로 사용된 유전자는 TGset ∩ DEGset(j) 였다. TGset ∩ DEGset에 대한 영향력을 기준으로 전사인자(TF)를 순위 매기기 위해 레이블 영향력 최대화(Labeled influence maximization) 알고리즘을 사용한다. 각 TF의 영향력은 시간에 따라 차별적으로 발현된 유전자들을 얼마나 많이 조절하는지에 따라 결정된다.

```python
# Function to calculate influence
def calculate_influence(time_specific_networks, tf_network, deg_sets, target_genes, de_levels, rounds=100):
    influence_results = {}
    
    for time_point, network in time_specific_networks.items():
        nodes = network['nodes']
        edges = network['edges']
        vj = list(set(nodes) & set(target_genes) & set(deg_sets[time_point]))
        
        if not vj:
            continue
        
        tf_set = list(set(tf_network['TF'].unique()) & set(nodes))
        
        # Initialize DE(s) and IL(t) using precomputed differential expression levels
        de = de_levels[time_point].to_dict()['differential_expression']
        il = {tf: 0 for tf in tf_set}
        
        for _ in range(rounds):
            g_prime_edges = [(u, v) for u, v in edges if np.random.rand() > 0.5]
            g_prime = {node: [] for node in nodes}
            
            for u, v in g_prime_edges:
                g_prime[u].append(v)
            
            for tf in tf_set:
                reachable = set()
                queue = [tf]
                
                while queue:
                    current = queue.pop(0)
                    for neighbor in g_prime[current]:
                        if neighbor not in reachable and neighbor not in tf_set:
                            reachable.add(neighbor)
                            queue.append(neighbor)
                
                if reachable:
                    il[tf] += sum(de[node] for node in reachable if node in de) / len(reachable)
        
        for tf in il:
            il[tf] /= rounds
        
        influence_results[time_point] = il
    
    return influence_results

influence_results = calculate_influence(time_specific_networks, tf_network, deg_sets, target_genes, de_levels)

influence_results
```
```
{'t2': {'tf7': 0.0,
  'tf4': 0.0,
  'tf6': 0.0,
  'tf2': -3.76,
  'tf3': 0.0,
  'tf1': 0.0,
  'tf5': 0.0},
 't3': {'tf7': 0.0,
  'tf4': -3.6,
  'tf6': -2.45,
  'tf2': -2.36,
  'tf3': 0.0,
  'tf1': -2.45,
  'tf5': 0.0},
 't4': {'tf7': 0.0,
  'tf4': -4.36,
  'tf6': -3.57,
  'tf2': 0.66,
  'tf3': 0.0,
  'tf1': -3.5,
  'tf5': 0.0},
 't5': {'tf7': 0.0,
  'tf4': -3.57,
  'tf6': -4.13,
  'tf2': -3.96,
  'tf3': 0.0,
  'tf1': -3.78,
  'tf5': 0.0},
 't6': {'tf7': 0.0,
  'tf4': 0.0,
  'tf6': 0.0,
  'tf2': -0.28,
  'tf3': 0.0,
  'tf1': 0.0,
  'tf5': 0.0},
 't7': {'tf7': 0.0,
  'tf4': -3.76,
  'tf6': -4.24,
  'tf2': 2.35,
  'tf3': 0.0,
  'tf1': -4.32,
  'tf5': 0.0}}
```

노드의 초기 가중치는 절대 차등 발현량으로, TF의 초기 영향력(IL)을 0으로 설정한다. 다음으로, 각 에지에 대해 확률이 1 - p인 에지를 선택하여 부분 그래프 G'를 생성한다. 여기서 p는 원본 그래프 G에서 에지의 가중치이며 초기 전사 인자 네트워크에서 가중치를 설정하지 않았기 때문에 p=0.5가 되었다. 다음으로, IL은 생성된 부분 그래프 G'에서 TF가 도달할 수 있는 노드의 가중치 평균(∑DE(s′)/|AllReachableNodesG′(t)|)에 따라 증가한다. 이 과정을 Round번 반복하고, 여러 라운드 동안의 평균을 계산하여 시점별로 IL을 업데이트함으로써 각 TF의 영향을 평가한다.

## NP 알고리즘으로 시점별 주요 조절 TF 식별

네트워크 전파를 사용하여 주요 조절 전사인자(TF)를 식별한다. 네트워크 전파 알고리즘으로는 랜덤 워크 재시작(RWR) 알고리즘이 사용되었다.

```python
def network_propagation(W, p0, alpha, iterations):
    p = p0
    for _ in range(iterations):
        p = alpha * p0 + (1 - alpha) * np.dot(W, p)
    return p

def identify_major_regulatory_tfs(time_specific_networks, tf_network, deg_sets, de_levels, alpha=0.7, iterations=100):
    major_tfs = {}
    
    for time_point, network in time_specific_networks.items():
        nodes = network['nodes']
        edges = network['edges']
        
        # Create adjacency matrix W
        node_index = {node: idx for idx, node in enumerate(nodes)}
        W = np.zeros((len(nodes), len(nodes)))
        for u, v in edges:
            W[node_index[u], node_index[v]] = 1
        
        # Normalize the adjacency matrix
        row_sums = W.sum(axis=1, keepdims=True)
        row_sums[row_sums == 0] = 1  # Avoid division by zero
        W = W / row_sums
        
        # Initialize influence scores
        tf_influence_scores = {tf: influence_results[time_point].get(tf, 0) for tf in nodes if tf.startswith('tf')}
        sorted_tfs = sorted(tf_influence_scores.keys(), key=lambda x: tf_influence_scores[x], reverse=True)
        
        # Differential expression values
        de = de_levels[time_point]['differential_expression'].to_dict()
        de_values = np.array([de[node] if node in de else 0 for node in nodes])
        
        # Network propagation and selection of major regulatory TFs
        p0 = np.zeros(len(nodes))
        selected_tfs = []
        max_scc = -1
        
        for tf in sorted_tfs:
            p0[node_index[tf]] = tf_influence_scores[tf]
            propagated_values = network_propagation(W, p0, alpha, iterations)
            
            # Check for constant input arrays
            if np.all(propagated_values == propagated_values[0]) or np.all(de_values == de_values[0]):
                scc = -1
            else:
                scc, _ = spearmanr(propagated_values, de_values)
            
            if scc > max_scc:
                max_scc = scc
                selected_tfs.append(tf)
            else:
                p0[node_index[tf]] = 0
        
        major_tfs[time_point] = selected_tfs
        
        #for tf in sorted_tfs:
        #    p0[node_index[tf]] = tf_influence_scores[tf]
        #    propagated_values = network_propagation(W, p0, alpha, iterations)
        #    scc, _ = spearmanr(propagated_values, de_values)
            
        #    if scc > max_scc:
        #        max_scc = scc
        #        selected_tfs.append(tf)
        #    else:
        #        p0[node_index[tf]] = 0
        
        #major_tfs[time_point] = selected_tfs
    
    return major_tfs

major_tfs = identify_major_regulatory_tfs(time_specific_networks, tf_network, deg_sets, de_levels)

major_tfs
```
```
{'t2': ['tf2'],
 't3': ['tf2', 'tf4'],
 't4': ['tf2', 'tf1'],
 't5': ['tf4', 'tf1', 'tf2'],
 't6': ['tf2'],
 't7': ['tf2', 'tf4']}
```
