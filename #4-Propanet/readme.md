# Outline

데이터 로드

시점별 차등 발현 유전자 식별

시점별 전사 인자 네트워크 생성

전사 인자 영향력 순위 측정

NP 알고리즘으로 시점별 주요 조절 TF 식별

# Code

관심 유전자셋이 주어졌을때, 관심(타겟) 유전자셋 TGset을 대상으로 하는 중요한 TF를 식별할 수 있다. 연구의 목적은 중요한 TF 리스트와 시점별 차등 발현 유전자 정보를 통해 시점별 TF 네트워크를 구축해서, 각 시점에서 주요 TF 및 그들의 대상 유전자 사이의 시간에 따라 변화하는 관계성을 네트워크를 통해 설명하는 것이다.

## 데이터 로드

유전자 발현량 데이터, 전사 인자 네트워크, 관심 유전자 리스트를 가져온다. 발현량 데이터의 row는 17개 유전자(g1-g10, tf1-tf7)이고, column은 10개 샘플이다. 전사 인자 네트워크는 node1이 전사인자, node2가 유전자이다. 유전자 발현량 데이터는 시계열 데이터이고 7개 시점(t1-t7)을 가진다. 따라서 row의 총 개수는 17*7=189 이다. 

```python
expression = load_expression()
tf_network = load_tf_network()
target_genes = load_target_genes()
```
```
expression
>>
	p1	p2	p3	p4	p5	p6	p7	p8	p9	p10	time
g1	12	7	8	5	12	10	12	11	11	13	t1
g2	6	11	7	9	5	4	4	6	4	7	t1
g3	7	16	6	9	4	8	7	6	6	7	t1
g4	5	17	15	16	4	6	5	5	4	4	t1
g5	6	11	11	21	6	4	4	6	6	5	t1
...	...	...	...	...	...	...	...	...	...	...	...
tf3	14	77	42	21	7	7	14	7	0	21	t7
tf4	9	6	13	11	8	8	9	9	10	9	t7
tf5	9	4	7	9	8	11	9	8	9	9	t7
tf6	13	9	6	9	13	15	12	12	13	14	t7
tf7	10	6	9	6	8	8	9	8	9	10	t7
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
['g2', 'g4', 'g6', 'g8', 'g10', 'g12', 'g14', 'g16', 'g18', 'g20']
```

## 시점별 차등 발현 유전자 식별

시계열 발현량 데이터로 시간 흐름에 따른 차등 발현 유전자를 식별한다. t2-t7 6개 시점에 대해, t1과의 발현량이 유의미하게 변화한 유전자를 식별하였다.

```python
deg_sets, de_levels = identify_degs(expression)

# counts of differentially expressed genes
for key, value in deg_sets.items():
    print(key, '-', len(value))
```
```
t2 - 18
t3 - 9
t4 - 17
t5 - 20
t6 - 24
t7 - 25
```

## 시점별 전사 인자 네트워크 생성

시점(j)별 차등 발현 유전자(DEG)와 타겟 유전자(TG)를 사용해서 시점별 초기 네트워크를 생성한다. 그래프의 노드는 TGset ∩ DEGset(j) + TF 이다. TGset과 TF는 전역적인 정보이고, DEGset(j)은 지역적인 시점별 정보이다. 

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

time_specific_networks = construct_time_specific_networks(tf_network, deg_sets, target_genes)

# Main pipeline function
def propa_net_pipeline():
    expression = load_expression()
    tf_network = load_tf_network()
    target_genes = load_target_genes()
    deg_sets, de_levels = identify_degs(expression)
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
```python
# time specific graphs
for key, value in time_specific_graphs.items():
    print(key, '\n', value)
```
```
t2 
    node1 node2
0    tf1    g2
1    tf1   g16
2    tf2    g8
3    tf2   g12
4    tf4    g2
5    tf4    g6
6    tf4   g18
7    tf4   g20
8    tf6    g2
9    tf7   g16
10   tf7   g18
t3 
   node1 node2
0   tf1    g2
1   tf2    g4
2   tf2   g12
3   tf4    g2
4   tf6    g2
t4 
   node1 node2
0   tf1    g2
1   tf1   g16
2   tf2    g4
3   tf2    g8
4   tf2   g12
5   tf3   g10
6   tf4    g2
7   tf4   g20
8   tf6    g2
9   tf7   g16
t5 
    node1 node2
0    tf1    g2
1    tf1   g16
2    tf2    g4
3    tf2    g8
4    tf2   g12
5    tf3   g10
6    tf4    g2
7    tf4    g6
8    tf4   g20
9    tf6    g2
10   tf7   g16
t6 
    node1 node2
0    tf1    g2
1    tf1   g16
2    tf2    g4
3    tf2    g8
4    tf2   g12
5    tf3   g10
6    tf4    g2
7    tf4    g6
8    tf4   g14
9    tf4   g18
10   tf4   g20
11   tf6    g2
12   tf7   g16
13   tf7   g18
t7 
    node1 node2
0    tf1    g2
1    tf1   g16
2    tf2    g4
3    tf2    g8
4    tf2   g12
5    tf3   g10
6    tf4    g2
7    tf4    g6
8    tf4   g14
9    tf4   g18
10   tf4   g20
11   tf6    g2
12   tf7   g16
13   tf7   g18
```

## 전사 인자 영향력 순위 측정

TGset ∩ DEGset 에 대한 영향력을 기준으로 전사 인자(TF)를 순위 매기기 위해 레이블 영향력 최대화(Labeled influence maximization) 알고리즘을 사용한다. 각 TF의 영향력은 시간에 따라 차별적으로 발현된 유전자들을 얼마나 많이 조절하는지에 따라 결정된다.

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

for key, value in influence_results.items():
    print(key, '-', value)
```
```
t2 - {'tf5': 0.0, 'tf6': 8.55, 'tf7': 18.2, 'tf3': 0.0, 'tf2': 21.0, 'tf4': 21.55333333333333, 'tf1': 16.195}
t3 - {'tf5': 0.0, 'tf6': 7.42, 'tf7': 0.0, 'tf3': 0.0, 'tf2': 27.49, 'tf4': 6.58, 'tf1': 7.56}
t4 - {'tf5': 0.0, 'tf6': 17.49, 'tf7': 14.31, 'tf3': 15.84, 'tf2': 27.786666666666665, 'tf4': 28.99, 'tf1': 22.62}
t5 - {'tf5': 0.0, 'tf6': 25.92, 'tf7': 25.92, 'tf3': 39.0, 'tf2': 38.50166666666666, 'tf4': 40.6, 'tf1': 35.04}
t6 - {'tf5': 0.0, 'tf6': 16.53, 'tf7': 23.19, 'tf3': 29.15, 'tf2': 47.95666666666667, 'tf4': 36.45833333333334, 'tf1': 29.535}
t7 - {'tf5': 0.0, 'tf6': 34.56, 'tf7': 31.28, 'tf3': 22.44, 'tf2': 58.01, 'tf4': 51.53783333333334, 'tf1': 41.335}
```

그래프 G에서, 각 노드의 초기 가중치는 차등 발현량으로 설정한다. TF의 초기 영향력(IL)을 0으로 설정한다. 

각 에지에 대해 확률이 1 - p인 에지를 선택하여 부분 그래프 G'를 생성한다. 여기서 p는 원본 그래프 G에서 에지의 가중치이며 초기 전사 인자 네트워크에서 가중치를 설정하지 않았기 때문에 p=0.5가 되었다. 다음으로, IL은 생성된 부분 그래프 G'에서 TF가 도달할 수 있는 노드의 가중치 평균(∑DE(s′)/|AllReachableNodesG′(t)|)에 따라 증가한다. 이 과정을 Round번 반복하고, 여러 라운드 동안의 평균을 계산하여 시점별로 IL을 업데이트함으로써 각 TF의 영향을 평가한다.

## 시점별 주요 조절 TF 식별

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

원문 링크: [PropaNet: Time-Varying Condition-Specific Transcriptional Network Construction by Network Propagation][1]


[1]: https://www.frontiersin.org/journals/plant-science/articles/10.3389/fpls.2019.00698/full
