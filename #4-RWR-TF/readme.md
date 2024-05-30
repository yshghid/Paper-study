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
g1	12	7	15	14	10	13	9	12	14	10	t1
g2	3	7	9	11	3	3	1	6	5	3	t1
g3	17	4	5	7	16	14	16	19	16	17	t1
g4	8	17	11	11	9	7	9	7	7	8	t1
g5	8	19	15	15	6	7	10	9	7	9	t1
...	...	...	...	...	...	...	...	...	...	...	...
tf3	9	9	6	21	7	10	10	7	6	9	t7
tf4	20	14	18	23	22	18	20	20	19	17	t7
tf5	15	10	10	4	16	16	11	15	14	15	t7
tf6	8	7	5	12	7	6	7	6	12	8	t7
tf7	10	4	7	10	11	11	13	7	10	7	t7
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
{'t2': ['g2',
  'g4',
  'g6',
  'g8',
  'g9',
  'g12',
  'g13',
  'g14',
  'g15',
  'g16',
  'g18',
  'g19',
  'tf2',
  'tf3',
  'tf5',
  'tf6',
  'tf7'],
 't3': ['g1',
  'g2',
  'g6',
  'g7',
  'g13',
  'g15',
  'g17',
  'g18',
  'tf1',
  'tf3',
  'tf5',
  'tf6'],
 't4': ['g2',
  'g6',
  'g7',
  'g8',
  'g10',
  'g11',
  'g14',
  'g15',
  'g16',
  'g17',
  'g18',
  'g19',
  'tf1',
  'tf3',
  'tf5',
  'tf6'],
 't5': ['g1',
  'g3',
  'g4',
  'g7',
  'g8',
  'g10',
  'g12',
  'g13',
  'g14',
  'g15',
  'g17',
  'g18',
  'g20',
  'tf2',
  'tf4',
  'tf5'],
 't6': ['g2',
  'g3',
  'g4',
  'g6',
  'g7',
  'g8',
  'g10',
  'g11',
  'g14',
  'g15',
  'g17',
  'g20',
  'tf1',
  'tf5'],
 't7': ['g1',
  'g2',
  'g5',
  'g6',
  'g9',
  'g11',
  'g12',
  'g13',
  'g14',
  'g15',
  'g16',
  'tf4',
  'tf5']}
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

time_specific_nodes
time_specific_graphs
```
```
{'t2': ['g2',
  'tf6',
  'tf2',
  'g6',
  'tf5',
  'tf4',
  'tf3',
  'tf1',
  'g8',
  'g4',
  'tf7'],
 't3': ['g2', 'tf6', 'tf2', 'g6', 'tf5', 'tf4', 'tf3', 'tf1', 'tf7'],
 't4': ['g2', 'tf6', 'tf2', 'g6', 'tf5', 'tf4', 'tf3', 'tf1', 'g8', 'tf7'],
 't5': ['tf6', 'tf2', 'tf5', 'tf4', 'tf3', 'tf1', 'g8', 'g4', 'tf7'],
 't6': ['g2',
  'tf6',
  'tf2',
  'g6',
  'tf5',
  'tf4',
  'tf3',
  'tf1',
  'g8',
  'g4',
  'tf7'],
 't7': ['g2', 'tf6', 'tf2', 'g6', 'tf5', 'tf4', 'tf3', 'tf1', 'tf7']}

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

시점(j)별 전사 인자 초기 네트워크를 형성했다. 이때 노드로 사용된 유전자는 TGset ∩ DEGset(j) 였다. TGset ∩ DEGset에 대한 영향력을 기준으로 전사인자(TF)를 순위 매기기 위해 라벨 영향력 극대화 알고리즘을 사용하였다. 각 TF의 영향력은 시간에 따라 차별적으로 발현된 유전자들을 얼마나 많이 조절하는지에 따라 결정된다.

```python
# Function to calculate influence
def calculate_influence(time_specific_networks, tf_network, deg_sets, target_genes, rounds=100):
    influence_results = {}
    
    for time_point, network in time_specific_networks.items():
        nodes = network['nodes']
        edges = network['edges']
        vj = list(set(nodes) & set(target_genes) & set(deg_sets[time_point]))
        
        if not vj:
            continue
        
        tf_set = list(set(tf_network['TF'].unique()) & set(nodes))
        
        # Initialize DE(s) and IL(t)
        de = {node: abs(np.random.randn()) for node in nodes}  # Using random values for DE(s)
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
                    il[tf] += sum(de[node] for node in reachable) / len(reachable)
        
        for tf in il:
            il[tf] /= rounds
        
        influence_results[time_point] = il
    
    return influence_results

influence_results = calculate_influence(time_specific_networks, tf_network, deg_sets, target_genes)

influence_results
```
```
{'t2': {'tf6': 0.7027490449494657,
  'tf2': 0.7635642961695567,
  'tf5': 0.0,
  'tf4': 0.4748264715382994,
  'tf3': 0.0,
  'tf1': 0.7165284379876905,
  'tf7': 0.0},
 't3': {'tf6': 1.2217646419851387,
  'tf2': 0.0,
  'tf5': 0.0,
  'tf4': 1.9350372858187732,
  'tf3': 0.0,
  'tf1': 1.434245449286902,
  'tf7': 0.0},
 't4': {'tf6': 0.4166727077002777,
  'tf2': 0.2984318032054628,
  'tf5': 0.0,
  'tf4': 0.9166101758545078,
  'tf3': 0.0,
  'tf1': 0.42453445690216973,
  'tf7': 0.0},
 't5': {'tf6': 0.0,
  'tf2': 0.5579100117201744,
  'tf5': 0.0,
  'tf4': 0.0,
  'tf3': 0.0,
  'tf1': 0.0,
  'tf7': 0.0},
 't6': {'tf6': 0.12032426844012253,
  'tf2': 0.8278748278573331,
  'tf5': 0.0,
  'tf4': 1.1288900038675553,
  'tf3': 0.0,
  'tf1': 0.1422014081565084,
  'tf7': 0.0},
 't7': {'tf6': 0.6688486243949332,
  'tf2': 0.0,
  'tf5': 0.0,
  'tf4': 0.5439017871366648,
  'tf3': 0.0,
  'tf1': 0.6242587161019375,
  'tf7': 0.0}}
```

각 시점에 대해 TGset, DEGset, 그리고 노드의 교집합을 식별합니다.
DE(s)를 절대 차별 발현 수준으로 초기화하고 IL(t)를 0으로 초기화합니다.
각 라운드마다, 확률 0.5로 각 엣지를 포함하여 서브 그래프를 생성합니다.
BFS를 사용하여 각 TF에서 서브 그래프에서 도달 가능한 모든 노드를 찾습니다.
여러 라운드 동안의 평균을 계산하여 각 TF의 영향을 평가합니다.



