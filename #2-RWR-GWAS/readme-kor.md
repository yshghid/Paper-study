# Outline

1) GWAS 카탈로그 로드, SNP 통계 테스트
2) 시드 유전자, 타겟 유전자 세팅
3) 그래프 로드
4) 초기 가중치 설정
5) RWR 알고리즘 수행

# Pipeline

## GWAS 카탈로그 로드, SNP 통계 테스트

질병 연관 SNP에 대한 정보는 GWAS 카탈로그로 다운로드할 수 있다. 카탈로그에는 각 SNP의 염색체 위치, 연관된 유전자, 통계적 유의성 (p-value), accession id 등이 제공된다. 해당 데이터에서 SNP는 1283개이고 유전자는 20개(g1-g20)이다.

```python
# Load gwas catalog
gwas_catalog = load_catalog()
gwas_catalog
```
```
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

질병 연관 SNP 정보를 유전자 단위 점수로 변환하는 여러 툴이 있다. 툴들은 통계 테스트 방법이 조금씩 다르며, 해당 논문에서는 PEGASUS를 사용하였다. PEGASUS는 GWAS 카탈로그를 입력으로 받아, 질병 연관 SNP 정보를 활용하여 유전자마다 질병과의 연관성에 대한 통계적 유의성(Pvalue)과 점수(Score)를 생성해준다. 아래 PEGASUS 데이터는 g1-g20에 대한 질병 연관성 분석 결과이다.

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

## 시드 유전자, 타겟 유전자 세팅

PEGASUS 데이터 상에서 p-value<0.05 인 유의미한 유전자들은 시드 유전자로 선택된다. GWAS 카탈로그 상에서 질병과의 연관성이 언급된 SNP와 연관된 유전자들은 타겟 유전자로 선택된다. PEGASUS 데이터는 유전자 20개 중 15개의 p-value가 <0.05 이므로 시드 유전자는 15개, 타겟 유전자는 20개가 된다.

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
시드와 타겟 유전자 이름이 정수인 이유는, 네트워크 전파에는 노드에 대한 명명법으로 gene symbol이 아닌 integer 형태인 ncbi id를 주로 사용하기 때문이다. gene symbol id g1-g20 각각에 대한 ncbi id는 g를 뺀 뒤의 숫자로 할당되었다. 

## 그래프 로드

ppi 네트워크는 소스 노드와 타겟 노드 정보로 구성돼있다. 아래 네트워크는 20개의 노드에 대해 80개의 에지를 갖는다. 본 논문에서 에지의 가중치는 고려되지 않았다.

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

## 초기 가중치 설정

로드한 네트워크에서 시드의 초기 가중치를 설정한다. 초기 가중치는 PEGASUS 데이터에서 얻은 Score 를 사용한다. 시드는 PEGASUS 데이터에서 얻은 15개 유전자였으며, 시드 유전자가 아닌 경우 초기 가중치는 0으로 설정된다.

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

## RWR 알고리즘 수행

랜덤 워크 재시작(RWR) 알고리즘으로는 pagerank 함수를 사용했다. PageRank는 웹 페이지의 중요성을 계산하기 위해 Google에서 개발한 알고리즘으로, 각 노드(예: 웹 페이지)는 연결된 다른 노드들로부터 가중치를 받는다. pagerank에서는 모든 노드에 초기 가중치를 균일하게 분배하고 모든 노드에서 시작하지만, RWR은 시드 노드에 특정 가중치를 할당하고 시드 노드에서 시작한다.

pagerank의 α는 임의의 노드로 이동할 확률을 나타낸다. 값이 클수록 노드 간의 전환을 허용한다. RWR의 α는 restart probability로서 시드 노드로 돌아올 확률을 나타내며, 시드 노드로부터의 관계를 강조한다. 본 논문에서는 α = 0.1로 RWR을 수행하면 본질적으로 질병 유전자 즉 시드 노드의 한-홉 이웃(one-hop neighbourhood)에만 전파되었고 α = 0.3은 두- 혹은 세-홉 이웃까지 영향이 전파되는 것으로 언급되어 있으며, α = 0.3 으로 수행하였다.

```
def perform_rwr_nx(alpha, graph, seeds):   
    rwr_scores = pagerank(graph, alpha=alpha, personalization=seeds)
    return rwr_scores

alpha = 0.3
rwr_scores = perform_rwr_nx(alpha, graph, pagerank_seeds)
rwr_scores
```
```
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

결과를 정리하면 다음과 같다.

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

초기 가중치가 0으로 할당되어 순위가 20위였던 유전자 g1, g18의 순위가 네트워크 전파 후 각각 14위, 15위로 재랭크된 것을 확인할 수 있다. 

원문: [Network propagation for GWAS analysis: a practical guide to leveraging molecular networks for disease gene discovery][1]

[1]: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10858647/
