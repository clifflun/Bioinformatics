from collections import defaultdict
import copy
import random


with open('ba3f.txt', 'r') as f:
	graph=defaultdict(list)
	lines=f.read().splitlines()
	for line in lines:
		k,v=line.split(' -> ')
		v_split=v.split(',')
		for vs in v_split:
			graph[k].append(vs)



def traverse(k, graph):
	v = graph[k]
	v_len=len(v)
	if v_len == 0:
		return 
	next_node=v.pop(random.randrange(len(v)))
	path.append(next_node)
	traverse(next_node, graph)

def main():
	exp_len = len([x for xs in list(graph.values()) for x in xs])+1
	path_len=0
	# print(exp_len)

	# print(list(graph.keys()))

	while path_len != exp_len:
		k=random.choice(list(graph.keys()))
		# print(k)
		path=[k]
		tmp_graph=copy.deepcopy(graph)
		traverse(k,tmp_graph)
		path_len=len(path)
		if path_len > 3000:
			print(path_len)
	print('->'.join(path))



def find_eu_cycle(graph, node):
	cycle = [node]
	while graph[node]:
		cycle = cycle[:1] + find_eu_cycle(
			graph,
			graph[node].pop(),
		) + cycle[1:]
	return cycle



eu_cycle = find_eu_cycle(graph, '0')
print('->'.join(eu_cycle))