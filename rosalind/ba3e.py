from collections import defaultdict

with open('ba3e.txt', 'r') as f:
	text=f.read().splitlines()
	k = len(text[0])


def kmer_composition(text,k):
	output=[]
	for i in range(len(text)-k+1):
		output.append(text[i:i+k])	
	return	output

def prefix(string):
	return string[:-1]

def suffix(string):
	return string[1:]

def debruijn_graph(dna):
	graph=defaultdict(list)
	for d in dna:
		graph[prefix(d)].append(suffix(d))
	return graph

def print_db_graph(graph):
	# print(graph)
	for k,v in graph.items():
		x=','.join(v)
		print(f'{k} -> {x}')

# dna=kmer_composition(text, k)
graph=debruijn_graph(text)
print_db_graph(graph)