from collections import defaultdict

with open('ba3d.txt', 'r') as f:
	k=int(f.readline().strip())
	text=f.readline().strip()

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

dna=kmer_composition(text, k)
graph=debruijn_graph(dna)
print_db_graph(graph)