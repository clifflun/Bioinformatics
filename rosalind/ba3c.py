from collections import defaultdict

with open('ba3c.txt', 'r') as f:
	dna=f.read().splitlines()

def prefix(string):
	return string[:-1]

def suffix(string):
	return string[1:]

def overlap_graph(dna):
	graph=[]
	for d in dna:
		tmp=[]
		tmp.append(d) #0
		tmp.append(prefix(d)) #1 prefix
		tmp.append(suffix(d)) #2 suffix
		graph.append(tmp)
	for i in range(len(graph)):
		for j in range(len(graph)):
			if i==j:
				continue
			# print(graph[i][1] , graph[j][2])
			if graph[i][1] == graph[j][2]:
				print(graph[j][0], '->',graph[i][0])

overlap_graph(dna)