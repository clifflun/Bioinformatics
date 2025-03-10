from collections import defaultdict
import itertools
import copy
import random

def get_db_graph():
	with open('ba3k.txt', 'r') as f:
		text=f.read().splitlines()		

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
	return debruijn_graph(text)

class CONTIGS:
	def __init__(self, graph):
		self.graph=graph
		self.in_degree=defaultdict(int)
		self.out_degree=defaultdict(int)
		self.delta_degree=defaultdict(int)
		self.stack=[]
		self.all_nodes=[]
		self.start_nodes=[]
		self.contigs=[]

	def get_degrees(self):		
		for k, values in self.graph.items():
			self.out_degree[k] = len(values)
			for v in values:
				self.in_degree[v] +=1
		all_in = set(self.in_degree.keys())
		all_out =  set(self.out_degree.keys())			
		self.all_nodes = all_in.union(all_out)
		self.delta_degree = defaultdict(int)
		for n in self.all_nodes:
			self.delta_degree[n] = self.in_degree[n]-self.out_degree[n]
		
	def get_start_nodes(self):
		self.start_nodes=[]
		for n in self.all_nodes:
			if self.out_degree[n]>0:
				if self.in_degree[n] != 1 and self.out_degree != 1:
					self.start_nodes.append(n)

	def dfs(self, key):
		next_v=self.graph[key].pop()
		self.stack.append(next_v)
		if self.in_degree[next_v] ==1 and self.out_degree[next_v]==1:
			self.dfs(next_v)
		else:
			return
			# self.stack=[]

	def get_genome_path(self):
		string=self.stack[0]
		for i in range(1,len(self.stack)):
			string+=self.stack[i][-1]
		self.contigs.append(string)
		self.stack=[]

	def get_contigs(self):
		while self.start_nodes:
			for sn in self.start_nodes: 
				self.out_degree[sn] -=1
				self.stack.append(sn)
				self.dfs(sn)
				self.get_genome_path()
				self.get_start_nodes()



def main():
	with open('ba3k_out.txt', 'r') as f:
		expected_out=f.read().splitlines()
	graph=get_db_graph()
	contig=CONTIGS(graph)
	contig.get_degrees()
	contig.get_start_nodes()
	contig.get_contigs()
	out=sorted(contig.contigs)
	res = list(set(expected_out) - set(out))
	print(res)
	print(len(res), len(expected_out), len(out))

if __name__ == '__main__':
	import sys
	sys.setrecursionlimit(15000)
	main()