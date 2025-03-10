from collections import defaultdict
import itertools
import copy
import random

def get_binary_k_db(k):
	bins=[]
	for pr in itertools.product(range(2),repeat=k):
		_str=''
		for p in pr:
			_str+=str(p)
		bins.append(_str)

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

	return debruijn_graph(bins)



class EU_PATH:
	def __init__(self, graph, k):
		self.graph=graph
		self.in_degree=defaultdict(int)
		self.out_degree=defaultdict(int)
		self.delta_degree=defaultdict(int)
		self.has_eu_path=False
		self.start=0
		self.end=0
		self.stack=[]
		self.eu_path=[]
		self.k=k

	def in_out_degree(self):		
		for k, values in self.graph.items():
			self.out_degree[k] = len(values)
			for v in values:
				self.in_degree[v] +=1
		all_in = set(self.in_degree.keys())
		all_out =  set(self.out_degree.keys())			
		all_nodes = all_in.union(all_out)
		self.delta_degree = defaultdict(int)
		for n in all_nodes:
			self.delta_degree[n] = self.in_degree[n]-self.out_degree[n]
		return self.in_degree, self.out_degree, self.delta_degree

	def verify_eu_path_exist(self):
		pos_count = 0
		neg_count = 0 
		for k,v in self.delta_degree.items():
			if v == 1:
				pos_count+=1
				self.end = k
			if v == -1: 
				neg_count+=1
				self.start = k
		if pos_count == 1 and neg_count == 1:
			self.has_eu_path=True


	def dfs(self, key):
		if self.out_degree[key] > 0:
			_next=self.graph[key].pop()	
			self.stack.append(_next)
			self.out_degree[key]-=1
			# print(key, self.out_degree[key], self.graph)
			self.dfs(_next)
		if self.out_degree[key] == 0:
			self.backtrack()

	def backtrack(self):
		_tmp = self.stack.pop()
		# print(_tmp)
		if self.out_degree[_tmp] == 0:
			self.eu_path.append(_tmp)
		if self.out_degree[_tmp] > 0:
			self.dfs(_tmp)

	def find_eu_path(self):
		self.in_out_degree()
		self.verify_eu_path_exist()
		if self.has_eu_path:
			self.stack.append(self.start)
			while self.stack:
				self.dfs(self.start)
			self.eu_path=self.eu_path[::-1]
		else: 
			print('no eu path')

	def find_eu_cycle(self):
		self.in_out_degree()
		self.verify_eu_path_exist()
		self.stack.append(list(self.graph.keys())[0])
		while self.stack:
			self.dfs(list(self.graph.keys())[0])
		self.eu_path=self.eu_path[::-1]

	def print_eu_path(self):
		# print(self.stack)
		# print(self.eu_path)
		print('->'.join(self.eu_path))

	def print_genome_cycle(self):
		if not self.stack:
			init=self.eu_path[0]
			for i in range(1,len(self.eu_path)):
				init+=self.eu_path[i][-1]
			print('0'+init[:-self.k])

	def print_genome_path(self):
		if not self.stack:
			init=self.eu_path[0]
			for i in range(1,len(self.eu_path)):
				init+=self.eu_path[i][-1]
			print(init)
		
def main():
	k=10
	graph=get_binary_k_db(k)
	eu_path=EU_PATH(graph, k)
	eu_path.find_eu_cycle()
	# eu_path.print_eu_path()
	eu_path.print_genome_cycle()



if __name__ == '__main__':
	import sys
	sys.setrecursionlimit(15000)
	main()