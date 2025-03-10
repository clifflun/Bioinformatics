from collections import defaultdict
import itertools
import copy
import random

def get_paired_db_graph():
	with open('ba3j.txt', 'r') as f:
		k,d=f.readline().strip().split(' ')
		text=f.read().splitlines()
		text_pair=[tuple(t.split('|')) for t in text]

	def prefix(string):
		return string[:-1]

	def suffix(string):
		return string[1:]

	def paried_debruijn_graph(dna):
		graph=defaultdict(list)
		for d in dna:
			# print(d)
			graph[prefix(d[0])+'|'+prefix(d[1])].append(suffix(d[0])+'|'+suffix(d[1]))
		return graph
	p_graph=paried_debruijn_graph(text_pair)

	return k,d,p_graph

class EU_PATH:
	def __init__(self, graph, k, d):
		self.graph=graph
		self.in_degree=defaultdict(int)
		self.out_degree=defaultdict(int)
		self.delta_degree=defaultdict(int)
		self.has_eu_path=False
		self.start=0
		self.end=0
		self.stack=[]
		self.eu_path=[]
		self.k=int(k)
		self.d=int(d)

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
		# print(key, self.out_degree[key], self.graph)	
		if self.out_degree[key] > 0:
			_next=self.graph[key].pop()	
			self.stack.append(_next)
			self.out_degree[key]-=1
			self.dfs(_next)
		if self.out_degree[key] == 0:
			self.backtrack()

	def backtrack(self):
		_tmp = self.stack.pop()
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

	def print_genome_path_paired(self):
		if not self.stack:
			string=self.eu_path[0][:self.k-1]
			for i in range(1,len(self.eu_path)):
				string+=self.eu_path[i].split('|')[0][-1]
			for i in range((len(self.eu_path)-self.k-self.d), len(self.eu_path)):
				string+=self.eu_path[i].split('|')[1][-1]
			return string

def main():
	k, d, graph=get_paired_db_graph()
	# print(k,d,graph)
	eu_path=EU_PATH(graph, k, d)
	eu_path.find_eu_path()
	# eu_path.print_eu_path()
	out=eu_path.print_genome_path_paired()
	print(out)


if __name__ == '__main__':
	import sys
	sys.setrecursionlimit(15000)
	main()