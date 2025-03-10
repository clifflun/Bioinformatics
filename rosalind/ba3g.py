from collections import defaultdict
import copy
import random

def load_graph():
	with open('ba3g.txt', 'r') as f:
		graph=defaultdict(list)
		lines=f.read().splitlines()
		for line in lines:
			k,v=line.split(' -> ')
			v_split=v.split(',')
			for vs in v_split:
				graph[k].append(vs)
	return graph

class EU_PATH:
	def __init__(self, graph):
		self.graph=graph
		self.in_degree=defaultdict(int)
		self.out_degree=defaultdict(int)
		self.delta_degree=defaultdict(int)
		self.has_eu_path=False
		self.start=0
		self.end=0
		self.stack=[]
		self.eu_path=[]

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
			return
		self.has_eu_path=False	


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
		print(self.delta_degree)
		print(self.start, self.end)
		print(self.has_eu_path)
		if self.has_eu_path:
			self.stack.append(self.start)
			while self.stack:
				self.dfs(self.start)
			self.eu_path=self.eu_path[::-1]

	def print_eu_path(self):
		# print(self.stack)
		# print(self.eu_path)
		print('->'.join(self.eu_path))
		
def main():
	graph=load_graph()
	eu_path=EU_PATH(graph)
	eu_path.find_eu_path()
	# eu_path.print_eu_path()


if __name__ == '__main__':
	main()