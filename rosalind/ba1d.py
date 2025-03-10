with open('ba1d.txt', 'r') as f:
	pattern=f.readline().strip()
	genome=f.readline().strip()

def find_idx(pattern, genome):
	for i in range(len(genome)-len(pattern)+1):
		if pattern == genome[i:i+len(pattern)]:
			print(i, end=' ')

find_idx(pattern, genome)