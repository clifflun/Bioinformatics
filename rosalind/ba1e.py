with open('ba1e.txt', 'r') as f:
	genome=f.readline().strip()
	k,L,t=f.readline().strip().split()
	k=int(k)
	L=int(L)
	t=int(t)


def kmer_freq(string, k):
	kmer={}
	for i in range(len(string)-k+1):
		substring=string[i:i+k]
		kmer[substring]=kmer.get(substring,0)+1
	return kmer

def find_idx(pattern, genome):
	idx=[]
	for i in range(len(genome)-len(pattern)+1):
		if pattern == genome[i:i+len(pattern)]:
			idx.append(i)
	return idx


def clumps(genome, k, L, t):
	kmer=kmer_freq(genome, k)
	clumps=[]
	patterns=[]
	for k,v in kmer.items():
		if v >= t:
			patterns.append(k)
	for p in patterns:
		idx=find_idx(p, genome)
		for i in range(len(idx)):
			for j in range(i+1,len(idx)):
				if idx[i]-idx[j] <= L and j-i >=t:
					clumps.append(p)
	clumps = set(clumps)
	return clumps
clumps=clumps(genome, k, L, t)
[print(c, end=' ') for c in clumps]