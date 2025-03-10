import itertools

with open('ba2a.txt', 'r') as f:
	k,d=f.readline().strip().split()
	k=int(k)
	d=int(d)
	text=f.read().splitlines()

def ham_dist(s1,s2):
	ham_dist=0
	for i,j in zip(s1,s2):
		if i!=j:
			ham_dist+=1
	return ham_dist

pattern_count={}
it=itertools.product('ACGT', repeat=k)
while True:
	kmer=next(it, 'end')
	if kmer=='end':
		break
	string=''.join(kmer)
	for t in text:
		for i in range(len(t)-k+1):
			pattern=t[i:i+k]
			if ham_dist(string, pattern) <=d:
				pattern_count[string]=pattern_count.get(string,0)+1
				break

max_count=max(pattern_count.values())
for k,v in pattern_count.items():
	if v==max_count:
		print(k, end=' ')