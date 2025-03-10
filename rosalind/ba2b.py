import itertools

with open('ba2b.txt', 'r') as f:
	k=f.readline().strip()
	k=int(k)
	text=f.read().splitlines()

def ham_dist(s1,s2):
	ham_dist=0
	for i,j in zip(s1,s2):
		if i!=j:
			ham_dist+=1
	return ham_dist

def min_ham_dist(pattern, text):
	cur_min=99999999
	for i in range(len(text)-len(pattern)+1):
		score=ham_dist(text[i:i+len(pattern)], pattern)
		if score<cur_min:
			cur_min=score
	return cur_min



count={}
it=itertools.product('ACGT', repeat=k)
while True:
	kmer=next(it, 'end')
	if kmer=='end':
		break
	string=''.join(kmer)
	for t in text:
		count[string]=count.get(string,0)+min_ham_dist(string, t)

print(count)

min_count=min(count.values())
for k,v in count.items():
	if v==min_count:
		print(k, end=' ')