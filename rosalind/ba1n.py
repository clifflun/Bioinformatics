import itertools

with open('ba1n.txt', 'r') as f:
	text=f.readline().strip()
	d=int(f.readline().strip())

def ham_dist(s1,s2):
	ham_dist=0
	for i,j in zip(s1,s2):
		if i!=j:
			ham_dist+=1
	return ham_dist

it=itertools.product('ACGT', repeat=len(text))
while True:
	kmer=next(it, 'end')
	if kmer=='end':
		break
	string=''.join(kmer)
	if ham_dist(string, text) <=d:
		print(string)