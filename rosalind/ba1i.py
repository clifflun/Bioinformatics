import itertools

with open('ba1i.txt', 'r') as f:
	text=f.readline().strip()
	k,d=f.readline().strip().split()
	k=int(k)
	d=int(d)

def ham_dist(s1,s2):
	ham_dist=0
	for i,j in zip(s1,s2):
		if i!=j:
			ham_dist+=1
	return ham_dist


count={}
it=itertools.product('ATCG', repeat=k)
while True:
	kmer=next(it, 'end')
	if kmer=='end':
		break
	string=''.join(kmer)
	for i in range(len(text)-k+1):
		if ham_dist(string,text[i:i+k])<=d:
			count[string]=count.get(string,0)+1

max_count=max(count.values())
for k,v in count.items():
	if v==max_count:
		print(k, end=' ')