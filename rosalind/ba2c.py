import itertools

with open('ba2c.txt', 'r') as f:
	text=f.readline().strip()
	k=int(f.readline().strip())
	profile={}
	profile['A']=[float(a) for a in f.readline().strip().split(' ')]
	profile['C']=[float(a) for a in f.readline().strip().split(' ')]
	profile['G']=[float(a) for a in f.readline().strip().split(' ')]
	profile['T']=[float(a) for a in f.readline().strip().split(' ')]
# print(profile)

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

probability={}
for i in range(len(text)-k+1):
	substring=text[i:i+k]
	prob=1
	for j in range(len(substring)):
		# print(substring,substring[j],j)
		prob=prob*profile[substring[j]][j]
	probability[substring]=probability.get(substring,0)+prob

max_count=max(probability.values())
for k,v in probability.items():
	if v==max_count:
		print(k, end=' ')