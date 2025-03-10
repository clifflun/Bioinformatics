import random
import copy

with open('ba2g.txt', 'r') as f: 
	x= f.readline().strip().split()
	k,t,N=map(int, x)
	dna=f.read().splitlines()


NUC=['A', 'C', 'G', 'T']

def score(motifs):
	score=[]
	for i in range(len(motifs[0])):
		tmp={}
		for m in motifs:
			tmp[m[i]]=tmp.get(m[i],0)+1
		max_count=max(tmp.values())
		score.append(len(motifs)-max_count)
	return	sum(score)


def profile(motifs):
	profile={'A':[], 'C':[], 'G':[], 'T':[]}
	for i in range(len(motifs[0])):
		tmp={}
		for m in motifs:
			tmp[m[i]]=tmp.get(m[i],0)+1
		for n in NUC:
			try:
				profile[n].append(tmp[n]/sum(tmp.values())+1)
			except:
				profile[n].append(1/sum(tmp.values())+1)
	return profile

def random_n(n):
	return random.randrange(n)

def gibbs_sampler(text, profile, k):
	probability={}
	for i in range(len(text)-k+1):
		substring=text[i:i+k]
		prob=1
		for j in range(len(substring)):
			# print(substring,substring[j],j)
			prob=prob*profile[substring[j]][j]
		probability[substring]=probability.get(substring,0)+prob
	_val=list(probability.values())
	_keys=list(probability.keys())
	_val_trans=[v/sum(_val) for v in _val]
	_kmer=random.choices(_keys, _val_trans)
	return ''.join(_kmer)
	


def profile_most_kmer(text, profile, k):
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
			return k

def gibbs(dna, k, N):
	init_motifs=[]
	best_motifs=[]
	for d in dna:
		random_start=random.randrange(len(dna[0])-k+1)
		init_motifs.append(d[random_start:random_start+k])
	for d in dna:
		random_start=random.randrange(len(dna[0])-k+1)
		best_motifs.append(d[random_start:random_start+k])

	for i in range(100):
		randn=random_n(len(dna))
		# print(randn, init_motifs)
		del init_motifs[randn]
		# print(init_motifs)
		tmp_profile=profile(init_motifs)
		print(tmp_profile)
		init_motifs.insert(randn, gibbs_sampler(dna[randn], tmp_profile, k))
		# print(init_motifs)
		if score(init_motifs) < score(best_motifs):
			best_motifs = init_motifs
		else:
			return best_motifs





cur_score=float('inf')
final_motifs=[]
for i in range(1000):
	tmp_motifs=gibbs(dna,k,N)
	tmp_score=score(tmp_motifs)
	if tmp_score < cur_score:
		# print(tmp_score)
		cur_score=tmp_score
		final_motifs=tmp_motifs

# print(final_motifs)
[print(f) for f in final_motifs]