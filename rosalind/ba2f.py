import random
import copy

with open('ba2f.txt', 'r') as f: 
	k,t= f.readline().strip().split()
	k,t=int(k),int(t)
	dna=f.read().splitlines()
# print(k,t,dna)

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
				profile[n].append(tmp[n]/len(motifs)+1)
			except:
				profile[n].append(1/len(motifs)+1)
	return profile


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

def randomized(dna, k):
	init_motifs=[]
	best_motifs=[]
	for d in dna:
		random_start=random.randrange(len(dna[0])-k+1)
		init_motifs.append(d[random_start:random_start+k])
	for d in dna:
		random_start=random.randrange(len(dna[0])-k+1)
		best_motifs.append(d[random_start:random_start+k])

	while True:
		tmp_profile=profile(init_motifs)
		motifs=[]
		for d in dna:
			motifs.append(profile_most_kmer(d,tmp_profile,k))
		if score(motifs) < score(best_motifs):
			best_motifs = motifs
			init_motifs =motifs
		else:
			return best_motifs

cur_score=float('inf')
final_motifs=[]
for i in range(1000):
	tmp_motifs=randomized(dna,k)
	tmp_score=score(tmp_motifs)
	if tmp_score < cur_score:
		# print(tmp_score)
		cur_score=tmp_score
		final_motifs=tmp_motifs

# print(final_motifs)
[print(f) for f in final_motifs]