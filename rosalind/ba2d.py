with open('ba2d.txt', 'r') as f: 
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
				profile[n].append(tmp[n]/len(motifs))
			except:
				profile[n].append(0)
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

def greedy(dna, k):
	best_motifs = [d[:k] for d in dna]
	motifs = [dna[0][i:i+k] for i in range(len(dna[0])-k+1)]

	for m in motifs:
		tmp_motifs=[]
		tmp_motifs.append(m)
		for i in range(1,len(dna)):
			tmp_profile=profile(tmp_motifs)
			tmp_motifs.append(profile_most_kmer(dna[i], tmp_profile, k))
		if score(tmp_motifs) < score(best_motifs):
			best_motifs = tmp_motifs
	return best_motifs

output=greedy(dna,k)
[print(o) for o in output]