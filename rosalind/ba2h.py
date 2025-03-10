with open('ba2h.txt', 'r') as f: 
	pattern=f.readline().strip()
	dna=f.readline().strip().split()



def ham_dist(s1,s2):
	ham_dist=0
	for i,j in zip(s1,s2):
		if i!=j:
			ham_dist+=1
	return ham_dist

def dbpas(pattern,dna):
	k=len(pattern)
	distance=0
	for d in dna:
		max_h_dist=float('inf')
		for i in range(len(d)-k+1):
			if max_h_dist > ham_dist(pattern,d[i:i+k]):
				max_h_dist=ham_dist(pattern,d[i:i+k])
		distance=distance+max_h_dist
	return distance


print(dbpas(pattern,dna))