with open('ba1f.txt', 'r') as f:
	genome=f.readline().strip()


def skew(genome):
	skew=[0]*(len(genome)+1)	
	for i in range(len(genome)):
		if genome[i] == 'C':
			skew[i+1]=skew[i]-1
		elif genome[i] == 'G':
			skew[i+1]=skew[i]+1
		else: 
			skew[i+1]=skew[i]
	return skew

skew=skew(genome)
min_skew=min(skew)
for i in range(len(skew)):
	if min_skew==skew[i]:
		print(i, end=' ')