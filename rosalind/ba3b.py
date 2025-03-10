with open('ba3b.txt', 'r') as f:
	dna=f.read().splitlines()

def genome_path(dna):
	init=dna[0]
	for i in range(1,len(dna)):
		init+=dna[i][-1]
	return init


print(genome_path(dna))