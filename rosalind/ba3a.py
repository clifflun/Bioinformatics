with open('ba3a.txt', 'r') as f:
	k=int(f.readline().strip())
	text=f.readline().strip()

def kmer_composition(text): -> list
	output=[]
	for i in range(len(text)-k+1):
		output.append(text[i:i+k])	
	return	output