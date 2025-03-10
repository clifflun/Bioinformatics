import itertools

with open('ba1k.txt', 'r') as f:
	text=f.readline().strip()
	k=f.readline().strip()
	k=int(k)

def get_freq(pattern, text):
	count=0
	for i in range(len(text)-len(pattern)+1):
		substring=text[i:i+len(pattern)]
		if substring == pattern:
			count+=1
	return count

freq=[]
it=itertools.product('ACGT', repeat=k)
while True:
	kmer=next(it, 'end')
	if kmer=='end':
		break
	string=''.join(kmer)
	freq.append(get_freq(string, text))

[print(f, end=' ') for f in freq]
