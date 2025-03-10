with open('ba1h.txt', 'r') as f:
	pattern=f.readline().strip()
	text=f.readline().strip()
	d=int(f.readline().strip())

def ham_dist(s1,s2):
	ham_dist=0
	for i,j in zip(s1,s2):
		if i!=j:
			ham_dist+=1
	return ham_dist

for i in range(len(text)-len(pattern)+1):
	substring=text[i:i+len(pattern)]
	if ham_dist(pattern, substring)<=d:
		print(i, end=' ')