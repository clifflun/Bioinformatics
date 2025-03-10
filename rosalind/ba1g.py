with open('ba1g.txt', 'r') as f:
	s1=f.readline().strip()
	s2=f.readline().strip()

dist=0
for a,b in zip(s1, s2):
	if a!=b: 
		dist+=1
print(dist)