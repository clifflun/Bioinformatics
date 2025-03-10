with open('ba10b.txt', 'r') as f:
	string=f.readline().strip()
	f.readline()
	alpha=f.readline().strip().split()
	f.readline()
	hpath=f.readline().strip()
	f.readline()
	states=f.readline().strip().split()


trans = {'Ax':0.242,
		'Ay':0.477,
		'Az':0.281,
		'Bx':0.595,
		'By':0.274,
		'Bz':0.131}


final_prob=1
for x,y in zip(hpath, string):
	sub_path=x+y
	print(sub_path)
	final_prob = final_prob * trans[sub_path]
print(final_prob)