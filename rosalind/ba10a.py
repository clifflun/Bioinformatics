prob = {'AA':0.479,
		'AB':0.521,
		'BA':0.515,
		'BB':0.485}

path = 'ABBBAABAAABBABBBBABAABABAABBABABBABBBBBAABAABAAABB'

final_prob=1
for i in range(len(path)-1):
	sub_path=path[i:i+2]
	print(sub_path)
	final_prob = final_prob * prob[sub_path]
print(final_prob*0.5)