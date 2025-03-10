import itertools

print('starting')
with open('ba1l.txt', 'r') as f:
	text=f.readline().strip()

# i=0
# it=itertools.product('ACGT', repeat=len(text))
# while True:
# 	kmer=next(it, 'end')
# 	if kmer=='end':
# 		break
# 	string=''.join(kmer)
# 	if string == text:
# 		print(string, i)
# 		break	
# 	i+=1
# 	print(i)
	

def PatternToNumber(pattern):
    seq = {'A':0, 'C':1, 'G':2, 'T':3}
    k = len(pattern)
    num = 0
    for i in range(k):
        num += seq[pattern[i]] * pow(4, k-i-1)
    return num

print(PatternToNumber(text))