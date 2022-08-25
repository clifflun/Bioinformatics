#!/usr/bin/env python
# coding: utf-8

# This is basically the most bare bones dynamic programming solution.
#
# 1. Create scoring matrix accroding to the rules
# 2. Backtrack the sequence
#
# Rules:
#
# if match: +1 to top-left diagonal
# if not match : max(top, left)
#
# Backtrack:
#
# find the matches == find the cell that was +1 from diagonal
#
# The property the matched cell is that ***BOTH*** top and left cells have to be -1 of current cell

# In[1]:


import numpy as np


# In[19]:

#
# ##getting fasta in list format
# fa_id = []
# temp = ''
# for line in open('LCSQ.txt'):
#     if line.startswith('>'):
#         fa_id.append(line.strip()[1:])
#         temp +='>'
#     else:
#         temp += line.strip()
# fa_seq = temp.split('>')[1:]
#
# print(fa_id)
# print(fa_seq)
#

# # In[20]:
#
#
# s1 = fa_seq[0]
# s2 = fa_seq[1]
# 
#
# # In[23]:


def backtrack(s1, s2):
    ## generating scoring matrix
    m = np.zeros((len(s1)+1, len(s2)+1))
    for i in range(len(s1)):
        for j in range(len(s2)):
            if s1[i] == s2[j]:
                m[i+1,j+1] = m[i,j] +1 # +1 from top-left, matched
            else:
                m[i+1,j+1] = max(m[i,j+1], m[i+1, j]) # max(top, left), not matched

    ## backtracking

    ## finding starting index
    start = np.where(m == m.max()) #start from highest score
    i = start[0][-1]
    j = start[1][-1]


    lcs = ''
    while i!=0 and j!=0:
        if m[i-1,j] == m[i,j]: #move to the next cell
            i -= 1
        elif m[i,j-1] == m[i,j]:
            j -= 1
        elif m[i-1,j] == m[i,j-1]: # apppend nucleotide if found diagonal cell
            lcs += s1[i-1]
            i -= 1
            j -= 1
    return(lcs[::-1])


# In[22]:




# In[ ]:
