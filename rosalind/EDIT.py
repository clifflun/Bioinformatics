#!/usr/bin/env python
# coding: utf-8

# In[27]:


import numpy as np


# Dynammic programming
#
# Very similar idea to longest common subsequence
#
# Any indel or substitution = +1
#
# Wagner-Fischer algo
#

# In[78]:


# ##getting fasta in list format
# fa_id = []
# temp = ''
# for line in open('EDIT.txt'):
#     if line.startswith('>'):
#         fa_id.append(line.strip()[1:])
#         temp +='>'
#     else:
#         temp += line.strip()
# fa_seq = temp.split('>')[1:]
#
# print(fa_id)
# # print(fa_seq)


# # In[79]:
#
#
# s1 = fa_seq[0]
# s2 = fa_seq[1]


# In[88]:


def wf(s1,s2):
    m = np.zeros((len(s1)+1, len(s2)+1))
    m[:,0] = range(len(s1)+1)
    m[0,:] = range(len(s2)+1)
    for i in range(len(s1)):
        for j in range(len(s2)):
            if s1[i] == s2 [j]:
                sub = 0
            else:
                sub = 1
            m[i+1,j+1] = min(m[i,j+1]+1, m[i+1,j]+1, m[i,j]+sub)
    return(m)


# # In[89]:
#
#
# get_ipython().run_cell_magic('time', '', 'wf(s1, s2)')
#
#
# # In[ ]:
#
