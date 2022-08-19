#!/usr/bin/env python
# coding: utf-8

# In[1]:


import itertools


# In[85]:

#
# with open('LEXF.txt', 'r') as f:
#     line = f.readline()
#     n = int(f.readline())
# x= line.strip().split()
# print(x,n)


# In[88]:


def main(x,n):
    tmp=[]
    for i in itertools.product(x, repeat = n):
        tmp.append(i)
    tmp = [list(i) for i in tmp]
    tmp = [''.join(i) for i in tmp]
    out = sorted([i for i in set(tmp)])
    return out




if __name__=="__main_":
    main()
