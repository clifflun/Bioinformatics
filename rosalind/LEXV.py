#!/usr/bin/env python
# coding: utf-8

# In[28]:


import itertools
import LEXF


# In[63]:


with open('LEXV.txt', 'r') as f: 
    x = f.readline().split()
    n = int(f.readline())

# x = ['D', 'N', 'A']
# n = 3


# In[64]:


def main(x,n):
    master = []
    for i in range(1, n+1):
        master.append(LEXF.main(x,i))
    out = [i for j in master for i in j]
    order = {key:i for i, key in enumerate(x)}
    out = sorted(out, key=lambda a: [order.get(c,ord(c)) for c in a])
    return out


# In[65]:


if __name__=='__main__':
    sort = main(x,n)
    for i in sort:
        print(i)


# In[48]:





# In[53]:





# In[ ]:




