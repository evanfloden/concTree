#! /usr/bin/env python
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning) 

import sys
sys.path.append("/software/cn/el6.3/python-3.0.0/site-packages/ete2-2.0rev111/")

from ete2 import Tree

t=Tree(sys.argv[1])
sup_sum=0
count=0
for node in t.traverse("postorder"):
  if not node.is_leaf():
    sup_sum+=node.support;
    count+=1;

print '%.4f'%(sup_sum/count);
