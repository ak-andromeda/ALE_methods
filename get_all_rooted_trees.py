import sys
import os
import ete3
import copy


def get_all_rooted_trees(tree_file):

    with open(tree_file) as f:
        t = ete3.Tree(f.readline().strip(), format=1)   
        
    # First I construct a dict that allows me to map easily the nodes later
      
    n2sps = dict()
    
    # Second, I root randomly the tree (this would prevent an error later if the tree is unrooted in the first place)
    
    t.set_outgroup(t.get_leaves()[0])
    
    for n in t.traverse():
        if n.is_leaf():
            continue
        else:
            
            c1, c2 = n.get_children()
            c1ln = (c1.get_leaves()[0]).name
            c2ln = (c2.get_leaves()[0]).name
        n2sps[n] = (c1ln, c2ln)
            
    for n in t.traverse():          
        if n.is_root():
            print(t.write(format=1))        
        elif (n.up).is_root():
            continue
        else:   
            mt = copy.deepcopy(t)
            
            if not n.is_leaf():                
                mn = mt.get_common_ancestor(n2sps[n])
                mt.set_outgroup(mn)                
            else:
                mt.set_outgroup(mt&(n.name))                
            print(mt.write(format=1))

    #print(t.write(format=1))

if __name__ == "__main__":

    if len(sys.argv) != 2:
        print("usage: python get_all_rooted_trees.py tree_file")
        print("It will reroot the tree in each branch and output to the stdout")        
        exit(0)
        
    scr, tree_file = sys.argv
    get_all_rooted_trees(tree_file)
