###** Adrian A. Davin 2021 **##

import sys
import os
import ete3

def root_tree_in_position(tree_file, sp1, sp2):

    with open(tree_file) as f:
        t = ete3.Tree(f.readline().strip(), format=1)
        
    if len(t.get_children()) != 2:
        print("Tree is not rooted. Please, use this tree as input and run this script again")
        t.set_outgroup(t.get_leaves()[0])
        print(t.write(format=1))
        return None
    
    try: 
        nsp1 = t&sp1
    except:
        print("Could not find %s in tree" % sp1)
        return None
    try: 
        nsp2 = t&sp2
    except:
        print("Could not find %s in tree" % sp2)
        return None
    
    ca = t.get_common_ancestor(nsp1, nsp2)
    try:
        t.set_outgroup(ca)
    except:
        print("Cannot root there! Is it already rooted in that branch?")
        return None
    
    print(t.write(format=1))

if __name__ == "__main__":

    if len(sys.argv) != 4:
        print("usage: python root_tree_in_position.py tree_file sp1 sp2")
        print("It will reroot the tree in the branch leading to the common ancestor of sp1 and sp2")
        print("Use a rooted tree in the first place!")
        print("It might not work with non-binary trees")
        print("sp1 and sp2 must be the names of two leaves in the species tree")
        exit(0)
        

    scr, tree_file, sp1, sp2 = sys.argv
    root_tree_in_position(tree_file, sp1, sp2)
