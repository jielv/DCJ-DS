#!/users/jl38/import/Python-2.6/python

#argv[1]:number of DCJ moves that will perform
#argv[2]:persentage of marked genes
#argv[3]:directory of starting genome, each chromosome is stored in one file as a string.

import sys
import math
import random
import re
import time
import copy
import shutil
from collections import deque

class Gene:
    """
    define a Gene class, for each Gene, storing 10 kinds of informations.
    """
    def __init__(self,id,pa,hc,tc,flip,height,size,circle,mark,marksignal):
        self.Id=id
        self.Parent=pa
        self.Headchild=hc
        self.Tailchild=tc
        self.Flip=flip
        self.Height=height
        self.Size=size
        self.Circle=circle
        self.Mark=mark
        self.Marksignal=marksignal


def get_genes(file):
    """
    Given a genome, store it in two different ways: as a gene list or a dictionary.
    For the gene list, index is in-order number of gene, value is all kinds of information belonging to this gene, which is stored in the form of a class called Gene.
    For the gene dictionary, key is gene id, value is the same as that of gene list.
    Also,set up a mapping between index of list and id of dictionary.
    """
    f=open(file,'r')
    genel=[]
    gened={}
    map={}
    list=f.read().split()
    j=0
    for i in list:
        genel.append(Gene(str(i),'Null','Null','Null',0,0,0,0,0,0))
        gened[str(i)]=Gene(str(i),'Null','Null','Null',0,0,0,0,0,0)
        if not map.has_key(j):
            map[j]=str(i)
        j+=1
    return [genel,gened,map,j]

class Genome:
    def __init__(self,genel,gened):
        """
        a class for chromosome,two attribure, Genes and Genome. Genes is list of all genes belonging to this chromosome in Gene class format, Genome is dictionary
        storing all genes in Gene class format.
        """
        self.Genes=genel
        self.Genome=gened

    def initGenome(self,map):
        """
        The genome will be initialized using gene list at first, then by the mapping information between in-order number and gene id,
        The genome can also be initilazid to gene dictionry.
        """
        #index of beginning of the gene list
        B=0
        #index of the end of the gene list
        E=len(self.Genes)-1
        #root for this chromosome, which is the middle gene from this list.
        root=self.Genes[int(math.ceil(float(E)/2))].Id 
        def Makechild(c,p):
            """
            given a possible child and its possible parent. If it is a real child, initialize child's parent, and return the child Id.
            If it is 'Null', just return 'Null' as the child.
            """
            if c!='Null':
                self.Genes[c].Parent=self.Genes[p].Id
                return self.Genes[c].Id
            else: 
                return c
            
        def Buildtree(b,e):
            """
            set child-parent relationships for each gene.
            """
            #curent parent is always in the middle of the chromosome
            m=int(math.ceil(float(b+e)/2))
            #if the index of beginning is still less than that of end,set the children of current node.
            if b<e:
                self.Genes[m].Headchild=Makechild(Buildtree(b,m-1),m)
                self.Genes[m].Tailchild=Makechild(Buildtree(m+1,e),m)
            if b<=e:
                return m
            else:
                return 'Null'
        Buildtree(B,E)
        #Buildtree initialize the parent-child relationship using Genes,but we want what we really want is to storint this information
        #in self.Genome. 
        for i in range(len(self.Genes)):
            self.Genome[map[i]]=self.Genes[i]
        self.Genes=[]
        return root

    def update_SH(self,node):
            
        left=self.Genome[node].Headchild
        right=self.Genome[node].Tailchild
            
        if (left=='Null')&(right=='Null'):
            self.Genome[node].Size=1
            self.Genome[node].Height=1
            self.Genome[node].Mark=self.Genome[node].Marksignal
        elif (left=='Null')&(right!='Null'):
            self.Genome[node].Size=self.Genome[right].Size+1
            self.Genome[node].Height=self.Genome[right].Height+1
            self.Genome[node].Mark=self.Genome[right].Mark+self.Genome[node].Marksignal
        elif (left!='Null')&(right=='Null'):
            self.Genome[node].Size=self.Genome[left].Size+1
            self.Genome[node].Height=self.Genome[left].Height+1
            self.Genome[node].Mark=self.Genome[left].Mark+self.Genome[node].Marksignal
        else:
            self.Genome[node].Size=self.Genome[left].Size+self.Genome[right].Size+1
            self.Genome[node].Height=max(self.Genome[left].Height,self.Genome[right].Height)+1
            self.Genome[node].Mark=self.Genome[left].Mark+self.Genome[right].Mark+self.Genome[node].Marksignal
# given one node, cleave the the chromosome before or after this node and decide the flip status of the resulting two parts.  
    def cleave(self,i,switch):
# get all the acestral nodes of the cut points, i.
        def visit(n,path):
            path.append(n)
            i=self.Genome[n].Parent
            if i != 'Null':
                visit(i,path)
# loop the list of acestral nodes(including i) and decide the flip status for each node.         
        def flipstatus(path):
            num=0
            flip={}
            path2=map(str,path)
            path2.reverse()
            for n in path2:
                num+=self.Genome[n].Flip
                if not flip.has_key(n):
                    flip[n]=0
                if (num%2)==0:
                    flip[n]=0
                else:
                    flip[n]=1
            return flip
# initialization of cleave, setting up the value of x, y, left and right respectively.
# x, current node; y, child node of x; left, root of the left subtree; right, root of the right subtree.
# initial value: x, parent of i; y, cutting point i; left, left child of i; right right child of i.
        path=[]
        visit(i,path)
        flip=flipstatus(path)
        y=path[0]
        left=''
        right=''
        if switch == 0:
            if flip[y] == 0:
                if self.Genome[y].Headchild != 'Null':
                    left=self.Genome[y].Headchild
                    self.Genome[y].Headchild='Null'
                    self.Genome[left].Parent='Null'
                    right=y
                    self.Genome[left].Flip ^=self.Genome[y].Flip
                else:
                    left=self.Genome[y].Headchild
                    right=y
            else:
                if self.Genome[y].Tailchild != 'Null':
                    left=self.Genome[y].Tailchild
                    self.Genome[y].Tailchild='Null'
                    self.Genome[left].Parent='Null'
                    right=y
                    self.Genome[left].Flip ^=self.Genome[y].Flip
                else:
                    left=self.Genome[y].Tailchild
                    right=y
        else:
            if flip[y] == 0:
                if self.Genome[y].Tailchild != 'Null':
                    right=self.Genome[y].Tailchild
                    self.Genome[y].Tailchild='Null'
                    self.Genome[right].Parent='Null'
                    left=y
                    self.Genome[right].Flip ^=self.Genome[y].Flip
                else:
                    right=self.Genome[y].Tailchild
                    left=y
            else:
                if self.Genome[y].Headchild != 'Null':
                    right=self.Genome[y].Headchild
                    self.Genome[y].Headchild='Null'
                    self.Genome[right].Parent='Null'
                    left=y
                    self.Genome[right].Flip ^=self.Genome[y].Flip
                else:
                    right=self.Genome[y].Headchild
                    left=y
        self.update_SH(y)

#loop trougth all ancestral nodes and update the left and right subtree each time.             
        for n in path:
            y=n
            x=self.Genome[y].Parent
            if x=='Null':
                break
            if ((y==self.Genome[x].Headchild)^(flip[x])):
                if right==y:
                    right=x
                else:
                    self.Genome[y].Parent='Null'
                    if flip[x]==0:
                        self.Genome[x].Headchild=right
                    else:
                        self.Genome[x].Tailchild=right
                    if right != 'Null':
                        self.Genome[right].Parent=x
                    right=x
                if left!='Null':
                    self.Genome[left].Flip ^=self.Genome[x].Flip
            #self.update_SH(x)
            else:
                if left==y:
                    left=x
                else:
                    self.Genome[y].Parent='Null'
                    if flip[x]==0:    
                        self.Genome[x].Tailchild=left
                    else:
                        self.Genome[x].Headchild=left
                    if left != 'Null':
                        self.Genome[left].Parent=x
                    left=x
                if right!='Null':
                    self.Genome[right].Flip ^=self.Genome[x].Flip
            self.update_SH(x)
        return (left,right)
        
#Check the invert state for each node    
    def checkInvert(self,node):
        print '#',node,self.Genome[node].Flip
        flip_all={}
        def visit(n,path,dic):
                print n
                if dic.has_key(n):
                    return dic[n] 
                path.append(n)
                i=self.Genome[n].Parent
                if i != 'Null':
                    return(visit(i,path,dic))
        def flipcheck(path,num,dic):
            path2=map(str,path)
            path2.reverse()
            for n in path2:
                if not dic.has_key(n):
                    dic[n]='Null'
                num^=self.Genome[n].Flip
                dic[n]=num
        
        if not flip_all.has_key(node):
            flip_all[node]=self.Genome[node].Flip
        for leaf in self.Genome.keys():
            if (self.Genome[leaf].Headchild=='Null')&(self.Genome[leaf].Tailchild=='Null'):
                path=[]
                a=visit(leaf,path,flip_all)
                print 'a',a
#print 'flip_all[node]',flip_all[node],type(node)
                flipcheck(path,a,flip_all)
        return flip_all
            
#tranvese the tree,print nodes in linear-order, and visulize the result by dot plot.    
    def treeWalker(self,root,dic,file):
        f=open(file,'w')
        f.write('digraph '+'id'+' '+'{'+'\n\tordering=out\n')
        order=[]
        def tranvese(node,invert,f,order):
            if invert[node]== 0:
                if self.Genome[node].Headchild != 'Null':
                    tranvese(self.Genome[node].Headchild,invert,f,order)
                order.append(self.Genome[node].Id)
                if (self.Genome[node].Headchild !='Null')|(self.Genome[node].Tailchild !='Null'):
                    if self.Genome[node].Headchild =='Null': 
                        f.write('\t'+self.Genome[node].Id+' -> '+'Hnil'+self.Genome[node].Id+'\n')
                        f.write('\t'+self.Genome[node].Id+' -> '+self.Genome[node].Tailchild+'\n')
                    elif self.Genome[node].Tailchild =='Null':
                        f.write('\t'+self.Genome[node].Id+' -> '+self.Genome[node].Headchild+'\n')
                        f.write('\t'+self.Genome[node].Id+' -> '+'Tnil'+self.Genome[node].Id+'\n')
                    else:
                        f.write('\t'+self.Genome[node].Id+' -> '+self.Genome[node].Headchild+'\n')
                        f.write('\t'+self.Genome[node].Id+' -> '+self.Genome[node].Tailchild+'\n')
                if self.Genome[node].Tailchild != 'Null':
                    tranvese(self.Genome[node].Tailchild,invert,f,order)
            if invert[node] != 0:
                if self.Genome[node].Tailchild != 'Null':
                    tranvese(self.Genome[node].Tailchild,invert,f,order)
                order.append(self.Genome[node].Id)
                if (self.Genome[node].Headchild !='Null')|(self.Genome[node].Tailchild !='Null'):
                    if self.Genome[node].Headchild =='Null': 
                        f.write('\t'+self.Genome[node].Id+' -> '+'Hnil'+self.Genome[node].Id+'\n')
                        f.write('\t'+self.Genome[node].Id+' -> '+self.Genome[node].Tailchild+'\n')
                        f.write('\t'+self.Genome[node].Id+' '+'[color=red]'+'\n')
                    elif self.Genome[node].Tailchild =='Null':
                        f.write('\t'+self.Genome[node].Id+' -> '+self.Genome[node].Headchild+'\n')
                        f.write('\t'+self.Genome[node].Id+' -> '+'Tnil'+self.Genome[node].Id+'\n')
                        f.write('\t'+self.Genome[node].Id+' '+'[color=red]'+'\n')
                    else:
                        f.write('\t'+self.Genome[node].Id+' -> '+self.Genome[node].Headchild+'\n')
                        f.write('\t'+self.Genome[node].Id+' -> '+self.Genome[node].Tailchild+'\n')
                        f.write('\t'+self.Genome[node].Id+' '+'[color=red]'+'\n')
                else:
                    f.write('\t'+self.Genome[node].Id+' '+'[color=red]'+'\n')
                if self.Genome[node].Headchild != 'Null':
                    tranvese(self.Genome[node].Headchild,invert,f,order)
        tranvese(root,dic,f,order)
        f.write('}')
        return order

    def checkInvert_test(self,node_list):
#print '#',node,self.Genome[node].Flip
#print '##',node_list
        flip_all={}
        def visit(n,path,dic):
#print n
                if dic.has_key(n):
                    return dic[n] 
                path.append(n)
                i=self.Genome[n].Parent
                if i != 'Null':
                    return(visit(i,path,dic))
        def flipcheck(path,num,dic):
            path2=map(str,path)
            path2.reverse()
            for n in path2:
                if not dic.has_key(n):
                    dic[n]='Null'
                num^=self.Genome[n].Flip
                dic[n]=num
        for node in node_list:
            if not flip_all.has_key(node):
                flip_all[node]=self.Genome[node].Flip
        for leaf in self.Genome.keys():
            if (self.Genome[leaf].Headchild=='Null')&(self.Genome[leaf].Tailchild=='Null'):
                path=[]
                a=visit(leaf,path,flip_all)
#print '@@',path
#               print 'a',a
                flipcheck(path,a,flip_all)
        return flip_all
# sys.exit()
    def treeWalker_test(self,root_list,dic,a=0,b=0,p=0,q=0):
        #f=open(file,'w')
#f.write('digraph '+'id'+' '+'{'+'\n\tordering=out\n')
#        f2=open(file2,'w')
#        f2.write('graph '+'id'+' '+'{'+'\n')
#        order=[]
#        print root_list
        def tranvese(node,invert,order):
            if invert[node]== 0:
                if self.Genome[node].Headchild != 'Null':
                    tranvese(self.Genome[node].Headchild,invert,order)
                order.append(self.Genome[node].Id)
#if (self.Genome[node].Headchild !='Null')|(self.Genome[node].Tailchild !='Null'):
#                    if self.Genome[node].Headchild =='Null': 
#                        f.write('\t'+self.Genome[node].Id+' -> '+'Hnil'+self.Genome[node].Id+'\n')
#                        f.write('\t'+self.Genome[node].Id+' -> '+self.Genome[node].Tailchild+'\n')
#                    elif self.Genome[node].Tailchild =='Null':
#                        f.write('\t'+self.Genome[node].Id+' -> '+self.Genome[node].Headchild+'\n')
#                        f.write('\t'+self.Genome[node].Id+' -> '+'Tnil'+self.Genome[node].Id+'\n')
#                   else:
#                        f.write('\t'+self.Genome[node].Id+' -> '+self.Genome[node].Headchild+'\n')
#                        f.write('\t'+self.Genome[node].Id+' -> '+self.Genome[node].Tailchild+'\n')
                if self.Genome[node].Tailchild != 'Null':
                    tranvese(self.Genome[node].Tailchild,invert,order)
            if invert[node] != 0:
                if self.Genome[node].Tailchild != 'Null':
                    tranvese(self.Genome[node].Tailchild,invert,order)
                order.append(self.Genome[node].Id)
#                if (self.Genome[node].Headchild !='Null')|(self.Genome[node].Tailchild !='Null'):
#                    if self.Genome[node].Headchild =='Null': 
#                        f.write('\t'+self.Genome[node].Id+' -> '+'Hnil'+self.Genome[node].Id+'\n')
#                        f.write('\t'+self.Genome[node].Id+' -> '+self.Genome[node].Tailchild+'\n')
#                        f.write('\t'+self.Genome[node].Id+' '+'[color=red]'+'\n')
#                    elif self.Genome[node].Tailchild =='Null':
#                        f.write('\t'+self.Genome[node].Id+' -> '+self.Genome[node].Headchild+'\n')
#                        f.write('\t'+self.Genome[node].Id+' -> '+'Tnil'+self.Genome[node].Id+'\n')
#                        f.write('\t'+self.Genome[node].Id+' '+'[color=red]'+'\n')
#                    else:
#                        f.write('\t'+self.Genome[node].Id+' -> '+self.Genome[node].Headchild+'\n')
#                        f.write('\t'+self.Genome[node].Id+' -> '+self.Genome[node].Tailchild+'\n')
#                        f.write('\t'+self.Genome[node].Id+' '+'[color=red]'+'\n')
#                else:
#                    f.write('\t'+self.Genome[node].Id+' '+'[color=red]'+'\n')
                if self.Genome[node].Headchild != 'Null':
                    tranvese(self.Genome[node].Headchild,invert,order)
        new_genome=[]
#        r_i=0
        
        for root in root_list:
#           output='/Users/jielv/genome_evolution/'+node+'/'+str(r_i)+'.txt'
#           r_i+=1
#           out=open(output,'w')
            order=[]
            tranvese(root,dic,order)
#            for each in order:
#               out.write(each+' ')
            order1=map(str,order)
            new_genome.append(order1)
#            if self.Genome[root].Circle==1:
#                f.write('\t'+str(root)+' '+'[style=filled]'+'\n')
#                f2.write('\t'+order[0]+' -- '+order[-1]+'\n')
#            while (1):
#                if len(order)<1:
#                    break
#                if len(order)<2:
#                    f2.write('\t'+order[0]+'\n')
#                else:
#                    f2.write('\t'+order[0]+' -- '+order[1]+'\n')
#                order.pop(0)
#        if a!=int(p):
#            f.write('\t'+'tel'+p+' '+'[shape=plaintext]'+'\n')
#            f2.write('\t'+'tel'+p+' '+'[shape=plaintext]'+'\n')
#        if q=='Null':
#            f.write('\t'+'tel'+q+' '+'[shape=plaintext]'+'\n')
#            f2.write('\t'+'tel'+q+' '+'[shape=plaintext]'+'\n')
#        elif b!=int(q):
#            f.write('\t'+'tel'+q+' '+'[shape=plaintext]'+'\n')
#            f2.write('\t'+'tel'+q+' '+'[shape=plaintext]'+'\n')
#        f.write('\t'+str(a)+' '+'[shape=box]'+'\n')
#        f.write('\t'+str(b)+' '+'[shape=diamond]'+'\n')
#        f2.write('\t'+str(a)+' '+'[shape=box]'+'\n')
#        f2.write('\t'+str(b)+' '+'[shape=diamond]'+'\n')
#                left=self.findLeftMostPoint(root)
#                right=self.findRightMostPoint(root)
#       f.write('\t'+left+' -> '+right+'\n')

#        f.write('}')
#        f2.write('}')
        return new_genome
        
#return order
# find the right most node of the left subtree
    def findRightMostPoint(self,root):
        if root =='Null':
            return root
        flag=[0]
        def find_right_most_point(root,flag):
            flag[0]^=self.Genome[root].Flip
            if flag[0]==0:
                next=self.Genome[root].Tailchild
                if next == 'Null':
                    return(root)
                else:    
                    return(find_right_most_point(next,flag))
            else:
                next=self.Genome[root].Headchild
                if next == 'Null':
                    return (root)
                else:    
                    return(find_right_most_point(next,flag))
        return(find_right_most_point(root, flag))
# find the left most node of the right subtree
    def findLeftMostPoint(self,root):
        if root == 'Null':
            return root
        flag=[0]
        def find_right_most_point(root,flag):
            flag[0]^=self.Genome[root].Flip
            if flag[0]==0:
                next=self.Genome[root].Headchild
                if next == 'Null':
                    return(root)
                else:    
                    return(find_right_most_point(next,flag))
            else:
                next=self.Genome[root].Tailchild
                if next == 'Null':
                    return (root)
                else:    
                    return(find_right_most_point(next,flag))
        return(find_right_most_point(root, flag))

# join the new root and two new subtrees together
    def rejoinSubtree(self,left,new_root,right):
        invert_bit=self.Genome[new_root].Flip
        if invert_bit == 0:
            self.Genome[new_root].Headchild=left
            self.Genome[new_root].Tailchild=right
        else:
            self.Genome[new_root].Headchild=right
            self.Genome[new_root].Tailchild=left
        if left != 'Null':
            self.Genome[left].Parent=new_root
            self.Genome[left].Flip^=invert_bit
        if right != 'Null':
            self.Genome[right].Parent=new_root
            self.Genome[right].Flip^=invert_bit
        self.update_SH(new_root)

# cut a circular chromosome
    def cutCircle(self,root,i):
        self.Genome[root].Circle=0
        (left,right)=self.cleave(i,0)
        (left,new_root,right)=self.getNewRoot(right,left)
        self.rejoinSubtree(left,new_root,right)
        return(left,new_root,right)

    def cutLine(self,i):
        (left,right)=self.cleave(i,0)
        (left,new_root,right)=self.getNewRoot(left,right)
        self.rejoinSubtree(left,new_root,right)
        return(left,new_root,right)
    def circulization(self,root):
        if root=='Null':
            return 'Null'
        def get_parent(n,path):
            path.append(n)
            i=self.Genome[n].Parent
            if i != 'Null':
                get_parent(i,path)
#print root
        if self.Genome[root].Parent != 'Null':
            path=[]
            get_parent(root,path)
            root=path[-1]
#print root
        self.Genome[root].Circle=1
        return root
#print '#',left,right

    def getParent(self,node):
        path=[]
        def get_parent(n,path):
                path.append(n)
                i=self.Genome[n].Parent
                if i != 'Null':
                    get_parent(i,path)
        if node != 'Null' and self.Genome[node].Parent != 'Null':
            get_parent(node,path)
            node=path[-1]
        return node


# before rejoin, get the root of left_subtree, new_root and root of the right_subtree    
    def getNewRoot(self,left,right):
#       print '@',left,right
        path=[]
        def get_parent(n,path):
                path.append(n)
                i=self.Genome[n].Parent
                if i != 'Null':
                    get_parent(i,path)

        if left != 'Null' and self.Genome[left].Parent != 'Null':
            get_parent(left,path)
            left=path[-1]
            path=[]
        if right != 'Null' and self.Genome[right].Parent != 'Null':
            get_parent(right,path)
            right=path[-1]
#print '#',left,right
        if left != 'Null' and right != 'Null':
            if self.Genome[left].Size>self.Genome[right].Size:
                new_root=self.findRightMostPoint(left)
                (left,new_root)=self.cleave(new_root,0)
                #print '%% cut at left'
            else:
                new_root=self.findLeftMostPoint(right)
                (new_root,right)=self.cleave(new_root,1)
                #print '&& cut at right'
        elif left == 'Null':
            new_root=self.findLeftMostPoint(right)
            (new_root,right)=self.cleave(new_root,1)
        elif right == 'Null':
            new_root=self.findRightMostPoint(left)
            (left,new_root)=self.cleave(new_root,0)
        
        return (left,new_root,right)

# given a gene, cut before the gene, and generate the resulting piceses(one if the original is circluar and two if the original is linear)
    def cut(self,i):
        (left,right)=self.cleave(i,0)
        if ((left != 'Null') and (self.Genome[left].Circle == 1)) or (self.Genome[right].Circle == 1):
            if (left != 'Null') and (self.Genome[left].Circle == 1):
                original_root=left
                self.Genome[left].Circle=0
            if self.Genome[right].Circle == 1:
                original_root=right
                self.Genome[right].Circle=0
            (left,new_root,right)=self.getNewRoot(right,left)
            self.rejoinSubtree(left,new_root,right)
            return (new_root,)
        else:
            return(left,right)
            
            



        
############################ Main program start#################
import os

def ini_genome(file):

    temp=get_genes(file)
    chrome=Genome(temp[0],temp[1])
    root=chrome.initGenome(temp[2])
    #visit(root,chrome)
    return (chrome,root)

def get_graph(chrome,root,dic,file):
    return(chrome.treeWalker(root,dic,file))

def visit(node,chrome):
    """
    initializing other information for each gene except child-parent relationship,like size,hight,Mark,and Marksignal.
    """
    left=chrome.Genome[node].Headchild
    right=chrome.Genome[node].Tailchild
    for i in [left,right]:
        if i != 'Null':
            visit(i,chrome)

    if (left=='Null')&(right=='Null'):
        chrome.Genome[node].Size=1
        chrome.Genome[node].Height=1
        chrome.Genome[node].Mark=chrome.Genome[node].Marksignal
    elif (left=='Null')&(right!='Null'):
        chrome.Genome[node].Size=chrome.Genome[right].Size+1
        chrome.Genome[node].Height=chrome.Genome[right].Height+1
        chrome.Genome[node].Mark=chrome.Genome[right].Mark+chrome.Genome[node].Marksignal
    elif (left!='Null')&(right=='Null'):
        chrome.Genome[node].Size=chrome.Genome[left].Size+1
        chrome.Genome[node].Height=chrome.Genome[left].Height+1
        chrome.Genome[node].Mark=chrome.Genome[left].Mark+chrome.Genome[node].Marksignal
    else:
        chrome.Genome[node].Size=chrome.Genome[left].Size+chrome.Genome[right].Size+1
        chrome.Genome[node].Height=max(chrome.Genome[left].Height,chrome.Genome[right].Height)+1
        chrome.Genome[node].Mark=chrome.Genome[left].Mark+chrome.Genome[right].Mark+chrome.Genome[node].Marksignal



def select_cut_test(j,n):
    if j<n:
        return(str(j),0)
    else:
        if j==n:
            return('1',1)
        else:
            return('4',1)
            
        #return(random.choice(['1','4']),1)
        
def select_cut(set,n,pool):
    #set includes all roots of linear chromosome
    j=random.randrange(n+len(set))
    # if choose a position which is not the end of chromosome(can be the head of the chromosome, because we choose to cut before gene)
    if j<n:
        return(pool[j],0)
    #if choose a end of a chomosome, randomly decide which.
    else:
        return(random.choice(list(set)),1)
def cut_chrom(i_list,genome):
    #if the cut point is at the end of a chromosome,return the root and ''
    if i_list[1]:
        return(i_list[0],'Null')
    else:
        return(genome.cut(i_list[0]))
def join_chrom(left,right,genome):
    if left=="Null" and right=="Null":
        return 'Null'
    (left,new_root,right)=genome.getNewRoot(left,right)
    genome.rejoinSubtree(left,new_root,right)
#print 'left new_root right',left,new_root,right
    return new_root

def DCJ_test(genome,n,pool,set1,set2,count,Ni,Nt,move=None,result=[]):

#print m,n
#roots of pieces need to be fliped
    flip=[]                                         
#pieces need to be joined together
    join=[]
#pieces need to be circulization
    cir=[]
#roots of pieces directly added to the new roots list
    add=[]
#list of pieces inverted, used for control the inversion size
    inversion=[]
# a list storing the number of maked genes for different original chromosome( 1 or 2)
    root_mark=[]
# pieces that will be seperated from it's original chromosome, for each seperation, just pick one, used for control DCJ with balanced gene.
    frag=[]
# pairs of pieces from different original chromomsomes that will be joined together in new chromosomes
    merge=[]
#    if move:
#        (m,n,flag)=move
#    print '#', (m,n,flag)
    a=select_cut(set2,n,pool)
    b=select_cut(set2,n,pool)
    #determine which chromosome the cut point belong
    a_root=genome.getParent(a[0])
    b_root=genome.getParent(b[0])
    #determine how many marked genes the chosen chromosome have
    a_mark=genome.Genome[a_root].Mark
    b_mark=genome.Genome[b_root].Mark
    #put this information in a list, if the two cut points come from the same chromosome, just storing one
    root_mark.append(a_mark)
    if a_root != b_root:
        root_mark.append(b_mark)
    
    piecesA=cut_chrom(a,genome)
    A=len(piecesA)
    if A==2:
        x=piecesA[0]
        y=piecesA[1]
#       print '**',b[0],x,y
        if (b[1] and (b[0]==x)):
            b=(y,1)
    piecesB=cut_chrom(b,genome)
    B=len(piecesB)
    flag=random.choice((0,1))
#print '#',a
#    print '#',b
#    print '#',flag
#    print '!',piecesA
#    print '!',piecesB
#   for each in genome.Genome:
#        print '@@',each,genome.Genome[each].Mark
    U3=len(set(piecesA+piecesB))
    U1=len(set(piecesA+piecesB).difference(['Null']))
    U2=0
    for i in piecesA+piecesB:
            if i == 'Null':
                U2+=1
    U=U1+U2
#print A,B,U,U1,U2,U3
    if (A==1 and B==1):                             # two cuts at different circluar chromosome
        if flag == 0:
            flip.append(piecesB[0])
        #    genome.Genome[piecesB[0]].Flip^=1
        join.append((piecesA[0],piecesB[0]))
        merge.append((a_mark,b_mark))
        if piecesA[0] != 'Null':
            cir.append(piecesA[0])
        else:
            cir.append(piecesB[0])
            
#        new_root=join_chrom(piecesA[0],piecesB[0])
#        genome.Genome[new_root].Circle=1
#        add_list=(new_root,)
    elif (A==1 and B==2 and U==3) or(A==2 and B==1):            # two cuts at one circluar, one linear chromosome
        if A>B:
            (piecesA,piecesB)=(piecesB,piecesA)
        if flag == 0:
            flip.append(piecesA[0])
#genome.Genome[piecesA[0]].Flip^=1
        join.append((piecesB[0],piecesA[0]))
        join.append((piecesA[0],piecesB[1]))
        merge.append((a_mark,b_mark))
#        new_root=join_chrom(piecesB[0],piecesA[0])
#new_root=join_chrom(new_root,piecesB[1])    
#        add_list=(new_root,)
    elif (A==1 and B==2 and U==2):                              #two cuts at the same circular chromosome
        if piecesB[0] == 'Null':                                #two cuts at the same point of the same circular chromosom
            add.append(piecesB[1])
#            add_list=(piecesB[1],)
        else:
            if flag == 0:
                flip.append(piecesB[1])
                if genome.Genome[piecesB[1]].Size>genome.Genome[piecesB[0]].Size:
                    inversion.append(piecesB[0])
                else:
                    inversion.append(piecesB[1])
#               genome.Genome[piecesB[1]].Flip^=1
                join.append((piecesB[0],piecesB[1]))
                cir.append(piecesB[1])
#               new_root=join_chrom(piecesB[0],piecesB[1])
#genome.Genome[new_root].Circle=1
#                add_list=(new_root,)
            else:
                frag.append(piecesB[0])
                cir.append(piecesB[0])
                cir.append(piecesB[1])

#                genome.Genome[piecesB[0]].Circle=1
#                genome.Genome[piecesB[1]].Circle=1
#                add_list=piecesB
    elif (A==2 and B==2 and U==4 and U3!=2):                #two cuts at two different linear chromosome,U3!=2 used to exclude the situation that two cuts at same tail
        def trans(a):
            if a != 'Null':
                a=genome.Genome[a].Mark
            else:
                a=0
            return a
        if flag == 0:
            join.append((piecesA[0],piecesB[1]))
            join.append((piecesB[0],piecesA[1]))
            merge.append((trans(piecesA[0]),trans(piecesB[1])))
            merge.append((trans(piecesB[0]),trans(piecesA[1])))
#            if piecesA[0] != 'Null':
#                join.append(piecesA[0])
#            else:
#                join.append(piecesB[1])

#            new_root1=join_chrom(piecesA[0],piecesB[1])
#            new_root2=join_chrom(piecesB[0],piecesA[1])
        else:
            if piecesB[0]!= 'Null':
                flip.append(piecesB[0])
            join.append((piecesA[0],piecesB[0]))
                
#                genome.Genome[piecesB[0]].Flip^=1
#            new_root1=join_chrom(piecesA[0],piecesB[0])
            if piecesA[1]!= 'Null':
                flip.append(piecesA[1])
            join.append((piecesA[1],piecesB[1]))
            merge.append((trans(piecesA[0]),trans(piecesB[0])))
            merge.append((trans(piecesA[1]),trans(piecesB[1])))
#merge.append((piecesA[0],piecesB[0]))
#            merge.append((piecesA[1],piecesB[1]))
        frag.append(piecesA[0])
        frag.append(piecesB[0]) 
#                genome.Genome[piecesA[1]].Flip^=1
#            new_root2=join_chrom(piecesA[1],piecesB[1])
#        add_list=(new_root1,new_root2)
    elif (A==2 and B==2 and U==3) or U3==2:             #two cuts at a same linear chromosome, including two cuts at the tail 
        #two cuts at the middle and is the same point of a linear chromosome
        if (piecesB[0]=='Null' and piecesA[1]==piecesB[1]) and (piecesA[0] != 'Null') and (piecesA[1] != 'Null'):
            add.append(piecesA[0])
            add.append(piecesA[1])
            frag.append(piecesA[0])
#            add_list=piecesA
        else:
            if (piecesA[0]==piecesB[0] or piecesA[0]==piecesB[1]) and piecesA[0] != 'Null':     # second cut happend on the left fragment of the first cut
                if flag == 0:
                    join.append((piecesB[0],piecesA[1]))
                    if piecesB[1] != 'Null':
                        cir.append(piecesB[1])
                        frag.append(piecesB[1])
                    
#                    new_root1=join_chrom(piecesB[0],piecesA[1])
#                    new_root2=piecesB[1]
#                    if new_root2 != 'Null':
#                        genome.Genome[new_root2].Circle=1
#                    add_list=(new_root1,new_root2)
                else:
                    if piecesB[1]!= 'Null':
                        flip.append(piecesB[1])
                        inversion.append(piecesB[1])
                    join.append((piecesB[0],piecesB[1]))
                    if piecesB[0] != 'Null':
                        join.append((piecesB[0],piecesA[1]))
                    else:
                        join.append((piecesB[1],piecesA[1]))
#                        genome.Genome[piecesB[1]].Flip^=1
#                    new_root1=join_chrom(piecesB[0],piecesB[1])
#                    new_root2=join_chrom(new_root1,piecesA[1])
#                    add_list=(new_root2,)
            else:                                                                           #second cut happend on the right fragment of the first cut
                if flag == 0:
                    join.append((piecesA[0],piecesB[1]))
                    if piecesB[0] != 'Null':
                        cir.append(piecesB[0])
                        frag.append(piecesB[0])
#                    new_root1=join_chrom(piecesA[0],piecesB[1])
#                    new_root2=piecesB[0]
#                    if new_root2 != 'Null':
#                        genome.Genome[new_root2].Circle=1
#                    add_list=(new_root1,new_root2)
                else:
                    if piecesB[0]!= 'Null':
                        flip.append(piecesB[0])
                        inversion.append(piecesB[0])
                    join.append((piecesA[0],piecesB[0]))
                    if piecesA[0] != 'Null':
                        join.append((piecesA[0],piecesB[1]))
                    else:
                        join.append((piecesB[0],piecesB[1]))
#                    if piecesB[0]!= 'Null':
#                        genome.Genome[piecesB[0]].Flip^=1
#                    new_root1=join_chrom(piecesA[0],piecesB[0])
#                    new_root2=join_chrom(new_root1,piecesB[1])
#                    add_list=(new_root2,)
#    print add_list
#    print set1,set2
#    update(set1,set2,add_list)
    signal=0
#    for i in inversion:
#        if genome.Genome[i].Size>4:
#            signal=1
#            inversion_list.append((a,b))
#            print 'invertion size too big!!'
    for i in range(len(frag)):
        if (frag[i] == 'Null') or (genome.Genome[frag[i]].Mark == 0) or (genome.Genome[frag[i]].Mark == root_mark[i]):
            continue
        else:
            signal=1
            break
#print '$',merge
    for i in merge:
#        if len(mer_list)>0:
#            break
        if (i[0]==0) or (i[1]==0):
#print 'continue',merge
            continue
        else:
            signal=1
            break
#            mer_list.append((a,b)
    if signal==1:
        undo_cut(piecesB,set1,set2,genome)
        undo_cut(piecesA,set1,set2,genome)
    else:
#        mer_list.append((a,b,flag))
        if len(cir)==0 and len(inversion)>0:
            Ni[0]+=1
        if len(cir)==0 and len(frag)>0:
            Nt[0]+=1
        do_join(flip,join,cir,add,set1,set2,genome)
        count[0]+=1
    if move:
        result.append((signal,a,b,flag))
        print '#',move
        print '#',result
    
    return (a[0],b[0],signal,a,b,flag)

def update(set1,set2,add_list,genome):
    set1.update(add_list)
    temp1=set1.copy()
    for i in temp1:
        if i=='Null' or genome.Genome[i].Parent != 'Null':
            set1.difference_update([i])
    temp1=set1.copy()
    set2.clear()
    for i in temp1:
        if genome.Genome[i].Circle == 1:
            continue
        set2.update([i])
#    print set1,set2

def undo_cut(pieces,set1,set2,genome):
    if len(pieces)==1:
        genome.Genome[pieces[0]].Circle=1
        undo_list=[pieces[0]]
    else:
        new_root=join_chrom(pieces[0],pieces[1],genome)
        undo_list=[new_root]
    update(set1,set2,undo_list,genome)

def do_join(f,j,c,a,set1,set2,genome):
    add_list=[]
    if len(f)>0:
        for i in f:
            genome.Genome[i].Flip^=1
    if len(j)>0:
        for i in j:
            new_root=join_chrom(i[0],i[1],genome)
            add_list.append(new_root)
    if len(c)>0:
        for i in c:
            genome.circulization(i)
            add_list.append(i)
    if len(a)>0:
        for i in a:
            add_list.append(i)
    update(set1,set2,add_list,genome)
#    return add_list
def get_chromosome_content(node,chrome,gene_set):
    """
    given a root, get IDs of all genes belong to the chromosome.
    """
    left=chrome.Genome[node].Headchild
    right=chrome.Genome[node].Tailchild
    for i in [left,right]:
        if i != 'Null':
            get_chromosome_content(i,chrome,gene_set)
    gene_set.add(node)


############################mian program###############
unit=int(float(sys.argv[1]))
#print unit
#sys.exit() 
    
genome1=Genome([],{})
def initial_genome(genome):            
    root_set=set()
    linear_root_set=set()
    for name in os.listdir(sys.argv[3]):
        file=sys.argv[3]+'/'+name
        (g1,g2)=ini_genome(file)
        genome.Genome=dict(genome.Genome,**g1.Genome)
        root_set.update([g2])
    pool=genome.Genome.keys()
    n=len(pool)
    mark_list=random.sample(pool,int(float(sys.argv[2])*n)/100)
    dic_mark={}
    #chrom_dic={}
    for w in mark_list:
        w=str(w)
        genome.Genome[w].Marksignal=1
        r=genome.getParent(w)
        if not dic_mark.has_key(r):
            dic_mark[r]=[]
        dic_mark[r].append(w)
    for node in root_set:
        visit(node,genome)
    for i in root_set:
        if genome.Genome[i].Circle == 0:
            linear_root_set.update([i])
    return pool,n,root_set,linear_root_set,dic_mark



def simulation(genome,n,pool,step,root_set,dic_mark):
    count=[0]
    Ni=[0]
    Nt=[0]
    Na=0
    linear_root_set=set()
    for i in root_set:
        if genome.Genome[i].Circle == 0:
            linear_root_set.update([i])
    while count[0]<step:
        (p,q,signal,a_test,b_test,flag)=DCJ_test(genome,n,pool,root_set,linear_root_set,count,Ni,Nt)
        Na+=1
#       print "&&&",p,q
#       print root_set,linear_root_set
#    print Na,Ni[0],Nt[0],
    dic={}
    dic=genome.checkInvert_test(root_set)
#    print len(dic)
    if len(dic) != n:
        print '##something wrong',len(dic),n
        print dic
        sys.exit()
#os.remove(node)
#    os.mkdir(node)
    new_chrom=genome.treeWalker_test(root_set,dic)
    for group in dic_mark.values():
        mark_root=genome.getParent(group[0])
        if genome.Genome[mark_root].Mark != len(group):
            print 'something wrong!',mark_root,genome.Genome[mark_root].Mark,group
            sys.exit()
        for each in group:
            mark_root2=genome.getParent(each)
            if mark_root != mark_root2:
                print 'something wrong!',mark_root,mark_root2
                sys.exit()
#    for i in new_chrom:
#print i
#    print count
    return new_chrom

pool1,n1,root_set_ori,linear_root_set_ori,dic_mark=initial_genome(genome1)

#sim_genome is the simulated genome stored in a list, whose element is again a list storing a whole chrmosome.
sim_genome=simulation(genome1,n1,pool1,unit,root_set_ori,dic_mark)


################generate synin file###############

def genome_dic(folder):
    """
    given a directory of genome file, generating a dictonary, with key be each gene and value be a tuple, the first element is chromosome position and second element be the    cordinate
    """
    i=0  
    j=0  
    dict={}
    for name in os.listdir(folder):
        file=folder+'/'+name
        f=open(file,'r')
        list=f.read().split()
        if len(list)==0:
            continue
        for each in list:
            if each not in dict:
                dict[each]=(i,j)
            j+=1 
        i+=1 
    return dict 

dict_original=genome_dic(sys.argv[3])

chn=0
posi=0
dict_simulated={}

for chr in sim_genome:
    if len(chr)==0:
        continue
    for each in chr: 
        if each not in dict_simulated:
            dict_simulated[each]=(chn,posi)
        posi+=1
    chn+=1

c=0
for i in dict_original:
    print 'g1.'+str(i)+'\t'+'g1.'+str(dict_original[i][0])+'\t'+str(dict_original[i][1])+'\t'+'g2.'+str(i)+'\t'+'g2.'+str(dict_simulated[i][0])+'\t'+str(dict_simulated[i][1])
    c+=1 

