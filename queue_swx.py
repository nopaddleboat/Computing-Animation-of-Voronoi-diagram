import numpy as np
import quicksort_swx as qsort

class qnode:#queue node,actually, a tree node
    def __init__(self,coord):
        self.Coord=coord
        self.lc=None#left child
        self.rc=None#right child
        self.pr=None#parent
        self.leaf=None#point to leaf of status tree for circle event, None for site event
        self.cnt=None#lowpnt of the circle event
class node_q:
    def __init__(self,qnode):
        self.pre=None
        self.next=None
        self.qnode=qnode

def genBST(sortedsites):#generate BST using sorted events
    len1=sortedsites.shape[0]#length of the sites list
    #print("len1: ",len1)
    if len1==1:
        root=qnode(sortedsites[0,:])
        return root
    elif len1==2:
        root=qnode(sortedsites[0,:])
        root.rc=qnode(sortedsites[1,:])
        root.rc.pr=root
        return root
    else:
        mid=int((len1+1)/2)
        #print("mid",mid)
        root=qnode(sortedsites[mid-1,:])
        root.lc=genBST(sortedsites[0:mid-1,:])
        root.rc=genBST(sortedsites[mid:len1,:])
        root.lc.pr=root
        root.rc.pr=root
        return root

class queue_tree:
    def __init__(self,node_q):
        self.head=node_q
        self.tail=node_q
    def insert(self,node_q):#when insert a node, the queue should not be empty
        node_q.next=self.tail
        self.tail.pre=node_q
        self.tail=node_q
    def popout(self):
        self.head=self.head.pre

def printtree(root):#print the tree rooted at "root": using breadth first search
    qt1=queue_tree(node_q(root))
    while qt1.head is not None:
        print(qt1.head.qnode.Coord)
        # if qt1.head.qnode.pr is not None:
        #     print("pr.lc? pr.rc?",(qt1.head.qnode.pr.lc==qt1.head.qnode)or(qt1.head.qnode.pr.rc==qt1.head.qnode))
        if qt1.head.qnode.lc is not None:
            qt1.insert(  node_q(qt1.head.qnode.lc)  )
        if qt1.head.qnode.rc is not None:
            qt1.insert(  node_q(qt1.head.qnode.rc)  )
        qt1.popout()


def cmpsites(site1,site2):#compare qnode
    if (site1.Coord[1] > site2.Coord[1]) or \
       ((site1.Coord[1] == site2.Coord[1]) and (site1.Coord[0] < site2.Coord[0])):
         return True #site1>site2
    else:
         return False #site1<site2
 
class queue_v:#queue for computing voronoi diagram
    def __init__(self,sites):#use n sites to initialize Q
        indexes=np.arange(0,len(sites))
        sortedsites,sortedindexes=qsort.quicksort(sites,indexes)
        newsites=np.zeros([len(sites),3])#third indexes show the index 0: highest site 1:second highest
        newsites[:,2]=np.arange(len(sites))
        newsites[:, 0:2]=sortedsites
        print("sortedsites: ",newsites)
        root=genBST( np.flip(newsites,0))#  sortedsites
        # print("Initial Q: ")
        printtree(root)
        # print("end of initial Q:")
        self.root=root  #pointer to the root
        #find the largest node
        self.sortedlist=np.copy(newsites)
        p=root
        while True:
            if p.rc is None:
                break
            p=p.rc
        print("outq of initial Q: ",p.Coord)
        self.outq=p    #point to the largest element
    def insert(self,event):#operate on qtree,event is q qnode
        current=self.root
        if current is None:
            print("current is None!")
            self.root=event
            self.outq=event
        else:
            while True :
                if cmpsites(current,event) is True:#event<current
                    if current.lc is None:
                        event.pr=current
                        current.lc=event
                        break
                    else:
                        current=current.lc
                else:  #event>current
                     if current.rc is None:
                        event.pr=current
                        current.rc=event
                        break
                     else:
                        current=current.rc
        #update the top element in the queue
            if self.outq.rc is not None:
                self.outq=self.outq.rc
    def delete(self,event):
        #check if event is root
        if event.pr is None:#event is self.root
            print("delete:event is root")
            if event.rc is None and event.lc is None:
                print("event.pr.rc is None and event.pr.lc is None")
            if  event.lc is None and event.rc is None: #event has no child
                self.root=None
                print("(event.lc==None) and (event.rc==None)")
            elif event.lc is None: #has one child:left child
                self.root=event.rc
                event.rc.pr = None
                print("delete line129")
            elif event.rc is None: #has one child:right child
                print("delete line131")
                self.root=event.lc
                event.lc.pr=None
            else:#has two children
                print("has two children")
            # #find the smallest element/leaf in the right subtree,this element will be the new root
                p=event.rc
                while True:
                    if p.lc is not None:
                        p=p.lc
                    else:#p is what we want to replace event,break
                        if p.rc is None:
                            if p.pr != event:
                                p.pr.lc=None
                        else:
                            if p.pr!=event:
                                p.pr.lc=p.rc
                        p.lc = event.lc#operations to replace event with p
                        event.lc.pr = p
                        if event.rc!=p:
                            p.rc = event.rc
                            event.rc.pr = p
                        p.pr = None
                        self.root = p
                        print("Found p: ",p.Coord)
                        break
        else:#event is not self.root
            print("event is not Root",event.Coord)
            if (event.lc is None) and (event.rc is None): #event has no child
                # left child or right
                if event.pr.lc==event:#left child
                    event.pr.lc=None
                else:
                    event.pr.rc=None
                    # print("line157")
            elif event.lc is None: #has one child
                event.rc.pr = event.pr
                if event.pr.lc==event:
                    event.pr.lc=event.rc
                else:
                    event.pr.rc=event.rc
                    # print("line164")
            elif event.rc is None: #has one child
                # if event.rc is None:
                #     if event.pr.lc is not None:
                #         print("line171",event.pr.lc.Coord)
                #     if event.pr.rc is not None:
                #         print("line171", event.pr.rc.Coord)
                #     print(event.pr.lc,event.pr.rc,event,event.lc.Coord)
                event.lc.pr=event.pr
                if event.pr.lc==event:
                    event.pr.lc=event.lc
                else:
                    event.pr.rc=event.lc
                if event.pr.rc is not None:
                    print("event.pr.rc.Coord",event.pr.rc.Coord)
            else:#has two children    #find the smallest element in the right subtree
                p=event.rc
                while True:
                    if p.lc is not None:
                        p=p.lc
                    else:
                        if p.rc is None:
                            if p.pr != event:
                                p.pr.lc=None
                        else:
                            if p.pr!=event:
                                p.pr.lc=p.rc
                        p.lc = event.lc  # operations to replace event with p
                        event.lc.pr = p
                        if event.rc!=p:
                            p.rc = event.rc
                            event.rc.pr = p
                        p.pr = event.pr
                        if event.pr.lc==event:
                            event.pr.lc=p
                        else:
                            event.pr.rc = p
                        break
                print("line187")
        #update outq:finding max node

        p = self.root
        n=1
        if p is not None:
            while n<100:
                n=n+1
                if p.rc is None:
                    break
                p = p.rc
                print(p.Coord)
        else:
            print("outq is None")
        if n>90:
            print("delete line206", self.root.Coord, p.rc, p.rc.rc)
        self.outq = p
        # print("p",p.Coord)

    def printall(self):
        printtree(self.root)

#    def search(self,event):
##    def printall(self):#print all nodes
