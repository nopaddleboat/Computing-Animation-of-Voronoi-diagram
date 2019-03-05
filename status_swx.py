import numpy as np
import math

def combreakpoint(sitei,sitej,ys):#Ys is y coordiante of sweep line
    xy = np.array([0.0, 0.0])
    if ys==sitei[1]:
        xy[0]=sitei[0]
        xy[1]=((sitei[0]-sitej[0])**2+sitej[1]**2-ys**2)/2.0/(sitej[1]-ys)
    elif ys==sitej[1]:
        xy[0] = sitej[0]
        xy[1] = ((sitej[0] - sitei[0]) ** 2 + sitei[1]**2 - ys ** 2) / 2.0 / (sitei[1] - ys)
    else:
        xi=sitei[0]
        yi=sitei[1]
        xj = sitej[0]
        yj = sitej[1]
        xyc=(sitei[0:2]+sitej[0:2])/2.0
        d=math.sqrt( (xi-xj)**2+(yi-yj)**2 )/2.0
        ny=(xj-xi)/d/2.0
        nx = -(yj - yi) /d/2.0
        lamda=np.roots( [ny**2-1,2*ny*(xyc[1]-ys),(xyc[1]-ys)**2-d**2] )
        x1=xyc[0]+lamda[0]*nx
        if (x1-xj)*(yi-ys)<(x1-xi)*(yj-ys):
            y1 = xyc[1] + lamda[0] * ny

        else:
            x1=xyc[0]+lamda[1]*nx
            y1=xyc[1]+lamda[1]*ny
        xy[0] = x1
        xy[1] = y1
    return xy


class leafnode:#leaf node of the status tree
    def __init__(self,site):
        self.site=site#site defining the arc:coordinate
        self.eventnode=None#usually circle event in Q indicating that this arc will disappear
        self.pr=None#parent node,No next node, since it's a leaf
        self.lsite = None  # left site defining the leafnode:always None
        self.rsite = None  # right site defining the rightnode:always None
        self.trinode1=None#record two other arcs of the in the circle event
        self.trinode2=None
        self.lc = None  # left child
        self.rc = None  # right child
    def setvalue(self,parent,circleevent):
        if parent is not None:
            self.pr=parent
        if circleevent is not None:
            self.eventnode=circleevent
class innode:#internal node of the status tree
    def __init__(self,lsite,rsite):
        self.lsite=lsite#left site defining the breakpoint:coords
        self.rsite=rsite#right site defining the break point:coords
        self.Enode=None#point to edge node in DCEL
        self.lc=None#left child
        self.rc=None#right child
        self.pr=None#parent node
        self.dir=np.array([0.0,0.0,0.0])#1 up;-1down
    def setvalue(self,Enode,lc,rc,pr):
        self.Enode = None  # point to edge node in DCEL
        if lc is not None:
            self.lc = lc  # left child
        if rc is not None:
            self.rc = rc  # right child
        if pr is not None:
            self.pr = pr  # parent node
    def setdir(self,val):
        self.dir=val

def rightRotate(rootnode):
    A=rootnode.lc.lc
    B=rootnode.lc.rc
    C=rootnode.rc
    rootnodeL=rootnode.lc
    if rootnode.pr is None:
        rootnodeL.pr=None
    elif rootnode.pr.lc==rootnode:
        rootnode.pr.lc=rootnodeL
        rootnodeL.pr=rootnode.pr
    else:# rootnode.pr.rc==rootnode:
        rootnode.pr.rc = rootnodeL
        rootnodeL.pr = rootnode.pr
    rootnodeL.rc = rootnode
    rootnode.pr = rootnodeL
    rootnodeL.lc = A
    A.pr = rootnodeL
    rootnode.lc = B
    B.pr = rootnode
    rootnode.rc = C
    C.pr = rootnode
    return rootnodeL

def leftRotate(rootnode):
    A = rootnode.lc
    B = rootnode.rc.lc
    C = rootnode.rc.rc
    rootnodeR = rootnode.rc
    if rootnode.pr is None:
        # print("rootnode.pr is None")
        rootnodeR.pr=None
    elif rootnode.pr.lc is rootnode:
        rootnode.pr.lc = rootnodeR
        rootnodeR.pr = rootnode.pr
    else:  # rootnode.pr.rc==rootnode:
        rootnode.pr.rc = rootnodeR
        rootnodeR.pr = rootnode.pr
    rootnodeR.lc = rootnode
    rootnode.pr = rootnodeR
    rootnodeR.rc = C
    C.pr = rootnodeR
    rootnode.lc = A
    A.pr = rootnode
    rootnode.rc = B
    B.pr = rootnode
    # print(A,B,C)
    return rootnodeR

def splay(Sroot,node1,ysweep):
    if Sroot is None or Sroot.lc is None:
        return Sroot
    SrootBp = combreakpoint(Sroot.lsite, Sroot.rsite, ysweep)
    InBp = combreakpoint(node1.lsite, node1.rsite, ysweep)
    print("line121",SrootBp[0], InBp[0])
    if SrootBp[0]==InBp[0]:
        return Sroot
    if SrootBp[0]>InBp[0]:#key lies in left subtree
        print("line124>")
        if Sroot.lc.lc is None:
            print("line126None")
            return Sroot
        SrootlcBp=combreakpoint(Sroot.lc.lsite, Sroot.lc.rsite, ysweep)
        if SrootlcBp[0]>InBp[0]:#zig-zig
            print("line129>")
            splay(Sroot.lc.lc,node1,ysweep)
            Sroot=rightRotate(Sroot)
        elif SrootlcBp[0]<InBp[0]:#zig-zag
            print("line133>")
            splay(Sroot.lc.rc,node1,ysweep)
            if Sroot.lc.rc.lc is not None:
                leftRotate(Sroot.lc)
        if Sroot.lc.lc is None:#do second rotate
            print("line139None")
            return Sroot
        else:
            print("line141None")
            return rightRotate(Sroot)
    else:
        print("line141<")
        if Sroot.rc.lc is None:
            print("line146None")
            return Sroot
        SrootrcBp=combreakpoint(Sroot.rc.lsite, Sroot.rc.rsite, ysweep)
        print(SrootrcBp[0],InBp[0])
        if SrootrcBp[0]>InBp[0]:#zag-zig
            print("line148<")
            splay(Sroot.rc.lc,node1,ysweep)
            if Sroot.rc.lc.lc is not None:
                rightRotate(Sroot.rc)
        elif SrootrcBp[0] < InBp[0]:  # zig-zag
            print("line152<")
            splay(Sroot.rc.rc, node1,ysweep)
            Sroot=leftRotate(Sroot)
        if Sroot.rc.lc is None:#do second rotate
            print("line159None")
            return Sroot
        else:
            print("line163None")
            return leftRotate(Sroot)


class statustree:
    def __init__(self):
        self.root=None#must be an innode
    def setroot(self,event):
        if self.root is None:
            self.root=event
        else:
            print("Error! a root already exists!")
    def replaceleaf(self,leaf,newroot):#replace the leaf with a root newroot
        if (self.root.lsite is None) and (self.root.rsite is None):#leaf is root
            self.root=newroot
            print("root replaced!")
        else:#leaf is not root
            newroot.pr=leaf.pr
            if leaf.pr.lc==leaf: #replace leaf is left child
                leaf.pr.lc=newroot
            elif leaf.pr.rc==leaf: #replace leaf is right child
                leaf.pr.rc=newroot
            leaf.pr=None
            print("leaf replaced!")
    def deleteleaf(self,leaf):
        print("deleteleaf: ",leaf.site)
        if leaf.pr.pr is None:#leaf is child of root
            if leaf.pr.lc==leaf:#leaf is a leaf child
                if leaf.pr.rc.lsite is None:#rchild is a leaf
                    self.root=leaf
                else:#innode
                    self.root=leaf.pr.rc
                    leaf.pr.rc.pr=None
            elif leaf.pr.rc==leaf:#leaf is a right child
                if leaf.pr.lc.lsite is None:#rchild is a leaf
                    self.root=leaf
                else:#innode
                    self.root=leaf.pr.lc
                    leaf.pr.lc.pr=None
        else:
            #leaf is not child of root
            flag_leafpr=-1#leaf's parent is a left child
            flag_leaf=-1#leaf is left or right
            if leaf.pr.pr.rc==leaf.pr:
                flag_leafpr=1
            if leaf.pr.rc==leaf:
                flag_leaf=1

            if flag_leaf==-1:
                if flag_leafpr == -1:
                    leaf.pr.pr.lc = leaf.pr.rc
                else:
                    leaf.pr.pr.rc = leaf.pr.rc
                leaf.pr.rc.pr=leaf.pr.pr
            elif flag_leaf==1:
                if flag_leafpr == -1:
                    leaf.pr.pr.lc = leaf.pr.lc
                else:
                    leaf.pr.pr.rc = leaf.pr.lc
                leaf.pr.lc.pr=leaf.pr.pr
        leaf.pr = None

    def depoutT(self,root):#output leaves from left to right
        if root.lc is None:#reach a leaf
            print(root.site)
            return None
        else:
            self.depoutT(root.lc)
            self.depoutT(root.rc)
    def widoutT(self,ysweep):#output breakpoints
        pall = self.root
        innodeList = []
        innodeList.append(pall)
        bplist = []
        leavescoord = []
        breaknode = []
        while len(innodeList) > 0:
            tarbppall = combreakpoint(innodeList[0].lsite, innodeList[0].rsite, ysweep)  # target breakpoint
            bplist.append(tarbppall)
            breaknode.append(innodeList[0])
            inDel = innodeList[0]
            if inDel.lc.lc is not None:
                innodeList.append(inDel.lc)
                # innodeList.insert(0, inDel.lc)
            if inDel.rc.lc is not None:
                innodeList.append(inDel.rc)
                # innodeList.insert(0, inDel.rc)

            if inDel.lc.lc is None:
                leavescoord.append(inDel.lc.site)
            if inDel.rc.lc is None:
                leavescoord.append(inDel.rc.site)
            #
            innodeList.remove(inDel)
        print("all breakpoints!")
        # print(np.asarray(bplist))
        for i in range(len(breaknode)):
            print(bplist[i], breaknode[i].lsite, breaknode[i].rsite)
    def outputTree(self,ysweep):
        pall = self.root
        innodeList = []
        innodeList.append(pall)
        bplist = []
        leavescoord = []
        breaknode = []
        while len(innodeList)>0:
            tarbppall = combreakpoint(innodeList[0].lsite, innodeList[0].rsite, ysweep)  # target breakpoint
            bplist.append(tarbppall)
            breaknode.append(innodeList[0])
            inDel=innodeList[0]
            if inDel.lc.lc is not None:
                innodeList.append(inDel.lc)
                # innodeList.insert(0, inDel.lc)
            if inDel.rc.lc is not None:
                innodeList.append(inDel.rc)
                # innodeList.insert(0, inDel.rc)


            if inDel.lc.lc is None:
                leavescoord.append(inDel.lc.site)
            if inDel.rc.lc is None:
                leavescoord.append(inDel.rc.site)
        #
            innodeList.remove(inDel)

        # print(np.asarray(bplist))
        print("status structure:")
        for i in range(len(breaknode)):
            string1=str(bplist[i])
            if breaknode[i].pr is None:
                string1=string1+"pr:None"
            else:
                tarbppall = combreakpoint(breaknode[i].pr.lsite, breaknode[i].pr.rsite, ysweep)  # target breakpoint
                string1=string1+"pr:"+str(tarbppall)
            if breaknode[i].lc.lc is None:
                string1=string1+"lc:leaf"+str(breaknode[i].lc.site)
            else:
                tarbppall = combreakpoint(breaknode[i].lc.lsite, breaknode[i].lc.rsite, ysweep)  # target breakpoint
                string1 = string1 + "lc:innode" + str(tarbppall)
            if breaknode[i].rc.lc is None:
                string1 = string1 + "rc:leaf" + str(breaknode[i].rc.site)
            else:
                tarbppall = combreakpoint(breaknode[i].rc.lsite, breaknode[i].rc.rsite, ysweep)  # target breakpoint
                string1 = string1 + "rc:innode" + str(tarbppall)
            string1=string1+"(lsite,rsite)"+str(breaknode[i].lsite)+ str(breaknode[i].rsite)
            print(string1)
        # print("all leaves")
        # print(np.sort(np.asarray(leavescoord), axis=0))
        # print(np.asarray(leavescoord))
    def searchinnode(self,lsite,rsite,ysweep):
        p = self.root  # start from the root of the tree
        tarbp = combreakpoint(lsite, rsite, ysweep)#target breakpoint
        print("tarbp: ",tarbp)
        #compute all the break points
        # self.widoutT(ysweep)
        # print("left to right leaves")
        # self.depoutT(self.root)
        # self.outputTree(ysweep)
        while True:
            bp = combreakpoint(p.lsite, p.rsite, ysweep)
            # print("bp: ",p.lsite, p.rsite,bp)
            if tarbp[0] < bp[0]:
                if (p.lc.lc is None) and (p.lc.rc is None):  # p.lc is a leafnode,alpha FOUND=!
                    return None
                else:  # p.lc is innode
                    p = p.lc
            elif tarbp[0] > bp[0]:
                if (p.rc.lc is None) and (p.rc.rc is None):  # p.rc is a leafnode,alpha FOUND=!
                    return None
                else:  # p.rc is innode
                    p = p.rc
            else:
                break
        return p  # found target innode




