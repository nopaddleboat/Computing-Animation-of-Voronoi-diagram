import numpy as np
import status_swx as status
import queue_swx as qs
import quicksort_swx as quicksort
from numpy import linalg as LA
import copy
CRED = '\033[91m'
CEND = '\033[0m'
CYEL = '\033[93m'#yellow
CG='\033[32m'#green

def test1(infor,key):
    if key is True:
        print(infor)
def findneighborleaf(arc,dir):#dir:1:find left neighbour;2:find right neighbour #arc is a leaf node
    exitnum=6
    testfind=False
    if dir==1:#search the left neighbour leaf in the status tree
        if arc.pr is not None:
            if arc.pr.lc!=arc:#has left leaf
                if arc.pr.lc.lsite is None:
                    test1("line12",testfind)
                    return arc.pr.lc
                else:
                    sp=arc.pr.lc
                    while True:
                        if sp.rc.lsite is None:  # leaf neighbour leaf found
                            break
                        sp = sp.rc
                    test1("line23", testfind)
                    return sp.rc
            else:
                sp = arc.pr  # sliding pointer
                while True:
                    if sp.pr is None:  # at the root
                        test1("line18", testfind)
                        return None
                    elif sp.pr.lc != sp:
                        break
                    else:
                        pass
                    sp = sp.pr
                sp = sp.pr.lc
                if sp.rc is None:  # sp is leafnode
                    test1("line27", testfind)
                    return sp
                while True:
                    if sp.rc.lsite is None:  # leaf neighbour leaf found
                        break
                    sp = sp.rc
                test1("line33", testfind)
                return sp.rc
        else:
            test1("line36",testfind)
            return None#no left neighbour found
    else:#search right neighbour
        if arc.pr is not None:
            if arc.pr.rc != arc:  # has right leaf
                if arc.pr.rc.lsite is None:
                    test1("line44", testfind)
                    return arc.pr.rc
                else:
                    sp=arc.pr.rc
                    while True:
                        if sp.lc.lsite is None:  # leaf neighbour leaf found
                            break
                        sp = sp.lc
                    test1("line56", testfind)
                    return sp.lc
            else:#arc is right child
                sp = arc.pr  # sliding pointer
                while True:
                    if sp.pr is None:  # at the root
                        test1("line48", testfind)
                        return None
                    elif sp.pr.rc != sp:
                        break
                    else:
                        pass
                    sp = sp.pr
                sp = sp.pr.rc
                if sp.lc is None:  # sp is leafnode
                    test1("line60", testfind)
                    return sp
                while True:
                    if sp.lc.lsite is None:  # leaf neighbour leaf found
                        break
                    sp = sp.lc
                test1("line66", testfind)
                return sp.lc
        else:
            test1("line69", testfind)
            return None  # no left neighbour found
def v3v4n3n4(pi,midleaf,thirdleaf,tritype):#compute  v3 v4 n3 n4
    n3=0
    n4=0
    v3=0
    v4=0
    if tritype is True:##<beta alpha pi>
        bp = status.combreakpoint(thirdleaf.site, midleaf.site, pi.site[1])
        v3 = bp[0]
        v4 = bp[1]
        # print("V3V4  1")
        if thirdleaf.site[0] == midleaf.site[0]:
            n4=0
            if thirdleaf.site[1]>midleaf.site[1]:
                n3 = -1
            else:
                n3=1
            # print("V3V4  2")
        elif thirdleaf.site[0] < midleaf.site[0]:  # downward
            n3 = (thirdleaf.site[1] - midleaf.site[1]) / (thirdleaf.site[0] - midleaf.site[0])
            n4 = -1
            # print("V3V4  3")
        else:  # upward
            n3 = -(thirdleaf.site[1] - midleaf.site[1]) / (thirdleaf.site[0] - midleaf.site[0])
            n4 = 1
            # print("V3V4  4")
    else:#<pi alpha eta>
        bp = status.combreakpoint(midleaf.site, thirdleaf.site, pi.site[1])
        v3 = bp[0]
        v4 = bp[1]

        if thirdleaf.site[0] == midleaf.site[0]:
            n4 = 0
            if thirdleaf.site[1] > midleaf.site[1]:
                n3 = 1
            else:
                n3 = -1
            # print("V3V4  1")
        elif thirdleaf.site[0] < midleaf.site[0]:  # upward
            n3 = -(thirdleaf.site[1] - midleaf.site[1]) / (thirdleaf.site[0] - midleaf.site[0])
            n4 = 1
            # print("V3V4  2")
        else:  # downward
            n3 = (thirdleaf.site[1] - midleaf.site[1]) / (thirdleaf.site[0] - midleaf.site[0])
            n4 = -1
            # print("V3V4  3")
    # print("pi,midleaf,thirdleaf,tritype",pi,midleaf,thirdleaf,tritype)
    # print("tritype v3v4n3n4 function v3 v4 n3 n4",tritype, v3, v4, n3, n4)
    return v3,v4,n3,n4
def CheckIsConverge(v1, v2, n1, n2, v3, v4, n3, n4,pi):
    print(CG,"CheckIsConverge v1,v2,n1,n2", v1, v2, n1, n2,CEND)
    print(CG,"CheckIsConverge v3,v4,n3,n4", v3, v4, n3, n4,CEND)
    if abs(n1*n4-n2*n3)<1e-8:
        print(CG,"skew:",n1*n4-n2*n3,CEND)
        isTri=False
        cnt=None
        lowpnt=None
        return isTri, cnt, lowpnt
    isTri=False
    lowpnt = np.array([0.0, 0.0])  # lowest point of triple
    cnt = np.array([0.0, 0.0])  # center of circle
    det=n1*n4-n3*n2
    if det!=0:
        lamda=(v3*n4-v4*n3-v1*n4+v2*n3)/det
        mu=-( v1*n2-v2*n1-v3*n2+v4*n1 )/det
        if lamda<0 or mu<0:
            isTri=False
            # print("CheckIsConverge  1")
        else:
            isTri =True
            cnt[0]=v1+lamda*n1
            cnt[1] = v2 + lamda * n2
            radius=((pi.site[0] - cnt[0]) ** 2 + (pi.site[1] - cnt[1]) ** 2) ** 0.5
            lowpnt[0] = cnt[0]
            lowpnt[1] = cnt[1] - radius
    else:
        isTri=False
        # print("CheckIsConverge  2")
    return isTri, cnt, lowpnt
def checktriple(pi,midleaf,thirdleaf,tritype):#pi is the new arc,return low point coord and disappearing arc
    # midleaf:alpha; thirdleaf:beta or eta
    isTri=False
    lowpnt = np.array([0.0,0.0])#lowest point of triple
    cnt = np.array([0.0,0.0])#center of circle
    if tritype is True:#<beta alpha pi>
        if pi.site[1]==midleaf.site[1]:#pi and alpha have same y coord
            if thirdleaf.site[1]==midleaf.site[1]:
                isTri=False
            else:
                #obtain (v1 v2) (n1 n2)
                v1 = (pi.site[0] + midleaf.site[0]) / 2.0
                #obtain (v3 v4) (n3 n4)
                v3, v4, n3, n4=v3v4n3n4(pi,midleaf,thirdleaf,tritype)
                mu=(v1-v3)/n3
                if mu>=0:
                    isTri=True
                    cnt[0]=v1
                    cnt[1]=v4+mu*n4
                    radius = ((pi.site[0] - cnt[0]) ** 2 + (pi.site[1] - cnt[1]) ** 2) ** 0.5
                    lowpnt[0] = cnt[0]
                    lowpnt[1] = cnt[1] - radius
                else:
                    isTri=False
        else:
            # print("line167")
            # obtain (v1 v2) (n1 n2)
            if pi.site[0]==midleaf.site[0]:
                n1=-1
                n2=0
                v1=pi.site[0]
                v2=(pi.site[1]+midleaf.site[1])/2.0
            elif pi.site[0]<midleaf.site[0]:#upward
                n1 = -(pi.site[1]-midleaf.site[1])/(pi.site[0]-midleaf.site[0])
                n2 = 1
                v1 = pi.site[0]
                v2 = (( v1-midleaf.site[0] )**2+midleaf.site[1]**2-pi.site[1]**2)/2.0/(midleaf.site[1]-pi.site[1])
            else:#downward
                n1 = (pi.site[1] - midleaf.site[1]) / (pi.site[0] - midleaf.site[0])
                n2 = -1
                v1 = pi.site[0]
                v2 = ((v1 - midleaf.site[0]) ** 2 + midleaf.site[1] ** 2 - pi.site[1] ** 2) / 2.0 / (midleaf.site[1] - pi.site[1])

            print("beta v1,v2,n1,n2",v1,v2,n1,n2)
            # obtain (v3 v4) (n3 n4)
            if thirdleaf.site[1] == midleaf.site[1]:
                v3=(thirdleaf.site[0] + midleaf.site[0])/2.0
                lamda=(v3-v1)/n1
                if lamda>=0:#can be zero!!!
                    isTri=True
                    cnt[0]=v3
                    cnt[1]=v2+lamda*n2
                    radius=((pi.site[0]-cnt[0])**2+(pi.site[1]-cnt[1])**2)**0.5
                    lowpnt[0]=cnt[0]
                    lowpnt[1]=cnt[1]-radius
                else:
                    isTri=False
            else:
                v3, v4, n3, n4 = v3v4n3n4(pi,midleaf,thirdleaf,tritype)
                print("beta: v1,v2,n1,n2", v1, v2, n1, n2)
                print("beta: v3,v4,n3,n4", v3, v4, n3, n4)
                isTri, cnt, lowpnt = CheckIsConverge(v1, v2, n1, n2, v3, v4, n3, n4,pi)

    else:#<pi alpha eta>
        if pi.site[1]==midleaf.site[1]:#pi and alpha have same y coord
            if thirdleaf.site[1]==midleaf.site[1]:
                isTri=False
            else:
                #obtain (v1 v2) (n1 n2)
                v1 = (pi.site[0] + midleaf.site[0]) / 2.0
                #obtain (v3 v4) (n3 n4)
                v3, v4, n3, n4=v3v4n3n4(pi,midleaf, thirdleaf, tritype)
                mu=(v1-v3)/n3
                if mu>=0:
                    isTri=True
                    cnt[0]=v1
                    cnt[1]=v4+mu*n4
                    radius = ((pi.site[0] - cnt[0]) ** 2 + (pi.site[1] - cnt[1]) ** 2) ** 0.5
                    lowpnt[0] = cnt[0]
                    lowpnt[1] = cnt[1] - radius
                else:
                    isTri=False
        else:
            # obtain (v1 v2) (n1 n2)
            if pi.site[0]==midleaf.site[0]:
                n1=1
                n2=0
                v1=pi.site[0]
                v2=(pi.site[1]+midleaf.site[1])/2.0
            elif pi.site[0]<midleaf.site[0]:#upward
                n1 = (pi.site[1]-midleaf.site[1])/(pi.site[0]-midleaf.site[0])
                n2 = -1
                v1 = pi.site[0]
                v2 = (( v1-midleaf.site[0] )**2+midleaf.site[1]**2-pi.site[1]**2)/2.0/(midleaf.site[1]-pi.site[1])
            else:#downward
                n1 = -(pi.site[1] - midleaf.site[1]) / (pi.site[0] - midleaf.site[0])
                n2 = 1
                v1 = pi.site[0]
                v2 = ((v1 - midleaf.site[0]) ** 2 + midleaf.site[1] ** 2 - pi.site[1] ** 2) / 2.0 / (midleaf.site[1] - pi.site[1])

            # obtain (v3 v4) (n3 n4)
            if thirdleaf.site[1] == midleaf.site[1]:
                v3=(thirdleaf.site[0] + midleaf.site[0])/2.0
                lamda=(v3-v1)/n1
                if lamda>=0:
                    isTri=True
                    cnt[0]=v3
                    cnt[1]=v2+lamda*n2
                    radius=((pi.site[0]-cnt[0])**2+(pi.site[1]-cnt[1])**2)**0.5
                    lowpnt[0]=cnt[0]
                    lowpnt[1]=cnt[1]-radius
                else:
                    isTri=False
            else:
                v3, v4, n3, n4 = v3v4n3n4(pi,midleaf, thirdleaf, tritype)
                # print("eta: v1,v2,n1,n2", v1, v2, n1, n2)
                # print("eta: v3,v4,n3,n4", v3, v4, n3, n4)
                isTri, cnt, lowpnt = CheckIsConverge(v1, v2, n1, n2, v3, v4, n3, n4,pi)

    disarc = midleaf
    if tritype is True:  # return midleaf leftleaf and right leaf
        print("isTri: ", isTri, "cnt", cnt, "lowpnt", lowpnt)
        print("pi", pi.site, "midleaf", midleaf.site, "thirdleaf", thirdleaf.site)
        return isTri, cnt, lowpnt, disarc, thirdleaf, pi  # if isTri is False,the other two terms are meaningless
    else:
        print("isTri: ", isTri, "cnt", cnt, "lowpnt", lowpnt)
        print("pi", pi.site, "midleaf", midleaf.site, "thirdleaf", thirdleaf.site)
        return isTri, cnt, lowpnt, disarc, pi, thirdleaf

def dirNewedge(leftcurve,rightcurve,nodetype):#direction of new edge:up,down,left,right
    #nodetype: 0 pointing to endpoint, 1 starting from endpoint
    dir=np.array([0.0,0.0,0.0])#direction of the edge
    dir[2]=nodetype
    if leftcurve.site[0]==rightcurve.site[0]:
        if leftcurve.site[1]>rightcurve.site[1]:
            dir[0]=-1
            dir[1]=0
        elif leftcurve.site[1]<rightcurve.site[1]:
            dir[0] = 1
            dir[1] = 0
    elif leftcurve.site[0]<rightcurve.site[0]:
        dir[0]=0
        dir[1]=-1
    else:#leftcurve.site[0]>rightcurve.site[0]
        dir[0] = 0
        dir[1] = 1
    return dir

'''--------------------------------
#deal with degeneracy where several points lie in same/first horizontal line
#return root of the T through recursive operations
#--------------------------------'''
def HandleFirstColine(indexes):#initial index np.array[0 1 2 3]
    global innodelist
    global leaflist
    len1=len(indexes)
    root=None
    if len1==1:#add leafs
        root=innodelist[indexes[0]]
        root.lc=leaflist[indexes[0]]
        leaflist[indexes[0]].pr=root
        root.rc = leaflist[indexes[0]+1]
        leaflist[indexes[0] + 1].pr=root
    elif len1==2:
        root=innodelist[indexes[0]]
        root.rc=innodelist[indexes[1]]
        innodelist[indexes[1]].pr=root
        root.lc=leaflist[indexes[0]]
        leaflist[indexes[0]].pr=root
        root.rc.lc=leaflist[indexes[1]]
        leaflist[indexes[1]].pr=root.rc
        root.rc.rc=leaflist[indexes[1]+1]
        leaflist[indexes[1] + 1].pr=root.rc
    else:
        mid=round(len1/2)
        root=innodelist[indexes[mid-1]]
        root.lc=HandleFirstColine(np.copy(indexes[0:mid-1]))
        root.lc.pr=root
        root.rc=HandleFirstColine(np.copy(indexes[mid:len1]))
        root.rc.pr=root
    return root

def HandleSiteEvent(D,Q,T,event):
    print("HandleSiteEvent: ")
    if T.root is None:#T is empty
        #handle a several point on the first scan line
        global innodelist
        global leaflist
        innodelist = []
        leaflist = []
        Enodelist=[]
        leaf1 = status.leafnode(event.Coord)
        leaflist.append(leaf1)
        leftcoord=np.copy(event.Coord)
        Q.delete(Q.outq)
        print("leftcoord",leftcoord)
        while True:
            if Q.outq.Coord[1] == leftcoord[1]:
                print("event.Coord",Q.outq.Coord[1])
                innode1 = status.innode(np.copy(leftcoord), np.copy(Q.outq.Coord))  # <pj,pi>
                leaf2 = status.leafnode(np.copy(Q.outq.Coord))
                leftcoord = np.copy(Q.outq.Coord)
                Q.delete(Q.outq)
                innodelist.append(innode1)
                leaflist.append(leaf2)
                Enode=D.addEdge(1)
                Enode_Tw = D.addEdge(2)
                Enode.Twin=Enode_Tw
                Enode_Tw.Twin=Enode
                innode1.Enode = Enode
                innode1.dir = dirNewedge(leaflist[len(innodelist)-2], leaf2, 0)
            else:
                break
        for i in range(len(innodelist)):
            print("innodelist[i]: ",i,innodelist[i].lsite,innodelist[i].lsite)
        if len(innodelist)==0:#just add event to the root
            T.setroot(leaflist[0])
            return None
        elif len(innodelist)==1:
            T.setroot(innodelist[0])
            innodelist[0].lc=leaflist[0]
            leaflist[0].pr=innodelist[0]
            innodelist[0].rc = leaflist[1]
            leaflist[1].pr=innodelist[0]
        elif len(innodelist)==2:
            T.setroot(innodelist[0])
            innodelist[0].lc = leaflist[0]
            leaflist[0].pr=innodelist[0]
            innodelist[0].rc = innodelist[1]
            innodelist[1].pr=innodelist[0]
            innodelist[1].lc = leaflist[1]
            leaflist[1].pr=innodelist[1]
            innodelist[1].rc = leaflist[2]
            leaflist[2].pr=innodelist[1]
        else:
            indexes=np.arange(len(innodelist))
            T.setroot(HandleFirstColine(indexes))
        return None
    else:#search in T for alpha which is above event
        if (T.root.lsite is None) and (T.root.rsite is None) :#root is a leafnode;the arc is alpha
            alpha=T.root
        else:#search alpha arc
            p=T.root#start from the root of the tree
            n=1
            while True:
                bp=status.combreakpoint(p.lsite,p.rsite,event.Coord[1])
                if event.Coord[0]<bp[0]:
                    # print("at left")
                    if (p.lc.lc is None) and (p.lc.rc is None):#p.lc is a leafnode,alpha FOUND=!
                        alpha=p.lc#pass leaf to alpha
                        break
                    else:#p.lc is innode
                        p=p.lc
                else:
                    # print("at right")
                    if (p.rc.lc is None) and (p.rc.rc is None):#p.rc is a leafnode,alpha FOUND=!
                        alpha=p.rc#pass leafnode to p
                        break
                    else:#p.rc is innode
                        p=p.rc
                n=n+1

           #check if alpha has circle event
            if alpha.eventnode is not None:
                Q.delete(alpha.eventnode)#false alarm, delete such circle event
                print("Q.delete(alpha.eventnode)")
        print("alpha founded: ",alpha.site)
        if alpha.pr is None:
            print("alpha.pr is None")
        # find left and right arc of alpha<pi alpha eta><beta alpha pi>
        beta = None
        eta = None
        beta = findneighborleaf(alpha, 1)#find left leaf
        if beta is not None:
            print("beta founded: ", beta.site)
        eta = findneighborleaf(alpha, 2)#find right leaf
        if eta is not None:
            print("eta founded: ", eta.site)

        # replace founded leaf p with subtree
        innode1 = status.innode(alpha.site, event.Coord)  # <pj,pi>
        innode2 = status.innode(event.Coord, alpha.site)  # <pi,pj>
        leafL = status.leafnode(alpha.site)
        leafM = status.leafnode(event.Coord)
        leafR = status.leafnode(alpha.site)
        innode1.setvalue(Enode=None, lc=leafL, rc=innode2, pr=None)
        innode2.setvalue(Enode=None, lc=leafM, rc=leafR, pr=innode1)
        leafL.setvalue(parent=innode1, circleevent=None)
        leafM.setvalue(parent=innode2, circleevent=None)
        leafR.setvalue(parent=innode2, circleevent=None)
        T.replaceleaf(alpha, innode1)

        #before replacing alpha with new inserted,find left and right neighbor leaf for later use
        #add records to DCEL:no vertex is needed to added at this point
        Enode=D.addEdge(1)#add origin node
        Enode_Tw = D.addEdge(2)
        Enode.Twin = Enode_Tw
        Enode_Tw.Twin = Enode
        innode1.Enode=Enode#point to DCEL
        innode2.Enode=Enode#point to DCEL
        #check triple consecutive arcs
        #find right triples: arc in the triple at highest position will disappear
        cir_eve = None
        cir_eve1 = None
        trinode11=None#left leaf node
        trinode12 = None#right leaf node
        trinode21 = None#left leaf node
        trinode22 = None#right leaf node
        cnt=None
        cnt1=None
        pi=leafM#make pi be leafnode
        #set direction of the traced edge
        #for innode1
        innode1.dir=dirNewedge(leafL, leafM,0)
        innode2.dir = dirNewedge(leafM, leafR,0)

        if beta is not None:
            isTri, cnt, lowpnt, disarc, trinode11, trinode12 = checktriple(pi,leafR,beta,True)#alpha should be replaced by leafR
            # print(" pi.Coord",pi.site)
            print(" beta:circle event check",isTri)
            if isTri is True:
                cir_eve = qs.qnode(lowpnt)  # circle event
                if trinode12.pr is None:
                    print("trinode12.pr is None")
        if eta is not None:
            isTri, cnt1, lowpnt1, disarc1, trinode21, trinode22 = checktriple(pi,leafL,eta,False)#alpha should be replaced by leafL
            print(" eta:circle event check", isTri)
            if isTri is True:
                cir_eve1 = qs.qnode(lowpnt1)  # circle event

        if cir_eve is not None:
            # print(" beta circle event!",lowpnt)
            Q.insert(cir_eve)
            cir_eve.leaf = leafL
            cir_eve.cnt = cnt

            if (leafL.eventnode is not None) and (leafL.eventnode.Coord[1]<cir_eve.Coord[1]):
                leafL.eventnode = cir_eve
                leafL.trinode1=trinode11
                leafL.trinode2=trinode12
            elif leafL.eventnode is None:
                leafL.eventnode = cir_eve
                leafL.trinode1 = trinode11
                leafL.trinode2 = trinode12
            else:
                pass

            if leafL.trinode2.pr is None:
                print("trinode12.pr is None")
            else:
                print("trinode12.pr is not None:",leafL.trinode2.site)
        if cir_eve1 is not None:
            # print(" eta circle event!",lowpnt1)
            Q.insert(cir_eve1)
            cir_eve1.leaf = leafR
            cir_eve1.cnt = cnt1

            if (leafR.eventnode is not None) and (leafR.eventnode.Coord[1]<cir_eve1.Coord[1]):
                leafR.eventnode = cir_eve1
                leafR.trinode1 =trinode21
                leafR.trinode2 =trinode22
            elif leafR.eventnode is None:
                leafR.eventnode = cir_eve1
                leafR.trinode1 = trinode21
                leafR.trinode2 = trinode22
            else:
                pass

def del_cir_Q(Q,cir_eve,gama):#delete circle events involving alpha from Q
    print("del_cir_Q: ")
    if cir_eve is not None:
        if cir_eve.eventnode is not None:
            print("cir_eve.trinode1 2,gamma.site",cir_eve.trinode1.site,cir_eve.trinode2.site,gama.site)
            # just use trinode1's site, it may be deleted
            #the following equivalent to (cir_eve.trinode1== gama) or (cir_eve.trinode2== gama)
            if LA.norm(cir_eve.trinode1.site-gama.site)==0 or LA.norm(cir_eve.trinode2.site-gama.site)==0:
                Q.delete(cir_eve.eventnode)
                print(" delete circle!")

def DirVec(Lleaf,Rleaf):#compute the direction vector of two sites
    dir = np.array([0.0, 0.0])  # direction of the edge
    dir1=Lleaf-Rleaf
    if Lleaf[0] == Rleaf[0]:
        if Lleaf[1] > Rleaf[1]:
            dir[0] = -1
            dir[1] = 0
        elif Lleaf[1] < Rleaf[1]:
            dir[0] = 1
            dir[1] = 0
    elif Lleaf[0] < Rleaf[0]:
        dir[0] = dir1[1]/dir1[0]
        dir[1] = -1
    else:  # leftcurve.site[0]>rightcurve.site[0]
        dir[0] = -dir1[1]/dir1[0]
        dir[1] = 1
    return dir


def HandleCircleEvent(D,Q,T,event):
    print("HandleCircleEvent: ", event.leaf.site)
    ysweep=event.Coord[1]
    cnt=event.cnt#center of the circle
    Lleaf1=None
    Lleaf2=None
    Rleaf1=None
    Rleaf2=None

    # print("event.leaf event.leaf.trinode1 2 ",event.leaf.site,event.leaf.trinode1.site,event.leaf.trinode2.site)
    Lleaf1=findneighborleaf(event.leaf, 1)
    Rleaf1=findneighborleaf(event.leaf, 2)#don't use event.leaf.trinode2, which may be destroyed in deleteleaf
    if Lleaf1 is not None and Rleaf1 is not None:
        print("Lleaf1,Rleaf1 ", Lleaf1.site, " ", Rleaf1.site)
    else:
        print("Lleaf1 is None OR Rleaf1 is None")

    if Lleaf1 is not None:
        Lleaf2 = findneighborleaf(Lleaf1, 1)
        if Lleaf2 is not None:
            print("Lleaf1,Lleaf2 ",Lleaf1.site," ",Lleaf2.site)
    else:
        print("Error line310, Lleaf1 should exist!")
    if Rleaf1 is not None:
        Rleaf2 = findneighborleaf(Rleaf1, 2)
        if Rleaf2 is not None:
            print("Rleaf1,Rleaf2 ", Rleaf1.site, " ", Rleaf2.site)
    else:
        print("Error line315, Rleaf1 should exist!")

    # step1.1 delete the leaf gama
    # step1.2 delete relevant circle events from Q
    del_cir_Q(Q, Lleaf1, event.leaf)
    del_cir_Q(Q, Lleaf2, event.leaf)
    del_cir_Q(Q, Rleaf1, event.leaf)
    del_cir_Q(Q, Rleaf2, event.leaf)

    #step2:before delete leaf, add edges
    #step2.1 add center as vertex
    Vnode = D.addVertex()#voronoi diagram node
    Vnode.setdata(cnt, IncidentEdge=None)
    #add new break point
    #add records to relevant edges
    innode_addE2 = event.leaf.pr  # innode edge to be deleted
    flagleft=True
    # print("event.leaf.pr.lsite",event.leaf.pr.lsite)
    if event.leaf.site[0]==event.leaf.pr.lsite[0] and event.leaf.site[1]==event.leaf.pr.lsite[1]:#leaf is left of breakpoint
        # print("LLL:",Lleaf1.site,event.leaf.site,event.Coord[1])
        innode_addE1=T.searchinnode(Lleaf1.site,event.leaf.site,event.Coord[1])#innode edge to be updated
        print("Lleaf1.site,event.leaf.site",Lleaf1.site,event.leaf.site)
        innode_addE1.rsite = Rleaf1.site
        innode_addE1.dir = dirNewedge(Lleaf1,Rleaf1,1)  # ray point downward,used to calculate intersection of ray and bounding box
        dir_vec=DirVec(Lleaf1.site,event.leaf.site)
        dir_vec1 = DirVec(event.leaf.site,Rleaf1.site)#for innode2
    else:#leaf is right of breakpoint
        flagleft=False
        innode_addE1 = T.searchinnode( event.leaf.site,Rleaf1.site,event.Coord[1])  # innode edge to be updated
        innode_addE1.lsite = Lleaf1.site
        innode_addE1.dir = dirNewedge(Lleaf1,Rleaf1,1)# ray point downward,used to calculate intersection of ray and bounding box
        dir_vec = DirVec(event.leaf.site,Rleaf1.site)
        dir_vec1 = DirVec(Lleaf1.site,event.leaf.site)  # for innode2
    print("dirvec innode1 innode2", dir_vec, dir_vec1,innode_addE2.lsite,innode_addE2.rsite)

    if innode_addE1.Enode.Origin is None:  # add Vnode as origin
        innode_addE1.Enode.Origin = Vnode
        #rotate dir_vec anticlockwise
        print("innode_addE1.Enode.Origin is None")
        dir_vec=np.array([-dir_vec[1],dir_vec[0]])

        #add face information
        if flagleft is True:
            if -np.dot(dir_vec,event.leaf.site[0:2]-Vnode.Coord[0:2])>0:#start from vnode, reverse dir_vec;  face corresponds to eventleaf
                innode_addE1.Enode.IncidentFace=D.Flist[int(event.leaf.site[2])]
                innode_addE1.Enode.Twin.IncidentFace=D.Flist[int(Lleaf1.site[2])]

            else:
                innode_addE1.Enode.IncidentFace=D.Flist[int(Lleaf1.site[2])]
                innode_addE1.Enode.Twin.IncidentFace=D.Flist[int(event.leaf.site[2])]
        else:
            if -np.dot(dir_vec,event.leaf.site[0:2]-Vnode.Coord[0:2])>0:#start from vnode, reverse dir_vec;  face corresponds to eventleaf
                innode_addE1.Enode.IncidentFace=D.Flist[int(event.leaf.site[2])]
                innode_addE1.Enode.Twin.IncidentFace=D.Flist[int(Rleaf1.site[2])]
            else:
                innode_addE1.Enode.IncidentFace=D.Flist[int(Rleaf1.site[2])]
                innode_addE1.Enode.Twin.IncidentFace=D.Flist[int(event.leaf.site[2])]
    else:  # create twin edge, and add pointers between existing and twin edge
        innode_addE1.Enode.Twin.setdata(Origin=Vnode, Twin=None, IncidentFace=None, Next=None, Prev=None)
    if innode_addE1.Enode.Origin==Vnode:
        node11=innode_addE1.Enode
    else:
        node11=innode_addE1.Enode.Twin


    dir_vec2 = DirVec(Lleaf1.site, Rleaf1.site)  # for new breakpoint
    Enode_newbp = D.addEdge(1)  # node for the new break point,store it in internal node!edge for new breakpoint
    Enode_tw = D.addEdge(2)  # create twin node to store voroni vertex
    innode_addE1.Enode = Enode_newbp
    print("dir_vec2, Lleaf1.site[0:2]", dir_vec2)
    dir_vec2 = np.array([-dir_vec2[1], dir_vec2[0]])#rotate 90 degrees before dot product

    if np.dot(dir_vec2, Lleaf1.site[0:2]-Vnode.Coord[0:2]) > 0:

        Enode_newbp.setdata(Origin=Vnode, Twin=Enode_tw, IncidentFace=D.Flist[int(Lleaf1.site[2])], Next=None, Prev=None)
        Enode_tw.setdata(Origin=None, Twin=Enode_newbp, IncidentFace=D.Flist[int(Rleaf1.site[2])], Next=None, Prev=None)
    else:
        Enode_newbp.setdata(Origin=Vnode, Twin=Enode_tw, IncidentFace=D.Flist[int(Rleaf1.site[2])], Next=None, Prev=None)
        Enode_tw.setdata(Origin=None, Twin=Enode_newbp, IncidentFace=D.Flist[int(Lleaf1.site[2])], Next=None, Prev=None)
    node31 = Enode_newbp


    if innode_addE2 is not None:
        if innode_addE2.Enode.Origin is None:#add Vnode as origin
            dir_vec1 = np.array([-dir_vec1[1], dir_vec1[0]])
            innode_addE2.Enode.Origin=Vnode
            # add face information (event.leaf.site,Rleaf1)
            print("innode_addE2 up",-np.dot(dir_vec1, event.leaf.site[0:2]-Vnode.Coord[0:2]),dir_vec1,Vnode.Coord[0:2])
            if flagleft is True:
                if -np.dot(dir_vec1, event.leaf.site[0:2]-Vnode.Coord[0:2]) > 0:  # start from vnode;  face corresponds to eventleaf

                    innode_addE2.Enode.IncidentFace=D.Flist[int(event.leaf.site[2])]
                    innode_addE2.Enode.Twin.IncidentFace=D.Flist[int(Rleaf1.site[2])]
                else:
                    innode_addE2.Enode.IncidentFace=D.Flist[int(Rleaf1.site[2])]
                    innode_addE2.Enode.Twin.IncidentFace=D.Flist[int(event.leaf.site[2])]
            else:
                if -np.dot(dir_vec1, event.leaf.site[0:2]-Vnode.Coord[0:2]) > 0:  # start from vnode, reverse dir_vec;  face corresponds to eventleaf
                    innode_addE2.Enode.IncidentFace=D.Flist[int(event.leaf.site[2])]
                    innode_addE2.Enode.Twin.IncidentFace=D.Flist[int(Lleaf1.site[2])]
                else:
                    innode_addE2.Enode.IncidentFace=D.Flist[int(Lleaf1.site[2])]
                    innode_addE2.Enode.Twin.IncidentFace=D.Flist[int(event.leaf.site[2])]
        else:#create twin edge, and add pointers between existing and twin edge
            innode_addE2.Enode.Twin.Origin=Vnode
    else:
        print("Error, line343 innode_addE1 should exist!")
    if innode_addE2.Enode.Origin==Vnode:
        node21 = innode_addE2.Enode
    else:
        node21=innode_addE2.Enode.Twin

    #set pre and next of nodes
    #update node31 information
    print("node11",node11.EdgeName,node11.IncidentFace.Face,node11.Twin.IncidentFace.Face)
    print("node21", node21.EdgeName, node21.IncidentFace.Face,node21.Twin.IncidentFace.Face)
    print("node31", node31.EdgeName, node31.IncidentFace.Face,node31.Twin.IncidentFace.Face)
    if node31.IncidentFace==node11.Twin.IncidentFace:
        node31.Prev=node11.Twin
    else:
        node31.Prev=node21.Twin
    if node31.Twin.IncidentFace == node11.IncidentFace:
        node31.Twin.Next = node11
    else:
        node31.Twin.Next = node21
    # update node11 information
    if node11.IncidentFace==node31.Twin.IncidentFace:
        node11.Prev=node31.Twin
    else:
        node11.Prev=node21.Twin
    if node11.Twin.IncidentFace == node31.IncidentFace:
        node11.Twin.Next = node31
    else:
        node11.Twin.Next = node21
    # update node21 information
    if node21.IncidentFace==node11.Twin.IncidentFace:
        node21.Prev=node11.Twin
    else:
        node21.Prev=node31.Twin
    if node21.Twin.IncidentFace == node11.IncidentFace:
        node21.Twin.Next = node11
    else:
        node21.Twin.Next = node31

    T.deleteleaf(event.leaf)
    #update tuples in internal node:at most two internal nodes are affected
    #step 3:check new circle triples
    cir_eve3 = None
    cir_eve4 = None
    trinode31 = None
    trinode32 = None
    trinode41 = None
    trinode42 = None
    cnt3 = None
    cnt4 = None
    isTri=None

    if (Lleaf2 is not None) and (Lleaf1 is not None) and (Rleaf1 is not None):
        min=-1#default beta1 lowest
        if Lleaf2.site[1]>Lleaf1.site[1]:
            min=0
            if Lleaf1.site[1]>Rleaf1.site[1]:
                min=1
        else:
            if Lleaf2.site[1]>Rleaf1.site[1]:
                min=1
        if min==0:
            pass
        elif min==-1:#equivalent to<pi alpha eta>
            isTri,cnt3, lowpnt,  disarc, trinode31, trinode32 = checktriple(Lleaf2, Lleaf1, Rleaf1,False)
            if isTri is True:
                cir_eve3 = qs.qnode(lowpnt)  # circle event found
        else:#equivalent to<beta alpha pi>
            isTri, cnt3, lowpnt, disarc, trinode31, trinode32 = checktriple(Rleaf1, Lleaf1, Lleaf2, True)
            if isTri is True:
                cir_eve3 = qs.qnode(lowpnt)  # circle event found
    if (Lleaf1 is not None) and (Rleaf1 is not None) and (Rleaf2 is not None):
        min = -1  # default beta1 lowest
        if Lleaf1.site[1] > Rleaf1.site[1]:
            min = 0
            if Rleaf1.site[1] > Rleaf2.site[1]:
                min = 1
        else:
            if Lleaf1.site[1] > Rleaf2.site[1]:
                min = 1
        if min == 0:
            pass
        elif min == -1:  # equivalent to<pi alpha eta>
            isTri, cnt4, lowpnt1, disarc1, trinode41, trinode42 = checktriple(Lleaf1, Rleaf1, Rleaf2, False)
            if isTri is True:
                cir_eve4 = qs.qnode(lowpnt1)  # circle event found
        else:  # equivalent to<beta alpha pi>
            isTri, cnt4, lowpnt1, disarc1, trinode41, trinode42 = checktriple(Rleaf2, Rleaf1, Lleaf1, True)
            if isTri is True:
                cir_eve4 = qs.qnode(lowpnt1)  # circle event found

    if cir_eve3 is not None:
        print("cir_eve3!")
        Q.insert(cir_eve3)
        cir_eve3.leaf = disarc
        cir_eve3.cnt = cnt3
        if disarc.pr is None:
            print("354 is None")
        if (disarc.eventnode is not None) and (disarc.eventnode.Coord[1] < cir_eve3.Coord[1]):
            disarc.eventnode = cir_eve3
            disarc.trinode1 = trinode31
            disarc.trinode2 = trinode32
        elif disarc.eventnode is None:
            disarc.eventnode = cir_eve3
            disarc.trinode1 = trinode31
            disarc.trinode2 = trinode32
        else:
            pass
    if cir_eve4 is not None:
        print("cir_eve4!")
        Q.insert(cir_eve4)
        cir_eve4.leaf = disarc1
        cir_eve4.cnt =cnt4
        if (disarc1.eventnode is not None) and (disarc1.eventnode.Coord[1] < cir_eve4.Coord[1]):
            disarc1.eventnode = cir_eve4
            disarc1.trinode1 = trinode41
            disarc1.trinode2 = trinode42
        elif disarc1.eventnode is None:
            disarc1.eventnode = cir_eve4
            disarc1.trinode1 = trinode41
            disarc1.trinode2 = trinode42
        else:
            pass
    # if T.root.rc.lc is not None:
        # T.root=status.leftRotate(T.root,ysweep)
    T.root = status.splay(T.root,innode_addE1,ysweep)

class node_qin:#node for storing innode
    def __init__(self,qnode):
        self.pre=None
        self.next=None
        self.qnode=qnode

class queue_innode:
    def __init__(self,node_q):
        self.head=node_q
        self.tail=node_q
    def insert(self,node_q):#when insert a node, the queue should not be empty
        node_q.next=self.tail
        self.tail.pre=node_q
        self.tail=node_q
    def popout(self):
        self.head=self.head.pre
class boxnode:
    def __init__(self,interpnt,Enode):
        self.interpnt = interpnt
        self.Enode = Enode
def AddBoundingBox(T,D,sites,Fzero,sortedlist):#Fzero bounding box
    #step 1 find bounds:Leftb,Rightb,Topb,Bottomb
    bdsize=5
    Leftb=np.amin(sites[:, 0])-bdsize
    Rightb= np.amax(sites[:, 0])+bdsize
    Bottomb= np.amin(sites[:, 1])-bdsize
    Topb=np.amax(sites[:, 1])+bdsize
    boxcoords=np.array([[Leftb,Topb],[Rightb,Topb],[Rightb,Bottomb],[Leftb,Bottomb]])
    #add 4 vertices and edges
    nodesList=[]
    BoxlinkL=[]#left
    BoxlinkR=[]#right
    BoxlinkT=[]#top
    BoxlinkB=[]#bottom
    interpnt=np.array([0.0,0.0])

    for i in range(0,4):
        nodesList.append( D.addVertex() )
        nodesList[i].setdata(boxcoords[i,:], IncidentEdge=None)

    #step 2 scan the DCEL , complete the DCEL
    #search T for internal nodes, which represents rays and will intersect with bounding box
    #traverse T,find all innodes,use a queue to store
    # compute intersection point i.e., vertex to be added
    p=D.EdgeHead
    while True:
        if p is None:
            break
        if p.EdgeName[1]==1 and p.Twin.Origin is None:#rays found
            lsite=sortedlist[p.Twin.IncidentFace.Face-1,:]
            rsite= sortedlist[p.IncidentFace.Face-1,:]
            cx = (lsite[0] + rsite[0] ) / 2.0  # center point of site1 and site2 <cx,cy>
            cy = (lsite[1] + rsite[1]) / 2.0
            direction=DirVec(lsite,rsite)
            if direction[1] == -1:  # point downward
                if lsite[1] == rsite[1]:  # vertical line
                    interpnt[0] = (lsite[0] + rsite[0]) / 2.0
                    interpnt[1] = Bottomb
                    BoxlinkB.append(boxnode(np.copy(interpnt), p))
                    print("BoxlinkB1")
                else:
                    sk = (rsite[1] - lsite[1]) / (rsite[0] - lsite[0])
                    interpnt[0] = -sk * (Bottomb - cy) + cx
                    if interpnt[0] > Rightb:
                        interpnt[0] = Rightb
                        interpnt[1] = -1.0 / sk * (Rightb - cx) + cy
                        bnode = boxnode(np.copy(interpnt), p)
                        BoxlinkR.append(bnode)
                        print("BoxlinkR2")
                    elif interpnt[0] < Leftb:
                        interpnt[0] = Leftb
                        interpnt[1] = -1.0 / sk * (Leftb - cx) + cy
                        BoxlinkL.append(boxnode(np.copy(interpnt), p))
                        print("BoxlinkL3")
                    else:
                        interpnt[1] = Bottomb
                        BoxlinkB.append(boxnode(np.copy(interpnt), p))
                        print("BoxlinkB4")
            elif direction[1] == 1:  # point upward
                if lsite[1] == rsite[1]:  # vertical line
                    interpnt[0] = (lsite[0] + rsite[0]) / 2.0
                    interpnt[1] = Topb
                    BoxlinkT.append(boxnode(np.copy(interpnt), p))
                    print("BoxlinkT5")
                else:
                    sk = (rsite[1] - lsite[1]) / (rsite[0] - lsite[0])
                    interpnt[0] = -sk * (Topb - cy) + cx
                    if interpnt[0] > Rightb:
                        interpnt[0] = Rightb
                        interpnt[1] = -1.0 / sk * (Rightb - cx) + cy
                        bnode = boxnode(np.copy(interpnt), p)
                        BoxlinkR.append(bnode)
                        print("BoxlinkR6")
                    elif interpnt[0] < Leftb:
                        interpnt[0] = Leftb
                        interpnt[1] = -1.0 / sk * (Leftb - cx) + cy
                        BoxlinkL.append(boxnode(np.copy(interpnt), p))
                        print("BoxlinkL7")
                    else:
                        interpnt[1] = Topb
                        BoxlinkT.append(boxnode(np.copy(interpnt), p))
                        print("BoxlinkT8")
            else:  # horizontal
                interpnt[1] = cy
                if direction[0] == 1:  # point to right
                    interpnt[0] = Rightb
                    bnode = boxnode(np.copy(interpnt), p)
                    BoxlinkR.append(bnode)
                    print("BoxlinkR9")
                else:  # point to left
                    interpnt[0] = Leftb
                    BoxlinkL.append(boxnode(np.copy(interpnt), p))
                    print("BoxlinkL10")


            if p.EdgeName[1] == 1:
                if p.Origin is not None:
                    out = "v" + str(p.EdgeName[0]) + " Origin" + str(p.Origin.Coord)
                else:
                    out = "v" + "No origin"
                if p.Twin is not None:
                    if p.Twin.Origin is not None:
                        out = out + "Twin-Origin" + str(p.Twin.Origin.Coord)
                    else:
                        out = out + "Twin(No origin)"
                print(out, "dir:", direction,p.IncidentFace.Face,p.Twin.IncidentFace.Face)
        p = p.next

        #add vertex and twin edge to D structure
        # Vnode = D.addVertex()  # voronoi diagram node
        # Vnode.setdata(interpnt, IncidentEdge=None)
    # print(D.edgenum)

    print("len(BottomT): ",len(BoxlinkT))
    print("len(BottomB): ", len(BoxlinkB))
    print("len(BottomL): ", len(BoxlinkL))
    print("len(BottomR): ", len(BoxlinkR))

    BoxlinkT=quicksort.quicksortBList(BoxlinkT, 1)
    BoxlinkB=quicksort.quicksortBList(BoxlinkB, 2)
    BoxlinkL=quicksort.quicksortBList(BoxlinkL, 3)
    BoxlinkR=quicksort.quicksortBList(BoxlinkR, 4)

    # Fnode = D.addFace(i + 1)
    B0=[]
    B1=[]
    B2=[]
    B3=[]
    B0.append(0)
    B1.append(0)
    B2.append(0)
    B3.append(0)
    B4links=BoxlinkT+B1+BoxlinkR[::-1]+B2+BoxlinkB[::-1]+B3+BoxlinkL+B0
    print("boxcoords",boxcoords)
    print("B4links[i].interpnt")
    for i in range(len(B4links)):
        if B4links[i]!=0:
            print(B4links[i].interpnt)
        else:
            print(B4links[i])
    indexcorner=0
    firstinter=0
    for i in range(4):
        if B4links[i]!=0:
            firstinter=i
            break
        else:
            indexcorner=indexcorner+1
    addedcorner=0

    for i in range(firstinter,len(B4links)):#start from
        print("i",i,firstinter)
        #print("firstinter indexcorner addedcorner", firstinter, indexcorner,addedcorner)
        if i==firstinter:
            E01 = D.addEdge(1)
            E02 = D.addEdge(2)
            print("case1:",E01.EdgeName,indexcorner)
            N1 = D.addVertex()
            N1.Coord = np.copy(np.array([B4links[i].interpnt[0], B4links[i].interpnt[1], 0]))
            B4links[i].Enode.Twin.Origin = N1
            B4links[i].Enode.Next = E01
            E01.setdata(Origin=N1, Twin=E02,
                        IncidentFace=B4links[i].Enode.IncidentFace,
                        Next=None, Prev=B4links[i].Enode)
            E02.setdata(Origin=nodesList[indexcorner], Twin=E01,
                        IncidentFace=Fzero,
                        Next=None, Prev=None)
            Epre = E01
            Npre = N1
            Eface=B4links[i].Enode.IncidentFace
        else:
            if B4links[i]==0:#corner node
                E11 = D.addEdge(1)
                E12 = D.addEdge(2)
                print("case2:", E11.EdgeName,E12.EdgeName)
                Epre.Twin.Next=E12
                if B4links[i-1]==0:#pre is corner node
                    Epre.Prev = E11
                    E11.setdata(Origin=nodesList[(indexcorner+addedcorner+1)%4], Twin=E12,
                            IncidentFace=Eface,
                            Next=Epre, Prev=None)
                else:
                    B4links[i-1].Enode.Twin.Prev=E11
                    E11.setdata(Origin=nodesList[(indexcorner+addedcorner+1)%4], Twin=E12,
                                IncidentFace=Epre.Prev.Twin.IncidentFace,
                                Next=Epre.Prev.Twin, Prev=None)
                E12.setdata(Origin=Npre, Twin=E11,
                            IncidentFace=Fzero,
                            Next=None, Prev=Epre.Twin)
                addedcorner=addedcorner+1
                Epre = E11
                Npre = E11.Origin
                Eface = E11.IncidentFace
            else:#not corner node
                E11 = D.addEdge(1)
                E12 = D.addEdge(2)
                N1 = D.addVertex()
                print("Case3: ", E11.EdgeName,B4links[i].Enode.EdgeName)
                N1.Coord = np.copy(np.array([B4links[i].interpnt[0], B4links[i].interpnt[1], 0]))
                Epre.Twin.Next = E12
                B4links[i].Enode.Next = E11

                B4links[i].Enode.Twin.Origin = N1

                if B4links[i - 1] == 0:  # pre is corner node
                    Epre.Prev = E11
                    E11.setdata(Origin=N1, Twin=E12,
                                IncidentFace=B4links[i].Enode.IncidentFace,
                                Next=Epre, Prev=B4links[i].Enode)
                else:
                    B4links[i - 1].Enode.Twin.Prev = E11
                    E11.setdata(Origin=N1, Twin=E12,
                                IncidentFace=B4links[i].Enode.IncidentFace,
                                Next=B4links[i-1].Enode.Twin, Prev=B4links[i].Enode)
                E12.setdata(Origin=Npre, Twin=E11,
                                IncidentFace=Fzero,
                                Next=None, Prev=Epre.Twin)
                Epre = E11
                Npre = N1
                Eface = E11.IncidentFace
        p11 = D.EdgeHead
        while True:
            if p11 is None:
                break
            if p11.EdgeName[0] == 36 and p11.EdgeName[1]==1:
                if p11.Next is not None:
                    print("p11.EdgeName[0]==36", p11.Next.EdgeName)
                else:
                    print("p11.Next is  None!")
            p11 = p11.next

    if addedcorner==1:
        E21 = D.addEdge(1)#1 to 0
        E22 = D.addEdge(2)
        E31 = D.addEdge(1)#2 to 1
        E32 = D.addEdge(2)
        E41 = D.addEdge(1)#3 to 2
        E42 = D.addEdge(2)
        E21.setdata(Origin=nodesList[1], Twin=E22,
                    IncidentFace=E01.IncidentFace,
                    Next=Epre, Prev=E31)
        E22.setdata(Origin=nodesList[0], Twin=E21,
                    IncidentFace=Fzero,
                    Next=E32, Prev=Epre.Twin)
        E31.setdata(Origin=nodesList[2], Twin=E32,
                    IncidentFace=E01.IncidentFace,
                    Next=E21, Prev=E41)
        E32.setdata(Origin=nodesList[1], Twin=E31,
                    IncidentFace=Fzero,
                    Next=E42, Prev=E42)
        E41.setdata(Origin=nodesList[3], Twin=E42,
                    IncidentFace=E01.IncidentFace,
                    Next=E01, Prev=E31)
        E42.setdata(Origin=nodesList[2], Twin=E41,
                    IncidentFace=Fzero,
                    Next=E01.Twin, Prev=E32)
        E01.Next=E41
        E02.Prev=E42
        Epre.Twin.Next=E22
        Epre.Prev =E21
        print("case1")
    elif addedcorner==2:
        print("case2")
        E21 = D.addEdge(1)  # 1 to 0
        E22 = D.addEdge(2)
        E31 = D.addEdge(1)  # 2 to 1
        E32 = D.addEdge(2)
        E21.setdata(Origin=nodesList[1], Twin=E22,
                    IncidentFace=E01.IncidentFace,
                    Next=Epre, Prev=E31)
        E22.setdata(Origin=nodesList[0], Twin=E21,
                    IncidentFace=Fzero,
                    Next=E32, Prev=Epre.Twin)
        E31.setdata(Origin=nodesList[2], Twin=E32,
                    IncidentFace=E01.IncidentFace,
                    Next=E01, Prev=E21)
        E32.setdata(Origin=nodesList[1], Twin=E31,
                    IncidentFace=Fzero,
                    Next=E01.Twin, Prev=E22)
        E01.Next = E31
        E02.Prev = E32
        Epre.Twin.Next = E22
        Epre.Prev = E21
    elif addedcorner==3:
        print("case3")
        E21 = D.addEdge(1)  # 1 to 0
        E22 = D.addEdge(2)
        E21.setdata(Origin=nodesList[1], Twin=E22,
                    IncidentFace=E01.IncidentFace,
                    Next=Epre, Prev=E01)
        E22.setdata(Origin=nodesList[0], Twin=E21,
                    IncidentFace=Fzero,
                    Next=E01.Twin, Prev=Epre.Twin)
        E01.Next = E21
        E02.Prev = E22
        Epre.Twin.Next = E22
        Epre.Prev = E21
    else:# addedcorner==4:
        print("case4")
        E01.Next = Epre
        E02.Prev = Epre.Twin
        Epre.Twin.Next = E02
        Epre.Prev = E01
    # print("E01",E01.Origin.Coord)
    return boxcoords



