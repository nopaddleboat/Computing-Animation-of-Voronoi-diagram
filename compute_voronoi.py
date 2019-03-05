import numpy as np
import re
import quicksort_swx as qs
import queue_swx
import status_swx
import DCEL
import handle_event
import matplotlib.pyplot as plt

CRED = '\033[91m'
CEND = '\033[0m'
CYEL = '\033[93m'#yellow
CHIGHBLUE='\033[0;94m'

sitesfile=open("sites.txt","r")
nums=[]
for line in sitesfile:
    nums=nums+[float(s) for s in re.findall(r'-?\d+\.?\d*',line)]#use regular expression
    #print(nums,len(nums))
sitesfile.close()
num_sites=int( len(nums)/2 )
sites=np.zeros( (num_sites,2) )#be careful! Don't forget the ()
for i in range(0,num_sites):
    sites[i]=nums[i*2:i*2+2]
print("sites coordinates:")
print(sites)

#initiailize Q
Q=queue_swx.queue_v(sites)
#Q.printall()
T=status_swx.statustree()
D=DCEL.DCEL()
print("Voronoi diagram computation process: ")

#create len(sites) faces to DCEL and store it to DCEL
Facelist=[]
for i in range(0,num_sites):
    Fnode=D.addFace(i+1)
    Facelist.append(Fnode)
D.Flist=Facelist
Fzero=D.addFace(0)
# print(D.Flist)
################test
# site1=np.array([1., 1.])
# site2=np.array([3.666, 5.444])
# site3=np.array([14. ,  4.5])
# handle_event.checktriple_test(site1,site2,site3)
##################
allVorVer=np.zeros([100,2])
numver=0
n=1
while Q.outq is not None:
    event=Q.outq
    if n!=1:
        # print(CYEL, "B$ delete Q", CEND)
        # Q.printall()
        # print(CYEL, "QEND", CEND)

        Q.delete(Q.outq)#delete processed data
        # if n==5:
        # print(CYEL,"After delete Q",CEND)
        # Q.printall()
        # print(CYEL, "QEND", CEND)
        # if Q.outq is not None:
        #     print(CYEL, n, "Q.outq", Q.outq.Coord,CEND)
        # else:
        #     print(CYEL, n, "Q.outq is NONE!!!!!",CEND)
    if event.leaf is None:#event is a site event
        print(CYEL,n, "th iteration: ", event.Coord,CEND)
        # T.depoutT(T.root)
        handle_event.HandleSiteEvent(D,Q,T,event)
    else:#event is a circle event
        allVorVer[numver,:]=event.cnt
        numver=numver+1
        print(CRED, n, "th iteration: ", event.Coord, "cnt:", event.cnt, CEND)
        handle_event.HandleCircleEvent(D,Q,T,event)
    n=n+1
    if n==1:
        break
print("Detected Voronoi Verteces: ")
for i in range(numver):
    print(allVorVer[i,:])
if numver==0:
    print(CHIGHBLUE,"All sites are colinear!",CEND)
if numver>0:
    box=handle_event.AddBoundingBox(T,D,allVorVer[0:numver,:],Fzero,Q.sortedlist)
    #plot the voronoi diagram from DCEL:
    D.mergeVertices()
    D.completeDCEL()
    print("---------------------------------")
    D.PrintEdges()
    print("---------------------------------")
    print("**************Voronoi Diagram**********************")
    Edgecoords=D.printall()
    Edgecoords,edgelist = np.asarray(Edgecoords)
    plt.axis([box[0,0]-2, box[1,0]+2, box[2,1]-2, box[0,1]+2])#[Leftb,Topb],[Rightb,Topb],[Rightb,Bottomb],[Leftb,Bottomb]
    Lcoords=np.zeros([3,4])
    for i in range(len(Edgecoords)):
        # k=i
        # i=i+25
        # if i==42:#test the face information
        #     te=edgelist[i]
        #     tepre=edgelist[i].Prev
        #     tenext=edgelist[i].Next
        #     Lcoords[0,0:2]=te.Origin.Coord[0:2]
        #     Lcoords[0,2:4] = te.Twin.Origin.Coord[0:2]
        #     Lcoords[1, 0:2] = tepre.Origin.Coord[0:2]
        #     Lcoords[1, 2:4] = tepre.Twin.Origin.Coord[0:2]
        #     Lcoords[2, 0:2] = tenext.Origin.Coord[0:2]
        #     Lcoords[2, 2:4] = tenext.Twin.Origin.Coord[0:2]
        #     plt.plot([Lcoords[0,0],Lcoords[0,2]], [Lcoords[0,1],Lcoords[0,3]], '--',color='red',linewidth=4.0)
        #     plt.plot(Lcoords[0,2],Lcoords[0,3],'*r',markersize=16)
        #     plt.plot([Lcoords[1, 0], Lcoords[1, 2]], [Lcoords[1, 1], Lcoords[1, 3]],'--', color='black',linewidth=4.0)
        #     plt.plot([Lcoords[2, 0], Lcoords[2, 2]], [Lcoords[2, 1], Lcoords[2, 3]],linestyle='--', marker='o', color='b',linewidth=4.0)
        #     print(Lcoords,te.EdgeName,tepre.EdgeName,tenext.EdgeName)
        # i=k
        plt.plot(Edgecoords[i][0,[0,2]], Edgecoords[i][0,[1,3]],color = 'green')
        plt.plot(Edgecoords[i][0,0], Edgecoords[i][0,1], marker='o', markersize=6,color = 'green')
    for i in range(len(sites)):
        plt.plot(Q.sortedlist[i, 0], Q.sortedlist[i, 1], marker='+', markersize=5,color = 'green')
        plt.text(Q.sortedlist[i, 0], Q.sortedlist[i, 1], str(Q.sortedlist[i, 2]+1), fontsize=12)


    Delau=DCEL.VorToDel(D,box,Q.sortedlist)

    p=Delau.EdgeHead
    Dnode=np.zeros([1,4])
    n=0
    while True:

        if p is None:
            break
        if p.EdgeName[1]==1:
            n = n + 1
            Dnode[0,0:2]=np.copy(p.Origin.Coord[0:2])
            Dnode[0, 2:4] = np.copy(p.Twin.Origin.Coord[0:2])
            # print('n:',n,Dnode,p.EdgeName)
            plt.plot(Dnode[0,[0,2]], Dnode[0,[1,3]], color='red')
        p=p.next

    Delau.completeDCEL()
    print("**************Delaunay Triangulation**********************\n")
    Edgecoords=Delau.printall()
    file = open("voronoi.txt", "a")
    file.write("**********Voronoi Diagram***************\n")
    file.close()
    D.ToTxt("Voronoi")
    file = open("voronoi.txt", "a")
    file.write("\n**********Delaunay Triangulation***************\n")
    file.close()
    Delau.ToTxt("Delaunay")


    print("Detected Voronoi Verteces: ")
    for i in range(numver):
        print(allVorVer[i,:])
    plt.show()
    print("Delau.PrintEdges()")
    Delau.PrintEdges()


