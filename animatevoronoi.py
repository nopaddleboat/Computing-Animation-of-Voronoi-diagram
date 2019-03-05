import numpy as np
import re
import quicksort_swx as qs
import queue_swx
import status_swx
import DCEL
import handle_event
import matplotlib.pyplot as plt
from matplotlib import collections as mc
# def combreakpoint(sitei,sitej,ys):#Ys is y coordiante of sweep line
def FindAllLeaves(root):#recursively
    global leaves
    global leaveIndex
    p=root
    if p.lc is not None:
        if p.lc.lc is None:#leaf found!
            leaves[leaveIndex,:]=np.copy(p.lc.site[0:2])
            leaveIndex=leaveIndex+1
        else:#lc is innode
            FindAllLeaves(p.lc)
    if p.rc is not None:
        if p.rc.lc is None:#leaf found!
            leaves[leaveIndex,:]=np.copy(p.rc.site[0:2])
            leaveIndex=leaveIndex+1
        else:#rc is innode
            FindAllLeaves(p.rc)
    if p.rc is None and p.lc is None:#root is a leaf
        leaves[leaveIndex,:]=np.copy(p.site[0:2])
        leaveIndex=leaveIndex+1
    # print("leaveIndex",leaveIndex)
    return None

def plotvor(T,SpL):#plot the computed voronoi vertices and beachlines,SpL:sweepline

    global traceoutlines
    global IndTrace
    global trLineInd

    global leaves
    global leaveIndex

    global sites
    global Updateleaves

    global EndPoint

    stind=1#startindex
    plt.plot(sites[:, 0], sites[:, 1], 'o', color='C1')

    if T.root is not None and Updateleaves is True:#tree is not eampty
        leaveIndex = 0
        leaves=leaves*0
        FindAllLeaves(T.root)
        IndTrace = trLineInd  # start index for inserting tracing out line for traceoutlines
        trLineInd=trLineInd+leaveIndex-1#add traceoutlines lines to be traced out for each breakpoint
        stind=0
        EndPoint=0
    #plot all the beach lines and traced out break points

    for i in range(leaveIndex):#starting from zero leaveIndex
        if i==0:#leftmost leaf
            if leaveIndex>1:
                traceoutlines[IndTrace+i,stind,:]=status_swx.combreakpoint(leaves[i,:],leaves[i+1,:],SpL)
                if leaves[i, 1] - SpL == 0:
                    x1 = np.ones([1,50])*leaves[i,0]
                    y1 = np.linspace(SpL, traceoutlines[IndTrace+i,stind,1], num=50)
                else:
                    x1 = np.linspace(traceoutlines[IndTrace+i,stind,0] - 5, traceoutlines[IndTrace+i,stind,0], num=50)
                    y1 = ((x1 - leaves[i, 0]) ** 2 + leaves[i, 1] ** 2 - SpL ** 2) / 2.0 / (leaves[i, 1] - SpL)
                # print('line1', traceoutlines[IndTrace+i,stind,:])
            else:
                x1 = np.linspace(leaves[i,0]-5, leaves[i,0]+5, num=50)
                y1 = ((x1 - leaves[i, 0]) ** 2 + leaves[i, 1] ** 2 - SpL ** 2) / 2.0 / (leaves[i, 1] - SpL)
            #plotbeachline
            plt.plot(x1,y1,'green')


        elif i==leaveIndex-1:#rightmost leaf
            x2 = np.linspace(traceoutlines[IndTrace+i-1,stind,0], traceoutlines[IndTrace+i-1,stind,0] + 10, num=50)
            if leaves[i, 1] == SpL :
                y2 = np.linspace(SpL, traceoutlines[IndTrace+i,stind,1], num=50)
            else:
                y2 = ((x2 - leaves[i, 0]) ** 2 + leaves[i, 1] ** 2 - SpL ** 2) / 2.0 / (leaves[i, 1] - SpL)
            plt.plot(x2, y2,'green')
        else:#middle leaves
            traceoutlines[IndTrace + i,stind,:] = status_swx.combreakpoint(leaves[i,:], leaves[i + 1,:], SpL)
            x3 = np.linspace(traceoutlines[IndTrace + i-1,stind,0], traceoutlines[IndTrace + i,stind,0], num=50)

            if leaves[i, 1] == SpL:
                y3 = np.linspace(SpL, traceoutlines[IndTrace + i,stind,1], num=50)
            else:
                y3 = ((x3 - leaves[i, 0]) ** 2 + leaves[i, 1] ** 2 - SpL ** 2) / 2.0 / (leaves[i, 1] - SpL)
            plt.plot(x3, y3,'green')
    if EndPoint==0:
        traceoutlines[IndTrace:IndTrace+leaveIndex-1,1,:]=traceoutlines[IndTrace:IndTrace+leaveIndex-1,0,:]
            # print('line3', traceoutlines[IndTrace + i, stind, :])
    #plot breakpoints
    # for i in range(trLineInd):
    #     plt.plot(traceoutlines[IndTrace + i,0,0],traceoutlines[IndTrace + i,0,1],'+')
    #     plt.plot(traceoutlines[IndTrace + i, 1, 0], traceoutlines[IndTrace + i, 1, 1], '+')
    # print('traceoutlines[0:trLineInd]', traceoutlines[0:trLineInd])
    EndPoint=EndPoint+1
    return None

CRED = '\033[91m'
CEND = '\033[0m'
CYEL = '\033[93m'#yellow

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
n=1
#create len(sites) faces to DCEL and store it to DCEL
Facelist=[]
for i in range(0,num_sites):
    Fnode=D.addFace(i+1)
    Facelist.append(Fnode)
D.Flist=Facelist
Fzero=D.addFace(0)
# print(D.Flist)

allVorVer=np.zeros([100,2])
numver=0
plt.ion()#Turn interactive mode on.

traceoutlines=np.zeros([1000,2,2])
#animationi view box
bdsize=10.0
Leftb=np.amin(sites[:, 0])-bdsize
Rightb= np.amax(sites[:, 0])+bdsize
Bottomb= np.amin(sites[:, 1])-bdsize
Topb=np.amax(sites[:, 1])+bdsize
sweeprange=np.linspace(Topb-bdsize-1, Bottomb-20, num=200)
SweepPos=0#position of sweep line
event_pre=None
leaves=np.zeros([len(sites)*10,2])#store all the leaves
leaveIndex=0
trLineInd=0#trace out line index
IndTrace=0
flag=True
testnum=400
pausetime=0.001
Updateleaves=True
EndPoint=0#whether endpoint of the line has been inserted
fig, ax = plt.subplots()
while Q.outq is not None and flag is True:
    leaves=leaves*0
    event=Q.outq
    event_pre=event.Coord
    print("event_pre",event_pre)
    if n!=1:
        Q.delete(Q.outq)#delete processed data


    if event.leaf is None:#event is a site event
        print(CYEL,n, "th iteration: ", event.Coord,CEND)
        handle_event.HandleSiteEvent(D,Q,T,event)
    else:#event is a circle event
        allVorVer[numver,:]=event.cnt
        numver=numver+1
        print(CRED, n, "th iteration: ", event.Coord, "cnt:", event.cnt, CEND)
        handle_event.HandleCircleEvent(D,Q,T,event)
    n=n+1

    if Q.outq is not None:
        if event_pre[1]!=Q.outq.Coord[1]:
            Updateleaves = True
            while True:
                if SweepPos==len(sweeprange)-1:
                    break
                plt.plot([Leftb,Rightb],[sweeprange[SweepPos],sweeprange[SweepPos]])
                plotvor(T,sweeprange[SweepPos])
                if trLineInd>1 :
                   lines = np.copy(traceoutlines[0:trLineInd])
                   len1=lines.shape[0]
                   tmp = np.copy(lines[0:len1, 0, 1])
                   lines[0:len1, 0, 1] = np.copy(lines[0:len1, 1, 0])
                   lines[0:len1, 1, 0] = tmp
                   a1 = np.copy(np.reshape(lines, [1, len1*2, 2]))
                   plt.plot(*a1[0],color='black')


                plt.axis([Leftb-20, Rightb+20, Bottomb-20, Topb+20])
                plt.draw()
                plt.show()
                plt.pause(pausetime)
                Updateleaves = False
                if SweepPos ==testnum :#len(sweeprange)
                    flag=False
                    break
                if SweepPos==len(sweeprange)-1:
                    pass
                else:
                    plt.clf()
                if sweeprange[SweepPos+1]<Q.outq.Coord[1]:
                    if sweeprange[SweepPos]==Q.outq.Coord[1]:
                        break
                    else:
                        sweeprange[SweepPos]=Q.outq.Coord[1]
                else:
                    SweepPos = SweepPos + 1

    else:
        Updateleaves = True
        while SweepPos<len(sweeprange):
            plotvor(T, sweeprange[SweepPos])
            if trLineInd > 1:
                lines = np.copy(traceoutlines[0:trLineInd])
                len1 = lines.shape[0]
                tmp = np.copy(lines[0:len1, 0, 1])
                lines[0:len1, 0, 1] = np.copy(lines[0:len1, 1, 0])
                lines[0:len1, 1, 0] = tmp
                a1 = np.copy(np.reshape(lines, [1, len1 * 2, 2]))
                plt.plot(*a1[0], color='black')

            plt.axis([Leftb-20, Rightb+20, Bottomb-20, Topb+20])
            plt.draw()
            plt.pause(pausetime)
            Updateleaves = False
            if SweepPos ==len(sweeprange) :#len(sweeprange)
                flag = False
                break
            if SweepPos == len(sweeprange)-1:
                pass
            else:
                plt.clf()
            SweepPos = SweepPos + 1

print("Detected Voronoi Verteces: ")
for i in range(numver):
    print(allVorVer[i,:])
# box=handle_event.AddBoundingBox(T,D,allVorVer[0:numver,:],Fzero)
#plot the voronoi diagram from DCEL:
Edgecoords=D.printall()
Edgecoords = np.asarray(Edgecoords)
plt.axis([box[0,0]-2, box[1,0]+2, box[2,1]-2, box[0,1]+2])#[Leftb,Topb],[Rightb,Topb],[Rightb,Bottomb],[Leftb,Bottomb]
for i in range(len(Edgecoords)):
    plt.plot(Edgecoords[i][0,[0,2]], Edgecoords[i][0,[1,3]])
    plt.plot(Edgecoords[i][0,0], Edgecoords[i][0,1], marker='o', markersize=10)
for i in range(len(sites)):
    plt.plot(sites[i, 0], sites[i, 1], marker='+', markersize=10,color = 'red')
plt.show()
