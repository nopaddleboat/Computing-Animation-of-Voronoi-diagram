import numpy as np
from numpy import linalg as LA
import quicksort_swx as qs


class VertexNode:
    def __init__(self, index):
        '''data'''
        self.Vertex = index  # vertex index
        self.Coord = np.array([0.0, 0.0, 0.0])  # coordinate, third will not be used
        self.IncidentEdge = None  # IncidentEdge: Edgenode
        '''pointers'''
        self.next = None
        self.prev = None

    def setdata(self, Coord=None, IncidentEdge=None):
        if Coord is not None:  # do not use Coord!=None, 'cause the result will be [true,true]
            self.Coord = Coord
        if IncidentEdge is not None:
            self.IncidentEdge = IncidentEdge


class FaceNode:
    def __init__(self, index):
        '''data'''
        self.Face = index
        self.site = np.array([0.0, 0.0, 0.0])  # store site coordinate of the voronoi cell
        self.OuterComponent = None  # EdgeNode
        self.InnerComponent = None  # edge node
        '''pointers'''
        self.next = None
        self.prev = None

    def setdata(self, OuterComponent=None, InnerComponent=None):
        if OuterComponent is not None:
            self.OuterComponent = OuterComponent
        if InnerComponent is not None:
            self.InnerComponent = InnerComponent


class EdgeNode:  # i.e., half edge node
    def __init__(self, data):
        '''data'''
        self.EdgeName = data  # is integer array, e.g., 1,2 means e1,2
        self.Origin = None  # vertex node
        self.Twin = None  # same as edge name
        self.IncidentFace = None  # FaceNode
        self.Next = None  # integer array
        self.Prev = None  # integer array
        '''pointers'''
        self.next = None
        self.prev = None

    def setdata(self, Origin=None, Twin=None, IncidentFace=None, Next=None, Prev=None):
        if Origin is not None:
            self.Origin = Origin
        if Twin is not None:
            self.Twin = Twin
        if IncidentFace is not None:
            self.IncidentFace = IncidentFace
        if Next is not None:
            self.Next = Next
        if Prev is not None:
            self.Prev = Prev


class DCEL:
    def __init__(self):
        self.VertexHead = None
        self.FaceHead = None
        self.EdgeHead = None  # contains three tables: vertext, face, edge, stored in three linked list
        self.vertexnum = 0
        self.edgenum = 0  # num of edge pair,i.e., total number of edges/2
        self.Flist = None

    def addVertex(self):
        self.vertexnum = self.vertexnum + 1
        Vnode = VertexNode(self.vertexnum)
        if self.VertexHead == None:
            self.VertexHead = Vnode
        else:
            Vnode.next = self.VertexHead
            Vnode.next.prev = Vnode
            self.VertexHead = Vnode
        return Vnode

    def addFace(self, Face):
        Fnode = FaceNode(Face)
        if self.FaceHead is None:
            self.FaceHead = Fnode
        else:
            Fnode.next = self.FaceHead
            Fnode.next.prev = Fnode
            self.FaceHead = Fnode
        print("Face added! Face,self.FaceHead.Face",Face,self.FaceHead.Face)
        return Fnode

    def addEdge(self, EdgeName):  # Edgename=1 original edge,2 twin edge
        if EdgeName == 1:
            self.edgenum = self.edgenum + 1
        Enode = EdgeNode(np.array([self.edgenum, EdgeName]))
        if self.EdgeHead is None:
            self.EdgeHead = Enode
        else:
            Enode.next = self.EdgeHead
            Enode.next.prev = Enode
            self.EdgeHead = Enode
        print("Edge", "type: ", EdgeName, self.EdgeHead.EdgeName, "added!!!!!!!!!!!!!!!!!!!!!!!!!!!!!")
        return Enode

    def deleteNode(self, Ndel):  # Edel to be deleted
        if Ndel is self.VertexHead:
            self.VertexHead = Ndel.prev
            Ndel.prev.Next = None
        else:
            Ndel.prev.next = Ndel.next
            Ndel.next.prev = Ndel.prev
        # print(Ndel.Coord,"Ndel deleted!")

    def deleteEdge(self, Edel):  # Edel to be deleted
        if Edel is self.EdgeHead:
            self.EdgeHead = Edel.prev
            Edel.prev.Next = None
        else:
            Edel.prev.next = Edel.next
            Edel.next.prev = Edel.prev
        # print(Edel.EdgeName, "Edge deleted!")

    def printall(self):

        '''print vertex table'''
        p = self.VertexHead
        print("Double Connected Edge List")
        print("Vertex    Coordinates   IncidentEdge")
        while True:
            if p.IncidentEdge is None:
                print("p.IncidentEdge is None",p.Coord)
            print("v%d" % p.Vertex, "  (%f,%f)" % (p.Coord[0], p.Coord[1]), " e%d,%d" % \
                  (p.IncidentEdge.EdgeName[0], p.IncidentEdge.EdgeName[1]))
            if p.next == None:
                break
            p = p.next
        #
        p = self.FaceHead
        print("Face  OuterComponent   InnerComponent")
        while True:
            if p.OuterComponent == None:
                Outer = "nil"
            else:
                Outer = "e" + str(p.OuterComponent.EdgeName[0]) + "," + str(p.OuterComponent.EdgeName[1])
            if p.InnerComponent == None:
                Inner = "nil"
            else:
                Inner = "e" + str(p.InnerComponent.EdgeName[0]) + "," + str(p.InnerComponent.EdgeName[1])
            print("f%d" % p.Face, Outer, " ", Inner)
            if p.next == None:
                break
            p = p.next

        p = self.EdgeHead
        Edgecoords = []
        Edges = []
        print("HalfEdge  Origin   Twin IncidentFace Next Prev")
        while True:
            Edges.append(p)
            # print("p.EdgeName p.Twin.Origin", p.EdgeName,p.Twin.Origin.Coord)
            # print("p.EdgeName p.Origin,p.Twin.Origin",p.EdgeName,p.Origin.Coord,p.Twin.Origin.Coord)
            outE = "e" + str(p.EdgeName[0]) + "," + str(p.EdgeName[1])
            outO = "v" + str(p.Origin.Vertex)
            outT = "e" + str(p.Twin.EdgeName[0]) + "," + str(p.Twin.EdgeName[1])
            outI = "f" + str(p.IncidentFace.Face)
            if p.Next is None:
                print(" p.Next is None", p.EdgeName,p.Origin.Coord,p.Twin.Origin.Coord)
            outN = "e" + str(p.Next.EdgeName[0]) + "," + str(p.Next.EdgeName[1])
            outP = "e" + str(p.Prev.EdgeName[0]) + "," + str(p.Prev.EdgeName[1])
            print(outE, outO, outT, outI, outN, outP)
            if p.next == None:
                break
            p = p.next
            tmp = np.zeros([1, 4])
            tmp[0, 0:2] = p.Origin.Coord[0:2]
            tmp[0, 2:4] = p.Twin.Origin.Coord[0:2]
            if p.EdgeName[1] == 1:
                Edgecoords.append(np.copy(tmp))
        return Edgecoords, Edges

    def ToTxt(self,graphtype):#output the result to txt file as required
        file = open("voronoi.txt", "a")
        VTxt=[]
        ETxt=[]
        FTxt=[]
        Ename='e'
        Fname='f'
        if graphtype=="Voronoi":#Voronoi
            Ename='e'
            Fname='c'
        else:#Delaunay
            Ename='d'
            Fname = 't'
        '''print vertex table'''
        p = self.VertexHead
        # print("Double Connected Edge List")
        # print("Vertex    Coordinates   IncidentEdge")
        while True:
            string1="v%d" % p.Vertex+ "  (%f,%f)" % (p.Coord[0], p.Coord[1])+" "+Ename+"%d,%d" % \
                  (p.IncidentEdge.Origin.Vertex, p.IncidentEdge.Twin.Origin.Vertex)
            VTxt.append(string1)
            if p.next == None:
                break
            p = p.next

        p = self.FaceHead
        # print("Face  OuterComponent   InnerComponent")
        while True:
            if p.OuterComponent == None:
                Outer = "nil"
            else:
                Outer = Ename + str(p.OuterComponent.Origin.Vertex) + "," + str(p.OuterComponent.Twin.Origin.Vertex)
            if p.InnerComponent == None:
                Inner = "nil "
            else:
                Inner = Ename + str(p.InnerComponent.Origin.Vertex) + "," + str(p.InnerComponent.Twin.Origin.Vertex)
            string1=Fname+"%d" % p.Face+" "+Outer+" "+Inner
            FTxt.append(string1)
            if p.next == None:
                break
            p = p.next

        p = self.EdgeHead
        Edgecoords = []
        Edges = []
        # print("HalfEdge  Origin   Twin IncidentFace Next Prev")
        while True:
            Edges.append(p)
            # print("p.EdgeName p.Twin.Origin", p.EdgeName,p.Twin.Origin.Coord)
            # print("p.EdgeName p.Origin,p.Twin.Origin",p.EdgeName,p.Origin.Coord,p.Twin.Origin.Coord)
            outE = Ename + str(p.Origin.Vertex) + "," + str(p.Twin.Origin.Vertex)
            outO = "v" + str(p.Origin.Vertex)
            outT = Ename + str(p.Twin.Origin.Vertex) + "," + str(p.Origin.Vertex)
            outI = Fname + str(p.IncidentFace.Face)
            outN = Ename + str(p.Next.Origin.Vertex) + "," + str(p.Next.Twin.Origin.Vertex)
            outP = Ename + str(p.Prev.Origin.Vertex) + "," + str(p.Prev.Twin.Origin.Vertex)
            string1=outE+" "+outO+" "+outT+" "+outI+" "+outN+" "+outP
            ETxt.append(string1)
            if p.next == None:
                break
            p = p.next
            tmp = np.zeros([1, 4])
            tmp[0, 0:2] = p.Origin.Coord[0:2]
            tmp[0, 2:4] = p.Twin.Origin.Coord[0:2]
            if p.EdgeName[1] == 1:
                Edgecoords.append(np.copy(tmp))
        #write to txt file
        #write vertex
        for i in range(0,len(VTxt)):
            file.write(VTxt[len(VTxt)-1-i]+"\n")
        # write face
        file.write("\n\n"+FTxt[0]+"\n")
        for i in range(0, len(FTxt) - 1):
            file.write(FTxt[len(FTxt) - 1 - i] + "\n")
        #write edge
        file.write("\n\n")
        for i in range(0, len(ETxt) ):
            file.write(ETxt[len(ETxt) - 1 - i] + "\n")
        file.close()
        return None

    def PrintEdges(self):
        p = self.EdgeHead
        print("HalfEdge  Origin   Twin IncidentFace Next Prev")
        out = None
        while True:
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
                if p.IncidentFace is None or p.Twin.IncidentFace is None:
                    out=out+"None"
                else:
                    out=out+"complete face"
                print(out)#, p.IncidentFace.Face, p.Twin.IncidentFace.Face, p.EdgeName
            if p.next == None:
                break
            p = p.next

    def completeDCEL(self):
        # complete the information of face and vertex through traverse all the edges
        p = self.EdgeHead
        n = 0
        while True:
            n=n+1
            if p is None:
                break
            p.Origin.IncidentEdge = p
            if p.IncidentFace.Face == 0:
                p.IncidentFace.InnerComponent = p
            else:
                p.IncidentFace.OuterComponent = p
            p = p.next

        return None

    # deal with degeneracies where co-circular points generate zero length edge
    # step1 merge two vertices of the edge
    # step2 delete one of the vertices
    def mergeVertices(self):  # compute the length of each edge
        self.PrintEdges()
        p = self.EdgeHead
        n = 0
        while True:
            # print(n,"n")
            if p is None:
                break
            if p.EdgeName[1] == 1:
                # print(LA.norm(p.Origin.Coord[0:2]-p.Twin.Origin.Coord[0:2]))
                if LA.norm(p.Origin.Coord[0:2] - p.Twin.Origin.Coord[0:2]) < 1e-3:
                    # print("detected!", p.Origin.Coord, p.Twin.Origin.Coord)
                    n = n + 1
                    mergenode = p.Origin  # vertex to be preserved
                    p.Next.Prev = p.Prev
                    p.Prev.Next = p.Next
                    p.Next.Origin = mergenode
                    p1 = p.Next.Twin
                    while True:
                        if p1.Next is p.Twin:
                            break
                        p1.Next.Origin = mergenode
                        p1 = p1.Next.Twin
                        # print("node changed!")
                    p1.Twin.Origin = mergenode
                    p.Twin.Next.Prev = p.Twin.Prev
                    p.Twin.Prev.Next = p.Twin.Next
                    # print("node changed!END")

                    self.deleteNode(p.Twin.Origin)
                    self.deleteEdge(p)
                    self.deleteEdge(p.Twin)

            p = p.next  # edgenode deleted


'''Recursivley construct Delaunay triangle'''
def AddDTriangle(AdjEdge1, AdjDir1, box, Delau, sortedsites, VerList,FList,DeVertexList):
    coord1 = AdjDir1.Twin.Origin.Coord
    p = AdjDir1
    global FN

    cond1 = coord1[0] != box[0, 0] and coord1[0] != box[1, 0] and coord1[1] != box[0, 1] and coord1[1] != box[2, 1]
    # print("AdjDir1 cond2",AdjDir1.Origin.Coord,AdjDir1.Twin.Origin.Coord, AdjEdge1.EdgeName, AdjEdge1.Origin.Coord, AdjEdge1.Twin.Origin.Coord)
    cond2 = VerList[AdjDir1.Twin.Origin.Vertex] is False
    if cond1 and cond2:  # inner vertex found
        V1 = AdjEdge1.Twin.Origin
        V2 = AdjEdge1.Origin
        V3 = DeVertexList[p.Next.Twin.IncidentFace.Face - 1]
        E01 = AdjEdge1.Twin
        E02 = AdjEdge1
        E11 = Delau.addEdge(1)
        E12 = Delau.addEdge(2)
        E21 = Delau.addEdge(1)
        E22 = Delau.addEdge(2)
        F1 = Delau.addFace(FN)
        FN = FN + 1
        FList.append(F1)
        AdjEdge1.Twin.setdata(Origin=None, Twin=None, IncidentFace=F1, Next=E11, Prev=E21)

        E11.setdata(Origin=V2, Twin=E12, IncidentFace=F1, Next=E21, Prev=E01)
        E12.setdata(Origin=V3, Twin=E11, IncidentFace=None, Next=None, Prev=None)
        E21.setdata(Origin=V3, Twin=E22, IncidentFace=F1, Next=E01, Prev=E11)
        E22.setdata(Origin=V1, Twin=E21, IncidentFace=None, Next=None, Prev=None)
        # print("sub E11 E21", V2.Coord, V3.Coord, V3.Coord, V1.Coord)

        AdjEdge = []
        AdjDirList = []


        AdjEdge.append(E21)
        Epre = E11
        Vpre = [V2, V3]

        AdjDirList.append(p.Next)
        p1 = p
        n=1
        while n<10:
            n=n+1
            # print("n",n)
            # print("p1.Next.Twin.Next is p.Twin.Prev.Twin:",p1.Next.Twin.Next.EdgeName,p.Twin.Prev.Twin.EdgeName)
            if p1.Next.Twin.Next is p.Twin.Prev.Twin:  # degree not more than 3
                break
            # V4 = Delau.addVertex()
            E31 = Delau.addEdge(1)
            E32 = Delau.addEdge(2)
            E41 = Delau.addEdge(1)
            E42 = Delau.addEdge(2)
            # V4.Coord = sortedsites[p1.Next.Twin.Next.Twin.IncidentFace.Face - 1]
            V4 = DeVertexList[p1.Next.Twin.Next.Twin.IncidentFace.Face - 1]

            AdjEdge.append(E41)
            F1 = Delau.addFace(FN)
            FList.append(F1)
            FN = FN + 1
            # print("len(FList),FN",len(FList),FN)
            Epre.Twin.setdata(Origin=None, Twin=None, IncidentFace=F1, Next=E31, Prev=E41)
            E31.setdata(Origin=Vpre[0], Twin=E32, IncidentFace=F1, Next=E41, Prev=Epre.Twin)
            E32.setdata(Origin=V4, Twin=E31, IncidentFace=None, Next=None, Prev=None)
            E41.setdata(Origin=V4, Twin=E42, IncidentFace=F1, Next=Epre.Twin, Prev=E31)
            E42.setdata(Origin=Vpre[1], Twin=E41, IncidentFace=None, Next=None, Prev=None)
            # print("sub E31 E41", Vpre[0].Coord, V4.Coord, V4.Coord, Vpre[1].Coord)
            Vpre[1] = V4
            Epre = E31
            AdjDirList.append(p1.Next.Twin.Next)
            p1 = p1.Next.Twin

        AdjEdge.append(Epre)
        AdjDirList.append(p.Twin.Prev.Twin)

        VerList[p.Twin.Origin.Vertex] = True  # mark the vertex dealt with
        for i in range(len(AdjEdge)):
            if AdjDirList[i].Twin.Origin:
                AddDTriangle(AdjEdge[i], AdjDirList[i], box, Delau, sortedsites, VerList,FList,DeVertexList)
    return None


'''Generete DCEL from Voronoi DCEL'''
'''
boxcoords=np.array([[Leftb,Topb],[Rightb,Topb],[Rightb,Bottomb],[Leftb,Bottomb]])
'''


def VorToDel(D, box, sortedsites):  # D represent DCEL of voronoi
    print("constructing Delaunay triangulation from Voronoi diagram")
    global FN
    FN = 1
    Delau = DCEL()  # DCEL for Delaunay triangulation
    Nver = D.VertexHead.Vertex + 1  # number of vertex,some may be deleted e.g. 1 2 4 5 ...
    VerList = []  # for missing replace with False
    DeVertexList=[]
    for i in range(Nver):
        VerList.append(False)  # once AddDTriangle act on one node, set the value to True
    for i in range(len(sortedsites)):
        newV=Delau.addVertex()
        newV.Coord = sortedsites[i]
        DeVertexList.append(newV)
    p = D.EdgeHead
    while True:
        if p.IncidentFace.Face != 0 and p.Twin.IncidentFace.Face != 0:  # inner edge found
            break
        p = p.next
    # find a inner vertex
    coord0 = p.Origin.Coord
    if coord0[0] != box[0, 0] and coord0[0] != box[1, 0] and coord0[1] != box[0, 1] and coord0[1] != box[
        2, 1]:  # inner vertex found
        p = p.Twin
    # deal with the first point to construct delaunay triangulation
    V1=DeVertexList[p.IncidentFace.Face - 1]
    V2 = DeVertexList[p.Twin.IncidentFace.Face - 1]
    V3 = DeVertexList[p.Next.Twin.IncidentFace.Face - 1]

    # print("sortedsites",sortedsites,p.Origin.Coord,p.Twin.Origin.Coord,p.IncidentFace.Face,p.Twin.IncidentFace.Face)
    # print("V1.Coord", V1.Coord, V2.Coord, V3.Coord)
    E01 = Delau.addEdge(1)  # intersect p
    E02 = Delau.addEdge(2)
    E11 = Delau.addEdge(1)
    E12 = Delau.addEdge(2)
    E21 = Delau.addEdge(1)
    E22 = Delau.addEdge(2)
    F1 = Delau.addFace(FN)
    E01.setdata(Origin=V1, Twin=E02, IncidentFace=F1, Next=E11, Prev=E21)
    E02.setdata(Origin=V2, Twin=E01, IncidentFace=None, Next=None, Prev=None)
    E11.setdata(Origin=V2, Twin=E12, IncidentFace=F1, Next=E21, Prev=E01)
    E12.setdata(Origin=V3, Twin=E11, IncidentFace=None, Next=None, Prev=None)
    E21.setdata(Origin=V3, Twin=E22, IncidentFace=F1, Next=E01, Prev=E11)
    E22.setdata(Origin=V1, Twin=E21, IncidentFace=None, Next=None, Prev=None)
    # print("main",V1.Coord,V2.Coord,V2.Coord,V3.Coord,V3.Coord,V1.Coord)

    AdjEdge = []
    AdjDirList = []
    FList = []
    FList.append(F1)
    FN = FN + 1
    AdjEdge.append(E01)
    AdjEdge.append(E21)
    Epre = E11
    Vpre = [V2, V3]
    AdjDirList.append(p.Twin)
    AdjDirList.append(p.Next)
    p1 = p
    n=1
    while n<10:
        n=n+1
        # print("Main:n",n)

        if p1.Next.Twin.Next is p.Twin.Prev.Twin:  # degree not more than 3
            break

        E31 = Delau.addEdge(1)
        E32 = Delau.addEdge(2)
        E41 = Delau.addEdge(1)
        E42 = Delau.addEdge(2)
        V4 = DeVertexList[p1.Next.Twin.Next.Twin.IncidentFace.Face - 1]
        AdjEdge.append(E41)
        F1 = Delau.addFace(FN)
        FList.append(F1)
        FN = FN + 1
        Epre.Twin.setdata(Origin=None, Twin=None, IncidentFace=F1, Next=E31, Prev=E41)
        E31.setdata(Origin=Vpre[0], Twin=E32, IncidentFace=F1, Next=E41, Prev=Epre.Twin)
        E32.setdata(Origin=V4, Twin=E31, IncidentFace=None, Next=None, Prev=None)
        E41.setdata(Origin=V4, Twin=E42, IncidentFace=F1, Next=Epre.Twin, Prev=E31)
        E42.setdata(Origin=Vpre[1], Twin=E41, IncidentFace=None, Next=None, Prev=None)
        # print("main E31 E41", Vpre[0].Coord, V4.Coord, V4.Coord, Vpre[1].Coord)
        Vpre[1] = V4
        Epre = E31
        AdjDirList.append(p1.Next.Twin.Next)
        p1 = p1.Next.Twin

    AdjEdge.append(Epre)
    AdjDirList.append(p.Twin.Prev.Twin)

    VerList[p.Twin.Origin.Vertex] = True  # mark the vertex dealt with
    for i in range(len(AdjEdge)):
        if AdjDirList[i].Twin.Origin:
            AddDTriangle(AdjEdge[i], AdjDirList[i], box, Delau, sortedsites, VerList,FList,DeVertexList)

    # print("Delau.PrintEdges()")
    # Delau.PrintEdges()
    # after triangulation completed, add face zero i.e., f0 to the edges
    # sort all the 1 edges and combine the duplicated edges
    AllEdges = []
    AllEdgeInd = []
    p = Delau.EdgeHead
    while True:
        if p is None:
            break
        if p.EdgeName[1] == 1:
            AllEdges.append(p)
            AllEdgeInd.append(np.array([p.Origin.Coord[2], p.Twin.Origin.Coord[2]]))
        p = p.next
    AllEdgeInd_np=np.zeros([len(AllEdgeInd),2])
    for i in range(len(AllEdgeInd)):
        if AllEdgeInd[i][0]<AllEdgeInd[i][1]:
            AllEdgeInd_np[i,:]=np.copy(AllEdgeInd[i])
        else:
            AllEdgeInd_np[i, 0] = AllEdgeInd[i][1]
            AllEdgeInd_np[i, 1] = AllEdgeInd[i][0]
        # print("AllEdgeInd", AllEdgeInd[i])

    newsites, newindexes=qs.quicksort(AllEdgeInd_np,np.arange(len(AllEdgeInd_np)))
    # print(newsites)
    # p11 = Delau.VertexHead
    # print("all vertices")
    # while True:
    #     if p11 is None:
    #         break
    #     print(p11.Vertex, p11.Coord)
    #     p11 = p11.next
    for i in range(len(AllEdgeInd)-1):
        # print("New AllEdgeInd", AllEdgeInd[newindexes[i]],newindexes[i],newsites[i])
        if newsites[i,0]==newsites[i+1,0] and newsites[i,1]==newsites[i+1,1]:
            Edel=AllEdges[newindexes[i+1]]#to be deleted
            AllEdges[newindexes[i]].Twin.setdata(Origin=None, Twin=None,IncidentFace=Edel.IncidentFace,
                                                 Next=Edel.Prev, Prev=Edel.Next)
            Edel.Next.Prev=AllEdges[newindexes[i]].Twin
            Edel.Prev.Next = AllEdges[newindexes[i]].Twin

            Delau.deleteEdge(Edel)
            Delau.deleteEdge(Edel.Twin)

    F0 = Delau.addFace(0)
    p = Delau.EdgeHead
    nloop2=0
    while nloop2<500:
        nloop2=nloop2+1
        if p is None:
            break
        if p.IncidentFace is None:
            p.IncidentFace=F0
            #also add Next and Prev information for such edges
            #Next
            p1=p.Twin.Prev
            # print("p.IncidentFace is None",p.Origin.Coord,p.Twin.Origin.Coord)
            while True:
                # print("while1")
                if p1.Twin.IncidentFace is None or p1.Twin.IncidentFace == F0:
                    p.Next=p1.Twin
                    break
                p1=p1.Twin.Prev
            #Prev
            p1 = p.Twin.Next
            while True:
                # print("while2")
                if p1.Twin.IncidentFace is None or p1.Twin.IncidentFace == F0:
                    p.Prev = p1.Twin
                    break
                p1 = p1.Twin.Next
        p=p.next

    # print("Delau.PrintEdges()")
    # Delau.PrintEdges()
    print("End of delaunay!")
    return Delau






