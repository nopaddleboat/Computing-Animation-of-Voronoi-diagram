import numpy as np
import copy
'''''quick sort: can be used to sort any list of objects by redefining the cmpsites functions'''
def cmpsites(twosites,indexes):#compare composite number
    flag1=(twosites[0,1]<twosites[1,1])
    flag2=(twosites[0,1]==twosites[1,1])and( twosites[0,0]>twosites[1,0] )
    flag=flag1 or flag2
    if flag:
        tmp=np.copy(twosites[0,:])
        twosites[0,:]=twosites[1,:]
        twosites[1,:]=tmp  
        tmp=np.copy(indexes)
        indexes=[tmp[1],tmp[0]]
    return twosites,indexes #deep copy,will not affect inputs

def quicksort(sites,indexes):
    len_sites= len(sites)
    if len_sites==1:
        return sites,indexes
    elif len_sites==2:#2 elements
        return cmpsites(sites,indexes)
    else:#more than 3 elements
        mid=int(len_sites/2)
        #divide and conquer
        sites[0:mid],indexes[0:mid]\
                                     =quicksort(sites[0:mid,:],indexes[0:mid])
        sites[mid:len_sites,:],indexes[mid:len_sites]\
                                    =quicksort(sites[mid:len_sites,:],indexes[mid:len_sites])
        newsites=np.zeros( sites.shape )
        newindexes=np.arange(0,len_sites)
        tail=[0,mid]
        mov_ind=0
        for i in range(0,len_sites):
            if tail[0]==mid:
                mov_ind=tail[1]
                tail[1]=tail[1]+1
            elif tail[1]==len_sites:
                mov_ind=tail[0]
                tail[0]=tail[0]+1
            else:
                twosites=sites[tail,:]
                tmp1,tmp2=cmpsites(twosites,tail)
                mov_ind=tmp2[0]
                if tmp2[0]>=mid:
                    tail[1]=tail[1]+1
                else:
                    tail[0]=tail[0]+1
            newsites[i]=sites[mov_ind]
            newindexes[i]=indexes[mov_ind]
        return newsites,newindexes

def comtwonodes(twonode,type):#true: 1 less than 2
    if type==1 or type==2:
        if twonode[0].interpnt[0]>twonode[1].interpnt[0]:
            return False
        return  True
    elif type==3 or type==4:
        if twonode[0].interpnt[1]>twonode[1].interpnt[1]:
            return  False
        return  True
    return  None

def SortTwoBnodes(twonode,type):#compare composite number1T 2B 3L 4R, output small to large
    newtwonode=copy.copy(twonode)
    if type==1 or type==2:
        if twonode[0].interpnt[0]>twonode[1].interpnt[0]:
            newtwonode[0]=twonode[1]
            newtwonode[1] = twonode[0]
    elif type==3 or type==4:
        if twonode[0].interpnt[1]>twonode[1].interpnt[1]:
            newtwonode[0]=twonode[1]
            newtwonode[1] = twonode[0]
    else:
        pass
    return newtwonode
def quicksortBList(iBList,type):
    # print("Initial")
    # for i in range(len(iBList)):
    #     print(iBList[i].interpnt)

    BList=copy.copy(iBList)
    len_sites= len(BList)
    print(len_sites)
    if len_sites==0:
        return iBList
    if len_sites==1:
        return copy.copy(BList)
    elif len_sites==2:#2 elements
        return copy.copy(SortTwoBnodes(BList,type))
    else:#more than 3 elements
        mid=int(len_sites/2)
        #divide and conquer
        BList[0:mid]=copy.copy(quicksortBList(BList[0:mid],type))
        BList[mid:len_sites]=copy.copy(quicksortBList(BList[mid:len_sites],type))
        newList=copy.copy(BList)
        pointer1=0
        pointer2=mid
        twosites = copy.copy(BList[0:2])
        # print("b4Final")
        # for i in range(len(newList)):
            # print(newList[i].interpnt)
        # print("--------")

        for i in range(0, len_sites):
            # print("p1 p2",pointer1,pointer2)
            if pointer1==mid:
                newList[i]=copy.copy(BList[pointer2])
                pointer2 = pointer2 + 1
                # print("line1")
            elif pointer2==len_sites:
                newList[i] = copy.copy(BList[pointer1])
                pointer1 = pointer1 + 1
                # print("line2")
            else:
                twosites[0] = copy.copy(BList[pointer1])
                twosites[1] = copy.copy(BList[pointer2])
                if comtwonodes(twosites, type) is True:
                    newList[i] = copy.copy(BList[pointer1])
                    pointer1=pointer1+1
                    # print("line3")
                else:
                    newList[i] = copy.copy(BList[pointer2])
                    pointer2 = pointer2 + 1
                    # print("line4")
        # print("Final")
        # for i in range(len(newList)):
        #     print(newList[i].interpnt)
        # print("--------")
        return newList
