import copy
import time
import matplotlib.pyplot as plt
from random import randint

class Neighber(object):
    def __init__(self, str):
        self.str=str
        self.idxs=[]

#@Query wmer
class Qwmer(object):
    def __init__(self, str, i):
        self.str=str
        self.idx=i

class MSP(object):
    def __init__(self, score, sQ, eQ, sDB, eDB):
        self.score= score
        self.startQuery= sQ
        self.endQuery=eQ
        self.startDB=sDB
        self.endDB=eDB

#PARAMETERS:
optimizaion_flag=1 #proning for extending neighbers.
d=1 #d=1 show time \ progress
ass2graphPlot=5 #values: 1-> plot avarege neighbers amount. 2-> plot running time (PRT2ITS=iterations).
                # 3-> 1+2.
                # 4->comparing with local alignment log
                # 0->default: part 1 with parameters w_val,T_val,X_val
                # 5-> report top N MSPs. (part 3)
N=5 #num of MSPs to report
runtimes_mode=0;#mode=0-> X random. 1-> T random. 2-> w random. 5-> all random
PART2ITRS=10 #part 2 iterations for runtimes plot
#parameters values (when 'ass2graphPlot' != 2, 3)
w_val=4
T_val=31
X_val=10

def main():
    runtime=-1
    db_file = open('db2.fasta', 'r')
    query_file = open('query2.fasta', 'r')
    name1, query=readFasta(query_file)
    name2, db=readFasta(db_file)

    if (ass2graphPlot==1 or ass2graphPlot==3): #For w=3 and different values of T,
        # plot a graph of the average number of neighbors for each w-mer.
        neighbers_avrg_plot(query, 3)

    print("Part 1:")
    iterations=1
    if (ass2graphPlot==2 or ass2graphPlot==3): #ploting runtimes by parameters (part 2)
        iterations=PART2ITRS
        Xvalues=[]; Tvalues=[]; Wvalues=[]; runtimes=[]; names=[]

    for i in range(iterations):
        t=time.time()

        if (ass2graphPlot == 2 or ass2graphPlot == 3):
            print("Iteration ",(i+1))
            w, T, X= rand_params(runtimes_mode)
        else:
            w=w_val; T=T_val; X=X_val

        print("Starting with w=",w,"T=",T,"X=",X)

        dbmap=preproc(db, w)

        if (d==1):
            print(time.time()-t, "(mapping neighbers)")
        neighbers_map=mapNeighbers(query, T, w)
        if (d==1):
            print(time.time()-t,"(updating idxs)")
        updateIdxs(neighbers_map, dbmap) #of db

    #   printQmersAndNeighbers(neighbers_map) #for debbuging
        if (d==1):
            print(time.time()-t,"(extending neighbers)")
        MSPs=extendNeighbers(neighbers_map, db, query, X, w)
        sortedMSPs=sortedReportMSPs(MSPs)

        runtime=time.time()-t
        print("Iteration ",(i+1)," runtime:",time.time()-t)
        if (ass2graphPlot == 2 or ass2graphPlot == 3):
            updateRuntimePlotParams(Xvalues,Tvalues,Wvalues,runtimes,names, X, T, w, runtime)

    if (ass2graphPlot==2 or ass2graphPlot==3):
        plot_by_parameteres(names, Xvalues, Tvalues, Wvalues, runtimes)

    if (ass2graphPlot==4): #comparing to local alignment (part 2)
        compareWithLAlog(sortedMSPs, runtime)

    if (ass2graphPlot==5):
        reportTopMSPs(sortedMSPs, N)

def updateRuntimePlotParams(Xvalues,Tvalues,Wvalues,runtimes,names, X, T, w, runtime):
    Xvalues.append(X)  # plot by Xs
    Tvalues.append(T)  # plot by X
    Wvalues.append(w)  # plot by X
    runtimes.append(runtime)
    names.append(" X=" + str(X) + ", T=" + str(T) + ", w=" + str(w) + " | ")

def rand_params(mode): #mode=0-> X random. 1-> T random. 2-> w random. 5-> all random
    if mode==0:
        w=w_val; T=T_val; X= randint(3, 8)
    if mode==1:
        w=w_val; T = randint(4, 7) * (w + 3); X = X_val
    if mode==2:
        w = randint(3, 8); T = T_val; X = X_val
    if mode==5:
        w = randint(3, 8)
        T = randint(4, 7) * (w + 2)
        X = randint(3, 8)
    return w,T,X

def plot_by_parameteres(names, Xvalues, Tvalues, Wvalues, runtimes):
    title=""
    for x in names:
        title=title+x
    plt.title(title+"\nRed by X. Blue by T. Green by W")
    plt.plot(Xvalues, runtimes, 'rs', Tvalues, runtimes, 'bs', Wvalues, runtimes, 'gs')
    plt.xlabel("X/T/w")
    plt.ylabel("Runtime")
    plt.show()

def neighbers_avrg_plot(query, w):
    TvaluesList=[]
    neighbersAmountsList=[]
    print("Calculating graph of Part 2")
    for i in range(8,20,1):
        if (d==1):
            print("calculatin for T=",i)
        T=i
        neighbers_map = mapNeighbers(query, T, w)
        neighAmount=0
        for wmer in neighbers_map:
            neighAmount += len(neighbers_map[wmer])
        avrg=neighAmount / len(neighbers_map.keys())
        TvaluesList.append(T)
        neighbersAmountsList.append(avrg)
    plt.plot(TvaluesList, neighbersAmountsList)
    plt.xlabel("T")
    plt.ylabel("Average neighbers Amount for each w-mer")
    plt.show()

def compareWithLAlog(sortedMSPs, runtime):
    print("----------------------------")
    print("Comparing to Local Alignment (query1, db1):::")
    print("Local Alignment Results: score = 6987\nquery: (16, 1382)\ndb: (2999, 4352)")
    print("Local Alignment runtime: 320 seconds")
    reportTopMSPs(sortedMSPs, 3)
    print("BLAST runtime: ",runtime, "seconds")
    print("----------------------------")

def reportTopMSPs(sortedList, x):
    top_msps=[]
    if (x > len(sortedList)):
        print ("There are ", len(sortedList),"( <",x,") MSPs")
        x = len(sortedList)
    for i in range(x):
        top_msps.append(sortedList[i])
    print("Top",x,"MSPs:::")
    for msp in top_msps:
        print("score: ", msp.score, ", Query: ", msp.startQuery, "-", msp.endQuery, ", DB:", msp.startDB, "-",
              msp.endDB)


def printQmersAndNeighbers(neighbers_map):
    for i in neighbers_map.keys():
        print("---")
        print ("Qmer idx:",i.idx,"str:",i.str)
        for j in range(len(neighbers_map[i])):
            print ("Neigh idxs:",neighbers_map[i][j].idxs,"str",neighbers_map[i][j].str)
            print ("score: ", score(i.str, neighbers_map[i][j].str))
        print("---")


def updateIdxs(neighs, dbmap):
    for wmer in neighs.keys():
        for neigh in neighs[wmer]:
            neighStr = ''.join(neigh.str)
            if neighStr in dbmap.keys():
          #      print(neigh.str)
              #  if (equals1(list(key), list(neigh.str)))==1:
                    neigh.idxs = dbmap[neighStr]
                    break; #key appears maximum once on db_map.key()

def mapNeighbers(query, T, w):
    window=""
    i=0
    Qlen=len(query)
    dict={}
    for i in range (w):
        window=window + query[i]
    i=0
    while i<Qlen-w:
        wmer=Qwmer(window, i)
        dict[wmer] = findNeighbers(window, T, w)
        window=window[1:]+query[i+w]
        i+=1
    return dict

def findNeighbers(window, T, w):
    amino_dict={'A':0, 'R':1, 'N':2, 'D':3, 'C':4, 'Q':5, 'E':6, 'G':7, 'H':8, 'I':9,
                'L':10, 'K':11, 'M':12, 'F':13, 'P':14, 'S':15, 'T':16, 'W':17, 'Y':18,
                'V':19, 'B':20, 'Z':21, 'X':22} #no space ('*':23)
    list1=[]
    temp=window
    firstNeigh=Neighber(list(window))
    if (score(window, window)>=T):
        list1.append(firstNeigh) #str is neigh of himself (the only one that doesn't checked on the recursion)
    efficient_find(amino_dict, window, list(temp), list1, 0, w, T)
    return list1

#to be checked again
#@returns list of neighbers of that windows
def efficient_find(amino_dict, wmer, temp, list1, i, w, T):
    if i==w:
        return;
    elif i<w:
        for x in amino_dict.keys():
            nextTemp = copy.deepcopy(temp)
            nextTemp[i] = x
            if (score(nextTemp, wmer)>=T and equals1(temp,nextTemp)==-1): #its a neighber (also thats optimum from this position)
                toAdd=copy.deepcopy(nextTemp)                             #second condition says its not the same string like last recursion
                neigh=Neighber(toAdd)
                list1.append(neigh)
                efficient_find(amino_dict, wmer, nextTemp, list1, i + 1, w, T)
            elif (score(nextTemp, wmer)>=T):
                efficient_find(amino_dict, wmer, nextTemp, list1, i + 1, w, T)
        return;

#checks content of l1 an l2.
def equals1(l1,l2):
    if len(l1)!=len (l2): return -1
    for i in range (len(l1)):
        if (l1[i]!=l2[i]): return -1
    return 1

def score(s1, s2): #same length
    amino_dict={'A':0, 'R':1, 'N':2, 'D':3, 'C':4, 'Q':5, 'E':6, 'G':7, 'H':8, 'I':9,
                'L':10, 'K':11, 'M':12, 'F':13, 'P':14, 'S':15, 'T':16, 'W':17, 'Y':18,
                'V':19, 'B':20, 'Z':21, 'X':22, '*':23}
    sum=0;
    for i in range(len(s1)):
        sum+=blosum[amino_dict[s1[i]]][amino_dict[s2[i]]]
    return sum

def opt(str):
    amino_dict={'A':0, 'R':1, 'N':2, 'D':3, 'C':4, 'Q':5, 'E':6, 'G':7, 'H':8, 'I':9,
                'L':10, 'K':11, 'M':12, 'F':13, 'P':14, 'S':15, 'T':16, 'W':17, 'Y':18,
                'V':19, 'B':20, 'Z':21, 'X':22, '*':23}
    sum=0;
    for i in range(len(str)):
        sum+=blosum[amino_dict[str[i]]][amino_dict[str[i]]]
    return sum

#@gets dict- keys=wmers. values: lists of Neighbers (type).
def extendNeighbers(map, db, query, X, w):
    optimizaion_dict = {}
    MSPs={}
    loading=0
    for wmer in map.keys():
    ##    loading += 1
    ##    print(loading/len(map.keys()))
        for neigh in map[wmer]: #for each neighber of the query wmer
            for j in range(len(neigh.idxs)): #j is instance of the neighber in DB
                extend_str(neigh.idxs[j], wmer.idx, db, query, X, w, MSPs, optimizaion_dict)
    return MSPs


#@tries to extend until score drops X below maximum score
#extends to the right until drops x under max, and then to the right.
def extend_str(DBidx, Qidx, db, query, X, w, MSPs, optimizaion_dict):
    q1=""
    db1=""
    for i in range(w):
        q1 = q1 + query[i+Qidx]
        db1 = db1 + db[i+DBidx]
    s=score(q1, db1)
    max=s
    max_query=q1
    max_db=db1
    maxQidx=Qidx
    maxDBidx=DBidx
    i=w
    proning_flag=0
    while (s>=max-X and (i+DBidx<len(db)) and (i+Qidx<len(query)) and proning_flag==0): #extend to the right
            if (optimizaion_flag==1):
                proning_flag=optimization(Qidx, DBidx, i, max, optimizaion_dict)
            q1=q1+query[i+Qidx]
            db1=db1+db[i+DBidx]
            s=score(q1, db1)
            if s>max:
                max=s
                maxQidx=Qidx
                maxDBidx=DBidx
                max_query=q1
                max_db=db1
            i += 1
    Qidx=maxQidx
    DBidx=maxDBidx
    q1=max_query
    db1=max_db
    s=max
    while (s>=max-X and DBidx>0 and Qidx>0 and proning_flag==0): #extends to the left
         q1=query[Qidx-1]+q1
         db1=db[DBidx-1]+db1
         DBidx -= 1
         Qidx -= 1
         s=score(q1, db1)
         if s>max:
             max = s
             maxQidx = Qidx
             maxDBidx = DBidx
             max_query = q1
             max_db = db1
    if (proning_flag==0):
        msp=MSP(max, maxQidx, maxQidx+len(max_query), maxDBidx, maxDBidx+len(max_db))
        if not (max_query in MSPs.keys()):
            MSPs[max_query]=[]
        MSPs[max_query].append(msp)

def optimization(Qidx, DBidx, i, max, optimizaion_dict):
    # note- the input wmer.idx order is increasing
    hit_params = [Qidx, DBidx, i, max] #all parameters for debbuging
    hit_str = str(Qidx + i - 1) + "-" + str(DBidx + i - 1)  # the end idxs of the MSP at that point
    if not (hit_str in optimizaion_dict.keys()):
        optimizaion_dict[hit_str] = [Qidx, DBidx, i, max]
    else:  # if we did extend already to that indexes in query and DB, its for sure contains the current. (see note)
        if optimizaion_dict[hit_str][3] > max:  # no point to extend anymore
            return 1
        else:
            optimizaion_dict[hit_str] = [Qidx, DBidx, i, max] #update the max params that extends to that idx in query and DB
    return 0

#@returns dict of data base. wmer: indexes
def preproc(db, w):
    window=""
    i=0
    dblen=len(db)
    dict={}
    for i in range(w):
        window=window + db[i] ##init with first wmer
    while i<(dblen-1):
        if not(window in dict.keys()):
            dict[window] = []
        dict[window].append(i-(w-1)) #w-1 because idx starts from 0
        i+=1
        window=window[1:]+db[i]
    if not(window in dict.keys()): #add last window
        dict[window] = []
    dict[window].append(i-(w-1))
    return dict

def reportMSP(map):
    for q in map.keys():
        checkduplicates(map[q])
        print ("<",q,"> :")
        for msp in map[q]:
            print("score: ",msp.score,", Query: ",msp.startQuery,"-",msp.endQuery,", DB:", msp.startDB,"-",msp.endDB)
        print("----------------------------")

def sortedReportMSPs(map):
    msps_list=[]
    for q in map.keys():
        checkduplicates(map[q])
    for q in map.keys():
        for msp in map[q]:
            msps_list.append(msp)
    sortedList = sorted(msps_list, key=lambda x: x.score, reverse=True)
    for msp in sortedList:
        print("score: ", msp.score, ", Query: ", msp.startQuery, "-", msp.endQuery, ", DB:", msp.startDB, "-",
              msp.endDB)
    print("----------------------------")
    return sortedList


def checkduplicates(list1):
    i=0
    list1strs=[]
    for j in range(len(list1)):
        list1strs.append(str(list1[j].score)+"-"+str(list1[j].startQuery)+"-"+str(list1[j].endQuery)+"-"+str(list1[j].startDB)+"-"+str(list1[j].endDB))
    while i<len(list1):
        if list1strs[i] in list1strs[(i+1):]:
            del list1[i]
            del list1strs[i]
            i -= 1
        i+=1

#@returns scores matrix
#@matrix columns and lines order: A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  *
def readBlosum(blosumFile):
    arr=blosumFile.read()
    arr=arr.split('\n')
    matrix=[[0 for x in range(24)] for y in range(24)]
    i=0
    for line in arr:
        if line[0]!='#' and line[0]!=' ':
            matrix[i]=line[2:]
            matrix[i]=matrix[i].strip(' ')
            matrix[i]=matrix[i].replace("  "," ")
            matrix[i]=matrix[i].split(' ')
            matrix[i] = list(map(int, matrix[i]))
            i+=1
    return matrix

#@returns name, seq.
def readFasta(fastaFile):
    arr=fastaFile.read()
    arr=arr.strip('>')
    arr = arr.rstrip()
    flag=0 # flag indicates if reads the name of the sequence
    seq=""
    name=""
    for i in range (len(arr)):
        if arr[i]=='\n' and flag==0:
            flag=1
            name=arr[:i]
            i += 1
        elif arr[i]=='\n' and flag==1:
            i+=1
        elif flag==1:
            seq=seq+arr[i]
    return name,seq

blosumfile = open('blosum_62.txt', 'r')
blosum=readBlosum(blosumfile)

if __name__ == '__main__':
	main()