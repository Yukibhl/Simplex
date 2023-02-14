import numpy as np

class pls :
    def __init__(self,maxz,nbrc,b,A):
        self.maxz=maxz
        self.nbrc=nbrc
        self.b=b
        self.A=A
    def makeT0(self):
        T0=np.zeros((self.nbrc+1,len(self.maxz)+self.nbrc+2))
        T0[:self.nbrc,len(self.maxz):self.nbrc+len(self.maxz)]=np.identity(self.nbrc)
        T0[:self.nbrc,:len(self.maxz)]=self.A     
        for n in range(len(self.maxz)):
            T0[self.nbrc,n]=-self.maxz[n]
        T0[self.nbrc,len(self.maxz)+self.nbrc]=1
        for j in range(self.nbrc):
                T0[j,len(self.maxz)+self.nbrc+1]=self.b[j]
        return T0

def pivot(pls,T0):
    entrant=None
    for i in range(len(pls.maxz)+pls.nbrc):
        if T0[pls.nbrc,i]<0 :
            entrant=i
            break
    s=[]
    
    if entrant==None :return None
    else :
        for i in range(pls.nbrc):
         if T0[i,entrant]>0 :
            s.append(T0[i,-1]/T0[i,entrant])
         else :
            s.append('*')
        o=set(s)
        if allpositif(s) and len(o)!=1 :
            lignepivot=s.index(min(s))
            return [lignepivot,entrant]
    
        if allpositif(s)==False and len(o)!=1:
         m=list(o)
         m.remove('*')
         lignepivot=s.index(min(m))
         return [lignepivot,entrant]
        if len(o)==1 and allpositif(s):
         u=[]
         if '*' in o :return None
         for i in range(pls.nbrc):
            u.append(T0[i,-1])
         lignepivot=u.index(min(u))
         return[lignepivot,entrant]

def allpositif(L) :
    for i in L :
        if type(i)!=str:
          if i < 0  :
            return False
        else :
            return False
    return True

def existzero(L):
    s=[]
    for i in L :
        if i==0 :
            s.append(i)
    if s==[] : return False
    else :
        return s


def simplexphase2(pls,T0):
        piv=pivot(pls,T0)
        while ( pivot(pls,T0)!=None) :
            piv=pivot(pls,T0)
            P=np.identity(pls.nbrc+1)
            for i in range(pls.nbrc+1):
                if i!=piv[0] :
                    P[i,piv[0]]=-T0[i,piv[1]]/T0[piv[0],piv[1]]
                else :
                    P[i,piv[0]]=1/T0[piv[0],piv[1]]
            T0=np.dot(P,T0)
        return T0
    

    
def ibase(pls,T0):
    l=[]
    I=np.identity(pls.nbrc)
    k=list(I)
    for i in range(len(k)) :
        k[i]=list(k[i])
    for j in range(len(pls.maxz)+pls.nbrc):
        if list(T0[:pls.nbrc,j])in k:
            l.append(j)
    return l
def simplexphase1(pls,T):
    K=[[0] for i in range(pls.nbrc+1)]
    Taux=np.append(T[:pls.nbrc+1,:pls.nbrc+len(pls.maxz)],np.array(K),axis=1)
    Taux=np.append(Taux,T[:pls.nbrc+1,len(pls.maxz)+pls.nbrc:len(pls.maxz)+pls.nbrc+2],axis=1)
    s=[]
    for i in range(pls.nbrc) :
        if T[i,pls.nbrc+len(pls.maxz)+1]<0 :
            s.append(i)
    A=[[0] for i in range(pls.nbrc+1)]
    for i in range(pls.nbrc):
        if i in s :
            A[i][0]=-1
    Taux=np.append(np.array(A),Taux,axis=1)
    l=[0 for i in range(len(s)+len(pls.maxz)+pls.nbrc+3)]
    for t in range(len(s)) :
        l[t]=1
    l[len(s)+len(pls.maxz)+pls.nbrc]=1
    Taux=np.append(Taux,np.array([l]),axis=0)
    entrant=0
    while entrant <len(s) or allpositif(Taux[pls.nbrc+1,:len(s)+len(pls.maxz)+pls.nbrc])==False:
        p=[]
        if entrant>=len(s) :
            for i in range(len(s)+len(pls.maxz)+pls.nbrc):
                if Taux[pls.nbrc+1,i]<0:
                    entrant=i
                    break
        for i in range(pls.nbrc):
            if Taux[i,entrant]>0 :
                p.append(Taux[i,len(s)+pls.nbrc+len(pls.maxz)+2]/Taux[i,entrant] )
            else :
                p.append(np.inf)
        if len(p)==p.count(np.inf):
            p=[]
            for i in range(pls.nbrc):
                p.append(Taux[i,len(s)+pls.nbrc+len(pls.maxz)+2] )
        lignepivot=p.index(min(i for i in p if i!=np.inf))
        piv=[lignepivot,entrant]
        P=np.identity(pls.nbrc+2)
        for i in range(pls.nbrc+2):
                if i!=piv[0] :
                    P[i,piv[0]]=-Taux[i,piv[1]]/Taux[piv[0],piv[1]]
                else :
                    P[i,piv[0]]=1/Taux[piv[0],piv[1]]
        Taux=np.dot(P,Taux)
        entrant+=1
    return Taux

maxz=input('Maxz :')
maxz=maxz.split(' ')
for w in range(len(maxz)):
    maxz[w]=float(maxz[w])   
nbrc=int(input('Number of ineqalities :'))
b=input('Matrix b :')
b=b.split(' ')
for w in range(len(b)):
    b[w]=float(b[w])
A=np.zeros((nbrc,len(maxz)))
for m in range(nbrc):
            y='Inequality n '+str(m+1)+' :'
            t=input(y)
            t=t.split(' ')
            for j in range(len(t)):
                    A[m,j]=float(t[j])
p=pls(maxz,nbrc,b,A)
T0=p.makeT0()
s=[]
for i in range(p.nbrc) :
        if T0[i,p.nbrc+len(p.maxz)+1]<0 :
            s.append(i)
if not allpositif(b) :
    Taux=simplexphase1(p,T0)
    if Taux[-1,-1] <0 :
        print('No admissible solutions')
    elif Taux[-1,-1]==0 :
        Aaux=Taux[:p.nbrc,len(s):len(p.maxz)+1]
        baux=Taux[:p.nbrc+1,-1]
        maxzaux=Taux[p.nbrc,len(s):len(p.maxz)+1]
        pne=pls(maxzaux,nbrc,baux,Aaux)
        T0aux=pne.makeT0()
        T0aux[-1,:len(maxzaux)]=-T0aux[-1,:len(maxzaux)]
        for i in range(pne.nbrc):
            T0aux[:pne.nbrc+1,i+len(maxzaux)]=Taux[:nbrc+1,i+len(s)+len(maxzaux)]
        O=simplexphase2(pne,T0aux)
        ivb=ibase(pne,T0aux)
        print('Possible solutions :')
        z=0
        for i in ivb :
            print('x'+str(i+1)+'= '+str(O[z,-1]))
            z=z+1
        print('Z optimal'+'= '+str(O[p.nbrc,-1]))
else :
    T=simplexphase2(p,T0)
    ivb=ibase(p,T)
    ivnb=[]
    for i in range(len(maxz)+nbrc):
        if i not in ivb :
            ivnb.append(i)
    n=[]
    for q in ivnb :
     n.append(T[p.nbrc,q])

    if not allpositif(n) :               
               print('No optimal solutions')
    else :
     if existzero(list(T[nbrc,i] for i in ivnb))==False:
       if existzero(list(T[:nbrc,-1]))==False :                   
        print('Only one optimal solution :')
        z=0
        for i in ivb :
            print('x'+str(i+1)+'= '+str(T[z,-1]))
            z=z+1
       else :
        print('A degenered solution : ')       
        z=0
        for i in ivb :
            print('x'+str(i+1)+'= '+str(T[z,-1]))
            z=z+1
     else :                                                           
        X1=T[:nbrc,-1]
        entrant=ivnb[0]
        f=[]
        for i in range(nbrc):
             f.append(T[i,-1]/T[i,entrant])
        mi=f.index(min(f))
        piv=(mi,entrant)
        P=np.identity(p.nbrc+1)
        for i in range(p.nbrc+1):
                if i!=piv[0] :
                    P[i,mi]=-T[i,entrant]/T[mi,entrant]
                else :
                    P[i,mi]=1/T[mi,entrant]
        T=np.dot(P,T)
        X2=T[:nbrc,-1]
        ch='Various solutions : \n a*'+str(X1)+'+(1-a)*'+str(X2)+'\n'+'where a in [0,1]'
        print(ch)
     print('Z optimal'+'= '+str(T[p.nbrc,-1]))

        
        
        
        


