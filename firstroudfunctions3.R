
remove(list=ls())
# K maximum  memory
K=5

########################################################
# function  defining  the  path  from  a  "nontrivial"  latent state L
#  to  the  next   "notrivial"  latent  state  with less  memory  
# nontrivial  are  those  states  such  that the next  transition is  nondeterministic
###########################################################

pathLM<-function(L){
a=-sum(L<=0)	
L[L<=0]=a	
Path=array(,c(2*K,K))
Path[1,]=L
if(a>(-K) ){
for(j in 1:(-a+1)){
	if((j+1)<=K){Path[j+1,]=c(L[(j+1):K],rep((a-1),j))}
	if((j+1)>K){Path[j+1,]=rep((a-1),j)}
	}
kkk=j+1
L=Path[j+1,]
while( sum(L[1:(-a+1)]==(a-1) )< (-a+1) ){
g1=L[1]
L=c(L[2:K],g1)	
kkk=kkk+1	
Path[kkk,]=L
	}
	}

Path= as.matrix(Path[1: sum(! is.na(Path[,1]) ), ])	
if(sum(dim(Path)==c(K,1))==2){Path=t(Path)}

	return(Path)
	}

#examples
L=c(-2,-2,1,2,1) #  the  random walker  remember  only  the  last  three  states (I am  using  -2  to  denote s2)
pathLM(L)
L=c(-2,-2,-1,-2,-1) 
pathLM(L)
L=c(1,2,1,2,1) 
pathLM(L)

########################################################
### function  defining  the  path  from  a  "nontrivial"  latent state
#  to  the  next   "notrivial"  latent  state  with  more  memory  
#  without  loss of  generalty  the RW  goes  to  0  
########################################################


pathMM<-function(L){
a=-sum(L<=0)	
L[L<=0]=a	
Path=array(,c(2*K,K))
Path[1,]=L
Path[2,]=c(L[2:K],0)
L=Path[2,]
if(a<(-1)){
for(j in 1:(-a-1)){Path[j+2,]=c(L[(j+1):K],rep((a+1),j))}

kkk=j+2
L=Path[j+2,]

while( sum(L[1:(-a-1)]==(a+1) )< (-a-1) ){
g1=L[1]
L=c(L[2:K],g1)	
kkk=kkk+1	
Path[kkk,]=L
	}
	}
Path= as.matrix(Path[1: sum(! is.na(Path[,1]) ), ])

	if(sum(dim(Path)==c(K,1))==2){Path=t(Path)}
	return(Path)	
	
	}


#examples
L=c(-3,-3,-3,2,1) #  the  random walker  remember  only  the  last 2  states (I am  using  -3  to  denote s3)
pathMM(L)
L=c(-2,-2,-1,-2,-1) 
pathMM(L)
L=c(1,2,1,2,1) 
pathMM(L)

#################################################################
# function  giving  all  possible  latent  paths  from  a "nontrivial"  latent
# state  to  a nontrivial  latent  state  with  the  addition   of   x  
#################################################################
   
allpl<-function(L,x){
	a=sum(L<=0)
	if(a>0){L[1:a]=-a}
	aa=K-a
	aaa=aa+1
	HG<-array(,c(aaa,(2*K)^2,K)) 
	for( fgt in  0:aa){		
		LL=L
		k=0
		if(fgt>0){for(j in 1:fgt){ m=pathLM(LL);LL=m[dim(m)[1],];kk=k+dim(m)[1];HG[fgt+1,(k+1):kk,]=m;k=kk}}
		m=pathMM(LL);kk=k+dim(m)[1];HG[fgt+1,(k+1):kk,]=m;k=kk
		HHG=HG[fgt+1,1:k,]
		R=rep(1,k)
		for(i in 2:k){ if (sum(HHG[i,]==HHG[i-1,])==K){R[i]=0}}
		HHG=HHG[R==1,]
		HHG=HHG+x*(HHG==0)
		HG[fgt+1,,]=rbind(HHG,array(,c((2*K)^2-dim(HHG)[1],K)))
		}
		return(HG)
		}

#examples
L=c(-3,-3,-3,1,2) #  the  random walker  remember  only  the  last 2  states (I am  using  -3  to  denote s3)
ex=allpl(L,7)
ex[1,1:7,]
ex[2,1:12,]
ex[3,1:17,]

L=c(-6,-6,-6,-6,-6) #  the  random walker  remember  only  the  last 2  states (I am  using  -3  to  denote s3)
ex=allpl(L,7)
ex[1,,]
dim(ex)

L=c(6,6,6,6,6) #  the  random walker  remember  only  the  last 2  states (I am  using  -3  to  denote s3)
ex=allpl(L,7)
ex[1,1:2,]
ex[2,1:7,]
ex[3,1:17,]
ex[4,1:22,]
ex[5,1:27,]
ex[6,1:32,]

dim(ex)

########################################################
# a  function that  (i)  when  zzz=0
# return  the  transition matrix  (prior) for  the latent states when  \mathcal{X}  is  singleton
# and (ii)  when zzz=1   return  the corresponding  list   of   latent states  
########################################################

  
tm<-function(alpha,zzz){
LS=array(0,c(2,K))
for(k  in  0:K){LS=rbind(LS,as.vector(c(rep(0,k), rep(1,K-k))))}
LS=LS[3:length(LS[,1]),]

LS2=array(0,c(2,K))
for(i  in 1:length(LS[,1])){ LS2=rbind(LS2,pathLM(LS[i,]));LS2=rbind(LS2,rep(-1,K))}

LS3=array(0,c(2,K))
for(i  in 1:length(LS[,1])){ LS3=rbind(LS3,pathMM(LS[i,]));LS3=rbind(LS3,rep(-1,K))}
LS2=LS2[3:length(LS2[,1]),]	
LS3=LS3[3:length(LS3[,1]),]	
LS3=LS3+(LS3==0) 
LS2=rbind(rep(-1,K),LS2)
LS3=rbind(rep(-1,K),LS3)
LS5=rbind(LS2,LS3)

for(i in 1:K){LS5=LS5[order(LS5[,i]),]}
R=rep(0, length(LS5[,1]))
for(i in 1:(length(LS5[,1])-1)){if(!(sum(LS5[i,]==LS5[i+1,])==K)){R[i]=1}}
R[i+1]=1;LS5=LS5[R==1,];R=rep(0, length(LS5[,1]))
for(i in  1:K){R=R+(!(LS5[,i]==-1))}	
LS5=LS5[R>0,]

TRAN=array(0,c(length(LS5[,1]),length(LS5[,1])))

LLS2=rep(0,length(LS2[,1]))
for(i in  1:length(LS2[,1])){
	R=c(1:length(LS5[,1]))
	for(j in 1:K){R=R*(LS5[,j]==LS2[i,j])}
	LLS2[i]=sum(R)}
	
LLS3=rep(0,length(LS3[,1]))
for(i in  1:length(LS3[,1])){
	R=c(1:length(LS5[,1]))
	for(j in 1:K){R=R*(LS5[,j]==LS3[i,j])}
	LLS3[i]=sum(R)}
	
for(i  in  1:(length(LLS2)-1)){
	if(LLS2[i]>0){if(LLS2[i+1]>0){TRAN[LLS2[i],LLS2[i+1]]=1}}}	
for(i  in  1:(length(LLS3)-1)){
	if(LLS3[i]>0){if(LLS3[i+1]>0){TRAN[LLS3[i],LLS3[i+1]]=1}}}
	
ids=(c(1:length(TRAN[,1])))[(TRAN%*%rep(1,length(TRAN[,1])))>1]

for(i  in  1: length(ids)){
	id2=(c(1:length(TRAN[,1])))[TRAN[ids[i],]>0]
	TRAN[ids[i],id2[1]]=alpha;TRAN[ids[i],id2[2]]=1-alpha
if( sum( LS5[id2[1],]>0)<sum( LS5[id2[2],]>0)){TRAN[ids[i],id2[1]]=1-alpha;TRAN[ids[i],id2[2]]=alpha}}

if(zzz==0){return(TRAN)}
if(zzz==1){return(LS5)}
}

#example 
K=3
alpha=0.2
zzz=1
tm(alpha,zzz)
zzz=0
tm(alpha,zzz)

#############################################
# function  giving  the stationary  distribution###  
#of  the latent states when  \mathcal{X}  is  singleton##
######################################################

st<-function(alpha){
zzz=0
B=tm(alpha,zzz)
a=eigen(t(B))
a=abs(t(as.vector(a$vectors[,1])/sum((a$vectors[,1]))))
return(a)
}

#example
K=3
alpha=0.2
zzz=0
B=tm(alpha,zzz)
a=st(alpha)
max(a%*%B-a)
min(a%*%B-a)

#################################################################
# iw: function giving  the initial  weight  of a  latent  state  L (vector  of  length K) 
#################################################################

iw<-function(L,BETA,N,Listelementary,weightselementary){
L1=L;L1[L1>0]=1
R=rep(1,length(Listelementary[,1]))
for(i in 1:(dim(Listelementary)[2])){R=R* (Listelementary[,i]==L1[i])}
w=sum(weightselementary*R)
w=w*BETA*(N^(-sum(L>0)))
return(w)}

#example
K=3
alpha=0.2;
BETA=2
L=c(1,2,-2)
N=5
zzz=1
Listelementary=tm(alpha,zzz)
zzz=0
B=tm(alpha,zzz)
weightselementary=st(alpha)
max(weightselementary%*%B-weightselementary)
min(weightselementary%*%B-weightselementary)
iw(L,BETA,N,Listelementary,weightselementary)

#################################################################
# function giving  the initial  weight  of an edge from a  latent  state  Lp  to  a  
#  latent state  La.  N  is  the  total  number  of  simbols 
# BETA  is  a prior  parameter   
#################################################################

ie<-function(Lp,La,BETA,N,Listelementary,weightselementary,B){
w=0	
L1=Lp;L1[L1>0]=1
R=rep(1,length(Listelementary[,1]))
for(i in 1:(dim(Listelementary)[2])){R=R* (Listelementary[,i]==L1[i])}
R1=sum(R*c(1:length(R)))

L1=La;L1[L1>0]=1
R=rep(1,length(Listelementary[,1]))
for(i in 1:(dim(Listelementary)[2])){R=R* (Listelementary[,i]==L1[i])}
R2=sum(R*c(1:length(R)))

if(R1>0 & R2>0){if(  B[R1,R2]>0  & sum(Lp[2:length(Lp)]==La[1:(length(La)-1)])==(-1+length(Lp))){
w=iw(Lp,BETA,N,Listelementary,weightselementary)*B[R1,R2]
ww=w
if(La[length(La)]>0){ww=w/N}
if(La[length(La)]>0){if(sum(B[R1,]>10^(-10))==1){
		if(!(sum(Lp<0)==length(Lp)) ){ 
			ww=w
			if(! (Lp[1]==La[length(La)])  ){ ww=0 }}}}
}}
return(ww)}

#example
K=3
alpha=0.2;
BETA=2
Lp=c(-2,-2,1)
La=c(-2,1,2)
N=3
zzz=1
Listelementary=tm(alpha,zzz)
zzz=0
B=tm(alpha,zzz)
weightselementary=st(alpha)
ie(Lp,La,BETA,N,Listelementary,weightselementary,B)

#double check (you can ignore):

L5=array(,c(10000,11))
N=3
kkkk=0
for(i in 1: length(Listelementary[,1])){
L1=	Listelementary[i,]
for(j in 1: length(Listelementary[,1])){
L2=Listelementary[j,]
if(B[i,j]>0){
L3=c(L1,L2[K])
for(i1 in 1:3){
for(i2 in 1:3){
for(i3 in 1:3){
for(i4 in 1:3){
L4=L3
if(L4[1]>0){L4[1]=i1}
if(L4[2]>0){L4[2]=i2}
if(L4[3]>0){L4[3]=i3}
if(L4[4]>0){L4[4]=i4}
kkkk=kkkk+1
L5[kkkk,1:4]=L4
}}}}
	}}
	}
L5=L5[1:kkkk,]
L5=L5[order(L5[,1]),]	
L5=L5[order(L5[,2]),]	
L5=L5[order(L5[,3]),]	
L5=L5[order(L5[,4]),]	
R=rep(1,length(L5[,1]))
for(i in 2:length(L5[,1])){
	if(sum(L5[i,1:4]==L5[i-1,1:4])==4){R[i]=0}}
	L5=L5[R>0,]
for(i in 1:length(L5[,1])){ 
	Lp=L5[i,1:3]
	La=L5[i,2:4]
	L5[i,5]= ie(Lp,La,BETA,N,Listelementary,weightselementary,B)
 }	
sum(L5[,5]>0   )
L5=L5[L5[,5]>0,]


for(i in 1:length(L5[,1])){ 
	Lp=c(L5[i,4],L5[i,3],L5[i,2])
	La=c(L5[i,3],L5[i,2],L5[i,1])
	L5[i,6]= ie(Lp,La,BETA,N,Listelementary,weightselementary,B)
 }	
L5[,5]/L5[,6]
sum(L5[,5])
max(L5[,5]/L5[,6])
min(L5[,5]/L5[,6])

for(i in 1:length(L5[,1])){ 
	Lp=L5[i,1:3]
	kzv=0
for(j in 1:length(L5[,1])){ 

if(sum(Lp==L5[j,1:3])==3){kzv=kzv+L5[j,5]}

}
L5[i,7]=kzv
}

for(i in 1:length(L5[,1])){ 
	Lp=L5[i,1:3]
	kzv=0
for(j in 1:length(L5[,1])){ 

if(sum(Lp==L5[j,2:4])==3){kzv=kzv+L5[j,5]}

}
L5[i,8]=kzv
}
L5[,7]/L5[,8]
max(L5[,7]/L5[,8])
min(L5[,7]/L5[,8])


for(i in 1:length(L5[,1])){ 
	Lp=L5[i,1:3]
   L5[i,9]=iw(Lp,BETA,N,Listelementary,weightselementary)
}
L5[,9]/L5[,8]
max(L5[,9]/L5[,8])
min(L5[,9]/L5[,8])

sort(L5[,5]/L5[,9])

for(i in 1:length(L5[,1])){ 
	Lp=L5[i,1:4]
	Lpp=c(Lp[4],Lp[3],Lp[2],Lp[1])
	kzv=0
for(j in 1:length(L5[,1])){ 

if(sum(Lpp==L5[j,1:4])==4){kzv=kzv+1}

}
L5[i,10]=kzv
}

L5[,10]
max(L5[,10])
min((L5[,10]))

L5[,11]=0
h=(c(-3,-3,-3))
for(i in 1:1000000){
  g=trunc(runif(1,1,1+length(L5[,10])))	
	if(sum(L5[g,1:3]==h)==3){ 
		h=L5[g,2:4]
		for(j  in 1:length(L5[,10])){		
			if(sum(L5[j,1:3]==h)==3){ L5[j,11]=1  }
			}}}
			min(L5[,11])

#################################################################
# fuunction that  identify   the   line   of  a  matrix  M (K\times watever)  that  
# mathches a  vector  V (K \times 1)
#  the  raws  in  M  are assumed  distinct  
#################################################################
match<-function(M,V){
	R=rep(1,length(M[,1]))
	for(i in 1:(dim(M)[2])){R=R*(M[,i]==V[i])}
	return(sum(R*c(1:length(R))))
	}

#example 
V=c(1,2,2)
M=rbind(c(1,1,1),c(2,2,2),c(3,3,3),c(1,2,2))
match(M,V)
M=rbind(c(1,1,1),c(2,2,2),c(3,3,3),c(1,2,6))
match(M,V)

########################################################
# reverse  a vector
########################################################
 rever<-function(L){
 	LL=L
 	for(i in 1:(length(L))){LL[(length(L))-i+1]=L[i] }
 	return(LL)
 	}
#example
rever(c(1,2,3))

##########################################################
# functions that  update  summary  statistics
###########################################################

updsss<-function(OldM,OldW,New){
	id1=match(OldM,New)
    if(id1>0){OldW[id1]=1+OldW[id1]}
	    if(id1==0){
	    	OldM=as.matrix(rbind(OldM,t(as.matrix(New))))
	    	OldW=c(OldW,1)
	    	}
  return(cbind(OldM, as.matrix(OldW)))
	}
#example
OldM=rbind( c(-3,-3,-3),c(-3,-3,1),c(-3,1,-2),c(1,-2,-2), c(-2,-2,1))
OldW=rep(1,5)
New=c(1,1,1)
updsss(OldM,OldW,New)
New=c(-3,-3,1)
updsss(OldM,OldW,New)

#################################################################
# function giving  the log-transition  probability   to  La  from the  
# state  Lp (wich  is  assumed non trivial)  given  the  sufficient
# summary statistics LE (listof  edges -all distinct-) LEW (list of  edges weights-all distinct-)
# LS (list of latent states-all distinct-) LSW  (list of latent states weghts-all distinct-)
#################################################################
prtlsle<-function(La,Lp,LE,LEW,LS,LSW,BETA,N,Listelementary,weightselementary,B){
denom1=0
nom1=0
ed=c(Lp,La[length(La)])
id1=match(LE,ed)
id2=match(LS,Lp)
if(id1>0){	nom1=LEW[id1]}
if(id2>0){	denom1=LSW[id2]}
denom2=iw(Lp,BETA,N,Listelementary,weightselementary)
nom2=ie(Lp,La,BETA,N,Listelementary,weightselementary,B)
return(log(nom1+nom2)-log(denom1+denom2))
	}

#example
K=3
alpha=0.2;
BETA=2

a1=allpl(c(-3,-3,-3),1)[1,1:5,]
a2=allpl(a1[5,],2)[1,1:5,]
a3=allpl(a2[5,],2)[2,1:8,]
PF=rbind(a1,a2[2:5,],a3[2:8,])

LS1=array(0,c(2,3));LSW1=c(0,0)
LE1=array(0,c(2,4));LEW1=c(0,0)

for(i in 2:16){
A=updsss(LS1, LSW1, PF[i-1,])
LS1=A[,1:(dim(A)[2]-1)]
LSW1=A[,(dim(A)[2])]
A=updsss(LS1, LSW1, rever(PF[i,]))
LS1=A[,1:(dim(A)[2]-1)]
LSW1=A[,(dim(A)[2])]
A=updsss(LE1, LEW1, c(PF[i-1,],PF[i,length(PF[i,])]))
LE1=A[,1:(dim(A)[2]-1)]
LEW1=A[,(dim(A)[2])]
A=updsss(LE1, LEW1, rever(c(PF[i-1,],PF[i,length(PF[i,])])))
LE1=A[,1:(dim(A)[2]-1)]
LEW1=A[,(dim(A)[2])]
}

Lp=PF[16,]
La=c(2,2,3)
LEW=LEW1;LE=LE1;LS=LS1;LSW=LSW1
prtlsle(La,Lp,LE,LEW,LS,LSW,BETA,N,Listelementary,weightselementary,B)
#################################################################
#function  giving    the log-probability,  given  the   sufficient
# summary statistics LE (listof  edges) LEW (list of  edges weights)
# LS (list of latent states) LSW  (list of latent states weghts)
#  of  a future  path PF assuming  the  current  state  is PF[1,]
#################################################################

prtplsle<-function(PF,LE,LEW,LS,LSW,BETA,N,Listelementary,weightselementary,B){		LE1=LE;LEW1=LEW;LS1=LS;LSW1=LSW;
		h=prtlsle(PF[2,],PF[1,],LE,LEW,LS,LSW,BETA,N,Listelementary,weightselementary,B)
if((dim(PF)[1])>2){
	for(i in 2:(dim(PF)[1]-1)){ 
#update
A=updsss(LS1, LSW1, PF[i-1,])
LS1=A[,1:(dim(A)[2]-1)]
LSW1=A[,(dim(A)[2])]
A=updsss(LS1, LSW1, rever(PF[i,]))
LS1=A[,1:(dim(A)[2]-1)]
LSW1=A[,(dim(A)[2])]
A=updsss(LE1, LEW1, c(PF[i-1,],PF[i,length(PF[i,])]))
LE1=A[,1:(dim(A)[2]-1)]
LEW1=A[,(dim(A)[2])]
A=updsss(LE1, LEW1, rever(c(PF[i-1,],PF[i,length(PF[i,])])))
LE1=A[,1:(dim(A)[2]-1)]
LEW1=A[,(dim(A)[2])]
#compute
h=h+prtlsle(PF[i+1,],PF[i,],LE1,LEW1,LS1,LSW1,BETA,N,Listelementary,weightselementary,B)
		}
			}
			return(h)
}

#example
K=3
alpha=0.2;
BETA=2

a1=allpl(c(-3,-3,-3),1)[1,1:5,]
a2=allpl(a1[5,],2)[1,1:5,]
a3=allpl(a2[5,],2)[2,1:8,]
PF=rbind(a1,a2[2:5,],a3[2:8,])

LS1=array(0,c(2,3));LSW1=c(0,0)
LE1=array(0,c(2,4));LEW1=c(0,0)

for(i in 2:16){
A=updsss(LS1, LSW1, PF[i-1,])
LS1=A[,1:(dim(A)[2]-1)]
LSW1=A[,(dim(A)[2])]
A=updsss(LS1, LSW1, rever(PF[i,]))
LS1=A[,1:(dim(A)[2]-1)]
LSW1=A[,(dim(A)[2])]
A=updsss(LE1, LEW1, c(PF[i-1,],PF[i,length(PF[i,])]))
LE1=A[,1:(dim(A)[2]-1)]
LEW1=A[,(dim(A)[2])]
A=updsss(LE1, LEW1, rever(c(PF[i-1,],PF[i,length(PF[i,])])))
LE1=A[,1:(dim(A)[2]-1)]
LEW1=A[,(dim(A)[2])]
}
LEW=LEW1;LE=LE1;LS=LS1;LSW=LSW1
past=PF
PF=allpl(c(-1,2,2),1)[2,1:8,]

prtplsle(PF,LE,LEW,LS,LSW,BETA,N,Listelementary,weightselementary,B)


#################################################################
#function  giving    the log-probabilities,  given  the   past  LE,LEW,LS,LSW,  with non trivial final state
# Lp of  all  the future  paths  that  terminates   in  x
#################################################################

prtab<-function(x,Lp,LE,LEW,LS,LSW,BETA,N,Listelementary,weightselementary,B){
L=Lp
U=allpl(L,x)
g=dim(U)
vec=rep(0,dim(U)[1])
for(i in  1:(dim(U)[1])){
h=sum(! is.na(U[i,,1]))	
vec[i]= prtplsle(U[i,1:h,],LE,LEW,LS,LSW,BETA,N,Listelementary,weightselementary,B)
	}
	return(vec)
	}

#x;Lp=CLS2[i,];LE;LEW=LEW2[,i];LS;LSW=LSW2[,i]

#example
K=3
alpha=0.2;
BETA=2

a1=allpl(c(-3,-3,-3),1)[1,1:5,]
a2=allpl(a1[5,],2)[1,1:5,]
a3=allpl(a2[5,],2)[2,1:8,]
PF=rbind(a1,a2[2:5,],a3[2:8,])

LS1=array(0,c(2,3));LSW1=c(0,0)
LE1=array(0,c(2,4));LEW1=c(0,0)

for(i in 2:16){
A=updsss(LS1, LSW1, PF[i-1,])
LS1=A[,1:(dim(A)[2]-1)]
LSW1=A[,(dim(A)[2])]
A=updsss(LS1, LSW1, rever(PF[i,]))
LS1=A[,1:(dim(A)[2]-1)]
LSW1=A[,(dim(A)[2])]
A=updsss(LE1, LEW1, c(PF[i-1,],PF[i,length(PF[i,])]))
LE1=A[,1:(dim(A)[2]-1)]
LEW1=A[,(dim(A)[2])]
A=updsss(LE1, LEW1, rever(c(PF[i-1,],PF[i,length(PF[i,])])))
LE1=A[,1:(dim(A)[2]-1)]
LEW1=A[,(dim(A)[2])]
}
LEW=LEW1;LE=LE1;LS=LS1;LSW=LSW1
x=1
prtab(x,Lp,LE,LEW,LS,LSW,BETA,N,Listelementary,weightselementary,B)
###########################################################################
##################################################################
#  function  for  particle  filter  
#  input:
#  1)  Data  a  sequence   of  simbols in  {1,....,N}  (these  are the  data) 
#  2)  step,  i.e. step of  the smc slgorithm, it  ranges  in  1,....,length(D)   
#  3)  as  in the  previous  functions: BETA,N,Listelementary,weightselementary,B  
#################################################################
###########################################################################
# LS (K\times many)  list  of  distinct  latent states
# LE ((K+1)\times many)  list  of  distinct  latent  edges 
# LEW  weights  of  LE  many \times #particles    
# LSW  weights  of  LS  many \times #particles 
# CLS  current  latent  state   many\times K
###########################################################################

pf<-function(Data, step,CLS,LE,LEW,LS,LSW,BETA,N,Listelementary,weightselementary,B){
	#resampling
	x=Data[step]
	w=rep(0,length(CLS[,1])) 
	for(i in  1: length(w)){
	w[i]=sum(exp(prtab(x,CLS[i,],LE,LEW[,i],LS,LSW[,i],BETA,N,Listelementary,weightselementary,B)))}    
   w=w/sum(w)	
   index=sample(length(w),length(w),replace=TRUE,prob=w)
   	
CLS2=CLS
LEW2=LEW
LSW2=LSW

for(i in 1: length(w)){ 
	CLS2[i,]= CLS[index[i],]
	LEW2[,i]= LEW[,index[i]]
		LSW2[,i]= LSW[,index[i]]
	   }

# length the pat
l=rep(0,length(w))
for(i in 1: length(w)){ 
pr=(exp(prtab(x,CLS2[i,],LE,LEW2[,i],LS,LSW2[,i],BETA,N,Listelementary,weightselementary,B)))
l[i]=sample(length(pr),1, prob=pr/(sum(pr)));
}

CLS3=CLS
# update  lists  and weights
for(i in 1: length(l)){ 
	PP=(allpl(CLS2[i,],x))[l[i],,];
	PP=PP[1: (sum(!is.na(PP[,1]))),];
	CLS3[i,]=PP[length(PP[,1]),]
#update
 LEW1=LEW2[,i]
 LSW1=LSW2[,i]
 LS1=LS
 LE1=LE 
 for(j in 2:length(PP[,1])){
A=updsss(LS1, LSW1, PF[j-1,])
LS1=A[,1:(dim(A)[2]-1)]
LSW1=A[,(dim(A)[2])]
A=updsss(LS1, LSW1, rever(PF[j,]))
LS1=A[,1:(dim(A)[2]-1)]
LSW1=A[,(dim(A)[2])]
A=updsss(LE1, LEW1, c(PF[j-1,],PF[j,length(PF[j,])]))
LE1=A[,1:(dim(A)[2]-1)]
LEW1=A[,(dim(A)[2])]
A=updsss(LE1, LEW1, rever(c(PF[j-1,],PF[j,length(PF[j,])])))
LE1=A[,1:(dim(A)[2]-1)]
LEW1=A[,(dim(A)[2])]
} 
h=-length(LEW2[,i])+length(LEW1)
LE=LE1
LEW2=rbind(LEW2,array(0,c(h,(dim(LEW2)[2]))))
LEW2[,i]=LEW1  
  
h=-length(LSW2[,i])+length(LSW1)
LS=LS1
LSW2=rbind(LSW2,array(0,c(h,(dim(LSW2)[2]))))
LSW2[,i]=LSW1  
   }

LSW<<-LSW2
LE<<-LE
LS<<-LS
LEW<<-LEW2
CLS<<-CLS3
}

#example  (run  all  the  code  above to  get  the  followig  objects:
Data=c(1,1,2,3)
CLS=array(-3,c(30,3))
LE=array(0,c(2,4))
LEW=array(0,c(2,30))
LS=array(0,c(2,3))
LSW=array(0,c(2,30))
for(step  in 1:4){
pf(Data, step,CLS,LE,LEW,LS,LSW,BETA,N,Listelementary,weightselementary,B)
}

CLS
LE
LEW
LS
LSW