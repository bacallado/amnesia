from numpy import *
from numpy.random import *

p = arange(1,9) 
m = []

for i in range(100000):
	m.append(sum(p>0))	
	#print m[-1]
	l = randint(0,7)
	a = p[l]
	b = p[l+1]
	p[l]=-b
	p[l+1]=-a


m = array(m)

tot2 = 0
tot22 = 0
tot222 = 0
tot2222 = 0
for i in range(len(m)-3):
	if m[i]==2:
		tot2 +=1
		if m[i+1]==2:
			tot22 +=1
			if m[i+2]==2:
				tot222 +=1
				if m[i+3]==2:
					tot2222 +=1

print sum(m==0)
print sum(m==2)
print sum(m==4)
print sum(m==6)
print sum(m==8)

f = open("SignedPermutations.dat","w")
savetxt(f,m,fmt='%i')


print float(tot2222)/float(tot222)

print float(tot222)/float(tot22)

print float(tot22)/float(tot2)


