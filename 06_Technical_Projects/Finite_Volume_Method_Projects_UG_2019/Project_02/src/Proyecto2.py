print('\033[H\033[J')
import numpy as np
import matplotlib.pyplot as plt
from Problema2 import Problema2
from scipy import interpolate


n=10
gamma=1
L=1
dx=dy=L/n

esquema='CDS'

if esquema=='UDS':
    alpha=1
    betha=0
if esquema=='CDS':
    alpha=0
    betha=1

Dcons=np.zeros(((n)),float)
Dphi=np.zeros(((n)),float)
x=np.zeros(((n)),float)


sol=Problema2(n,dx,np,gamma,alpha,betha)
for k in range(0,n):
    Dcons[k]=sol[k*(n)+k]
    Dphi[k]=sol[-k*n-n+k]
    x[k]=dx/2+dx*k

n=50
dx=dy=L/n

Dcons2=np.zeros(((n)),float)
Dphi2=np.zeros(((n)),float)
x2=np.zeros(((n)),float)

sol2=Problema2(n,dx,np,gamma,alpha,betha)
for k in range(0,n):
    Dcons2[k]=sol2[k*(n)+k]
    Dphi2[k]=sol2[-k*n-n+k]
    x2[k]=dx/2+dx*k
    
n=100
dx=dy=L/n

Dcons3=np.zeros(((n)),float)
Dphi3=np.zeros(((n)),float)
x3=np.zeros(((n)),float)

sol3=Problema2(n,dx,np,gamma,alpha,betha)
for k in range(0,n):
    Dcons3[k]=sol3[k*(n)+k]
    Dphi3[k]=sol3[-k*n-n+k]
    x3[k]=dx/2+dx*k
#    
#
A=np.zeros(((n,n)),float) 
arr = sol3.reshape((n,n))
for m in range (0,n):
    A[m]=arr[n-1-m]    
    
plt.figure (1)
plt.plot(x,Dphi, label = "n=10", color = 'red')
plt.plot(x2,Dphi2, label = "n=50", color = 'blue')
plt.plot(x3,Dphi3, label = "n=100", color = 'green')
plt.legend(loc="upper right")
plt.hold(True)
plt.grid(True)
plt.grid(color = '0.5', linestyle = '--', linewidth = 1)
plt.xlabel(r"$X$", fontsize = 12, color = 'black')
plt.ylabel(r"$\phi$", fontsize = 15, color = 'black')
plt.title('$CDS\ \ \Gamma=1$',fontsize = 12, color = 'black', verticalalignment = 'baseline', horizontalalignment = 'center')
plt.savefig('CDS_Gamma=1.pdf')
plt.show()


xm=np.arange(dx/2,L,dx)
ym=np.arange(dx/2,L,dx)
xx,yy=np.meshgrid(xm,ym)
f=interpolate.interp2d(xm,ym,sol)
im = plt.imshow(A,cmap='jet')
plt.colorbar()
plt.hsv
plt.xlabel(r"$X\ [\%]$", fontsize = 12, color = 'black')
plt.ylabel(r"$Y\ [\%]$", fontsize = 15, color = 'black')
plt.title('$CDS\ \ \Gamma=1$',fontsize = 12, color = 'black', verticalalignment = 'baseline', horizontalalignment = 'center')
plt.savefig('CDS_Gamma=1_Dis.pdf')











