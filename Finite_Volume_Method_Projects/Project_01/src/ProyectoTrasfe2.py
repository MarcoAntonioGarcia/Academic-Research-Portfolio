print('\033[H\033[J')
import numpy as np
import matplotlib.pyplot as plt
from gauss_seidel import gauss_seidel
n=100 #orden del sistema (numero de nodos)
n_t=800 #nodos temporales
error=.01 #error propuesto
A=np.zeros(((n,n)),float)
T=70*np.ones(((n)),float)
B=np.zeros(((n)),float)
Ba=np.zeros(((n)),float)
Xc=np.zeros(((n)),float)
E=np.zeros(((n)),float)
Ln=np.zeros(((n)),float)
Lnt=np.zeros(((8)),float)

Tn0=np.ones(((8)),float)
Tnm=np.ones(((8)),float)
Tnf=np.ones(((8)),float)

Tp2=np.ones(((n)),float)
Tp4=np.ones(((n)),float)
Tp6=np.ones(((n)),float)

L=.5/12
t=8/3600
dx=L/n
dt=t/n_t
k_unam=13
k_fimee=22.9
alpha_unam=.1775
alpha_fimee=.298
Tfusion_unam=2570
Tfusion_fimee=2770
T_0h=4500
T_0c=70
h_h=1000
h_c=300
Bh_u=h_h*dx/k_unam
Bh_f=h_h*dx/k_fimee
F_u=alpha_unam*dt/(dx**2)
F_f=alpha_fimee*dx/(dx**2)
Bc_u=h_c*dx/k_unam
Bc_f=h_c*dx/k_fimee


fimee=1
unam=2
material=fimee

if material==fimee:
    alpha=alpha_fimee
    Tfusion=Tfusion_fimee
    Bh=Bh_f
    Bc=Bc_f
    k=k_fimee
    F=F_f

if material==unam:
    alpha=alpha_unam
    Tfusion=Tfusion_unam
    Bh=Bh_u
    Bc=Bc_u
    k=k_unam
    F=F_u
    
explicito=1
implicito=2
CrankN=3

case=implicito
if case==1:
   a=1
   b=0
if case==2:
   a=0
   b=1
if case==3:
   a=.5
   b=.5


for r in range(0,8):   
   Ba=T 
   Tn0[r]=T[0]
   Tnm[r]=T[50]
   Tnf[r]=T[99]
   Lnt[r]=dt*r*3600*100
   if r==2:
      Tp2=T
   if r==4:
      Tp4=T
   if r==6:
      Tp6=T
   
   for i in range(0,n):  
       aux1=i+1
       aux2=i-1
       Ln[i]=dx*i
       
          
       if i==0:
          A[i,i]=a*(1)+b*(1+2*Bh*F+2*F)
          A[i,aux1]=a*(0)-b*(2*F)
          B[i]=a*(Ba[i]*(1-2*F-2*F*Bh)+2*F*(Ba[aux1]+Bh*T_0h))+b*(Ba[i]+2*Bh*F*T_0h)
       if i==(n-1):
          A[i,i]=a*(1)+b*(1+2*Bc*F+2*F)
          A[i,aux2]=a*(0)-b*(2*F)
          B[i]=a*(Ba[i]*(1-2*F-2*F*Bc)+2*F*(Ba[aux2]+Bc*T_0c))+b*(Ba[i]+2*Bc*F*T_0c)
       if i!=0 and i!=n-1:
          A[i,i]=a*1+b*(1+2*F)
          A[i,aux1]=-a*(0)-b*(F)
          A[i,aux2]=-a*(0)-b*(F)
          B[i]=a*(F*(Ba[aux1]+Ba[aux2])+(1-2*F)*Ba[i])+b*(Ba[i])  
   X=gauss_seidel(A,B,T,Xc,E,n,error,np)


plt.figure()
plt.plot(Lnt,Tn0,label = "x=0 ", color = 'red')
plt.hold(True)
plt.plot(Lnt,Tnm,label = "x=L/2 ", color = 'blue')
plt.plot(Lnt,Tnf,label = "x=l", color = 'green')
plt.grid(color = '0.5', linestyle = '--', linewidth = 1)

plt.savefig('GraficaT-t.pdf')
plt.show()


plt.figure()

plt.hold(True)
plt.plot(Ln,Tp2,label = "t=2 s", color = 'blue')
plt.plot(Ln,Tp4,label = "t=4 s", color = 'green')
plt.plot(Ln,Tp6,label = "t=6 s", color = 'yellow')
plt.plot(Ln,T,label = "t=8 s", color = 'red')
#plt.grid(True)
plt.grid(color = '0.5', linestyle = '--', linewidth = 1)

plt.savefig('GraficaT-x.pdf')
plt.show()

Xg=np.zeros(((10*n)),float)
for k in range(1,11):
   for j in range(0,11*n):
       if j<k*n and j>=(k-1)*n:
          Xg[j]=T[j-(k-1)*n]    
          
plt.close('all')
arr = Xg.reshape((10,n))
fig = plt.figure(figsize=(7, 3))
im = plt.imshow(arr,cmap='jet')
plt.colorbar()
plt.savefig('perfil.pdf')