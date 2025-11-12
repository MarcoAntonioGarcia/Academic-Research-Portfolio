print('\033[H\033[J')
import numpy as np
import matplotlib.pyplot as plt
from gauss_seidel import gauss_seidel
n=100 #orden del sistema (numero de nodos)
n_time=110 #orden del tiempo (numero de nodos temporales)
time=1.1 #segundos
error=.001 #error propuesto

A=np.zeros(((n,n)),float)
X=np.zeros(((n)),float)
B=np.zeros(((n)),float)
Xc=np.zeros(((n)),float)
E=np.zeros(((n)),float)
Ln=np.zeros(((n)),float)
B0=np.zeros(((n)),float)

#vectores para gracicar cada 25% de avance del tiempo de solucion 
X25=np.zeros(((n)),float)
X50=np.zeros(((n)),float)
X75=np.zeros(((n)),float)
X100=np.zeros(((n)),float)
X110=np.zeros(((n)),float)


theta=1
T_a=100
T_b=100
L=1
dx=L/n
dt=time/n_time
De=Dw=.5
a=dx/dt
ae=De/dx
aw=Dw/dx

for j in range(0,n_time):
   for i in range(0,n): 
       aux1=i+1
       aux2=i-1
       Ln[i]=aux1*dx-(dx/2)
       if i==0:
          A[i,i]=a+(theta)*(ae+2*aw)
          A[i,aux1]=-theta*ae
          B[i]=(a+(1-theta)*(-ae-2*aw))*B0[i]+(ae*(1-theta))*B0[aux1]+(2*aw*(1-theta))*T_a+(2*aw*theta)*T_a
       if i==(n-1):
          A[i,i]=a+(theta)*(2*ae+aw)
          A[i,aux2]=-theta*aw
          B[i]=(a+(1-theta)*(-2*ae-aw))*B0[i]+(ae*(1-theta))*B0[aux2]+(2*ae*(1-theta))*T_b+(2*ae*theta)*T_b
       if i!=0 and i!=n-1:
          A[i,i]=a+theta*(ae+aw)
          A[i,aux1]=-theta*ae
          A[i,aux2]=-theta*aw
          B[i]=(a+(1-theta)*(-ae-aw))*B0[i]+(ae*(1-theta))*B0[aux1]+(aw*(1-theta))*B0[aux2]
   AA=np.linalg.inv(A)  
   sol=np.dot(AA,B) 
   #X=gauss_seidel(A,B,X,Xc,E,n,error,np)
   B0=sol
   if j==25:
       X25=sol
   if j==50:
       X50=sol
   if j==75:
       X75=sol
   if j==100:
       X100=sol
   if j==109:
       X110=sol
       
plt.figure()
plt.plot(Ln,X25,label = "t=.25 s", color = 'red')
plt.hold(True)
plt.plot(Ln,X50,label = "t=.50 s", color = 'blue')
plt.plot(Ln,X75, label = "t=.75 s", color = 'green')
plt.plot(Ln,X100, label = "t=1 s", color = 'orange')
plt.plot(Ln,X110, label = "t=1.1 s", color = 'black')
plt.grid(True)
plt.legend(loc="lower left")
plt.xlabel(r"$Longitud\ [m]$", fontsize = 12, color ='black')
plt.ylabel(r"$Temperatura\ [^{\circ} C]$", fontsize = 12, color = 'black')
plt.title(r"$Distribucion\ de\ Temperatura$",fontsize = 12, color = 'black', verticalalignment = 'baseline', horizontalalignment = 'center')
plt.savefig('Grafica1.pdf')
plt.show()

espesor=int(n/10)
X25_1=np.zeros(((espesor*n)),float)
X50_1=np.zeros(((espesor*n)),float)
X75_1=np.zeros(((espesor*n)),float)
X100_1=np.zeros(((espesor*n)),float)
X110_1=np.zeros(((espesor*n)),float)
for k in range(1,espesor+1):
   for l in range(0,(espesor+1)*n):
       if l<k*n and l>=(k-1)*n:
          X25_1[l]=X25[l-(k-1)*n]   
          X50_1[l]=X50[l-(k-1)*n]
          X75_1[l]=X75[l-(k-1)*n] 
          X100_1[l]=X100[l-(k-1)*n]
          X110_1[l]=X100[l-(k-1)*n]

plt.figure(1)
arr = X25.reshape((1,n))
fig = plt.figure(figsize=(15, 3))
im1 = plt.imshow(arr, cmap='jet')
plt.colorbar(extend='both')
plt.clim(70, 100);
plt.xlabel(r"$Longitud\ [\%]$", fontsize = 12, color = 'black')
plt.title(r"$Distribucion\ de\ Temperatura\ [^{\circ} C]\ t=0.25\ s$",fontsize = 12, color = 'black', verticalalignment = 'baseline', horizontalalignment = 'center')
plt.savefig('Grafica2_a.pdf')

plt.figure(2) 
arr2 = X50.reshape((1,n))
fig = plt.figure(figsize=(15, 3))
im2 = plt.imshow(arr2, cmap='jet')
plt.colorbar(extend='both')
plt.clim(70, 100);
plt.xlabel(r"$Longitud\ [\%]$", fontsize = 12, color = 'black')
plt.title(r"$Distribucion\ de\ Temperatura\ [^{\circ} C]\ t=0.50\ s$",fontsize = 12, color = 'black', verticalalignment = 'baseline', horizontalalignment = 'center')
plt.savefig('Grafica2_b.pdf')

plt.figure(3)
arr3 = X75.reshape((1,n))
fig = plt.figure(figsize=(15, 3))
im3 = plt.imshow(arr3, cmap='jet')
plt.colorbar(extend='both')
plt.clim(70, 100);
plt.xlabel(r"$Longitud\ [\%]$", fontsize = 12, color = 'black')
plt.title(r"$Distribucion\ de\ Temperatura\ [^{\circ} C]\ t=0.75\ s$",fontsize = 12, color = 'black', verticalalignment = 'baseline', horizontalalignment = 'center')
plt.savefig('Grafica2_c.pdf')  

plt.figure(4)
arr4 = X100.reshape((1,n))
fig = plt.figure(figsize=(15, 3))
im4 = plt.imshow(arr4, cmap='jet')
plt.colorbar(extend='both')
plt.clim(70, 100);
plt.xlabel(r"$Longitud\ [\%]$", fontsize = 12, color = 'black')
plt.title(r"$Distribucion\ de\ Temperatura\ [^{\circ} C]\ t=1\ s$",fontsize = 12, color = 'black', verticalalignment = 'baseline', horizontalalignment = 'center')
plt.savefig('Grafica2_d.pdf')

plt.figure(4)
arr5 = X110.reshape((1,n))
fig = plt.figure(figsize=(15, 3))
im4 = plt.imshow(arr5, cmap='jet')
plt.colorbar(extend='both')
plt.clim(99, 100);
plt.xlabel(r"$Longitud\ [\%]$", fontsize = 12, color = 'black')
plt.title(r"$Distribucion\ de\ Temperatura\ [^{\circ} C]\ t=1.1\ s$",fontsize = 12, color = 'black', verticalalignment = 'baseline', horizontalalignment = 'center')
plt.savefig('Grafica2_e.pdf')



