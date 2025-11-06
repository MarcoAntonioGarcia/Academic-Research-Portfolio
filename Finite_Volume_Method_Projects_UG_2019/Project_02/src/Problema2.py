
def Problema2(n,dx,np,gamma,alpha,betha):
   #nodosX=nodosY
   nt=n**2 #nodos totales del sistema
   error=.001 #error propuesto
   A=np.zeros(((nt,nt)),float)
   X=np.zeros(((nt)),float)
   B=np.zeros(((nt)),float)
   Xc=np.zeros(((nt)),float)
   E=np.zeros(((nt)),float)
   Ln=np.zeros(((nt)),float)
   A2=np.zeros(((n,n)),float)     #matriz de visualizacion de malla

   rho=1
   u=v=2
   Fe=Fw=rho*u
   Fn=Fs=rho*v
   De=Dw=Ds=Dn=gamma/dx
   print(Dn)

   Cs=0
   Ce=0
   Cn=100
   Cw=100


   for i in range(0,nt): 
       #nodos de las esquinas
       if i==0:               #nodo inferior izquierda
           A[i,i]=alpha*(Fn+Fe+Fw+Fs+De+2*Dw+Dn+2*Ds)+betha*(.5*(Fn+Fe-Fs-Fw+Fs+Fw)+De+Dn+2*Ds+2*Dw)
           A[i,n]=alpha*(-Dn)+betha*(.5*Fn-Dn)
           A[i,i+1]=alpha*(-De)+betha*(.5*Fe-De)
           B[i]=alpha*((Fw+Dw)*2*Cw+(Fs+Ds)*2*Cs)+betha*(2*Cs*(.5*Fs+Ds)+2*Cw*(.5*Fw+Dw))
       if i==(n-1):           #nodo inferior derecha
           A[i,i]=alpha*(Fn+Fe+Fs+2*De+Dw+2*Ds+Dn)+betha*(.5*(Fn+Fe-Fs-Fw+Fs-Fe)+2*De+Dn+Dw+2*Ds)
           A[i,2*n-1]=alpha*(-Dn)+betha*(.5*Fn-Dn)
           A[i,i-1]=alpha*(-Fw-Dw)+betha*(-.5*Fw-Dw)
           B[i]=alpha*((De)*2*Ce+(Fs+Ds)*2*Cs)+betha*(2*Cs*(.5*Fs+Ds)+2*Ce*(-.5*Fe+De))

       if i==((n-1)*n):       #nodo superior izquierda
           A[i,i]=alpha*(Fe+Fn+Fw+De+2*Dw+2*Dn+Ds)+betha*(.5*(Fn+Fe-Fw-Fs-Fn+Fw)+De+2*Dn+Ds+2*Dw)
           A[i,i+1]=alpha*(-De)+betha*(.5*Fe-De)
           A[i,i-n]=alpha*(-Fs-Ds)+betha*(-.5*Fs-Ds)
           B[i]=alpha*((Fw+Dw)*2*Cw+(Dn)*2*Cn)+betha*(2*Cs*(-.5*Fe+Dn)+2*Cw*(.5*Fw+Dw))

       if i==(nt-1):          #nodo superior derecha
           A[i,i]=alpha*(Fe+Fn+2*De+Dw+2*Dn+Ds)+betha*(.5*(Fn+Fe-Fs-Fw-Fn-Fe)+2*De+2*Dn+Ds+Dw)
           A[i,i-1]=alpha*(-Dw-Fw)+betha*(-.5*Fw-Dw)
           A[i,i-n]=alpha*(-Ds-Fs)+betha*(-.5*Fs-Ds)
           B[i]=alpha*((Dn)*2*Cn+(De)*2*Ce)+betha*(2*Cn*(-.5*Fn+Dw)+2*Ce*(-.5*Fe+De))

        
       #nodos de las fronteras 
       if i>0 and i<n-1:        #frontera sur 
           A[i,i]=alpha*(Fn+Fe+Fs+De+Dw+Dn+2*Ds)+betha*(.5*(Fn+Fe-Fs-Fw+Fs)+De+Dn+Dw+2*Ds)
           A[i,i+1]=alpha*(-De)+betha*(.5*Fe-De)
           A[i,i+(n)]=alpha*(-Dn)+betha*(.5*Fn-Dn)
           A[i,i-1]=alpha*(-Fw-Dw)+betha*(-.5*Fw-Dw)
           B[i]=alpha*((Fs+Ds)*2*Cs)+betha*(2*Cs*(.5*Fs+Ds))

          #nodos centrales...
           for j in range(1,(n-1)):
               A[j+n*i,j+n*i]=alpha*(Fn+Fe+De+Dw+Dn+Ds)+betha*(.5*(Fn+Fe-Fs-Fw)+De+Dn+Ds+Dw)
               A[j+n*i,j+n*i+1]=alpha*(-De)+betha*(.5*Fe-De)
               A[j+n*i,j+n*i-1]=alpha*(-Fw-Dw)+betha*(-.5*Fw-Dw)
               A[j+n*i,j+n*i-n]=alpha*(-Fs-Ds)+betha*(-.5*Fs-Ds)
               A[j+n*i,j+n*i+n]=alpha*(-Dn)+betha*(.5*Fn-Dn)

       if i<n-1 and i>0:        #forntera oeste
           A[i*n,i*n]=alpha*(Fn+Fe+Fw+De+2*Dw+Dn+Ds)+betha*(.5*(Fn+Fe-Fs-Fw+Fw)+De+Dn+Ds+2*Dw)
           A[i*n,i*n+n]=alpha*(-Dn)+betha*(.5*Fn-Dw)
           A[i*n,i*n+1]=alpha*(-De)+betha*(.5*Fe-De)
           A[i*n,i*n-n]=alpha*(-Fs-Ds)+betha*(-.5*Fs-Ds)
           B[i*n]=alpha*(2*Cw*(Fw+Dw))+betha*(2*Cw*(.5*Fw+Dw))

       if i<n-1 and i>0:         #frontera este 
           A[i*n+n-1,i*n+n-1]=alpha*(Fn+Fe+2*De+Dw+Dn+Ds)+betha*(.5*(Fn+Fe-Fs-Fw-Fe)+Dn+Ds+Dw+2*De)
           A[i*n+n-1,i*n+2*n-1]=alpha*(-Dn)+betha*(.5*Fn-Dn)
           A[i*n+n-1,i*n-1]=alpha*(-Ds-Fs)+betha*(-.5*Fs-Ds)
           A[i*n+n-1,i*n+n-2]=alpha*(-Fw-Dw)+betha*(-.5*Fw-Dw)
           B[i*n+n-1]=alpha*(2*Ce*(De))+betha*(2*Ce*(-.5*Fe+De))

       if i>(nt-n) and i<nt-1:   #frontera norte
           A[i,i]=alpha*(Fe+Fn+De+Dw+2*Dn+Ds)+betha*(.5*(Fn+Fe-Fs-Fw-Fn)+Dn+Ds+Dw+2*Dn)
           A[i,i-1]=alpha*(-Fw-Dw)+betha*(-.5*Fw-Dw)
           A[i,i+1]=alpha*(-De)+betha*(.5*Fe-De)
           A[i,i-n]=alpha*(-Ds-Fs)+betha*(-.5*Fs-Ds)
           B[i]=alpha*(2*Cn*(Dn))+betha*(2*Cn*(-.5*Fn+Dn))

#X=gauss_seidel(A,B,X,Xc,E,n,error,np)
   AA=np.linalg.inv(A)  
   sol=np.dot(AA,B)    
   print(A)
   print(B)
   return(sol)


    

    
    



        
        
        
        
        
        
    





