def gauss_seidel(A,B,X,Xc,E,n,error,np):
   ite=0
   Error_prom=1
   while error<Error_prom and ite<10000:
      ite=ite+1
      for i in range(0,n):
         Xc[i]=X[i] 
         suma=0
         for j in range(0,n):
            if j!=i:
               suma=suma+(A[i,j]*X[j])
         X[i]=(B[i]-suma)/A[i,i]
         E[i]=abs(X[i])-abs(Xc[i])
      #E[i]=abs(X[i])-abs(Xc[i])
      Error_prom=abs(np.sum(E,axis=0)/n) 
    
   print('Matriz A')
   print(A)
   print('Vector B')
   print(B)
   print('X')
   print(X)
   print('Xc')
   print(Xc)
   print('iteracion')
   print(ite)
   print('Error Promedio')
   print(Error_prom)
   return(X)

