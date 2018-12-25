

L = 1
tm = 1
dt = 0.001
M = tm/dt
N = 50
X = seq(0,L,len=N+1)
T = seq(0,tm,len=M+1)
m = length(X); n=length(X);
x = matrix(rep(X,each=n),nrow=n);
y = matrix(rep(X,m),nrow=n)
x = c(x)
y = c(y)
NNODES = (N+1)^2
U = matrix(1, NNODES, 1)
NumTRI = 2*N^2
LNODES = matrix(0,NumTRI,3)

for (i in 1:N){
  for (j in 1:N){
    LNODES[i+2*(j-1)*N,1] = i+(j-1)*(N+1)
    LNODES[i+2*(j-1)*N,2] = i+j*(N+1)
    LNODES[i+2*(j-1)*N,3] = (i+1)+(j-1)*(N+1)
    LNODES[i+N+2*(j-1)*N,1] = i+1+j*(N+1)
    LNODES[i+N+2*(j-1)*N,2] = (i+1)+(j-1)*(N+1)
    LNODES[i+N+2*(j-1)*N,3] = i+j*(N+1)
  }
}

SP_Stiff <- matrix(0, NNODES, NNODES)
SP_Mass <- matrix(0, NNODES, NNODES)
LV = matrix(0, NNODES, 1)


for (n in 1:NumTRI){
  r1 = matrix(c(x[LNODES[n,1]],y[LNODES[n,1]]),nrow=2, byrow=FALSE)
  r2 = matrix(c(x[LNODES[n,2]],y[LNODES[n,2]]),nrow=2, byrow=FALSE);
  r3 = matrix(c(x[LNODES[n,3]],y[LNODES[n,3]]),nrow=2, byrow=FALSE);
  J = matrix(c(r2[1]-r1[1],r2[2]-r1[2]
               ,r3[1]-r1[1],r3[2]-r1[2]), nrow=2, byrow=TRUE) 
  Astiff = (1/(2*det(J)))*matrix(c(sum((r2-r3)*(r2-r3))
                                   ,sum((r2-r3)*(r3-r1)),sum((r2-r3)*(r1-r2))
                                   ,sum((r2-r3)*(r3-r1)),sum((r3-r1)*(r3-r1))
                                   ,sum((r3-r1)*(r1-r2)),sum((r2-r3)*(r1-r2))
                                   ,sum((r3-r1)*(r1-r2)),sum((r1-r2)*(r1-r2))), nrow=3, byrow=TRUE);
  Amass = (det(J))*matrix(c(1/12,1/24,1/24,1/24,1/12,1/24,1/24,1/24,1/12), nrow=3, byrow=TRUE);

  for (i in 1:3){
    for (j in 1:3){
      SP_Stiff[LNODES[n,i],LNODES[n,j]]=SP_Stiff[LNODES[n,i],LNODES[n,j]]+Astiff[i,j]
      SP_Mass[LNODES[n,i],LNODES[n,j]]=SP_Mass[LNODES[n,i],LNODES[n,j]]+Amass[i,j]
    }
  }
}


TMatrix = dt * SP_Stiff + SP_Mass
MatrixU <- matrix(0,length(T),length(U));

for (i in 1:NNODES){
  if (x[i]==0 | y[i]==0 | x[i]==L | y[i]==L){
    TMatrix[i,] = 0
    SP_Mass[i,] = 0
    TMatrix[i,i] = 1
  }
}

for (j in 1:M+1)
{
  RHS = SP_Mass %*% U
  U = solve(TMatrix,RHS)
  MatrixU[j,] = U
}

write.table(x,file="./coordx.txt",row.names=FALSE,col.names=FALSE)
write.table(y,file="./coordy.txt",row.names=FALSE,col.names=FALSE)
write.table(MatrixU,file="./SolutionMatrix.txt",row.names=FALSE,col.names=FALSE)
write.table(LNODES,file="./triangles.txt",row.names=FALSE,col.names=FALSE)
write.table(T,file="./timesteps.txt",row.names=FALSE,col.names=FALSE)


