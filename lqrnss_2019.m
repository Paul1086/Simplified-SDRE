function [x,u,K]=lqrnss(As,Bs,Fs,Qs,Rs,x0,tspan)

%Revision Date 04/20/2007
%
%This m-file calculates and plots the outputs for a Linear Quadratic Regulator Problem
%Based on provided linear state space matrices for A and B and Performance Index
%matrices for F, Q and R.  This function takes these inputs, and using the analytical 
%solution to the Riccati equation, formulates the optimal states and inputs.
%
%
%      SYNTAX:  [x,u,K]=lqrnss(A,B,F,Q,R,x0,tspan)
%
%      INPUTS (All numeric):
%      A,B         Matrices from xdot=Ax+Bu
%      F,Q,R       Performance Index Parameters; terminal cost, error and control weighting
%      x0          State variable initial condition
%      tspan       Vector containing time span [t0 tf]
%
%      OUTPUTS:
%      x           is the state variable vector
%      u           is the input vector
%      K           is the steady-state matrix inv(R)*B'*P
%
%      The system plots the x vector, and u vector in combinations of 3.  
%      The system plots the Riccati coefficients in combinations of 4.  
%

%Define variables to use in external functions

global A E F Md tf W11 W12 W21 W22 n 

%Check for correct number of inputs

if nargin<7
   error('Incorrect number of inputs specified')
   return
end   

%Convert Variables to normal symbology to prevent problems with global statement

A=As;
B=Bs;
F=Fs;
Q=Qs;
R=Rs;
plotflag=0;                 %set plotflag to a 1 to avoid plotting of data on figures

%Define secondary variables for global passing to ode-related functions and determine matrice size

[n,m]=size(A);             %Find dimensions of A
[nb,mb]=size(B);           %Find dimensions of B
[nq,mq]=size(Q);           %Find dimensions of Q
[nr,mr]=size(R);           %Find Dimensions of R
[nf,mf]=size(F);           %Find Dimensions of F
if n~=m                    %Verify A is square
   error('A must be square')
else
   [n,n]=size(A);
end   

%Data Checks for proper setup
if length(A)>rank(ctrb(A,B))                %Check for controllability
   error('System Not Controllable')
   return
end
if (n ~= nq) | (n ~= mq)            %Check that A and Q are the same size
   error('A and Q must be the same size');
   return
end
if ~(mf==1&nf==1)
   if (nq ~= nf) | (mq ~= mf)            %Check that Q and F are the same size 
      error('Q and F must be the same size');
      return
   end
end
if ~(mr==1&nr==1)
   if (mr ~= nr) | (mb ~= nr)
      error('R must be consistent with B');
      return
   end
end   
if (Q~=zeros(mq,nq))
   mq = norm(Q,1);                   % Check if Q is positive semi-definite and symmetric
   if any(eig(Q) < -eps*mq) | (norm(Q'-Q,1)/mq > eps)
	  disp('Warning: Q is not symmetric and positive semi-definite');
   end
end   
mr = norm(R,1);                   % Check if R is positive definite and symmetric
if any(eig(R) <= -eps*mr) | (norm(R'-R,1)/mr > eps)
	disp('Warning: R is not symmetric and positive definite');
end

%Define Initial Conditions for numerical solution of x states 

t0=tspan(1);
tf=tspan(2);
tspan=[tf t0];

%Define Calculated Matrices and Vectors

E=B*inv(R)*B';      %E Matrix E=B*(1/R)*B'

%Find Hamiltonian matrix needed to calculate Analytical Solution to Riccati Equation

Z=[A,-E;-Q,-A'];

%Find Eigenvectors

%Here D is the set of eigen values of hamiltonian Z in
%diagonal form

% W - each column represents the eigenvectors of Z corresponding
% eigenvalues. 

[W,D]=eig(Z);

%Find the diagonals from D and pick the negative diagonals to create a new matrix M

%with the following program I am doing m2 as [-0.8590 -4.0326]
%and index as [2 1 3 4]

j=n;
[m1,index1]=sort(real(diag(D)));
   for i=1:1:n
      m2(i)=m1(j);
      index2(i)=index1(j);
      index2(i+n)=index1(i+n);
      j=j-1;
   end
   
Md=-diag(m2);

%Rearrange W so that it corresponds to the sort of the eigenvalues
%Rearranging is necessary because we want that first two column will
%represent the eigenvectors of all the negative eigenvalues
%second two column will
%represent the eigenvectors of all the negative eigenvalues

for i=1:2*n
   w2(:,i)=W(:,index2(i));
end
W=w2;

%Define the Modal Matrix for D and Split it into Parts


%define modal matrix such that the eigen vectors corresponding to 
%negative eigenvalues are at the W11 and W21 column and the eigenvectors 
%corresponding to ositive eigenvalues are at W12 and W22 columns. 

W11=zeros(n);
W12=zeros(n);   
W21=zeros(n);
W22=zeros(n);
j=1;
  for i=1:2*n:(2*n*n-2*n+1);
      i;
    W11(j:j+n-1)=W(i:i+n-1);
    W21(j:j+n-1)=W(i+n:i+2*n-1);
    W12(j:j+n-1)=W(2*n*n+i:2*n*n+i+n-1);
    W22(j:j+n-1)=W(2*n*n+i+n:2*n*n+i+2*n-1);
    j=j+n;
  end   
  
%Define other initial conditions for calculation of P, g, x and u

t1=0.;
tx=0;                      %time array for x
tu=0;                      %time array for u
x=0.;                       %state vector

%Calculation of optimized x

   [tx,x]=ode45('lqrnssf',fliplr(tspan),x0,odeset('refine',2,'RelTol',1e-4,'AbsTol',1e-6));

%Find u vector

j=1;
us=0;                          %Initialize computational variable 
[ntx,mtx]=size(tx);
for i=1:1:mb
  for k=1:1:ntx
    Tt=-inv(W22-F*W12)*(W21-F*W11);
    P=(W21+W22*expm(-Md*(tf-tx(k)))*Tt*expm(-Md*(tf-tx(k))))*inv(W11+W12*expm(-Md*(tf-tx(k)))*Tt*expm(-Md*(tf-tx(k))));
    K=inv(R)*B'*P;
    us1=real(-K*x(k,:)');
    us(j)=us1(i);
    tu(j)=tx(k);
    j=j+1;
  end
  u(:,i)=us';
  us=0;
  j=1;
end

% Provide final steady-state K

P=W21/W11;
K=real(inv(R)*B'*P);



end

