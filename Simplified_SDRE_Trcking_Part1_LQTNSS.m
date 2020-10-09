function [x,y,u,trkerr]=lqtnss(As,Bs,Cs,Fs,Qs,Rs,x0,zs,tspan)

%Revision Date 06/15/2019
%Edited by Sudipta Paul, MSEE and DS Naidu
%This m-file calculates and plots the outputs for a Linear Quadratic Tracking Problem
%Based on provided linear state space matrices for A, B, and C and Performance Index
%matrices for F, Q and R.  If desired, a tracking input z may be input.  This function
%takes these inputs, and using the analytical solution to the Riccati equation, formulates
%the optimal states and inputs.
%
%
%      SYNTAX:  [x,u,y,trkerr]=lqtnss(A,B,C,F,Q,R,x0,z,tspan)
%
%      INPUTS (All numeric except z):
%      A,B,C       Matrices from xdot=Ax+Bu; y=Cx
%      F,Q,R       Performance Index Parameters; terminal cost, error and control weighting
%      x0          State variable initial condition
%      z           Tracking input, which must be a string in terms of t, otherwise input '0'
%      tspan       Vector containing time span [t0 tf]
%
%      OUTPUTS:
%      x           is the state variable vector
%      y           is the output vector
%      u           is the input vector
%      trkerr      is the tracking error, if tracking used
%
%      The system plots the diagonal Riccati coefficients, x vector, and u vector in 
%      combinations of 3.  The g vector is plotted if a tracking input is given.
%

%Check for correct number of inputs

if nargin<9
   error('Incorrect number of inputs specified')
   return
end   

%Convert Variables to normal symbology 

A=As;
B=Bs;
C=Cs;
F=Fs;
Q=Qs;
R=Rs;
z=zs;
plotflag=0;                 %set plotflag to a 1 to avoid plotting of data on figures

%Define secondary variables for global passing to ode-related functions and determine matrice size

[n,m]=size(A);             %Find dimensions of A
[nb,mb]=size(B);           %Find dimensions of B
[nc,mc]=size(C);           %Find dimensions of C
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
if length(A)>rank(obsv(A,C))                %Check for observability
   error('System Not Completely Observable')
   return
end
if (nc ~= nq) | (nc ~= mq)            %Check that C and Q have the appropriate size
   error('Number of C rows and the order of Q must be the same');
   return
end
if ~(mf==1&nf==1)
   if (nq ~= nf) | (mq ~= mf)            %Check that Q and F are the same size 
      error('Q and F must be the same size');
      return
   end
end
if ~(mr==1&nr==1)                       %Check that R and B are consistent
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

t=1;
znum=double(eval(z));
[nz,mz]=size(znum);
if length(znum)==1&znum==0
   zflag=0;                        %Flag to tell program tracking is used
else
   if mz~=1|nz~=nc               %Verify good z (and C) input
      error('Incorrect z vector or C matrix')
      return
   else
      zflag=1;                   %z input okay
   end
end
clear t

%Define Initial Conditions for numerical solution of g and x states 

t0=tspan(1);
tf=tspan(2);
t=tf;
gf=(C'*F*double(eval(z)))';
clear t
tspan=[tf t0];

%Define Calculated Matrices and Vectors

Wv=C'*Q;              %W matrix W=C'*Q
Vv=C'*Q*C;            %V matrix V=C'*Q*C
E=B*inv(R)*B';        %E Matrix E=B*(1/R)*B'
F=C'*F*C;             %Create a New F if C not the identity matrix.

%Find Hamiltonian matrix needed to calculate Analytical Solution to Riccati Equation

Z=[A,-E;-Vv,-A'];

%Find Eigenvectors

[W,D]=eig(Z);

%Find the diagonals from D and pick the negative diagonals to create a new matrix M

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

for i=1:2*n
   w2(:,i)=W(:,index2(i));
end
W=w2;

%Define the Modal Matrix for D and Split it into Parts

W11=zeros(n);
W12=zeros(n);   
W21=zeros(n);
W22=zeros(n);
j=1;
  for i=1:2*n:(2*n*n-2*n+1)
    W11(j:j+n-1)=W(i:i+n-1);
    W21(j:j+n-1)=W(i+n:i+2*n-1);
    W12(j:j+n-1)=W(2*n*n+i:2*n*n+i+n-1);
    W22(j:j+n-1)=W(2*n*n+i+n:2*n*n+i+2*n-1);
    j=j+n;
  end   
  
%Define other initial conditions for calculation of P, g, x and u

t1=0.;
tg=0.;                      %time array for g
ttrk=0.;                    %time array for error track variable trkerr
tx=0.;                      %time array for x
tu=0.;                      %time array for u
g=0.;                       %g vector
x=0.;                       %state vector

%Calculation of g 

if zflag~=0                      %Check to see if tracking used
   [tg,g]=ode15s(@lqtnssa,tspan,gf,odeset('refine',2,'RelTol',1e-4,'AbsTol',1e-6),A,E,F,Md,tf,W11,W12,W21,W22,Wv,z,n);
end   

%Calculation of optimized x

if zflag~=0                      %Check to see if tracking used
   [tx,x]=ode15s(@lqtnssb,flipud(tg),x0,odeset('refine',2,'RelTol',1e-4,'AbsTol',1e-6),tg,A,E,F,Md,tf,W11,W12,W21,W22,g,n,zflag);
   [nx,mx]=size(x);
   for i=1:1:nx
      y(i,:)=(C*x(i,:)')';
   end   
else
    y=0.;                       %output vector not used
   [tx,x]=ode15s(@lqtnssb,fliplr(tspan),x0,odeset('refine',2,'RelTol',1e-4,'AbsTol',1e-6),tg,A,E,F,Md,tf,W11,W12,W21,W22,g,n,zflag);
end

%Find Tracking array

i=1;
for ttrka=t0:.1:tf
   t=ttrka;
   trkerr(i,:)=double(eval(z))';
   ttrk(i)=ttrka;
   i=i+1;
end
clear t

%Find u vector

j=1;
us=0.;                          %Initialize computational variable 
[ntx,mtx]=size(tx);
gn=flipud(g);
for i=1:1:mb
  for k=1:1:ntx
    Tt=-inv(W22-F*W12)*(W21-F*W11);
    P=(W21+W22*expm(-Md*(tf-tx(k)))*Tt*expm(-Md*(tf-tx(k))))*inv(W11+W12*expm(-Md*(tf-tx(k)))*Tt*expm(-Md*(tf-tx(k))));
    K=inv(R)*B'*P;
    if zflag~=0                      %Check to see if tracking used
       K1=real(inv(R)*B');
       us1=real(-K*x(k,:)'+K1*gn(k,:)');
       us(j)=us1(i);
    else
       us1=real(-K*x(k,:)');
       us(j)=us1(i);
    end   
    tu(j)=tx(k);
    j=j+1;
  end
  u(:,i)=us';
  us=0;
  j=1;
end


%
%====================================================================================================
%====================================================================================================
function dg = lqtnssa(t,g,A,E,F,Md,tf,W11,W12,W21,W22,Wv,z,n)
% Function for g
%Revision Date 12/29/03
%

%Calculation of P, Riccati Analytical Solution

Tt=-inv(W22-F*W12)*(W21-F*W11);
P=(W21+W22*expm(-Md*(tf-t))*Tt*expm(-Md*(tf-t)))*inv(W11+W12*expm(-Md*(tf-t))*Tt*expm(-Md*(tf-t)));

%Other simplifications

ga=-[A-E*P]';
gb=-Wv*double(eval(z));
t=double(t);

%Definition of differential equations

dg=[ga*g+gb];

%
%====================================================================================================
%====================================================================================================
function dx = lqtnssb(t,x,tg,A,E,F,Md,tf,W11,W12,W21,W22,g,n,zflag)
% Function for x
%Revision Date 8/25/00
%

%Calculation of P, Riccati Analytical Solution
Tt=-inv(W22-F*W12)*(W21-F*W11);
P=(W21+W22*expm(-Md*(tf-t))*Tt*expm(-Md*(tf-t)))*inv(W11+W12*expm(-Md*(tf-t))*Tt*expm(-Md*(tf-t)));

%Other simplifications
if zflag~=0
   gs=[interp1(tg,g,t)];
   xb=[E*gs'];
else  
   xb=zeros(n,1);
end   

xa=[A-E*P];

%Definition of differential equations

dx=[xa*x+xb];

