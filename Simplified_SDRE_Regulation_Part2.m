tic;
clc

R1 = 2.5;

x0 = [2;1;3;4;4;3];
B = [1;0;0;0.5];

a7 = 6.54*10^(-7);
a6 = -4.54*10^(-7);
a5 = 1.3027*10^(-5);
a4 = -6.5416*10^(-3);
a3 = -9.7477*10^(-4);
a2 = 0.0081;
a1 = -0.0013;
a0 = 0.0061;

% a1=0.22;a2=116;a3=0.4;a4=0;a5=0;a6=5;a7=12.5;a0=0.08;c9=0.035;


%main Q and R
Q = [1000 0 0 0 0 0;
     0 3500 0 0 0 0;
     0 0 1000 0 0 0;
     0 0 0 100 0 0;
     0 0 0 0 1000 0;
     0 0 0 0 0 1000];

R = [15 5 0 0;
    5 15 0 0;
    0 0 10 0;
    0 0 0 500];
% 
% %testing R
% R = [15 5 0 0;
%     5 15 0 0;
%     0 0 10 0;
%     0 0 0 500];





%same Q and R as VSVP
% Q = [1 0 0 0 0 0;
%      0 1 0 0 0 0;
%      0 0 1 0 0 0;
%      0 0 0 1 0 0;
%      0 0 0 0 1 0;
%      0 0 0 0 0 1000];
% 
% R = [2 0 0 0;
%     0 0.25 0 0;
%     0 0 0.25 0;
%     0 0 0 0.35];





% Q = [1 0 0 0 0 0;
%      0 1 0 0 0 0;
%      0 0 1 0 0 0;
%      0 0 0 1 0 0;
%      0 0 0 0 1 0;
%      0 0 0 0 0 1];
% 
% R = [0.25 0 0 0;
%     0 0.25 0 0;
%     0 0 0.25 0;
%     0 0 0 0.25];



% F = 1*eye(6);

q = 1;

t1=0;
t2=0.001;
tf=8;
tspan = [t1 t2];
delta = 0.001;




F=[1.2535 0.5609 0.0087 0.0267 0.5657 0.0015;
0.2687 0.1038 0.0133 0.0019 0.0015 0.5587;
0 0 0 0 0 0;
15.8572 0.9702 0.5783 2.2907 0.0111 0.0008;
0 0 0 0 0 0;
0 0 0 0 0 0];





    while (t2<=4)
        
%        l  = ((x0(1))*R1)/u(1,1);
%        li_inverse = 1 / (l + 0.08*(x0(4))) - 0.035 / ((x0(4))^3 +1);
%        cp = 0.5176 * (116 * li_inverse - 0.4 * (x0(4)) - 5) * exp(-21 * li_inverse) + 0.0068 / li_inverse;

       
          q

       l1_lqrnss(q) = x0(1);
       l2_lqrnss(q) = x0(2);
       l3_lqrnss(q) = x0(3);
       l4_lqrnss(q) = x0(4);
       l5_lqrnss(q) = x0(5);
       l6_lqrnss(q) = x0(6);
       q;
       t2;
       z(q)=t2;
       l = [l1_lqrnss(q);l2_lqrnss(q);l3_lqrnss(q);l4_lqrnss(q);l5_lqrnss(q);l6_lqrnss(q)];

       
       % elements of A matrix 

       A11 = 10.65*(a5*R1*x0(4)*B(1,1)+a6*R1*((x0(4))^2)*B(1,1)+a7*R1*((x0(4))^3)*B(1,1));
       A12 = 0;
       A13 = -2.0833;
       A14 = 10.65*(a1*((B(1,1))^2) + a2*x0(4)*((B(1,1))^2) + a3*((x0(4))^2)*((B(1,1))^2));
       A15 = 0;
       A16 = 0;

       A21 = 0;
       A22 = 0;
       A23 = 4.5454;
       A24 = 0;
       A25 = 0;
       A26 = -5.9754;

       A31 = 450+19.17*(a5*R1*x0(4)*B(1,1)+a6*R1*((x0(4))^2)*B(1,1)+a7*R1*((x0(4))^3)*B(1,1));
       A32 = -75;
       A33 = -5.1136;
       A34 = 19.17*(a1*((B(1,1))^2) + a2*x0(4)*((B(1,1))^2) + a3*((x0(4))^2)*((B(1,1))^2));
       A35 = 0;
       A36 = 1.79;
       
       A41 = 0;
       A42 = 0;
       A43 = 0;
       A44 = -10;
       A45 = 0;
       A46 = 0;

       A51 = 0;
       A52 = 3*x0(6);
       A53 = 0;
       A54 = 0;
       A55 = -79.4033;
       A56 = 0;

       A61 = 0;
       A62 = -72.1848*(0.04156*x0(5)-0.4382);
       A63 = 0;
       A64 = 0;
       A65 = 0; 
       A66 = -79.4033;

     
       %n_final = n1+n2+n3+n4;
       % elements of matrix B 

       B11 = 10.65*(a0*B(1,1)+a4*x0(1)*(R1));
       B12 = 0;
       B13 = 0;
       B14 = 0;

       B21 = 0;
       B22 = 0;
       B23 = 0;
       B24 = 0;

       B31 = 19.17*(a0*B(1,1)+a4*x0(1)*(R1));
       B32 = 0;
       B33 = 0;
       B34 = 0;

       B41 = 0;
       B42 = 0;
       B43 = 0;
       B44 = 10;

       B51 = 0;
       B52 = -24.0616;
       B53 = 0;
       B54 = 0;
       
       B61 = 0;
       B62 = 0;
       B63 = -24.0616;
       B64 = 0;

       A = [A11 A12 A13 A14 A15 A16;
           A21 A22 A23 A24 A25 A26;
           A31 A32 A33 A34 A35 A36;
           A41 A42 A43 A44 A45 A46;
           A51 A52 A53 A54 A55 A56;
           A61 A62 A63 A64 A65 A66];

       B = [B11 B12 B13 B14;
            B21 B22 B23 B24;
            B31 B32 B33 B34;
            B41 B42 B43 B44;
            B51 B52 B53 B54;
            B61 B62 B63 B64];
        
  
        [x,B,K]=Simp_SDRE_Regulation_P1_lqrnss(A,B,F,Q,R,x0,tspan);

        xr = x(end,:);
        

        x1(q) = x(1);
        x2(q) = x(2);
        x3(q) = x(3);
        x4(q) = x(4);
        x5(q) = x(5);
        x6(q) = x(6);

        B;

        %control_u = -inv(R)*B'*P*x; % Commanded input vector
        c1_lqrnss(q) = B(1,1);
        c2_lqrnss(q) = B(1,2);
        c3_lqrnss(q) = B(1,3);
        c4_lqrnss(q) = B(1,4);
        c=[c1_lqrnss(q);c2_lqrnss(q);c3_lqrnss(q);c4_lqrnss(q)];
        %d_c(q) = control_u(1); % Commanded input for rudder angle
        %x0(1);
        x1(q+1) = xr(1);
        x2(q+1) = xr(2);
        x3(q+1) = xr(3);
        x4(q+1) = xr(4);
        x5(q+1) = xr(5);
        x6(q+1) = xr(6);
        
        
        x0 = [x1(q+1);x2(q+1);x3(q+1);x4(q+1);x5(q+1);x6(q+1)];
        
        t1 = t1+delta; 
        t2 = t2+delta;
        tspan = [t1 t2];
        q = q + 1;
        B;
 
  end
    



figure,plot(z(1:4000),l1_lqrnss(1:4000),'--b',z(1:4000),l2_lqrnss(1:4000),'r',z(1:4000),l3_lqrnss(1:4000),'g',z(1:4000),l4_lqrnss(1:4000),'y',z(1:4000),l5_lqrnss(1:4000),'m',z(1:4000),l6_lqrnss(1:4000),'k','LineWidth',1) 
legend('state-1','state-2','state-3','state-4','state-5','state-6','0')
figure,plot(z,l1_lqrnss,'r','LineWidth',1)
legend('state 1',0)
figure,plot(z,l2_lqrnss,'b','LineWidth',1)
legend('state 2',0)
figure,plot(z,l3_lqrnss,'r','LineWidth',1)
legend('state 3',0)
figure,plot(z,l4_lqrnss,'b','LineWidth',1)
legend('state 4',0)
figure,plot(z,l5_lqrnss,'b','LineWidth',1)
legend('state 5',0)
figure,plot(z,l6_lqrnss,'b','LineWidth',1)
legend('state 6',0)    
  

 

figure,plot(z(1:4000),c1_lqrnss(1:4000),'--b',z(1:4000),c2_lqrnss(1:4000),'r',z(1:4000),c3_lqrnss(1:4000),'g',z(1:4000),c4_lqrnss(1:4000),'y','LineWidth',2)
legend('input-1','input-2','input-3','input-4','0')
figure,plot(z(1:4000),c1_lqrnss(1:4000),'b','LineWidth',1)
legend('input-1',0)
figure,plot(z,c2_lqrnss,'b','LineWidth',1)
legend('input-1',0)
figure,plot(z,c3_lqrnss,'g','LineWidth',1)
legend('input-3',0)
figure,plot(z,c4_lqrnss,'y','LineWidth',1)
legend('input-4',0)
toc;
