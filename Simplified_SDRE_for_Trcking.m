tic
clc;
% clear all;
R1 = 2.5;


x0 = [20;120;3;1;5;3];
C=eye(6);

Q = [1000 0 0 0 0 0;
     0 150000 0 0 0 0;
     0 0 100 0 0 0;
     0 0 0 1 0 0;
     0 0 0 0 10 0;
     0 0 0 0 0 1000];

R = [0.1 0 0;
    0 1 0;
    0 0 1];


F = 0.25*eye(6)

t = 0;
q = 1;
u = [0;0;0.1734];

t1=0;
t2=0.005;
% tf=8;
tspan = [t1 t2];
delta = 0.005;


lamda = 8.1;


while (t2<=23)
    
     % v=7 || 22.68 || 136.08 || 
     % v=8 || 25.92 || 155.52 ||
     % v=9 || 29.16 || 174.96 ||

z = '[25;150;0;0;0;0]';     
      if t2<=4
          v=7;
          m1(q)=22.68;
          z = '[22.68;136.08;0;0;0;0]'
      elseif t2 <=8
          v=8;
          m1(q)=25.92;
          z = '[25.92;155.52;0;0;0;0]'
      elseif t2 <=12 
          v=9;
          m1(q)=29.16;
          z = '[29.16;174.96;0;0;0;0]'
      elseif t2 <=16 
          v=8;
          m1(q)=25.92;
          z = '[25.92;155.52;0;0;0;0]'
      elseif t2<=23
          v=7;
          m1(q)=22.68;
          z = '[22.68;136.08;0;0;0;0]'
          
      end
      
          
    
      jr = 2.88;
      rho = 1.25;
      R1 = 2.5;
      cp =0.47;
      eta = 1;
      jg = 0.22;
      p=3;
      fm = 0.4382;
      bg = 0.3;
      i=6;
      kg = 75;
      Tg(q)=p*fm*x0(6);
      
      a1 = 0.5109; a2 = 116; a3 = 0.4; a4 = 5; a5 = 21;
      a6 = 0.0068; a7 = 0.08; a8 = 0.035;


       lamda(q) = x0(1)*R1/v;
       nn = lamda(q);
%       
      
      beta1 = 0;    
      lamda_bar_1 = 1/(1/(nn+a7*beta1)-a8/(beta1^3+1));
      cp_1(q) = a1*(a2/lamda_bar_1-a3*beta1-a4)*exp(-a5/lamda_bar_1)+a6*nn;

      h_start = q;
    
       l1_lqtnss(q)=x0(1);
       l2_lqtnss(q)=x0(2);
       l3_lqtnss(q)=x0(3);
       l4_lqtnss(q)=x0(4);
       l5_lqtnss(q)=x0(5);
       l6_lqtnss(q)=x0(6);
       q;
       
       error2(q) = m1(q) - l1_lqtnss(q);
       
       l = [l1_lqtnss(q);l2_lqtnss(q);l3_lqtnss(q);l4_lqtnss(q);l5_lqtnss(q);l6_lqtnss(q)];

       z1(q) = t2;

       A11 = ((1/(2*jr))*(pi*rho*R1^2*cp*v^3)/((x0(1))^2));
       A12 = 0;
       %A13 = -2.0833;
       A13 = -i/(eta*jr);
       A14 = 0;
       A15 = 0;
       A16 = 0;

       A21 = 0;
       A22 = 0;
       %A23 = 4.5454;
       A23 = (1/jg);
       A24 = 0;
       A25 = 0;
       %A26 = -5.9754;
       A26 = -(1/jg)*p*fm;
       

       %A31 = 450+((3.60*v^3)/((x(1))^2));
       A31 = 450 + (((i*bg)/(2*jr))*((pi*rho*R1^2*cp*v^3)/((x0(1))^2)));
       A32 = -kg;
       A33 = -bg*(1/jg + (i^2)/(eta*jr));
       A34 = 0;
       A35 = 0;
       %A36 = 1.79;
       A36 = (bg/jr)*p*fm;
       
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

       B12 = 0;
       B13 = 0;
       B14 = 0;

       B22 = 0;
       B23 = 0;
       B24 = 0;

       B32 = 0;
       B33 = 0;
       B34 = 0;

       B42 = 0;
       B43 = 0;
       B44 = 10;

       B52 = -24.0616;
       B53 = 0;
       B54 = 0;

       B62 = 0;
       B63 = -24.0616;
       B64 = 0;

       A = [A11 A12 A13 A14 A15 A16;
           A21 A22 A23 A24 A25 A26;
           A31 A32 A33 A34 A35 A36;
           A41 A42 A43 A44 A45 A46;
           A51 A52 A53 A54 A55 A56;
           A61 A62 A63 A64 A65 A66];

       B = [B12 B13 B14;
            B22 B23 B24;
            B32 B33 B34;
            B42 B43 B44;
            B52 B53 B54;
            B62 B63 B64];
        
        
    
        [x,y,u,trkerr] = lqtnss(A,B,C,F,Q,R,x0,z,tspan);
        xr = x(end,:);
        ur = u(end,:);
      
        
        x1(q) = x(1);
        x2(q) = x(2);
        x3(q) = x(3);
        x4(q) = x(4);
        x5(q) = x(5);
        x6(q) = x(6);
        
        c1_lqtnss(q) = ur(1);
        c2_lqtnss(q) = ur(2);
        c3_lqtnss(q) = ur(3);

        x1(q+1) = xr(1);
        x2(q+1) = xr(2);
        x3(q+1) = xr(3);
        x4(q+1) = xr(4);
        x5(q+1) = xr(5);
        x6(q+1) = xr(6);
        
        x0 = [x1(q+1);x2(q+1);x3(q+1);x4(q+1);x5(q+1);x6(q+1)]
      
        t1 = t1 + delta;
        t2=t2+delta
        tspan = [t1 t2];
        q = q + 1

        
end

plot(z1(1:4000),l1_lqtnss(1:4000))

figure,plot(z1,l1_lqtnss,z1,m1,'r-.','LineWidth',3);
legend('actual','desired')
axis ([0 10 15 35])
figure,plot(z1,c1_lqtnss,z1,c2_lqtnss,z1,c3_lqtnss)
legend('u_q','u_q','\beta_d')
axis ([0 20 -400 700])
figure,plot(z1,error2)



figure,plot(z1,lamda)
legend('lamda')
figure,plot(z1,Tg)
legend('Gen Torque')
figure,plot(z1,cp_1)
% legend('power coefficient')

figure,plot(z1,state2,z1,m2,'r-.','LineWidth',3);
legend('actual','desired')
figure,plot(z1,state4,z1,m4,'r-.','LineWidth',3);
legend('actual','desired')
figure,plot(z1(1:4000),error2(1:4000))
legend('error signal')
figure,plot(z1,c1,'g',z1,c2,'b',z1,c3,'y','LineWidth',3);
legend('ud','uq','\beta_d','0')
