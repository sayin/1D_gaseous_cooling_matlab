%% Harsha Vaddireddy
%% Main Code for Gasesous cooling 
clc
clear all
Tg  =  884;   %%Temp of gas
TC  =  308;   %%Coolant of gas
mug  =  viscind(TC,809.1,2.641,2)*10^(-5); %%viscosity funtion
ug  =  182.6;   %%Velocities
uc  =  273.9;   %%Velocities
xi  =  0.017;   %%postion length injection
i  =  1;        %%Counter
[x,D,dx]  =  nozgeo(0.20);  
G  =  223.685; %%mass flow rate of gas
Mc  =  0.0278; %%mass flow rate of coolant
c  =  G*mug^(0.25)*Mc^(-1.28); %%
s  =  0.004;                   %%
xf  =  (((0.325*G^(0.8)*mug^(-0.2)*xi^(0.8)+Mc)/(0.325*G^(0.8)*mug^(-0.2)))^(-0.8))-xi;
%%Fictious leading edge

%%%%%%%%%%%%%%%%%%-----------------------------%%%%%%%%%%%%%%%%%%%%%%
X(i,1)  =  x(i,1)+xf;
Re  =  1.225*ug*X(i,1)*c/mug;
um  =  2.0437*(x(i+1,1)/s)^(-0.088)*ug;  %%Peak velcoity
delt  =  0.385*X(i,1)/(Re)^(0.2);
delm  =  0.026*x(i,1)+0.16*s;             
delh  =  0.028*x(i,1)+0.7*s;
delr  =  delt-delh;
mt  =  7*1.225*11.10*delt/8; 
p  =  delm*um/2;
q  =  (delh-delm)*(3*um+ug)/2;
r  =  (delr*(um+3*ug)/2);
k  =  mt/(p+q+r);
ma(1,1)  =  k*r;
mb(1,1)  =  k*q;
mc(1,1)  =  k*p;
if (mc(1,1)<Mc)
    Tc(1,1) = TC;  
end
if (mb(1,1)<(Mc-mc(1,1)))
    Tb(1,1) = TC;
end
if (mb(1,1)>(Mc-mc(1,1)))
    Tb(1,1) = ((Mc-mc(1,1))*TC+(mb(1,1)-(Mc-mc(1,1)))*Tg)/mb(1,1);
end
% if (ma(1,1)<(Mc-mb(1,1)-mc(1,1)))
%     Ta(1,1) = TC;
% end
% if (ma(1,1)>(Mc-mb(1,1)-mc(1,1)))
%     j = Mc-mc(1,1)-mb(1,1);
%     Ta(1,1) = (j*TC+(ma(1,1)-j)*Tg)/ma(1,1)
% end
% Ta(1,1) = Tg;
    Tc(1,1) = 308;
    Tb(1,1) = 318;
    Ta(1,1) = 580;

for i = 1:numel(x)-1
    %x(i+1) = x(i)+dx;
    X(i+1,1) = x(i+1,1)+xf;
    %mug = viscmul(Ta(i));
    mug = viscind(Tc(i),71.4,3.798,28)*10^(-5);
    Re = G*X(i+1,1)*c/mug;
    um = 2.0437*(x(i+1,1)/s)^(-0.088)*ug;
    delt = 0.385*X(i+1,1)/(Re)^(0.2);
    delm = 0.026*x(i+1,1)+0.16*s   ;         %Delta at peak velocity
    delh = 0.028*x(i+1,1)+0.7*s   ;          %Delta at half 
    delr = delt-delh;
    mt = 7*1.225*11.10*delt/8;
    p = delm*um/2;
    q = (delh-delm)*(3*um+ug)/2;
    r = (delr*(um+3*ug)/2.2);
    k = mt/(p+q+r);
    ma(i+1,1) = k*r;
    mb(i+1,1) = k*q;
    mc(i+1,1) = k*p;
    dm(i+1,1) = 0.2611*1.225*ug*Re^(-0.2)*dx(i,1);
    dTr = rad(Ta(i,1),dx(i),ma(i,1));
    Ta(i+1,1) = (dm(i+1,1)*Tg*cp(Tg)+ma(i,1)*Ta(i,1)*cp(Ta(i,1)))/(ma(i+1,1)*cp(Ta(i,1)))-dTr;
    Tb(i+1,1) = ((mb(i+1,1)-mb(i,1))*Ta(i+1,1)*cp(Ta(i+1,1))+mb(i,1)*Tb(i,1)*cp(Tb(i,1)))/(mb(i+1,1)*cp(Tb(i,1)));
    Tc(i+1,1) = ((mc(i+1,1)-mc(i,1))*Tb(i+1,1)*cp(Tb(i+1,1))+mc(i,1)*Tc(i,1)*cp(Tc(i,1)))/(mc(i+1,1)*cp(Tc(i,1)));
end
plot(x,Tc,'k');
xlabel('x (m)');
ylabel('T (K)')
hold on

filename = 'ewan102.txt';
A = importdata(filename);
plot(A(:,1),A(:,2),'r*');

% hold on
% b = 0.017:0.0004:0.10;
% y = 1730.76*b+284.717;
% plot(b,y);

