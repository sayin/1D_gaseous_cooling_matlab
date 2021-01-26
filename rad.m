%% Harsha Vaddireddy
%% Graduate Research Assitant-CFD Lab OSU
%% Function for accounting for radaitive heat transfer

function [ dtr ] = rad(Tb,dx,M)
Tg=3080;
%cpg=1416;
D=0.05;
t=Tg/1000;
p1=(11.85+0.28*0.1*11.85);
a1=(0.1/(t^1.45))+1;
b1=0.23;
c1=1.47;
eps1=0.157*(1-((a1-1)*(1-p1)/(a1+b1-1+p1))*exp(-c1*(log10(0.225*t^2))));
p2=(11.85+(2.56*0.2*11.85)/(sqrt(t)));
a2=1.888-2.053*log10(t);
b2=1.1/(t)^1.4;
c2=0.5;
eps2=0.359*(1-((a2-1)*(1-p2)/(a2+b2-1+p2))*exp(-c2*(log10(0.225*t^2))));
de=0.00822*(log10(0.3*1185))^(2.76);
eps=eps1+eps2-de;
dQr=eps*5.67*10^-8*(Tg^4-Tb^4);
Qr=dQr*pi*D*(dx);
cpg=cp(Tg);
dtr=-(Qr/(cpg*M));
end

