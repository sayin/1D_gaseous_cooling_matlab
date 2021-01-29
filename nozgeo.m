%% Harsha Vaddireddy
%% Function for specif nozzle geomnetry

function [x,D,delx] = nozgeo(lm)
%filename = 'nozzle.txt';
%A = importdata(filename);
Dch = 31*10^-3;
l1 = 88.9*10^-3;
r1 = 15.5*10^-3;
thc = 30;
Dt = 15.5*10^-3;
r2 = 7.8*10^-3;
thd = 15;
l2 = l1+r1*sin(thc);
D2 = Dch-2*r1*(1-cosd(thc));
D3 = Dt+2*r2*(1-cosd(thc));
l3 = l2+((D2-D3)/(2*tand(thc)));
lt = l3+r2*sind(thc);
l5 = lt+r2*sind(thd);
D5 = Dt+2*r2*(1-cosd(thd));
l = 0:0.0004:lm;
D = zeros(numel(l),1);
delx = zeros(numel(l)-1,1);
delx(1) = 0.0001;
x = D;
x(1) = 0.017;
for i = 1:numel(l)-1
    if l(i)<l1+.000001
        D(i) = Dch;
    end
    if l(i)<l2 && l(i)>l1
        D(i) = Dch-2*(r1-(r1^2-(l(i)-l1)^2)^0.5);
    end
     if l(i)<l3 && l(i)>l2
        D(i) = D2-2*(l(i)-l2)*tand(thc);
     end
     if l(i)<l5 && l(i)>l3
        D(i) = Dt+2*(r2-(r2^2-(l(i)-lt)^2)^0.5);
     end
     if l(i)>l5
        D(i) = D5+2*(l(i)-l5)*tand(thd);
     end
     if D(i)<Dt
         D(i) = D(i-1);
     end
     if i>1
     delx(i) = ((l(i+1)-l(i))^2+(D(i-1)-D(i))^2)^0.5;
     end
    x(i+1) = x(i)+delx(i); %downstream position
end
%D(i+1) = D(i);
end

