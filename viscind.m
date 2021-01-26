%% Harsha Vaddireddy
%% Graduate Research Assitant-CFD Lab OSU
%% Function for polynomial viscosity 

function [nu]  =  viscind(T,epk,sigma,M)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Ts = T/epk;
A = 1.16145;
B = 0.14874;
C = 0.52487;
D = 0.77320;
E = 2.16178;
F = 2.43787;
omgv = (A/Ts^B)+(C/exp(D*Ts))+(E/exp(Ts*F));
nu = 26.69*(M*T)^0.5/(sigma^2*omgv);
end

