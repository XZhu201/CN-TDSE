% dt=0.1; %0.04;
% dx=0.1;  dy=dx;
% 
% Ks=2^11;   
% x=(-Ks/2:Ks/2-1)*dx;  y=x;
% L0=Ks*dx;
L0=409.4;
Ks=2^11;
x=linspace(-L0/2,L0/2,Ks);

y=x;
dx=L0/(Ks-1); 
dy=dx;