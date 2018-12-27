 clear all;
 clc;
%  matlabpool local 16
%  atomparameters;
 Gridparameters;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
c=2.997925e8;
e=1.60217733e-19;
% me=9.1093897e-31;
Epsilon0=8.854187818e-12;
% hbar=1.05457266e-34;

Tau0=2.41888129e-17;       % 1a.u.=Tau0秒（时间）
a0=5.29177249e-11;         % 1a.u.=a0米（长度）
Eh=4.35974819e-18;         % 1a.u.=Eh焦（能量）=27.2113961ev 
% EE=5.14220824e11;          % 1a.u.=EE伏/米（电场）
laserI=0.4;       % PW/cm^2
laserlambda=800;   % nm

lambda=laserlambda*1.0e-9;
I=laserI*1.0e15;

E0=sqrt(2*I*10^4/Epsilon0/c)/Eh*e*a0;

T0=lambda/c/Tau0;
omega=2*pi/T0; 

% E0=0.06;

dt=0.02;
% n1=1;n2=3;n3=1;
N1=8;
N=12;     %N>=n1-n2-n3
T1=N1*T0;
T=N*T0;
t=0:dt:T;

eps=0.8;    %0<=eps<=1
% Ex=1.0./sqrt(1.0+eps.^2)*E0*cos(omega*t).*((t<=n1*T0).*t/(n1*T0)+(t>n1*T0).*(t<=(n1+n2)*T0)+(t>(n1+n2)*T0).*(t<(n1+n2+n3)*T0).*((n1+n2+n3)*T0-t)/(n3*T0));
% Ey=eps./sqrt(1.0+eps.^2)*E0*sin(omega*t).*((t<=n1*T0).*t/(n1*T0)+(t>n1*T0).*(t<=(n1+n2)*T0)+(t>(n1+n2)*T0).*(t<(n1+n2+n3)*T0).*((n1+n2+n3)*T0-t)/(n3*T0));
Ex=1.0./sqrt(1.0+eps.^2)*E0*(sin(pi*t/T1)).^2.*cos(omega*(t-T1/2)).*(t<T1);
Ey=eps./sqrt(1.0+eps.^2)*E0*(sin(pi*t/T1)).^2.*sin(omega*(t-T1/2)).*(t<T1);
Ax=-cumsum(Ex)*dt; 
Ay=-cumsum(Ey)*dt; 
% figure
% plot3(t,Ex,Ey)
% figure
% plot(t/T0,Ex,'r',t/T0,Ey,'b');
% figure
% plot(t/T0,Ax,'r',t/T0,Ay,'b');
% error
K=length(t);
K0=floor(T0./dt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%…………………… 空间，时间，动量 ……………

% load CO ub2
load He_2px
upx=-ubx;
load He_2py
upy=-uby;
up_plus=(upx+1i*upy)/sqrt(2);
up_minus=(upx-1i*upy)/sqrt(2);
% figure;imagesc(abs(up_plus).^2)
Ksm=length(up_minus(1,:));
indexm=((Ks-Ksm)/2+1):((Ks-Ksm)/2+Ksm);

 YY=zeros(Ks,Ks);
 u0=zeros(Ks,Ks);
 u0(indexm,indexm)=up_minus;
 u=u0;

[xx,yy]=meshgrid(x,y);
% V=-Z./sqrt(xx.^2+yy.^2+a);
V0=-(1+9*exp(-(xx.^2+yy.^2)))./sqrt(xx.^2+yy.^2+2.88172);
R_c=120;
deta_R=8;


pmx=2*pi/dx;
pmy=2*pi/dy;
px=linspace(-pmx/2,pmx/2,Ks);
py=linspace(-pmy/2,pmy/2,Ks);
[ppx,ppy]=meshgrid(px,py);
dpx=2*pi/(L0);
dpy=2*pi/(L0);



absorber = ( 1-1./(1+exp(-(sqrt(xx.^2+yy.^2)-R_c)/deta_R)) );

absorber1 = ( 1-1./(1+exp(-(sqrt(xx.^2+yy.^2)-180)/8)) );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AB=-4*1i*(dx)^2/dt;
%   YY=0;
   for j=1:K
%     tic  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    u=u.*exp(-1i*dt/2*(V0+Ex(j).*xx+Ey(j).*yy));    
    parfor jj=1:Ks
        u(jj,:)=difLap(u(jj,:),AB);
%         jj
    end    
    parfor jj=1:Ks
        u(:,jj)=difLap(u(:,jj),AB);
%         jj
    end        
    u=u.*exp(-1i*dt/2*(V0+Ex(j).*xx+Ey(j).*yy)); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    u=u.*absorber1;   
%     toc
       
       if mod(j,100)==0;
        j
       end
       
      
%       if (mod(j,floor(K0/10))==0)
         
         Uex=(1-absorber).*u;
         u=absorber.*u;
         Uexv=exp( -1i*( (Ax(j)-Ax(1))*xx+(Ay(j)-Ay(1))*yy  ) ).*Uex;                                           %%%        
         clear Uex;
        
         Uexp= fftshift( fft2(Uexv) ) *dx*dy/(2*pi);
         
         ph1=(ppx.^2/2+ppy.^2/2)*(t(K)-j*dt) + ppx*sum(Ax(j:K))*dt + sum(Ax(j:K).^2)/2*dt+ ppy*sum(Ay(j:K))*dt + sum(Ay(j:K).^2)/2*dt;    %%%
         phase=exp(-1i*ph1);
         clear ph1;
        
        Y=Uexp.*phase;
        YY=YY+Y;
        clear phase Uexp Y;
%       end
    
     
   end
   

save('2pplusATI800nm','ppx','ppy', 'YY', 'u');
% matlabpool close
% imagesc(px,py,log10(abs(YY).^2)) 
% N_plot=(Ks/2-300):(Ks/2+300);
% pxm=ppx(N_plot,N_plot);
% pym=ppy(N_plot,N_plot);
% YYp=YY(N_plot,N_plot);
% 
% 
% 
% imagesc(px(N_plot),py(N_plot),log10(abs(YYp).^2));axis xy;









