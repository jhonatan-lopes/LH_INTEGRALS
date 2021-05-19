
clc
clear all

% TEST FOR Z PATH CUT

c=1;
a=1;

A=-1;
B=-1;

mu=1;

nu=0.3;
kap=3-4*nu; % Kolosov's constant


N=10;

auxi=1:N;
s=cos(pi*auxi/(N+1));
W=(1-s.^2)./((N+1));
%s=(1+s)/2;
clear auxi

% t(k)
auxk=1:N+1;
t=cos(pi*(auxk-1/2)./(N+1));
%t=(1+t)/2;
clear auxk

[s,t]=meshgrid(s,t);

% zeta=s*c/a;
% delta=t*c/a;

zeta=c/(2*a)*(t+1);
delta=c/(2*a)*(s+1);

lh='221';

r=1;
par=1;


z1=zeta-delta;
z2=-(zeta+delta);


km=sqrt(4*r./((r+1).^2+z1.^2));
kp=sqrt(4*r./((r+1).^2+z2.^2));

[Km,Em]=ellipke(km.^2);
[Kp,Ep]=ellipke(kp.^2);


[nl,nc]=size(km);

JOLD=zeros(nl,nc);

%for k=1:niter
for i=1:nl
    for j=1:nc
        JOLD(i,j)=LH_INTEGRALS_OLD(par,r,z1(i,j),0,Km(i,j),Em(i,j),lh);
    end
end
%end

J=LH_INTEGRALS(par,r,z1,0,Km,Em,lh);

J -  JOLD

% toc
% 
% tic
% 
% for k=1:niter
%     
%     J2=LH_INTEGRALS(r,z,0,'213byr');
%     
% end
% 
% J2
% 
% toc
% 
% J2-J