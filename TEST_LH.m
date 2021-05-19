
clc
clear all

niter=1;

rm=0.2:.2:2.2;
zm=0:.2:2;

% rm=[2 1 .5 0]

[r,z]=meshgrid(rm,zm)

tic

k=sqrt(4.*r./((r+1).^2+z.^2));
[K,E]=ellipke(k.^2);

[nl,nc]=size(k);

J=zeros(nl,nc);

for k=1:niter
for i=1:nl
    for j=1:nc
        J(i,j)=LH_INTEGRALS(1,r(i,j),z(i,j),0,K(i,j),E(i,j),'01m1');
    end
end
end

J

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