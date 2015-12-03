close all
Exact='oui';
Neumann='non';
Periodique='non';
xi=0.5;
num=4;
visualisation='oui';
validation='oui';
% Tracer des courbes de convergence pour Dirichlet
if strcmp(Exact,'oui')
ErreursExact=zeros(4,2);
h=[0.1,0.05,0.02,0.01];
% [ErreursDirichlet(1,1),ErreursDirichlet(1,2)]=principal_dirichlet('geomCarre_h0pt1.msh','oui','non');
[ErreursExact(2,1),ErreursExact(2,2)]=prob_exact('geomCarre_h0pt05.msh',xi,num,validation,visualisation);
[ErreursExact(3,1),ErreursExact(3,2)]=prob_exact('geomCarre_h0pt02.msh',xi,num,validation,visualisation);
[ErreursExact(4,1),ErreursExact(4,2)]=prob_exact('geomCarre_h0pt01.msh',xi,num,validation,visualisation);

figure
plot(log(h),log(ErreursExact(:,1)),'*',log(h),log(ErreursExact(:,2)),'*',log(h),2*log(h)+1,log(h),log(h))
legend('Erreur L2','Erreur H1','ordre 2','ordre 1')
title('Evolution de l erreur pour le probleme exact')
end

% Tracer des courbes de convergence pour Neumann
if strcmp(Neumann,'oui')
ErreursNeumann=zeros(3,2);
h=[0.1,0.05,0.02,0.01];
[ErreursNeumann(1,1),ErreursNeumann(1,2)]=principal_neumann('geomCarre_h0pt1.msh',1,'oui','non','non');
[ErreursNeumann(2,1),ErreursNeumann(2,2)]=principal_neumann('geomCarre_h0pt05.msh',1,'oui','non','non');
[ErreursNeumann(3,1),ErreursNeumann(3,2)]=principal_neumann('geomCarre_h0pt02.msh',1,'oui','non','non');
[ErreursNeumann(4,1),ErreursNeumann(4,2)]=principal_neumann('geomCarre_h0pt01.msh',1,'oui','non','non');
figure
plot(log(h),log(ErreursNeumann(:,1)),'*',log(h),log(ErreursNeumann(:,2)),'*',log(h),2*log(h)+1,log(h),log(h))
legend('Erreur L2','Erreur H1','ordre 2','ordre 1')
title('Evolution de l erreur pour des conditions de Neumann')
end

% Tracer des courbes de convergence pour périodique
if strcmp(Periodique,'oui')
ErreursPeriodique=zeros(3,2);
h=[0.1,0.05,0.02,0.01];
[ErreursPeriodique(1,1),ErreursPeriodique(1,2)]=principal_periodique('geomCarre_h0pt1.msh','oui','non','non');
[ErreursPeriodique(2,1),ErreursPeriodique(2,2)]=principal_periodique('geomCarre_h0pt05.msh','oui','non','non');
[ErreursPeriodique(3,1),ErreursPeriodique(3,2)]=principal_periodique('geomCarre_h0pt02.msh','oui','non','non');
[ErreursPeriodique(4,1),ErreursPeriodique(4,2)]=principal_periodique('geomCarre_h0pt01.msh','oui','non','non');
figure
plot(log(h),log(ErreursPeriodique(:,1)),'*',log(h),log(ErreursPeriodique(:,2)),'*',log(h),2*log(h)+1,log(h),log(h))
legend('Erreur L2','Erreur H1','ordre 2','ordre 1')
title('Evolution de l erreur pour des conditions de Periodique')
end