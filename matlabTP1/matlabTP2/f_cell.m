function [Lel1,Lel2] = f_cell(S1,S2,S3,num)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f_neumann :
% Evaluation de la fonction second membre.
%
% SYNOPSIS val = f_neumann(x,y)
%          
% INPUT * x,y : les 2 coordonnees du point ou on veut evaluer la fonction.
%
% OUTPUT - Lel1, Lel2 vecteurs
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% preliminaires, pour faciliter la lecture:
x1 = S1(1); y1 = S1(2);
x2 = S2(1); y2 = S2(2);
x3 = S3(1); y3 = S3(2);

% Deformation par rapport au triangle de reference Fl(x)=Bx+a
a=zeros(2,1);
B=zeros(2,2);
a=[x1;y1];
B(:,1)=[x2;y2]-a;
B(:,2)=[x3;y3]-a;
C=inv(B);

% Points de Gauss Lobatto
s=[1/3,(6-sqrt(15))/21,(6+sqrt(15))/21,(9+2*sqrt(15))/21,(9-2*sqrt(15))/21];
Mq=[s(1),s(1);s(2),s(2);s(2),s(4);s(4),s(2);s(3),s(3);s(3),s(5);s(5),s(3)];
wq=[9/80;(155-sqrt(15))/2400;(155-sqrt(15))/2400;(155-sqrt(15))/2400;(155+sqrt(15))/2400;(155+sqrt(15))/2400;(155+sqrt(15))/2400];

% les 3 normales a l'arete opposees (de la longueur de l'arete)
norm = zeros(3, 2);
x1 = 0; y1 = 0;
x2 = 1; y2 = 0;
x3 = 0; y3 = 1;
norm(1, :) = [y2-y3, x3-x2];
norm(2, :) = [y3-y1, x1-x3];
norm(3, :) = [y1-y2, x2-x1];

% D est, au signe pres, deux fois l'aire du triangle
D = ((x2-x1)*(y3-y1) - (y2-y1)*(x3-x1));
if (abs(D) <= eps) 
  error('l aire d un triangle est nulle!!!'); 
end;

% calcul des seconds membres
% -----------------------------
Lel1=zeros(3,1);
Lel2=zeros(3,1);
for j=1:3
    for k=1:7
        Y=B*Mq(k,:)'+a;
        Lel1(j)=Lel1(j)+1/abs(D)*wq(k)*(mat_A(Y(1),Y(2),num)*[1;0])'*C'*norm(j,:)'*abs(det(B));
        Lel2(j)=Lel2(j)+1/abs(D)*wq(k)*(mat_A(Y(1),Y(2),num)*[0;1])'*C'*norm(j,:)'*abs(det(B));
    end
end; % i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
