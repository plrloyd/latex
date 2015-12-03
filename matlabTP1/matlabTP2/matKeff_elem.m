function [Kel] = matKeff_elem(S1, S2, S3, Aeff)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mat_elem :
% calcul la matrices de raideur elementaire en P1 lagrange
%
% SYNOPSIS [Kel] = mat_elem(S1, S2, S3)
%          
% INPUT * S1, S2, S3 : les 2 coordonnees des 3 sommets du triangle 
%                      (vecteurs reels 1x2)
%
% OUTPUT - Kel matrice de raideur elementaire (matrice 3x3)
%
% NOTE (1) le calcul est exacte (pas de condensation de masse)
%      (2) calcul direct a partir des formules donnees par 
%          les coordonnees barycentriques 
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


% calcul de la matrice de raideur
% -------------------------------
Kel = zeros(3,3);
for i=1:3
 for j=1:3
        Kel(i,j) = Kel(i,j) + 1/(abs(D)^2)*(0.5*Aeff*C'*norm(i,:)')'*C'*norm(j,:)'*abs(det(B));
 end; % j
end; % i
end
% Kel=1/(2*abs(D)) * norm*norm';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
