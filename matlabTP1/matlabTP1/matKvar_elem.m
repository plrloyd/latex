function [Kel] = matKvar_elem(S1, S2, S3)
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

A=zeros(2,2);
A(1,1)=x2-x1;
A(2,1)=y2-y1;
A(1,2)=x3-x1;
A(2,2)=y3-y1;

x1 = 0; y1 = 0;
x2 = 1; y2 = 0;
x3 = 0; y3 = 1;
% les 3 normales a l'arete opposees (de la longueur de l'arete)
norm = zeros(3, 2);
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
w0=9/80;
w1=(155-sqrt(15))/2400;
w2=(155+sqrt(15))/2400;
w=[w0,w1,w1,w1,w2,w2,w2];
s0=1/3;
s1=(6-sqrt(15))/21;
s2=(6+sqrt(15))/21;
s3=(9+2*sqrt(15))/21;
s4=(9-2*sqrt(15))/21;
x=[s0,s1,s1,s3,s2,s2,s4];
y=[s0,s1,s3,s1,s2,s4,s2];
Kel = zeros(3,3);
detA=det(A);
for i=1:3
  for j=1:3
      for q=1:7
      xtest=A*[x(q);y(q)];
      t1=(coeff(xtest(1),xtest(2))*(inv(A)'*norm(i,:)'));
      t2=(inv(A)'*norm(j,:)');
      Kel(i,j)=Kel(i,j)+w(q)*t1'*t2*abs(detA);
      end
  end; % j
end; % i

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Coeff] = coeff(x,y)
Coeff=zeros(2,2);
Coeff(1,1)=1;
Coeff(2,2)=1;
end