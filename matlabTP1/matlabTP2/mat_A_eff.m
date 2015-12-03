function [AA] = mat_A_eff(num,nom_maillage)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mat_A :
% Evaluation de la matrice A .
%
% SYNOPSIS val = mat_A(x,y)
%          
% INPUT * x,y : les 2 coordonnees du point ou on veut evaluer la fonction.
%
% OUTPUT - val: valeur de la matrice sur ce point.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes]=lecture_msh(nom_maillage);
[UU1,UU2,EL2,EH1]=prob_cell(nom_maillage,num,'oui','non');
AA=zeros(2,2);
% boucle sur les triangles
% ------------------------
for l=1:Nbtri
    % Coordonnees des sommets du triangles
    S1=Coorneu(Numtri(l,1),:);
    S2=Coorneu(Numtri(l,2),:);
    S3=Coorneu(Numtri(l,3),:);
  
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

    J1=[1;0]+C'*(UU1(Numtri(l,1))*norm(1,:)'+UU1(Numtri(l,2))*norm(2,:)'+UU1(Numtri(l,3))*norm(3,:)');
    J2=[0;1]+C'*(UU2(Numtri(l,1))*norm(1,:)'+UU2(Numtri(l,2))*norm(2,:)'+UU2(Numtri(l,3))*norm(3,:)');

     for k=1:7
        Y=B*Mq(k,:)'+a;
        mat=mat_A(Y(1),Y(2),num);
        AA(1,1) = AA(1,1) + wq(k)*(mat*J1)'*J1*abs(det(B));
        AA(2,2) = AA(2,2) + wq(k)*(mat*J2)'*J2*abs(det(B));
        AA(1,2) = AA(1,2) + wq(k)*(mat*J1)'*J2*abs(det(B));
        AA(2,1) = AA(2,1) + wq(k)*(mat*J2)'*J1*abs(det(B));
        
     end; % k


end


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
