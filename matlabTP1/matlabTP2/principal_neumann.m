function [ EL2, EH1 ] = principal_neumann( nom_maillage, xi, validation, visualisation,verification )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% principal_dirichlet :
% Resolution du probleme de Poisson avec conditions de Neumann
%
% SYNOPSIS [ EL2, EH1 ] = principal_neumann(nom_maillage, verification, validation,visualisation)
%          
% INPUT * nom_maillage, xi, verification, validation  : nom du maillage utilise
%                                                       micro periode
%                                                       options de tests et
%                                                       affichage
%                                                       visualisation
%
% OUTPUT - [EL2, EH1] : erreur L2 et H1
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% =====================================================
%
%
% une routine pour la mise en oeuvre des EF P1 Lagrange
% pour l'equation de Laplace suivante, avec conditions de
% Neumann sur le maillage nom_maillage.msh
%
% | -Delta  A(x/xi)u = f,   dans \Omega
% |         du/dn = 0,   sur le bord
%
% =====================================================


% lecture du maillage et affichage
% ---------------------------------
% nom_maillage = 'geomCarre.msh';
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes]=lecture_msh(nom_maillage);

% ----------------------
% calcul des matrices EF
% ----------------------

% declarations
% ------------
KK = sparse(Nbpt,Nbpt); % matrice de rigidite
MM = sparse(Nbpt,Nbpt); % matrice de rigidite
LL = zeros(Nbpt,1);     % vecteur second membre

% boucle sur les triangles
% ------------------------
for l=1:Nbtri
  % Coordonnees des sommets du triangles
%   S1=zeros(1,2);
%   S2=zeros(1,2);
%   S3=zeros(1,2);
  
  S1=Coorneu(Numtri(l,1),:);
  S2=Coorneu(Numtri(l,2),:);
  S3=Coorneu(Numtri(l,3),:);
  % calcul des matrices elementaires du triangle l 
  
   Kel=matKvar_elem(S1, S2, S3, xi);
           
   Mel=matM_elem(S1, S2, S3);
    
  % On fait l'assemmblage de la matrice globale et du second membre
  
  for i=1:3
     I=Numtri(l,i);
        for j=1:3
            J=Numtri(l,j);
            MM(I,J)=MM(I,J)+Mel(i,j);
            KK(I,J)=KK(I,J)+Kel(i,j);
        end
     
  end
end % for l

% Verification des matrices
% -------------------------
% verification = 'non';
if strcmp(verification,'oui')
    % Symetrique
    Test1=MM-MM';
    Test2=KK-KK';
    test1=max(max(abs(Test1)))
    test2=max(max(abs(Test2)))
    
    % Definie positivite de MM
    test=eig(MM);
    test=min(test)
    % Noyau de KK
    testb=ones(Nbpt,1);
    testb=testb'*KK*testb
end


% Calcul du second membre L
% -------------------------
	% utiliser la routine f.m

FF = f_neumann(Coorneu(:,1),Coorneu(:,2));
LL = MM*FF;


% inversion
% ----------
UU = KK\LL;

%visualisation
%-------------
if strcmp(visualisation,'oui')
affiche(UU, Numtri, Coorneu, sprintf('Neumann - %s', nom_maillage));
end

% validation = 'oui';
% validation
% ----------
if strcmp(validation,'oui')
UU_exact = cos(pi*Coorneu(:,1)).*cos(pi*Coorneu(:,2));
if strcmp(visualisation,'oui')
affiche(UU_exact, Numtri, Coorneu, sprintf('Neumann - %s', nom_maillage));
end
% Calcul de l erreur
NL2=UU_exact'*MM*UU_exact;
NH1=NL2+UU_exact'*KK*UU_exact;

EL2=(UU-UU_exact)'*MM*(UU-UU_exact);
EH1=sqrt(EL2+(UU-UU_exact)'*KK*(UU-UU_exact))/NH1;
EL2=sqrt(EL2)/NL2;
% attention de bien changer le terme source (dans FF)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

