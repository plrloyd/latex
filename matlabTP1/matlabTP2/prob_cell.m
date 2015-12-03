function [ UU1,UU2, EL2, EH1 ] = prob_cell( nom_maillage,num, validation, visualisation )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prob_cell :
% Resolution des problemes de cellule avec conditions periodiques
%
% SYNOPSIS [ EL2, EH1 ] = prob_cell(nom_maillage, validation,visualisation)
%          
% INPUT * nom_maillage, verification, validation  : nom du maillage utilise
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
% pour l'equation de Laplace suivante, avec conditions periodiques
% sur le maillage nom_maillage.msh
%
% | -div(A grad u ) u + u= f,   dans \Omega
% |         u periodique,   sur le bord
%
% =====================================================


% lecture du maillage et affichage
% ---------------------------------
% nom_maillage = 'geomCarre_per.msh';
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes]=lecture_msh(nom_maillage);
Corres=correspondance(Refneu,Coorneu);

% ----------------------
% calcul des matrices EF
% ----------------------

% declarations
% ------------
KK = sparse(Nbpt,Nbpt); % matrice de rigidite
MM = sparse(Nbpt,Nbpt); % matrice de rigidite
LL1 = zeros(Nbpt,1);     % vecteur second membre
LL2 = zeros(Nbpt,1);     % vecteur second membre
eta=0.00001;

% boucle sur les triangles
% ------------------------
for l=1:Nbtri
  % Coordonnees des sommets du triangles
  S1=Coorneu(Numtri(l,1),:);
  S2=Coorneu(Numtri(l,2),:);
  S3=Coorneu(Numtri(l,3),:);
  % calcul des matrices elementaires du triangle l 
  
   Kel=matKvar_elem(S1, S2, S3, 1, num);
           
   Mel=matM_elem(S1, S2, S3);
   
   [Lel1,Lel2]=f_cell(S1,S2,S3,num);
  % On fait l'assemblage de la matrice globale et du second membre
  for i=1:3
     I=Numtri(l,i);
        for j=1:3
            J=Numtri(l,j);
            MM(I,J)=MM(I,J)+Mel(i,j);
            KK(I,J)=KK(I,J)+Kel(i,j);
        end
     LL1(I)=LL1(I)-Lel1(i);
     LL2(I)=LL2(I)-Lel2(i);
  end

end % for l

% % Calcul du second membre L
% % -------------------------
% [FF1,FF2] = f_cell(Coorneu(:,1),Coorneu(:,2));
% LL1 = MM*FF1;
% LL2 = MM*FF2;
% Pseudo elimination
% ------------------
AA=eta*MM+KK;
[AA,LL1,LL2]=elimine_periodique(AA,LL1,LL2,Refneu,Corres);

% inversion
% ----------
UU1 = AA\LL1;
UU2 = AA\LL2;

% Periodisation
% -------------

UU1=Periodisation(UU1,Refneu,Corres);
UU2=Periodisation(UU2,Refneu,Corres);

% visualisation
% -------------
if strcmp(visualisation,'oui')
affiche(UU1, Numtri, Coorneu, sprintf('Cellule - 1 - %s', nom_maillage));
affiche(UU2, Numtri, Coorneu, sprintf('Cellule - 2 - %s', nom_maillage));
end
% validation = 'oui';
% validation
% ----------
if strcmp(validation,'oui')
UU_exact = zeros(length(Refneu),1);
if strcmp(visualisation,'oui')
affiche(UU_exact, Numtri, Coorneu, sprintf('Cellule - ref - %s', nom_maillage));
end
% Calcul de l erreur
NL2=UU_exact'*MM*UU_exact;
NH1=NL2+UU_exact'*KK*UU_exact;

EL2=(UU1-UU_exact)'*MM*(UU1-UU_exact);
EH1=sqrt(EL2+(UU1-UU_exact)'*KK*(UU1-UU_exact))/NH1;
EL2=sqrt(EL2)/NL2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

