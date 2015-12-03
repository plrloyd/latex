
function [ EL2, EH1 ] = prob_homogene( nom_maillage, num, validation, visualisation )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% principal_dirichlet :
% Resolution du probleme de Poisson avec conditions de Dirichlet
%
% SYNOPSIS [ EL2, EH1 ] = principal_dirichlet(nom_maillage, verification, validation,visualisation)
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
% pour l'equation de Laplace suivante, avec conditions de
% Dirichlet sur le maillage nom_maillage.msh
%
% | -\Delta u + u= f,   dans \Omega
% |         u = 0,   sur le bord
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
AA = sparse(Nbpt,Nbpt); % matrice de rigidite
LL = zeros(Nbpt,1);     % vecteur second membre

% boucle sur les triangles
% ------------------------
for l=1:Nbtri
  % Coordonnees des sommets du triangles
  S1=Coorneu(Numtri(l,1),:);
  S2=Coorneu(Numtri(l,2),:);
  S3=Coorneu(Numtri(l,3),:);
  % calcul des matrices elementaires du triangle l 
  
   Kel=matK_elem(S1, S2, S3);
           
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

% Calcul du second membre L
% -------------------------

FF = f_dirichlet(Coorneu(:,1),Coorneu(:,2));
LL = MM*FF;

% Pseudo elimination
% ------------------

[tilde_MM,tilde_LL]=elimine(MM,LL,Refneu);
[tilde_KK,tilde_LL]=elimine(KK,LL,Refneu);
tilde_AA=tilde_MM+tilde_KK;
% inversion
% ----------
UU = tilde_AA\tilde_LL;

% visualisation
% -------------
if strcmp(visualisation,'oui')
    affiche(UU, Numtri, Coorneu, sprintf('Neumann - %s', nom_maillage));
end

% validation
% ----------
if strcmp(validation,'oui')
UU_exact = sin(pi*Coorneu(:,1)).*sin(2*pi*Coorneu(:,2));
if strcmp(visualisation,'oui')
affiche(UU_exact, Numtri, Coorneu, sprintf('Neumann - %s', nom_maillage));
end
% Calcul de l erreur
NL2=UU_exact'*tilde_MM*UU_exact;
NH1=NL2+UU_exact'*tilde_KK*UU_exact;

EL2=(UU-UU_exact)'*tilde_MM*(UU-UU_exact);
EH1=sqrt(EL2+(UU-UU_exact)'*KK*(UU-UU_exact))/NH1;
EL2=sqrt(EL2)/NL2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

