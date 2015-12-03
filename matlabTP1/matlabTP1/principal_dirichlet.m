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
nom_maillage = 'geomCarre.msh';
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes]=lecture_msh(nom_maillage);

% ----------------------
% calcul des matrices EF
% ----------------------

% declarations
% ------------
KK = sparse(Nbpt,Nbpt); % matrice de rigidite
MM = sparse(Nbpt,Nbpt); % matrice de rigidite
LL = zeros(Nbpt,1);     % vecteur second membre

for l=1:Nbtri
  % Coordonnees des sommets du triangles
  % A COMPLETER
  S1=Coorneu(Numtri(l,1),:);
  S2=Coorneu(Numtri(l,2),:);
  S3=Coorneu(Numtri(l,3),:);
  % calcul des matrices elementaires du triangle l 
  
  Kel=matK_elem(S1, S2, S3);         
  Mel=matM_elem(S1, S2, S3);
    
  % On fait l'assemmblage de la matrice globale et du second membre
  for i=1:3
      I = Numtri(l,i);
      for j=1:3
          J = Numtri(l,j);
          MM(I,J) = MM(I,J) + Mel(i,j);
          KK(I,J) = KK(I,J) + Kel(i,j);
      end
  end

end % for l

% Calcul du second membre L
% -------------------------
	% A COMPLETER
	% utiliser la routine f.m
FF = f(Coorneu(:,1),Coorneu(:,2));
LL = MM*FF;
% inversion
% ----------
[tilde_MM,tilde_LL]=elimine(MM,LL,Refneu);
[tilde_KK,tilde_LL]=elimine(KK,LL,Refneu);
tilde_KK=tilde_KK-diag(Refneu>0);
%UU = (MM+KK)\LL;
tilde_AA=tilde_MM+tilde_KK;
UU=tilde_AA\tilde_LL;
% visualisation
% -------------
affiche(UU, Numtri, Coorneu, sprintf('Neumann - %s', nom_maillage));

validation = 'oui';
% validation
% ----------
if strcmp(validation,'oui')
UU_exact = sin(pi*Coorneu(:,1)).*sin(pi*Coorneu(:,2));
% Calcul de l erreur L2
ErreurL2=(UU-UU_exact)'*(tilde_MM-diag(Refneu>0))*(UU-UU_exact)/((UU_exact)'*(tilde_MM-diag(Refneu)>0)*(UU_exact))
% Calcul de l erreur H1
ErreurH1=(UU-UU_exact)'*tilde_KK*(UU-UU_exact)/((UU_exact)'*tilde_KK*(UU_exact))
% attention de bien changer le terme source (dans FF)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

