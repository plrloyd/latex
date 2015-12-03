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
nom_maillage = 'geomCarre_per.msh';
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
    % A COMPLETER
    S1=Coorneu(Numtri(l,1),:);
    S2=Coorneu(Numtri(l,2),:);
    S3=Coorneu(Numtri(l,3),:);
    % calcul des matrices elementaires du triangle l
    
    Kel=matKvar_elem(S1, S2, S3);
    %Kel=matK_elem(S1, S2, S3);
    Mel=matM_elem(S1, S2, S3);
    
    % On fait l'assemblage de la matrice globale et du second membre
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
AA=KK+MM;
% inversion
% ----------
corres=corres_bord(Refneu,Coorneu);
[tilde_MM,tilde_LL]=elimine_periodique2(MM,LL,Refneu,corres,true);
[tilde_KK,tilde_LL]=elimine_periodique2(KK,LL,Refneu,corres,false);
%[tilde_AA,tilde_LL]=elimine_periodique2(AA,LL,Refneu,corres,true);
%UU = (MM+KK)\LL;
tilde_AA=tilde_MM+tilde_KK;
%[tilde_AA,tilde_LL]=elimine_periodique(AA,LL,Refneu);
UU=tilde_AA\tilde_LL;
UU_final=UU;
UU_exact_t = cos(pi*Coorneu(:,1)).*cos(pi*Coorneu(:,2));
for i=1:length(Refneu)
    if(corres(i)~=0 && corres(i)~=i)
        UU_final(i)=UU(corres(i));
        UU_exact_t(i)=0;
    end
end
% visualisation
% -------------
affiche(UU_final, Numtri, Coorneu, sprintf('Periodique - %s', nom_maillage));

validation = 'oui';
% validation
% ----------
if strcmp(validation,'oui')
    
    affiche(UU_exact, Numtri, Coorneu, sprintf('Periodique - exact %s', nom_maillage));
    % Calcul de l erreur L2
    ErreurL2=sqrt((UU-UU_exact_t)'*(tilde_MM)*(UU-UU_exact_t)/((UU_exact_t)'*(tilde_MM)*(UU_exact_t)))
    % Calcul de l erreur H1
    ErreurH1=sqrt((UU-UU_exact_t)'*tilde_KK*(UU-UU_exact_t)/((UU_exact_t)'*tilde_KK*(UU_exact_t)))
    % attention de bien changer le terme source (dans FF)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

