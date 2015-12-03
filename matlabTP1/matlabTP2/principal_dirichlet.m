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
nom_maillage = 'geomCarre_h0pt02.msh';
[Nbpt,Nbtri,Coorneu,Refneu,Numtri,Reftri,Nbaretes,Numaretes,Refaretes]=lecture_msh(nom_maillage);

% ----------------------
% calcul des matrices EF
% ----------------------

% declarations
% ------------
xitest=[0.1,0.05,0.01];
err_l2=zeros(1,length(xitest));
err_h1=zeros(1,length(xitest));
num=4;

KK = sparse(Nbpt,Nbpt); % matrice de rigidite

KK3 = sparse(Nbpt,Nbpt); % matrice de rigidite
MM = sparse(Nbpt,Nbpt); % matrice de rigidite
LL = zeros(Nbpt,1);     % vecteur second membre

Aeff=mat_A_eff(num,'geomCell.msh');
for l=1:Nbtri
    % Coordonnees des sommets du triangles
    % A COMPLETER
    S1=Coorneu(Numtri(l,1),:);
    S2=Coorneu(Numtri(l,2),:);
    S3=Coorneu(Numtri(l,3),:);
    % calcul des matrices elementaires du triangle l
    
    Kel=matKeff_elem(S1, S2, S3,Aeff);
    Kel3=matK_elem(S1, S2, S3);
    Mel=matM_elem(S1, S2, S3);
    Kel2=matKvar_elem(S1, S2, S3,xi,num);
    % On fait l'assemmblage de la matrice globale et du second membre
    for i=1:3
        I = Numtri(l,i);
        for j=1:3
            J = Numtri(l,j);
            MM(I,J) = MM(I,J) + Mel(i,j);
            KK(I,J) = KK(I,J) + Kel(i,j);
            KK3(I,J) = KK3(I,J) + Kel3(i,j);
        end
    end
    
end % for l

% Calcul du second membre L
% -------------------------
% A COMPLETER
% utiliser la routine f.m
FF = f_exact(Coorneu(:,1),Coorneu(:,2),0,1);
LL = MM*(FF*0+1);
% inversion
% ----------
[tilde_MM,tilde_LL]=elimine(MM,LL,Refneu);
[tilde_KK,tilde_LL]=elimine(KK,LL,Refneu);
[tilde_KK3,tilde_LL]=elimine(KK3,LL,Refneu);
%tilde_KK=tilde_KK-diag(Refneu>0);
%UU = (MM+KK)\LL;
tilde_AA=tilde_KK;
UU=tilde_AA\tilde_LL;
affiche(UU, Numtri, Coorneu, sprintf('sol homogene - %s', nom_maillage));
UUtest=zeros(length(UU),length(xitest));
for p=1:length(xitest)
    KK2 = sparse(Nbpt,Nbpt); % matrice de rigidite
    xi=xitest(p);
    for l=1:Nbtri
        % Coordonnees des sommets du triangles
        % A COMPLETER
        S1=Coorneu(Numtri(l,1),:);
        S2=Coorneu(Numtri(l,2),:);
        S3=Coorneu(Numtri(l,3),:);
        % calcul des matrices elementaires du triangle l
        Kel2=matKvar_elem(S1, S2, S3,xi,num);
        % On fait l'assemmblage de la matrice globale et du second membre
        for i=1:3
            I = Numtri(l,i);
            for j=1:3
                J = Numtri(l,j);
                KK2(I,J) = KK2(I,J) + Kel2(i,j);
            end
        end
        
    end % for l
    [tilde_KK2,tilde_LL]=elimine(KK2,LL,Refneu);
    UUxi=tilde_KK2\tilde_LL;
    % visualisation
    % -------------
    affiche(UUxi, Numtri, Coorneu, sprintf('sol xi - %s', nom_maillage));
    validation = 'non';
    % validation
    % ----------
    % Calcul de l erreur L2
    L2=((UU)'*(tilde_MM)*(UU));
    ErreurL2=sqrt(((UUxi-UU)'*(tilde_MM)*(UUxi-UU))/L2);
    err_l2(p)=ErreurL2;
    % Calcul de l erreur H1
    H1=L2+((UU)'*tilde_KK3*(UU));
    ErreurH1=sqrt((ErreurL2^2+(UUxi-UU)'*tilde_KK3*(UUxi-UU))/H1);
    err_h1(p)=ErreurH1;
    UUtest(:,p)=UUxi;
    % attention de bien changer le terme source (dans FF)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                        fin de la routine
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

