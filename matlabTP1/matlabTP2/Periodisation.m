function [ UU_tilde ] = Periodisation( UU, Refneu,Corres )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Periodise la solution :
% Effectue la pseudo-elimination
%
% SYNOPSIS [ tilde_UU ] = elimine( UU, Refneu )
%          
% INPUT * UU, Refneu : le vecteur a periodiser et les références aux 
% noeuds sur le bord
%
% OUTPUT - [ tilde_UU ]: le vecteur après traitement
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

UU_tilde=UU;
Index=find((Refneu>0)==1);
for i=1:length(Index)
    if Corres(Index(i))>0
        UU_tilde(Index(i))=UU(Corres(Index(i)));
    end
end
