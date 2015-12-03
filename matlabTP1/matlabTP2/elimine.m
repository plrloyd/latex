function [ tilde_AA, tilde_LL ] = elimine( AA, LL , Refneu )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% elimine :
% Effectue la pseudo-elimination
%
% SYNOPSIS [ tilde_AA, tilde_LL ] = elimine( AA, LL , Refneu )
%          
% INPUT * AA, LL , Refneu : les 2 matrices � g�rer et les r�f�rences aux 
% noeuds sur le bord
%
% OUTPUT - [ tilde_AA, tilde_LL ]: les 2 matrices apr�s traitement
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tilde_AA=AA;
tilde_LL=LL;
for i=1:length(Refneu)
    if Refneu(i)~=0
        tilde_LL(i)=0;
        tilde_AA(i,:)=0;
        tilde_AA(:,i)=0;
        tilde_AA(i,i)=1;
    end

end

