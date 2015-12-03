function [ tilde_AA, tilde_LL1,tilde_LL2 ] = elimine_periodique( AA, LL1, LL2 , Refneu , Corres)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% elimine :
% Effectue la pseudo-elimination
%
% SYNOPSIS [ tilde_AA, tilde_LL ] = elimine( AA, LL , Refneu )
%          
% INPUT * AA, LL , Refneu : les 2 matrices à gérer et les références aux 
% noeuds sur le bord
%
% OUTPUT - [ tilde_AA, tilde_LL ]: les 2 matrices après traitement
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Index=find((Refneu>0)==1);
for i=1:length(Index)
    if Corres(Index(i))>0
        LL1(Corres(Index(i)))=LL1(Corres(Index(i)))+LL1(Index(i));
        LL2(Corres(Index(i)))=LL2(Corres(Index(i)))+LL2(Index(i));
        AA(Corres(Index(i)),:)=AA(Corres(Index(i)),:)+AA(Index(i),:);
        AA(:,Corres(Index(i)))=AA(:,Corres(Index(i)))+AA(:,Index(i));
        LL1(Index(i))=0;
        LL2(Index(i))=0;
        AA(Index(i),:)=0;
        AA(:,Index(i))=0;
        AA(Index(i),Index(i))=1;
    end
end




% Nombre de noeuds sur le bord hors coins
% N=0;
% for i=1:length(Refneu)
%     if Refneu(i)==1
%         % Pour les coins
%         if i>1 && i<=4 
%             LL(1)=LL(1)+LL(i);
%             AA(1,:)=AA(1,:)+AA(i,:);
%             AA(:,1)=AA(:,1)+AA(:,i);
%             LL(i)=0;
%             AA(i,:)=0;
%             AA(:,i)=0;
%             AA(i,i)=1;
%         end
%     end
%     if Refneu(i)==2
%         N=+1;
%     end
%     if Refneu(i)==3 || Refneu(i)==4
%         LL(i-N*2)=LL(i-N*2)+LL(i);
%         AA(i-N*2,:)=AA(i-N*2,:)+AA(i,:);
%         AA(:,i-N*2)=AA(:,i-N*2)+AA(:,i);
%         LL(i)=0;
%         AA(i,:)=0;
%         AA(:,i)=0;
%         AA(i,i)=1;
%     end
% end
% 
tilde_AA=AA;
tilde_LL1=LL1;
tilde_LL2=LL2;
end

