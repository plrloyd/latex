function [ Corres ] = correspondance( Refneu, Coorneu )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

Index=find((Refneu>0)==1);
Index1=find((Refneu==1)==1);
Index2=find((Refneu==2)==1);
Index3=find((Refneu==3)==1);
Index4=find((Refneu==4)==1);

h=(max(max(Coorneu))-min(min(Coorneu)))/(sum((Refneu==2))+2);
h=h/100;

Corres=Refneu*0;
% On cherche les coins et le bords 1
for i=1:length(Index1)
    % Coin en bas a gauche
    if -h<Coorneu(Index1(i),1) && Coorneu(Index1(i),1)<h && -h<Coorneu(Index1(i),2) && Coorneu(Index1(i),2)<h
        for j=1:length(Index1)
            % Coin en bas a droite 
            if max(Coorneu(:,1))-h<Coorneu(Index1(j),1) && Coorneu(Index1(j),1)<max(Coorneu(:,1))+h && -h<Coorneu(Index1(j),2) && Coorneu(Index1(j),2)<h
                Corres(Index1(j))=Index1(i);
            end
            % Coin en haut a gauche 
            if max(Coorneu(:,1))-h<Coorneu(Index1(j),1) && Coorneu(Index1(j),1)<max(Coorneu(:,1))+h && max(Coorneu(:,2))-h<Coorneu(Index1(j),2) && Coorneu(Index1(j),2)<max(Coorneu(:,2))+h
                Corres(Index1(j))=Index1(i);
            end
            % Coin en haut a droite 
            if -h<Coorneu(Index1(j),1) && Coorneu(Index1(j),1)<h && max(Coorneu(:,2))-h<Coorneu(Index1(j),2) && Coorneu(Index1(j),2)<max(Coorneu(:,2))+h
                Corres(Index1(j))=Index1(i);
            end
        end 
    end
    % Bords 1
    if -h<Coorneu(Index1(i),2) && Coorneu(Index1(i),2)<h && ( Coorneu(Index1(i),1)>h && (Coorneu(Index1(i),1)>max(Coorneu(:,1))+h || Coorneu(Index1(i),1)<max(Coorneu(:,1))-h)) 
        for j=1:length(Index3)
            if Coorneu(Index3(j),1)-h<Coorneu(Index1(i),1) && Coorneu(Index1(i),1)<Coorneu(Index3(j),1)+h
                Corres(Index3(j))=Index1(i);
            end
        end   
    end
end
% Bords 2
for i=1:length(Index2)
    if  Coorneu(Index2(i),2)>h && (Coorneu(Index2(i),2)>h+max(Coorneu(:,2)) || Coorneu(Index2(i),2)<max(Coorneu(:,2))-h)
        for j=1:length(Index4)
            if Coorneu(Index4(j),2)-h<Coorneu(Index2(i),2) && Coorneu(Index2(i),2)<Coorneu(Index4(j),2)+h 
                Corres(Index4(j))=Index2(i);
            end
        end   
    end
end
end

