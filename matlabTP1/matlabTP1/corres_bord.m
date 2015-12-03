function [ corres ] = corres_bord( Refneu,Coorneu )
index=find((Refneu>0)>0);
corres=zeros(1,length(Refneu));
maxi1=max(Coorneu(:,1));
mini1=min(Coorneu(:,1));
maxi2=max(Coorneu(:,2));
mini2=min(Coorneu(:,2));
list=zeros(1,length(index));
n=0;
check=1-Refneu>0;
new_Coorneu=Coorneu;
test=find(Coorneu(:,1)==maxi1);
new_Coorneu(test,1)=ones(1,length(new_Coorneu(test,1)))*mini1;
test=find(Coorneu(:,2)==maxi2);
new_Coorneu(test,2)=ones(1,length(new_Coorneu(test,2)))*mini2;
for i=1:length(index)
    if(sum(check)~=length(Refneu) && check(index(i))~=1)
        ind1=find(abs(new_Coorneu(:,1)-new_Coorneu(index(i),1))<10^-10);
        ind2=find(abs(new_Coorneu(:,2)-new_Coorneu(index(i),2))<10^-10);
        a=intersect(ind1,ind2);
        n=n+length(a);
        for j=1:length(a)
            corres(a(j))=a(1);
            check(a(j))=1;
        end
    end
end
end

