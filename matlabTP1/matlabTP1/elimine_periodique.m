function [ tilde_AA,tilde_LL ] = elimine_periodique( AA,LL, Refneu )
nb1=sum(Refneu==1);
nb2=sum(Refneu==2);
tilde_LL=LL;
tilde_AA=AA;
j=0;
m=0;
for i=2:4
tilde_AA(1,:)=tilde_AA(1,:)+tilde_AA(i,:);
tilde_AA(:,1)=tilde_AA(:,i)+tilde_AA(:,1);
tilde_LL(1)=LL(1)+LL(i);  
tilde_AA(i,:)=0;
tilde_AA(:,i)=0;
tilde_AA(i,i)=1;
tilde_LL(i)=0;
end
for i=1:length(Refneu)
    if (Refneu(i)==3 )
        j=j+1;
        ind=nb1+1-j;
        tilde_AA(ind,:)=tilde_AA(i,:)+tilde_AA(ind,:);
        tilde_AA(:,ind)=tilde_AA(:,i)+tilde_AA(:,ind);
        tilde_LL(ind)=LL(ind)+LL(i);
        tilde_AA(i,:)=0;
        tilde_AA(:,i)=0;
        tilde_AA(i,i)=1;
        tilde_LL(i)=0;
    else if (Refneu(i)==4)
            m=m+1;
            ind=nb1+nb2+1-m;
            tilde_LL(ind)=LL(ind)+LL(i);
            tilde_AA(ind,:)=tilde_AA(i,:)+tilde_AA(ind,:);
            tilde_AA(:,ind)=tilde_AA(:,i)+tilde_AA(:,ind);
            tilde_LL(i)=0;
            tilde_AA(i,:)=0;
            tilde_AA(:,i)=0;
            tilde_AA(i,i)=1;
        end
    end
end
end

