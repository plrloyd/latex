function [ tilde_AA,tilde_LL ] = elimine_periodique2( AA,LL, Refneu ,corres,boolk)
tilde_LL=LL;
tilde_AA=AA;
for i=1:length(Refneu)
    if(corres(i)~=0 && corres(i)~=i)
        tilde_AA(corres(i),:)=tilde_AA(corres(i),:)+tilde_AA(i,:);
        tilde_AA(:,corres(i))=tilde_AA(:,corres(i))+tilde_AA(:,i);
        tilde_LL(corres(i))=LL(corres(i))+LL(i);
        tilde_AA(i,:)=0;
        tilde_AA(:,i)=0;
        if(boolk)
        tilde_AA(i,i)=1;
        end
        tilde_LL(i)=0;
    end
end

end

