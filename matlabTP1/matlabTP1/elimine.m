function [ tilde_AA,tilde_LL ] = elimine( AA,LL, Refneu )
tilde_LL=LL.*(1-Refneu>0);
tilde_AA=AA;
for i=1:length(Refneu)
    if (Refneu(i)>=1)
        tilde_AA(i,:)=0;
        tilde_AA(:,i)=0;
        tilde_AA(i,i)=1;
    end
end
end

