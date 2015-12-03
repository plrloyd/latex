function val = f(x,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f :
% Evaluation de la fonction second membre.
%
% SYNOPSIS val = f(x,y)
%          
% INPUT * x,y : les 2 coordonnees du point ou on veut evaluer la fonction.
%
% OUTPUT - val: valeur de la fonction sur ce point.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%A COMPLETER
%val = ones(length(x),1); %test1
%val= (5*pi*pi+1)*cos(pi*x).*cos(2*pi*y); %principal neumann
%val= (4*pi*pi)*cos(pi*x).*cos(pi*y); %principal dirichlet
val= (2*pi*pi+1)*sin(pi*x+pi/2).*sin(pi*y+pi/2);%principal periodique
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
