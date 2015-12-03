function val = f_neumann(x,y)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f_neumann :
% Evaluation de la fonction second membre.
%
% SYNOPSIS val = f_neumann(x,y)
%          
% INPUT * x,y : les 2 coordonnees du point ou on veut evaluer la fonction.
%
% OUTPUT - val: valeur de la fonction sur ce point.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

val = (2*pi^2)*(cos(pi*x).*cos(pi*y));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
