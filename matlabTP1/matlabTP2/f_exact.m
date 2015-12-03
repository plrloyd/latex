function val = f_exact(x,y,xi,num)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f_dirichlet :
% Evaluation de la fonction second membre.
%
% SYNOPSIS val = f_dirichlet(x,y)
%          
% INPUT * x,y : les 2 coordonnees du point ou on veut evaluer la fonction.
%
% OUTPUT - val: valeur de la fonction sur ce point.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch num
    case 0
        val = (2*pi^2)*(sin(pi*x).*sin(pi*y));
    case 1
        val = (3*pi^2)*(sin(pi*x).*sin(pi*y));
    case 2
        val = pi^2*(6+sin(2*pi*x/xi)).*(sin(pi*x).*sin(pi*y))-2*pi^2*cos(2*pi*x/xi).*(cos(pi*x).*sin(pi*y))/xi;
    case 3
        val = pi^2*(6+sin(2*pi*x/xi)+sin(2*pi*y/xi)).*(sin(pi*x).*sin(pi*y))-(2*pi^2/xi)*(cos(2*pi*x/xi).*cos(pi*x).*sin(pi*y)+cos(2*pi*y/xi).*cos(pi*y).*sin(pi*x));
    case 4
        val = pi^2*(2+sin(2*pi*x/xi)).*(4+sin(2*pi*y/xi)).*(sin(pi*x).*sin(pi*y))-(2*pi^2/xi)*(cos(pi*x).*sin(pi*y)+sin(pi*x).*cos(pi*y)).*((2+sin(2*pi*x/xi)).*cos(2*pi*y/xi)+(4+sin(2*pi*y/xi)).*cos(2*pi*x/xi));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
