function val = mat_A(x,y, num)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mat_A :
% Evaluation de la matrice A .
%
% SYNOPSIS val = mat_A(x,y)
%          
% INPUT * x,y : les 2 coordonnees du point ou on veut evaluer la fonction.
%
% OUTPUT - val: valeur de la matrice sur ce point.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch num
    case 0
        val = [1 0;0 1];
    case 1
        val = [1 0;0 2];
    case 2
        val = [2+sin(2*pi*x) 0;0 4];
    case 3
        val = [2+sin(2*pi*x) 0;0 4+sin(2*pi*x)];
    case 4
        val = (2+sin(2*pi*x))*(4+sin(2*pi*y))*[1 0;0 1];
    otherwise
        'FAIL'
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                     fin de la fonction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
