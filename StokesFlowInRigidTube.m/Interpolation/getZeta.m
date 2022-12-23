function [zeta_1, zeta_2] = getZeta(ind)
%%
% Intrinsic coordinates of the element nodes

switch(ind)
    case(1)
        zeta_1 = -1; zeta_2 = -1;
    case(2)
        zeta_1 = 1; zeta_2 = -1;
    case(3)
        zeta_1 = 1; zeta_2 = 1;
    case(4)
        zeta_1 = -1; zeta_2 = 1;
    case(5)
        zeta_1 = 0; zeta_2 = -1;
    case(6)
        zeta_1 = 1; zeta_2 = 0;
    case(7)
        zeta_1 = 0; zeta_2 = 1;
    case(8)
        zeta_1 = -1; zeta_2 = 0;
    case(9)
        zeta_1 = 0; zeta_2 = 0;
end
