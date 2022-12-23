function GxpN = NearlySingularIntegrals_GxpN_M(xi_e, x, mu, s01, s02, ...
                                               thet_min, thet_max, rmax_fn,...
                                               grx, grw, gtx, gtw)
%%
% xi_e                  matrix of element node coord (3x9)
% x                     coordinates of target point (source point)
% mu                    viscosity
% s01, s02              singular point
% thet_min, thet_max    integ limits
% rmax_fn               index of the rmax function
% grx, grw, gtx, gtw	gauss integration locations (in [-1,1]) and weights
% GxpN                  BEM matrices (3x27)

% Ntilde                singularity-cancelled interpolation (9,1)
%%
ngr = length(grx);      % number of integration points in r
ngt = length(gtx);      % number of integration points n theta
eye = [1;0;0;0;1;0;0;0;1];  % vector representation of identity matrix

Jthet = (thet_max - thet_min)/2;

GxpN = zeros(3*27,1);
for mm = 1:ngt
    thet = thet_min + (gtx(mm)+1)*Jthet;
    rm = rmax(thet, s01, s02, rmax_fn);
    Jr = rm/2;
    for nn = 1:ngr
        rr = (grx(nn)+1)*Jr;

        s1 = s01 + rr*cos(thet);
        s2 = s02 + rr*sin(thet);

        N = interpolate_9nodequad(s1, s2);  % Interpolation function
        xi = xi_e*N;                        % gauss point in element
        r = x - xi;
        normr = norm(r);
        rhat = r/normr;

        Prhat = [rhat(1)*rhat(1);
                 rhat(1)*rhat(2);
                 rhat(1)*rhat(3);
                 rhat(2)*rhat(1);
                 rhat(2)*rhat(2);
                 rhat(2)*rhat(3);
                 rhat(3)*rhat(1);
                 rhat(3)*rhat(2);
                 rhat(3)*rhat(3)];
        [~, JJ] = jacobian_9nodequad(s1, s2, xi_e);

        G = (eye + Prhat)/(8*pi*mu);

        wt = ((rr*JJ)/normr)*Jr*Jthet*grw(nn)*gtw(mm);
        for kk = 1:9
            GxpN(9*kk-8:9*kk) = GxpN(9*kk-8:9*kk) + G*(N(kk)*wt);
        end
    end
end
GxpN = reshape(GxpN,3,27);
end
% ------------------------------------------------------------------------------
function rmax = rmax(thet, s01, s02, fn_ind)

switch(fn_ind)
    case(1)
        rmax = 2/cos(thet);
    case(2)
        rmax = 2/cos(pi/2-thet);
    case(3)
        rmax = 2/cos(thet-pi/2);
    case(4)
        rmax = 2/cos(pi-thet);
    case(5)
        rmax = 2/cos(thet-pi);
    case(6)
        rmax = 2/cos(3*pi/2-thet);
    case(7)
        rmax = 2/cos(thet-3*pi/2);
    case(8)
        rmax = 2/cos(thet);
    case(9)
        rmax = (1-s01)/cos(thet);
    case(10)
        rmax = 2/cos(pi/2-thet);
    case(11)
        rmax = (s01+1)/cos(pi-thet);
    case(12)
        rmax = (1-s02)/cos(thet-pi/2);
    case(13)
        rmax = 2/cos(pi-thet);
    case(14)
        rmax = (s02+1)/cos(3*pi/2-thet);
    case(15)
        rmax = (s01+1)/cos(thet-pi);
    case(16)
        rmax = 2/cos(3*pi/2-thet);
    case(17)
        rmax = (1-s01)/cos(thet);
    case(18)
        rmax = (s02+1)/cos(pi/2+thet);
    case(19)
        rmax = 2/cos(thet);
    case(20)
        rmax = (1-s02)/cos(pi/2-thet);
    case(21)
        rmax = (1-s01)/cos(thet);
    case(22)
        rmax = (1-s02)/cos(pi/2-thet);
    case(23)
        rmax = (s01+1)/cos(pi-thet);
    case(24)
        rmax = (s02+1)/cos(3*pi/2-thet);
end
end

