function K = WeaklySingularIntegrals_KxN_M(xi_e, x, s01, s02, ...
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
% Kx                    BEM matrices (3x27)

% Nbar                  singularity-cancelled interpolation (9,1)
%%
ngr = length(grx);      % number of integration points in r
ngt = length(gtx);      % number of integration points n theta

Jthet = (thet_max - thet_min)/2;

Kx = zeros(3*27,1);
for mm = 1:ngt
    thet = thet_min + (gtx(mm)+1)*Jthet;
    rm = rmax(thet, rmax_fn);
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
        [nhat, JJ] = jacobian_9nodequad(s1, s2, xi_e);

        K = Prhat*(3*sum(rhat.*nhat)/(4*pi));

        Nbar = interpolate_polar_9nodequad(rr, thet, s01, s02);
        denom_vec = xi_e*Nbar;
        denom = norm(denom_vec);

        wt = JJ/denom*Jr*Jthet*grw(nn)*gtw(mm);
        for kk = 1:9
            Kx(9*kk-8:9*kk) = Kx(9*kk-8:9*kk) + K*(Nbar(kk)*wt/denom);
        end
    end
end
K = reshape(Kx,3,27);
end
%%
function rmax = rmax(thet, fn_ind)

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
        rmax = 2/cos(2*pi-thet);
    case(9)
        rmax = 1/cos(thet);
    case(10)
        rmax = 2/cos(pi/2-thet);
    case(11)
        rmax = 1/cos(pi-thet);
    case(12)
        rmax = 1/cos(thet-pi/2);
    case(13)
        rmax = 2/cos(pi-thet);
    case(14)
        rmax = 1/cos(3*pi/2-thet);
    case(15)
        rmax = 1/cos(thet-pi);
    case(16)
        rmax = 2/cos(3*pi/2-thet);
    case(17)
        rmax = 1/cos(thet);
    case(18)
        rmax = 1/cos(pi/2+thet);
    case(19)
        rmax = 2/cos(thet);
    case(20)
        rmax = 1/cos(pi/2-thet);
    case(21)
        rmax = 1/cos(thet);
    case(22)
        rmax = 1/cos(pi/2-thet);
    case(23)
        rmax = 1/cos(pi-thet);
    case(24)
        rmax = 1/cos(3*pi/2-thet);

end
end

