function GKx = WeaklySingularIntegrals_GxN_KxN_M(xi_e, chi, mu,...
                                                 zetachi1, zetachi2, ...
                                                 thet_min, thet_max, ...
                                                 rmax_fn,...
                                                 grx, grw, gtx, gtw)
%%
% xi_e = matrix of element node coord (3x9)
% chi = coordinates of target point on element
% mu = viscosity
% zetachi1, zetachi2, thet_min, thet_max = singular point and integ limits
% rmax_fn = index of the rmax function
% grx, grw, gtx, gtw = gauss integration locations (in [-1,1]) and weights
% GKxs = BEM matrices (3x54)

% Nbar = singularity-cancelled interpolation (9,1)
%%
ngr = length(grx); % number of integration points in r
ngt = length(gtx); % number of integration points in theta
eye = [1;0;0;0;1;0;0;0;1]; % vector representation of identity matrix

Jthet = (thet_max - thet_min)/2;

GKx = zeros(3*54,1);
for mm = 1:ngt
    thet = thet_min + (gtx(mm)+1)*Jthet;
    rm = rmax(thet, rmax_fn);
    Jr = rm/2;
    for nn = 1:ngr
        rr = (grx(nn)+1)*Jr;

        zeta1 = zetachi1 + rr*cos(thet);
        zeta2 = zetachi2 + rr*sin(thet);

%       Interpolation function
        N = interpolate_9nodequad(zeta1, zeta2);
        xi = xi_e*N; % gauss point in element
        r = chi - xi;
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
        [nhat, JJ] = jacobian_9nodequad(zeta1, zeta2, xi_e);

        G = (eye + Prhat)/(8*pi*mu);
        K = Prhat*(3*sum(rhat.*nhat)/(4*pi));

        Nbar = interpolate_polar_9nodequad(rr, thet, zetachi1, zetachi2);
        denom_vec = xi_e*Nbar;
        denom = norm(denom_vec);

        wt = JJ/denom*Jr*Jthet*grw(nn)*gtw(mm);
        for kk = 1:9
          GKx(9*kk-8:9*kk) = GKx(9*kk-8:9*kk) + G*(N(kk)*wt);
          GKx(9*kk+73:9*kk+81) = GKx(9*kk+73:9*kk+81) + K*(Nbar(kk)*wt/denom);
        end
    end
end

GKx = reshape(GKx,3,54);

end
% ------------------------------------------------------------------------------
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

