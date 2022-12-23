function kappabrev = curvatureMetric(EE, FF, GG, LL, MM, NN, WW, ll, mm, nn)
%% Curvature metric
kappabrev(:,:,1) = (GG.*(ll-LL) - FF.*(mm-MM))./WW; % kappabrev_thet_thet
kappabrev(:,:,2) = (EE.*(mm-MM) - FF.*(ll-LL))./WW; % kappabrev_phi_thet
kappabrev(:,:,3) = (GG.*(mm-MM) - FF.*(nn-NN))./WW; % kappabrev_thet_phi
kappabrev(:,:,4) = (EE.*(nn-NN) - FF.*(mm-MM))./WW; % kappabrev_phi_phi