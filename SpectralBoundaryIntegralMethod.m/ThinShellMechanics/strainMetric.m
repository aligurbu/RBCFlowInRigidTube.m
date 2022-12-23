function epsilbrev = strainMetric(EE, GG, FF, WW, ee, gg, ff)
%% Strain metric
epsilbrev(:,:,1) = (GG.*ee - FF.*ff - WW)./(2*WW); % epsilbrev_thet_thet
epsilbrev(:,:,2) = (EE.*ff - FF.*ee)./(2*WW);      % epsilbrev_phi_thet
epsilbrev(:,:,3) = (GG.*ff - FF.*gg)./(2*WW);      % epsilbrev_thet_phi
epsilbrev(:,:,4) = (EE.*gg - FF.*ff - WW)./(2*WW); % epsilbrev_phi_phi