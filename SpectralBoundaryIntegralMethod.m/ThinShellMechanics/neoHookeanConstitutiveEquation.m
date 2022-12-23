function nbrev = neoHookeanConstitutiveEquation(ED, ES, epsilbrev)
%% Neo-Hookean model
det_U2 = (2*epsilbrev(:,:,1)+1).*(2*epsilbrev(:,:,4)+1) - ...
          4*epsilbrev(:,:,3).*epsilbrev(:,:,2);
I2 = log(det_U2);
% nbrev(:,:,1) = 0.5*(ES + (ED*I2 - ES).*(2*epsilbrev(:,:,4)+1)./det_U2); % nbrev_thet_thet
nbrev(:,:,1) = (ES + (ED*I2 - ES).*(2*epsilbrev(:,:,4)+1)./det_U2); % nbrev_thet_thet
nbrev(:,:,2) = -(ED*I2 - ES).*epsilbrev(:,:,3)./det_U2;                 % nbrev_phi_thet
nbrev(:,:,3) = -(ED*I2 - ES).*epsilbrev(:,:,2)./det_U2;                 % nbrev_thet_phi
% nbrev(:,:,4) = 0.5*(ES + (ED*I2 - ES).*(2*epsilbrev(:,:,1)+1)./det_U2); % nbrev_phi_phi
nbrev(:,:,4) = (ES + (ED*I2 - ES).*(2*epsilbrev(:,:,1)+1)./det_U2); % nbrev_phi_phi