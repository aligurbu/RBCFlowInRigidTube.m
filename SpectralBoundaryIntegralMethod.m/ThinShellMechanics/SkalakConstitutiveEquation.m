function nbrev = SkalakConstitutiveEquation(ED, ES, epsilbrev)
%% Skalak constitutive equation
I1 = 2*(epsilbrev(:,:,1) + epsilbrev(:,:,4));
I2 = I1 + 4*(epsilbrev(:,:,1).*epsilbrev(:,:,4) - ...
             epsilbrev(:,:,3).*epsilbrev(:,:,2));
nbrev(:,:,1) = 0.5*(ED*I2.*(2*epsilbrev(:,:,4) + 1) + ...
                    ES*(I1 - 2*epsilbrev(:,:,4)));  % nbrev_thet_thet
nbrev(:,:,2) = -(ED*I2 - ES).*epsilbrev(:,:,3);     % nbrev_phi_thet
nbrev(:,:,3) = -(ED*I2 - ES).*epsilbrev(:,:,2);     % nbrev_thet_phi
nbrev(:,:,4) = 0.5*(ED*I2.*(2*epsilbrev(:,:,1) + 1) + ...
                    ES*(I1 - 2*epsilbrev(:,:,1)));  % nbrev_phi_phi