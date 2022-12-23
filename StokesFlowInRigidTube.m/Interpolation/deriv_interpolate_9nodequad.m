function DN = deriv_interpolate_9nodequad(s1, s2)
%%
% Compute derivatives of intepolation functions for 9-node quadrilateral 
% element at an intrinsic point (s1, s2) \in [-1, 1] x [-1, 1]
%% 
% DN = [    dN(s1,s2)/ds1               dN(s1,s2)/ds2          ];
DN = zeros(9,2);
DN(1,1) = ((2*s1 - 1)*s2*(s2 - 1))/4; DN(1,2) = (s1*(s1 - 1)*(2*s2 - 1))/4;
DN(2,1) = ((2*s1 + 1)*s2*(s2 - 1))/4; DN(2,2) = (s1*(s1 + 1)*(2*s2 - 1))/4;
DN(3,1) = ((2*s1 + 1)*s2*(s2 + 1))/4; DN(3,2) = (s1*(s1 + 1)*(2*s2 + 1))/4;
DN(4,1) = ((2*s1 - 1)*s2*(s2 + 1))/4; DN(4,2) = (s1*(s1 - 1)*(2*s2 + 1))/4;
DN(5,1) = -s1*s2*(s2 - 1);            DN(5,2) = ((1 - s1^2)*(2*s2 - 1))/2;
DN(6,1) = ((2*s1 + 1)*(1 - s2^2))/2;  DN(6,2) = -s1*(1 + s1)*s2;
DN(7,1) = -s1*s2*(s2 + 1);            DN(7,2) = ((1 - s1^2)*(2*s2 + 1))/2;
DN(8,1) = ((2*s1 - 1)*(1 - s2^2))/2;  DN(8,2) = -s1*(s1 - 1)*s2;
DN(9,1) = -2*s1*(1 - s2^2);           DN(9,2) = -2*(1 - s1^2)*s2;