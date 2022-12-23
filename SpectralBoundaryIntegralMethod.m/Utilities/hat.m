function X = hat(x)
%% Produce a skew-symmetric matrix using the components of x vector
X = zeros(3);
X(1,2) = -x(3);
X(1,3) = x(2);
X(2,1) = x(3);
X(2,3) = -x(1);
X(3,1) = -x(2);
X(3,2) = x(1);
