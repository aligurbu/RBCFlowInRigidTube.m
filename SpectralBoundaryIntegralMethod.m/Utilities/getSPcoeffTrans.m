function P = getSPcoeffTrans(N)

% Given SH coefficients from Spherepack arranged in vector format as
% [a00,a01,a11,a02,a12,a22,...,a0N,a1N,...aNN,b11,b12,b22,...,b1N,...,bNN]'
% get the matrix to transform it to
% [a00*sqrt(pi/2), ...
%  -b11*sqrt(pi),a01*sqrt(pi/2),a11*sqrt(pi), ...
%  -b22*sqrt(pi),-b12*sqrt(pi),a02*sqrt(pi/2),a12*sqrt(pi),a22*sqrt(pi),...
%  -bNN*sqrt(pi),...,-b1N*sqrt(pi),a0N*sqrt(pi/2),a1N*sqrt(pi),...,aNN*sqrt(pi)]';
% so that rotation can be applied to it.
% P is constructed in sparse matrix format
% After rotation, the inverse of this transform has to be applied to get it
% back to Spherepack format.
% This is a permutation and scaling of the coefficients.

nnzP = (N+1)^2; % number of non-zero elements in P
len_a = (N+1)*(N+2)/2; % length of the "a" part of the coeff vector

IP = 1:(N+1)^2; % row indices of non-zero elements in P
JP = zeros(nnzP,1); % column indices of non-zero elements in P
PP = zeros(nnzP,1); % nonzeros elements of P

rowind = 0; % row index in JP and PP
for n = 0:N
    pos_a = n*(n+1)/2; % position in "a" part of the coeff vector
    pos_b = len_a + n*(n-1)/2; % position in "b" part of the coeff vector
                               % useful only for n > 0
    for m = -n:n
        rowind = rowind + 1;
        if (m<0)
            JP(rowind) = pos_b - m;
            PP(rowind) = -sqrt(pi);
        elseif (m==0)
            JP(rowind) = pos_a + 1;
            PP(rowind) = sqrt(pi/2);
        else % m>0
            JP(rowind) = pos_a + m + 1;
            PP(rowind) = sqrt(pi);
        end
    end
end

P = sparse(IP, JP, PP);
