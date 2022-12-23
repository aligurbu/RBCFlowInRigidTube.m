#include "fintrf.h"

subroutine mexFunction(nlhs, plhs, nrhs, prhs)

! call: [KN, K] = RegularIntegrals_KN_K(chi,xie,nhat,wJ,BasisFn);

implicit none

mwPointer plhs(*), prhs(*)
integer nlhs, nrhs

mwPointer mxCreateDoubleMatrix, mxGetPr
double precision mxGetScalar
mwPointer mxGetM, mxGetN

mwPointer chi, xie, nhat, wJ, BasisFn, KN, K
mwPointer NG2

integer*4 :: complexflag
complexflag = 0

chi = mxGetPr(prhs(1))
xie = mxGetPr(prhs(2))
nhat = mxGetPr(prhs(3))
wJ = mxGetPr(prhs(4))
BasisFn = mxGetPr(prhs(5))
NG2 = mxGetN(prhs(2))

plhs(1) = mxCreateDoubleMatrix(3_8, 27_8, complexFlag)
plhs(2) = mxCreateDoubleMatrix(3_8, 3_8, complexFlag)
KN = mxGetPr(plhs(1))
K = mxGetPr(plhs(2))

call RegularIntegrals_KN_K(%val(chi), %val(xie), %val(nhat), %val(wJ), & 
                           %val(BasisFn), int(NG2), %val(KN), %val(K))

end subroutine mexFunction

!------------------------------------------------------------------------
    
subroutine RegularIntegrals_KN_K(chi, xie, nhat, wJ, BasisFn, NG2, KN, K)

implicit none

integer :: NG2 ! number of gauss integration points
double precision :: chi(3) ! source point
double precision :: xie(3,NG2) ! element node coordinates
double precision :: nhat(3,NG2) ! element normal vector
double precision :: wJ(1,NG2) ! element weights and Jacobian
double precision :: BasisFn(9,NG2) ! shape functions
double precision :: KN(9,9), K(9,1)

double precision :: xi_(3) ! gauss point in element
double precision :: r(3), rhat(3), normr, nhat_(3), wJ_, Fn_(9,1), rrhat(6)
double precision :: KK_(6,1), KN_(6,9),  K_(6,1)

double precision, parameter :: pi = 3.141592653589793d0

integer :: mm ! gauss point index

KN_ = 0.d0
K_ = 0.d0
do mm = 1,NG2
    xi_ = xie(1:3,mm)
    r = chi - xi_
    normr = dsqrt(dot_product(r,r))
    rhat = r/normr

    nhat_ =  nhat(1:3,mm)
    wJ_ = wJ(1,mm)
    Fn_(:,1) = BasisFn(1:9,mm)

    rrhat = (/rhat(1)*rhat(1),rhat(1)*rhat(2),rhat(1)*rhat(3), &
              rhat(2)*rhat(2),rhat(2)*rhat(3), &
              rhat(3)*rhat(3) /)

    KK_(:,1) = rrhat* &
                 ((3.d0*dot_product(rhat,nhat_)*wJ_)/(4.d0*pi*normr*normr))

    KN_ = KN_ + matmul(KK_,transpose(Fn_)) 
    K_ = K_ + KK_
enddo

KN(1,:) = KN_(1,:)
KN(2,:) = KN_(2,:)
KN(3,:) = KN_(3,:)
KN(4,:) = KN_(2,:)
KN(5,:) = KN_(4,:)
KN(6,:) = KN_(5,:)
KN(7,:) = KN_(3,:)
KN(8,:) = KN_(5,:)
KN(9,:) = KN_(6,:)

K(1,1) = K_(1,1)
K(2,1) = K_(2,1)
K(3,1) = K_(3,1)
K(4,1) = K_(2,1)
K(5,1) = K_(4,1)
K(6,1) = K_(5,1)
K(7,1) = K_(3,1)
K(8,1) = K_(5,1)
K(9,1) = K_(6,1)

end subroutine RegularIntegrals_KN_K
    