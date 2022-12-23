#include "fintrf.h"

subroutine mexFunction(nlhs, plhs, nrhs, prhs)

! call: [GN, KN, K] = RegularIntegrals_GN_KN_K(chi,xie,nhat,wJ,BasisFn,mu)

implicit none

mwPointer plhs(*), prhs(*)
integer nlhs, nrhs

mwPointer mxCreateDoubleMatrix, mxGetPr
double precision mxGetScalar
mwPointer mxGetM, mxGetN

mwPointer chi, xie, nhat, wJ, BasisFn, GN, KN, K
mwPointer NG2

double precision :: mu

integer*4 :: complexflag
complexflag = 0

chi = mxGetPr(prhs(1))
xie = mxGetPr(prhs(2))
nhat = mxGetPr(prhs(3))
wJ = mxGetPr(prhs(4))
BasisFn = mxGetPr(prhs(5))
mu = mxGetScalar(prhs(6))
NG2 = mxGetN(prhs(2))

plhs(1) = mxCreateDoubleMatrix(3_8, 27_8, complexFlag)
plhs(2) = mxCreateDoubleMatrix(3_8, 27_8, complexFlag)
plhs(3) = mxCreateDoubleMatrix(3_8, 3_8, complexFlag)
GN = mxGetPr(plhs(1))
KN = mxGetPr(plhs(2))
K = mxGetPr(plhs(3))

call RegularIntegrals_GN_KN_K(%val(chi), %val(xie), %val(nhat), %val(wJ), & 
                              %val(BasisFn), mu, int(NG2), & 
                              %val(GN), %val(KN), %val(K))

end subroutine mexFunction

!------------------------------------------------------------------------
    
subroutine RegularIntegrals_GN_KN_K(chi, xie, nhat, wJ, BasisFn, & 
                                    mu, NG2, GN, KN, K)

implicit none

integer :: NG2 ! number of gauss integration points
double precision :: chi(3) ! source point
double precision :: xie(3,NG2) ! element node coordinates
double precision :: nhat(3,NG2) ! element normal vector
double precision :: wJ(1,NG2) ! element weights and Jacobian
double precision :: BasisFn(9,NG2) ! Basis functions
double precision :: mu ! viscosity

double precision :: GN(9,9), KN(9,9), K(9,1)

double precision :: xi_(3) ! gauss point in element
double precision :: r(3), rhat(3), normr, nhat_(3), wJ_, Fn_(9,1), rrhat(6)
double precision :: GG_(6,1), GN_(6,9), KK_(6,1), KN_(6,9), K_(6,1) 

double precision, parameter :: pi = 3.141592653589793d0
double precision, parameter, dimension(6) :: eye = &
                    (/1.d0,0.d0,0.d0,1.d0,0.d0,1.d0/)
                    ! vector representation of identity matrix

integer :: mm ! gauss point index

GN_ = 0.d0
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

    GG_(:,1) = (eye + rrhat)*(wJ_/(8.d0*pi*mu*normr))
    KK_(:,1) = rrhat* & 
                 ((3.d0*dot_product(rhat,nhat_)*wJ_)/(4.d0*pi*normr*normr))

    GN_ = GN_ + matmul(GG_,transpose(Fn_))
    KN_ = KN_ + matmul(KK_,transpose(Fn_)) 
    K_ = K_ + KK_
enddo

GN(1,:) = GN_(1,:)
GN(2,:) = GN_(2,:)
GN(3,:) = GN_(3,:)
GN(4,:) = GN_(2,:)
GN(5,:) = GN_(4,:)
GN(6,:) = GN_(5,:)
GN(7,:) = GN_(3,:)
GN(8,:) = GN_(5,:)
GN(9,:) = GN_(6,:)

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

end subroutine RegularIntegrals_GN_KN_K
    