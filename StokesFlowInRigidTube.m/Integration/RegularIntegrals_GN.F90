#include "fintrf.h"

subroutine mexFunction(nlhs, plhs, nrhs, prhs)

! call: GN = RegularIntegrals_GN(chi,xie,wJ,BasisFn,mu) 

implicit none

mwPointer plhs(*), prhs(*)
integer nlhs, nrhs

mwPointer mxCreateDoubleMatrix, mxGetPr
double precision mxGetScalar
mwPointer mxGetM, mxGetN

mwPointer chi, xie, wJ, BasisFn, GN
mwPointer NG2

double precision :: mu

integer*4 :: complexflag
complexflag = 0

chi = mxGetPr(prhs(1))
xie = mxGetPr(prhs(2))
wJ = mxGetPr(prhs(3))
BasisFn = mxGetPr(prhs(4))
mu = mxGetScalar(prhs(5))
NG2 = mxGetN(prhs(2))

plhs(1) = mxCreateDoubleMatrix(3_8, 27_8, complexFlag)
GN = mxGetPr(plhs(1))

call RegularIntegrals_GN(%val(chi), %val(xie), %val(wJ), %val(BasisFn), mu, &
                         int(NG2), %val(GN))

end subroutine mexFunction

!------------------------------------------------------------------------
    
subroutine RegularIntegrals_GN(chi, xie, wJ, BasisFn, mu, NG2, GN)

implicit none

double precision :: chi(3) ! source point
integer :: NG2 ! number of gauss integration points
double precision :: xie(3,NG2) ! element node coordinates
double precision :: wJ(1,NG2) ! element weights and Jacobian
double precision :: BasisFn(9,NG2) ! shape functions
double precision :: mu ! viscosity

double precision :: GN(9,9)

double precision :: xi_(3) ! gauss point in element
double precision :: r(3), rhat(3), normr, wJ_, Fn_(9,1), rrhat(6)
double precision :: GG_(6,1), GN_(6,9)

double precision, parameter :: pi = 3.141592653589793d0

integer :: mm ! gauss point index
double precision, parameter, dimension(6) :: eye = &
                    (/1.d0,0.d0,0.d0,1.d0,0.d0,1.d0/)
                    ! vector representation of identity matrix

GN_ = 0.d0
do mm = 1,NG2
    xi_ = xie(1:3,mm)
    r = chi - xi_
    normr = dsqrt(dot_product(r,r))
    rhat = r/normr

    wJ_ = wJ(1,mm)
    Fn_(:,1) = BasisFn(1:9,mm)

    rrhat = (/rhat(1)*rhat(1),rhat(1)*rhat(2),rhat(1)*rhat(3), &
              rhat(2)*rhat(2),rhat(2)*rhat(3), &
              rhat(3)*rhat(3) /)

    GG_(:,1) = (eye + rrhat)*(wJ_/(8.d0*pi*mu*normr))
    
    GN_ = GN_ + matmul(GG_,transpose(Fn_))
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

end subroutine RegularIntegrals_GN
    