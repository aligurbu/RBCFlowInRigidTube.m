#include "fintrf.h"

subroutine mexFunction(nlhs, plhs, nrhs, prhs)

! Matlab version of this. 
! call: GKx = WeaklySingularIntegrals_GxN_KxN_M(xi_e, chi, mu,...
!                                               zetachi1, zetachi2, ...
!                                               thet_min, thet_max, ...
!                                               rmax_fn,...
!                                               grx, grw, gtx, gtw)

implicit none

mwPointer plhs(*), prhs(*)
integer nlhs, nrhs

mwPointer mxCreateDoubleMatrix, mxGetPr
double precision mxGetScalar
mwPointer mxGetM, mxGetN

mwPointer xi_e, x, grx, grw, gtx, gtw, GKx
mwPointer ngr, ngt

double precision :: mu
double precision :: s01, s02, thet_min, thet_max
integer :: rmax_fn

integer*4 :: complexflag
complexflag = 0

xi_e = mxGetPr(prhs(1))
x = mxGetPr(prhs(2))
mu = mxGetScalar(prhs(3))
s01 = mxGetScalar(prhs(4))
s02 = mxGetScalar(prhs(5))
thet_min = mxGetScalar(prhs(6))
thet_max = mxGetScalar(prhs(7))
rmax_fn = int(mxGetScalar(prhs(8)))
grx = mxGetPr(prhs(9))
grw = mxGetPr(prhs(10))
ngr = max0(mxGetM(prhs(9)),mxGetN(prhs(9)))
gtx = mxGetPr(prhs(11))
gtw = mxGetPr(prhs(12))
ngt = max0(mxGetM(prhs(11)),mxGetN(prhs(11)))

plhs(1) = mxCreateDoubleMatrix(3_8, 54_8, complexFlag)
GKx = mxGetPr(plhs(1))

call get_9NodeQuad_GKx_singular(%val(xi_e), %val(x), mu, &
                                s01, s02, thet_min, thet_max, rmax_fn, &
                                int(ngr), %val(grx), %val(grw), &
                                int(ngt), %val(gtx), %val(gtw), &
                                %val(GKx))

end subroutine mexFunction
    
!------------------------------------------------------------------------
    
subroutine get_9NodeQuad_GKx_singular(xi_e, x, mu, &
                                      s01, s02, thet_min, thet_max, rmax_fn, &
                                      ngr, grx, grw, ngt, gtx, gtw, &
                                      GKx)

implicit none

double precision :: xi_e(3,9) ! element node coordinates
double precision :: x(3) ! source point
double precision :: mu ! viscosity
double precision :: s01, s02, thet_min, thet_max ! singular point and integ limits
integer :: rmax_fn ! index of the rmax function
integer :: ngr, ngt ! number of gauss points in r and theta directions
double precision :: grx(ngr), grw(ngr), gtx(ngt), gtw(ngt) ! gauss pt locs and wts
double precision :: GKx(162) ! 3*54

double precision :: N(9), DN(9,2) ! interpolation fn and its derivative
double precision :: Nbar(9) ! singularity-cancelled interpolation

double precision :: xi(3) ! gauss point in element
double precision :: r(3), rhat(3), Prhat(9), normr, nhat(3), JJ, dxids(3,2)
double precision :: G(9), K(9) ! vector representations of 3x3 matrices

double precision, parameter :: pi = 3.141592653589793d0
double precision, parameter, dimension(9) :: eye = &
                    (/1.d0,0.d0,0.d0,0.d0,1.d0,0.d0,0.d0,0.d0,1.d0/)
                    ! vector representation of identity matrix

double precision :: thet, Jthet, rm, Jr, rr, denom_vec(3), denom
                    
integer :: mm, nn, kk ! gauss point index
double precision :: s1, s2, wt

double precision :: rmax

Jthet = (thet_max - thet_min)/2.d0

GKx = 0.d0
do mm = 1,ngt
    thet = thet_min + (gtx(mm)+1.d0)*Jthet
    rm = rmax(thet, rmax_fn)
    Jr = rm/2.d0
    do nn = 1,ngr
        rr = (grx(nn)+1.d0)*Jr
        
        s1 = s01 + rr*dcos(thet);
        s2 = s02 + rr*dsin(thet);

        ! Interpolation function
        N = (/ (s1*(s1 - 1.d0)*s2*(s2 - 1.d0))/4.d0, &
               (s1*(s1 + 1.d0)*s2*(s2 - 1.d0))/4.d0, &
               (s1*(s1 + 1.d0)*s2*(s2 + 1.d0))/4.d0, &
               (s1*(s1 - 1.d0)*s2*(s2 + 1.d0))/4.d0, &
               ((1.d0 - s1)*(1.d0 + s1)*s2*(s2 - 1.d0))/2.d0, &
               (s1*(s1 + 1.d0)*(1.d0 - s2)*(1.d0 + s2))/2.d0, &
               ((1.d0 - s1)*(1.d0 + s1)*s2*(s2 + 1.d0))/2.d0, &
               (s1*(s1 - 1.d0)*(1.d0 - s2)*(1.d0 + s2))/2.d0, &
               (1.d0 - s1)*(1.d0 + s1)*(1.d0 - s2)*(1.d0 + s2) /)
        
        xi = matmul(xi_e,N)
        
        r = x - xi
        normr = dsqrt(dot_product(r,r))
        rhat = r/normr
        
        Prhat = (/rhat(1)*rhat(1),rhat(1)*rhat(2),rhat(1)*rhat(3), &
                  rhat(2)*rhat(1),rhat(2)*rhat(2),rhat(2)*rhat(3), &
                  rhat(3)*rhat(1),rhat(3)*rhat(2),rhat(3)*rhat(3) /)
        
        DN(1:9,1) = (/ (((2.d0*s1 - 1.d0)*s2*(s2 - 1.d0))/4.d0), &
                       (((2.d0*s1 + 1.d0)*s2*(s2 - 1.d0))/4.d0), & 
                       (((2.d0*s1 + 1.d0)*s2*(s2 + 1.d0))/4.d0), &
                       (((2.d0*s1 - 1.d0)*s2*(s2 + 1.d0))/4.d0), &
                       (-s1*s2*(s2 - 1.d0)), &
                       (((2.d0*s1 + 1.d0)*(1.d0 - s2)*(1.d0 + s2))/2.d0), &
                       (-s1*s2*(s2 + 1.d0)), &
                       (((2.d0*s1 - 1.d0)*(1.d0 - s2)*(1.d0 + s2))/2.d0), &
                       (-2.d0*s1*(1.d0 - s2)*(1.d0 + s2)) /)
        DN(1:9,2) = (/ ((s1*(s1 - 1.d0)*(2.d0*s2 - 1.d0))/4.d0), &
                       ((s1*(s1 + 1.d0)*(2.d0*s2 - 1.d0))/4.d0), &
                       ((s1*(s1 + 1.d0)*(2.d0*s2 + 1.d0))/4.d0), &
                       ((s1*(s1 - 1.d0)*(2.d0*s2 + 1.d0))/4.d0), &
                       (((1.d0 - s1)*(1.d0 + s1)*(2.d0*s2 - 1.d0))/2.d0), &
                       (-s1*(1.d0 + s1)*s2), &
                       (((1.d0 - s1)*(1.d0 + s1)*(2.d0*s2 + 1.d0))/2.d0), &
                       (-s1*(s1 - 1.d0)*s2), &
                       (-2.d0*(1.d0 - s1)*(1.d0 + s1)*s2) /)
        
        dxids = matmul(xi_e,DN)
        
        nhat = (/ (dxids(2,1)*dxids(3,2)-dxids(3,1)*dxids(2,2)), &
                  (dxids(3,1)*dxids(1,2)-dxids(1,1)*dxids(3,2)), &
                  (dxids(1,1)*dxids(2,2)-dxids(2,1)*dxids(1,2)) /)
        JJ = dsqrt(dot_product(nhat,nhat))
        nhat = nhat/JJ;
        
        G = (eye + Prhat)/(8.d0*pi*mu)
        K = Prhat*(3.d0*dot_product(rhat,nhat)/(4.d0*pi))
        
        call interpolate_polar_9nodequad(rr, thet, s01, s02, Nbar)
        denom_vec = matmul(xi_e,Nbar)
        denom = dsqrt(dot_product(denom_vec,denom_vec))
        
        wt = JJ/denom*Jr*Jthet*grw(nn)*gtw(mm)
        do kk = 1,9
            GKx(9*kk-8:9*kk) = GKx(9*kk-8:9*kk) + G*(N(kk)*wt)
            GKx(9*kk+73:9*kk+81) = GKx(9*kk+73:9*kk+81) + K*(Nbar(kk)*wt/denom)
        enddo
    enddo
enddo
    
end subroutine get_9NodeQuad_GKx_singular

!------------------------------------------------------------------------------
    
function rmax(thet, fn_ind)

implicit none

double precision :: thet, rmax
integer :: fn_ind

double precision, parameter :: pi = 3.141592653589793d0

select case(fn_ind)
case(1)
    rmax = 2.d0/dcos(thet)
case(2)
    rmax = 2.d0/dcos(pi/2.d0-thet)
case(3)
    rmax = 2.d0/dcos(thet-pi/2.d0)
case(4)
    rmax = 2.d0/dcos(pi-thet)
case(5)
    rmax = 2.d0/dcos(thet-pi)
case(6)
    rmax = 2.d0/dcos(3.d0*pi/2.d0-thet)
case(7)
    rmax = 2.d0/dcos(thet-3.d0*pi/2.d0)
case(8)
    rmax = 2.d0/dcos(2.d0*pi-thet)
case(9)
    rmax = 1.d0/dcos(thet)
case(10)
    rmax = 2.d0/dcos(pi/2.d0-thet)
case(11)
    rmax = 1.d0/dcos(pi-thet)
case(12)
    rmax = 1.d0/dcos(thet-pi/2.d0)
case(13)
    rmax = 2.d0/dcos(pi-thet)
case(14)
    rmax = 1.d0/dcos(3.d0*pi/2.d0-thet)
case(15)
    rmax = 1.d0/dcos(thet-pi)
case(16)
    rmax = 2.d0/dcos(3.d0*pi/2.d0-thet)
case(17)
    rmax = 1.d0/dcos(thet)
case(18)
    rmax = 1.d0/dcos(pi/2.d0+thet)
case(19)
    rmax = 2.d0/dcos(thet)
case(20)
    rmax = 1.d0/dcos(pi/2.d0-thet)
case(21)
    rmax = 1.d0/dcos(thet)
case(22)
    rmax = 1.d0/dcos(pi/2.d0-thet)
case(23)
    rmax = 1.d0/dcos(pi-thet)
case(24)
    rmax = 1.d0/dcos(3.d0*pi/2.d0-thet)
end select

end function rmax

!-----------------------------------------------------------------------
    
subroutine interpolate_polar_9nodequad(r, thet, s01, s02, Nbar)

! computes (N(s01+r*cos(theta),s02+r*sin(theta)) - N(s01,s02))/r

implicit none

double precision :: r, thet, s01, s02, Nbar(9)

double precision :: ct, st, c2t, s2t, r3c2ts2t, &
                    r2c2tst, r2s2tct, rc2t, rs2t, rstct

ct = dcos(thet)
st = dsin(thet)

c2t = ct*ct
s2t = st*st

r3c2ts2t = (r*r*r) * c2t * s2t

r2c2tst = (r*r) * c2t * st
r2s2tct = (r*r) * s2t * ct

rc2t = r * c2t
rs2t = r * s2t
rstct = r * st * ct

Nbar(1) = (r3c2ts2t &
      + (2.d0*s02-1.d0)*r2c2tst + (2.d0*s01-1.d0)*r2s2tct &
      + s02*(s02-1.d0)*rc2t + s01*(s01-1.d0)*rs2t + (2.d0*s01-1.d0)*(2.d0*s02-1.d0)*rstct &
      + (2.d0*s01-1.d0)*s02*(s02-1.d0)*ct + s01*(s01-1.d0)*(2.d0*s02-1.d0)*st)/4.d0
  
Nbar(2) = (r3c2ts2t &
      + (2.d0*s02-1.d0)*r2c2tst + (2.d0*s01+1.d0)*r2s2tct &
      + s02*(s02-1.d0)*rc2t + s01*(s01+1.d0)*rs2t + (2.d0*s01+1.d0)*(2.d0*s02-1.d0)*rstct &
      + (2.d0*s01+1.d0)*s02*(s02-1.d0)*ct + s01*(s01+1.d0)*(2.d0*s02-1.d0)*st)/4.d0
  
Nbar(3) = (r3c2ts2t &
      + (2.d0*s02+1.d0)*r2c2tst + (2.d0*s01+1.d0)*r2s2tct &
      + s02*(s02+1.d0)*rc2t + s01*(s01+1.d0)*rs2t + (2.d0*s01+1.d0)*(2.d0*s02+1.d0)*rstct &
      + (2.d0*s01+1.d0)*s02*(s02+1.d0)*ct + s01*(s01+1.d0)*(2.d0*s02+1.d0)*st)/4.d0

Nbar(4) = (r3c2ts2t &
      + (2.d0*s02+1.d0)*r2c2tst + (2.d0*s01-1.d0)*r2s2tct &
      + s02*(s02+1.d0)*rc2t + s01*(s01-1.d0)*rs2t + (2.d0*s01-1.d0)*(2.d0*s02+1.d0)*rstct &
      + (2.d0*s01-1.d0)*s02*(s02+1.d0)*ct + s01*(s01-1.d0)*(2.d0*s02+1.d0)*st)/4.d0

Nbar(5) = (-r3c2ts2t &
      - (2.d0*s02-1.d0)*r2c2tst - 2.d0*s01*r2s2tct &
      - s02*(s02-1.d0)*rc2t + (1.d0-s01)*(1.d0+s01)*rs2t - 2.d0*s01*(2.d0*s02-1.d0)*rstct &
      - 2.d0*s01*s02*(s02-1.d0)*ct + (1.d0-s01)*(1.d0+s01)*(2.d0*s02-1.d0)*st)/2.d0

Nbar(6) = (-r3c2ts2t &
      - 2.d0*s02*r2c2tst - (2.d0*s01+1.d0)*r2s2tct &
      + (1.d0-s02)*(1.d0+s02)*rc2t - s01*(s01+1.d0)*rs2t - 2.d0*(2.d0*s01+1.d0)*s02*rstct &
      + (2.d0*s01+1.d0)*(1.d0-s02)*(1.d0+s02)*ct - 2.d0*s01*(s01+1.d0)*s02*st)/2.d0

Nbar(7) = (-r3c2ts2t &
      - (2.d0*s02+1.d0)*r2c2tst - 2.d0*s01*r2s2tct &
      - s02*(s02+1.d0)*rc2t + (1.d0-s01)*(1.d0+s01)*rs2t - 2.d0*s01*(2.d0*s02+1.d0)*rstct &
      - 2.d0*s01*s02*(s02+1.d0)*ct + (1.d0-s01)*(1.d0+s01)*(2.d0*s02+1.d0)*st)/2.d0

Nbar(8) = (-r3c2ts2t &
      - 2.d0*s02*r2c2tst - (2.d0*s01-1.d0)*r2s2tct &
      + (1.d0-s02)*(1.d0+s02)*rc2t - s01*(s01-1.d0)*rs2t - 2.d0*(2.d0*s01-1.d0)*s02*rstct &
      + (2.d0*s01-1.d0)*(1.d0-s02)*(1.d0+s02)*ct - 2.d0*s01*(s01-1.d0)*s02*st)/2.d0

Nbar(9) = r3c2ts2t &
      + 2.d0*s02*r2c2tst + 2.d0*s01*r2s2tct &
      - (1.d0-s02)*(1.d0+s02)*rc2t - (1.d0-s01)*(1.d0+s01)*rs2t + 4.d0*s01*s02*rstct &
      - 2.d0*s01*(1.d0-s02)*(1.d0+s02)*ct - 2.d0*(1.d0-s01)*(1.d0+s01)*s02*st
      
end subroutine interpolate_polar_9nodequad
