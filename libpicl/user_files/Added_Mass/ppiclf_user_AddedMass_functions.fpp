!------------------------------------------------------------------------
!
! Created May 20, 2024
!
! Subroutine for added mass
!   also called the quasi-unsteady force,
!   or the inviscid unsteady force in case of the Euler equations
!
! Contains functions:
!    real*8 function B11_11(d, alpha)
!    real*8 function B11_22(d, alpha)
!    real*8 function B12_11(d, alpha)
!    real*8 function B12_22(d, alpha)
!    real*8 function IA_analytical(rmax, rad, alpha)
!    real*8 function II_analytical(rmax, rad, alpha)
!    real*8 function rdf_analytical(r, alpha)
!    real*8 function IA_numerical_integrand(d, rad, alpha)
!    real*8 function II_numerical_integrand(d, rad, alpha)
!    real*8 function IA_numerical(rmax, rad, alpha)
!    real*8 function II_numerical(rmax, rad, alpha)
!    
! Contains subroutines:
!     rotation_matrix(x, y, z, Q)
!     resistance_pair(x, y, z, alpha, rad, R)
!
! Implementing Added Mass Algorithm from S.Briney (2024)
!  
! n       = number of points
! alpha   = volume fraction
! rad     = particle radius
! d       = center-to-center distance
! rmax    = center-to-center max neighbor distance
! R       = resistance matrix (output)
! x       = x_2 - x_1
! y       = y_2 - y_1
! z       = z_2 - z_1
! dr_max  = max interaction distance between particles considered 
! poins   = 3xn array of points x, y, z. Initialized as points(3,n)
!
! correction_analytical_always 
!    = if true, always use the analytical distant neighbor correction
!    > if false, use numerical when dr_max/rad < 3.49
!
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! Functions for calculating the four curves that define the binary model
! These the the vol fraction-corrected binary added mass terms
!   taken from Appendix C of Briney et at (2024).
!
! parallel added mass
module ppiclf_user_AM_functions

    implicit none
    contains

    real*8 function B11_11(d, alpha)
    
        ! input
        real*8 d, alpha
        
        ! local vars
        ! asymptotic solution (Beguin et al. 2016)
        ! empirical higher order terms
        ! empirical volume fraction correction
        real*8 asym, hot, alpha_corr
    
        ! Beguin et al. 2016
        asym = 0.5*(3.0/64.0 * (2.0/d)**6 + 9.0/256.0*(2.0/d)**8) 

        hot = 0.059/((d-1.9098)*(d-0.4782)**2 * (d**3))
        alpha_corr = (alpha*alpha - 0.0902*alpha) * 70.7731/((d+2.0936)**6)
        
        B11_11 = asym + hot + alpha_corr
        
        return
    end function B11_11

    ! perpendicular added mass
    real*8 function B11_22(d, alpha)
    
        ! input
        real*8 d, alpha
        
        ! local vars
        ! hot combined with alpha correction in this case
        ! See comment in calc_B11_11
        real*8 asym, hot
        real*8 A, B
        
        ! Beguin et al. 2016
        asym = 0.5*(3.0/256.0 * (2.0/d)**6 + 3.0/256.0 *(2.0/d)**8) 
        
        A = 0.0003 + 0.0262*alpha*alpha - 0.0012*alpha
        B = 1.3127 + 1.0401*alpha*alpha - 1.2519*alpha
        hot = A / ((d-B)**6)
        
        B11_22 = asym + hot
        
        return
    end function B11_22
    
    ! parallel induced added mass
    real*8 function B12_11(d, alpha)
        
        ! input
        real*8 d, alpha
        
        ! local vars
        real*8 asym, hot, alpha_corr
    
        ! Beguin et al. 2016
        asym = 0.5*(-3.0/8.0*(2.0/d)**3 - 3.0/512.0 * (2.0/d)**9) 

        hot = -0.0006/((d-1.5428)**5)
        alpha_corr = -0.7913*alpha *exp(-(0.9801 - 0.1075*alpha)*d)
        
        B12_11 = asym + hot + alpha_corr

        return
    end function B12_11
      
    ! perpendicular induced added mass
    real*8 function B12_22(d, alpha)
        
        ! input
        real*8 d, alpha
        
        ! local vars
        real*8 asym, hot
        
        ! Beguin et al. 2016
        asym = 0.5*(3.0/16.0 * (2.0/d)**3 + 3.0/4096.0 * (2.0/d)**9) 
        
        ! includes alpha correction
        hot = (-0.1985 + 16.7372*alpha)/(d**6) + (1.3907 - 48.2604*alpha)/(d**8) 
        
        B12_22 = asym + hot
        
        return
    end function B12_22
    !
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    !
    !     Analytical correction functions for distant neighbors. 
    !     Assumes g(r) = 1. Good for maxr >= 3.5.
    !      
    !     These can be found in Appendix D of Briney etal (2024).
    !
    !     added mass
    !     rmax = center-to-center maximum neighbor distance
    !     rad = particle radius
    !     alpha = volume fraction
    real*8 function IA_analytical(rmax, rad, alpha)

            
        ! input
        real*8 rmax, rad, alpha
        
        ! local vars
        real*8 term1, term2, term3, c, numerator, B11, B22, maxr
    
        maxr = rmax / rad
        
        ! B11
        term1 = (5.0*maxr**2 + 9.0)/(10.0*maxr**5)
        
        term2 = ((0.00720826 - 0.0150737 * maxr) * log(maxr - 1.9098)       &
                - 0.120023 * maxr * log(maxr - 0.4782)                      &
                + 0.057395 * log(maxr - 0.4782)                             &
                + (0.135097 * maxr - 0.0646033) * log(maxr) - 0.0861828)    &
                                                  /(maxr - 0.4782)
            
        term3 = (alpha**2 - 0.0902*alpha)                                   &
                * (0.333333 * maxr**2 + 0.348933 * maxr + 0.146105)         &
                /(maxr + 2.0936)**5
            
        B11 = term1 + term2 + term3
            
        ! B22
        term1 = (5.0*maxr**2 + 12.0)/(40.0*maxr**5)
        
        numerator = 0.0262*alpha**2 - 0.0012*alpha + 0.0003
        c = -1.0401*alpha**2 + 1.2519*alpha - 1.3127
        term2 = numerator * (10.0*maxr**2 + 5.0*maxr*c + c**2)              &
                / (30.0*(maxr+c)**5)
      
        B22 = term1 + term2
        
        IA_analytical = alpha*(B11 + 2.0*B22)
        
        !     write(1,*) IA_analytical, rmax, rad, alpha,
        !    >                maxr, term1, term2, term3, B11, B22

        return
    end function IA_analytical

    ! induced added mass
    real*8 function II_analytical(rmax, rad, alpha)
      
        ! input
        real*8 rmax, rad, alpha
        
        ! local vars
        real*8 term1, term2, term3, B11, B22, c, maxr, prefactor
      
        maxr = rmax / rad ! scale the filter width
        
        ! B11
        term1 = -1.0/(4.0*maxr**6) 
        
        term2 = -6.0e-4 * (0.5*maxr**2 - 0.516067*maxr + 0.199744)/(1.5482 - maxr)**4
      
        prefactor = -0.7913*alpha
        c = 0.9801 - 0.1075*alpha
        term3 = prefactor * (maxr *c *(maxr *c + 2.0) + 2.0) * exp((-maxr *c))/c**3
      
        B11 = term1 + term2 + term3
        
        ! B22
        term1 = 1.0 / (32.0 * maxr**6)
        
        term2 = (-0.1985 + 16.7372*alpha) / (3.0*maxr**3)
        
        term3 = (1.3907 - 48.2604*alpha) / (5.0*maxr**5)
        
        B22 = term1 + term2 + term3
        
        II_analytical = alpha*(B11 + 2.0*B22)
            
        !          write(2,*) II_analytical, rmax, rad, alpha,
        !     >                maxr, term1, term2, term3, B11, B22
      
        return
    end function II_analytical
    !
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    !
    !     Numerical distant neighbor corrections
    !      
    !     Uses the rdf function
    !      
    !     r = dist between particles / diameter 
    !          - normalization different than my convention!
    !     phi = particle volume fraction
    !     Trokhymchuk et al. 2005. The Journal of Chemical Physics.
    !     https://doi.org/10.1063/1.1979488
    !     erratum: https://doi.org/10.1063/1.2188941

    real*8 function rdf_analytical(r, alpha)
      
        ! inputs
        real*8 r, alpha
      
        ! local vars
        real*8  eta, d, mu, alpha0, beta0, denom_gamma, gamma, omega, kappa
        real*8  alpha1, beta, rstar, gm, gsig, Brad, Arad, delta, Crad, g2
      
        eta = alpha
      
        d   = (2.0*eta*(eta*eta - 3.0*eta - 3.0     &
            + sqrt(3.0*(eta**4-2.0*(eta**3)         &
            +eta**2+6.0*eta+3.0))))**(1.0/3.0)
      
        !mu    = (2.0*eta/(1.0-eta))*(-1.0-(d/(2.0*eta))-(eta/d))
        mu    = (2.0*eta/(1.0-eta))*(-1.0-(d/(2.0*eta))+(eta/d)) ! see erratum
      
        alpha0 = (2.0*eta/(1.0-eta))*(-1.0+(d/(4.0*eta))-(eta/(2.0*d)))
      
        beta0 = (2.0*eta/(1.0-eta))*sqrt(3.0)*(-(d/(4.0*eta))-(eta/(2.0*d)))
      
        denom_gamma = (alpha0**2 + beta0**2 - 2.0*mu*alpha0)*(1.0 + 0.5*eta) - mu*(1.0 + 2.0*eta) ! erratum
      
        gamma = atan(-(1.0/beta0)*((alpha0*(alpha0**2+beta0**2)     &
                    -mu*(alpha0**2-beta0**2))*(1.0+0.5*eta)         &
                    +(alpha0**2+beta0**2-mu*alpha0)*(1.0+2.0*eta))  &
                    /(denom_gamma)) ! erratum
      
        !gamma = atan(-(1.0/beta0)*((alpha0*(alpha0**2+beta0**2)
        ! > -mu*(alpha0**2+beta0**2))*(1.0+0.5*eta)+
        ! > (alpha0**2+beta0**2-mu*alpha0)*(1.0+2.0*eta)));
      
        omega = -0.682*exp(-24.697*eta)+4.72+4.45*eta
        
        kappa = 4.674*exp(-3.935*eta)+3.536*exp(-56.27*eta)
        
        alpha1 = 44.554+79.868*eta+116.432*eta*eta-44.652*exp(2.0*eta)
        
        beta  = -5.022+5.857*eta+5.089*exp(-4.0*eta)
        
        rstar = 2.0116-1.0647*eta+0.0538*eta*eta
      
        gm    = 1.0286-0.6095*eta+3.5781*(eta**2)-21.3651*(eta**3)+42.6344*(eta**4)-33.8485*(eta**5)
        
        gsig  = (1.0/(4.0*eta))*(((1.0+eta+(eta**2)-(2.0/3.0)*(eta**3)-(2.0/3.0)*(eta**4))/((1.0-eta)**3))-1.0)
        
        Brad = (gm-(gsig/rstar)*exp(mu*(rstar - 1.0)))      &
                    *rstar/(cos(beta*(rstar-1.0)+gamma)     &
                    *exp(alpha1*(rstar-1.0))                &
                    -cos(gamma)*exp(mu*(rstar-1.0)))
      
        Arad = gsig - Brad*cos(gamma)
        
        delta = -omega*rstar - atan((kappa*rstar+1.0)/(omega*rstar))
        
        Crad = rstar*(gm-1.0)*exp(kappa*rstar)/cos(omega*rstar+delta)
        
        g2 = (Arad/r)*exp(mu*(r-1.0)) + (Brad/r)*cos(beta*(r-1.0)+gamma)*exp(alpha1*(r-1.0))
      
        if (r > rstar) then
            g2 = 1.0 +(Crad/r)*cos(omega*r+delta)*exp(-kappa*r)
        else if (r < 1) then
            g2 = 0.0
        end if
      
        rdf_analytical = g2
        
        return
    end function rdf_analytical
    !
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    !
    !     d = center to center distance
    !     rad = radius
    !     alpha = volume fraction

    real*8 function IA_numerical_integrand(d, rad, alpha)
      
        ! input
        real*8 d, rad, alpha
        
        ! local vars
        real*8 g, f11, f22, rho, pi, w

        ! ! declare funcitons used
        ! real*8 rdf_analytical, B11_11, B11_22
      
        pi = 4.0*ATAN(1.0)
        
        g = rdf_analytical((d/(2.0*rad)), alpha)
        
        f11 = B11_11(d/rad, alpha)
        f22 = B11_22(d/rad, alpha)
        
        rho = alpha / (4.0/3.0*pi) ! number density
        w = 1.0/3.0 * rho * 4.0*pi*(d/rad)**2
        
        IA_numerical_integrand = w*g*(f11 + 2.0*f22)

        return
    end function IA_numerical_integrand
    !
    !-----------------------------------------------------------------------
    !-----------------------------------------------------------------------
    !
    !     d = center to center distance
    !     rad = radius
    !     alpha = volume fraction

    real*8 function II_numerical_integrand(d, rad, alpha)
      
        ! input
        real*8 d, rad, alpha
        
        ! local vars
        real*8 g, f11, f22, rho, pi, w
            
        ! ! declare funcitons used
        ! real*8 rdf_analytical, B12_11, B12_22
      
        pi = 4.0*ATAN(1.0)
        
        g = rdf_analytical((d/(2.0*rad)), alpha)
        
        f11 = B12_11(d/rad, alpha)
        f22 = B12_22(d/rad, alpha)
        
        rho = alpha / (4.0/3.0*pi) ! number density
        w = 1.0/3.0 * rho * 4.0*pi*(d/rad)**2
        
        II_numerical_integrand = w*g*(f11 + 2.0*f22)

        return
    end function II_numerical_integrand
    !
    !-----------------------------------------------------------------------
    !
      
    real*8 function IA_numerical(rmax, rad, alpha)
      
        ! input
        real*8 rmax, rad, alpha
        
        ! local vars
        real*8 maxr, dr, r, integral, f, coef(2)
        integer npts, i, j
            
        ! ! declare funcitons used
        ! real*8 IA_numerical_integrand
        
        npts = 1001 ! must be odd
        
        maxr = rad*20.0 ! essentially infinity
        
        coef(1) = 4.0
        coef(2) = 2.0
        
        if (maxr > rmax) then
            dr = (maxr - rmax) / (npts + 1)
        
            integral = 0.0
            r = rmax
        
            ! Simpson's rule
            ! endpoints first
            f = IA_numerical_integrand(r, rad, alpha)
            integral = integral + f
        
            f = IA_numerical_integrand(maxr, rad, alpha)
            integral = integral + f
        
            ! middle points
            do i=2,npts-1
                r = r + dr
                f = IA_numerical_integrand(r, rad, alpha)
        
                j = mod(i, 2) + 1
                integral = integral + coef(j)*f ! 4, 2, 4, 2, ...
            end do
        
            IA_numerical = dr * integral / 3.0
        
        else
            IA_numerical = 0.0
        end if
        
        return
    end function IA_numerical
      
    !
    !-----------------------------------------------------------------------
    !
      
    real*8 function II_numerical(rmax, rad, alpha)
      
        ! input
        real*8 rmax, rad, alpha
        
        ! local vars
        real*8 maxr, dr, r, integral, f, coef(2)
        integer npts, i, j
        
        ! ! declare funcitons used
        ! real*8 II_numerical_integrand
          
        npts = 1001 ! must be odd
        
        maxr = rad*20.0 ! essentially infinity
        
        coef(1) = 4.0
        coef(2) = 2.0
        
        if (maxr > rmax) then
            dr = (maxr - rmax) / (npts + 1)
        
            integral = 0.0
            r = rmax
        
            ! Simpson's rule
            ! endpoints first
            f = II_numerical_integrand(r, rad, alpha)
            integral = integral + f
        
            f = II_numerical_integrand(maxr, rad, alpha)
            integral = integral + f
        
            ! middle points
            do i=2,npts-1
                r = r + dr
                f = II_numerical_integrand(r, rad, alpha)
        
                j = mod(i, 2) + 1
                integral = integral + coef(j)*f ! 4, 2, 4, 2, ...
            end do
        
            II_numerical = dr * integral / 3.0
        
        else
            II_numerical = 0.0
        end if
        
        return 
    end function II_numerical


    !-----------------------------------------------------------------------
    !
    ! Calculates the rotation matrix Q to align the coordinate 
    !    system such that the point (x, y, z) is on the x-axis

    subroutine rotation_matrix(x, y, z, Q)
       
        ! inputs
        real*8 x, y, z
        
        ! output rotation matrix
        real*8 Q(3, 3)
        
        ! local vars
        real*8 gamma, beta
           
        gamma =  atan2(y, x)
        beta  = -atan2(z, sqrt(x*x + y*y))
        
        Q(1, 1) = cos(beta)*cos(gamma)
        Q(1, 2) = cos(beta)*sin(gamma)
        Q(1, 3) = -sin(beta)
        
        Q(2, 1) = -sin(gamma)
        Q(2, 2) = cos(gamma)
        Q(2, 3) = 0.0
        
        Q(3, 1) = sin(beta)*cos(gamma)
        Q(3, 2) = sin(beta)*sin(gamma)
        Q(3, 3) = cos(beta)
        
        return
    end subroutine rotation_matrix
    !
    !-----------------------------------------------------------------------
    !
    ! This is the one of the functions the user should typically call.
    ! Returns the resistance matrix for a 2 particle system 
    !    of arbitrary orientation
    ! x = x2 - x1, y = y2 - y1, z = z2 - z1 (relative position)
    ! rad = particle radius
    ! R = resistance matrix (output)

    subroutine resistance_pair(x, y, z, alpha, rad, R)
       
        ! input: x, y, z of second particle relative to the second particle
        ! x = x2 - x1, etc.
        ! rad = particle radius (monodisperse)
        real*8 x, y, z, alpha, rad
       
        ! output: resistance matrix
        real*8 R(6, 6)
       
        ! local vars
        real*8, dimension(3, 3) :: Q, Qt, B11, B12
        
        real*8 dist
        
        ! ! declare functions
        ! real*8 B11_11, B11_22, B12_11, B12_22
       
        ! get rotation matrix Q
        call rotation_matrix(x, y, z, Q)
        
        ! normalize the distance by the particle radius
        dist = sqrt(x*x + y*y + z*z) / rad 
        
        ! initially set coefficients to 0
        B11(1:3, 1:3) = 0.0
        B12(1:3, 1:3) = 0.0
        
        ! self acceleration
        B11(1, 1) = B11_11(dist, alpha)
        B11(2, 2) = B11_22(dist, alpha)
        B11(3, 3) = B11(2, 2)
       
        ! neighbor acceleration
        B12(1, 1) = B12_11(dist, alpha)
        B12(2, 2) = B12_22(dist, alpha)
        B12(3, 3) = B12(2, 2)
        
        ! rotate
        Qt  = transpose(Q)
        B11 = matmul(Qt, matmul(B11, Q))
        B12 = matmul(Qt, matmul(B12, Q))
            
        ! store resitance matrix
        R(1:3, 1:3) = B11
        R(1:3, 4:6) = B12
        
        R(4:6, 1:3) = B12
        R(4:6, 4:6) = B11
          
        ! write(7059,*) x, y, z, alpha, rad, dist 
        !
        ! write(7060,*) B11(1,1:3), B11(2,1:3), B11(3,1:3)
        ! write(7061,*) B12(1,1:3), B12(2,1:3), B12(3,1:3)
        ! write(7062,*) C11(1,1:3), C11(2,1:3), C11(3,1:3)
        ! write(7063,*) C12(1,1:3), C12(2,1:3), C12(3,1:3)
        ! write(7064,*) Qt(1,1:3), Qt(2,1:3), Qt(3,1:3)
        !
        ! write(7069,*) R(1,1:6)
        ! write(7069,*) R(2,1:6)
        ! write(7069,*) R(3,1:6)
        ! write(7069,*) R(4,1:6)
        ! write(7069,*) R(5,1:6)
        ! write(7069,*) R(6,1:6)
        ! write(7069,*) "------------------------------------"

        return
    end subroutine resistance_pair

end module ppiclf_user_AM_functions