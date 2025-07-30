subroutine potem(t, s, p, theta)
 
!=======================================================================
!     this subroutine calculates potential temperature as a function
!     of in-situ temperature, salinity, and pressure.
!
!     input [units]:
!       in-situ temperature (t): [degrees centigrade]
!       salinity (s): [per mil]
!       pressure (p): [decibars, approx. as meters of depth]
!     output [units]:
!       potential temperature (theta): [degrees centigrade]
!
!     references:
!        based on Fofonoff and Froese (1958) as shown in ...
!        Fofonoff, N., The Sea: Vol 1, (ed. M. Hill). Interscience,
!          New York, 1962, page 17, table iv.
!
!=======================================================================
 
    implicit double precision(a-h, o-z)
 
    b1    = -1.60d-5 * p
    b2    = 1.014d-5 * p * t
    t2    = t * t
    t3    = t2 * t
    b3    = -1.27d-7 * p * t2
    b4    = 2.7d-9 * p * t3
    b5    = 1.322d-6 * p * s
    b6    = -2.62d-8 * p * s * t
    s2    = s * s
    p2    = p * p
    b7    = 4.1d-9 * p * s2
    b8    = 9.14d-9 * p2
    b9    = -2.77d-10 * p2 * t
    b10   = 9.5d-13 * p2 * t2
    b11   = -1.557d-13 * p2 * p
    potmp = b1 + b2 + b3 + b4 + b5 + b6 + b7 + b8 + b9 + b10 + b11
    theta = t - potmp
 
    return
end subroutine potem




subroutine unesco(t, s, pin, rho)
 
!=======================================================================
!     this subroutine calculates the density of seawater using the
!     standard equation of state recommended by unesco(1981).
!
!     input [units]:
!       in-situ temperature (t): [degrees centigrade]
!       salinity (s): [practical salinity units]
!       pressure (pin): [decibars, approx. as meters of depth]
!     output [units]:
!       density(rho): kilograms per cubic meter
!
!     references:
!        Gill, A., Atmosphere-Ocean Dynamics: International Geophysical
!         Series No. 30. Academic Press, London, 1982, pp 599-600.
!        UNESCO, 10th report of the joint panel on oceanographic tables
!          and standards. UNESCO Tech. Papers in Marine Sci. No. 36,
!          Paris, 1981.
!
!=======================================================================
 
    implicit double precision(a-h, o-z)
 
    c1p5 = 1.5d0
 
    ! convert from depth [m] (decibars) to bars
    p = pin * 1.0d-1
 
    rw = 9.99842594d2 + 6.793952d-2 * t - 9.095290d-3 * t**2  &
       + 1.001685d-4 * t**3 - 1.120083d-6 * t**4 + 6.536332d-9 * t**5
 
    rsto = rw + (8.24493d-1 - 4.0899d-3 * t + 7.6438d-5 * t**2 - 8.2467d-7 * t**3 + 5.3875d-9 * t**4) * s  &
         + (-5.72466d-3 + 1.0227d-4 * t - 1.6546d-6 * t**2) * s**c1p5 + 4.8314d-4 * s**2
 
    xkw = 1.965221d4 + 1.484206d2 * t - 2.327105d0 * t**2 + 1.360477d-2 * t**3 - 5.155288d-5 * t**4
 
    xksto = xkw + (5.46746d1 - 6.03459d-1 * t + 1.09987d-2 * t**2 - 6.1670d-5*t**3) * s  &
          + (7.944d-2 + 1.6483d-2 * t - 5.3009d-4 * t**2) * s**c1p5
 
    xkstp = xksto + (3.239908d0 + 1.43713d-3 * t + 1.16092d-4 * t**2 - 5.77905d-7 * t**3) * p  &
          + (2.2838d-3 - 1.0981d-5 * t - 1.6078d-6 * t**2) * p * s   &
          + 1.91075d-4 * p * s**c1p5                                 &
          + (8.50935d-5 - 6.12293d-6 * t + 5.2787d-8 * t**2) * p**2  &
          + (-9.9348d-7 + 2.0816d-8*t + 9.1697d-10*t**2) * p**2 * s
 
    rho = rsto / (1.0d0 - p/xkstp)
 
    return
end subroutine unesco
