MODULE container_mod

! File: container_mod.f90
! Associated Program(s): container_main.f90
! Description: Subroutines for container modeling
! Dan Steinhoff
! 13 March 2013

! Revision history:
!   2023-08-08, CNY: added manual emptying (mempty and empty_intv)
!   2023-08-24, CNY: fixed addition of manual emptying in read_params()
!   2023-09-06, CNY: fixed typos in comments where "water" was used instead of "air"
!                    (referring to the temperature values used in calc_tg_init() and calc_tg_wave())
!   2023-09-11, CNY: edited read_params() to handle longer filenames (formerly up to 100 char, now up to 300 char)
!                     - edited variable declarations for data_file, outfile (CHARACTER(len=100) -> CHARACTER(len=300))
!                     - edited string formatting for data_file, outfile (A100 -> A300)
!   2023-11-17, CNY: added subroutine calc_vpd


! Procedure List:
! read_params
! calc_area
! calc_tg_init
! calc_tg_wave
! calc_z
! calc_tau
! calc_svp
! calc_vp
! calc_vpd  (ADDED ON 2023-11-17, CNY)
! calc_swdown
! calc_swside
! calc_lwdown
! calc_lwup
! calc_lwsidein
! calc_lwsideout
! calc_htop
! calc_hside
! calc_lh
! calc_ghc
! calc_ghw
! calc_cond
! calc_evap
! update_wh
! update_tw
! update_tcon

IMPLICIT NONE

!!! Constant Declarations
! The thermodynamic and physics constants here come from a combination of sources, including Monteith and Unsworth (2008), Stull (1988), and Oke (1987)
! Thermodynamic
REAL, PARAMETER :: RHO_W = 1.0E03        ! Density of water at 20C (kg m-3)
REAL, PARAMETER :: KA = 2.06E-05         ! Thermal diffusivity of air (m2 s-1)
REAL, PARAMETER :: KG = 3.40E-07         ! Thermal diffusivity of clay soil, 50% saturated (m2 s-1)
REAL, PARAMETER :: K_G = 0.915           ! Thermal conductivity of soil, assuming clay 50% saturated (W m-1 K-1).
REAL, PARAMETER :: K_W = 0.57            ! Thermal conductivity of water (W m-1 K-1)
REAL, PARAMETER :: PSY = 0.067           ! Psychrometer constant (kPa K-1)
REAL, PARAMETER :: LV = 2.45E06          ! Latent heat of vaporization at 20C (J kg-1)
REAL, PARAMETER :: RV = 461.             ! Gas constant for water vapor (J K-1 kg-1)
REAL, PARAMETER :: HC_A = 1.23E03        ! Volumetric heat capacity of air (J m-3 K-1)
REAL, PARAMETER :: HC_W = 4.18E06        ! Volumetric heat capacity of water (J m-3 K-1)
REAL, PARAMETER :: HC_C = 2.20E06        ! Volumetric heat capacity of container (J m-3 K-1)
! Radiation
REAL, PARAMETER :: SB = 5.67E-08         ! Stefan-Boltzmann constant (W m-2 K-4)
REAL, PARAMETER :: REF_G = 0.2           ! Reflectivity of ground surface
REAL, PARAMETER :: R_W = 1.33            ! Refractive index of water
REAL, PARAMETER :: SC = 1366.            ! Solar constant (W m-2)
REAL, PARAMETER :: E_W = 0.95            ! Emissivity of water
REAL, PARAMETER :: E_G = 0.95            ! Emissivity of ground surface (soil, grass, etc.)
REAL, PARAMETER :: E_C = 0.95            ! Emissivity of container
REAL, PARAMETER :: E_L = 0.95            ! Emissivity of cloud
REAL, PARAMETER :: E_B = 0.90            ! Emissivity of shelter, could be tree (Em~0.98) or concrete roof (Em~0.82).  Use some middle value.
! Heat Transfer
REAL, PARAMETER :: B1 = 0.5              ! Nusselt number coefficient, horizontal surface
REAL, PARAMETER :: A1 = 0.25             ! Nusselt number exponent, horizontal surface
REAL, PARAMETER :: B2 = 0.58             ! Nusselt number coefficient, vertical surface, laminar flow
REAL, PARAMETER :: A2 = 0.25             ! Nusselt number exponent, vertical surface, laminar flow
! Experiment-specific
REAL, PARAMETER :: WHMIN = 0.015         ! Minimum water level allowed (m)
REAL, PARAMETER :: WHFLOOR = 0.          ! Water height floor (m)
REAL, PARAMETER :: P = 86400.            ! Period (seconds in one day)
REAL, PARAMETER :: DZG = 0.05            ! Depth of soil temperature observation/estimate (m)
REAL, PARAMETER :: EVAP_LOW = 0.02       ! Evaporation rate when water level too low (mm/interval of input data)
! Misc.
REAL, PARAMETER :: MISS = -9999.         ! Missing value
REAL, PARAMETER :: PI = ACOS(-1.)        ! pi
SAVE

CONTAINS

    SUBROUTINE read_params (param_file, data_file, city, station, lat, lon, elev, cshape, rt1, rt2, rb1, rb2, h, refl, k_c, thick, wh, tw, tg, time_step, t_diff, shade, cflag, sflag, mfill, fill_intv, mempty, empty_intv, outfile, cld, cmd, chd, cln, cmn, chn)  ! ADDED mempty and empty_intv (2023-08-08, CNY)
    ! Description: Reads experiment parameters from a text file into program
    ! Input: Parameters file from calling program (param_file)
    ! Output: Parameters described below
    IMPLICIT NONE

    ! Variable Declarations
    CHARACTER(len=*), INTENT(IN) :: param_file          ! Input file name
    CHARACTER(len=300), INTENT(OUT) :: data_file        ! Data file name
    CHARACTER(len=300), INTENT(OUT) :: outfile          ! Output file
    CHARACTER(len=30), INTENT(OUT) :: city              ! City
    CHARACTER(len=30), INTENT(OUT) :: station           ! Station
    REAL, INTENT(OUT) :: lat                            ! Latitude
    REAL, INTENT(OUT) :: lon                            ! Longitude
    REAL, INTENT(OUT) :: elev                           ! Elevation
    CHARACTER(len=5), INTENT(OUT) :: cshape             ! Container shape
    REAL, INTENT(OUT) :: rt1                            ! Container radius 1 (Top)
    REAL, INTENT(OUT) :: rt2                            ! Container radius 2 (Top)
    REAL, INTENT(OUT) :: rb1                            ! Container radius 1 (Body)
    REAL, INTENT(OUT) :: rb2                            ! Container radius 2 (Body)
    REAL, INTENT(OUT) :: h                              ! Container height
    REAL, INTENT(OUT) :: refl                           ! Container reflectivity
    REAL, INTENT(OUT) :: k_c                            ! Container thermal conductivity
    REAL, INTENT(OUT) :: thick                          ! Container thickness
    REAL, INTENT(OUT) :: wh                             ! Initial water height
    REAL, INTENT(OUT) :: tw                             ! Initial water temperature
    REAL, INTENT(OUT) :: tg                             ! Ground temperature
    REAL, INTENT(OUT) :: time_step                      ! Model time step
    REAL, INTENT(OUT) :: t_diff                         ! Ratio of input data interval to time step
    REAL, INTENT(OUT) :: shade                          ! Shade fraction
    LOGICAL, INTENT(OUT) :: cflag                       ! Cloud data input
    LOGICAL, INTENT(OUT) :: sflag                       ! Soil temperature input
    LOGICAL, INTENT(OUT) :: mfill                       ! Manual fill flag
    INTEGER, INTENT(OUT) :: fill_intv                   ! Fill interval
    LOGICAL, INTENT(OUT) :: mempty                      ! Manual empty flag   (ADDED ON 2023-08-08, CNY)
    INTEGER, INTENT(OUT) :: empty_intv                  ! Empty interval      (ADDED ON 2023-08-08, CNY)
    REAL, INTENT(OUT), OPTIONAL :: cld                  ! Low-cloud fraction - Day
    REAL, INTENT(OUT), OPTIONAL :: cmd                  ! Middle-cloud fraction - Day
    REAL, INTENT(OUT), OPTIONAL :: chd                  ! High-cloud fraction - Day
    REAL, INTENT(OUT), OPTIONAL :: cln                  ! Low-cloud fraction - Night
    REAL, INTENT(OUT), OPTIONAL :: cmn                  ! Middle-cloud fraction - Night
    REAL, INTENT(OUT), OPTIONAL :: chn                  ! High-cloud fraction - Night
    INTEGER :: iostat, status                           ! Status variables

    ! Open parameters file
    OPEN (UNIT=100, FILE=param_file, STATUS='OLD', ACTION='READ', IOSTAT=iostat)

    ! If open successful, then continue
    fileopen: IF (iostat == 0) THEN
        READ (100,1000,IOSTAT=status) data_file, city, station, lat, lon, elev, cshape, rt1, rt2, rb1, rb2, h, refl, k_c, thick, wh, tw, tg, time_step, t_diff, shade, cflag, sflag, mfill, fill_intv, mempty, empty_intv, outfile, cld, cmd, chd, cln, cmn, chn     ! ADDED mempty and empty_intv (2023-08-08, CNY)
        1000 FORMAT (A300,/,A30,/,A30,/,F10.4,/,F10.4,/,F10.4,/,A5,/,F10.4,/,F10.4,/,F10.4,/,F10.4,/,F10.4,/,F10.4,/,F10.4,/,F10.4,/,F10.4,/,F10.4,/,F10.4,/,F10.4,/,F10.4,/,F10.4,/,L,/,L,/,L,/,I6,/,L,/,I6,/,A300,/,F10.4,/,F10.4,/,F10.4,/,F10.4,/,F10.4,/,F10.4) ! ADDED formatting for mempty and empty_intv (2023-08-24, CNY)
        readif: IF (status > 0) THEN
            WRITE (*,*) 'An error occurred reading array'
        END IF readif

    ELSE fileopen
        WRITE (*,1010) iostat
        1010 FORMAT (1X, 'File open failed for values: ', I6)
    ENDIF fileopen

    CLOSE (UNIT=100)

    END SUBROUTINE read_params

    SUBROUTINE calc_area (rt1, rt2, rb1, rb2, h, thick, cshape, area_top, area_bot, area_side, area_beam, vol_cont)
    ! Description: Calculates container areas for top, sides, and a direct beam
    ! Input: Radius/length values and height
    ! Output: Areas
    IMPLICIT NONE

    ! Variable Declarations
    REAL, INTENT(IN) :: rt1                              ! Top radius / length 1
    REAL, INTENT(IN) :: rt2                              ! Top radius / length 2
    REAL, INTENT(IN) :: rb1                              ! Body radius / length 1
    REAL, INTENT(IN) :: rb2                              ! Body radius / length 2
    REAL, INTENT(IN) :: h                                ! Height
    REAL, INTENT(IN) :: thick
    CHARACTER(len=*), INTENT(IN) :: cshape               ! Container shape (ROUND, RECT, TIRE)
    REAL, INTENT(OUT) :: area_top                        ! Area of top of container
    REAL, INTENT(OUT) :: area_bot                        ! Area of bottom of container
    REAL, INTENT(OUT) :: area_side                       ! Area of side of container
    REAL, INTENT(OUT) :: area_beam                       ! Area of a direct beam of container
    REAL, INTENT(OUT) :: vol_cont                        ! Volume of container itself

    IF (cshape == 'ROUND') THEN
        area_top = PI*rt1*rt2
        area_bot = PI*rb1*rb2
        area_side = 2.*PI*rb1*h
        area_beam = 2.*rb1*h
        vol_cont = (PI*h*thick*((2.*rb1)+thick))+(PI*rb1*rb2*thick) ! First term is sides, second term is bottom
    ELSE IF (cshape == 'RECT') THEN
        area_top = rt1*rt2
        area_bot = rb1*rb2
        area_side = 2.*(rb1+rb2)*h
        area_beam = (rb1+rb2)*h
        vol_cont =(h*thick*(rb1+rb2+thick))+(rb1*rb2*thick) ! First term is sides, second term is bottom
    ELSE IF (cshape == 'TIRE') THEN
        area_top = 2.*PI*((rt1**2)-(rt2**2))
        area_bot = 2.*PI*((rb1**2)-(rb2**2))
        area_side = 2.*PI*rb1*h
        area_beam = 2.*rb1*h
        vol_cont = 2.*PI*h*thick*((2.*rb1)+thick) ! No bottom area of tire
    END IF

    END SUBROUTINE calc_area

    SUBROUTINE calc_tg_init (ta, ta_avg, tg_clim, cnt, ta_array, tg)
    ! Description: Calculates ground temperature estimate using recent air temperature average
    ! Input: Existing ground temperature, air temperature, number of times read
    ! Output: Updated ground temperature
    IMPLICIT NONE

    ! Variable Declarations
    REAL, INTENT(IN) :: ta                               ! Air temperature (K)
    REAL, INTENT(INOUT) :: ta_avg                        ! Running average of air temperature
    REAL, INTENT(IN) :: tg_clim                          ! Ground temperature climatology for first part of record
    INTEGER, INTENT(IN) :: cnt                           ! Count of number of times read
    REAL, INTENT(INOUT), DIMENSION(:) :: ta_array        ! Array of air temperatures over past 24 hours, running  ! Fixed typo, "water"->"air" (2023-09-06, CNY)
    REAL, INTENT(OUT) :: tg                              ! Ground Temperature (K)

    ta_avg = ((ta_avg*REAL(cnt-1.))+ta)/REAL(cnt)
    tg = ((tg_clim*(24.-REAL(cnt)))+(ta_avg*REAL(cnt)))/24.
    ta_array(INT(cnt)) = ta

    END SUBROUTINE calc_tg_init

    SUBROUTINE calc_tg_wave (ta, cnt, tod, time_step, ta_array, tg)
    ! Description: Calculates ground temperature estimate using ground temperature wave - offset in time, and dampened compared to air temperature cycle
    ! Input: Air temperature, counter, hour of day, hourly temperatures over previous day
    ! Output: Updated ground temperature
    IMPLICIT NONE

    ! Variable Declarations
    REAL, INTENT(IN) :: ta                               ! Air temperature (K)     ! Fixed typo, "water"->"air" (2023-09-06, CNY)
    INTEGER, INTENT(IN) :: cnt                           ! Count of number of times read
    REAL, INTENT(IN) :: tod                              ! Hour of day
    REAL, INTENT(IN) :: time_step                        ! Model time step
    REAL, INTENT(INOUT), DIMENSION(:) :: ta_array        ! Array of air temperatures over past 24 hours, running  ! Fixed typo, "water"->"air" (2023-09-06, CNY)
    REAL, INTENT(OUT) :: tg                              ! Ground temperature (K)
    INTEGER :: ind                                       ! Array index
    REAL :: d                                            ! Depth term in temperature equation (m)
    REAL :: t                                            ! Time variable in temperature equation (s)
    REAL :: tg_bar                                       ! Average ground temperature over a day, use air temperature, assume equilibrium  ! Fixed typo, "water"->"air" (2023-09-06, CNY)
    REAL :: as                                           ! Amplitude of ground temperature over a day, again use air temperature           ! Fixed typo, "water"->"air" (2023-09-06, CNY)

    ! Update our running 24 hour air temperature array  ! Fixed typo, "water"->"air" (2023-09-06, CNY)
    ind = MOD(REAL(cnt),24.)
    IF (ind == 0.) THEN
        ind = 24.
    END IF
    ta_array(INT(ind)) = ta

    ! Get daily average
    tg_bar = SUM(ta_array)/24.

    ! Get daily amplitude
    as = (MAXVAL(ta_array)-MINVAL(ta_array))/2.

    ! Get updated tg, this is based on sinusoidal daily ground temperature wave
    ! See Arya (2001), p.53 for derivation
    d = SQRT(P*KG/PI)
    t = (tod*time_step)-39600.    ! Taking time relative to 11 AM, which we assume is when tg=tg_bar
    tg = tg_bar+as*EXP(-1.*DZG/d)*SIN((2.*PI*t/P)-(DZG/d))

    END SUBROUTINE calc_tg_wave

    SUBROUTINE calc_z (lat, doy, tod, z)
    ! Description: Calculates solar zenith angle
    ! Input: Latitude, day of year, time of day (local hours)
    ! Output: Solar zenith angle
    IMPLICIT NONE

    ! Variable Declarations
    REAL, INTENT(IN) :: lat                               ! Latitude
    REAL, INTENT(IN) :: doy                               ! Day of year
    REAL, INTENT(IN) :: tod                               ! Time of day (local hours) from 1200
    REAL, INTENT(OUT) :: z                                ! Solar zenith angle (radians)
    REAL :: dec                                           ! Solar declination
    REAL :: ha                                            ! Hour angle

    ! Calculate solar declination
    dec = -0.4093*COS(2.*PI*(doy+10.)/365.) ! Cooper (1969)

    ! Calculate hour angle
    ha = ABS(2.*PI*((12.-tod)/24.))

    ! Calculate zenith angle
    z = ACOS(SIN(dec)*SIN(2.*PI*lat/360.)+COS(dec)*COS(2.*PI*lat/360.)*COS(ha))

    IF (z > (PI/2.)) THEN
        z = PI/2.
    END IF

    IF (z < 0.) THEN
        z = 0.
    END IF

    END SUBROUTINE calc_z

    SUBROUTINE calc_tau (cl, cm, ch, z, tau, tau_cs)
    ! Description: Calculates an estimate of optical thickness for use in SW and LW calculations.  Uses a simple parameterization from Burridge and Gadd (1974)
    ! Input: Cloud fractions for low, middle, and high levels, and solar zenith angle
    ! Output: Optical thickness (tau)
    IMPLICIT NONE

    ! Variable Declarations
    REAL, INTENT(IN) :: cl                              ! Low-cloud fraction
    REAL, INTENT(IN) :: cm                              ! Middle-cloud fraction
    REAL, INTENT(IN) :: ch                              ! High-cloud fraction
    REAL, INTENT(IN) :: z                               ! Solar zenith angle
    REAL, INTENT(OUT) :: tau                            ! Optical thickness
    REAL, INTENT(OUT) :: tau_cs                         ! Clear-sky optical thickness

    tau = -1*ALOG((0.6+(0.2*COS(z)))*(1.-(0.4*ch))*(1.-(0.7*cm))*(1.-(0.4*cl)))
    tau_cs = -1*ALOG((0.6+(0.2*COS(z)))) ! Clear-sky tau value

    END SUBROUTINE calc_tau

    SUBROUTINE calc_svp (ta, svp)
    ! Description: Calculates saturation vapor pressure
    ! Input: Air temperature (K)
    ! Output: Saturation vapor pressure (kPa)
    IMPLICIT NONE

    ! Variable Declarations
    REAL, INTENT(IN) :: ta                              ! Air temperature (K)
    REAL, INTENT(OUT) :: svp                            ! Saturation vapor pressure (kPa)

    svp = 0.611*EXP((LV/RV)*((1/273.)-(1./ta)))

    END SUBROUTINE calc_svp

    SUBROUTINE calc_vp (svp, rh, vp)
    ! Description: Calculates vapor pressure
    ! Input: Saturation vapor pressure (kPa), relative humidity (%)
    ! Output: Vapor pressure (kPa)
    IMPLICIT NONE

    ! Variable Declarations
    REAL, INTENT(IN) :: svp                             ! Saturation vapor pressure (kPa)
    REAL, INTENT(IN) :: rh                              ! Relative humidity (%)
    REAL, INTENT(OUT) :: vp                             ! Vapor pressure (kPa)

    vp = svp*(rh/100.)

    END SUBROUTINE calc_vp

    SUBROUTINE calc_vpd (svp, vp, vpd)   ! (ADDED ON 2023-11-17, CNY)
    ! Description: Calculates vapor pressure deficit
    ! Input: Saturation vapor pressure (kPa), Vapor pressure (kPa)
    ! Output: Vapor pressure deficit (kPa)
    IMPLICIT NONE

    ! Variable Declarations
    REAL, INTENT(IN) :: svp                             ! Saturation vapor pressure (kPa)
    REAL, INTENT(IN) :: vp                              ! Vapor pressure (kPa)
    REAL, INTENT(OUT) :: vpd                            ! Vapor pressure deficit (kPa)

    vpd = svp-vp

    END SUBROUTINE calc_vpd

    SUBROUTINE calc_swdown (area_top, tau, tau_cs, z, doy, shade, clf, st, swdown)
    ! Description: Calculates downward shortwave power
    ! Input: Container area, optical thickness, solar zenith angle, day of year
    ! Output: Downward shortwave radiative power (and flux for use in swside)
    IMPLICIT NONE

    ! Variable Declarations
    REAL, INTENT(IN) :: area_top                        ! Container area - Top
    REAL, INTENT(IN) :: tau                             ! Optical thickness
    REAL, INTENT(IN) :: tau_cs                          ! Clear-sky optical thickness
    REAL, INTENT(IN) :: z                               ! Solar zenith angle
    REAL, INTENT(IN) :: doy                             ! Day of year
    REAL, INTENT(IN) :: shade                           ! Shade fraction
    REAL, INTENT(OUT) :: clf                            ! Cloud fraction function
    REAL, INTENT(OUT) :: st                             ! Total downward shortwave flux (W m-2)
    REAL, INTENT(OUT) :: swdown                         ! Downward shortwave power (W)
    REAL :: tr                                          ! Atmospheric transmissivity
    REAL :: tr_cs                                       ! Clear-sky transmissivity
    REAL :: zd                                          ! Zenith angle in degrees
    REAL :: am                                          ! Air mass number
    REAL :: df                                          ! Distance factor
    REAL :: r                                           ! Angle of refraction of water
    REAL :: ref_w                                       ! Reflectivity of water

    ! Calculate the transmissivity from the optical thickness
    tr = EXP(-1.*tau)
    tr_cs = EXP(-1.*tau_cs)

    ! Calculate air mass number, using equation from Kasten and Young (1989)
    zd = z*(360./(2.*PI))
    am = (1./(ABS(COS(z))+0.50572*(96.07995-zd)**(-1.6364)))

    ! Calculate distance factor
    df = (1.+0.01673*COS(2.*PI*(doy-1.)/365.))**2

    ! Calculate an estimate of the total downward shortwave flux (W m-2)
    st = SC*df*(tr**am)*COS(z)

    ! Calculate angle of refraction of water
    r = ASIN(SIN(z)/R_W)

    ! Calculate an estimate of the albedo of water, based on solar zenith angle (i.e., Fresnel albedo, see Cogley 1979)
    ref_w = 0.5*((((SIN(z-r))**2)/((SIN(z+r))**2))+(((TAN(z-r))**2)/((TAN(z+r))**2)))

    ! Calculate downward shortwave radiative power, by getting the area of the water surface and taking into account reflectivity of water.
    swdown = st*area_top*(1.-ref_w)*(1.-shade)

    IF (swdown < 0.) THEN
        swdown = 0.
    END IF

    ! For longwave calculations, calculate clf
    clf = 1.-(tr/tr_cs)

    END SUBROUTINE calc_swdown

    SUBROUTINE calc_swside (refl, area_side, area_beam, z, tau, st, shade, ssd, ssb, swside)
    ! Description: Calculates sideways shortwave power, by separating the total shortwave into direct and diffuse components.
    ! Input: Container reflectivity, container area, optical thickness, solar zenith angle, total shortwave, shade
    ! Output: Sideways shortwave radiative power (along with direct and diffuse components)
    IMPLICIT NONE

    ! Variable Declarations
    REAL, INTENT(IN) :: refl                            ! Container reflectivity
    REAL, INTENT(IN) :: area_side                       ! Container area - Side
    REAL, INTENT(IN) :: area_beam                       ! Container area - Beam
    REAL, INTENT(IN) :: z                               ! Solar zenith angle
    REAL, INTENT(IN) :: tau                             ! Optical thickness
    REAL, INTENT(IN) :: st                              ! Total downward shortwave flux (W m-2)
    REAL, INTENT(IN) :: shade                           ! Shade
    REAL, INTENT(OUT) :: ssd                            ! Diffuse solar component
    REAL, INTENT(OUT) :: ssb                            ! Direct solar component
    REAL, INTENT(OUT) :: swside                         ! Sideways shortwave power (W)
    REAL :: sd                                          ! Diffuse shortwave flux (W m-2)
    REAL :: sb                                          ! Direct shortwave flux (W m-2)

    ! Calculate diffuse component, using parameterization in Monteith and Unsworth (2008)
    sd = st*((0.68*tau)+0.10)

    ! Calculate diffuse shortwave power.  This is a combination of the diffuse radiation from atmosphere (scattering, forward reflection, etc.), and reflection from ground.
    ! For both of these components, only half of the contribution is assumed to reach the container (the other half goes away from the container).
    ssd = (1.-refl)*(area_side/2.)*((sd*(1.-(shade/2.)))+(REF_G*st*(1.-shade)))
    IF (ssd < 0.) THEN
        ssd = 0.
    END IF

    ! Calculate the direct shortwave component.  This is just the total flux minus the diffuse component.
    sb = st-sd

    ! Calculate direct shortwave power.  This is the direct beam shortwave on the sides of the container.  Dependent on sun angle and shade.
    ssb = (1.-refl)*area_beam*sb*TAN(z)*(1.-shade)
    IF ((TAN(z) < 0.).OR.(ssb < 0.)) THEN
        ssb = 0.
    END IF

    ! Calculate the total sideways shortwave power.
    swside = ssb+ssd

    END SUBROUTINE calc_swside

    SUBROUTINE calc_lwdown (area_top, vp, ta, shade, clf, e_a, lwdown)
    ! Description: Calculates downward longwave radiative power.
    ! Input: Container area, vapor pressure, air temperature, cloud temperature, total cloud fraction, shade
    ! Output: Downward longwave radiative power
    IMPLICIT NONE

    ! Variable Declarations
    REAL, INTENT(IN) :: area_top                        ! Container area - Top
    REAL, INTENT(IN) :: vp                              ! Vapor pressure (kPa)
    REAL, INTENT(IN) :: ta                              ! Air temperature (K)
    REAL, INTENT(IN) :: shade                           ! Shade fraction
    REAL, INTENT(IN) :: clf                             ! Cloud fraction function
    REAL, INTENT(OUT) :: e_a                            ! Atmospheric emissivity
    REAL, INTENT(OUT) :: lwdown                         ! Downward longwave power (W)a
    REAL :: et                                          ! LW parameter

    ! Get effective emissivity of cloudless atmosphere, modified for clouds by Crawford and Duchon (1999)
    !e_a = (clf+(1.-clf)*(1.24*(((vp*10.)/ta)**(1./7.)))) ! Brutsaert (1975)
    et = 46.5*((vp*10.)/ta) ! Prata (1996)
    e_a = (clf+(1.-clf))*(1.-(1.+et)*EXP(-1.*(1.2+(3.0*et))**0.5)) ! Prata (1996)

    ! Get downward longwave power
    ! First term is for atmosphere above container, second term is for the overlying shading surface (modified by atmospheric absorption, and assume same temperature as air).
    lwdown = area_top*(((1.-shade)*(e_a*SB*(ta**4)))+(shade*E_B*SB*(ta**4)))

    END SUBROUTINE calc_lwdown

    SUBROUTINE calc_lwup (area_top, tw, lwup)
    ! Description: Calculates upward longwave radiative power.
    ! Input: Container area, water temperature
    ! Output: Upward longwave radiative power
    IMPLICIT NONE

    ! Variable Declarations
    REAL, INTENT(IN) :: area_top                        ! Container area - Top
    REAL, INTENT(IN) :: tw                              ! Water temperature (K)
    REAL, INTENT(OUT) :: lwup                           ! Upward longwave power (W)

    lwup = area_top*E_W*SB*(tw**4)

    END SUBROUTINE calc_lwup

    SUBROUTINE calc_lwsidein (area_side, tg, ta, e_a, lwsidein)
    ! Description: Calculates sideways inbound longwave radiative power.
    ! Input: Container area, ground temperature
    ! Output: Sideways inbound longwave radiative power.  Assume half from ground, and half from surrounding air.
    IMPLICIT NONE

    ! Variable Declarations
    REAL, INTENT(IN) :: area_side                       ! Container area - Side
    REAL, INTENT(IN) :: tg                              ! Ground temperature (K)
    REAL, INTENT(IN) :: ta                              ! Air temperature (K)
    REAL, INTENT(IN) :: e_a                             ! Atmospheric emissivity
    REAL, INTENT(OUT) :: lwsidein                       ! Sideways inbound longwave power (W)

    lwsidein = (area_side/2.)*((E_G*SB*(tg**4))+(e_a*SB*(ta**4)))

    END SUBROUTINE calc_lwsidein

    SUBROUTINE calc_lwsideout (area_side, tcon, lwsideout)
    ! Description: Calculates sideways outbound longwave radiative power.
    ! Input: Container area, container temperature
    ! Output: Sideways outbound longwave radiative power
    IMPLICIT NONE

    ! Variable Declarations
    REAL, INTENT(IN) :: area_side                       ! Container area - Side
    REAL, INTENT(IN) :: tcon                            ! Container temperature (K)
    REAL, INTENT(OUT) :: lwsideout                      ! Sideways outbound longwave power (W)

    lwsideout = area_side*E_C*SB*(tcon**4)

    END SUBROUTINE calc_lwsideout

    SUBROUTINE calc_htop (area_top, r, ta, tw, rheat, htop)
    ! Description: Calculates sensible heat power.
    ! Input: Container area, radius, air temperature, water temperature
    ! Output: Sensible heat power to/from top
    IMPLICIT NONE

    ! Variable Declarations
    REAL, INTENT(IN) :: area_top                        ! Container area - Top
    REAL, INTENT(IN) :: r                               ! Container radius (for use in diameter calc)
    REAL, INTENT(IN) :: ta                              ! Air temperature (K)
    REAL, INTENT(IN) :: tw                              ! Water temperature (K)
    REAL, INTENT(OUT) :: rheat                          ! Heat resistance (s m-1)
    REAL, INTENT(OUT) :: htop                           ! Sensible heat to/from top (W)
    REAL :: d                                           ! Container diameter (m)
    REAL :: dcm                                         ! Container diameter (cm)
    REAL :: gr                                          ! Grashof Number
    REAL :: nu                                          ! Nusselt Number

    ! Get container diameter in m and in cm
    d = 2.*r
    dcm = 100.*d

    ! Calculate Grashof Number.  See Monteith and Unsworth (2008).
    gr = 158.*(dcm**3)*ABS(tw-ta)

    ! Calculate Nusselt Number
    nu = B1*(gr**A1)

    ! Calculate heat resistance
    rheat = d/(nu*KA)

    ! Calculate sensible heat
    htop = (area_top*HC_A*(tw-ta))/rheat

    END SUBROUTINE calc_htop

    SUBROUTINE calc_hside (area_side, r, ta, tcon, hside)
    ! Description: Calculates sensible heat power.
    ! Input: Container area, container radius, container height, air temperature, container temperature
    ! Output: Sensible heat power to/from sides
    IMPLICIT NONE

    ! Variable Declarations
    REAL, INTENT(IN) :: area_side                       ! Container area - Side
    REAL, INTENT(IN) :: r                               ! Container radius (m)
    REAL, INTENT(IN) :: ta                              ! Air temperature (K)
    REAL, INTENT(IN) :: tcon                            ! Container temperature (K)
    REAL, INTENT(OUT) :: hside                          ! Sensible heat to/from side (W)
    REAL :: d                                           ! Container diameter (m)
    REAL :: dcm                                         ! Container diameter (cm)
    REAL :: gr                                          ! Grashof Number
    REAL :: nu                                          ! Nusselt Number
    REAL :: rheat                                       ! Heat resistance (s m-1)

    ! Get container diameter in m and in cm
    d = 2.*r
    dcm = 100.*d

    ! Calculate Grashof Number.  See Monteith and Unsworth (2008).
    gr = 158.*(dcm**3)*ABS(tcon-ta)

    ! Calculate Nusselt Number
    nu = B2*(gr**A2)

    ! Calculate heat resistance
    rheat = d/(nu*KA)

    ! Calculate sensible heat
    hside = (area_side*HC_A*(tcon-ta))/rheat

    END SUBROUTINE calc_hside

    SUBROUTINE calc_lh (area_top, vp, svp, rheat, lh)
    ! Description: Calculates latent heat power.
    ! Input: Container radius, vapor pressure, saturation vapor pressure
    ! Output: Latent heat from top
    IMPLICIT NONE

    ! Variable Declarations
    REAL, INTENT(IN) :: area_top                        ! Container area - Top
    REAL, INTENT(IN) :: vp                              ! Vapor pressure (kPa)
    REAL, INTENT(IN) :: svp                             ! Saturation vapor pressure (kPa)
    REAL, INTENT(IN) :: rheat                           ! Heat resistance from sensible heat flux calculation (s m-1)
    REAL, INTENT(OUT) :: lh                             ! Latent heat power (W)

    lh = (area_top*HC_A*(svp-vp))/(PSY*rheat)

    END SUBROUTINE calc_lh

    SUBROUTINE calc_ghc (area_bot, tg, tcon, k_c, thick, ghc)
    ! Description: Calculates ground heat flux.  This is by conduction between the container and the ground.
    ! Input: Container area, ground temperature, container temperature, container thermal conductivity, container thickness
    ! Output: Ground heat power to/from bottom of container
    IMPLICIT NONE

    ! Variable Declarations
    REAL, INTENT(IN) :: area_bot                        ! Container area - Bottom
    REAL, INTENT(IN) :: tg                              ! Ground temperature (K)
    REAL, INTENT(IN) :: tcon                            ! Container temperature (K)
    REAL, INTENT(IN) :: k_c                             ! Thermal conductivity of container (W m-1 K-1)a
    REAL, INTENT(IN) :: thick                           ! Thickness of container (m)
    REAL, INTENT(OUT) :: ghc                            ! Conduction (W)
    REAL :: dt                                          ! Temperature difference (K)
    REAL :: sum_dzk                                     ! Sum of dz/k terms

    ! Calculate the temperature difference between water and ground.
    dt = tcon-tg

    ! Summing together the depths / thermal conductivities for ground and container
    sum_dzk = (DZG/K_G)+(thick/k_c)

    ! Calculate conduction between ground and container
    ghc = area_bot*(dt/sum_dzk)

    END SUBROUTINE calc_ghc

    SUBROUTINE calc_ghw (area_bot, tcon, tw, wh, k_c, thick, ghw)
    ! Description: Calculates conduction between the container bottom and the water.
    ! Input: Container bottom area, container temperature, water temperature, water height, thermal conductivity of container, container thickness
    ! Output: Conduction to/from bottom of container and water
    IMPLICIT NONE

    ! Variable Declarations
    REAL, INTENT(IN) :: area_bot                        ! Container area - Bottom
    REAL, INTENT(IN) :: tcon                            ! Container temperature (K)
    REAL, INTENT(IN) :: tw                              ! Water temperature (K)
    REAL, INTENT(IN) :: wh                              ! Water height (m)
    REAL, INTENT(IN) :: k_c                             ! Thermal conductivity of container (W m-1 K-1)
    REAL, INTENT(IN) :: thick                           ! Thickness of container (m)
    REAL, INTENT(OUT) :: ghw                            ! Conduction (W)
    REAL :: dt                                          ! Temperature difference (K)
    REAL :: dzw                                         ! Distance to take the temperature gradient over - Water
    REAL :: sum_dzk                                     ! Sum of dz/k terms

    ! Calculate the temperature difference between water and ground.
    dt = tw-tcon

    ! Calculate the distance to take the temperature gradient over for water.
    dzw = wh/2.

    ! Summing together the depths / thermal conductivities for ground, container, and water
    sum_dzk = (dzw/K_W)+(thick/k_c)

    ! Calculate conduction between container and water (through bottom)
    ghw = area_bot*(dt/sum_dzk)

    END SUBROUTINE calc_ghw

    SUBROUTINE calc_cond (area_side, tcon, tw, k_c, thick, r, cond)
    ! Description: Calculates conduction between container side walls and water.
    ! Input: Container area, container temperature, water temperature, container thermal conductivity, container thickness, container radius
    ! Output: Conduction to/from side walls of container and water
    IMPLICIT NONE

    ! Variable Declarations
    REAL, INTENT(IN) :: area_side                       ! Container area - Sides (m)
    REAL, INTENT(IN) :: tcon                            ! Container temperature (K)
    REAL, INTENT(IN) :: tw                              ! Water temperature (K)
    REAL, INTENT(IN) :: k_c                             ! Thermal conductivity of container (W m-1 K-1)
    REAL, INTENT(IN) :: thick                           ! Thickness of container (m)
    REAL, INTENT(IN) :: r                               ! Container radius (m)
    REAL, INTENT(OUT) :: cond                           ! Conduction (W)
    REAL :: dt                                          ! Temperature difference (K)
    REAL :: sum_dzk                                     ! Sum of dz/k terms

    ! Calculate the temperature difference between water and ground.
    dt = tw-tcon

    ! Summing together the depths / thermal conductivities for ground, container, and water
    sum_dzk = (thick/k_c)+(r/K_W)

    ! Calculate conduction between container and water (through sides)
    cond = area_side*(dt/sum_dzk)

    END SUBROUTINE calc_cond

    SUBROUTINE calc_evap (area_top, lh, time_step, evap)
    ! Description: Calculates evaporation
    ! Input: Container radius, latent heat
    ! Output: Evaporation
    IMPLICIT NONE

    ! Variable Declarations
    REAL, INTENT(IN) :: area_top                        ! Container area - Top
    REAL, INTENT(IN) :: lh                              ! Latent heat (W)
    REAL, INTENT(IN) :: time_step                       ! Model time step
    REAL, INTENT(OUT) :: evap                           ! Evaporation (mm)

    evap = (lh*time_step*1000.)/(area_top*LV*RHO_W)

    END SUBROUTINE calc_evap

    SUBROUTINE update_wh (precip, area_top, area_bot, evap, h, shade, mfill, fill_intv, mempty, empty_intv, cnt, t_diff, precip_eff, wh)  ! ADDED mempty and empty_intv (2023-08-08, CNY)
    ! Description: Calculates updated water height in container
    ! Input: Precipitation (mm h-1), evaporation (mm h-1), container height (m), shade, water height (m)
    ! Output: Updated water height
    IMPLICIT NONE

    ! Variable Declarations
    REAL, INTENT(IN) :: precip                          ! Precipitation (mm)
    REAL, INTENT(IN) :: area_top                        ! Container area - Top
    REAL, INTENT(IN) :: area_bot                        ! Container area - Bottom
    REAL, INTENT(INOUT) :: evap                         ! Evaporation (mm)
    REAL, INTENT(IN) :: h                               ! Container height (m)
    REAL, INTENT(IN) :: shade                           ! Shade fraction
    LOGICAL, INTENT(IN) :: mfill                        ! Manual fill flag
    INTEGER, INTENT(IN) :: fill_intv                    ! Manual fill interval in hours
    LOGICAL, INTENT(IN) :: mempty                       ! Manual empty flag                (ADDED ON 2023-08-08, CNY)
    INTEGER, INTENT(IN) :: empty_intv                   ! Manual empty interval in hours   (ADDED ON 2023-08-08, CNY)
    INTEGER, INTENT(IN) :: cnt                          ! Running count of time periods done
    REAL, INTENT(IN) :: t_diff                          ! Ratio of input data interval to time step
    REAL, INTENT(OUT) :: precip_eff                     ! Effective precipitation (reaching container, mm)
    REAL, INTENT(INOUT) :: wh                           ! Water height (m)

    ! Modify precipitation based on shade.  Assume that rain only reaches container for fraction not shaded.
    precip_eff = precip*(1.-shade)*(area_top/area_bot)*(1./t_diff)

    ! Calculate water height
    wh = wh+((precip_eff-evap)/1000.)

    ! Manual empty  (ADDED ON 2023-08-08, CNY)
    ! Empty every empty_intv hours (EXCLUDING the very first timestep).
    IF (mempty .AND. (cnt-1) /= 0 .AND. (MODULO(cnt-1,empty_intv) == 0)) THEN
        wh = 0.
    END IF

    ! Manual fill
    ! Fill every fill_intv hours (INCLUDING the very first timestep).
    IF (mfill .AND. (MODULO(cnt-1,fill_intv) == 0)) THEN
        wh = h
    END IF

    ! Don't let the water height exceed the container height
    IF (wh > h) THEN
        wh = h
    END IF

    ! Don't let water height get below the specified water height floor
    IF (wh <= WHFLOOR) THEN
        wh = WHFLOOR
        evap = 0.
    END IF

    END SUBROUTINE update_wh

    SUBROUTINE update_tw (bal_w, area_top, wh, time_step, tw)
    ! Description: Calculates updated water temperature
    ! Input: Energy balance, container area, water height in container, time step water temperature
    ! Output: Updated water temperature
    IMPLICIT NONE

    ! Variable Declarations
    REAL, INTENT(IN) :: bal_w                           ! Energy balance (W)
    REAL, INTENT(IN) :: area_top                        ! Container area - Top
    REAL, INTENT(IN) :: wh                              ! Updated height of water in container (m)
    REAL, INTENT(IN) :: time_step                       ! Model time step
    REAL, INTENT(INOUT) :: tw                           ! Updated water temperature (K to C)
    REAL :: dt                                          ! Temperature change (K)

    dt = (bal_w*time_step)/(area_top*wh*HC_W)

    tw = tw+dt

    END SUBROUTINE update_tw

    SUBROUTINE update_tcon (bal_c, vol_cont, time_step, tcon)
    ! Description: Calculates updated container temperature
    ! Input: Energy balance, container area, time step, container temperature
    ! Output: Updated container temperature
    IMPLICIT NONE

    ! Variable Declarations
    REAL, INTENT(IN) :: bal_c                           ! Energy balance (W)
    REAL, INTENT(IN) :: vol_cont                        ! Volume of container itself
    REAL, INTENT(IN) :: time_step                       ! Model time step
    REAL, INTENT(INOUT) :: tcon                         ! Updated container temperature (K to C)
    REAL :: dt                                          ! Temperature change (K)

    dt = (bal_c*time_step)/(vol_cont*HC_C)

    tcon = tcon+dt

    END SUBROUTINE update_tcon

END MODULE container_mod
