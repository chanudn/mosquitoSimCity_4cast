PROGRAM container

! File: container_main.f90
! Description: Main program for modeling surface energy budget for water containers.
! Input: Module file for subroutines, T and RH data, rainfall, optional soil and cloud data, initialization parameters
! Output: Modeled water temperature and height in container
! Dan Steinhoff
! 13 March 2013

! Revision history:
!   2023-08-08, CNY: added manual emptying (mempty and empty_intv)
!   2023-08-24, CNY: added use of io_emsg to debug failure in opening input file
!	2023-09-11, CNY: edited to handle longer filenames (formerly up to 100 char, now up to 300 char)
!                     - edited variable declarations for data_file, outfile (CHARACTER(len=100) -> CHARACTER(len=300))
!                     - edited variable declaration for io_emsg (CHARACTER(len=100) -> CHARACTER(len=350))
!   2023-11-17, CNY: edited to calculate VPD and print it to output file
!   2023-11-30, CNY: changed clouds metadata in header of output file
!                     - formatting now matches that of other text in header, simplifying programmatic parsing of this text
!                     - changed text indicating that cloud data was provided by user (from "Provided" to "Observed")
!                     - added text indicating that cloud data was specified as constant day/night values ("Specified")

USE container_mod
IMPLICIT NONE

!!! Constant Declarations are in the module container_mod

!!! Variable Declarations
! Files
CHARACTER(len=100) :: param_file         ! Parameters file
CHARACTER(len=300) :: data_file          ! Data file name
CHARACTER(len=300) :: outfile            ! Output file name
! Location parameters
CHARACTER(len=30) :: city                ! City
CHARACTER(len=30) :: station             ! Station
REAL :: lat                              ! Latitude
REAL :: lon                              ! Longitude
REAL :: elev                             ! Elevation (m)
! Time parameters
REAL :: doy                              ! Day of year
REAL :: tod                              ! Time of day
CHARACTER(len=8) :: date                 ! Date being read (YYYYMMDDHH)
CHARACTER(len=4) :: year                 ! Year
CHARACTER(len=2) :: month                ! Month
CHARACTER(len=2) :: day                  ! Day
REAL :: time_step                        ! Time step of the model
REAL :: t_diff                           ! Ratio of input data interval to time step
! Container parameters
REAL :: rt1                              ! Container radius 1 (Top) (m)
REAL :: rt2                              ! Container radius 2 (Top) (m)
REAL :: rb1                              ! Container radius 1 (Body) (m)
REAL :: rb2                              ! Container radius 2 (Body) (m)
REAL :: h                                ! Container height (m)
REAL :: refl                             ! Container reflectivity
REAL :: wh                               ! Initial water height (m)
REAL :: tw                               ! Initial water temperature (C)
REAL :: tcon                             ! Container temperature (c)
REAL :: shade                            ! Shade fraction
REAL :: tw_init                          ! Initial water temperature (C)
REAL :: k_c                              ! Thermal conductivity of container (W m-1 K-1)
REAL :: thick                            ! Container thickness (m)
REAL :: area_top                         ! Container area - top (m2)
REAL :: area_bot                         ! Container area - bottom (m2)
REAL :: area_side                        ! Container area - side (m2)
REAL :: area_beam                        ! Container area - beam (m2)
REAL :: vol_cont                         ! Volume of container itself
INTEGER :: fill_intv                     ! Manual fill interval (hours)
INTEGER :: empty_intv                    ! Manual empty interval (hours)  (ADDED ON 2023-08-08, CNY)
CHARACTER(len=5) :: cshape               ! Shape of container (currently 'ROUND', 'RECT', or 'TIRE')
! Atm. variables
REAL :: ta                               ! Observed air temperature (C)
REAL :: rh                               ! Observed relative humidity (%)
REAL :: z                                ! Solar zenith angle (rad)
REAL :: tau                              ! Optical thickness
REAL :: tau_cs                           ! Clear-sky optical thickness
REAL :: cld                              ! Low-cloud fraction - day
REAL :: cmd                              ! Middle-cloud fraction - day
REAL :: chd                              ! High-cloud fraction - day
REAL :: cln                              ! Low-cloud fraction - night
REAL :: cmn                              ! Middle-cloud fraction - night
REAL :: chn                              ! High-cloud fraction - night
REAL :: cl                               ! Low-cloud fraction
REAL :: cm                               ! Middle-cloud fraction
REAL :: ch                               ! High-cloud fraction
REAL :: ct                               ! Total cloud fraction
REAL :: svp                              ! Saturation vapor pressure (kPa)
REAL :: vp                               ! Vapor pressure (kPa)
REAL :: vpd                              ! Vapor pressure deficit (kPa)   (ADDED ON 2023-11-17, CNY)
REAL :: rheat                            ! Heat resistance (s m-1)
REAL :: evap                             ! Evaporation (mm)
REAL :: precip                           ! Precipitation (mm)
REAL :: precip_eff                       ! Effective precipitation (mm)
REAL :: ta_avg                           ! Running average of air temperature
REAL :: e_a                              ! Atmospheric emissivity
REAL :: clf                              ! Cloud fraction function
REAL, DIMENSION(24) :: ta_array          ! Array of air temperature values
! Ground variables
REAL :: tg                               ! Ground temperature (C)
REAL :: tg_clim                          ! Ground temperature climatology (C)
! Energy Bal. variables
REAL :: st                               ! Total shortwave flux (W m-2)
REAL :: ssb                              ! Direct solar component (W)
REAL :: ssd                              ! Diffuse solar component (W)
REAL :: swdown                           ! Downward shortwave power (W)
REAL :: swside                           ! Sideways shortwave power (W)
REAL :: lwdown                           ! Downward longwave power (W)
REAL :: lwup                             ! Upward longwave power (W)
REAL :: lwsidein                         ! Sideways inbound longwave power (W)
REAL :: lwsideout                        ! Sideways outbound longwave power (W)
REAL :: htop                             ! Sensible heat to/from top (W)
REAL :: hside                            ! Sensible heat to/from sides (W)
REAL :: lh                               ! Latent heat to/from top (W)
REAL :: ghc                              ! Ground heat power  - Ground/Container (W)
REAL :: ghw                              ! Ground heat power  - Container/Water (W)
REAL :: cond                             ! Conduction between container side walls and water (W)
REAL :: bal_w                            ! Energy balance - Water (W)
REAL :: bal_c                            ! Energy balance - Container (W)
! Flags
LOGICAL :: cflag                         ! Whether cloud info provided in input file (TRUE yes, FALSE no)
LOGICAL :: sflag                         ! Whether soil temperature info provided in input file (TRUE yes, FALSE no)
LOGICAL :: mfill                         ! Manual fill flag, every week (TRUE yes, FALSE no)
LOGICAL :: mempty                        ! Manual empty flag (TRUE yes, FALSE no)  (ADDED ON 2023-08-08, CNY)
! Misc.
CHARACTER(len=200) :: header             ! File header
CHARACTER(len=350) :: io_emsg            ! Status message    ! ADDED for debugging (2023-08-24, CNY)
INTEGER :: iostat, status                ! Status variables
INTEGER :: cnt                           ! Counter of times read
INTEGER :: i                             ! Do-loop counter

!!! Get contents of parameter file
param_file = 'param_file.txt'
CALL read_params (TRIM(param_file), data_file, city, station, lat, lon, elev, cshape, rt1, rt2, rb1, rb2, h, refl, k_c, thick, wh, tw, tg_clim, time_step, t_diff, shade, cflag, sflag, mfill, fill_intv, mempty, empty_intv, outfile, cld, cmd, chd, cln, cmn, chn) ! ADDED mempty and empty_intv (2023-08-08, CNY)

!!! Open output file and print out header
OPEN (UNIT=200, FILE=TRIM(outfile), STATUS='UNKNOWN', ACTION='WRITE', POSITION='REWIND', IOSTAT=iostat)
IF (iostat /= 0) THEN
    WRITE (*,2000) iostat
    2000 FORMAT (1X, 'Problem opening output file: ', I6)
    STOP
END IF
2010 FORMAT ('Location: City: ', A, 2X, 'Station Name: ', A, 2X, 'Lat: ', F8.4, 2X, 'Lon: ', F8.4, 2X, 'Elev: ', F7.2, ' m')
WRITE (200,2010) TRIM(city), TRIM(station), lat, lon, elev
2020 FORMAT ('Container Specs: Shape: ', A5, 2X, 'Lengths (Top): ', F6.4, ',', F6.4, 2X, 'Lengths (Body): ', F6.4, ',', F6.4, 2X, 'Height: ', F6.4, 2X, 'Refl: ', F6.2, 2X, 'Conduc: ', F6.2, 2X, 'Thickness: ', F6.4)
WRITE (200,2020) cshape, rt1, rt2, rb1, rb2, h, refl, k_c, thick
2025 FORMAT ('     Initial Water Lev: ', F6.4, 2X, 'Initial Water Temp: ', F6.2, 2X, 'Manual Fill: ', L, 2X, 'Fill Intv: ', I5, 2X, 'Manual Empty: ', L, 2X, 'Empty Intv: ', I5) ! ADDED mempty and empty_intv (2023-08-08, CNY)
WRITE (200,2025) wh, tw, mfill, fill_intv, mempty, empty_intv ! ADDED mempty and empty_intv (2023-08-08, CNY)
IF (cflag) THEN
    2030 FORMAT ('Physical Specs: Ground T: ', F6.2, 2X, 'Clouds: Observed', 2X, 'Time Step: ', F6.2, 2X, 'Shade: ', F4.2)  ! CHANGED formatting of clouds text to be consistent with other text (makes programmatic parsing of this info easier) (2023-11-30, CNY)
WRITE (200,2030) tg_clim, time_step, shade
ELSE
    2035 FORMAT ('Physical Specs: Ground T: ', F6.2, 2X, 'Clouds: Specified', 2X, 'Clouds - Day (Low, Middle, High): ', F4.2, 2X, F4.2, 2X, F4.2, 2X, 'Clouds - Night (Low, Middle, High): ', F4.2, 2X, F4.2, 2X, F4.2, 2X, 'Time Step : ', F6.2, 2X, 'Shade: ', F4.2) ! CHANGED formatting of clouds text to be consistent with other text (makes programmatic parsing of this info easier) (2023-11-30, CNY)
    WRITE (200,2035) tg_clim, cld, cmd, chd, cln, cmn, chn, time_step, shade
END IF
2040 FORMAT ('-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------') ! ADDED more dashes (2023-11-17 and 2024-02-05, CNY)
WRITE (200,2040)
2050 FORMAT ('Year Month Day Hour   SWDOWN   SWSIDE   LWDOWN     LWUP     LWIN    LWOUT     SHUP    SHOUT     LHUP      GHW      GHC     COND    BAL_W    BAL_C     PRCP  PRCPEFF     EVAP       WH       TG       TA     TCON       TW      VPD')  ! ADDED VPD (2023-11-17, CNY) AND RAW PRECIP (2024-02-05, CNY)
WRITE (200,2050)
2060 FORMAT ('                         (W)      (W)      (W)      (W)      (W)      (W)      (W)      (W)      (W)      (W)      (W)      (W)      (W)      (W)     (mm)     (mm)     (mm)     (mm)      (C)      (C)      (C)      (C)    (kPa)')  ! ADDED (kPa) (2023-11-17, CNY) AND RAW PRECIP (2024-02-05, CNY)
WRITE (200,2060)

!!! Read the input data, and go through energy budget calculations at each time
! Open input data files: Must have T, RH, and precipitation.  Cloud fraction and ground temperature optional.
OPEN (UNIT=100, FILE=TRIM(data_file), STATUS='OLD', ACTION='READ', IOSTAT=iostat, IOMSG=io_emsg) ! ADDED ",IOMSG=io_emsg" for debugging (2023-08-24, CNY)

! If open successful, then continue
fileopen: IF (iostat == 0) THEN
    ! Read header
    READ (100,*) header

    ! Misc. parameter setup
    cnt = 1                  ! Initialize counter
    ta_avg = 0.              ! Initialize running temperature average
    tg_clim = tg_clim+273.15 ! Convert climatological tg to K
    tw = tw+273.15           ! Convert water temperature to K
    tcon = tw                ! Initialize container temperature to water temperature
    tw_init = tw             ! Set initial water temperature

    ! Loop through the data file
    readloop: DO
        ! Read a line of data
        IF (cflag) THEN ! We have cloud cover data
            IF (sflag) THEN ! We have soil temperature data
                READ (100,*,IOSTAT=status) year, month, day, tod, doy, ta, rh, precip, cl, cm, ch, tg
                ! Convert tg from C to K
                tg = tg+273.15
            ELSE ! We don't have soil temperature data
                READ (100,*,IOSTAT=status) year, month, day, tod, doy, ta, rh, precip, cl, cm, ch
            END IF
        ELSE ! We don't have cloud cover data
            IF (sflag) THEN ! We have soil temperature data
                READ (100,*,IOSTAT=status) year, month, day, tod, doy, ta, rh, precip, tg
                ! Convert tg from C to K
                tg = tg+273.15
            ELSE ! We don't have soil temperature data
                READ (100,*,IOSTAT=status) year, month, day, tod, doy, ta, rh, precip
            END IF
        END IF
        IF (status /= 0) EXIT

        !!! Get some preliminary parameters
        ! Convert air temperature from C to K
        ta = ta+273.15

        ! If cloud fraction data not provided (instead specified), determine which cloud fraction values to use - day or night
        IF (.NOT.cflag) THEN
            IF ((tod >= 13.).AND.(tod <= 16.)) THEN
                cl=cld
                cm=cmd
                ch=chd
            ELSE IF ((tod >= 21.).OR.(tod <= 8.)) THEN
                cl=cln
                cm=cmn
                ch=chn
            ELSE !IF (((tod >= 9).AND.(tod <= 12)).OR.((tod >= 17).AND.(tod <=20))) THEN
                cl=(cld+cln)/2.
                cm=(cmd+cmn)/2.
                ch=(chd+chn)/2.
            END IF
        END IF

        ! Get total cloud fraction.  Add up the three components, but keep at or below 1.
        ct = 0.5*(cl+cm+ch)
        IF (ct > 1.) THEN
            ct = 1.
        END IF

        ! Handle missing precip data
        IF (precip < 0.) THEN
            precip = 0.
        END IF

        ! Saturation vapor pressure
        CALL calc_svp (ta, svp)

        ! Vapor pressure
        CALL calc_vp (svp, rh, vp)

        ! Vapor pressure deficit  (ADDED ON 2023-11-17, CNY)
        CALL calc_vpd (svp, vp, vpd)

        ! Soil temperature, if we need it.  Initialized (first 24 hours) as averaged air temperature.  After that, follows a diurnal sinusoidal wave.
        IF (.NOT.sflag) THEN
            IF (cnt < 25) THEN
                CALL calc_tg_init (ta, ta_avg, tg_clim, cnt, ta_array, tg)
            ELSE
                CALL calc_tg_wave (ta, cnt, tod, time_step, ta_array, tg)
            END IF
        END IF

        ! Container surface areas: Top water surface, container sides, approximation of solar beam on side walls.
        CALL calc_area (rt1, rt2, rb1, rb2, h, thick, cshape, area_top, area_bot, area_side, area_beam, vol_cont)

        ! Loop through each time step, in between reading input data
        DO i = 1, INT(t_diff)

            ! If water temperature is marked as missing (for example, rainfall has been received to restore water height to minimum level),
            ! then set it back to the initial value.
            IF (tw < 0.) THEN
                tw = tw_init
                tcon = tw
            END IF

            ! Adjust tod for minutes past the hour
            tod = tod+(1./60.)

            ! Solar zenith angle
            CALL calc_z (lat, doy, tod, z)

            ! Optical thickness
            CALL calc_tau (cl, cm, ch, z, tau, tau_cs)

            ! Check the water level in the container
            IF (wh > WHMIN) THEN
                ! OK.  Proceed.

                !!! Calculate energy balance terms
                ! Shortwave down
                CALL calc_swdown (area_top, tau, tau_cs, z, doy, shade, clf, st, swdown)

                ! Shortwave side
                CALL calc_swside (refl, area_side, area_beam, z, tau, st, shade, ssd, ssb, swside)

                ! Longwave down
                CALL calc_lwdown (area_top, vp, ta, shade, clf, e_a, lwdown)

                ! Longwave up
                CALL calc_lwup (area_top, tw, lwup)

                ! Longwave side in
                CALL calc_lwsidein (area_side, tg, ta, e_a, lwsidein)

                ! Longwave side out
                CALL calc_lwsideout (area_side, tcon, lwsideout)

                ! Sensible heat top
                CALL calc_htop (area_top, rt1, ta, tw, rheat, htop)

                ! Sensible heat sides
                CALL calc_hside (area_side, rb1, ta, tcon, hside)

                ! Latent heat top
                CALL calc_lh (area_top, vp, svp, rheat, lh)

                ! Ground heat - Conduction from ground to/from container
                CALL calc_ghc (area_bot, tg, tcon, k_c, thick, ghc)

                ! Ground heat - Conduction from container to/from water
                CALL calc_ghw (area_bot, tcon, tw, wh, k_c, thick, ghw)

                ! Container - Water Conduction
                CALL calc_cond (area_side, tcon, tw, k_c, thick, rb1, cond)

                !!! Calculate energy balances
                bal_w = swdown+lwdown-lwup-htop-lh-ghw-cond
                bal_c = swside+lwsidein-lwsideout-hside-ghc+ghw+cond

                !! Update terms
                ! Calculate evaporation
                CALL calc_evap (area_top, lh, time_step, evap)

                ! Update container water volume
                CALL update_wh (precip, area_top, area_bot, evap, h, shade, mfill, fill_intv, mempty, empty_intv, cnt, t_diff, precip_eff, wh)  ! ADDED mempty and empty_intv (2023-08-08, CNY)

                ! Update water temperature
                CALL update_tw (bal_w, area_top, wh, time_step, tw)

                ! Update container temperature
                CALL update_tcon (bal_c, vol_cont, time_step, tcon)

            ELSE
                ! Water level too low.  Abort calculations.  Update water level, based on rainfall
                swdown=MISS
                ssb=MISS
                ssd=MISS
                swside=MISS
                lwdown=MISS
                lwup=MISS
                lwsidein=MISS
                lwsideout=MISS
                htop=MISS
                hside=MISS
                lh=MISS
                ghc=MISS
                ghw=MISS
                cond=miss
                bal_w=MISS
                bal_c=MISS
                evap=EVAP_LOW/t_diff
                tw=MISS+273.15
                tcon=MISS+273.15
                CALL update_wh (precip, area_top, area_bot, evap, h, shade, mfill, fill_intv, mempty, empty_intv, cnt, t_diff, precip_eff, wh)  ! ADDED mempty and empty_intv (2023-08-08, CNY)
            END IF

            IF (i == 1) THEN
                ! Print everything to output file
                2070 FORMAT (A4, 1X, A5, 1X, A3, 1X, I4, 1X, F8.2, 1X, F8.2, 1X, F8.2, 1X, F8.2, 1X, F8.2, 1X, F8.2, 1X, F8.2, 1X, F8.2, 1X, F8.2, 1X, F8.2, 1X, F8.2, 1X, F8.2, 1X, F8.2, 1X, F8.2, &
                             1X, F8.2, 1X, F8.2, 1X, F8.2, 1X, F8.2, 1X, F8.2, 1X, F8.2, 1X, F8.2, 1X, F8.2, 1X, F8.2) ! ADDED VPD (2023-11-17, CNY) AND RAW PRECIP (2024-02-05, CNY)
                WRITE (200,2070) year, month, day, INT(tod), swdown, swside, lwdown, lwup, lwsidein, lwsideout, htop, hside, lh, ghw, ghc, cond, bal_w, bal_c, precip, precip_eff*t_diff, &
                                 evap*t_diff, wh*1000., tg-273.15, ta-273.15, tcon-273.15, tw-273.15, vpd ! ADDED VPD (2023-11-17, CNY) AND RAW PRECIP (2024-02-05, CNY)
            END IF

        END DO

        ! Increment counter
        cnt=cnt+1

    END DO readloop

    ! Either exited from error or from finishing data file
    readif: IF (status > 0) THEN
        WRITE (*,*) 'An error occurred reading data'
    END IF readif

! If open was not successful...
ELSE fileopen
    WRITE (*,1010) iostat
    1010 FORMAT (1X, 'File open failed for values: ', I6)
    WRITE (*,1020) io_emsg                     ! ADDED for debugging (2023-08-24, CNY)
    1020 FORMAT (1X, 'IO_emsg: ', A150)        ! ADDED for debugging (2023-08-24, CNY)
END IF fileopen

! Close files
CLOSE (UNIT=100)
CLOSE (UNIT=200)

WRITE (*,*) 'Finished!'

END PROGRAM container
