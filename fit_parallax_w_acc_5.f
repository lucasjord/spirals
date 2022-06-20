      program fit_parallax_w_acc

c----------------------------
c     Fits parallax and proper motion to VLBA data for one source
c     allowing for accelerated motions

c     Parameters: par, x0, x_dot, x_acc, y0, y_dot, y_acc

c     Input files:
c        parallax.dat.             contains (t,x, errx, y, erry) data
c        fit_par_acc_control.inp   various controls and parameters

c     Version 3: Assumes RA controls parallax when assigning degrees
c     of freedom for X and Y sigma_pdf calculations

c     Version 4: Replaced approximate distance of the Earth from
c                the solar-sytem barycenter with the NOVAS.
c                Now allow input observations times as one of three types:
c                1) decimal years (use day_number/365 or /366 as needed),
c                2) Julian Day Number
c                3) Modified JDN

c     Version 5: fixed up units for issues going from years to days for
c                input data


      implicit real*8 (a-h,o-z)

      character*48  data_file, out_file
      real*8     dates(1000), x(1000), xerr(1000), y(1000), yerr(1000)
      real*8     times(1000), data(1000), model(1000), resids(1000), 
     +           res_err(1000), partls(1000,401)
      real*8     params(401), new_params(401), param_sigmas(401),
     +           scaled_sigmas(401)
      real*8     model_slope_rem
      character*16  parnames(401)
      integer    paramids(401)
      logical    print_cor, print_resids

      integer    pos_type(1000)
      common /model_info/ t_mid,pos_type,ra_rad,dec_rad

c     Set some maximum control values, switches, and logical unit numbers...
      idebug       = 0                       ! used by least-sq routine
      itermax      = 1
      max_num_data = 1000

      print_cor    = .false.
      print_resids = .false.

      pi  = 4.d0*atan(1.d0)

      lu_control   = 9
      lu_print     = 6
      lu_out       = 7
      lu_out_2     = 10
      lu_out_3     = 11
      lu_out_4     = 12
      lu_data      = 8

c     =============================================================
c     Get input control parameters...
      call control_parameters ( lu_control, lu_out,
     +                  ra, dec, t_0, err_min_ra, err_min_dec,
     +                  num_params, num_solved_for,
     +                  params, paramids, parnames )

      t_mid = t_0

c     Convert RA, Dec to radians
      ra_rad  = ra * pi/12.d0                      ! radians
      dec_rad = dec * pi/180.d0

c     =============================================================
c     Input data for 1 source ...  
      data_file = 'parallax.dat'      ! hardwired data file name
      call read_data( lu_data, data_file, lu_out,
     +                t_0, n_pts, dates, x, xerr, y, yerr )

c     Check if t_mid was specified in input file, if not use center time
      if ( t_mid .eq. 0.d0 ) t_mid = t_0
      write (lu_print,1000) t_mid
 1000 format(/' Center time for proper motion equation:',f17.5)

c     Transfer 2-D data to 1-D working arrays...
      i_d = 0
      do n = 1, n_pts

         i_d = i_d + 1
         data(i_d)    = x(n)
         times(i_d)   = dates(n)
         res_err(i_d)=sqrt(err_min_ra**2 + xerr(n)**2)
         pos_type(i_d)= 1               ! indicates East offset

         i_d = i_d + 1
         data(i_d)    = y(n)
         times(i_d)   = dates(n)
         res_err(i_d)=sqrt(err_min_dec**2 + yerr(n)**2)
         pos_type(i_d)= 2               ! indicates North offset

       enddo

       num_data = i_d

C       ------------------------------------------------------------
C                       Iterate Least-Squares Fitting

        sigsq     = 9999.d0
        num_deg_freedom = num_data - num_solved_for
        print_cor = .false.
        iterpass  = 0

        do while ( iterpass.le.itermax )

           if ( iterpass .ge. 1 ) then

                call calc_partials ( num_params, params, paramids,
     +                 num_data, times, 
     +                 partls )

                if ( iterpass.eq.itermax ) print_cor = .true.
                call least_squares_fit( idebug,print_cor,max_num_data,
     +                         num_params,num_data,
     +                         parnames,params,paramids,resids,
     +                         partls,res_err,new_params,param_sigmas )

                call update_params ( num_params, paramids, 
     +                               params, new_params )

           endif

                call calc_residuals ( params, num_data, times, data,
     +                 model, resids ) 

                if ( iterpass .gt. 0 ) then
                   call calc_sigsq ( num_data, num_solved_for,
     +                  resids, res_err, 
     +                  sigsq_old, sigsq, sigma_pdf)
                endif

                iterpass = iterpass + 1

        enddo

        write (lu_print,2500) sigma_pdf, num_deg_freedom
 2500   format (/' sigma_pdf =',f10.3,
     +           ' for',i4,' degrees of freedom')

c       Print out scaled uncertainties...
        write (lu_print,2398)
 2398   format(/' Uncertainties scaled by "sigma_pdf":')
        do n_p=1, num_params
           scaled_sigmas(n_p) = sigma_pdf * param_sigmas(n_p)
           write (lu_print,2399) parnames(n_p), scaled_sigmas(n_p)
 2399      format(3x,a16,f12.6)
        enddo

c       Print out distance from scaled parallax
        par     = params(1) 
        par_err = scaled_sigmas(1)
        d_par   = 1.d0 / par
        d_low   = 1.d0 / (par + par_err)
        d_high  = 1.d0 / (par - par_err)
        err_low = d_par  - d_low
        err_high= d_high - d_par
        write (lu_print,2400) d_par, err_high, err_low
 2400   format(/' Distance =',f10.4,'; +err=',f8.4,', -err=',f8.4,
     +          ' kpc')

C       ------------------------------------------------------------
C                               PRINTOUT
C       Printout residuals?
        if ( print_resids ) then
           write (lu_print,3000) 
3000       format (//'         Data       Model    Residual      ',
     +               ' Error     Res/Err'/)
           do n = 1, num_data
                 wtres = resids(n) / res_err(n)
                 write (lu_print,3100) data(n), model(n), resids(n), 
     +                                          res_err(n), wtres
3100             format (1x,5f12.4)
           enddo
        endif

C       ------------------------------------------------------------
C                               PLOT FILES

        i0 = 1           ! params index before first motion parameter

c       Make ascii file containing the velocity, data, model, and resids 
c       over the velocity range used.
        do i_pos = 1, 2        ! separate files for RA an DEC models

         if ( i_pos .eq. 1 ) then
            write (out_file,4050) 
 4050       format('par_fit_results_ra.dat')
         endif
         if ( i_pos .eq. 2 ) then
            write (out_file,4060)
 4060       format('par_fit_results_dec.dat')
         endif

         open (unit=lu_out, file=out_file, form='formatted')

         write (lu_out,4100) 
 4100    format( '!    Times     Data      Model     Resids  ScaledErr')
         do n = 1, num_data, 2
            nn = (n-1) + i_pos
            scaled_err = res_err(nn) * sigma_pdf

            xJDN_1980 = 2444239.5d0
            time_yrs = 1980.d0 + (times(nn) - xJDN_1980)/365.2422d0
           
            write (lu_out,4200) time_yrs,data(nn),model(nn),
     +                          resids(nn),scaled_err
 4200       format(1x,f10.3,4f10.4)
         enddo

         close ( unit=lu_out )

        enddo

c       Make another output file of the model evaluated on a finer
c       grid of velocities for smoother plotting...
        do i_pos = 1, 2        ! separate files for RA an DEC models  

           if ( i_pos .eq. 1 ) then
                 out_file = 'par_fit_model_ra.dat'
              else
                 out_file = 'par_fit_model_dec.dat'
           endif

           open (unit=lu_out_2, file=out_file, form='formatted')

           n_mod_pts = 101
           t_span = times(num_data) - times(1)
           dt     = t_span / 100.d0
           do i = 1, n_mod_pts
              data_time = times(1) + (i-1)*dt
              call calc_model ( params, data_time, i_pos,
     +                          position ) 

              xJDN_1980 = 2444239.5d0
              time_yrs = 1980.d0 + (data_time-xJDN_1980)/365.2422d0

              write (lu_out_2,4200) time_yrs, position
	   enddo

           close ( unit=lu_out_2 )

        enddo

c       -------------------
c       Make even more output files of data and model with slope removed...
        do i_pos = 1, 2         ! separate files for RA an DEC models

           if ( i_pos .eq. 1 ) then
                 out_file = 'par_fit_results_desloped_ra.dat'
              else
                 out_file = 'par_fit_results_desloped_dec.dat'
           endif

           open (unit=lu_out_3, file=out_file, form='formatted')

           write (lu_out_3,4101) 
 4101      format('!   Velocity   Data      Model     Resids',
     +            '  ScaledErr   with baseline removed.')

           i0 = 1
           do n = 1, num_data, 2
              nn = (n-1) + i_pos
              del_t = times(nn) - t_mid           ! days
              del_t_yrs = del_t/365.2422d0        ! years
c             slope contribution...
c             (NB: slope pivots about t_mid)
              if ( i_pos .eq. 1 ) then
c                 RA data point...
                  x_off   = params( i0 + 1 )
                  x_slope = params( i0 + 2 )
                  x_acc   = params( i0 + 3 )
                  pm = x_off + x_slope*del_t_yrs 
     +                       + 0.5d0*x_acc*del_t_yrs**2
                else
c                 Dec data point...
                  y_off   = params( i0 + 4 )
                  y_slope = params( i0 + 5 )
                  y_acc   = params( i0 + 6 )
                  pm = y_off + y_slope*del_t_yrs 
     +                       + 0.5d0*y_acc*del_t_yrs**2
              endif

              data_slope_rem = data(nn) - pm
              model_slope_rem= model(nn)- pm
              scaled_err = res_err(nn) * sigma_pdf

              xJDN_1980 = 2444239.5d0
              time_yrs = 1980.d0 + (times(nn)-xJDN_1980)/365.2422d0

              write (lu_out_3,4200) time_yrs, data_slope_rem, 
     +                            model_slope_rem, resids(nn),
     +                            scaled_err
           enddo

           close ( unit=lu_out_3 )

        enddo

c       ---------------------------
        do i_pos = 1, 2           ! separate files for RA an DEC models
           if ( i_pos .eq. 1 ) then
                 out_file = 'par_fit_model_desloped_ra.dat'
              else
                 out_file = 'par_fit_model_desloped_dec.dat'
           endif

           open (unit=lu_out_4, file=out_file, form='formatted')
         
           i0 = 1
           n_mod_pts = 101
           do i = 1, n_mod_pts

              data_time = times(1) + (i-1)*dt
              del_t     = data_time - t_mid
              del_t_yrs = del_t/365.2422d0
              call calc_model ( params, data_time, i_pos,
     +                          position ) 

              if ( i_pos .eq. 1 ) then
c                 RA data point...
                  x_off   = params( i0 + 1 )
                  x_slope = params( i0 + 2 )
                  x_acc   = params( i0 + 3 )
                  pm = x_off + x_slope*del_t_yrs 
     +                       + 0.5d0*x_acc*del_t_yrs**2
                else
c                 Dec data point...
                  y_off   = params( i0 + 4 )
                  y_slope = params( i0 + 5 )
                  y_acc   = params( i0 + 6 )
                  pm = y_off + y_slope*del_t_yrs 
     +                       + 0.5d0*y_acc*del_t_yrs**2
              endif

              xJDN_1980 = 2444239.5d0
              time_yrs = 1980.d0 + (data_time-xJDN_1980)/365.2422d0

              position_pm_removed = position - pm
              write (lu_out_4,4200) time_yrs, position_pm_removed

	   enddo
           close ( unit=lu_out_4 )
        enddo

      stop
      end
c
c======================================================================
	subroutine open_ascii_file ( lu_in, ascii_file, lu_out )

C	Opens ascii file and reads and prints comments (if first character is "!") 
C       Leaves file open, ready to read data.

	character*48       ascii_file

	character*90       comment
	character*1        c_1
	equivalence       (c_1, comment)

        logical            first_comment

C	Open input ascii file...

        first_comment = .true.

	write (6,1000) ascii_file
 1000	format(/' Openning input ascii file "',a32,'"')

	open (unit=lu_in, file=ascii_file, status='old')

C       Read all comment lines (but stop at data)
	write (6,1010) 
 1010	format(' Comment lines follow:')

	c_1 = '!'
	do while ( c_1 .eq. '!' )

	   read  (lu_in,1020) comment
 1020	   format(a90)

	       if ( c_1 .eq. '!' ) then
		   write (6,1030) comment
 1030		   format(1x,a90)
		 else
		   backspace (unit=lu_in)
	       endif

	enddo

	write (6,1040)                ! print a blank line at end
 1040	format(' End of Comment lines',/1x)

	return
	end
c
c===================================================================
      subroutine read_data ( lu_data, data_file, lu_out, 
     +                       t_0, n_pts, dates, x, xerr, y, yerr )

      implicit real*8 (a-h,o-z)

      character*48  data_file
      character*1   c_1
      
      real*8   dates(1000), x(1000), xerr(1000), y(1000), yerr(1000)

      call open_ascii_file ( lu_data, data_file, lu_out )

      i    = 0
      ieof = 0

      t_0  = 0.d0

C     Now read in parallax data...
      do while ( ieof .eq. 0 )

c        Check for flagged line(s)...
         c_1 = '!'
         do while ( ieof.eq.0 .and. c_1.eq.'!' )

            read(lu_data,1000,iostat=ieof) c_1
 1000       format(a1)

         enddo
         if ( ieof.eq.0 ) then 

            backspace (unit=lu_data)

c           Now read data line...
            read(lu_data,*) date, offx, offxerr, offy, offyerr

            i = i + 1
            dates(i) = date            ! years
            x(i)     = offx            ! mas
            xerr(i)  = offxerr
            y(i)     = offy            ! mas
            yerr(i)  = offyerr

            write (6,2000) date, offx, offxerr, offy, offyerr 
 2000       format(1x,f15.5,3x,4f10.3)

            t_0 = t_0 + date

         endif

      enddo

      n_pts = i
      print*,' Number of un-flagged data epochs read in: ', n_pts

      if ( n_pts .gt. 0 ) then
         t_0 = t_0 / n_pts       ! unweighted center time
      endif

      close ( unit=lu_data )

      return
      end
c
c===================================================================
        subroutine control_parameters ( lu_control, lu_out, 
     +                  ra, dec, t_0, err_min_ra, err_min_dec,
     +                  num_params, num_solved_for,
     +                  params, paramids, parnames )


        implicit real*8 (a-h,o-z)

        real*8          params(401), x_paramids(401)
        character*16    parnames(401)
        integer         paramids(401)
        character*48    control_file

        control_file = 'fit_par_acc_control.inp'   ! hardwired input control file name

        call open_ascii_file ( lu_control, control_file, lu_out )

        read (lu_control,8010) ra                  ! R.A. (decimal hours)
        read (lu_control,8010) dec                 ! Dec  (decimal degrees)
        read (lu_control,8010) t_0                 ! mid-time (years)
        read (lu_control,8010) err_min_ra          ! data error floor for Right Ascension (mas)
	read (lu_control,8010) err_min_dec         ! data error floor for Declination (mas)

        read (lu_control,8010) params(1), x_paramids(1) ! parallax

        i0 = 1                       ! starting index for motion parameters
        read (lu_control,8010) params(i0+1), x_paramids(i0+1) ! X offset 
        read (lu_control,8010) params(i0+2), x_paramids(i0+2) ! X slope
        read (lu_control,8010) params(i0+3), x_paramids(i0+3) ! X acc
        read (lu_control,8010) params(i0+4), x_paramids(i0+4) ! Y offset
        read (lu_control,8010) params(i0+5), x_paramids(i0+5) ! Y slope
        read (lu_control,8010) params(i0+6), x_paramids(i0+6) ! Y acc

8010    format(2f10.3)

        write (6,8010) ra                       ! RA (hr)
        write (6,8010) dec                      ! Dec (deg)
        write (6,8010) t_0                      ! t-mid (years)
        write (6,8010) err_min_ra               ! data error floor (mas)
	write (6,8010) err_min_dec              ! data error floor (mas)

        write (6,8010) params(1), x_paramids(1) ! parallax

        i0 = 1                       ! starting index for motion parameters
        write(6,8010) params(i0+1), x_paramids(i0+1) ! X offset 
        write(6,8010) params(i0+2), x_paramids(i0+2) ! X slope
        write(6,8010) params(i0+3), x_paramids(i0+3) ! X acc
        write(6,8010) params(i0+4), x_paramids(i0+4) ! Y off
        write(6,8010) params(i0+5), x_paramids(i0+5) ! Y off
        write(6,8010) params(i0+6), x_paramids(i0+6) ! Y off

C       Set solve-for codes...
        num_params = 7                          ! number of parms
        num_solved_for = 0
        do n = 1, num_params
           paramids(n) = int ( x_paramids(n) )
           if (  paramids(n) .eq. 1 ) then
                 num_solved_for = num_solved_for + 1
           endif
        enddo

C       Enter parameter names...
        parnames(1) = 'Parallax  (mas) '
        parnames(2) = 'X off     (mas) '
        parnames(3) = 'X slope (mas/y) '
        parnames(4) = 'X acc (mas/y^2) '
        parnames(5) = 'Y off     (mas) '
        parnames(6) = 'Y slope (mas/y) '
        parnames(7) = 'Y acc (mas/y^2) '

        close(unit=lu_control)

        return
        end
c
	subroutine calc_residuals ( params, num_data, times, data,
     +                 model, resids )
 
C	Calculates models and resduals for entire data set...

	implicit real*8 ( a-h,o-z )

	real*8	params(401)

	real*8	times(1000), data(1000), model(1000), resids(1000)

        integer    pos_type(1000)
        common /model_info/ t_mid,pos_type,ra_rad,dec_rad

        do n = 1, num_data

           data_time = times(n)
           i_pos     = pos_type(n)
           call calc_model ( params, data_time, i_pos, 
     +                       position ) 

	   model(n) = position
	   resids(n) = data(n) - model(n)

	enddo
	
	return
	end
c
	subroutine calc_model ( params, data_time, i_pos, 
     +                          position ) 

c       Calculates a position (either RA or Dec given the data_time
c       and the model parameters of the parallax and 2 accelerating lines
c       (one each for the eastward and northward motions). 

	implicit real*8 ( a-h, o-z )

        real*8     params(401)

        integer    pos_type(1000)
        common /model_info/ t_mid,pos_type,ra_rad,dec_rad


c       Slope contribution...
c       (NB: slope pivots about t_mid)
        del_t = data_time - t_mid                      ! days
        del_t_yrs = del_t /365.2422d0                  ! years

c       Get proper motion parameters for the correct source/chan
        i0 = 1         ! params index before 1'st motion param

        if ( i_pos .eq. 1 ) then
c             RAcos(Dec) data... 
              x_off   = params( i0 + 1 )
              x_slope = params( i0 + 2 )
              x_acc   = params( i0 + 3 )
              pm = x_off + x_slope*del_t_yrs  
     +                   + 0.5d0*x_acc*del_t_yrs**2
           else
c             Dec data...
              y_off   = params( i0 + 4 )
              y_slope = params( i0 + 5 )
              y_acc   = params( i0 + 6 )
              pm = y_off + y_slope*del_t_yrs 
     +                   + 0.5d0*y_acc*del_t_yrs**2
        endif

c       Parallax offsets...
        call calc_parallax ( params, data_time, i_pos,
     +                       parallax_offset )

c       Sum contributions...
        position = pm + parallax_offset

	return 
	end
c
	subroutine calc_parallax ( params, t_yr, i_pos,
     +                       parallax_offset )

c	Calculates position offset expected from the parallax
c       of a source

        implicit real*8 (a-h, o-z)

        real*8     params(401), pos(3), vel(3)

        integer    pos_type(1000)
        common /model_info/ t_mid,pos_type,ra_rad,dec_rad

        pi    = 4.d0*atan(1.d0)
        twopi = 2.d0 * pi
        deg_to_rad = pi / 180.d0

        parallax = params(1)                    ! mas

        cos_ra = cos(ra_rad)
        sin_ra = sin(ra_rad)
        cos_dec= cos(dec_rad)
        sin_dec= sin(dec_rad)

c        Code used prior to Version 4
c        Longitude of sun is 0 at vernal equinox (~0.22 yrs)
c        sun_long = twopi * ( t_yr - 0.22d0 )    ! radians
c        cos_sun = cos( sun_long )
c        sin_sun = sin( sun_long )
c        obliquity = 23.4d0 * deg_to_rad         ! radians
c        cos_obl = cos( obliquity )
c        sin_obl = sin( obliquity )
c        Sun's position in cartesian coordinates in AU units
c        for a circular orbit
c        X = cos_sun
c        Y = sin_sun * cos_obl
c        Z = sin_sun * sin_obl
c        Correct for eccentricity of the Earth's orbit
c        earth_long = twopi * ( t_yr - 0.257d0 )    ! radians
c        factor = 1.0d0 + 0.0167d0 * sin( earth_long )

c       Get Julian Day Number (xJDN)...
        xJDN = 0.d0
        if ( t_yr.ge.1980.d0  .and. t_yr.lt.2050.d0 ) then
c          Have decimal year...convert to JDN
           call year_to_JDN(t_yr,xJDN)
        endif
        if ( t_yr.ge.44239.d0 .and. t_yr.lt.69807.d0) then
c          Have Modified Julian Day (MJD)...convert to JDN
           xJDN = 2400000.5d0 + t_yr
        endif
        if ( t_yr.ge.2444239.d0.and.t_yr.lt.2469807.d0) then
c          Have Julian Day (JDN) already...
           xJDN = t_yr
        endif
        if ( xJDN .eq. 0 ) then
            print *,' Invalid observation time; STOP'
            stop
        endif

c       Version 4 replaced above with the following from NOVAS routine
        iearth = 2        ! Earth's offset flag
        ik     = 1        ! Center of Sun (not barycenter) flag

        ik     = 0    !*************testing barycenter***********************

        call solsys (xJDN, iearth, ik, pos, vel, ierr)
        if ( ierr .ne. 0 ) then
            print *,' SOLSYS returned IERR.ne.0; STOP'
            stop
        endif

        X = -pos(1)                                       ! AU
        Y = -pos(2)
        Z = -pos(3)
        factor = 1.d0

        proj_east = ( Y*cos_ra - X*sin_ra )
        offset_east  = factor * parallax * proj_east     ! mas
 
        proj_north = Z*cos_dec - X*cos_ra*sin_dec - Y*sin_ra*sin_dec
        offset_north = factor * parallax * proj_north    ! mas

        if ( i_pos .eq. 1 ) then
c             RA component
              parallax_offset = offset_east
           else
c             Dec component
              parallax_offset = offset_north
        endif
 
	return
	end
c
        subroutine year_to_JDN ( t_yr, xJDN )

c       Given decimal year (eg, t_yr=2013.75) calculate Julian Day Number (xJDN)
c       Assumes t_yr = year + day_number/365 (or /366) for leap year

c       Works for 1980 and following, but note that it assumes a leap year 
c       every 4 years (so it will ultimately be wrong!

	implicit real*8 ( a-h,o-z )

        i_yr = int(t_yr)
        leap_year = mod(i_yr,4)             ! 0 on leap years; else 1, 2 or 3 
        num_years = i_yr - 1980             ! years since 1980
        num_leap_years = (num_years+3)/4    ! number of leap year days since 2000

        xJDN_1980 = 2444239.5d0             ! Jan 1, 1980 at 0hr UT

        xJDN_yr_jan1 = xJDN_1980 + num_years*365 + num_leap_years 
        
        if ( leap_year .ne. 0 ) then
              xJDN = xJDN_yr_jan1 + (t_yr-i_yr)*365.d0
           else
              xJDN = xJDN_yr_jan1 + (t_yr-i_yr)*366.d0
        endif

        return
        end
c
	subroutine calc_partials ( num_params, params, paramids,
     +                 num_data, times, 
     +                 partls )

	implicit real*8 ( a-h,o-z )

	real*8     params(401), params_wiggled(401)
	real*8	   partls(1000,401), times(1000)

	integer	   paramids(401)

        integer    pos_type(1000)
        common /model_info/ t_mid,pos_type,ra_rad,dec_rad

c	Set secondary parameter array and zero partials array...
	do i_p = 1, num_params

		params_wiggled(i_p) = params(i_p)

		do i_d = 1, num_data
			partls(i_d,i_p) = 0.d0
		enddo

	enddo

c       Calculate numerical partials...
	do i_d = 1, num_data

            data_time = times(i_d)
            i_pos     = pos_type(i_d)
            call calc_model ( params, data_time, i_pos, 
     +                        position )

            do i_p = 1, num_params

c               Calc partials for ID=1 or 2...
		if ( paramids(i_p) .ge. 1 ) then

c		   Change parameters by 1 part in 10^4; but never
c			by less than a value of 10^-6 (Jy or kms)...

                   del_param = abs( params(i_p) ) * 1.d-04
                   if ( del_param .lt. 1.d-06 ) del_param = 1.d-06

                   params_wiggled(i_p) = params(i_p) + del_param
			
                   call calc_model(params_wiggled, data_time, i_pos, 
     +                             position_wiggled)

                   partls(i_d,i_p)  = (position_wiggled - position) /
     +                                 del_param

                   params_wiggled(i_p) = params(i_p) 

                endif

             enddo

	enddo

	return
	end
c
	subroutine update_params(num_params, paramids,
     +                           params, new_params)

c       Takes parameter values before (params) and after (new_params)
c       least squares fit; adjusts the changes by the "gain";
c       and replaces the adjusted new parameters in the "params" array.

	implicit real*8 (a-h, o-z)

	real*8		params(401), new_params(401)
 
        integer         paramids(401)

        gain = 1.0                     ! hardwired gain parameter

	do j = 1, num_params

		delta_param = new_params(j) - params(j) 
		trial_param = params(j) + gain * delta_param
		params(j) = trial_param

	enddo

	return
	end
c
	subroutine calc_sigsq (num_data, num_solved_for,
     +			resids, res_err, 
     +			sigsq_old, sigsq, sigma_pdf)

C	Calculate new "sigsq"...

c       In version 3, changed d.o.f. calculation to assume that
c       RA controls parallax.

	implicit real*8	(a-h,o-z)

	real*8		resids(1000), res_err(1000)


c       First only RA data
	sigsq_old = sigsq
	sigsq     = 0.d0

	do i_d = 1, num_data, 2
		sigsq = sigsq + ( resids(i_d)/res_err(i_d) )**2
	enddo

	num_deg_freedom = num_data/2 - (num_solved_for+1)/2

	if ( num_deg_freedom .gt. 1 ) then
              sigsq_pdf = sigsq / float(num_deg_freedom) 
           else
              sigsq_pdf = 0.d0
        endif

	sigma_pdf = sqrt( sigsq_pdf )

        print *,' '
        print *,'sigma_pdf_RA = ',sigma_pdf

c       Next only Dec data
	sigsq_old = sigsq
	sigsq     = 0.d0

	do i_d = 2, num_data, 2
		sigsq = sigsq + ( resids(i_d)/res_err(i_d) )**2
	enddo

	num_deg_freedom = num_data/2 - (num_solved_for-1)/2

	if ( num_deg_freedom .gt. 1 ) then
              sigsq_pdf = sigsq / float(num_deg_freedom) 
           else
              sigsq_pdf = 0.d0
        endif

	sigma_pdf = sqrt( sigsq_pdf )

        print *,'sigma_pdf_Dec= ',sigma_pdf

c       Now all data data together
	sigsq_old = sigsq
	sigsq     = 0.d0

	do i_d = 1, num_data
		sigsq = sigsq + ( resids(i_d)/res_err(i_d) )**2
	enddo

	num_deg_freedom = num_data - num_solved_for

	if ( num_deg_freedom .ge. 1 ) then
              sigsq_pdf = sigsq / float(num_deg_freedom) 
           else
              sigsq_pdf = 0.d0
        endif

	sigma_pdf = sqrt( sigsq_pdf )

	return
	end
c
      SUBROUTINE LEAST_SQUARES_FIT ( IDEBUG,PRINT_COR,MAXDAT,NUMXES,
     +				N,VNAME,X,ID,S,C,E,XHAT2,ERROR )
C     SUBROUTINE FOR LINEAR LEAST SQUARE FITTING...
C     IT REQUIRES INPUT OF THE FOLLOWING:
C        IDEBUG    = DEBUGGING PARAMETER
C                    0 = SUPRESS PRINT OUT BEYOND SOLUTION AND CORRELATI
C                    1 = PRINT OUT ALL DATA AND MATRICES
C	 PRINT_COR = LOGICAL FLAG: PRINT CORRELATION MATRIX
C	 MAXDAT    = MAX NUM DATA POINTS.  DIMENSION LIMITS IN MAIN PROG.
C        NUMXES    = NUMBER OF PARAMETERS IN MODEL (WHETHER SOLVED FOR O
C        N         = NUMBER OF EQUATIONS
C        VNAME     = ALPHANUMERIC 'NAME' FOR THE PARAMETERS
C        X         = 1-D ARRAY OF INITIAL PARAMETER VALUES
C        ID        = 1-D ARRAY OF 'SOLVE FOR' CODES
C                  0 = HOLD PARAMETER CONSTANT
C                  1 = SOLVE FOR PARAMETER
C        S         = 1-D ARRAY OF DATA (E.G. VLB DIFFERENCED PHASES)
C        C         = 2-D ARRAY OF PARTIAL DERIVATIVES OF MODEL
C        E         = 1-D ARRAY OF ERRORS ON DATA POINTS
C
C     THE OUTPUT INCLUDES:
C        XHAT2     = 1-D ARRAY OF FINAL LEAST SQUARE VALUES FOR ALL PARA
C        THE SOLUTION AND CORRELATIONS ARE AUTOMATICALLY PRINTED OUT
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION C(MAXDAT,1)
      DIMENSION S(1),E(1),X(1),ID(1),XHAT2(1)

      character*16      VNAME(1)

        REAL*8          B(401,401),BSTORE(401,401),XIDENT(401,401),
     +                  COR(401,401),ERROR(401),STDEV(401),PR(401),
     +                  XHAT(401),SNEW(1000),PRODCT(1000)

        DIMENSION LL(401),MM(401)
        DATA      MSIZE /401/

C
      LOGICAL PRINT_COR
      character*8      quest(2)
      DATA QUEST/'YES     ','NO      '/

 9999 FORMAT (1X,10(1PD13.5))
 9998 FORMAT (1X,20F6.2)
 9417 FORMAT (/1X)

      IF (IDEBUG.EQ.0) GO TO 502
C
      WRITE (6,5522)
 5522 FORMAT (////50X,' C MATRIX')
      DO 500 I=1,N
      WRITE (6,9999) (C(I,J),J=1,NUMXES)
  500 CONTINUE
      WRITE (6,5523)
 5523 FORMAT (////40X,'        S(I)           E(I)')
      DO 501 I=1,N
      WRITE (6,9999) S(I),E(I)
  501 CONTINUE
  502 CONTINUE
      M=0
      do i = 1, numxes
         if ( id(i) .eq. 1 ) m = m + 1
      enddo
C        M IS THE NUMBER OF 'SOLVE FOR' PARAMETERS
      IF ( M.GT.N )  GO TO 998
C
      JNEW=0
      DO 12  J=1,NUMXES
	      IF (ID(J).ne.1) GO TO 12
	      JNEW=JNEW+1
	      DO 11 I=1,N
   11		 C(I,JNEW)=C(I,J)
   12 CONTINUE
C
C        WEIGHT EQUATIONS BY DIVIDING EACH BY THE ERROR-E(I)-(SEE SECT 2
C        OF TEXT)
      DO 13 I=1,N
	      SNEW(I)=S(I)/E(I)
	      DO 13 J=1,M
   13		 C(I,J)=C(I,J)/E(I)
      IF (IDEBUG.EQ.0) GO TO 21
      WRITE (6,2006)
 2006 FORMAT ('1'/////50X,' CSTAR MATRIX')
      DO 14 I=1,N
   14 WRITE (6,9999) (C(I,J),J=1,M)
      WRITE (6,2009)
 2009 FORMAT (/////50X,'     SSTAR')
      DO 15 I=1,N
   15 WRITE (6,9999) SNEW(I)
   21 CONTINUE
      ITEST2=0
      BMIN10=0.D0
      BMAX10=0.D0
      DO 22 I=1,M
      DO 22 J=1,M
   22 	B(I,J)=0.D0
      DO 24 I=1,M
      DO 24 J=1,M
	      DO 23 L=1,N
   23		 B(I,J)=B(I,J) + C(L,I)*C(L,J)
      IF (B(I,J).EQ.0.D0) GO TO 24
	      B1=DLOG10( DABS( B(I,J) ) )
	      IF (BMIN10.GT.B1) BMIN10=B1
	      IF (BMAX10.LT.B1) BMAX10=B1
   24 BSTORE(I,J)=B(I,J)
      IF (IDEBUG.EQ.0) GO TO 761
      WRITE (6,2010)
 2010 FORMAT ('1'/////50X,' C*TRANSPOSE C')
      DO 25 I=1,M
   25 WRITE (6,9999) (B(I,J),J=1,M)
  761 CONTINUE
C
C        THE SUBROUTINE 'MINV' IS A DOUBLE PRECISION MATRIX INVERSION
C        IT INVERTS THE MATRIX 'BB' AND STORES IT AS 'BB' (I.E. THE ORIGIN
C        MATIRX IS DESTROYED
C        *******************************
C
C
      IF (DABS(BMAX10-BMIN10).LT.18.D0) GO TO 962
      WRITE (6,9622) BMIN10,BMAX10
 9622 FORMAT (//1X,'********************     BMIN10 = ',F6.1,'
     1 BMAX10 = ',F6.1,'     ********************',/////1X)
  962 BADJST=10.D0**( (BMAX10+BMIN10)/2.D0 )
      DO 963 I=1,M
      DO 963 J=1,M
  963	 B(I,J)=B(I,J)/BADJST
      CALL MINV(B,M,DETERM,LL,MM,MSIZE)
      DO 964 I=1,M
      DO 964 J=1,M
  964	 B(I,J)=B(I,J)/BADJST
C
C
C        *******************************
C
      IF (IDEBUG.EQ.0) GO TO 226
      WRITE (6,2011) DETERM
 2011 FORMAT ('1'////'  THE DETERMINANT IS',1PD13.5)
      WRITE (6,2022)
 2022  FORMAT (////45X,' (C*TRANSPOSE C*) INVERSE')
      DO 26 I=1,M
   26 WRITE (6,9999) (B(I,J),J=1,M)
  226 CONTINUE
      DO 27 I=1,M
      DO 27 J=1,M
   27 	XIDENT(I,J)=0.D0
      DO 30 I=1,M
      DO 30 J=1,M
	      DO 28 L=1,M
   28		 XIDENT(I,J)=XIDENT(I,J) + B(I,L)*BSTORE(L,J)
C        XIDENT = B INVERSE  B     WHICH SHOULD BE THE IDENTITY MATRIX I
C        INVERSION ROUTINE WORKED.  ITEST2 CHECKS THAT XIDENT IS THE IDE
C        MATRIX
	      IF (I.EQ.J) GO TO 29
		      IF (DABS(XIDENT(I,J)).GT.1.D-06) ITEST2=1
		      GO TO 30
   29	      CONTINUE
	      IF (DABS(XIDENT(I,I)-1.D0).GT.1.D-06) ITEST2=-1
   30 CONTINUE
      IF (ITEST2.NE.0) GO TO 50
	      DO 33 I=1,M
		      XHAT(I)=0.D0
		      DO 32 J=1,N
		      PRODCT(J)=0.D0
			      DO 31 L=1,M
   31				 PRODCT(J)=PRODCT(J) + B(I,L)*C(J,L)
   32		      XHAT(I)=XHAT(I)+PRODCT(J)*SNEW(J)
   33	      CONTINUE
C        XHAT'S ARE (IF ADDED TO THE X'S) THE UNCONSTRAINED LEAST SQUARE
C        SOLUTION OF THE 'SOLVE FOR' PARAMETERS ONLY.
	      IN=0
C        XHAT2(J) = THE LEAST SQUARES SOLTIONS (INCLUDING NON SOLVE FOR
C        PARAMETERS
	      DO 43 J=1,NUMXES
		      IF (ID(J).ne.1) GO TO 42
		      IN=IN+1
		      XHAT2(J)=XHAT(IN) + X(J)
		      GO TO 43
   42 			XHAT2(J)=X(J)
   43	      CONTINUE
C
   50 CONTINUE

      DO 51 I=1,M
      DO 51 J=1,M
   51	 COR(I,J)=-BSTORE(I,J)/DSQRT(BSTORE(I,I)*BSTORE(J,J))
C
C
      IF (ITEST2.NE.0) GO TO 999
      INOTE=0
      DO 53 I=1,M
	      IF (B(I,I).GT.0.D0) GO TO 53
	      B(I,I)=-B(I,I)
	      INOTE = INOTE + 1
   53 CONTINUE
      IF (INOTE.EQ.0) GO TO 71
	WRITE (6,2071) INOTE
 2071	FORMAT (/////' ***** THERE ARE ',I2,' NEGATIVE VARIANCES'/'
     1             YOUR SOLUTION IS UNPHYSICAL')
   71 CONTINUE
      DO 54 J=1,M
   54	 STDEV(J)=DSQRT(B(J,J))
C        REARRANGE CALCULATED 1 SIGMA ERRORS TO APPLY TO THE SOLVED FOR
C        PARAMETERS
      IN=0
      DO 59 J=1,NUMXES
	      IF (ID(J).ne.1) GO TO 58
		      IN=IN+1
		      ERROR(J)=STDEV(IN)
		      GO TO 59
   58	      ERROR(J)=0.D0
   59 CONTINUE

      IF ( PRINT_COR )  THEN
C        OUTPUT FOR THE UNCONSTRAINED LEAST SQUARES SOLUTION
         WRITE (6,2040)
C
 2040    FORMAT (/1X,'    PARAMETER        ORIGINAL VALUE   ',
     1        'LEAST SQUARE VALUE  1 SIGMA ERRORS   SOLVED FOR?')
C
C
          DO 44 J=1,NUMXES
	      L=1
	      IF (ID(J).ne.1) L=2
  44	  WRITE (6,2041) VNAME(J), X(J), XHAT2(J), ERROR(J), QUEST(L)
 2041	  FORMAT (2X,A16,5X,3(F13.6,5X),4X,A8)
C
C
C	        CALCULATE THE CORRELATION COEFFICIENTS (COR)... 
	      DO 55 I=1,M
	      DO 55 J=1,M
   55		 COR(I,J)=B(I,J)/DSQRT(B(I,I)*B(J,J))
C
	      WRITE (6,2056)
 2056	      FORMAT(/10X,' THE CORRELATION COEFFICIENTS')
	      DO 57 I=1,M
   57		      WRITE (6,9998) (COR(I,J),J=1,M)

C  	      THE MULTIPLE CORRELATION COEFFICIENTS (PR) ARE CALCULATED
	      DO 60 I=1,M
   60		 PR(I)= 1.D0 - (1.d0/(BSTORE(I,I)*B(I,I)))
	      WRITE (6,2060)
 2060         FORMAT (//10X,'THE MULTIPLE CORRELATION COEFFICIENTS')

C			not sure if the true PR is the sqrt(PR)??
              I = 0
              DO 61 J=1,NUMXES
              	IF ( ID(J).ne.1 ) GO TO 61
              		I = I + 1
                        WRITE (6,2061) VNAME(J),PR(I)
 2061                   FORMAT (10X,A16,2X,F10.5)
   61         CONTINUE

  648 END IF
C
C
      GO TO 989
  998	WRITE (6,2094)  N,M
 2094	FORMAT (////' LEAST_SQUARES_FIT: # DATA POINTS (',I4,
     +              ') < # PARAMETERS (',I2,').  STOP' )
	STOP
C
  999 	WRITE (6,2095)  ITEST2
 2095 	FORMAT (////'  ***** ITEST2 =',I2,' ***** '
     +             /'       MATRIX INVERSION FAILED.')
	STOP
C
  989 RETURN
      END
c
      SUBROUTINE MINV(A,N,D,L,M,MSIZE)

      IMPLICIT REAL*8 (A-H,O-Z)
C     IBM SCIENTIFIC SUBROUTINE PACAKAGE PAGE 118
C     ******************************************************************
C
C     PURPOSE  INVERT A MATRIX
C
C     DESCRIPTION OF PARAMETERS
C        A  INPUT MATRIX, DESTROYED IN COMPUTATION AND REPLACED BY INVER
C        N  ORDER OF MATRIX A
C        D  RESULTANT DETERMINANT
C        L  WORK VECTOR OF LENGTH N
C        M  WORK VECTOR OF LENGTH N
C     MSIZE ORDER OF TWO DIMENSIONAL MATRIX A IN MAIN PROGRAM
C
C     METHOD
C        THE STANDARD GAUSS-JORDAN METHOD IS USED. THE DETERMINANT IS A
C        CALCULATED. A DETERMINANT OF ZERO INDICATES THAT THE MATRIX IS
C        SINGULAR.
C
C     ******************************************************************
C
      DIMENSION A(1),L(1),M(1)
C
C     STORE MATRIX A IN VECTOR FORM (COLUMN BY COLUMN)
C     SEE SSP PROGRAM ARRAY  P.98
C
      IF(MSIZE.EQ.N) GO TO 2
      NI=MSIZE-N
      IJ=0
      NM=0
      DO 230 K=1,N
      DO 225 LL=1,N
      IJ=IJ+1
      NM=NM+1
  225 A(IJ)=A(NM)
  230 NM=NM+NI
    2 CONTINUE
C
C     SEARCH FOR LARGEST ELEMENT
C
      D=1.0D0
      NK=-N
      DO 80 K=1,N
      NK=NK+N
      L(K)=K
      M(K)=K
      KK=NK+K
      BIGA=A(KK)
      DO 20 J=K,N
      IZ=N*(J-1)
      DO 20 I=K,N
      IJ=IZ+I
   10 IF(ABS(BIGA)-ABS(A(IJ))) 15,20,20
   15 BIGA=A(IJ)
      L(K)=I
      M(K)=J
   20 CONTINUE
C
C     INTERCHANGE ROWS
C
      J=L(K)
      IF(J-K) 35,35,25
   25 KI=K-N
      DO 30 I=1,N
      KI=KI+N
      HOLD=-A(KI)
      JI=KI-K+J
      A(KI)=A(JI)
   30 A(JI)=HOLD
C
C     INTERCHANGE COLUMNS
C
   35 I=M(K)
      IF(I-K) 45,45,38
   38 JP=N*(I-1)
      DO 40 J=1,N
      JK=NK+J
      JI=JP+J
      HOLD=-A(JK)
      A(JK)=A(JI)
   40 A(JI)=HOLD
C
C     DIVIDE COLUMN BY MINUS PIVOT (VALUE OF PIVOT ELEMENT IS CONTAINED
C     BIGA)
C
   45 IF(BIGA) 48,46,48
   46 D=0.0D0
      RETURN
   48 DO 55 I=1,N
      IF(I-K) 50,55,50
   50 IK=NK+I
      A(IK)=A(IK)/(-BIGA)
   55 CONTINUE
C
C     REDUCE MATRIX
C
      DO 65 I=1,N
      IK=NK+I
      HOLD=A(IK)
      IJ=I-N
      DO 65 J=1,N
      IJ=IJ+N
      IF(I-K) 60,65,60
   60 IF(J-K) 62,65,62
   62 KJ=IJ-I+K
      A(IJ)=HOLD*A(KJ)+A(IJ)
   65 CONTINUE
C
C     DIVIDE ROW BY PIVOT
C
      KJ=K-N
      DO 75 J=1,N
      KJ=KJ+N
      IF(J-K) 70,75,70
   70 A(KJ)=A(KJ)/BIGA
   75 CONTINUE
C
C     PRODUCT OF PIVOTS
C
      D=D*BIGA
C
C     REPLACE PIVOT BY RECIPROCAL
C
      A(KK)=1.0D0/BIGA
   80 CONTINUE
C
C     FINAL ROW AND COLUMN INTERCHANGE
C
      K=N
  100 K=K-1
      IF(K) 150,150,105
  105 I=L(K)
      IF(I-K) 120,120,108
  108 JQ=N*(K-1)
      JR=N*(I-1)
      DO 110 J=1,N
      JK=JQ+J
      HOLD=A(JK)
      JI=JR+J
      A(JK)=-A(JI)
  110 A(JI)=HOLD
  120 J=M(K)
      IF(J-K) 100,100,125
  125 KI=K-N
      DO 130 I=1,N
      KI=KI+N
      HOLD=A(KI)
      JI=KI-K+J
      A(KI)=-A(JI)
  130 A(JI)=HOLD
      GO TO 100
  150 CONTINUE
C
C     PUT MATRIX BACK INTO SQUARE FORM
C
      IF(MSIZE.EQ.N) GO TO 4
      IJ=N*N+1
      NM=N*MSIZE+1
      DO 210 K=1,N
      NM=NM-NI
      DO 210 LL=1,N
      IJ=IJ-1
      NM=NM-1
  210 A(NM)=A(IJ)
    4 CONTINUE
      RETURN
      END

c============================================================
c         Below is from NOVAS package
c============================================================
*** SOLSYS VERSION 3 PACKAGE: SOLSYS, SUN, IDSS ***

      SUBROUTINE SOLSYS (TJD,M,K,POS,VEL,IERR)
*
*     SUBROUTINE SOLSYS VERSION 3.
*     THIS SUBROUTINE PROVIDES THE POSITION AND VELOCITY OF THE
*     EARTH AT EPOCH TJD BY EVALUATING A CLOSED-FORM THEORY WITHOUT
*     REFERENCE TO AN EXTERNAL FILE.  THIS ROUTINE CAN ALSO PROVIDE
*     THE POSITION AND VELOCITY OF THE SUN.
*
*          TJD  = TDB JULIAN DATE OF DESIRED EPOCH (IN)
*          M    = BODY IDENTIFICATION NUMBER (IN)
*                 SET M=0 OR M=1 FOR THE SUN
*                 SET M=2 OR M=3 FOR THE EARTH
*          K    = ORIGIN SELECTION CODE (IN)
*                 SET K=0 FOR ORIGIN AT SOLAR SYSTEM BARYCENTER
*                 SET K=1 FOR ORIGIN AT CENTER OF SUN
*          POS  = POSITION VECTOR, EQUATORIAL RECTANGULAR
*                 COORDINATES, REFERRED TO MEAN EQUATOR AND EQUINOX
*                 OF J2000.0, COMPONENTS IN AU (OUT)
*          VEL  = VELOCITY VECTOR, EQUATORIAL RECTANGULAR
*                 COORDINATES, REFERRED TO MEAN EQUATOR AND EQUINOX
*                 OF J2000.0, COMPONENTS IN AU/DAY (OUT)
*          IERR = ERROR INDICATOR (OUT)
*                 IERR=0 MEANS EVERYTHING OK
*                 IERR=1 MEANS TJD BEFORE FIRST ALLOWED DATE
*                 IERR=2 MEANS TJD AFTER LAST ALLOWED DATE
*
*
       DOUBLE PRECISION TJD,POS,VEL,PI,TWOPI,T0,OBL,EL,C,P,TLAST,
     .     PM,PA,PE,PJ,PO,PW,PL,PN,       
     .     TMASS,SE,CE,SI,CI,SN,CN,SW,CW,P1,P2,P3,Q1,Q2,Q3,ROOTE,A,B,
     .     QJD,E,MLON,MA,U,SINU,COSU,ANR,PPLAN,VPLAN,F,PBARY,VBARY,
     .     DFLOAT,DABS,DMOD,DSIN,DCOS,DSQRT
      DIMENSION POS(3), VEL(3), EL(21), C(13), P(3,3),
     .     PM(4), PA(4), PE(4), PJ(4), PO(4), PW(4), PL(4), PN(4), 
     .     A(3,4), B(3,4), PPLAN(3), VPLAN(3), PBARY(3), VBARY(3)
      SAVE
      
      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( TWOPI  = 2.D0 * PI             )
      PARAMETER ( T0     = 2451545.0D0           )
      PARAMETER ( OBL    = 23.43927944D0         )
*     T0 = TDB JULIAN DATE OF EPOCH J2000.0
*     OBL = OBLIQUITY OF ECLIPTIC AT EPOCH J2000.0
      
      DATA EL, C, P / 43*0.D0 /,   TLAST / 0.D0 /

*     ARRAYS BELOW CONTAIN MASSES AND ORBITAL ELEMENTS OF THE FOUR
*     LARGEST PLANETS (SEE EXPLANATORY SUPPLEMENT (1992), P. 316)
*     WITH ANGLES IN RADIANS
*     THIS DATA USED FOR BARYCENTER COMPUTATIONS ONLY
*                 JUPITER        SATURN        URANUS       NEPTUNE
      DATA PM /  1047.349D 0,  3497.898D 0,   22903.0D 0,   19412.2D 0 /
      DATA PA /  5.203363D 0,  9.537070D 0, 19.191264D 0, 30.068963D 0 /
      DATA PE /  0.048393D 0,  0.054151D 0,  0.047168D 0,  0.008586D 0 /
      DATA PJ /  0.022782D 0,  0.043362D 0,  0.013437D 0,  0.030878D 0 /
      DATA PO /  1.755036D 0,  1.984702D 0,  1.295556D 0,  2.298977D 0 / 
      DATA PW /  0.257503D 0,  1.613242D 0,  2.983889D 0,  0.784898D 0 /
      DATA PL /  0.600470D 0,  0.871693D 0,  5.466933D 0,  5.321160D 0 /
      DATA PN /  1.450138D-3,  5.841727D-4,  2.047497D-4,  1.043891D-4 /
     
      IF ( TLAST .LT. 1.D0 ) THEN
*         FIRST TIME COMPUTATIONS
*         MASS OF SUN PLUS FOUR INNER PLANETS
          TMASS = 1.D0 + 5.977D-6
          SE = DSIN ( OBL * PI / 180.D0 )
          CE = DCOS ( OBL * PI / 180.D0 )
          DO 15 I = 1, 4
              TMASS = TMASS + 1.D0 / PM(I)
*             COMPUTE SINE AND COSINE OF ORBITAL ANGLES
			  SI = DSIN ( PJ(I) )
			  CI = DCOS ( PJ(I) )
			  SN = DSIN ( PO(I) )
			  CN = DCOS ( PO(I) )
			  SW = DSIN ( PW(I) - PO(I) )
			  CW = DCOS ( PW(I) - PO(I) )
*             COMPUTE P AND Q VECTORS (SEE BROUWER & CLEMENCE (1961), 
*             METHODS OF CELESTIAL MECHANICS, PP. 35-36.)
			  P1 =    CW * CN - SW * SN * CI
			  P2 = (  CW * SN + SW * CN * CI ) * CE - SW * SI * SE
			  P3 = (  CW * SN + SW * CN * CI ) * SE + SW * SI * CE
			  Q1 =   -SW * CN - CW * SN * CI
			  Q2 = ( -SW * SN + CW * CN * CI ) * CE - CW * SI * SE
			  Q3 = ( -SW * SN + CW * CN * CI ) * SE + CW * SI * CE
              ROOTE = DSQRT ( 1.D0 - PE(I)**2 )
              A(1,I) = PA(I) * P1
              A(2,I) = PA(I) * P2
              A(3,I) = PA(I) * P3
              B(1,I) = PA(I) * ROOTE * Q1
              B(2,I) = PA(I) * ROOTE * Q2
              B(3,I) = PA(I) * ROOTE * Q3
  15      CONTINUE
          TLAST = 1.D0
      END IF 
      
      IERR = 0
*     VALID DATES ARE WITHIN 3 CENTURIES OF J2000, ALTHOUGH RESULTS
*     DETERIORATE GRADUALLY      
      IF ( TJD .LT. 2340000.5D0 ) IERR = 1
      IF ( TJD .GT. 2560000.5D0 ) IERR = 2
      IF ( IERR .NE. 0 ) GO TO 110
      IF ( M .GE. 2 ) GO TO 30

*     FORM HELIOCENTRIC COORDINATES OF SUN
  20  DO 25 J=1,3
          POS(J) = 0.D0
          VEL(J) = 0.D0
  25  CONTINUE
      IF ( K .GE. 1 ) GO TO 110
      GO TO 90
      
*     FORM HELIOCENTRIC COORDINATES OF EARTH
*     VELOCITIES ARE OBTAINED FROM CRUDE NUMERICAL DIFFERENTIATION
  30  DO 35 I = 1, 3
          QJD = TJD + DFLOAT(I-2) * 0.1D0
C         SUBROUTINE SUN COMPUTES EARTH-SUN VECTOR
          CALL SUN ( QJD, EL, C )
          CALL PRECES ( QJD, C(11), T0, POS )
          P(I,1) = -POS(1)
          P(I,2) = -POS(2)
          P(I,3) = -POS(3)
  35  CONTINUE
      DO 40 J=1,3
          POS(J) =   P(2,J)
          VEL(J) = ( P(3,J) - P(1,J) ) / 0.2D0  
  40  CONTINUE
      IF ( K .GE. 1 ) GO TO 110

*     IF K=0, MOVE ORIGIN TO SOLAR SYSTEM BARYCENTER
*     SOLAR SYSTEM BARYCENTER COORDINATES COMPUTED FROM KEPLERIAN
*     APPROXIMATIONS OF THE COORDINATES OF THE FOUR LARGEST PLANETS
  90  IF ( DABS ( TJD - TLAST ) .LT. 1.D-6 ) GO TO 99
      DO 92 J = 1, 3
          PBARY(J) = 0.D0
          VBARY(J) = 0.D0
  92  CONTINUE 
*     THE FOLLOWING LOOP CYCLES ONCE FOR EACH OF THE FOUR LARGE PLANETS
      DO 98 I = 1, 4
*         COMPUTE MEAN LONGITUDE, MEAN ANOMALY, AND ECCENTRIC ANOMOLY
          E = PE(I)
          MLON = PL(I) + PN(I) * ( TJD - T0 )
          MA = DMOD ( MLON - PW(I), TWOPI )
          U = MA + E * DSIN ( MA ) + 0.5D0 * E * E * DSIN ( 2.D0 * MA )
          SINU = DSIN ( U )
          COSU = DCOS ( U )
*         COMPUTE VELOCITY FACTOR     
          ANR = PN(I) / ( 1.D0 - E * COSU )
*         COMPUTE PLANET'S POSITION AND VELOCITY WRT EQ & EQ J2000
          PPLAN(1) = A(1,I) * ( COSU - E ) + B(1,I) * SINU 
          PPLAN(2) = A(2,I) * ( COSU - E ) + B(2,I) * SINU
          PPLAN(3) = A(3,I) * ( COSU - E ) + B(3,I) * SINU
          VPLAN(1) = ANR * ( -A(1,I) * SINU + B(1,I) * COSU )
          VPLAN(2) = ANR * ( -A(2,I) * SINU + B(2,I) * COSU )
          VPLAN(3) = ANR * ( -A(3,I) * SINU + B(3,I) * COSU )
*         COMPUTE MASS FACTOR AND ADD IN TO TOTAL DISPLACEMENT      
          F = 1.D0 / ( PM(I) * TMASS )
          PBARY(1) = PBARY(1) + PPLAN(1) * F
          PBARY(2) = PBARY(2) + PPLAN(2) * F
          PBARY(3) = PBARY(3) + PPLAN(3) * F
          VBARY(1) = VBARY(1) + VPLAN(1) * F
          VBARY(2) = VBARY(2) + VPLAN(2) * F
          VBARY(3) = VBARY(3) + VPLAN(3) * F
  98  CONTINUE
      TLAST = TJD
  99  DO 100 J=1,3
          POS(J) = POS(J) - PBARY(J)
          VEL(J) = VEL(J) - VBARY(J)
 100  CONTINUE 

 110  RETURN

      END



      SUBROUTINE SUN (DJ,EL,C)
C
C     FOR USE WITH SUBROUTINE SOLSYS VERSION 3.
C     THIS SUBROUTINE COMPUTES THE COORDINATES OF THE EARTH-SUN
C     POSITION VECTOR WITH RESPECT TO THE ECLIPTIC AND EQUATOR
C     OF DATE.  A MODIFIED FORM OF NEWCOMB'S THEORY ('TABLES OF THE
C     SUN', 1898) IS USED.  ONLY THE LARGEST PERIODIC PERTURBATIONS
C     ARE EVALUATED, AND VAN FLANDERN'S EXPRESSIONS FOR THE FUNDAMENTAL
C     ARGUMENTS ('IMPROVED MEAN ELEMENTS FOR THE EARTH AND MOON', 1981)
C     ARE USED.  THE ABSOLUTE ACCURACY IS NO WORSE THAN 1 ARCSECOND
C     (AVERAGE ERROR ABOUT 0.2 ARCSECOND) OVER 1800-2200.
C     (ADAPTED FROM SUBROUTINE IAUSUN BY P. M. JANICZEK, USNO.)
C
C          DJ   = TDB JULIAN DATE OF DESIRED EPOCH (IN)
C          EL   = ARRAY OF ORBITAL ELEMENTS (SEE BELOW) FOR
C                 EPOCH DJ (OUT)
C          C    = ARRAY OF COORDINATES (SEE BELOW) FOR
C                 EPOCH DJ (OUT)
C
C
      DOUBLE PRECISION DJ,EL,C,T,TP,T20,RO,GV,GM,GJ,GS,DL,DR,DB,DG,
     1 DBLARG,D,TWOPI,STR,RTD,R,TR,
     2 SINO,COSO,SINL,COSL,SINB,COSB,
     3 DSIN,DCOS,DMOD
C
      DIMENSION EL(21)
C
C     EL( 1)= SEMI-MAJOR AXIS, AU
C     EL( 2)= ORBITAL ECCENTRICITY
C     EL( 5)= LONGITUDE OF PERIGEE, RADIANS
C     EL( 9)= UNPERTURBED MEAN LONGITUDE, RADIANS
C     EL(10)= MEAN ANOMALY, AFFECTED BY LONG-PD PERTURBATIONS, RADIANS
C     EL(11)= UNPERTURBED RADIUS, AU
C     EL(12)= EQUATION OF THE CENTER, RADIANS
C     EL(13)= MEAN OBLIQUITY OF ECLIPTIC, RADIANS
C     EL(14)= MEAN LONGITUDE OF MOON, RADIANS
C     EL(15)= MEAN ANOMALY OF MOON, RADIANS
C     EL(16)= LUNAR MEAN ARGUMENT OF LATITUDE, RADIANS
C     EL(17)= MEAN LONGITUDE OF LUNAR ASCENDING NODE, RADIANS
C     EL(21)= MEAN LONGITUDE OF MOON'S PERIGEE, RADIANS
C             (REMAINING ELEMENTS OF ARRAY EL NOT USED)
C
      DIMENSION C(13)
C
C     C( 1) = PERTURBED RADIUS VECTOR, AU
C     C( 2) = SAME AS C(4), DEGREES
C     C( 3) = SAME AS C(5), DEGREES
C     C( 4) = ECLIPTIC LONGITUDE WRT MEAN ECL & EQUX OF DATE, RADIANS
C     C( 5) = ECLIPTIC LATITUDE  WRT MEAN ECL        OF DATE, RADIANS
C     C(11) = EQUATORIAL X WRT MEAN EQU & EQUX OF DATE, AU
C     C(12) = EQUATORIAL Y WRT MEAN EQU & EQUX OF DATE, AU
C     C(13) = EQUATORIAL Z WRT MEAN EQU & EQUX OF DATE, AU
C             (REMAINING ELEMENTS OF ARRAY C NOT USED)
C
C
C***********************************************************************
C
C     PART I    TABLES OF THE PERTURBATIONS
C
      DIMENSION X(8,46), X1(80), X2(80), X3(80), X4(80), X5(48)
      EQUIVALENCE (X(1, 1),X1(1))
      EQUIVALENCE (X(1,11),X2(1))
      EQUIVALENCE (X(1,21),X3(1))
      EQUIVALENCE (X(1,31),X4(1))
      EQUIVALENCE (X(1,41),X5(1))
C
C     PERTURBATIONS BY VENUS
C                  J    I     VC      VS    RHOC    RHOS      BC     BS
      DATA X1 /  - 1.,  0., +  33.,-  67., -  85.,-  39., +  24.,-  17.,
     2           - 1.,+ 1., +2353.,-4228., -2062.,-1146., -   4.,+   3.,
     3           - 1.,+ 2., -  65.,-  34., +  68.,-  14., +   6.,-  92.,
     4           - 2.,+ 1., -  99.,+  60., +  84.,+ 136., +  23.,-   3.,
     5           - 2.,+ 2., -4702.,+2903., +3593.,+5822., +  10.,-   6.,
     6           - 2.,+ 3., +1795.,-1737., - 596.,- 632., +  37.,-  56.,
     7           - 3.,+ 3., - 666.,+  27., +  44.,+1044., +   8.,+   1.,
     8           - 3.,+ 4., +1508.,- 397., - 381.,-1448., + 185.,- 100.,
     9           - 3.,+ 5., + 763.,- 684., + 126.,+ 148., +   6.,-   3.,
     *           - 4.,+ 4., - 188.,-  93., - 166.,+ 337.,     0.,    0./
      DATA X2 /  - 4.,+ 5., - 139.,-  38., -  51.,+ 189., -  31.,-   1.,
     2           - 4.,+ 6., + 146.,-  42., -  25.,-  91., +  12.,    0.,
     3           - 5.,+ 5., -  47.,-  69., - 134.,+  93.,     0.,    0.,
     4           - 5.,+ 7., - 119.,-  33., -  37.,+ 136., -  18.,-   6.,
     5           - 5.,+ 8., + 154.,    0.,     0.,-  26.,     0.,    0.,
     6           - 6.,+ 6., -   4.,-  38., -  80.,+   8.,     0.,    0.,
C
C     PERTURBATIONS BY MARS
C                  J    I     VC      VS    RHOC    RHOS      BC     BS
     7           + 1.,- 1., - 216.,- 167., -  92.,+ 119.,     0.,    0.,
     8           + 2.,- 2., +1963.,- 567., - 573.,-1976.,     0.,-   8.,
     9           + 2.,- 1., -1659.,- 617., +  64.,- 137.,     0.,    0.,
     *           + 3.,- 3., +  53.,- 118., - 154.,-  67.,     0.,    0./
      DATA X3 /  + 3.,- 2., + 396.,- 153., -  77.,- 201.,     0.,    0.,
     2           + 4.,- 3., - 131.,+ 483., + 461.,+ 125., +   7.,+   1.,
     3           + 4.,- 2., + 526.,- 256., +  43.,+  96.,     0.,    0.,
     4           + 5.,- 4., +  49.,+  69., +  87.,-  62.,     0.,    0.,
     5           + 5.,- 3., -  38.,+ 200., +  87.,+  17.,     0.,    0.,
     6           + 6.,- 4., - 104.,- 113., - 102.,+  94.,     0.,    0.,
     7           + 6.,- 3., -  11.,+ 100., -  27.,-   4.,     0.,    0.,
     8           + 7.,- 4., -  78.,-  72., -  26.,+  28.,     0.,    0.,
     9           + 9.,- 5., +  60.,-  15., -   4.,-  17.,     0.,    0.,
     *           +15.,- 8., + 200.,-  30., -   1.,-   6.,     0.,    0./
C
C     PERTURBATIONS BY JUPITER
C                  J    I     VC      VS    RHOC    RHOS      BC     BS
      DATA X4 /  + 1.,- 2., - 155.,-  52., -  78.,+ 193., +   7.,    0.,
     2           + 1.,- 1., -7208.,+  59., +  56.,+7067., -   1.,+  17.,
     3           + 1.,  0., - 307.,-2582., + 227.,-  89., +  16.,    0.,
     4           + 1.,+ 1., +   8.,-  73., +  79.,+   9., +   1.,+  23.,
     5           + 2.,- 3., +  11.,+  68., + 102.,-  17.,     0.,    0.,
     6           + 2.,- 2., + 136.,+2728., +4021.,- 203.,     0.,    0.,
     7           + 2.,- 1., - 537.,+1518., +1376.,+ 486., +  13.,+ 166.,
     8           + 3.,- 3., - 162.,+  27., +  43.,+ 278.,     0.,    0.,
     9           + 3.,- 2., +  71.,+ 551., + 796.,- 104., +   6.,-   1.,
     *           + 3.,- 1., -  31.,+ 208., + 172.,+  26., +   1.,+  18./
      DATA X5 /  + 4.,- 3., -  43.,+   9., +  13.,+  73.,     0.,    0.,
     2           + 4.,- 2., +  17.,+  78., + 110.,-  24.,     0.,    0.,
C
C     PERTURBATIONS BY SATURN
C                  J    I     VC      VS    RHOC    RHOS      BC     BS
     3           + 1.,- 1., -  77.,+ 412., + 422.,+  79., +   1.,+   6.,
     4           + 1.,  0., -   3.,- 320., +   8.,-   1.,     0.,    0.,
     5           + 2.,- 2., +  38.,- 101., - 152.,-  57.,     0.,    0.,
     6           + 2.,- 1., +  45.,- 103., - 103.,-  44.,     0.,    0./
C
C
C***********************************************************************
C
C     PART II   NECESSARY PRELIMINARIES
C
      DATA TWOPI /6.283185307179586D0/
      DATA STR   /206264806.2470964D0/
      DATA RTD   /57.295779513082321D0/
      DATA R     /1296000.0D0/
      TR = 1000.0D0 / STR
C
C     T  = TIME IN JULIAN CENTURIES FROM 1900 JANUARY 0
      T  = (DJ - 2415020.D0)/36525.D0
C
C     TP = TIME IN JULIAN YEARS     FROM 1850 JANUARY 0
      TP = (DJ - 2396758.D0)/365.25D0
C
C     T20= TIME IN JULIAN CENTURIES FROM J2000.0
      T20= (DJ - 2451545.D0)/36525.D0
C
C
C***********************************************************************
C
C     PART III  COMPUTATION OF ELLIPTIC ELEMENTS AND SECULAR TERMS
C
C     VAN FLANDERN'S EXPRESSIONS FOR MEAN ELEMENTS
      EL( 1) = 1.00000030007166D0
      EL( 2) = 0.016708320D0 + (-0.42229D-04 - 0.126D-06 * T20) * T20
      EL( 5) = 1018578.046D0 + (6190.046D0 +
     1                (1.666D0 + 0.012D0 * T20) * T20) * T20
      EL( 5) = EL( 5) * TR
      EL( 9) = 1009677.850D0 + (100.0D0 * R + 2771.27D0 +
     1                1.089D0 * T20) * T20
      EL( 9) = DMOD (EL( 9) * TR, TWOPI)
      EL(10) = 1287099.804D0 + (99.0D0 * R + 1292581.224D0 +
     1                (-0.577D0 - 0.012D0 * T20) * T20) * T20
      EL(10) = DMOD (EL(10) * TR, TWOPI)
      
C     EXPRESSION FOR OBLIQUITY FROM P03 (IAU 2006) PRECESSION       
      EL(13) = 84381.406D0 + (-46.836769D0 +
     1               (-0.0001831D0 + 0.00200340D0 * T20) * T20) * T20
      EL(13) = EL(13) * TR

C     KAPLAN CORRECTION TO SUN'S MEAN LONGITUDE TO FIT DE405 OVER
C     INTERVAL 1800-2200, USING P03 (IAU 2006) PRECESSION
      EL(9) = EL(9)
     1      + ( 0.1320D0 - 0.1355D0 * T20 ) * TR

C
C***********************************************************************
C
C     PART IV   LUNAR TERMS
C
C     VAN FLANDERN'S EXPRESSIONS FOR MEAN ELEMENTS
      EL(14) = 785939.157D0 + (1336.0D0 * R + 1108372.598D0
     1                + (-5.802D0 + 0.019D0 * T20) * T20) * T20
      EL(14) = DMOD (EL(14) * TR, TWOPI)
      EL(17) = 450160.280D0 + (-5.0D0 * R - 482890.539D0 +
     1                (7.455D0 + 0.008D0 * T20) * T20) * T20
      EL(17) = DMOD (EL(17) * TR, TWOPI)
      EL(21) = 300072.424D0 + (11.0D0 * R + 392449.965D0 +
     1                (-37.112D0 - 0.045D0 * T20) * T20) * T20
      EL(21) = DMOD (EL(21) * TR, TWOPI)
C
C     DERIVED ARGUMENTS
      EL(15) = EL(14) - EL(21)
      EL(16) = EL(14) - EL(17)
      EL(15) = DMOD (EL(15),TWOPI)
      EL(16) = DMOD (EL(16),TWOPI)
C     MEAN ELONGATION
      D      = EL(14) - EL(9)
C
C     COMBINATIONS OF ARGUMENTS AND THE PERTURBATIONS
      D = DMOD (D,TWOPI)
      ARG = D
      DL =    +  6469.*SIN(ARG) +  13.*SIN(3.*ARG)
      DR =    + 13390.*COS(ARG) +  30.*COS(3.*ARG)
C
      DBLARG = D + EL(15)
      DBLARG = DMOD (DBLARG,TWOPI)
      ARG = DBLARG
      DL = DL +  177.*SIN(ARG)
      DR = DR +  370.*COS(ARG)
C
      DBLARG = D - EL(15)
      DBLARG = DMOD (DBLARG,TWOPI)
      ARG = DBLARG
      DL = DL -  424.*SIN(ARG)
      DR = DR - 1330.*COS(ARG)
C
      DBLARG = 3.D0*D - EL(15)
      DBLARG = DMOD (DBLARG,TWOPI)
      ARG = DBLARG
      DL = DL +   39.*SIN(ARG)
      DR = DR +   80.*COS(ARG)
C
      DBLARG = D + EL(10)
      DBLARG = DMOD (DBLARG,TWOPI)
      ARG = DBLARG
      DL = DL -   64.*SIN(ARG)
      DR = DR -  140.*COS(ARG)
C
      DBLARG = D - EL(10)
      DBLARG = DMOD (DBLARG,TWOPI)
      ARG = DBLARG
      DL = DL +  172.*SIN(ARG)
      DR = DR +  360.*COS(ARG)
C
      EL(16) = DMOD (EL(16),TWOPI)
      ARG = EL(16)
      DB =    + 576.*SIN(ARG)
C
C
C***********************************************************************
C
C     PART V    COMPUTATION OF PERIODIC PERTURBATIONS
C
C     THE PERTURBING MEAN ANOMALIES
C
      GV  = 0.19984020D+01 + .1021322923D+02*TP
      GM  = 0.19173489D+01 + .3340556174D+01*TP
      GJ  = 0.25836283D+01 + .5296346478D+00*TP
      GS  = 0.49692316D+01 + .2132432808D+00*TP
      GV  = DMOD (GV,TWOPI)
      GM  = DMOD (GM,TWOPI)
      GJ  = DMOD (GJ,TWOPI)
      GS  = DMOD (GS,TWOPI)
C
C
C     MODIFICATION OF FUNDAMENTAL ARGUMENTS
C
C     APPLICATION OF THE JUPITER-SATURN GREAT INEQUALITY
C     TO JUPITER'S MEAN ANOMALY
C
      GJ = GJ + 0.579904067D-02 * DSIN (5.D0*GS - 2.D0*GJ
     1                 + 1.1719644977D0 - 0.397401726D-03*TP)
      GJ = DMOD (GJ,TWOPI)
C
C     LONG PERIOD PERTURBATIONS OF MEAN ANOMALY
C
      ST = T
C                ARGUMENT IS ( 4 MARS - 7 EARTH + 3 VENUS )
      DG = 266.* SIN (0.555015 + 2.076942*ST)
C                ARGUMENT IS ( 3 JUPITER - 8 MARS + 4 EARTH )
     1    + 6400.* SIN (4.035027 + 0.3525565*ST)
C                ARGUMENT IS ( 13 EARTH - 8 VENUS )
     2    + (1882.-16.*ST) * SIN (0.9990265 + 2.622706*ST)
C
C
C     COMPUTATION OF THE EQUATION OF THE CENTER
C
C     FORM PERTURBED MEAN ANOMALY
      EL(10) = DG/STR + EL(10)
      EL(10) = DMOD (EL(10),TWOPI)
      EL(12) =   DSIN(     EL(10)) * (6910057.D0 -(17240.D0+52.D0*T)*T)
     1         + DSIN(2.D0*EL(10)) * (  72338.D0 -    361.D0*T)
     2         + DSIN(3.D0*EL(10)) * (   1054.D0 -      1.D0*T)
C
C     THE UNPERTURBED RADIUS VECTOR
      RO     =                          30570.D0 -    150.D0*T
     1         - DCOS(     EL(10)) * (7274120.D0 - (18140.D0+50.D0*T)*T)
     2         - DCOS(2.D0*EL(10)) * (  91380.D0 -    460.D0*T)
     3         - DCOS(3.D0*EL(10)) * (   1450.D0 -     10.D0*T)
      EL(11) = 10.D0**(RO*1.D-09)
C
C
C     SELECTED PLANETARY PERTURBATIONS FROM NEWCOMB'S THEORY FOLLOW
C
C     PERTURBATIONS BY VENUS
      DO 20 K=1,16
C     ARGUMENT J * VENUS +   I * EARTH
      DBLARG = X(1,K)*GV + X(2,K)*EL(10)
      DBLARG = DMOD (DBLARG,TWOPI)
      ARG = DBLARG
      CS  = COS(ARG)
      SS  = SIN(ARG)
      DL  =(X(3,K)*CS  + X(4,K)*SS )+ DL
      DR  =(X(5,K)*CS  + X(6,K)*SS )+ DR
      DB  =(X(7,K)*CS  + X(8,K)*SS )+ DB
   20 CONTINUE
C
C     PERTURBATIONS BY MARS
      DO 30 K=17,30
C     ARGUMENT  J * MARS +   I * EARTH
      DBLARG = X(1,K)*GM + X(2,K)*EL(10)
      DBLARG = DMOD (DBLARG,TWOPI)
      ARG = DBLARG
      CS  = COS(ARG)
      SS  = SIN(ARG)
      DL  =(X(3,K)*CS  + X(4,K)*SS )+ DL
      DR  =(X(5,K)*CS  + X(6,K)*SS )+ DR
      DB  =(X(7,K)*CS  + X(8,K)*SS )+ DB
   30 CONTINUE
C
C     PERTURBATIONS BY JUPITER
      DO 40 K=31,42
C     ARGUMENT J*JUPITER +   I * EARTH
      DBLARG = X(1,K)*GJ + X(2,K)*EL(10)
      DBLARG = DMOD (DBLARG,TWOPI)
      ARG = DBLARG
      CS  = COS(ARG)
      SS  = SIN(ARG)
      DL  =(X(3,K)*CS  + X(4,K)*SS )+ DL
      DR  =(X(5,K)*CS  + X(6,K)*SS )+ DR
      DB  =(X(7,K)*CS  + X(8,K)*SS )+ DB
   40 CONTINUE
C
C     PERTURBATIONS BY SATURN
      DO 50 K=43,46
C     ARGUMENT J*SATURN  +   I * EARTH
      DBLARG = X(1,K)*GS + X(2,K)*EL(10)
      DBLARG = DMOD (DBLARG,TWOPI)
      ARG = DBLARG
      CS  = COS(ARG)
      SS  = SIN(ARG)
      DL  =(X(3,K)*CS  + X(4,K)*SS )+ DL
      DR  =(X(5,K)*CS  + X(6,K)*SS )+ DR
      DB  =(X(7,K)*CS  + X(8,K)*SS )+ DB
   50 CONTINUE
C
C
C***********************************************************************
C
C     PART VI   COMPUTATION OF ECLIPTIC AND EQUATORIAL COORDINATES
C
      C(1) = EL(11)*10.D0**(DR*1.D-09)
      C(4) = (DL + DG + EL(12))/STR + EL(9)
      C(4) = DMOD (C(4),TWOPI)
      C(5) = DB/STR
      C(2) = C(4)*RTD
      C(3) = C(5)*RTD
      SINO = DSIN (EL(13))
      COSO = DCOS (EL(13))
      SINL = DSIN (C(4))
      COSL = DCOS (C(4))
      SINB = DSIN (C(5))
      COSB = DCOS (C(5))
      C(11) = C(1) * (COSB * COSL)
      C(12) = C(1) * (COSB * SINL * COSO - SINB * SINO)
      C(13) = C(1) * (COSB * SINL * SINO + SINB * COSO)
C
C
C***********************************************************************
C
C
      RETURN
C
      END



      INTEGER FUNCTION IDSS ( NAME )
*
*     THIS FUNCTION RETURNS THE ID NUMBER OF A SOLAR SYSTEM BODY
*     FOR THE VERSION OF SOLSYS (OR SOLSYS-AUXPOS COMBINATION) IN USE.
*
*         NAME   = NAME OF BODY WHOSE ID NUMBER IS DESIRED, E.G.,
*                  'SUN', 'MOON, 'MERCURY', ETC., EXPRESSED AS ALL
*                  UPPER-CASE LETTERS (IN)
*         IDSS   = ID NUMBER OF BODY, FOR USE IN CALLS TO SOLSYS
*                  (FUNCTION VALUE RETURNED)
*
*     NOTE 1: IN THIS VERSION, ONLY THE FIRST THREE LETTERS OF THE
*     BODY'S NAME ARE USED FOR IDENTIFICATION.  ALTERNATIVE VERSIONS
*     MIGHT USE MORE LETTERS.
*
*     NOTE 2: IF NAME IS 'JD', IDSS RETURNS IDSS=1, SINCE SOLSYS 
*     VERSION 3 DOES NOT PROCESS SPLIT JULIAN DATES.    
*
*     NOTE 3: ALL VERSIONS OF IDSS MUST RETURN IDSS=-9999 FOR OBJECTS
*     THAT IT CANNOT IDENTIFY OR ARE UNSUPPORTED BY SOLSYS.
*
*
      CHARACTER NAME*(*), NAMEIN*3, NAMES*3
      DIMENSION NAMES(35), IDS(35)

      DATA NAMES / 'SUN', 'EAR', '---', '---', '---', '---', '---',
     .             '---', '---', '---', '---', '---', '---', '---',
     .             '---', '---', '---', '---', '---', '---', '---',
     .             '---', '---', '---', '---', '---', '---', '---',
     .             '---', '---', '---', '---', '---', '---', '---'  /
      DATA IDS   /     0,     3,     0,     0,     0,     0,     0,
     .                 0,     0,     0,     0,     0,     0,     0,
     .                 0,     0,     0,     0,     0,     0,     0,
     .                 0,     0,     0,     0,     0,     0,     0,
     .                 0,     0,     0,     0,     0,     0,     0  /
      DATA NUM   / 2 /

   3  FORMAT ( ' IDSS ERROR: NO BODY ID NUMBER FOUND FOR ', A )

      IDSS = -9999
      NAMEIN = NAME

*     LOOK THROUGH LIST OF BODY NAMES TO FIND MATCH
      DO 20 I = 1, NUM
          IF ( NAMEIN .EQ. NAMES(I) ) THEN
              IDSS = IDS(I)
              GO TO 30
          END IF
  20  CONTINUE
  
*     IF NO MATCH, CHECK FOR INQUIRY ABOUT SPLIT JULIAN DATES   
      IF ( NAMEIN .EQ. 'JD ' ) THEN
*         IN THIS CASE, SET IDSS=2 IF SOLSYS PROCESSES SPLIT
*         JULIAN DATES (IN SUCCESSIVE CALLS), IDSS=1 OTHERWISE 
          IDSS = 1
          GO TO 30
      END IF    

      WRITE ( *, 3 ) NAME

  30  RETURN

      END

c=================================================================

*  NOVAS FORTRAN VERS F3.1 of 2011 MARCH 21                       
                                   
************************************************************************
*                                                                      *
*                              N O V A S                               *
*           NAVAL OBSERVATORY VECTOR ASTROMETRY SOFTWARE               *
*                                                                      *
*                            G. H. KAPLAN                              *
*                        U.S. NAVAL OBSERVATORY                        *
*                                                                      *
************************************************************************



      SUBROUTINE PLACE ( TJD, OBJECT, LOCATN, ICOORD, STAR, OBSERV,
     .                   SKYPOS )
*
*     THIS SUBROUTINE COMPUTES THE APPARENT DIRECTION OF A STAR OR SOLAR
*     SYSTEM BODY AT A SPECIFIED TIME AND IN A SPECIFIED COORDINATE
*     SYSTEM.  BASED ON KAPLAN, ET AL. (1989), ASTRONOMICAL JOURNAL 97,
*     1197-1210, WITH SOME ENHANCEMENTS FROM KLIONER (2003),
*     ASTRONOMICAL JOURNAL 125, 1580-1597.
*
*          TJD    = TT JULIAN DATE FOR PLACE (IN)
*          OBJECT = CHARACTER STRING IDENTIFYING OBJECT OF INTEREST (IN)
*                   FOR SOLAR SYSTEM
*                   BODY,             SPECIFY THE NAME USING ALL UPPER-
*                                     CASE LETTERS ('SUN', 'MOON',
*                                     'JUPITER', ETC.),
*                                     - OR -
*                                     SPECIFY THE BODY ID NUMBER
*                                     IN A 4-CHARACTER STRING OF THE
*                                     FORM '=NNN', WHERE NNN IS THE
*                                     BODY ID NUMBER
*                   FOR STAR,         PROVIDE A BLANK STRING, THE WORD
*                                     'STAR', OR ANY STRING BEGINNING
*                                     WITH '*'
*          LOCATN = INTEGER CODE SPECIFYING LOCATION OF OBSERVER (IN)
*                   SET LOCATN=0 FOR OBSERVER AT GEOCENTER
*                   SET LOCATN=1 FOR OBSERVER ON SURFACE OF EARTH
*                   SET LOCATN=2 FOR OBSERVER ON NEAR-EARTH SPACECRAFT
*          ICOORD = INTEGER CODE SPECIFYING COORDINATE SYSTEM OF OUTPUT
*                   POSITION (IN)
*                   SET ICOORD=0 FOR GCRS (OR 'LOCAL GCRS')
*                   SET ICOORD=1 FOR TRUE EQUATOR AND EQUINOX OF DATE
*                   SET ICOORD=2 FOR TRUE EQUATOR AND CIO OF DATE
*                   SET ICOORD=3 FOR ASTROMETRIC COORDINATES, I.E.,
*                                WITHOUT LIGHT DEFLECTION OR ABERRATION
*          STAR   = ARRAY OF CATALOG DATA FOR STAR (IN)
*                   (NOT USED IF SOLAR SYSTEM BODY REQUESTED)
*                   STAR(1) = ICRS RIGHT ASCENSION IN HOURS
*                   STAR(2) = ICRS DECLINATION IN DEGREES
*                   STAR(3) = ICRS PROPER MOTION IN RA IN
*                             MILLIARCSECONDS/YEAR
*                   STAR(4) = ICRS PROPER MOTION IN DEC IN
*                             MILLIARCSECONDS/YEAR
*                   STAR(5) = PARALLAX IN MILLIARCSECONDS
*                   STAR(6) = RADIAL VELOCITY IN KILOMETERS/SECOND
*                   FURTHER STAR ARRAY ELEMENTS ARE NOT USED HERE
*                   BUT ARE RESERVED FOR FUTURE USE
*          OBSERV = ARRAY OF DATA SPECIFYING LOCATION OF OBSERVER (IN)
*                   (NOT USED IF LOCATN=0)
*                   FOR LOCATN=1,
*                   OBSERV(1) = GEODETIC LONGITUDE (WGS-84) OF OBSERVER
*                               (EAST +) IN DEGREES 
*                   OBSERV(2) = GEODETIC LATITUDE (WGS-84) OF OBSERVER
*                               (NORTH +) IN DEGREES
*                   OBSERV(3) = HEIGHT OF OBSERVER ABOVE ELLIPSOID
*                               IN METERS
*                   OBSERV(4) = VALUE OF DELTA-T IN SECONDS
*                               (DELTA-T=TT-UT1)
*                   OBSERV(5) = (NOT USED, RESERVED FOR FUTURE USE)
*                   OBSERV(6) = (NOT USED, RESERVED FOR FUTURE USE)
*                   FOR LOCATN=2,
*                   OBSERV(1) = GEOCENTRIC X IN KILOMETERS
*                   OBSERV(2) = GEOCENTRIC Y IN KILOMETERS
*                   OBSERV(3) = GEOCENTRIC Z IN KILOMETERS
*                   OBSERV(4) = GEOCENTRIC X-DOT IN KILOMETERS/SECOND
*                   OBSERV(5) = GEOCENTRIC Y-DOT IN KILOMETERS/SECOND
*                   OBSERV(6) = GEOCENTRIC Z-DOT IN KILOMETERS/SECOND
*                   WITH RESPECT TO TRUE EQUATOR AND EQUINOX OF DATE
*          SKYPOS = ARRAY OF OUTPUT DATA SPECIFYING OBJECT'S PLACE
*                   ON THE SKY AT TIME TJD, WITH RESPECT TO THE 
*                   SPECIFIED OUTPUT COORDINATE SYSTEM (OUT)
*                   SKYPOS(1) = X, DIMENSIONLESS      UNIT VECTOR
*                   SKYPOS(2) = Y, DIMENSIONLESS      TOWARD OBJECT
*                   SKYPOS(3) = Z, DIMENSIONLESS
*                   SKYPOS(4) = APPARENT, TOPOCENTRIC, OR ASTROMETRIC
*                               RIGHT ASCENSION IN HOURS
*                   SKYPOS(5) = APPARENT, TOPOCENTRIC, OR ASTROMETRIC
*                               DECLINATION IN DEGREES
*                   SKYPOS(6) = TRUE (GEOMETRIC, EUCLIDIAN) DISTANCE
*                               TO SOLAR SYSTEM BODY IN AU AT TIME TJD,
*                               OR 0.D0 FOR STAR
*                   SKYPOS(7) = RADIAL VELOCITY IN KILOMETERS/SECOND 
*                   FURTHER SKYPOS ARRAY ELEMENTS ARE NOT USED HERE
*                   BUT ARE RESERVED FOR FUTURE USE
* 
*     NOTE 1: VALUES OF LOCATN AND ICOORD FOR VARIOUS STANDARD KINDS
*     OF PLACE:
*     LOCATN=0 AND ICOORD=1 APPARENT PLACE
*     LOCATN=1 AND ICOORD=1 TOPOCENTRIC PLACE
*     LOCATN=0 AND ICOORD=0 VIRTUAL PLACE
*     LOCATN=1 AND ICOORD=0 LOCAL PLACE
*     LOCATN=0 AND ICOORD=3 ASTROMETRIC PLACE
*     LOCATN=1 AND ICOORD=3 TOPOCENTRIC ASTROMETRIC PLACE
*
*     NOTE 2: ARRAYS STAR AND SKYPOS MAY BE EXPANDED IN THE FUTURE, AND
*     THIS CAN BE ALLOWED FOR IN THE CALLING CODE BY DIMENSIONING
*     THESE ARRAYS WITH 20 AND 10 ELEMENTS, RESPECTIVELY, EVEN THOUGH
*     ELEMENTS BEYOND STAR(6) AND SKYPOS(7) ARE NOT NOW REFERRED TO IN
*     THIS SUBROUTINE.
*
*     NOTE 3: IF LOCATN=1 AND OBSERV(4)=0.D0, THE VALUE OF DELTA-T WILL
*     BE OBTAINED FROM GETDT, WHICH PROVIDES THE LAST VALUE OF DELTA-T
*     DEFINED BY THE USER VIA CALL TO SETDT.
*
*     NOTE 4: SKYPOS(7), THE RADIAL VELOCITY, IS THE PREDICTED 
*     RADIAL VELOCITY MEASURE (Z) TIMES THE SPEED OF LIGHT, AN
*     INHERENTLY SPECTROSCOPIC MEASURE.  FOR A STAR, IT
*     INCLUDES ALL EFFECTS, SUCH AS GRAVITATIONAL RED SHIFT,
*     CONTAINED IN THE CATALOG BARYCENTRIC RADIAL VELOCITY MEASURE,
*     WHICH IS ASSUMED GIVEN IN STAR(6).  FOR A SOLAR SYSTEM
*     BODY, IT APPLIES TO A FICTITIOUS EMITTER AT THE CENTER OF THE
*     OBSERVED OBJECT, ASSUMED MASSLESS (NO GRAVITATIONAL RED SHIFT),
*     AND DOES NOT IN GENERAL APPLY TO REFLECTED LIGHT.
*
*

* --- INITIAL DECLARATIONS---------------------------------------------

      IMPLICIT NONE
      INTEGER LOCATN,ICOORD,NTIMES,IEARTH,ISUN,IDBODY,IERR,LOC,J,KCIO,
     .     IDSS 
      DOUBLE PRECISION TJD,STAR,OBSERV,SKYPOS,
     .     T0,TLAST1,TLAST2,TTJD,TDBJD,C,X,SECDIF,TLIGHT,DIS,DT,
     .     FRLIMB,RCIO,PEB,VEB,PSB,VSB,POG,VOG,POB,VOB,
     .     POS1,VEL1,POS2,POS3,POS4,POS5,POS6,POS7,POS8,
     .     PX,PY,PZ,RA,DEC,RVS,RVD,RV,DABS,DSQRT
      CHARACTER*(*) OBJECT

      DIMENSION STAR(*), OBSERV(6), SKYPOS(*), 
     .     PEB(3), VEB(3), PSB(3), VSB(3),
     .     POG(3), VOG(3), POB(3), VOB(3),
     .     POS1(3), VEL1(3), POS2(3), POS3(3), POS4(3), POS5(3),
     .     POS6(3), POS7(3), POS8(3),
     .     PX(3), PY(3), PZ(3), RVS(3), RVD(3)

      SAVE

*     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
      DATA T0 / 2451545.00000000D0 /
      DATA TLAST1, TLAST2 / 0.D0, 0.D0 /,   NTIMES / 0 /

   3  FORMAT ( ' PLACE: CANNOT OBTAIN COORDINATES OF ', A, ' AT JD ',
     .     F10.1 )
   4  FORMAT ( ' PLACE: WILL NOT PROCESS EARTH AS OBSERVED OBJECT ',
     .     'EXCEPT WHEN LOCATN=2' )

* --- GET CONSTANTS, FIRST TIME ONLY ----------------------------------

      NTIMES = NTIMES + 1

      IF ( NTIMES .EQ. 1 ) THEN
          IEARTH = IDSS ( 'EARTH' )
          ISUN = IDSS ( 'SUN' )
*         GET C, THE SPEED OF LIGHT IN AU/DAY
          CALL ASTCON ( 'C(AU/DAY)', 1.D0,   C )
      END IF

* --- CHECK ON EARTH AS AN OBSERVED OBJECT ----------------------------

      IF ( OBJECT .EQ. 'EARTH' .AND. LOCATN .NE. 2 ) THEN
          WRITE ( *, 4 )
          GO TO 70
      END IF

* --- GET POSITION AND VELOCITY OF EARTH (GEOCENTER) AND SUN ----------

      IF ( DABS ( TJD - TLAST1 ) .GT. 1.D-8 ) THEN

*         COMPUTE TDBJD, THE TDB JULIAN DATE CORRESPONDING TO TTJD
          TTJD = TJD
          TDBJD = TJD
          CALL TIMES ( TDBJD, X,   SECDIF )
          TDBJD = TTJD + SECDIF / 86400.D0

*         GET POSITION AND VELOCITY OF THE EARTH WRT BARYCENTER OF
*         SOLAR SYSTEM, IN ICRS
          CALL SOLSYS ( TDBJD, IEARTH, 0,   PEB, VEB, IERR )
          IF ( IERR .NE. 0 ) THEN
              WRITE ( *, 3 ) 'EARTH', TJD
              GO TO 70
          END IF

*         GET POSITION AND VELOCITY OF THE SUN WRT BARYCENTER OF
*         SOLAR SYSTEM, IN ICRS
          CALL SOLSYS ( TDBJD, ISUN, 0,   PSB, VSB, IERR )
          IF ( IERR .NE. 0 ) THEN
              WRITE ( *, 3 ) 'SUN', TJD
              GO TO 70
          END IF

          TLAST1 = TJD

      END IF

* --- GET POSITION AND VELOCITY OF OBSERVER ---------------------------

      IF ( LOCATN .EQ. 1 .OR. LOCATN .EQ. 2 ) THEN

*         FOR TOPOCENTRIC PLACE, GET GEOCENTRIC POSITION AND VELOCITY
*         VECTORS OF OBSERVER
          CALL GEOPOS ( TTJD, LOCATN, OBSERV,   POG, VOG )
          LOC = 1

      ELSE

*         FOR GEOCENTRIC PLACE, THERE IS NOTHING TO DO
          DO 25 J = 1, 3
             POG(J) = 0.D0
             VOG(J) = 0.D0
  25      CONTINUE
          LOC = 0

      END IF

*     COMPUTE POSITION AND VELOCITY OF OBSERVER WRT BARYCENTER OF
*     SOLAR SYSTEM (GALILEAN TRANSFORMATION FROM GCRS to BCRS)
      DO 30 J=1,3
          POB(J) = PEB(J) + POG(J)
          VOB(J) = VEB(J) + VOG(J)
  30  CONTINUE

* --- FIND GEOMETRIC POSITION OF OBSERVED OBJECT ----------------------

      IF ( OBJECT .EQ. 'STAR' .OR. OBJECT .EQ. ' ' .OR.
     .     OBJECT(1:1) .EQ. '*' ) THEN

*         OBSERVED OBJECT IS STAR
          IDBODY = -9999

*         GET POSITION OF STAR UPDATED FOR ITS SPACE MOTION
          CALL VECTRS ( STAR(1), STAR(2), STAR(3), STAR(4),
     .                  STAR(5), STAR(6),   POS1, VEL1 )
          CALL DLIGHT ( POS1, POB,   DT )
          CALL PROPMO ( T0, POS1, VEL1, TDBJD + DT,   POS2 )
*         GET POSITION OF STAR WRT OBSERVER (CORRECTED FOR PARALLAX)
          CALL GEOCEN ( POS2, POB,   POS3, TLIGHT )
          DIS = 0.D0

      ELSE

*         OBSERVED OBJECT IS SOLAR SYSTEM BODY

*         GET ID NUMBER OF BODY
          IF ( OBJECT(1:1) .EQ. '=' ) THEN
              READ ( OBJECT, '(1X,I3)' ) IDBODY
          ELSE
              IDBODY = IDSS ( OBJECT )
              IF ( IDBODY .EQ. -9999 ) THEN
                  WRITE ( *, 3 ) OBJECT, TJD
                  GO TO 70
              END IF
          END IF

*         GET POSITION OF BODY WRT BARYCENTER OF SOLAR SYSTEM
          CALL SOLSYS ( TDBJD, IDBODY, 0,   POS1, VEL1, IERR )
          IF ( IERR .NE. 0 ) THEN
              WRITE ( *, 3 ) OBJECT, TJD
              GO TO 70
          END IF

*         GET POSITION OF BODY WRT OBSERVER, AND TRUE (EUCLIDIAN)
*         DISTANCE
          CALL GEOCEN ( POS1, POB,   POS2, TLIGHT )
          DIS = TLIGHT * C

*         GET POSITION OF BODY WRT OBSERVER, ANTEDATED FOR LIGHT-TIME
          CALL LITTIM ( TDBJD, IDBODY, POB, TLIGHT,   POS3, TLIGHT )

      END IF

* --- APPLY GRAVITATIONAL DEFLECTION OF LIGHT AND ABERRATION ----------

      IF ( ICOORD .EQ. 3 ) THEN

*         THESE CALCULATIONS ARE SKIPPED FOR ASTROMETRIC PLACE
          DO 40 J = 1, 3
              POS5(J) = POS3(J)
  40      CONTINUE

      ELSE

*         VARIABLE LOC DETERMINES WHETHER EARTH DEFLECTION IS INCLUDED
          IF ( LOC .EQ. 1 ) THEN
              CALL LIMANG ( POS3, POG,   X, FRLIMB )
              IF ( FRLIMB .LT. 0.8D0 ) LOC = 0
          END IF

*         COMPUTE GRAVITATIONAL DEFLECTION AND ABERRATION
          CALL GRVDEF ( TDBJD, LOC, POS3, POB,   POS4 )
          CALL ABERAT ( POS4, VOB, TLIGHT,   POS5 )
*         POSITION VECTOR IS NOW IN GCRS 

      END IF

* --- TRANSFORM, IF NECESSARY, TO OUTPUT COORDINATE SYSTEM ------------

      IF ( ICOORD .EQ. 1 ) THEN

*         TRANSFORM TO EQUATOR AND EQUINOX OF DATE
          CALL FRAME  ( POS5, 1,   POS6 )
          CALL PRECES ( T0, POS6, TDBJD,   POS7 )
          CALL NUTATE ( TDBJD, POS7,   POS8 )

      ELSE IF ( ICOORD .EQ. 2 ) THEN

*         TRANSFORM TO EQUATOR AND CIO OF DATE
          IF ( DABS ( TDBJD - TLAST2 ) .GT. 1.D-8 ) THEN
*             OBTAIN THE BASIS VECTORS, IN THE GCRS, OF THE CELESTIAL
*             INTERMEDIATE SYSTEM
              CALL CIOLOC ( TDBJD,   RCIO, KCIO )
              CALL CIOBAS ( TDBJD, RCIO, KCIO,   PX, PY, PZ )
              TLAST2 = TDBJD
          END IF
*         TRANSFORM POSITION VECTOR TO CELESTIAL INTERMEDIATE SYSTEM
          POS8(1) = PX(1) * POS5(1) + PX(2) * POS5(2) + PX(3) * POS5(3)
          POS8(2) = PY(1) * POS5(1) + PY(2) * POS5(2) + PY(3) * POS5(3)
          POS8(3) = PZ(1) * POS5(1) + PZ(2) * POS5(2) + PZ(3) * POS5(3)

      ELSE

*         NO TRANSFORMATION -- KEEP COORDINATES IN GCRS
*         (OR ICRS FOR ASTROMETRIC PLACE)
          DO 50 J = 1, 3
              POS8(J) = POS5(J)
  50      CONTINUE

      END IF
      
* --- GET RADIAL VELOCITY ---------------------------------------------      

*     SET UP STAR DATA, IF APPLICABLE      
      IF ( IDBODY .EQ. -9999 ) THEN
          RVS(1) = STAR(1)
          RVS(2) = STAR(2)
          RVS(3) = STAR(6)
          IF ( STAR(5) .LE. 0.D0 ) THEN
              VEL1(1) = 0.D0
              VEL1(2) = 0.D0
              VEL1(3) = 0.D0 
          END IF
      ELSE
          RVS(1) = 0.D0
          RVS(2) = 0.D0
          RVS(3) = 0.D0
      END IF

*     COMPUTE DISTANCES: OBSERVER-GEOCENTER, OBSERVER-SUN, OBJECT-SUN
      RVD(1) = DSQRT ( ( POB(1)  - PEB(1) ) ** 2
     .               + ( POB(2)  - PEB(2) ) ** 2
     .               + ( POB(3)  - PEB(3) ) ** 2 )
      RVD(2) = DSQRT ( ( POB(1)  - PSB(1) ) ** 2
     .               + ( POB(2)  - PSB(2) ) ** 2
     .               + ( POB(3)  - PSB(3) ) ** 2 )
      RVD(3) = DSQRT ( ( POS1(1) - PSB(1) ) ** 2
     .               + ( POS1(2) - PSB(2) ) ** 2
     .               + ( POS1(3) - PSB(3) ) ** 2 )
      
      CALL RADVL ( POS3, VEL1, VOB, RVS, RVD,   RV )
      
* --- FINISH UP -------------------------------------------------------

      CALL ANGLES ( POS8,   RA, DEC )

      X = DSQRT ( POS8(1)**2 + POS8(2)**2 + POS8(3)**2 )

      DO 60 J = 1, 3
          SKYPOS(J) = POS8(J) / X
  60  CONTINUE
  
      SKYPOS(4) = RA
      SKYPOS(5) = DEC
      SKYPOS(6) = DIS
      SKYPOS(7) = RV
      
      CALL SETVEC ( POS8 )  
  
  70  RETURN

      END



      SUBROUTINE PLACES
*
*     THE ENTRIES TO THIS SUBROUTINE PROVIDE 'FRONT ENDS' TO
*     SUBROUTINE PLACE, TAILORED TO SPECIFIC PLACE TYPES.  THEY
*     PROVIDE COMPATIBILITY WITH PREVIOUSLY SUPPORTED CALLING
*     SEQUENCES.
*
*
      IMPLICIT NONE
      INTEGER L,N,LOCATN,ICOORD
      DOUBLE PRECISION TJD,RAI,DECI,PMRA,PMDEC,PARLAX,RADVEL,
     .     UJD,GLON,GLAT,HT,RA,DEC,DIS,
     .     TTJD,GAST,DELTAT,STAR,OBSERV,SKYPOS,DMOD
      CHARACTER*4 OBJECT

      DIMENSION STAR(20), OBSERV(6), SKYPOS(10)

      SAVE   

      DATA STAR, OBSERV, SKYPOS / 36 * 0.D0 /
   
  
* --- GEOCENTRIC PLACES OF STARS --------------------------------------
*
*     ARGUMENTS COMMON TO THESE ENTRIES:
*
*          TJD    = TT JULIAN DATE FOR APPARENT PLACE (IN)
*          N      = BODY IDENTIFICATION NUMBER FOR THE EARTH (IN)
*                   (NO LONGER USED)
*          RAI    = ICRS RIGHT ASCENSION IN HOURS (IN)
*          DECI   = ICRS DECLINATION IN DEGREES (IN)
*          PMRA   = ICRS PROPER MOTION IN RA IN MILLIARCSECONDS/YEAR
*                   (IN)
*          PMDEC  = ICRS PROPER MOTION IN DEC IN MILLIARCSECONDS/YEAR
*                   (IN)
*          PARLAX = PARALLAX IN MILLIARCSECONDS (IN)
*          RADVEL = RADIAL VELOCITY IN KILOMETERS/SECOND (IN)
*          RA     = APPARENT RIGHT ASCENSION IN HOURS (OUT)
*          DEC    = APPARENT DECLINATION IN DEGREES (OUT)
*
*     NOTE: COORDINATE SYSTEM FOR OUTPUT RA AND DEC IS GCRS FOR
*     VPSTAR, ICRS FOR ASSTAR, AND EQUATOR AND EQUINOX OF DATE FOR
*     APSTAR.
*
*
*     THIS ENTRY PROVIDES THE APPARENT PLACE OF A STAR.
      ENTRY APSTAR (TJD,N,RAI,DECI,PMRA,PMDEC,PARLAX,RADVEL, RA,DEC)
      LOCATN = 0
      ICOORD = 1
      GO TO 10

*     THIS ENTRY PROVIDES THE VIRTUAL PLACE OF A STAR.
      ENTRY VPSTAR (TJD,N,RAI,DECI,PMRA,PMDEC,PARLAX,RADVEL, RA,DEC)
      LOCATN = 0
      ICOORD = 0
      GO TO 10

*     THIS ENTRY PROVIDES THE ASTROMETRIC PLACE OF A STAR.
      ENTRY ASSTAR (TJD,N,RAI,DECI,PMRA,PMDEC,PARLAX,RADVEL, RA,DEC)
      LOCATN = 0
      ICOORD = 3

  10  TTJD = TJD
      OBJECT = 'STAR'
      STAR(1) = RAI
      STAR(2) = DECI
      STAR(3) = PMRA
      STAR(4) = PMDEC
      STAR(5) = PARLAX
      STAR(6) = RADVEL

      CALL PLACE (TTJD,OBJECT,LOCATN,ICOORD,STAR,OBSERV, SKYPOS)
      
      RA  = SKYPOS(4)
      DEC = SKYPOS(5)
      
      RETURN


* --- TOPOCENTRIC PLACES OF STARS -------------------------------------
*
*     EACH OF THESE ENTRIES CAN BE CALLED ONLY AFTER A CALL TO THE
*     CORRESPONDING GEOCENTRIC ENTRY, WHICH SUPPLIES SOME REQUIRED DATA.
*     ARGUMENTS COMMON TO THESE ENTRIES:
*
*          UJD    = UT1 JULIAN DATE FOR TOPOCENTRIC PLACE (IN)
*          GLON   = GEODETIC (ITRS) LONGITUDE (EAST +) OF OBSERVER
*                   IN DEGREES (IN)
*          GLAT   = GEODETIC (ITRS) LATITUDE (NORTH +) OF OBSERVER
*                   IN DEGREES (IN)
*          HT     = HEIGHT OF OBSERVER IN METERS (IN)
*          RA     = TOPOCENTRIC RIGHT ASCENSION IN HOURS (OUT)
*          DEC    = TOPOCENTRIC DECLINATION IN DEGREES (OUT)
*
*     NOTE 1: COORDINATE SYSTEM FOR OUTPUT RA AND DEC IS 'LOCAL GCRS'
*     FOR LPSTAR, AND EQUATOR AND EQUINOX OF DATE FOR TPSTAR.
*
*     NOTE 2: UJD CAN ALSO BE GREENWICH APPARENT SIDEREAL TIME IN HOURS,
*     EQUIVALENT TO UT1 JULIAN DATE, BUT THIS OPTION WILL NOT BE
*     SUPPORTED INDEFINITELY.  ADVISE USING UJD = UT1 JULIAN DATE ONLY.
*
*
*     THIS ENTRY PROVIDES THE TOPOCENTRIC PLACE OF A STAR.
      ENTRY TPSTAR (UJD,GLON,GLAT,HT, RA,DEC)
      LOCATN = 1
      ICOORD = 1
      GO TO 20

*     THIS ENTRY PROVIDES THE LOCAL PLACE OF A STAR.
      ENTRY LPSTAR (UJD,GLON,GLAT,HT, RA,DEC)
      LOCATN = 1
      ICOORD = 0

  20  IF ( UJD .GT. 100.D0 ) THEN
          DELTAT = ( TTJD - UJD ) * 86400.D0
      ELSE
          GAST = DMOD ( UJD, 24.D0 )
          IF ( GAST .LT. 0.D0 ) GAST = GAST + 24.D0
          CALL PLACST ( GAST )
          DELTAT = 0.D0
      END IF
      OBSERV(1) = GLON
      OBSERV(2) = GLAT
      OBSERV(3) = HT
      OBSERV(4) = DELTAT

      CALL PLACE (TTJD,OBJECT,LOCATN,ICOORD,STAR,OBSERV, SKYPOS)
      
      RA  = SKYPOS(4)
      DEC = SKYPOS(5)
      
      RETURN


* --- GEOCENTRIC PLACES OF PLANETS ------------------------------------
*
*     ARGUMENTS COMMON TO THESE ENTRIES:
*
*          TJD    = TT JULIAN DATE FOR APPARENT PLACE (IN)
*          L      = BODY IDENTIFICATION NUMBER FOR DESIRED PLANET (IN)
*          N      = BODY IDENTIFICATION NUMBER FOR THE EARTH (IN)
*                   (NO LONGER USED)
*          RA     = APPARENT RIGHT ASCENSION IN HOURS (OUT)
*          DEC    = APPARENT DECLINATION IN DEGREES (OUT)
*          DIS    = TRUE DISTANCE FROM EARTH TO PLANET IN AU (OUT)
*
*     NOTE: COORDINATE SYSTEM FOR OUTPUT RA AND DEC IS GCRS FOR
*     VPPLAN, ICRS FOR ASPLAN, AND EQUATOR AND EQUINOX OF DATE FOR
*     APPLAN.
*
*     NOTE: 'PLANET' IS USED GENERICALLY FOR ANY SOLAR SYSTEM BODY. 
*
*     THIS ENTRY PROVIDES THE APPARENT PLACE OF A PLANET.
      ENTRY APPLAN (TJD,L,N, RA,DEC,DIS)
      LOCATN = 0
      ICOORD = 1
      GO TO 30

*     THIS ENTRY PROVIDES THE VIRTUAL PLACE OF A PLANET.
      ENTRY VPPLAN (TJD,L,N, RA,DEC,DIS)
      LOCATN = 0
      ICOORD = 0
      GO TO 30

*     THIS ENTRY PROVIDES THE ASTROMETRIC PLACE OF A PLANET.
      ENTRY ASPLAN (TJD,L,N, RA,DEC,DIS)
      LOCATN = 0
      ICOORD = 3

  30  TTJD = TJD
      WRITE ( OBJECT, '(''='',I3)' ) L

      CALL PLACE (TTJD,OBJECT,LOCATN,ICOORD,STAR,OBSERV, SKYPOS)
      
      RA  = SKYPOS(4)
      DEC = SKYPOS(5)
      DIS = SKYPOS(6)

      RETURN


* --- TOPOCENTRIC PLACES OF PLANETS -----------------------------------
*
*     EACH OF THESE ENTRIES CAN BE CALLED ONLY AFTER A CALL TO THE
*     CORRESPONDING GEOCENTRIC ENTRY, WHICH SUPPLIES SOME REQUIRED DATA.
*     ARGUMENTS COMMON TO THESE ENTRIES:
*
*          UJD    = UT1 JULIAN DATE FOR TOPOCENTRIC PLACE (IN)
*          GLON   = GEODETIC (ITRS) LONGITUDE (EAST +) OF OBSERVER
*                   IN DEGREES (IN)
*          GLAT   = GEODETIC (ITRS) LATITUDE (NORTH +) OF OBSERVER
*                   IN DEGREES (IN)
*          HT     = HEIGHT OF OBSERVER IN METERS (IN)
*          RA     = TOPOCENTRIC RIGHT ASCENSION IN HOURS (OUT)
*          DEC    = TOPOCENTRIC DECLINATION IN DEGREES (OUT)
*          DIS    = TRUE DISTANCE FROM OBSERVER TO PLANET IN AU (OUT)
*
*     NOTE 1: COORDINATE SYSTEM FOR OUTPUT RA AND DEC IS 'LOCAL GCRS'
*     FOR LPPLAN, AND EQUATOR AND EQUINOX OF DATE FOR TPPLAN.
*
*     NOTE 2: UJD CAN ALSO BE GREENWICH APPARENT SIDEREAL TIME IN HOURS,
*     EQUIVALENT TO UT1 JULIAN DATE, BUT THIS OPTION WILL NOT BE
*     SUPPORTED INDEFINITELY.  ADVISE USING UJD = UT1 JULIAN DATE ONLY.
*
*
*     THIS ENTRY PROVIDES THE TOPOCENTRIC PLACE OF A PLANET.
      ENTRY TPPLAN (UJD,GLON,GLAT,HT, RA,DEC,DIS)
      LOCATN = 1
      ICOORD = 1
      GO TO 40

*     THIS ENTRY PROVIDES THE LOCAL PLACE OF A PLANET.
      ENTRY LPPLAN (UJD,GLON,GLAT,HT, RA,DEC,DIS)
      LOCATN = 1
      ICOORD = 0

  40  IF ( UJD .GT. 100.D0 ) THEN
          DELTAT = ( TTJD - UJD ) * 86400.D0
      ELSE
          GAST = DMOD ( UJD, 24.D0 )
          IF ( GAST .LT. 0.D0 ) GAST = GAST + 24.D0
          CALL PLACST ( GAST )
          DELTAT = 0.D0
      END IF
      OBSERV(1) = GLON
      OBSERV(2) = GLAT
      OBSERV(3) = HT
      OBSERV(4) = DELTAT

      CALL PLACE (TTJD,OBJECT,LOCATN,ICOORD,STAR,OBSERV, SKYPOS)
      
      RA  = SKYPOS(4)
      DEC = SKYPOS(5)
      DIS = SKYPOS(6)

      RETURN

      END



      SUBROUTINE MPSTAR (TJD,N,RA,DEC, RAI,DECI)
*
*     THIS SUBROUTINE COMPUTES THE ICRS POSITION OF A STAR,
*     GIVEN ITS APPARENT PLACE AT DATE TJD.  PROPER MOTION, PARALLAX,
*     AND RADIAL VELOCITY ARE ASSUMED TO BE ZERO.
*
*          TJD    = TT JULIAN DATE OF APPARENT PLACE (IN)
*          N      = BODY IDENTIFICATION NUMBER FOR THE EARTH (IN)
*                   (NO LONGER USED)
*          RA     = APPARENT RIGHT ASCENSION IN HOURS, REFERRED TO
*                   TRUE EQUATOR AND EQUINOX OF DATE (IN)
*          DEC    = APPARENT DECLINATION IN DEGREES, REFERRED TO
*                   TRUE EQUATOR OF DATE (IN)
*          RAI    = ICRS RIGHT ASCENSION IN HOURS (OUT)
*          DECI   = ICRS DECLINATION IN DEGREES (OUT)
*
*
      DOUBLE PRECISION TJD,RA,DEC,RAI,DECI,T0,T1,RAINEW,DCINEW,
     .     RAIOLD,DCIOLD,STAR,OBSERV,SKYPOS,R,D,P,V,DELRA,DELDEC,
     .     DABS
      DIMENSION STAR(20), OBSERV(6), SKYPOS(10), P(3), V(3)
      SAVE

*     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
      DATA T0 / 2451545.00000000D0 /
      DATA STAR, OBSERV, SKYPOS / 36 * 0.D0 /

   3  FORMAT ( ' MPSTAR: NO CONVERGENCE AT COORDINATES ',
     .     F9.5, 1X, SP,F9.4, ' AT JD ', SS, F10.1 )

      T1 = TJD

*     GET INITIAL APPROXIMATION
      ITER = 0
      CALL VECTRS (RA,DEC,0.D0,0.D0,0.D0,0.D0,P,V)
      CALL PRECES (T1,P,T0,V)
      CALL ANGLES (V,RAINEW,DCINEW)

*     ITERATIVELY FIND ICRS COORDINATES THAT PRODUCE INPUT
*     APPARENT PLACE OF STAR AT DATE TJD
   20 ITER = ITER + 1
      RAIOLD = RAINEW
      DCIOLD = DCINEW
      STAR(1) = RAIOLD
      STAR(2) = DCIOLD
      STAR(3) = 0.D0
      STAR(4) = 0.D0
      STAR(5) = 0.D0
      STAR(6) = 0.D0
      CALL PLACE (T1,'STAR',0,1,STAR,OBSERV, SKYPOS)
      R = SKYPOS(4)
      D = SKYPOS(5)
      DELRA = R - RA
      DELDEC = D - DEC
      IF (DELRA.LT.-12.D0) DELRA = DELRA + 24.D0
      IF (DELRA.GT.+12.D0) DELRA = DELRA - 24.D0
      RAINEW = RAIOLD - DELRA
      DCINEW = DCIOLD - DELDEC
      IF (ITER.GT.30) THEN
          WRITE ( *, 3 ) RA, DEC, TJD
          GO TO 40
      END IF
      IF (DABS(RAINEW-RAIOLD).GT.1.D-12) GO TO 20
      IF (DABS(DCINEW-DCIOLD).GT.1.D-11) GO TO 20

   40 RAI = RAINEW
      DECI = DCINEW
      IF (RAI.LT. 0.D0) RAI = RAI + 24.D0
      IF (RAI.GE.24.D0) RAI = RAI - 24.D0

*     STORE COMPUTED POSITION VECTOR FOR POSSIBLE LATER RETRIEVAL
      CALL VECTRS (RAI,DECI,0.D0,0.D0,0.D0,0.D0,P,V)
      CALL SETVEC (P)

   50 RETURN

      END



      SUBROUTINE SIDTIM ( TJDH, TJDL, K,   GST )
*
*     THIS SUBROUTINE COMPUTES THE GREENWICH SIDEREAL TIME
*     (EITHER MEAN OR APPARENT) AT JULIAN DATE TJDH + TJDL.
*
*          TJDH   = UT1 JULIAN DATE, HIGH-ORDER PART (IN)
*          TJDL   = UT1 JULIAN DATE, LOW-ORDER PART (IN)
*                   THE JULIAN DATE MAY BE SPLIT AT ANY POINT, BUT
*                   FOR HIGHEST PRECISION, SET TJDH TO BE THE INTEGRAL
*                   PART OF THE JULIAN DATE, AND SET TJDL TO BE THE
*                   FRACTIONAL PART
*          K      = TIME SELECTION CODE (IN)
*                   SET K=0 FOR GREENWICH MEAN SIDEREAL TIME
*                   SET K=1 FOR GREENWICH APPARENT SIDEREAL TIME
*          GST    = GREENWICH (MEAN OR APPARENT) SIDEREAL TIME
*                   IN HOURS (OUT)
*
*     NOTE:  SEE ALSO SUBROUTINE SETDT TO SET THE VALUE OF DELTA-T
*     (DELTA-T = TT - UT1) TO BE USED HERE.
*
*
      DOUBLE PRECISION TJDH,TJDL,GST,PI,DEGCON,DELTAT,
     .     T0,UTJD,TTJD,TDBJD,SECDIF,A,THETA,RCIO,
     .     UNITX,W1,W2,X,Y,Z,EQ,HAEQ,EE,DMOD,DATAN2
      DIMENSION UNITX(3), W1(3), W2(3), X(3), Y(3), Z(3), EQ(3)
      SAVE 

      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( DEGCON = 180.D0 / PI           )

*     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
      DATA T0 / 2451545.00000000D0 /
      DATA UNITX / 1.D0, 0.D0, 0.D0 /

   3  FORMAT ( ' SIDTIM ERROR: CANNOT RETURN SIDEREAL TIME FOR ',
     .     'JD ', F10.1 )

      CALL GETDT ( DELTAT )

*     TIME ARGUMENT FOR PRECESSION AND NUTATION COMPONENTS OF SIDEREAL
*     TIME IS TDB
      UTJD = TJDH + TJDL
      TTJD = UTJD + DELTAT
      TDBJD = TTJD
      CALL TIMES ( TDBJD, A,   SECDIF )
      TDBJD = TTJD + SECDIF / 86400.D0

*     GET METHOD/ACCURACY MODE
      CALL GETMOD ( MODE )

      IF ( MODE .GE. 2 ) THEN
*          EQUINOX-BASED MODE
*          SEE USNO CIRCULAR 179, SECTION 2.6.2   
     
*         GET -1 TIMES THE MEAN OR TRUE RIGHT ASCENSION OF THE CIO     
          CALL EQXRA ( TDBJD, K,   RCIO )
*         GET EARTH ROTATION ANGLE
          CALL EROT ( TJDH, TJDL,   THETA )
*         COMBINE TO OBTAIN SIDEREAL TIME       
          GST = DMOD ( THETA / 15.D0 - RCIO, 24.D0 )
          IF ( GST .LT. 0.D0 ) GST = GST + 24.D0

      ELSE
*         CIO-BASED MODE
*         SEE USNO CIRCULAR 179, SECTION 6.5.4

*         GET EARTH ROTATION ANGLE
          CALL EROT ( TJDH, TJDL,   THETA )
*         OBTAIN THE BASIS VECTORS, IN THE GCRS, OF THE CELESTIAL
*         INTERMEDIATE SYSTEM
          CALL CIOLOC ( TDBJD,   RCIO, KCIO )
          IF ( RCIO .EQ. 99.D0 ) THEN
              WRITE ( *, 3 ) TDBJD
              GST = 99.D0
              GO TO 50
          END IF
          CALL CIOBAS ( TDBJD, RCIO, KCIO,   X, Y, Z )
*         COMPUTE THE DIRECTION OF THE TRUE EQUINOX IN THE GCRS
          CALL NUTATE ( -TDBJD, UNITX,   W1 )
          CALL PRECES ( TDBJD, W1, T0,   W2 )
          CALL FRAME ( W2, -1,    EQ )
*         COMPUTE THE HOUR ANGLE OF THE EQUINOX WRT THE TIO MERIDIAN
*         (NEAR GREENWICH, BUT PASSES THROUGH THE CIP AND TIO)
          HAEQ = THETA - DATAN2 ( EQ(1)*Y(1) + EQ(2)*Y(2) + EQ(3)*Y(3),
     .                            EQ(1)*X(1) + EQ(2)*X(2) + EQ(3)*X(3) )
     .                   * DEGCON

*         FOR MEAN SIDEREAL TIME, OBTAIN THE EQUATION OF THE EQUINOXES
*         AND SUBTRACT IT
          IF ( K .EQ. 0 ) THEN
              CALL ETILT ( TDBJD,   A, A, EE, A, A )
              HAEQ = HAEQ - EE / 240.D0
          END IF

          HAEQ = DMOD ( HAEQ, 360.D0 ) / 15.D0
          IF ( HAEQ .LT. 0.D0 ) HAEQ = HAEQ + 24.D0
          GST = HAEQ

      END IF

  50  RETURN

      END



      SUBROUTINE CIORA ( TJD,   RACIO )
*
*     THIS SUBROUTINE COMPUTES THE TRUE RIGHT ASCENSION OF THE CELESTIAL
*     INTERMEDIATE ORIGIN (CIO) AT A GIVEN TT JULIAN DATE.  THIS IS
*     -(EQUATION OF THE ORIGINS).
*
*          TJD    = TT JULIAN DATE (IN)
*          RACIO  = RIGHT ASCENSION OF THE CIO, WITH RESPECT TO THE
*                   TRUE EQUINOX OF DATE, IN HOURS (+ OR -) (OUT)
*                 
*
      DOUBLE PRECISION TJD,RACIO,PI,DEGCON,T0,T,SECDIF,TDBJD,RCIO,
     .     UNITX,W1,W2,X,Y,Z,EQ,AZ,DATAN2
      DIMENSION UNITX(3), W1(3), W2(3), X(3), Y(3), Z(3), EQ(3)
      SAVE

      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( DEGCON = 180.D0 / PI           )

*     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
      DATA T0 / 2451545.00000000D0 /
      DATA UNITX / 1.D0, 0.D0, 0.D0 /

   3  FORMAT ( ' CIORA ERROR: CANNOT RETURN CIO RA VALUE FOR JD ',
     .     F10.1 )

*     TDBJD IS THE TDB JULIAN DATE
      TDBJD = TJD
      CALL TIMES ( TDBJD, T,   SECDIF )
      TDBJD = TJD + SECDIF / 86400.D0

*     OBTAIN THE BASIS VECTORS, IN THE GCRS, FOR THE CELESTIAL
*     INTERMEDIATE SYSTEM DEFINED BY THE CIP (IN THE Z DIRECTION) AND
*     THE CIO (IN THE X DIRECTION)
      CALL CIOLOC ( TDBJD,   RCIO, KCIO )
      IF ( RCIO .EQ. 99.D0 ) THEN
          WRITE ( *, 3 ) TDBJD
          RACIO = 99.D0
          GO TO 50
      END IF
      CALL CIOBAS ( TDBJD, RCIO, KCIO,   X, Y, Z )

*     COMPUTE THE DIRECTION OF THE TRUE EQUINOX IN THE GCRS
      CALL NUTATE ( -TDBJD, UNITX,   W1 )
      CALL PRECES ( TDBJD, W1, T0,   W2 )
      CALL FRAME ( W2, -1,    EQ )

*     COMPUTE THE INTERMEDIATE RA OF THE TRUE EQUINOX (EQUATION OF
*     THE ORIGINS)
      AZ = DATAN2 ( EQ(1)*Y(1) + EQ(2)*Y(2) + EQ(3)*Y(3),
     .              EQ(1)*X(1) + EQ(2)*X(2) + EQ(3)*X(3) ) * DEGCON

*     THE RA OF THE CIO IS MINUS THIS COORDINATE
      RACIO = -AZ / 15.D0

  50  RETURN

      END



      SUBROUTINE TERCEL ( TJDH, TJDL, XP, YP, VEC1,   VEC2 )
*
*     THIS SUBROUTINE ROTATES A VECTOR FROM THE TERRESTRIAL TO THE
*     CELESTIAL SYSTEM.  SPECIFICALLY, IT TRANSFORMS A VECTOR IN THE
*     ITRS (A ROTATING EARTH-FIXED SYSTEM) TO THE GCRS (A LOCAL
*     SPACE-FIXED SYSTEM) BY APPLYING ROTATIONS FOR POLAR MOTION,
*     EARTH ROTATION, NUTATION, PRECESSION, AND THE DYNAMICAL-TO-GCRS
*     FRAME TIE.
*
*          TJDH   = UT1 JULIAN DATE, HIGH-ORDER PART (IN)
*          TJDL   = UT1 JULIAN DATE, LOW-ORDER PART (IN)
*                   THE JULIAN DATE MAY BE SPLIT AT ANY POINT, BUT
*                   FOR HIGHEST PRECISION, SET TJDH TO BE THE INTEGRAL
*                   PART OF THE JULIAN DATE, AND SET TJDL TO BE THE
*                   FRACTIONAL PART
*          XP     = CONVENTIONALLY-DEFINED X COORDINATE OF CELESTIAL
*                   INTERMEDIATE POLE WITH RESPECT TO ITRS POLE,
*                   IN ARCSECONDS (IN)
*          YP     = CONVENTIONALLY-DEFINED Y COORDINATE OF CELESTIAL
*                   INTERMEDIATE POLE WITH RESPECT TO ITRS POLE,
*                   IN ARCSECONDS (IN)
*          VEC1   = POSITION VECTOR, GEOCENTRIC EQUATORIAL RECTANGULAR
*                   COORDINATES, REFERRED TO ITRS AXES (TERRESTRIAL
*                   SYSTEM) (IN)
*          VEC2   = POSITION VECTOR, GEOCENTRIC EQUATORIAL RECTANGULAR
*                   COORDINATES, REFERRED TO GCRS AXES (CELESTIAL
*                   SYSTEM) (OUT)
*
*     NOTE 1:  SET XP=YP=0.D0 TO ELIMINATE POLAR MOTION ROTATION.
*
*     NOTE 2:  SEE ALSO SUBROUTINE SETDT TO SET THE VALUE OF DELTA-T
*     (DELTA-T = TT - UT1) TO BE USED HERE.
*
*     NOTE 3:  BOTH TJDH AND TJDL SHOULD BE NON-NEGATIVE FOR NORMAL USE
*     (TJDL=0.D0 IS OK).  A NEGATIVE VALUE OF TJDH IS USED TO INVOKE A
*     SPECIAL OPTION WHERE THE OUTPUT VECTOR IS PRODUCED WITH RESPECT
*     TO THE EQUATOR AND EQUINOX OF DATE, AND THE DATE FOR WHICH THE
*     TRANSFORMATION APPLIES IS TAKEN FROM TJDL ONLY.  THIS OPTION
*     WORKS ONLY IN 'EQUINOX' MODE.
*
*     NOTE 4: INPUT PARAMETERS XP, YP WERE XPOLE, YPOLE IN NOVAS F3.0.
*     THE NAMES WERE CHANGED FOR CONSISTANCY THROUGHOUT NOVAS AND WITH
*     IERS CONVENTIONS.
*
      DOUBLE PRECISION TJDH,TJDL,XP,YP,VEC1,VEC2,
     .     T0,DELTAT,UTJDH,UTJDL,UTJD,TTJD,TDBJD,T,SECDIF,
     .     GAST,RCIO,THETA,V1,V2,V3,V4,X,Y,Z
      DIMENSION VEC1(3), VEC2(3), V1(3), V2(3), V3(3), V4(3),
     .     X(3), Y(3), Z(3)

*     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
      DATA T0 / 2451545.00000000D0 /

      CALL GETDT ( DELTAT )

      IF ( TJDH .GE. 0.D0 ) THEN
          UTJDH = TJDH
          UTJDL = TJDL
      ELSE
          UTJDH = TJDL
          UTJDL = 0.D0
      END IF
      UTJD = UTJDH + UTJDL

*     TIME ARGUMENT FOR PRECESSION AND NUTATION IS TDB
      TTJD = UTJD + DELTAT
      TDBJD = TTJD
      CALL TIMES ( TDBJD, T,   SECDIF )
      TDBJD = TTJD + SECDIF / 86400.D0

*     GET METHOD/ACCURACY MODE
      CALL GETMOD ( MODE )

      IF ( MODE .GE. 2 ) THEN
*         'EQUINOX' MODE

*         APPLY POLAR MOTION
          IF ( XP .EQ. 0.D0 .AND. YP .EQ. 0.D0 ) THEN
              V1(1) = VEC1(1)
              V1(2) = VEC1(2)
              V1(3) = VEC1(3)
          ELSE
              CALL WOBBLE ( TDBJD, XP, YP, VEC1,   V1 )
          END IF

*         APPLY EARTH ROTATION
          CALL SIDTIM ( UTJDH, UTJDL, 1,   GAST )
          CALL SPIN ( -GAST * 15.D0, V1,   V2 )

*         SPECIAL OPTION SKIPS REMAINING TRANSFORMATIONS
          IF ( TJDH .LT. 0.D0 ) THEN
              VEC2(1) = V2(1)
              VEC2(2) = V2(2)
              VEC2(3) = V2(3)
          ELSE

*         APPLY NUTATION AND PRECESSION
          CALL NUTATE ( -TDBJD, V2,   V3 )
          CALL PRECES ( TDBJD, V3, T0,   V4 )

*         APPLY FRAME-TIE MATRIX
          CALL FRAME ( V4, -1, VEC2 )

          END IF

      ELSE
*         'CIO-TIO-THETA' MODE
*         SEE G. KAPLAN (2003), 'ANOTHER LOOK AT NON-ROTATING ORIGINS',
*         PROCEEDINGS OF IAU XXV JOINT DISCUSSION 16 (PREPRINT),
*         EQ. (3) AND (4).

*         APPLY POLAR MOTION, TRANSFORMING THE VECTOR TO THE TERRESTRIAL
*         INTERMEDIATE SYSTEM
          IF ( XP .EQ. 0.D0 .AND. YP .EQ. 0.D0 ) THEN
              V1(1) = VEC1(1)
              V1(2) = VEC1(2)
              V1(3) = VEC1(3)
          ELSE
              CALL WOBBLE ( TDBJD, XP, YP, VEC1,   V1 )
          END IF

*         OBTAIN THE BASIS VECTORS, IN THE GCRS, OF THE CELESTIAL
*         INTERMEDIATE SYSTEM
          CALL CIOLOC ( TDBJD,   RCIO, KCIO )
          CALL CIOBAS ( TDBJD, RCIO, KCIO,   X, Y, Z )

*         COMPUTE AND APPLY THE EARTH ROTATION ANGLE THETA, TRANSFORMING
*         THE VECTOR TO THE CELESTIAL INTERMEDIATE SYSTEM
          CALL EROT ( UTJDH, UTJDL,   THETA )
          CALL SPIN ( -THETA, V1,   V2 )

*         TRANSFORM THE VECTOR FROM THE CELESTIAL INTERMEDIATE SYSTEM
*         TO THE GCRS
          VEC2(1) = X(1) * V2(1) + Y(1) * V2(2) + Z(1) * V2(3)
          VEC2(2) = X(2) * V2(1) + Y(2) * V2(2) + Z(2) * V2(3)
          VEC2(3) = X(3) * V2(1) + Y(3) * V2(2) + Z(3) * V2(3)

      END IF

*     STORE COMPUTED POSITION VECTOR FOR POSSIBLE LATER RETRIEVAL
  50  CALL SETVEC ( VEC2 )

      RETURN

      END


      SUBROUTINE CELTER ( TJDH, TJDL, XP, YP, VEC1,   VEC2 )
*
*     THIS SUBROUTINE ROTATES A VECTOR FROM THE CELESTIAL TO THE
*     TERRESTRIAL SYSTEM.  SPECIFICALLY, IT TRANSFORMS A VECTOR IN THE
*     GCRS (A LOCAL SPACE-FIXED SYSTEM) TO THE ITRS (A ROTATING
*     EARTH-FIXED SYSTEM) BY APPLYING ROTATIONS FOR THE GCRS-TO-
*     DYNAMICAL FRAME TIE, PRECESSION, NUTATION, EARTH ROTATION,
*     AND POLAR MOTION.
*
*          TJDH   = UT1 JULIAN DATE, HIGH-ORDER PART (IN)
*          TJDL   = UT1 JULIAN DATE, LOW-ORDER PART (IN)
*                   THE JULIAN DATE MAY BE SPLIT AT ANY POINT, BUT
*                   FOR HIGHEST PRECISION, SET TJDH TO BE THE INTEGRAL
*                   PART OF THE JULIAN DATE, AND SET TJDL TO BE THE
*                   FRACTIONAL PART
*          XP     = CONVENTIONALLY-DEFINED X COORDINATE OF CELESTIAL
*                   INTERMEDIATE POLE WITH RESPECT TO ITRS POLE,
*                   IN ARCSECONDS (IN)
*          YP     = CONVENTIONALLY-DEFINED Y COORDINATE OF CELESTIAL
*                   INTERMEDIATE POLE WITH RESPECT TO ITRS POLE,
*                   IN ARCSECONDS (IN)
*          VEC1   = POSITION VECTOR, GEOCENTRIC EQUATORIAL RECTANGULAR
*                   COORDINATES, REFERRED TO GCRS AXES (CELESTIAL
*                   SYSTEM) (IN)
*          VEC2   = POSITION VECTOR, GEOCENTRIC EQUATORIAL RECTANGULAR
*                   COORDINATES, REFERRED TO ITRS AXES (TERRESTRIAL
*                   SYSTEM) (OUT)
*
*     NOTE 1:  SET XP=YP=0.D0 TO ELIMINATE POLAR MOTION ROTATION.
*
*     NOTE 2:  SEE ALSO SUBROUTINE SETDT TO SET THE VALUE OF DELTA-T
*     (DELTA-T = TT - UT1) TO BE USED HERE.
*
*     NOTE 3:  BOTH TJDH AND TJDL SHOULD BE NON-NEGATIVE FOR NORMAL USE
*     (TJDL=0.D0 IS OK).  A NEGATIVE VALUE OF TJDH IS USED TO INVOKE A
*     SPECIAL OPTION WHERE THE INPUT VECTOR IS ASSUMED TO BE WITH
*     RESPECT TO THE EQUATOR AND EQUINOX OF DATE, AND THE DATE FOR WHICH
*     THE TRANSFORMATION APPLIES IS TAKEN FROM TJDL ONLY.  THIS OPTION
*     WORKS ONLY IN 'EQUINOX' MODE.
*
*
      DOUBLE PRECISION TJDH,TJDL,XP,YP,VEC1,VEC2,
     .     T0,DELTAT,UTJDH,UTJDL,UTJD,TTJD,TDBJD,T,SECDIF,
     .     GAST,RCIO,THETA,V1,V2,V3,V4,X,Y,Z
      DIMENSION VEC1(3), VEC2(3), V1(3), V2(3), V3(3), V4(3),
     .     X(3), Y(3), Z(3)

*     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
      DATA T0 / 2451545.00000000D0 /

      CALL GETDT ( DELTAT )

      IF ( TJDH .GE. 0.D0 ) THEN
          UTJDH = TJDH
          UTJDL = TJDL
      ELSE
          UTJDH = TJDL
          UTJDL = 0.D0
      END IF
      UTJD = UTJDH + UTJDL

*     TIME ARGUMENT FOR PRECESSION AND NUTATION IS TDB
      TTJD = UTJD + DELTAT
      TDBJD = TTJD
      CALL TIMES ( TDBJD, T,   SECDIF )
      TDBJD = TTJD + SECDIF / 86400.D0

*     GET METHOD/ACCURACY MODE
      CALL GETMOD ( MODE )

      IF ( MODE .GE. 2 ) THEN
*         'EQUINOX' MODE

*         SPECIAL OPTION SKIPS INITIAL TRANSFORMATIONS
          IF ( TJDH .LT. 0.D0 ) THEN
              V3(1) = VEC1(1)
              V3(2) = VEC1(2)
              V3(3) = VEC1(3)
          ELSE

*         APPLY FRAME-TIE MATRIX
          CALL FRAME ( VEC1, 1, V1 )

*         APPLY PRECESSION AND NUTATION
          CALL PRECES ( T0, V1, TDBJD,   V2 )
          CALL NUTATE ( TDBJD, V2,   V3 )

          END IF

*         APPLY EARTH ROTATION
          CALL SIDTIM ( UTJDH, UTJDL, 1,   GAST )
          CALL SPIN ( GAST * 15.D0, V3,   V4 )

*         APPLY POLAR MOTION
          IF ( XP .EQ. 0.D0 .AND. YP .EQ. 0.D0 ) THEN
              VEC2(1) = V4(1)
              VEC2(2) = V4(2)
              VEC2(3) = V4(3)
          ELSE
              CALL WOBBLE ( -TDBJD, XP, YP, V4,   VEC2 )
          END IF

      ELSE
*         'CIO-TIO-THETA' MODE
*         SEE G. KAPLAN (2003), 'ANOTHER LOOK AT NON-ROTATING ORIGINS',
*         PROCEEDINGS OF IAU XXV JOINT DISCUSSION 16 (PREPRINT),
*         EQ. (3) AND (4).

*         OBTAIN THE BASIS VECTORS, IN THE GCRS, OF THE CELESTIAL
*         INTERMEDIATE SYSTEM
          CALL CIOLOC ( TDBJD,   RCIO, KCIO )
          CALL CIOBAS ( TDBJD, RCIO, KCIO,   X, Y, Z )

*         TRANSFORM THE VECTOR FROM THE GCRS TO THE
*         CELESTIAL INTERMEDIATE SYSTEM
          V1(1) = X(1) * VEC1(1) + X(2) * VEC1(2) + X(3) * VEC1(3)
          V1(2) = Y(1) * VEC1(1) + Y(2) * VEC1(2) + Y(3) * VEC1(3)
          V1(3) = Z(1) * VEC1(1) + Z(2) * VEC1(2) + Z(3) * VEC1(3)

*         COMPUTE AND APPLY THE EARTH ROTATION ANGLE THETA, TRANSFORMING
*         THE VECTOR TO THE TERRESTRIAL INTERMEDIATE SYSTEM
          CALL EROT ( UTJDH, UTJDL,   THETA )
          CALL SPIN ( THETA, V1,   V2 )

*         APPLY POLAR MOTION, TRANSFORMING THE VECTOR TO THE ITRS
          IF ( XP .EQ. 0.D0 .AND. YP .EQ. 0.D0 ) THEN
              VEC2(1) = V2(1)
              VEC2(2) = V2(2)
              VEC2(3) = V2(3)
          ELSE
              CALL WOBBLE ( -TDBJD, XP, YP, V2,   VEC2 )
          END IF

      END IF

*     STORE COMPUTED POSITION VECTOR FOR POSSIBLE LATER RETRIEVAL
  50  CALL SETVEC ( VEC2 )

      RETURN

      END




      SUBROUTINE GETHIP ( RAH, DECH, PMRAH, PMDECH, PARXH, RVH,
     .                    RA2, DEC2, PMRA2, PMDEC2, PARX2, RV2 )
*
*     THIS SUBROUTINE CONVERTS HIPPARCOS DATA AT EPOCH J1991.25
*     TO EPOCH J2000.0.  TO BE USED ONLY FOR HIPPARCOS OR TYCHO STARS
*     WITH LINEAR SPACE MOTION.  BOTH INPUT AND OUTPUT DATA IS IN THE
*     ICRS.
*
*          RAH    = HIPPARCOS RIGHT ASCENSION IN DEGREES (IN)
*          DECH   = HIPPARCOS DECLINATION IN DEGREES (IN)
*          PMRAH  = HIPPARCOS PROPER MOTION IN RA
*                   IN MILLIARCSECONDS/YEAR (IN)
*          PMDECH = HIPPARCOS PROPER MOTION IN DEC
*                   IN MILLIARCSECONDS/YEAR (IN)
*          PARXH  = HIPPARCOS PARALLAX IN MILLIARCSECONDS (IN)
*          RVH    = RADIAL VELOCITY AT HIPPARCOS EPOCH
*                   IN KILOMETERS/SECOND (IN)
*          RA2    = RIGHT ASCENSION AT J2000.0 IN HOURS (OUT)
*          DEC2   = DECLINATION AT J2000.0 IN DEGREES (OUT)
*          PMRA2  = PROPER MOTION IN RA AT J2000.0
*                   IN MILLIARCSECONDS/YEAR (OUT)
*          PMDEC2 = PROPER MOTION IN DEC AT J2000.0
*                   IN MILLIARCSECONDS/YEAR (OUT)
*          PARX2  = PARALLAX AT J2000.0 IN MILLIARCSECONDS (OUT)
*          RV2    = RADIAL VELOCITY AT J2000.0 IN KILOMETERS/SECOND
*                   (OUT)
*
*     NOTE:  INPUT RA IS IN DEGREES, AS PER HIPPARCOS, BUT OUTPUT RA
*     IS IN HOURS.
*
*
      DOUBLE PRECISION RAH,DECH,PMRAH,PMDECH,PARXH,RVH,
     .     RA2,DEC2,PMRA2,PMDEC2,PARX2,RV2,EPOCH1,EPOCH2

      DATA EPOCH1, EPOCH2 / 2448349.0625D0, 2451545.0000D0 /

      CALL CATRAN ( 1,EPOCH1,RAH/15.D0,DECH,PMRAH,PMDECH,PARXH,RVH,
     .                EPOCH2,RA2,      DEC2,PMRA2,PMDEC2,PARX2,RV2 )

      RETURN

      END



      SUBROUTINE CATRAN ( IT,
     .                    DATE1, RA1, DEC1, PMRA1, PMDEC1, PARX1, RV1,
     .                    DATE2, RA2, DEC2, PMRA2, PMDEC2, PARX2, RV2 )
*
*     THIS SUBROUTINE TRANSFORMS A STAR'S CATALOG QUANTITIES FOR
*     A CHANGE OF EPOCH AND/OR EQUATOR AND EQUINOX.  IT CAN ALSO BE
*     USED TO ROTATE CATALOG QUANTITIES ON THE DYNAMICAL EQUATOR AND
*     EQUINOX OF J2000.0 TO THE ICRS OR VICE VERSA.
*
*          IT     = TRANSFORMATION OPTION (IN)
*                   SET IT=1 TO CHANGE EPOCH (SAME EQUATOR AND EQUINOX)
*                   SET IT=2 TO CHANGE EQUATOR AND EQUINOX (SAME EPOCH)
*                   SET IT=3 TO CHANGE EQUATOR AND EQUINOX AND EPOCH
*                   SET IT=4 TO CHANGE EQUATOR AND EQUINOX OF J2000.0
*                            TO ICRS
*                   SET IT=5 TO CHANGE ICRS TO EQUATOR AND EQUINOX OF
*                            J2000.0
*          DATE1  = TT JULIAN DATE, OR YEAR, OF ORIGINAL CATALOG
*                   DATA (THE FOLLOWING SIX ARGUMENTS) (IN)
*          RA1    = ORIGINAL MEAN RIGHT ASCENSION IN HOURS (IN)
*          DEC1   = ORIGINAL MEAN DECLINATION IN DEGREES (IN)
*          PMRA1  = ORIGINAL PROPER MOTION IN RA
*                   IN MILLIARCSECONDS/YEAR (IN)
*          PMDEC1 = ORIGINAL PROPER MOTION IN DEC
*                   IN MILLIARCSECONDS/YEAR (IN)
*          PARX1  = ORIGINAL PARALLAX IN MILLIARCSECONDS (IN)
*          RV1    = ORIGINAL RADIAL VELOCITY IN KILOMETERS/SECOND
*                   (IN)
*          DATE2  = TT JULIAN DATE, OR YEAR, FOR TRANSFORMED
*                   OUTPUT DATA (THE FOLLOWING SIX ARGUMENTS) (IN)
*          RA2    = TRANSFORMED MEAN RIGHT ASCENSION IN HOURS (OUT)
*          DEC2   = TRANSFORMED MEAN DECLINATION IN DEGREES (OUT)
*          PMRA2  = TRANSFORMED PROPER MOTION IN RA
*                   IN MILLIARCSECONDS/YEAR (OUT)
*          PMDEC2 = TRANSFORMED PROPER MOTION IN DEC
*                   IN MILLIARCSECONDS/YEAR (OUT)
*          PARX2  = TRANSFORMED PARALLAX IN MILLIARCSECONDS (OUT)
*          RV2    = TRANSFORMED RADIAL VELOCITY IN KILOMETERS/SECOND
*                   (OUT)
*
*     NOTE 1:  DATE1 AND DATE2 MAY BE SPECIFIED EITHER AS A JULIAN
*     DATE (E.G., 2433282.5D0) OR A JULIAN YEAR AND FRACTION
*     (E.G., 1950.0D0).  VALUES LESS THAN 10000 ARE ASSUMED TO
*     BE YEARS.  FOR IT=2 OR IT=3, EITHER DATE1 OR DATE2 MUST BE
*     2451545.0 OR 2000.0 (J2000.0).  FOR IT=4 AND IT=5, DATE1 AND
*     DATE2 ARE IGNORED.
*
*     NOTE 2:  IT=1 UPDATES THE STAR'S DATA TO ACCOUNT FOR
*     THE STAR'S SPACE MOTION BETWEEN THE FIRST AND SECOND DATES,
*     WITHIN A FIXED REFERENCE SYSTEM.  IT=2 APPLIES A ROTATION
*     OF THE REFERENCE SYSTEM CORRESPONDING TO PRECESSION BETWEEN
*     THE FIRST AND SECOND DATES, BUT LEAVES THE STAR FIXED IN SPACE.
*     IT=3 PROVIDES BOTH TRANSFORMATIONS.  IT=4 AND IT=5 PROVIDE A
*     A FIXED ROTATION ABOUT VERY SMALL ANGLES (<0.1 ARCSECOND) TO
*     TAKE DATA FROM THE DYNAMICAL SYSTEM OF J2000.0 TO THE ICRS (IT=4)
*     OR VICE VERSA (IT=5).
*
*     NOTE 3:  FOR IT=1, INPUT DATA CAN BE IN ANY FIXED REFERENCE
*     SYSTEM. FOR IT=2 OR IT=3, THIS SUBROUTINE ASSUMES THE INPUT DATA
*     IS IN THE DYNAMICAL SYSTEM AND PRODUCES OUTPUT IN THE DYNAMICAL
*     SYSTEM.  FOR IT=4, THE INPUT DATA MUST BE ON THE DYNAMICAL EQUATOR
*     AND EQUINOX OF J2000.0.  FOR IT=5, THE INPUT DATA MUST BE IN THE
*     ICRS.
*
*     NOTE 4:  THIS SUBROUTINE CANNOT BE PROPERLY USED TO BRING DATA
*     FROM OLD STAR CATALOGS INTO THE MODERN SYSTEM, BECAUSE
*     OLD CATALOGS WERE COMPILED USING A SET OF CONSTANTS THAT ARE
*     INCOMPATIBLE WITH MODERN VALUES.  IN PARTICULAR, IT SHOULD NOT
*     BE USED FOR CATALOGS WHOSE POSITIONS AND PROPER MOTIONS WERE
*     DERIVED BY ASSUMING A PRECESSION CONSTANT SIGNIFICANTLY DIFFERENT
*     FROM THE VALUE IMPLICIT IN SUBROUTINE PRECES.
*
*
      DOUBLE PRECISION DATE1,RA1,DEC1,PMRA1,PMDEC1,PARX1,RV1,
     .     DATE2,RA2,DEC2,PMRA2,PMDEC2,PARX2,RV2,
     .     PI,SECCON,AUKM,C,TJD1,POS1,VEL1,TJD2,POS2,VEL2,
     .     PARALX,DIST,R,D,CRA,SRA,CDC,SDC,K,PMR,PMD,RVL,
     .     XYPROJ,DCOS,DSIN,DATAN2
      INTEGER IT,J

      DIMENSION POS1(3), VEL1(3), POS2(3), VEL2(3)

      SAVE

      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( SECCON = 180.D0 * 3600.D0 / PI )

      DATA NTIMES / 0 /

      NTIMES = NTIMES + 1
      IF ( NTIMES .EQ. 1 ) THEN
*         GET LENGTH OF AU IN KILOMETERS
          CALL ASTCON ( 'AU', 1.D-3, AUKM )
*         GET C, THE SPEED OF LIGHT IN KILOMETERS/SECOND
          CALL ASTCON ( 'C', 1.D-3,   C )
      END IF

* --- IF NECESSARY, COMPUTE JULIAN DATES ------------------------------

*     SUBROUTINE USES TDB JULIAN DATES INTERNALLY, BUT NO
*     DISTINCTION BETWEEN TDB AND TT IS NECESSARY

      IF ( DATE1 .LT. 10000.D0 ) THEN
           TJD1 = 2451545.0D0 + ( DATE1 - 2000.D0 ) * 365.25D0
      ELSE
           TJD1 = DATE1
      END IF
      IF ( DATE2 .LT. 10000.D0 ) THEN
           TJD2 = 2451545.0D0 + ( DATE2 - 2000.D0 ) * 365.25D0
      ELSE
           TJD2 = DATE2
      END IF

* --- CONVERT INPUT ANGULAR COMPONENTS TO VECTORS ---------------------

*     IF PARALLAX IS UNKNOWN, UNDETERMINED, OR ZERO, SET IT TO 1E-6
*     MILLIARCSECOND, CORRESPONDING TO A DISTANCE OF 1 GIGAPARSEC
      PARALX = PARX1
      IF ( PARALX .LE. 0.D0 ) PARALX = 1.D-6

*     CONVERT RIGHT ASCENSION, DECLINATION, AND PARALLAX TO POSITION
*     VECTOR IN EQUATORIAL SYSTEM WITH UNITS OF AU
      DIST = 1.D0 / DSIN ( PARALX * 1.D-3 / SECCON )
      R = RA1 * 54000.D0 / SECCON
      D = DEC1 * 3600.D0 / SECCON
      CRA = DCOS ( R )
      SRA = DSIN ( R )
      CDC = DCOS ( D )
      SDC = DSIN ( D )
      POS1(1) = DIST * CDC * CRA
      POS1(2) = DIST * CDC * SRA
      POS1(3) = DIST * SDC

*     COMPUTE DOPPLER FACTOR, WHICH ACCOUNTS FOR CHANGE IN
*     LIGHT TRAVEL TIME TO STAR
      K = 1.D0 / ( 1.D0 - RV1 / C )

*     CONVERT PROPER MOTION AND RADIAL VELOCITY TO ORTHOGONAL
*     COMPONENTS OF MOTION, IN SPHERICAL POLAR SYSTEM AT STAR'S
*     ORIGINAL POSITION, WITH UNITS OF AU/DAY
      PMR = PMRA1  / ( PARALX * 365.25D0 ) * K
      PMD = PMDEC1 / ( PARALX * 365.25D0 ) * K
      RVL = RV1 * 86400.D0 / AUKM          * K

*     TRANSFORM MOTION VECTOR TO EQUATORIAL SYSTEM
      VEL1(1) = - PMR * SRA - PMD * SDC * CRA + RVL * CDC * CRA
      VEL1(2) =   PMR * CRA - PMD * SDC * SRA + RVL * CDC * SRA
      VEL1(3) =               PMD * CDC       + RVL * SDC

* --- UPDATE STAR'S POSITION VECTOR FOR SPACE MOTION ------------------
*     (ONLY IF IT=1 OR IT=3)

      IF ( IT .EQ. 1 .OR. IT .EQ. 3 ) THEN
          DO 22 J=1,3
              POS2(J) = POS1(J) + VEL1(J) * ( TJD2 - TJD1 )
              VEL2(J) = VEL1(J)
   22     CONTINUE
      ELSE
          DO 24 J=1,3
              POS2(J) = POS1(J)
              VEL2(J) = VEL1(J)
   24     CONTINUE
      END IF

* --- PRECESS POSITION AND VELOCITY VECTORS ---------------------------
*     (ONLY IF IT=2 OR IT=3)

      IF ( IT .EQ. 2 .OR. IT .EQ. 3 ) THEN
          DO 32 J=1,3
              POS1(J) = POS2(J)
              VEL1(J) = VEL2(J)
   32     CONTINUE
          CALL PRECES ( TJD1, POS1, TJD2,    POS2 )
          CALL PRECES ( TJD1, VEL1, TJD2,    VEL2 )
      END IF

* --- ROTATE DYNAMICAL J2000.0 POSITION AND VELOCITY VECTORS TO ICRS --
*     (ONLY IF IT=4)

      IF ( IT .EQ. 4 ) THEN
          CALL FRAME ( POS1, -1,    POS2 )
          CALL FRAME ( VEL1, -1,    VEL2 )
      END IF

* --- ROTATE ICRS POSITION AND VELOCITY VECTORS TO DYNAMICAL J2000.0 --
*     (ONLY IF IT=5)

      IF ( IT .EQ. 5 ) THEN
          CALL FRAME ( POS1, 1,    POS2 )
          CALL FRAME ( VEL1, 1,    VEL2 )
      END IF

* --- CONVERT VECTORS BACK TO ANGULAR COMPONENTS FOR OUTPUT -----------

*     FROM UPDATED POSITION VECTOR, OBTAIN STAR'S NEW POSITION
*     EXPRESSED AS ANGULAR QUANTITIES
      XYPROJ = DSQRT ( POS2(1)**2 + POS2(2)**2 )
      R = 0.D0
      IF ( XYPROJ .GT. 0.D0 ) R = DATAN2 ( POS2(2), POS2(1) )
      RA2 = R * SECCON / 54000.D0
      IF ( RA2 .LT.  0.D0 ) RA2 = RA2 + 24.D0
      IF ( RA2 .GE. 24.D0 ) RA2 = RA2 - 24.D0
      D = DATAN2 ( POS2(3), XYPROJ  )
      DEC2 = D * SECCON / 3600.D0
      DIST = DSQRT ( POS2(1)**2 + POS2(2)**2 + POS2(3)**2 )
      PARALX = DASIN ( 1.D0 / DIST ) * SECCON * 1.D3
      PARX2 = PARALX

*     TRANSFORM MOTION VECTOR BACK TO SPHERICAL POLAR SYSTEM AT STAR'S
*     NEW POSITION
      CRA = DCOS ( R )
      SRA = DSIN ( R )
      CDC = DCOS ( D )
      SDC = DSIN ( D )
      PMR = - VEL2(1) * SRA       + VEL2(2) * CRA
      PMD = - VEL2(1) * CRA * SDC - VEL2(2) * SRA * SDC + VEL2(3) * CDC
      RVL =   VEL2(1) * CRA * CDC + VEL2(2) * SRA * CDC + VEL2(3) * SDC

*     CONVERT COMPONENTS OF MOTION FROM AU/DAY TO NORMAL
*     CATALOG UNITS
      PMRA2  = PMR * PARALX * 365.25D0   / K
      PMDEC2 = PMD * PARALX * 365.25D0   / K
      RV2    = RVL * ( AUKM / 86400.D0 ) / K

*     TAKE CARE OF ZERO-PARALLAX CASE
      IF ( PARX2 .LE. 1.01D-6 ) THEN
          PARX2 = 0.D0
          RV2 = RV1
      END IF

*     STORE COMPUTED POSITION VECTOR FOR POSSIBLE LATER RETRIEVAL
      CALL SETVEC ( POS2 )

      RETURN

      END



      SUBROUTINE ZDAZ ( UJD, XP, YP, GLON, GLAT, HT, RA, DEC, IREFR,
     .                  ZD, AZ, RAR, DECR )
*
*     THIS SUBROUTINE TRANSFORMS TOPOCENTRIC RIGHT ASCENSION AND
*     DECLINATION TO ZENITH DISTANCE AND AZIMUTH.  THIS ROUTINE USES
*     A METHOD THAT PROPERLY ACCOUNTS FOR POLAR MOTION, WHICH IS
*     SIGNIFICANT AT THE SUB-ARCSECOND LEVEL.  THIS SUBROUTINE
*     CAN ALSO ADJUST COORDINATES FOR ATMOSPHERIC REFRACTION.
*
*          UJD    = UT1 JULIAN DATE (IN)
*          XP     = CONVENTIONALLY-DEFINED X COORDINATE OF CELESTIAL
*                   INTERMEDIATE POLE WITH RESPECT TO ITRS POLE, IN
*                   ARCSECONDS (IN)
*          YP     = CONVENTIONALLY-DEFINED Y COORDINATE OF CELESTIAL
*                   INTERMEDIATE POLE WITH RESPECT TO ITRS POLE, IN
*                   ARCSECONDS (IN)
*          GLON   = GEODETIC (ITRS) LONGITUDE (EAST +) OF OBSERVER
*                   IN DEGREES (IN)
*          GLAT   = GEODETIC (ITRS) LATITUDE (NORTH +) OF OBSERVER
*                   IN DEGREES (IN)
*          HT     = HEIGHT OF OBSERVER IN METERS (IN)
*          RA     = TOPOCENTRIC RIGHT ASCENSION OF OBJECT OF INTEREST,
*                   IN HOURS, REFERRED TO TRUE EQUATOR AND EQUINOX
*                   OF DATE (IN)
*          DEC    = TOPOCENTRIC DECLINATION OF OBJECT OF INTEREST,
*                   IN DEGREES, REFERRED TO TRUE EQUATOR OF DATE (IN)
*          IREFR  = ATMOSPHERIC REFRACTION OPTION (IN)
*                   SET IREFR=0 FOR NO REFRACTION
*                   SET IREFR=1 TO INCLUDE REFRACTION
*          ZD     = TOPOCENTRIC ZENITH DISTANCE IN DEGREES,
*                   AFFECTED BY REFRACTION IF IREFR=1 (OUT)
*          AZ     = TOPOCENTRIC AZIMUTH (MEASURED EAST FROM NORTH)
*                   IN DEGREES (OUT)
*          RAR    = TOPOCENTRIC RIGHT ASCENSION OF OBJECT OF INTEREST,
*                   IN HOURS, REFERRED TO TRUE EQUATOR AND EQUINOX
*                   OF DATE, AFFECTED BY REFRACTION IF IREFR=1 (OUT)
*          DECR   = TOPOCENTRIC DECLINATION OF OBJECT OF INTEREST,
*                   IN DEGREES, REFERRED TO TRUE EQUATOR OF DATE,
*                   AFFECTED BY REFRACTION IF IREFR=1 (OUT)
*
*     NOTE 1:  XP AND YP CAN BE SET TO ZERO IF SUB-ARCSECOND ACCURACY IS
*     NOT NEEDED.  HT IS USED ONLY FOR REFRACTION, IF IREFR=1.  RA AND
*     DEC CAN BE OBTAINED FROM TPSTAR, TPPLAN, OR PLACE.
*
*     NOTE 2:  THE DIRECTONS ZD=0 (ZENITH) AND AZ=0 (NORTH) ARE
*     HERE CONSIDERED FIXED IN THE TERRESTRIAL SYSTEM.  SPECIFICALLY,
*     THE ZENITH IS ALONG THE GEODETIC NORMAL, AND NORTH IS TOWARD
*     THE ITRS POLE.
*
*     NOTE 3:  IF IREFR=0, THEN RAR=RA AND DECR=DEC.
*
*     NOTE 4: INPUT PARAMETERS XP, YP WERE X, Y IN NOVAS F3.0.
*     THE NAMES WERE CHANGED FOR CONSISTANCY THROUGHOUT NOVAS AND WITH
*     IERS CONVENTIONS.
*
      DOUBLE PRECISION UJD,XP,YP,GLON,GLAT,HT,RA,DEC,ZD,AZ,RAR,DECR,
     .     PI,DEGRAD,RADDEG,
     .     SINLAT,COSLAT,SINLON,COSLON,SINDC,COSDC,SINRA,COSRA,
     .     UZE,UNE,UWE,UZ,UN,UW,P,PR,PZ,PN,PW,PROJ,
     .     ZD0,ZD1,REFR,SINZD,COSZD,SINZD0,COSZD0,
     .     DSIN,DCOS,DSQRT,DATAN2
      DIMENSION UZE(3), UNE(3), UWE(3), UZ(3), UN(3), UW(3),
     .     P(3), PR(3)

      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( DEGRAD = PI / 180.D0           )
      PARAMETER ( RADDEG = 180.D0 / PI           )

      RAR    = RA
      DECR   = DEC
      SINLAT = DSIN ( GLAT * DEGRAD )
      COSLAT = DCOS ( GLAT * DEGRAD )
      SINLON = DSIN ( GLON * DEGRAD )
      COSLON = DCOS ( GLON * DEGRAD )
      SINDC  = DSIN ( DEC * DEGRAD )
      COSDC  = DCOS ( DEC * DEGRAD )
      SINRA  = DSIN ( RA * 15.D0 * DEGRAD )
      COSRA  = DCOS ( RA * 15.D0 * DEGRAD )

* --- SET UP ORTHONORMAL BASIS VECTORS IN LOCAL EARTH-FIXED SYSTEM ----

*     DEFINE VECTOR TOWARD LOCAL ZENITH IN EARTH-FIXED SYSTEM (Z AXIS)
      UZE(1) =  COSLAT * COSLON
      UZE(2) =  COSLAT * SINLON
      UZE(3) =  SINLAT

*     DEFINE VECTOR TOWARD LOCAL NORTH IN EARTH-FIXED SYSTEM (X AXIS)
      UNE(1) = -SINLAT * COSLON
      UNE(2) = -SINLAT * SINLON
      UNE(3) =  COSLAT

*     DEFINE VECTOR TOWARD LOCAL WEST IN EARTH-FIXED SYSTEM (Y AXIS)
      UWE(1) =  SINLON
      UWE(2) = -COSLON
      UWE(3) =  0.D0

* --- OBTAIN VECTORS IN CELESTIAL SYSTEM ------------------------------

*     ROTATE EARTH-FIXED ORTHONORMAL BASIS VECTORS TO CELESTIAL SYSTEM
*     (WRT EQUATOR AND EQUINOX OF DATE)
      CALL EQINOX
      CALL TERCEL ( -1.D0, UJD, XP, YP, UZE,   UZ )
      CALL TERCEL ( -1.D0, UJD, XP, YP, UNE,   UN )
      CALL TERCEL ( -1.D0, UJD, XP, YP, UWE,   UW )
      CALL RESUME

*     DEFINE UNIT VECTOR P TOWARD OBJECT IN CELESTIAL SYSTEM
*     (WRT EQUATOR AND EQUINOX OF DATE)
      P(1) = COSDC * COSRA
      P(2) = COSDC * SINRA
      P(3) = SINDC

* --- COMPUTE COORDINATES OF OBJECT WRT ORTHONORMAL BASIS -------------

*     COMPUTE COMPONENTS OF P -- PROJECTIONS OF P ONTO ROTATED
*     EARTH-FIXED BASIS VECTORS
      PZ = P(1) * UZ(1) + P(2) * UZ(2) + P(3) * UZ(3)
      PN = P(1) * UN(1) + P(2) * UN(2) + P(3) * UN(3)
      PW = P(1) * UW(1) + P(2) * UW(2) + P(3) * UW(3)

*     COMPUTE AZIMUTH AND ZENITH DISTANCE
      PROJ = DSQRT ( PN**2 + PW**2 )
      AZ = 0.D0
      IF ( PROJ .GT. 0.D0 ) AZ = -DATAN2 ( PW, PN ) * RADDEG
      IF ( AZ .LT.   0.D0 ) AZ = AZ + 360.D0
      IF ( AZ .GE. 360.D0 ) AZ = AZ - 360.D0
      ZD = DATAN2 ( PROJ, PZ ) * RADDEG

* --- APPLY ATMOSPHERIC REFRACTION IF REQUESTED -----------------------

      IF ( IREFR .EQ. 1 ) THEN

*         GET REFRACTION IN ZENITH DISTANCE
*         ITERATIVE PROCESS REQUIRED BECAUSE REFRACTION ALGORITHMS ARE
*         ALWAYS A FUNCTION OF OBSERVED (NOT COMPUTED) ZENITH DISTANCE
          ZD0 = ZD
  40      ZD1 = ZD
          CALL REFRAC ( HT, ZD,   REFR )
          ZD = ZD0 - REFR
*         REQUIRE CONVERGENCE TO 0.1 ARCSEC (ACTUAL ACCURACY LESS)
          IF ( DABS ( ZD - ZD1 ) .GT. 3.D-5 ) GO TO 40

*         APPLY REFRACTION TO CELESTIAL COORDINATES OF OBJECT
          IF ( REFR .GT. 0.D0 .AND. ZD .GT. 3.D-4 ) THEN

*             SHIFT POSITION VECTOR OF OBJECT IN CELESTIAL SYSTEM
*             TO ACCOUNT FOR FOR REFRACTION (SEE USNO/AA TECHNICAL
*             NOTE 1998-09)
              SINZD  = DSIN ( ZD * DEGRAD )
              COSZD  = DCOS ( ZD * DEGRAD ) 
              SINZD0 = DSIN ( ZD0 * DEGRAD )
              COSZD0 = DCOS ( ZD0 * DEGRAD )
*             COMPUTE REFRACTED POSITION VECTOR
              DO 50 J = 1, 3
  50          PR(J) = ( ( P(J) - COSZD0 * UZ(J) ) / SINZD0 ) * SINZD
     .                 +                  UZ(J)              * COSZD   

*             COMPUTE REFRACTED RIGHT ASCENSION AND DECLINATION
              PROJ = DSQRT ( PR(1)**2 + PR(2)**2 )
              RAR = 0.D0
              IF ( PROJ .GT. 0.D0 ) RAR = DATAN2 ( PR(2), PR(1) )
     .                                    * RADDEG / 15.D0
              IF ( RAR .LT.  0.D0 ) RAR = RAR + 24.D0
              IF ( RAR .GE. 24.D0 ) RAR = RAR - 24.D0
              DECR = DATAN2 ( PR(3), PROJ ) * RADDEG

          END IF

      END IF

* ---------------------------------------------------------------------

      RETURN

      END



      SUBROUTINE GCRSEQ ( TJD, ICOORD, RAG, DECG,   RA, DEC )
*
*     THIS SUBROUTINE CONVERTS GCRS RIGHT ASCENSION AND DECLINATION
*     TO COORDINATES WITH RESPECT TO THE EQUATOR OF DATE (MEAN OR TRUE).
*     FOR COORDINATES WITH RESPECT TO THE TRUE EQUATOR OF DATE, THE
*     ORIGIN OF RIGHT ASCENSION CAN BE EITHER THE TRUE EQUINOX OR THE
*     CELESTIAL INTERMEDIATE ORIGIN (CIO).
*
*          TJD    = TT JULIAN DATE OF EQUATOR TO BE USED FOR
*                   OUTPUT COORDINATES (IN)
*          ICOORD = COORDINATE SYSTEM SELECTION FOR OUTPUT
*                   COORDINATES (IN)
*                   SET ICOORD=0 FOR MEAN EQUATOR AND EQUINOX OF DATE
*                   SET ICOORD=1 FOR TRUE EQUATOR AND EQUINOX OF DATE
*                   SET ICOORD=2 FOR TRUE EQUATOR AND CIO OF DATE
*          RAG    = GCRS RIGHT ASCENSION IN HOURS (IN)
*          DECG   = GCRS DECLINATION IN DEGREES (IN)
*          RA     = RIGHT ASCENSION IN HOURS, REFERRED TO SPECIFIED
*                   EQUATOR AND RIGHT ASCENSION ORIGIN OF DATE (OUT)
*          DEC    = DECLINATION IN DEGREES, REFERRED TO SPECIFIED
*                   EQUATOR OF DATE (OUT)
*
*
      DOUBLE PRECISION TJD,RAG,DECG,RA,DEC,PI,RADCON,T0,T1,T,SECDIF,R,D,
     .     POS1,POS2,POS3,POS4,RCIO,X,Y,Z,DSIN,DCOS
      DIMENSION POS1(3), POS2(3), POS3(3), POS4(3), X(3), Y(3), Z(3)
      SAVE

      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( RADCON = PI / 180.D0           )

*     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
      DATA T0 / 2451545.00000000D0 /

*     T1 IS THE TDB JULIAN DATE
      CALL TIMES ( TJD, T,   SECDIF )
      T1 = TJD + SECDIF / 86400.D0

*     FORM POSITION VECTOR IN EQUATORIAL SYSTEM FROM INPUT COORDINATES
      R = RAG * 15.D0 * RADCON
      D = DECG * RADCON
      POS1(1) = DCOS ( D ) * DCOS ( R )
      POS1(2) = DCOS ( D ) * DSIN ( R )
      POS1(3) = DSIN ( D )

      IF ( ICOORD .LE. 1 ) THEN

*         TRANSFORM THE POSITION VECTOR FROM GCRS TO MEAN EQUATOR AND
*         EQUINOX OF DATE
          CALL FRAME  ( POS1, 1,   POS2 )
          CALL PRECES ( T0, POS2, T1,   POS3 )
*         IF REQUESTED, TRANSFORM FURTHER TO TRUE EQUATOR AND EQUINOX
*         OF DATE
          IF ( ICOORD .EQ. 1 ) THEN
              CALL NUTATE ( T1, POS3,   POS4)
          ELSE
              POS4(1) = POS3(1)
              POS4(2) = POS3(2)
              POS4(3) = POS3(3)
          END IF

      ELSE

*         OBTAIN THE BASIS VECTORS, IN THE GCRS, OF THE CELESTIAL
*         INTERMEDIATE SYSTEM
          CALL CIOLOC ( T1,   RCIO, KCIO )
          CALL CIOBAS ( T1, RCIO, KCIO,   X, Y, Z )

*         TRANSFORM POSITION VECTOR TO THE CELESTIAL INTERMEDIATE SYSTEM
*         (WHICH HAS THE CIO AS ITS ORIGIN OF RIGHT ASCENSION)
          POS4(1) = X(1) * POS1(1) + X(2) * POS1(2) + X(3) * POS1(3)
          POS4(2) = Y(1) * POS1(1) + Y(2) * POS1(2) + Y(3) * POS1(3)
          POS4(3) = Z(1) * POS1(1) + Z(2) * POS1(2) + Z(3) * POS1(3)

      END IF

      CALL ANGLES ( POS4,   RA, DEC )

      RETURN

      END



      SUBROUTINE EQECL ( TJD, ICOORD, RA, DEC,   ELON, ELAT )
*
*     THIS SUBROUTINE CONVERTS RIGHT ASCENSION AND DECLINATION
*     TO ECLIPTIC LONGITUDE AND LATITUDE.
*
*          TJD    = TT JULIAN DATE OF EQUATOR, EQUINOX, AND ECLIPTIC
*                   USED FOR COORDINATES (IN)
*          ICOORD = COORDINATE SYSTEM SELECTION (IN)
*                   SET ICOORD=0 FOR MEAN EQUATOR AND EQUINOX OF DATE
*                   SET ICOORD=1 FOR TRUE EQUATOR AND EQUINOX OF DATE
*                   (ECLIPTIC IS ALWAYS THE MEAN PLANE)
*          RA     = RIGHT ASCENSION IN HOURS, REFERRED TO SPECIFIED
*                   EQUATOR AND EQUINOX OF DATE (IN)
*          DEC    = DECLINATION IN DEGREES, REFERRED TO SPECIFIED
*                   EQUATOR AND EQUINOX OF DATE (IN)
*          ELON   = ECLIPTIC LONGITUDE IN DEGREES, REFERRED TO SPECIFIED
*                   ECLIPTIC AND EQUINOX OF DATE (OUT)
*          ELAT   = ECLIPTIC LATITUDE IN DEGREES, REFERRED TO SPECIFIED
*                   ECLIPTIC AND EQUINOX OF DATE (OUT)
*
*     NOTE:  TO CONVERT ICRS RA AND DEC TO ECLIPTIC COORDINATES (MEAN
*     ECLIPTIC AND EQUINOX OF J2000.0), SET TJD = 0.D0 AND ICOORD = 0.
*     EXCEPT FOR THE INPUT TO THIS CASE, ALL COORDINATES ARE DYNAMICAL.
*
*
      DOUBLE PRECISION TJD,RA,DEC,ELON,ELAT,PI,RADCON,R,D,XYPROJ,E,
     .     POS1,POS2,DSIN,DCOS,DSQRT,DATAN2
      DIMENSION POS1(3), POS2(3)
      SAVE

      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( RADCON = PI / 180.D0           )

*     FORM POSITION VECTOR IN EQUATORIAL SYSTEM FROM INPUT COORDINATES
      R = RA * 15.D0 * RADCON
      D = DEC * RADCON
      POS1(1) = DCOS ( D ) * DCOS ( R )
      POS1(2) = DCOS ( D ) * DSIN ( R )
      POS1(3) = DSIN ( D )

*     CONVERT THE VECTOR FROM EQUATORIAL TO ECLIPTIC SYSTEM
      CALL EQEC ( TJD, ICOORD, POS1,   POS2 )

*     DECOMPOSE ECLIPTIC VECTOR INTO ECLIPTIC LONGITUDE AND LATITUDE
      XYPROJ = DSQRT ( POS2(1)**2 + POS2(2)**2 )
      E = 0.D0
      IF ( XYPROJ .GT. 0.D0 ) E = DATAN2 ( POS2(2), POS2(1) )
      ELON = E / RADCON
      IF ( ELON .LT. 0.D0 ) ELON = ELON + 360.D0
      E = DATAN2 ( POS2(3), XYPROJ )
      ELAT = E / RADCON

      RETURN

      END
      
      
      
      SUBROUTINE EQEC ( TJD, ICOORD, POS1,   POS2 )
*
*     THIS SUBROUTINE CONVERTS AN EQUATORIAL POSITION VECTOR TO
*     AN ECLIPTIC POSITION VECTOR.
*
*          TJD    = TT JULIAN DATE OF EQUATOR, EQUINOX, AND ECLIPTIC
*                   USED FOR COORDINATES (IN)
*          ICOORD = COORDINATE SYSTEM SELECTION (IN)
*                   SET ICOORD=0 FOR MEAN EQUATOR AND EQUINOX OF DATE
*                   SET ICOORD=1 FOR TRUE EQUATOR AND EQUINOX OF DATE
*                   (ECLIPTIC IS ALWAYS THE MEAN PLANE)
*          POS1   = POSITION VECTOR, REFERRED TO SPECIFIED
*                   EQUATOR AND EQUINOX OF DATE (IN)
*          POS2   = POSITION VECTOR, REFERRED TO SPECIFIED
*                   ECLIPTIC AND EQUINOX OF DATE (OUT)
*
*     NOTE:  TO CONVERT ICRS VECTORS TO ECLIPTIC VECTORS (MEAN ECLIPTIC       
*     AND EQUINOX OF J2000.0 ONLY), SET TJD = 0.D0 AND ICOORD = 0.
*     EXCEPT FOR THE INPUT TO THIS CASE, ALL VECTORS ARE ASSUMED TO
*     BE WITH RESPECT TO A DYNAMICAL SYSTEM.
*
*
      DOUBLE PRECISION TJD,POS1,POS2,POS0,PI,RADCON,T0,T1,T,SECDIF,
     .     TLAST,OB2000,OBLM,OBLT,OBL,X,DSIN,DCOS
      DIMENSION POS1(3), POS2(3), POS0(3)
      SAVE

      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( RADCON = PI / 180.D0           )

*     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
      DATA T0 / 2451545.00000000D0 /
      DATA TLAST / 0.D0 /,   OB2000 / 0.D0 /

*     T1 IS THE TDB JULIAN DATE
      CALL TIMES ( TJD, T,   SECDIF )
      T1 = TJD + SECDIF / 86400.D0

      IF ( TJD .EQ. 0.D0 ) THEN
*         CASE WHERE INPUT VECTOR IS IN ICRS SYSTEM
          CALL FRAME ( POS1, 1,   POS0 )
*         GET MEAN OBLIQUITY AT J2000.0 IF NECESSARY
          IF ( OB2000 .EQ. 0.D0 ) CALL ETILT ( T0,  OB2000, X, X, X, X )
          OBL = OB2000 * RADCON
      ELSE 
*         CASE WHERE INPUT VECTOR IS IN EQUATOR OF DATE SYSTEM      
          POS0(1) = POS1(1)
          POS0(2) = POS1(2)
          POS0(3) = POS1(3)
*         GET MEAN AND TRUE OBLIQUITY
          IF ( DABS ( TJD - TLAST ) .GT. 1.D-8 ) THEN
              CALL ETILT ( T1,   OBLM, OBLT, X, X, X )
              TLAST = TJD
          END IF
*         SELECT MEAN OR TRUE OBLIQUITY         
          OBL = OBLM * RADCON
          IF ( ICOORD .EQ. 1 ) OBL = OBLT * RADCON 
      END IF

*     ROTATE EQUATORIAL POSITION VECTOR TO ECLIPTIC SYSTEM
      POS2(1) =  POS0(1)
      POS2(2) =  POS0(2) * DCOS ( OBL ) + POS0(3) * DSIN ( OBL )
      POS2(3) = -POS0(2) * DSIN ( OBL ) + POS0(3) * DCOS ( OBL )

*     STORE COMPUTED POSITION VECTOR FOR POSSIBLE LATER RETRIEVAL
      CALL SETVEC ( POS2 )

      RETURN

      END
      
      
      
      SUBROUTINE ECEQ ( TJD, ICOORD, POS1,   POS2 )
*
*     THIS SUBROUTINE CONVERTS AN ECLIPTIC POSITION VECTOR TO
*     AN EQUATORIAL POSITION VECTOR.
*
*          TJD    = TT JULIAN DATE OF EQUATOR, EQUINOX, AND ECLIPTIC
*                   USED FOR COORDINATES (IN)
*          ICOORD = COORDINATE SYSTEM SELECTION (IN)
*                   SET ICOORD=0 FOR MEAN EQUATOR AND EQUINOX OF DATE
*                   SET ICOORD=1 FOR TRUE EQUATOR AND EQUINOX OF DATE
*                   (ECLIPTIC IS ALWAYS THE MEAN PLANE)
*          POS1   = POSITION VECTOR, REFERRED TO SPECIFIED
*                   ECLIPTIC AND EQUINOX OF DATE (IN)
*          POS2   = POSITION VECTOR, REFERRED TO SPECIFIED
*                   EQUATOR AND EQUINOX OF DATE (OUT)
*
*     NOTE:  TO CONVERT ECLIPTIC VECTORS (MEAN ECLIPTIC AND EQUINOX OF
*     OF J2000.0 ONLY) TO ICRS VECTORS, SET TJD = 0.D0 AND ICOORD = 0.
*     EXCEPT FOR THE OUTPUT FROM THIS CASE, ALL VECTORS ARE ASSUMED TO
*     BE WITH RESPECT TO A DYNAMICAL SYSTEM.
*
*
      DOUBLE PRECISION TJD,POS1,POS2,POS0,PI,RADCON,T0,T1,T,SECDIF,
     .     TLAST,OB2000,OBLM,OBLT,OBL,X,DSIN,DCOS
      DIMENSION POS1(3), POS2(3), POS0(3)
      SAVE

      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( RADCON = PI / 180.D0           )

*     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
      DATA T0 / 2451545.00000000D0 /
      DATA TLAST / 0.D0 /,   OB2000 / 0.D0 /

*     T1 IS THE TDB JULIAN DATE
      CALL TIMES ( TJD, T,   SECDIF )
      T1 = TJD + SECDIF / 86400.D0

      IF ( TJD .EQ. 0.D0 ) THEN
*         CASE WHERE OUTPUT VECTOR IS TO BE IN ICRS SYSTEM
*         GET MEAN OBLIQUITY AT J2000.0 IF NECESSARY 
          IF ( OB2000 .EQ. 0.D0 ) CALL ETILT ( T0,  OB2000, X, X, X, X )
          OBL = OB2000 * RADCON
      ELSE 
*         CASE WHERE OUTPUT VECTOR IS TO BE IN EQUATOR OF DATE SYSTEM      
*         GET MEAN AND TRUE OBLIQUITY
          IF ( DABS ( TJD - TLAST ) .GT. 1.D-8 ) THEN
              CALL ETILT ( T1,   OBLM, OBLT, X, X, X )
              TLAST = TJD
          END IF
*         SELECT MEAN OR TRUE OBLIQUITY         
          OBL = OBLM * RADCON
          IF ( ICOORD .EQ. 1 ) OBL = OBLT * RADCON 
      END IF     

*     ROTATE ECLIPTIC POSITION VECTOR TO EQUATORIAL SYSTEM
      POS2(1) =  POS1(1)
      POS2(2) =  POS1(2) * DCOS ( OBL ) - POS1(3) * DSIN ( OBL )
      POS2(3) =  POS1(2) * DSIN ( OBL ) + POS1(3) * DCOS ( OBL )

      IF ( TJD .EQ. 0.D0 ) THEN
*         CASE WHERE OUTPUT VECTOR IS TO BE IN ICRS SYSTEM
          POS0(1) = POS2(1) 
          POS0(2) = POS2(2)
          POS0(3) = POS2(3)
          CALL FRAME ( POS0, -1,   POS2 )
      END IF    

*     STORE COMPUTED POSITION VECTOR FOR POSSIBLE LATER RETRIEVAL
      CALL SETVEC ( POS2 )

      RETURN

      END
      


      SUBROUTINE EQGAL ( RA, DEC,   GLON, GLAT )
*
*     THIS SUBROUTINE CONVERTS ICRS RIGHT ASCENSION AND DECLINATION
*     TO GALACTIC LONGITUDE AND LATITUDE.  IT USES THE TRANSFORMATION
*     GIVEN IN THE HIPPARCOS AND TYCHO CATALOGUES, VOL. 1,
*     SECTION 1.5.3.
*
*          RA     = ICRS RIGHT ASCENSION IN HOURS (IN)
*          DEC    = ICRS DECLINATION IN DEGREES (IN)
*          GLON   = GALACTIC LONGITUDE IN DEGREES (OUT)
*          GLAT   = GALACTIC LATITUDE IN DEGREES (OUT)
*
*
      DOUBLE PRECISION RA,DEC,GLON,GLAT,PI,RADCON,AG,R,D,XYPROJ,G,
     .     POS1,POS2,DSIN,DCOS,DSQRT,DATAN2
      DIMENSION POS1(3), POS2(3), AG(3,3)
      SAVE

      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( RADCON = PI / 180.D0           )

*     ROTATION MATRIX A_G FROM HIPPARCOS DOCUMENTATION EQ. 1.5.11
      DATA AG /
     .     -0.0548755604D0, +0.4941094279D0, -0.8676661490D0,
     .     -0.8734370902D0, -0.4448296300D0, -0.1980763734D0,
     .     -0.4838350155D0, +0.7469822445D0, +0.4559837762D0 /

*     FORM POSITION VECTOR IN EQUATORIAL SYSTEM FROM INPUT COORDINATES
      R = RA * 15.D0 * RADCON
      D = DEC * RADCON
      POS1(1) = DCOS ( D ) * DCOS ( R )
      POS1(2) = DCOS ( D ) * DSIN ( R )
      POS1(3) = DSIN ( D )

*     ROTATE POSITION VECTOR TO GALACTIC SYSTEM, USING HIPPARCOS
*     DOCUMENTATION EQ. 1.5.13
      POS2(1) = AG(1,1)*POS1(1) + AG(1,2)*POS1(2) + AG(1,3)*POS1(3)
      POS2(2) = AG(2,1)*POS1(1) + AG(2,2)*POS1(2) + AG(2,3)*POS1(3)
      POS2(3) = AG(3,1)*POS1(1) + AG(3,2)*POS1(2) + AG(3,3)*POS1(3)

*     DECOMPOSE GALACTIC VECTOR INTO LONGITUDE AND LATITUDE
      XYPROJ = DSQRT ( POS2(1)**2 + POS2(2)**2 )
      G = 0.D0
      IF ( XYPROJ .GT. 0.D0 ) G = DATAN2 ( POS2(2), POS2(1) )
      GLON = G / RADCON
      IF ( GLON .LT. 0.D0 ) GLON = GLON + 360.D0
      G = DATAN2 ( POS2(3), XYPROJ )
      GLAT = G / RADCON

*     STORE COMPUTED POSITION VECTOR FOR POSSIBLE LATER RETRIEVAL
      CALL SETVEC ( POS2 )

      RETURN

      END



      SUBROUTINE VECTRS (RA,DEC,PMRA,PMDEC,PARLLX,RV,POS,VEL)
*
*     THIS SUBROUTINE CONVERTS ANGULAR QUANTITIES RELATED TO A STAR'S
*     POSITION AND MOTION TO VECTORS.
*
*          RA     = RIGHT ASCENSION IN HOURS (IN)
*          DEC    = DECLINATION IN DEGREES (IN)
*          PMRA   = PROPER MOTION IN RA IN MILLIARCSECONDS PER YEAR
*                   (IN)
*          PMDEC  = PROPER MOTION IN DEC IN MILLIARCSECONDS PER YEAR
*                   (IN)
*          PARLLX = PARALLAX IN MILLIARCSECONDS (IN)
*          RV     = RADIAL VELOCITY IN KILOMETERS/SECOND (IN)
*          POS    = POSITION VECTOR, EQUATORIAL RECTANGULAR COORDINATES,
*                   WITH RESPECT TO SOLAR SYSTEM BARYCENTER, COMPONENTS
*                   IN AU (OUT)
*          VEL    = VELOCITY VECTOR, EQUATORIAL RECTANGULAR COORDINATES,
*                   WITH RESPECT TO SOLAR SYSTEM BARYCENTER, COMPONENTS
*                   IN AU/DAY (OUT)
*
*
      DOUBLE PRECISION RA,DEC,PMRA,PMDEC,PARLLX,RV,POS,VEL,
     .     PI,SECCON,AUKM,C,PARALX,DIST,R,D,CRA,SRA,CDC,SDC,K,
     .     PMR,PMD,RVL,DCOS,DSIN
      DIMENSION POS(3), VEL(3)
      SAVE

      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( SECCON = 180.D0 * 3600.D0 / PI )

      DATA NTIMES / 0 /

      NTIMES = NTIMES + 1
      IF ( NTIMES .EQ. 1 ) THEN
*         GET LENGTH OF AU IN KILOMETERS
          CALL ASTCON ( 'AU', 1.D-3,   AUKM )
*         GET C, THE SPEED OF LIGHT IN KILOMETERS/SECOND
          CALL ASTCON ( 'C', 1.D-3,   C )
      END IF

*     IF PARALLAX IS UNKNOWN, UNDETERMINED, OR ZERO, SET IT TO 1E-6
*     MILLIARCSECOND, CORRESPONDING TO A DISTANCE OF 1 GIGAPARSEC
      PARALX = PARLLX
      IF ( PARALX .LE. 0.D0 ) PARALX = 1.D-6

*     CONVERT RIGHT ASCENSION, DECLINATION, AND PARALLAX TO POSITION
*     VECTOR IN EQUATORIAL SYSTEM WITH UNITS OF AU
      DIST = 1.D0 / DSIN ( PARALX * 1.D-3 / SECCON )
      R = RA * 54000.D0 / SECCON
      D = DEC * 3600.D0 / SECCON
      CRA = DCOS ( R )
      SRA = DSIN ( R )
      CDC = DCOS ( D )
      SDC = DSIN ( D )
      POS(1) = DIST * CDC * CRA
      POS(2) = DIST * CDC * SRA
      POS(3) = DIST * SDC

*     COMPUTE DOPPLER FACTOR, WHICH ACCOUNTS FOR CHANGE IN
*     LIGHT TRAVEL TIME TO STAR
      K = 1.D0 / ( 1.D0 - RV / C )

*     CONVERT PROPER MOTION AND RADIAL VELOCITY TO ORTHOGONAL COMPONENTS
*     OF MOTION WITH UNITS OF AU/DAY
      PMR = PMRA  / ( PARALX * 365.25D0 ) * K
      PMD = PMDEC / ( PARALX * 365.25D0 ) * K
      RVL = RV * 86400.D0 / AUKM          * K

*     TRANSFORM MOTION VECTOR TO EQUATORIAL SYSTEM
      VEL(1) = - PMR * SRA - PMD * SDC * CRA + RVL * CDC * CRA
      VEL(2) =   PMR * CRA - PMD * SDC * SRA + RVL * CDC * SRA
      VEL(3) =               PMD * CDC       + RVL * SDC

      RETURN

      END



      SUBROUTINE ANGLES (POS,RA,DEC)
*
*     THIS SUBROUTINE CONVERTS A VECTOR TO ANGULAR QUANTITIES.
*
*          POS = POSITION VECTOR, EQUATORIAL RECTANGULAR
*                COORDINATES (IN)
*          RA  = RIGHT ASCENSION IN HOURS (OUT)
*          DEC = DECLINATION IN DEGREES (OUT)
*
*
      DOUBLE PRECISION POS,RA,DEC,PI,SECCON,XYPROJ,R,D,DSQRT,DATAN2
      DIMENSION POS(3)

      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( SECCON = 180.D0 * 3600.D0 / PI )

      XYPROJ = DSQRT(POS(1)**2 + POS(2)**2)
      R = 0.D0
      IF (XYPROJ.GT.0.D0) R = DATAN2(POS(2),POS(1))
      RA = R * SECCON / 54000.D0
      IF (RA.LT. 0.D0) RA = RA + 24.D0
      IF (RA.GE.24.D0) RA = RA - 24.D0
      D = DATAN2(POS(3),XYPROJ)
      DEC = D * SECCON / 3600.D0

*     STORE COMPUTED POSITION VECTOR FOR POSSIBLE LATER RETRIEVAL
      CALL SETVEC (POS)

      RETURN

      END



      SUBROUTINE PROPMO (TJD1,POS1,VEL1,TJD2,POS2)
*
*     THIS SUBROUTINE APPLIES PROPER MOTION, INCLUDING FORESHORTENING
*     EFFECTS, TO A STAR'S POSITION.
*
*          TJD1 = TDB JULIAN DATE OF FIRST EPOCH (IN)
*          POS1 = POSITION VECTOR OF STAR AT FIRST EPOCH (IN)
*          VEL1 = VELOCITY VECTOR OF STAR AT FIRST EPOCH (IN)
*          TJD2 = TDB JULIAN DATE OF SECOND EPOCH (IN)
*          POS2 = POSITION VECTOR OF STAR AT SECOND EPOCH (OUT)
*
*
      DOUBLE PRECISION TJD1,POS1,VEL1,TJD2,POS2
      DIMENSION POS1(3), VEL1(3), POS2(3)

      DO 20 J=1,3
   20 POS2(J) = POS1(J) + VEL1(J) * (TJD2 - TJD1)

      RETURN

      END



      SUBROUTINE GEOCEN (POS1,PE,POS2,TLIGHT)
*
*     THIS SUBROUTINE MOVES THE ORIGIN OF COORDINATES FROM THE
*     BARYCENTER OF THE SOLAR SYSTEM TO THE OBSERVER (OR THE
*     GEOCENTER).  I.E., THIS SUBROUTINE ACCOUNTS FOR PARALLAX
*     (ANNUAL+GEOCENTRIC OR JUST ANNUAL).
*
*          POS1   = POSITION VECTOR OF STAR OR PLANET, WITH RESPECT TO
*                   ORIGIN AT SOLAR SYSTEM BARYCENTER, COMPONENTS
*                   IN AU (IN)
*          PE     = POSITION VECTOR OF OBSERVER (OR THE GEOCENTER),
*                   WITH RESPECT TO ORIGIN AT SOLAR SYSTEM BARYCENTER,
*                   COMPONENTS IN AU (IN)
*          POS2   = POSITION VECTOR OF STAR OR PLANET, WITH RESPECT TO
*                   ORIGIN AT OBSERVER (OR THE GEOCENTER), COMPONENTS
*                   IN AU (OUT)
*          TLIGHT = LIGHT-TIME FROM STAR OR PLANET TO OBSERVER (OR THE
*                   GEOCENTER) IN DAYS (OUT)
*
*     NOTE: STAR AND PLANET ARE USED GENERICALLY FOR BODIES OUTSIDE AND
*           INSIDE THE SOLAR SYSTEM, RESPECTIVELY.
*
      DOUBLE PRECISION POS1,PE,POS2,TLIGHT,C,DSQRT
      DIMENSION POS1(3), PE(3), POS2(3)
      SAVE

      DATA NTIMES / 0 /

      NTIMES = NTIMES + 1
      IF (NTIMES.EQ.1) THEN
*         GET C, THE SPEED OF LIGHT IN AU/DAY
          CALL ASTCON ('C(AU/DAY)',1.D0,C)
      END IF

      DO 20 J=1,3
   20 POS2(J) = POS1(J) - PE(J)
      TLIGHT = DSQRT(POS2(1)**2 + POS2(2)**2 + POS2(3)**2) / C

      RETURN

      END



      SUBROUTINE GEOPOS (TJD,LOCATN,OBSERV,POS,VEL)
*
*     THIS SUBROUTINE COMPUTES THE GEOCENTRIC POSITION AND VELOCITY
*     OF AN OBSERVER ON THE SURFACE OF THE EARTH OR ON A NEAR-EARTH
*     SPACECRAFT.  THE FINAL VECTORS ARE EXPRESSED IN THE GCRS.
*
*          TJD    = TT JULIAN DATE (IN)
*          LOCATN = INTEGER CODE SPECIFYING LOCATION OF OBSERVER (IN)
*                   SET LOCATN=0 FOR OBSERVER AT GEOCENTER
*                   SET LOCATN=1 FOR OBSERVER ON SURFACE OF EARTH
*                   SET LOCATN=2 FOR OBSERVER ON NEAR-EARTH SPACECRAFT
*          OBSERV = ARRAY OF DATA SPECIFYING LOCATION OF OBSERVER (IN)
*                   FOR LOCATN=0, THIS ARRAY NOT USED
*                   FOR LOCATN=1,
*                   OBSERV(1) = GEODETIC LONGITUDE (WGS-84) OF OBSERVER
*                               (EAST +) IN DEGREES (IN)
*                   OBSERV(2) = GEODETIC LATITUDE (WGS-84) OF OBSERVER
*                               (NORTH +) IN DEGREES (IN)
*                   OBSERV(3) = HEIGHT OF OBSERVER ABOVE ELLIPSOID
*                               IN METERS (IN)
*                   OBSERV(4) = VALUE OF DELTA-T IN SECONDS (IN)
*                               (DELTA-T=TT-UT1)
*                   OBSERV(5) = (NOT USED, RESERVED FOR FUTURE USE)
*                   OBSERV(6) = (NOT USED, RESERVED FOR FUTURE USE)
*                   FOR LOCATN=2,
*                   OBSERV(1) = GEOCENTRIC X IN KILOMETERS          
*                   OBSERV(2) = GEOCENTRIC Y IN KILOMETERS          
*                   OBSERV(3) = GEOCENTRIC Z IN KILOMETERS          
*                   OBSERV(4) = GEOCENTRIC X-DOT IN KILOMETERS/SECOND  
*                   OBSERV(5) = GEOCENTRIC Y-DOT IN KILOMETERS/SECOND  
*                   OBSERV(6) = GEOCENTRIC Z-DOT IN KILOMETERS/SECOND
*                   WITH RESPECT TO TRUE EQUATOR AND EQUINOX OF DATE
*          POS    = POSITION VECTOR OF OBSERVER, WITH RESPECT TO ORIGIN
*                   AT GEOCENTER, REFERRED TO GCRS AXES, COMPONENTS
*                   IN AU (OUT)
*          VEL    = VELOCITY VECTOR OF OBSERVER, WITH RESPECT TO ORIGIN
*                   AT GEOCENTER, REFERRED TO GCRS AXES, COMPONENTS
*                   IN AU/DAY (OUT)
*
*     NOTE 1: IF LOCATN=1 AND OBSERV(4)=0.D0, THE VALUE OF DELTA-T WILL
*     BE OBTAINED FROM GETDT, WHICH PROVIDES THE LAST VALUE OF DELTA-T
*     DEFINED BY USER VIA CALL TO SETDT.
*
*     NOTE 2: THIS SUBROUTINE CALLS SUBROUTINE TERRA FOR AN OBSERVER
*     ON THE SURFACE OF THE EARTH.  TERRA NEGLECTS POLAR MOTION, AN
*     APPROXIMATION WHICH MAY YIELD UP TO 15 METERS ERROR IN POSITION
*     AND SEVERAL MILLIMETERS/SEC ERROR IN VELOCITY.  
*
*
      DOUBLE PRECISION TJD,OBSERV,POS,VEL,T0,TLAST,
     .     AU,DELTAT,TTJD,TDBJD,UT1JD,ST,GST,GMST,GAST,EQEQ,X,
     .     POS1,VEL1,POS2,VEL2,POS3,VEL3

      DIMENSION OBSERV(6), POS(3), VEL(3),
     .     POS1(3), VEL1(3), POS2(3), VEL2(3), POS3(3), VEL3(3)

      SAVE

*     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
      DATA T0 / 2451545.00000000D0 /
      DATA TLAST / 0.D0 /,   GST / -99.D0 /,   NTIMES / 0 /

      NTIMES = NTIMES + 1
      IF ( NTIMES .EQ. 1 ) THEN
*         GET AU, THE LENGTH OF THE ASTRONOMICAL UNIT IN KILOMETERS
          CALL ASTCON ( 'AU', 1.D-3,   AU )
      END IF

      IF ( LOCATN .EQ. 0 ) THEN
          POS(1) = 0.D0
          POS(2) = 0.D0
          POS(3) = 0.D0
          VEL(1) = 0.D0
          VEL(2) = 0.D0
          VEL(3) = 0.D0
          GO TO 70
      END IF

      TTJD  = TJD
*     TDB IS APPROXIMATED BY TT      
      TDBJD = TJD

*     GET GEOCENTRIC POSITION AND VELOCITY VECTORS OF OBSERVER WRT
*     EQUATOR AND EQUINOX OF DATE

      IF ( LOCATN .EQ. 1 ) THEN

*         OBSERVER ON SURFACE OF EARTH

*         TEMPORARY CODE TO USE SIDEREAL TIME PREVIOUSLY PROVIDED
          IF ( GST .NE. -99.D0 ) THEN
              GAST = GST
              GST = -99.D0
              GO TO 20
          END IF
*         END OF TEMPROARY CODE

*         GET DELTA-T VALUE
          IF ( OBSERV(4) .NE. 0.D0 ) THEN
              DELTAT = OBSERV(4) / 86400.D0
          ELSE
              CALL GETDT ( DELTAT )
          END IF

*         USING DELTA-T VALUE, COMPUTE UT1 AND SIDEREAL TIME
          IF ( TTJD .EQ. 0.D0 ) THEN
              UT1JD = TDBJD - DELTAT
          ELSE
              UT1JD = TTJD - DELTAT
          END IF
          IF ( DABS ( UT1JD - TLAST ) .GT. 1.D-8 ) THEN
              CALL EQINOX
              CALL SIDTIM ( UT1JD, 0.D0, 0,   GMST )
              CALL ETILT ( TDBJD,   X, X, EQEQ, X, X )
              CALL RESUME
              TLAST = UT1JD
          END IF
          GAST = GMST + EQEQ / 3600.D0

*         SUBROUTINE TERRA DOES THE HARD WORK, GIVEN SIDEREAL TIME
  20      CALL TERRA ( OBSERV(1), OBSERV(2), OBSERV(3), GAST,
     .                 POS1, VEL1 )


      ELSE IF ( LOCATN .EQ. 2 ) THEN

*         OBSERVER ON NEAR-EARTH SPACECRAFT

*         CONVERT UNITS TO AU AND AU/DAY
          DO 30 J = 1, 3
              POS1(J) = OBSERV(J)   / AU
              VEL1(J) = OBSERV(J+3) / AU * 86400.D0
  30      CONTINUE

      END IF

*     TRANSFORM GEOCENTRIC POSITION VECTOR OF OBSERVER TO GCRS
      CALL NUTATE ( -TDBJD, POS1,   POS2 )
      CALL PRECES ( TDBJD, POS2, T0,   POS3 )
      CALL FRAME ( POS3, -1,   POS )

*     TRANSFORM GEOCENTRIC VELOCITY VECTOR OF OBSERVER TO GCRS
      CALL NUTATE ( -TDBJD, VEL1,   VEL2 )
      CALL PRECES ( TDBJD, VEL2, T0,   VEL3 )
      CALL FRAME ( VEL3, -1,   VEL )

  70  RETURN


*     TEMPORARY CODE FOR COMPATIBILITY WITH OLD ROUTINES
      ENTRY PLACST ( ST )
      GST = ST
      RETURN

      END



      SUBROUTINE LITTIM (TJD,IDBODY,POSE,TLITE,POS,TLIGHT)
*
*     THIS SUBROUTINE COMPUTES THE POSITION OF A SOLAR SYSTEM BODY,
*     AS ANTEDATED FOR LIGHT-TIME.
*
*          TJD    = TDB JULIAN DATE OF OBSERVATION (IN)
*          IDBODY = ID NUMBER OF BODY, USED IN CALLS TO SOLSYS (IN)
*          POSE   = POSITION VECTOR OF OBSERVER (OR THE GEOCENTER),
*                   WITH RESPECT TO ORIGIN AT SOLAR SYSTEM BARYCENTER,
*                   REFERRED TO ICRS AXES, COMPONENTS IN AU (IN)
*          TLITE  = FIRST APPROXIMATION TO LIGHT-TIME, IN DAYS (IN)
*                   (CAN BE SET TO 0.D0 IF UNKNOWN)
*          POS    = POSITION VECTOR OF BODY, WITH RESPECT TO ORIGIN AT
*                   OBSERVER (OR THE GEOCENTER), REFERRED TO ICRS AXES,
*                   COMPONENTS IN AU (OUT)
*          TLIGHT = FINAL LIGHT-TIME, IN DAYS (OUT)
*
*
      DOUBLE PRECISION TJD,POSE,TLITE,POS,TLIGHT,T0,T1,T2,T3,TOL,
     .     POS1,VEL1,DINT,DABS 
      LOGICAL SPLIT 

      DIMENSION POSE(3), POS(3), POS1(3), VEL1(3)
      
      SAVE NTIMES, SPLIT
      
      DATA NTIMES / 0 /,   SPLIT / .FALSE. /

   3  FORMAT ( ' LITTIM: PROBLEM WITH BODY NUMBER ', I3, ' AT JD ',
     .     F10.1 )

      NTIMES = NTIMES + 1

*     ON FIRST CALL, CHECK WHETHER SOLSYS SUPPORTS SPLIT JULIAN DATES
      IF ( NTIMES .EQ. 1 ) SPLIT = IDSS ('JD') .EQ. 2
      
*     SET LIGHT-TIME CONVERGENCE TOLERANCE
      TOL = 1.D-9
      IF ( SPLIT .AND. TLITE .LT. 0.01D0 ) TOL = 1.D-12 
      
*     IF SOLSYS SUPPORTS SPLIT JULIAN DATES, SPLIT THE JULIAN DATE
*     INTO WHOLE DAYS + FRACTION OF DAY
      T0 = 0.D0
      IF ( SPLIT ) T0 = DINT ( TJD )
      T1 = TJD - T0
      T2 = T1 - TLITE
      IF ( SPLIT ) CALL SOLSYS ( T0, IDBODY, 0,   POS1, VEL1, IERR)
      ITER = 0

*     ITERATE TO OBTAIN CORRECT LIGHT-TIME (USUALLY CONVERGES RAPIDLY)
  40  CALL SOLSYS ( T2, IDBODY, 0,   POS1, VEL1, IERR )
      CALL GEOCEN ( POS1, POSE,   POS, TLIGHT )
      IF ( IERR .NE. 0 ) THEN
          WRITE ( *, 3 ) IDBODY, T0 + T2
          GO TO 70
      END IF
      T3 = T1 - TLIGHT
      IF ( DABS ( T3 - T2 ) .GT. TOL ) THEN
          ITER = ITER + 1
          IF ( ITER .GT. 10 ) THEN
              WRITE ( *, 3 ) IDBODY, T0 + T3
              GO TO 70
          END IF
          T2 = T3
          GO TO 40
      END IF

  70  RETURN

      END
         


      SUBROUTINE DLIGHT (POS1,PE,DIFLT)
*
*     THIS SUBROUTINE RETURNS THE DIFFERENCE IN LIGHT-TIME, FOR A STAR,
*     BETWEEN THE BARYCENTER OF THE SOLAR SYSTEM AND THE OBSERVER (OR
*     THE GEOCENTER).
*
*          POS1   = POSITION VECTOR OF STAR, WITH RESPECT TO ORIGIN AT
*                   SOLAR SYSTEM BARYCENTER (IN)
*          PE     = POSITION VECTOR OF OBSERVER (OR THE GEOCENTER),
*                   WITH RESPECT TO ORIGIN AT SOLAR SYSTEM BARYCENTER,
*                   COMPONENTS IN AU (IN)
*          DIFLT  = DIFFERENCE IN LIGHT TIME, IN THE SENSE STAR TO
*                   BARYCENTER MINUS STAR TO EARTH, IN DAYS (OUT)
*
*     -OR-
*
*     THIS SUBROUTINE RETURNS THE LIGHT-TIME FROM THE OBSERVER (OR THE
*     GEOCENTER) TO A POINT ON A LIGHT RAY THAT IS CLOSEST TO A
*     SPECIFIC SOLAR SYSTEM BODY.
*
*          POS1   = POSITION VECTOR TOWARD OBSERVED OBJECT, WITH RESPECT
*                   TO ORIGIN AT OBSERVER (OR THE GEOCENTER) (IN)
*          PE     = POSITION VECTOR OF SOLAR SYSTEM BODY, WITH RESPECT
*                   TO ORIGIN AT OBSERVER (OR THE GEOCENTER), COMPONENTS
*                   IN AU (IN)
*          DIFLT  = LIGHT TIME TO POINT ON LINE DEFINED BY POS1 THAT IS
*                   CLOSEST TO SOLAR SYSTEM BODY (POSITIVE IF LIGHT
*                   PASSES BODY BEFORE HITTING OBSERVER, I.E., IF
*                   POS1 IS WITHIN 90 DEGREES OF PE)(OUT)
*
*
      DOUBLE PRECISION POS1,PE,DIFLT,C,DIS,U1,DSQRT
      DIMENSION POS1(3), PE(3), U1(3)
      SAVE

      DATA NTIMES / 0 /

      NTIMES = NTIMES + 1
      IF (NTIMES.EQ.1) THEN
*         GET C, THE SPEED OF LIGHT IN AU/DAY
          CALL ASTCON ('C(AU/DAY)',1.D0,C)
      END IF

*     FROM POS1, FORM UNIT VECTOR U1 IN DIRECTION OF STAR OR
*     LIGHT SOURCE
      DIS = DSQRT ( POS1(1)**2 + POS1(2)**2 + POS1(3)**2 )
      DO 20 J=1,3
   20 U1(J) = POS1(J) / DIS

*     LIGHT-TIME RETURNED IS THE PROJECTION OF VECTOR PE ONTO THE UNIT
*     VECTOR U1 (FORMED FROM POS1), DIVIDED BY THE SPEED OF LIGHT
      DIFLT = ( PE(1)*U1(1) + PE(2)*U1(2) + PE(3)*U1(3) ) / C

      RETURN

      END



      SUBROUTINE GRVDEF (TJD,LOC,POS1,POBS,POS2)
*
*     THIS SUBROUTINE COMPUTES THE TOTAL GRAVITATIONAL DEFLECTION OF
*     LIGHT FOR THE OBSERVED OBJECT DUE TO THE MAJOR GRAVITATING BODIES
*     IN THE SOLAR SYSTEM.  THIS SUBROUTINE VALID FOR AN OBSERVED BODY
*     WITHIN THE SOLAR SYSTEM AS WELL AS FOR A STAR.  SEE KLIONER
*     (2003), ASTRONOMICAL JOURNAL 125, 1580-1597, SECTION 6.
*
*          TJD    = TDB JULIAN DATE OF OBSERVATION
*          LOC    = CODE FOR LOCATION OF OBSERVER, DETERMINING
*                   WHETHER THE GRAVITATIONAL DEFLECTION DUE TO THE
*                   EARTH ITSELF IS APPLIED (IN)
*                   SET LOC=0 FOR NO EARTH DEFLECTION (NORMALLY MEANS
*                             OBSERVER IS AT GEOCENTER)
*                   SET LOC=1 TO ADD IN EARTH DEFLECTION (NORMALLY
*                             MEANS OBSERVER IS ON OR ABOVE SURFACE
*                             OF EARTH, INCLUDING EARTH ORBIT)
*          POS1   = POSITION VECTOR OF OBSERVED OBJECT, WITH RESPECT TO
*                   ORIGIN AT OBSERVER (OR THE GEOCENTER), REFERRED
*                   TO ICRS AXES, COMPONENTS IN AU (IN)
*          POBS   = POSITION VECTOR OF OBSERVER (OR THE GEOCENTER),
*                   WITH RESPECT TO ORIGIN AT SOLAR SYSTEM BARYCENTER,
*                   REFERRED TO ICRS AXES, COMPONENTS IN AU (IN)
*          POS2   = POSITION VECTOR OF OBSERVED OBJECT, WITH RESPECT TO
*                   ORIGIN AT OBSERVER (OR THE GEOCENTER), REFERRED
*                   TO ICRS AXES, CORRECTED FOR GRAVITATIONAL
*                   DEFLECTION, COMPONENTS IN AU (OUT)
*
*
      DOUBLE PRECISION TJD,POS1,POBS,POS2,C,RMASS,RMASSE,PBODY,VBODY,
     .     PBODYO,X,TLT,DLT,TCLOSE,DSQRT
      CHARACTER*3 NAME
      DIMENSION POS1(3), POBS(3), POS2(3), NAME(10), ID(10), RMASS(10),
     .     PBODY(3), VBODY(3), PBODYO(3)
      SAVE

*     THE FOLLOWING LIST OF NAMES IDENTIFIES WHICH GRAVITATING BODIES
*     (ASIDE FROM THE EARTH) ARE POTENTIALLY USED -- LIST IS TAKEN FROM
*     KLIONER'S TABLE 1, THE ORDER BASED ON AREA OF SKY AFFECTED (COL 2)
      DATA NAME / 'SUN', 'JUP', 'SAT', 'MOO', 'VEN', 'URA', 'NEP',
     .     3*'   ' /
*     CHANGE VALUE OF NBODY TO INCLUDE OR EXCLUDE GRAVITATING BODIES
*     (NBODY=0 MEANS NO DEFLECTION CALCULATED, NBODY=1 MEANS SUN ONLY,
*     NBODY=2 MEANS SUN + JUPITER, ETC.)
      DATA NBODY / 3 /

      DATA NTIMES / 0 /

      NTIMES = NTIMES + 1
      IF ( NTIMES .EQ. 1 ) THEN
*         GET C, THE SPEED OF LIGHT IN AU/DAY
          CALL ASTCON ( 'C(AU/DAY)', 1.D0, C )
*         GET ID NUMBERS AND RECIPROCAL MASSES OF GRAVITATING BODIES
          DO 20 I = 1, NBODY
              ID(I) = IDSS ( NAME(I) )
              CALL ASTCON ( 'MASS_'//NAME(I), 1.D0,   RMASS(I) )
  20      CONTINUE
          IDE = IDSS ( 'EARTH' )
          CALL ASTCON ( 'MASS_EARTH', 1.D0,   RMASSE )
      END IF

*     INITIALIZE OUTPUT VECTOR OF OBSERVED OBJECT TO EQUAL INPUT VECTOR
      DO 30 J = 1, 3
          POS2(J) = POS1(J)
  30  CONTINUE
*     OPTION FOR NO DEFLECTION
      IF ( NBODY .LE. 0 ) GO TO 50

*     COMPUTE LIGHT-TIME TO OBSERVED OBJECT
      TLT = DSQRT ( POS1(1)**2 + POS1(2)**2 + POS1(3)**2 ) / C

*     CYCLE THROUGH GRAVITATING BODIES
      DO 40 I = 1, NBODY

          IF ( ID(I) .EQ. -9999 ) GO TO 40

*         GET POSITION OF GRAVITATING BODY WRT SS BARYCENTER AT TIME TJD
          CALL SOLSYS ( TJD, ID(I), 0,   PBODY, VBODY, IERR )

*         GET POSITION OF GRAVITATING BODY WRT OBSERVER AT TIME TJD
          CALL GEOCEN ( PBODY, POBS,   PBODYO, X )

*         COMPUTE LIGHT-TIME FROM POINT ON INCOMING LIGHT RAY THAT
*         IS CLOSEST TO GRAVITATING BODY
          CALL DLIGHT ( POS2, PBODYO,   DLT )

*         GET POSITION OF GRAVITATING BODY WRT SS BARYCENTER AT TIME
*         WHEN INCOMING PHOTONS WERE CLOSEST TO IT
          TCLOSE = TJD
          IF ( DLT .GT. 0.D0 ) TCLOSE = TJD - DLT
          IF ( TLT .LT. DLT  ) TCLOSE = TJD - TLT
          CALL SOLSYS ( TCLOSE, ID(I), 0,   PBODY, VBODY, IERR )

*         COMPUTE DEFLECTION DUE TO GRAVITATING BODY
          CALL GRVD ( POS2, POBS, PBODY, RMASS(I),   POS2 )

  40  CONTINUE

*     IF OBSERVER IS NOT AT GEOCENTER, ADD IN DEFLECTION DUE TO EARTH
      IF ( LOC .NE. 0 ) THEN

*         GET POSITION OF EARTH WRT SS BARYCENTER AT TIME TJD
          CALL SOLSYS ( TJD, IDE, 0,   PBODY, VBODY, IERR )

*         COMPUTE DEFLECTION DUE TO EARTH
          CALL GRVD ( POS2, POBS, PBODY, RMASSE,   POS2 )

      END IF

  50  RETURN

      END



      SUBROUTINE GRVD (POS1,POBS,PBODY,RMASS,POS2)
*
*     THIS SUBROUTINE CORRECTS POSITION VECTOR FOR THE DEFLECTION
*     OF LIGHT IN THE GRAVITATIONAL FIELD OF AN ARBITRARY BODY.  ADAPTED
*     FROM MURRAY (1981) MON. NOTICES ROYAL AST. SOCIETY 195, 639-648.
*     SEE ALSO FORMULAE IN SECTION B OF THE ASTRONOMICAL ALMANAC, OR
*     KAPLAN ET AL. (1989) ASTRONOMICAL JOURNAL 97, 1197-1210, SECTION
*     III F.  THIS SUBROUTINE VALID FOR AN OBSERVED BODY WITHIN THE
*     SOLAR SYSTEM AS WELL AS FOR A STAR.
*
*          POS1   = POSITION VECTOR OF OBSERVED OBJECT, WITH RESPECT TO
*                   ORIGIN AT OBSERVER (OR THE GEOCENTER), COMPONENTS
*                   IN AU (IN)
*          POBS   = POSITION VECTOR OF OBSERVER (OR THE GEOCENTER),
*                   WITH RESPECT TO ORIGIN AT SOLAR SYSTEM BARYCENTER,
*                   COMPONENTS IN AU (IN)
*          PBODY  = POSITION VECTOR OF GRAVITATING BODY, WITH RESPECT TO
*                   ORIGIN AT SOLAR SYSTEM BARYCENTER, COMPONENTS
*                   IN AU (IN)
*          RMASS  = RECIPROCAL MASS OF GRAVITATING BODY IN SOLAR MASS
*                   UNITS, THAT IS, SUN MASS / BODY MASS (IN)
*          POS2   = POSITION VECTOR OF OBSERVED OBJECT, WITH RESPECT TO
*                   ORIGIN AT OBSERVER (OR THE GEOCENTER), CORRECTED FOR
*                   GRAVITATIONAL DEFLECTION, COMPONENTS IN AU (OUT)
*
*
      DOUBLE PRECISION POS1,POBS,PBODY,RMASS,POS2,C,MAU,GS,PQ,PE,
     .     PMAG,EMAG,QMAG,PHAT,EHAT,QHAT,PDOTQ,EDOTP,QDOTE,
     .     FAC1,FAC2,P2J,DABS,DSQRT
      DIMENSION POS1(3), POBS(3), PBODY(3), POS2(3), PQ(3), PE(3),
     .     PHAT(3), EHAT(3), QHAT(3)
      SAVE

      DATA NTIMES / 0 /

      NTIMES = NTIMES + 1
      IF (NTIMES.EQ.1) THEN
*         GET C, THE SPEED OF LIGHT IN METERS/SECOND
          CALL ASTCON ( 'C', 1.D0, C )
*         GET MAU, THE LENGTH OF THE AU IN METERS
          CALL ASTCON ( 'AU', 1.D0, MAU )
*         GET GS, THE HELIOCENTRIC GRAVITATIONAL CONSTANT
          CALL ASTCON ( 'GS', 1.D0, GS )
      END IF

*     CONSTRUCT VECTOR PQ FROM GRAVITATING BODY TO OBSERVED OBJECT AND
*     CONSTRUCT VECTOR PE FROM GRAVITATING BODY TO OBSERVER
      DO 20 J=1,3
          PQ(J) = POBS(J) + POS1(J) - PBODY(J)
          PE(J) = POBS(J) - PBODY(J)
  20  CONTINUE

*     COMPUTE VECTOR MAGNITUDES AND UNIT VECTORS
      PMAG = DSQRT (POS1(1)**2 + POS1(2)**2 + POS1(3)**2)
      EMAG = DSQRT (  PE(1)**2 +   PE(2)**2 +   PE(3)**2)
      QMAG = DSQRT (  PQ(1)**2 +   PQ(2)**2 +   PQ(3)**2)
      DO 30 J = 1, 3
          PHAT(J) = POS1(J) / PMAG
          EHAT(J) =   PE(J) / EMAG
          QHAT(J) =   PQ(J) / QMAG
  30  CONTINUE

*     COMPUTE DOT PRODUCTS OF VECTORS
      PDOTQ = PHAT(1)*QHAT(1) + PHAT(2)*QHAT(2) + PHAT(3)*QHAT(3)
      EDOTP = EHAT(1)*PHAT(1) + EHAT(2)*PHAT(2) + EHAT(3)*PHAT(3)
      QDOTE = QHAT(1)*EHAT(1) + QHAT(2)*EHAT(2) + QHAT(3)*EHAT(3)

*     IF GRAVITATING BODY IS OBSERVED OBJECT, OR IS ON A STRAIGHT LINE
*     TOWARD OR AWAY FROM OBSERVED OBJECT TO WITHIN 1 ARCSEC,
*     DEFLECTION IS SET TO ZERO
      IF ( DABS ( EDOTP ) .GT. 0.99999999999D0 ) THEN
          DO 35 J=1,3
              POS2(J) = POS1(J)
  35      CONTINUE
          GO TO 50
      END IF

*     COMPUTE SCALAR FACTORS
      FAC1 = 2.0D0 * GS / (C * C * EMAG * MAU * RMASS)
      FAC2 = 1.0D0 + QDOTE

*     CONSTRUCT CORRECTED POSITION VECTOR POS2
      DO 40 J = 1, 3
          P2J = PHAT(J) + FAC1 * (PDOTQ*EHAT(J) - EDOTP*QHAT(J)) / FAC2
          POS2(J) = P2J * PMAG
  40  CONTINUE

  50  RETURN

      END



      SUBROUTINE ABERAT (POS1,VE,TLIGHT,POS2)
*
*     THIS SUBROUTINE CORRECTS POSITION VECTOR FOR ABERRATION OF LIGHT.
*     ALGORITHM INCLUDES RELATIVISTIC TERMS.  ADAPTED FROM MURRAY (1981)
*     MON. NOTICES ROYAL AST. SOCIETY 195, 639-648.
*
*          POS1   = POSITION VECTOR OF OBSERVED OBJECT, WITH REESPECT TO
*                   ORIGIN AT OBSERVER (OR THE GEOCENTER), COMPONENTS
*                   IN AU (IN)
*          VE     = VELOCITY VECTOR OF OBSERVER (OR THE GEOCENTER),
*                   WITH RESPECT TO ORIGIN AT SOLAR SYSTEM BARYCENTER,
*                   COMPONENTS IN AU/DAY (IN)
*          TLIGHT = LIGHT TIME FROM BODY TO OBSERVER (OR THE GEOCENTER)
*                   IN DAYS (IN)
*                   IF TLIGHT = 0.D0, THIS SUBROUTINE WILL COMPUTE
*          POS2   = POSITION VECTOR OF OBSERVED OBJECT, WITH RESPECT TO
*                   ORIGIN AT OBSERVER (OR THE GEOCENTER), CORRECTED
*                   FOR ABERRATION, COMPONENTS IN AU (OUT)
*
*
      DOUBLE PRECISION POS1,VE,TLIGHT,POS2,C,TL,P1MAG,VEMAG,
     .     BETA,DOT,COSD,GAMMAI,P,Q,R,DSQRT
      DIMENSION POS1(3), VE(3), POS2(3)
      SAVE

      DATA NTIMES / 0 /

      NTIMES = NTIMES + 1
      IF (NTIMES.EQ.1) THEN
*         GET C, THE SPEED OF LIGHT IN AU/DAY
          CALL ASTCON ('C(AU/DAY)',1.D0,C)
      END IF

      TL = TLIGHT
      P1MAG = TL * C
      IF (TL.NE.0.D0) GO TO 20
      P1MAG = DSQRT(POS1(1)**2 + POS1(2)**2 + POS1(3)**2)
      TL = P1MAG / C
   20 VEMAG = DSQRT(VE(1)**2 + VE(2)**2 + VE(3)**2)
      BETA = VEMAG / C
      DOT = POS1(1)*VE(1) + POS1(2)*VE(2) + POS1(3)*VE(3)
      COSD = DOT / (P1MAG * VEMAG)
      GAMMAI = DSQRT(1.D0 - BETA**2)
      P = BETA * COSD
      Q = (1.D0 + P / (1.D0 + GAMMAI)) * TL
      R = 1.D0 + P

      DO 30 J=1,3
   30 POS2(J) = (GAMMAI * POS1(J) + Q * VE(J)) / R

      RETURN

      END

   
      
      SUBROUTINE RADVL ( POS, VEL, VELOBS, STAR, DIST,   RV )
*
*     THIS SUBROUTINE PREDICTS THE RADIAL VELOCITY OF THE OBSERVED
*     OBJECT AS IT WOULD BE MEASURED BY SPECTROSCOPIC MEANS.  RADIAL
*     VELOCITY IS HERE DEFINED AS THE RADIAL VELOCITY MEASURE (Z)
*     TIMES THE SPEED OF LIGHT.  FOR A SOLAR SYSTEM BODY, IT APPLIES
*     TO A FICTITIOUS EMITTER AT THE CENTER OF THE OBSERVED OBJECT,
*     ASSUMED MASSLESS (NO GRAVITATIONAL RED SHIFT), AND DOES NOT
*     IN GENERAL APPLY TO REFLECTED LIGHT.  FOR STARS, IT INCLUDES
*     ALL EFFECTS, SUCH AS GRAVITATIONAL RED SHIFT, CONTAINED
*     IN THE CATALOG BARYCENTRIC RADIAL VELOCITY MEASURE, A SCALAR
*     DERIVED FROM SPECTROSCOPY.  NEARBY STARS WITH A KNOWN KINEMATIC
*     VELOCITY VECTOR (OBTAINED INDEPENDENTLY OF SPECTROSCOPY) CAN BE
*     TREATED LIKE SOLAR SYSTEM OBJECTS.  SEE LINDEGREN & DRAVINS
*     (2003), ASTRONOMY & ASTROPHYSICS 401, 1185-1201.
*     
*          POS    = GEOMETRIC POSITION VECTOR OF OBJECT WITH RESPECT TO
*                   OBSERVER, CORRECTED FOR LIGHT-TIME, IN AU (IN)
*          VEL    = VELOCITY VECTOR OF OBJECT WITH RESPECT TO SOLAR 
*                   SYSTEM BARYCENTER, COMPONENTS IN AU/DAY (IN)
*          VELOBS = VELOCITY VECTOR OF OBSERVER WITH RESPECT TO SOLAR
*                   SYSTEM BARYCENTER, COMPONENTS IN AU/DAY (IN)
*          STAR   = 3-ELEMENT ARRAY OF CATALOG DATA FOR A STAR, TO BE
*                   NON-ZERO IF OBSERVED OBJECT IS A STAR FOR WHICH THE
*                   CATALOG RADIAL VELOCITY IS CONSISTENT WITH
*                   THE IAU DEFINITION OF BARYCENTRIC RADIAL VELOCITY
*                   MEASURE (OTHERWISE ALL ELEMENTS SHOULD BE SET TO 
*                   0.D0 EXACTLY) (IN)
*                   STAR(1) = CATALOG RA IN HOURS 
*                   STAR(2) = CATALOG DEC IN DEGREES
*                   STAR(3) = Z*C, THE CATALOG BARYCENTRIC RADIAL
*                             VELOCITY MEASURE TIMES THE SPEED OF LIGHT,
*                             IN KILOMETERS/SECOND
*                   ALL THREE DATA ELEMENTS MUST APPLY TO THE SAME 
*                   EPOCH (USUALLY J2000.0 = JD 2451545.0 TT)
*          DIST   = 3-ELEMENT ARRAY OF DISTANCES IN AU (IN)
*                   DIST(1) = DISTANCE OF OBSERVER FROM THE GEOCENTER
*                   DIST(2) = DISTANCE OF OBSERVER FROM THE SUN
*                   DIST(3) = DISTANCE OF OBJECT FROM THE SUN
*          RV     = THE OBSERVED RADIAL VELOCITY MEASURE TIMES
*                   THE SPEED OF LIGHT, IN KILOMETERS/SECOND (OUT)
*
*     NOTE 1:  ALL THE INPUT ARGUMENTS ARE BCRS QUANTITIES, EXPRESSED
*     WITH RESPECT TO THE ICRS AXES.  VEL AND VELOBS ARE KINEMATIC
*     VELOCITIES -- DERIVED FROM GEOMETRY OR DYNAMICS, NOT SPECTROSCOPY.
*
*     NOTE 2:  IF ANY ELEMENT OF ARRAY STAR IS NON-ZERO, THE ALGORITHM
*     USED WILL BE CONSISTENT WITH THE IAU DEFINITION OF STELLAR
*     RADIAL VELOCITY, SPECIFICALLY, THE BARYCENTRIC RADIAL VELOCITY
*     MEASURE, WHICH IS DERIVED FROM SPECTROSCOPY.  IN THAT CASE,
*     THE VECTOR VEL CAN BE VERY APPROXIMATE -- OR, FOR DISTANT STARS
*     OR GALAXIES, ZERO -- AS IT WILL BE USED ONLY FOR A SMALL GEOMETRIC
*     CORRECTION THAT IS PROPORTIONAL TO PROPER MOTION.
*
*     NOTE 3:  ANY OF THE DISTANCES IN ARRAY DIST CAN BE SET TO ZERO
*     (0.D0) IF THE CORRESPONDING GENERAL RELATIVISTIC GRAVITATIONAL
*     POTENTIAL TERM IS NOT TO BE EVALUATED.  THESE TERMS 
*     GENERALLY ARE IMPORTANT ONLY AT THE METER/SECOND LEVEL.  IF
*     THE FIRST TWO DISTANCES ARE BOTH ZERO, AN AVERAGE VALUE 
*     WILL BE USED FOR THE RELATIVISTIC TERM FOR THE OBSERVER, 
*     APPROPRIATE FOR AN OBSERVER ON THE SURFACE OF THE EARTH.  THE
*     THIRD DISTANCE, IF GIVEN, IS USED ONLY FOR SOLAR SYSTEM OBJECTS.
*
*     
      DOUBLE PRECISION POS,VEL,VELOBS,STAR,DIST,RV,PI,RADCON,
     .     AU,C,GS,GE,C2,TOMS,TOMS2,
     .     POSMAG,UK,V2,VO2,R,PHIGEO,PHISUN,REL,ZC,RA,DC,DU,
     .     ZB1,KVOBS,KV,ZOBS1,DSQRT,DCOS,DSIN
     
      LOGICAL DOSTAR
      
      DIMENSION POS(3), VEL(3), VELOBS(3), STAR(3), DIST(3), UK(3),
     .     DU(3) 
      
      SAVE
      
      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( RADCON = PI / 180.D0           )

      DATA NTIMES / 0 /
      
      NTIMES = NTIMES + 1

      IF ( NTIMES .EQ. 1 ) THEN
*         GET AU, LENGTH OF ASTRONOMICAL UNIT IN METERS          
          CALL ASTCON ( 'AU', 1.D0,   AU )
*         GET C, THE SPEED OF LIGHT IN METERS/SECOND
          CALL ASTCON ( 'C', 1.D0,   C )
*         GET GS, HELIOCENTRIC GRAVITATIONAL CONSTANT
          CALL ASTCON ( 'GS', 1.D0,   GS )
*         GET GE, GEOCENTRIC GRAVITATIONAL CONSTANT
          CALL ASTCON ( 'GE', 1.D0,   GE )
*         (GS AND GE ARE IN METERS**3/SECOND**2)
          C2 = C**2
          TOMS = AU / 86400.D0 
          TOMS2 = TOMS**2
      END IF
      
      RV = 0.D0

*     COMPUTE LENGTH OF POSITION VECTOR = DISTANCE TO OBJECT IN AU     
      POSMAG = DSQRT ( POS(1)**2 + POS(2)**2 + POS(3)**2 )
      IF ( POSMAG .LT. 1.D-8 ) GO TO 70

*     DETERMINE HOW OBJECT IS TO BE PROCESSED
      DOSTAR = STAR(1) .NE. 0.D0 .OR.
     .         STAR(2) .NE. 0.D0 .OR.
     .         STAR(3) .NE. 0.D0 

*     COMPUTE UNIT VECTOR TOWARD OBJECT      
      DO 40 J = 1, 3
          UK(J) = POS(J) / POSMAG
  40  CONTINUE    

*     COMPUTE VELOCITY-SQUARED FACTORS
      V2  = ( VEL(1)   **2 + VEL(2)   **2 + VEL(3)   **2 ) * TOMS2
      VO2 = ( VELOBS(1)**2 + VELOBS(2)**2 + VELOBS(3)**2 ) * TOMS2

*     COMPUTE GEOPOTENTIAL AT OBSERVER, UNLESS OBSERVER IS GEOCENTRIC
      R = DIST(1) * AU
      PHIGEO = 0.D0
      IF ( R .GT. 1.D6 ) PHIGEO = GE / R

*     COMPUTE SOLAR POTENTIAL AT OBSERVER
      R = DIST(2) * AU
      PHISUN = 0.D0
      IF ( R .GT. 1.D8 ) PHISUN = GS / R
     
*     COMPUTE RELATIVISTIC POTENTIAL AND VELOCITY FACTOR FOR OBSERVER
      IF ( DIST(1) .NE. 0.D0 .OR. DIST(2) .NE. 0.D0 ) THEN 
*         LINDEGREN & DRAVINS EQ. (41), SECOND FACTOR IN PARENTHESES
          REL = 1.D0 - ( PHIGEO + PHISUN ) / C2 - 0.5D0 * VO2 / C2 
      ELSE    
*         LINDEGREN & DRAVINS EQ. (42), INVERSE
          REL = 1.D0 - 1.550D-8
      END IF
      
      IF ( DOSTAR ) THEN

*         FOR STARS, UPDATE BARYCENTRIC RADIAL VELOCITY MEASURE FOR
*         CHANGE IN VIEW ANGLE
          RA = STAR(1) * 15.D0 * RADCON
          DC = STAR(2) * RADCON
          DU(1) = UK(1) - ( DCOS ( DC ) * DCOS ( RA ) )
          DU(2) = UK(2) - ( DCOS ( DC ) * DSIN ( RA ) )
          DU(3) = UK(3) - ( DSIN ( DC )               ) 
          ZC = STAR(3) * 1.D3 + 
     .       ( VEL(1) * DU(1) + VEL(2) * DU(2) + VEL(3) * DU(3) ) * TOMS
          
*         COMPUTE OBSERVED RADIAL VELOCITY MEASURE OF A STAR (INVERSE OF
*         LINDEGREN & DRAVINS EQ. (41))
          ZB1 = 1.D0 + ZC / C
          KVOBS = ( UK(1) * VELOBS(1) 
     .            + UK(2) * VELOBS(2)  
     .            + UK(3) * VELOBS(3) ) * TOMS     
          ZOBS1 = ZB1 * REL / ( 1.D0 + KVOBS / C )

      ELSE

*         COMPUTE SOLAR POTENTIAL AT OBJECT, IF WITHIN SOLAR SYSTEM
          R = DIST(3) * AU
          PHISUN = 0.D0
          IF ( R .GT. 1.D8 .AND. R .LT. 1.D16 ) PHISUN = GS / R

*         COMPUTE OBSERVED RADIAL VELOCITY MEASURE OF A PLANET OR OTHER
*         OBJECT -- INCLUDING A NEARBY STAR -- WHERE KINEMATIC 
*         BARYCENTRIC VELOCITY VECTOR IS KNOWN AND GRAVITATIONAL
*         RED SHIFT IS NEGLIGIBLE (LINDEGREN & DRAVINS EQ. (40),
*         APPLIED AS PER S. KLIONER PRIVATE COMMUNICATION (2006))
          KV = ( UK(1) * VEL(1)
     .         + UK(2) * VEL(2)
     .         + UK(3) * VEL(3) ) * TOMS
          ZB1 = ( 1.D0 + KV / C ) /
     .          ( 1.D0 - PHISUN / C2 - 0.5D0 * V2 / C2 )
          KVOBS = ( UK(1) * VELOBS(1)
     .            + UK(2) * VELOBS(2)
     .            + UK(3) * VELOBS(3) ) * TOMS 
          ZOBS1 = ZB1 * REL / ( 1.D0 + KVOBS / C )

      END IF          

*     CONVERT OBSERVED RADIAL VELOCITY MEASURE TO KILOMETERS/SECOND 
      RV = ( ZOBS1 - 1.D0 ) * C / 1000.D0
     
  70  RETURN
      
      END
      
      
      
      SUBROUTINE PRECES (TJD1,POS1,TJD2,POS2)
*
*     THIS SUBROUTINE PRECESSES EQUATORIAL RECTANGULAR COORDINATES FROM
*     ONE EPOCH TO ANOTHER.  THE COORDINATES ARE REFERRED TO THE MEAN
*     DYNAMICAL EQUATOR AND EQUINOX OF THE TWO RESPECTIVE EPOCHS.  SEE
*     EXPLANATORY SUPPLEMENT TO THE ASTRONOMICAL ALMANAC, PP. 103-104,
*     AND CAPITAINE ET AL. (2003), ASTRONOMY AND ASTROPHYSICS 412,
*     567-586.
*
*          TJD1 = TDB JULIAN DATE OF FIRST EPOCH (IN)
*          POS1 = POSITION VECTOR, GEOCENTRIC EQUATORIAL RECTANGULAR
*                 COORDINATES, REFERRED TO MEAN DYNAMICAL EQUATOR AND
*                 EQUINOX OF FIRST EPOCH (IN)
*          TJD2 = TDB JULIAN DATE OF SECOND EPOCH (IN)
*          POS2 = POSITION VECTOR, GEOCENTRIC EQUATORIAL RECTANGULAR
*                 COORDINATES, REFERRED TO MEAN DYNAMICAL EQUATOR AND
*                 EQUINOX OF SECOND EPOCH (OUT)
*
*     NOTE:  EITHER TJD1 OR TJD2 MUST BE 2451545.0 (J2000.0) TDB.
*
*
      DOUBLE PRECISION TJD1,TJD2,POS1,POS2,PI,SECCON,T0,TLAST,T,
     .     EPS0,PSIA,OMEGAA,CHIA,SA,CA,SB,CB,SC,CC,SD,CD,
     .     XX,YX,ZX,XY,YY,ZY,XZ,YZ,ZZ,DABS,DCOS,DSIN
      DIMENSION POS1(3), POS2(3)
      SAVE

      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( SECCON = 180.D0 * 3600.D0 / PI )

*     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
      DATA T0 / 2451545.00000000D0 /
      DATA TLAST / 0.D0 /
      
*     INITIALIZE PRECESSION ROTATION MATRIX AS IDENTITY MATRIX
      DATA XX, XY, XZ / 1.D0, 0.D0, 0.D0 /
      DATA YX, YY, YZ / 0.D0, 1.D0, 0.D0 /
      DATA ZX, ZY, ZZ / 0.D0, 0.D0, 1.D0 /

   3  FORMAT ( ' PRECES ERROR: PRECESSION FROM JD ', F10.1, ' TO ',
     .     F10.1, ' NOT TO/FROM J2000' )

      IF ( TJD1 .NE. T0 .AND. TJD2 .NE. T0 ) THEN
          WRITE ( *, 3 ) TJD1, TJD2
          GO TO 50
      END IF

*     T IS TIME IN TDB CENTURIES BETWEEN THE TWO EPOCHS
      T = ( TJD2 - TJD1 ) / 36525.D0
      IF ( TJD2 .EQ. T0 ) T = -T
      IF ( DABS ( T - TLAST ) .LT. 1.D-15 ) GO TO 20

*     NUMERICAL COEFFICIENTS OF PSI_A, OMEGA_A, AND CHI_A, ALONG WITH
*     EPSILON_0, THE OBLIQUITY AT J2000.0, ARE 4-ANGLE FORMULATION
*     FROM CAPITAINE ET AL. (2003), EQS. (4), (37), & (39)
      EPS0   = 84381.406D0
      PSIA   = ( ( ( ( -    0.0000000951D0   * T
     .                 +    0.000132851D0  ) * T
     .                 -    0.00114045D0   ) * T
     .                 -    1.0790069D0    ) * T
     .                 + 5038.481507D0     ) * T
      OMEGAA = ( ( ( ( +    0.0000003337D0   * T
     .                 -    0.000000467D0  ) * T
     .                 -    0.00772503D0   ) * T
     .                 +    0.0512623D0    ) * T
     .                 -    0.025754D0     ) * T + EPS0
      CHIA   = ( ( ( ( -    0.0000000560D0   * T
     .                 +    0.000170663D0  ) * T
     .                 -    0.00121197D0   ) * T
     .                 -    2.3814292D0    ) * T
     .                 +   10.556403D0     ) * T
      EPS0 = EPS0 / SECCON
      PSIA = PSIA / SECCON
      OMEGAA = OMEGAA / SECCON
      CHIA = CHIA / SECCON
      SA = DSIN ( EPS0 )
      CA = DCOS ( EPS0 )
      SB = DSIN ( -PSIA )
      CB = DCOS ( -PSIA )
      SC = DSIN ( -OMEGAA )
      CC = DCOS ( -OMEGAA )
      SD = DSIN ( CHIA )
      CD = DCOS ( CHIA )

*     COMPUTE ELEMENTS OF PRECESSION ROTATION MATRIX
*     EQUIVALENT TO R3(CHI_A)R1(-OMEGA_A)R3(-PSI_A)R1(EPSILON_0)
      XX =  CD * CB - SB * SD * CC
      YX =  CD * SB * CA + SD * CC * CB * CA - SA * SD * SC
      ZX =  CD * SB * SA + SD * CC * CB * SA + CA * SD * SC
      XY = -SD * CB - SB * CD * CC
      YY = -SD * SB * CA + CD * CC * CB * CA - SA * CD * SC
      ZY = -SD * SB * SA + CD * CC * CB * SA + CA * CD * SC
      XZ =  SB * SC
      YZ = -SC * CB * CA - SA * CC
      ZZ = -SC * CB * SA + CC * CA

      TLAST = T

   20 IF ( TJD2 .EQ. T0 ) GO TO 30

*     PERFORM ROTATION FROM J2000.0 TO EPOCH
      POS2(1) = XX * POS1(1) + YX * POS1(2) + ZX * POS1(3)
      POS2(2) = XY * POS1(1) + YY * POS1(2) + ZY * POS1(3)
      POS2(3) = XZ * POS1(1) + YZ * POS1(2) + ZZ * POS1(3)
      GO TO 50

*     PERFORM ROTATION FROM EPOCH TO J2000.0
   30 POS2(1) = XX * POS1(1) + XY * POS1(2) + XZ * POS1(3)
      POS2(2) = YX * POS1(1) + YY * POS1(2) + YZ * POS1(3)
      POS2(3) = ZX * POS1(1) + ZY * POS1(2) + ZZ * POS1(3)

   50 RETURN

      END



      SUBROUTINE NUTATE (TJD,POS1,POS2)
*
*     THIS SUBROUTINE NUTATES EQUATORIAL RECTANGULAR COORDINATES FROM
*     THE MEAN DYNAMICAL EQUATOR AND EQUINOX OF EPOCH TO THE TRUE
*     EQUATOR AND EQUINOX OF EPOCH.  SEE EXPLANATORY SUPPLEMENT TO THE
*     ASTRONOMICAL ALMANAC, PP. 114-115.
*
*          TJD    = TDB JULIAN DATE OF EPOCH (IN)
*          POS1   = POSITION VECTOR, GEOCENTRIC EQUATORIAL RECTANGULAR
*                   COORDINATES, REFERRED TO MEAN DYNAMICAL EQUATOR AND
*                   EQUINOX OF EPOCH (IN)
*          POS2   = POSITION VECTOR, GEOCENTRIC EQUATORIAL RECTANGULAR
*                   COORDINATES, REFERRED TO TRUE EQUATOR AND EQUINOX
*                   OF EPOCH (OUT)
*
*     NOTE:  IF TJD IS NEGATIVE, INVERSE NUTATION (TRUE TO MEAN)
*     IS APPLIED.
*
*
      DOUBLE PRECISION TJD,POS1,POS2,TJD1,PI,SECCON,OBLM,OBLT,EQEQ,
     .     DPSI,DEPS,COBM,SOBM,COBT,SOBT,CPSI,SPSI,
     .     XX,YX,ZX,XY,YY,ZY,XZ,YZ,ZZ,DABS,DCOS,DSIN
      DIMENSION POS1(3), POS2(3)

      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( SECCON = 180.D0 * 3600.D0 / PI )

      TJD1 = DABS(TJD)

      CALL ETILT ( TJD1,   OBLM, OBLT, EQEQ, DPSI, DEPS )
      OBLM = OBLM * 3600.D0 / SECCON
      OBLT = OBLT * 3600.D0 / SECCON
      DPSI = DPSI / SECCON
      DEPS = DEPS / SECCON
      COBM = DCOS ( OBLM )
      SOBM = DSIN ( OBLM )
      COBT = DCOS ( OBLT )
      SOBT = DSIN ( OBLT )
      CPSI = DCOS ( DPSI )
      SPSI = DSIN ( DPSI )

*     COMPUTE ELEMENTS OF NUTATION ROTATION MATRIX
      XX =  CPSI
      YX = -SPSI * COBM
      ZX = -SPSI * SOBM
      XY =  SPSI * COBT
      YY =  CPSI * COBM * COBT + SOBM * SOBT
      ZY =  CPSI * SOBM * COBT - COBM * SOBT
      XZ =  SPSI * SOBT
      YZ =  CPSI * COBM * SOBT - SOBM * COBT
      ZZ =  CPSI * SOBM * SOBT + COBM * COBT
   10 IF ( TJD .LT. 0.D0 ) GO TO 30

*     PERFORM ROTATION FROM MEAN TO TRUE
   20 POS2(1) = XX * POS1(1) + YX * POS1(2) + ZX * POS1(3)
      POS2(2) = XY * POS1(1) + YY * POS1(2) + ZY * POS1(3)
      POS2(3) = XZ * POS1(1) + YZ * POS1(2) + ZZ * POS1(3)
      GO TO 50

*     PERFORM ROTATION FROM TRUE TO MEAN
   30 POS2(1) = XX * POS1(1) + XY * POS1(2) + XZ * POS1(3)
      POS2(2) = YX * POS1(1) + YY * POS1(2) + YZ * POS1(3)
      POS2(3) = ZX * POS1(1) + ZY * POS1(2) + ZZ * POS1(3)

   50 RETURN

      END



      SUBROUTINE SPIN (ANGL,POS1,POS2)
*
*     THIS SUBROUTINE TRANSFORMS A VECTOR FROM ONE COORDINATE SYSTEM
*     TO ANOTHER WITH SAME ORIGIN AND AXES ROTATED ABOUT THE
*     Z AXIS.
*
*          ANGL   = ANGLE OF COORDINATE SYSTEM ROTATION, POSITIVE
*                   COUNTERCLOCKWISE WHEN VIEWED FROM +Z,
*                   IN DEGREES (IN)
*          POS1   = POSITION VECTOR (IN)
*          POS2   = POSITION VECTOR EXPRESSED IN NEW COORDINATE
*                   SYSTEM ROTATED ABOUT Z BY ANGLE ANG (OUT)
*
*
      DOUBLE PRECISION ANGL,POS1,POS2,PI,ALAST,ANG,COSANG,SINANG,
     .     XX,YX,ZX,XY,YY,ZY,XZ,YZ,ZZ,DABS,DCOS,DSIN
      DIMENSION POS1(3), POS2(3)
      SAVE

      PARAMETER ( PI     = 3.14159265358979324D0 )

      DATA ALAST / -999.D0 /

      IF ( DABS ( ANGL - ALAST ) .GT. 1.D-12 ) THEN

          ANG = ANGL / 180.D0 * PI
          COSANG = DCOS ( ANG )
          SINANG = DSIN ( ANG )

*         ROTATION MATRIX FOLLOWS
          XX =  COSANG
          YX =  SINANG
          ZX =  0.D0
          XY = -SINANG
          YY =  COSANG
          ZY =  0.D0
          XZ =  0.D0
          YZ =  0.D0
          ZZ =  1.D0

          ALAST = ANGL

      END IF

*     PERFORM ROTATION
      POS2(1) = XX * POS1(1) + YX * POS1(2) + ZX * POS1(3)
      POS2(2) = XY * POS1(1) + YY * POS1(2) + ZY * POS1(3)
      POS2(3) = XZ * POS1(1) + YZ * POS1(2) + ZZ * POS1(3)

      RETURN

      END



      SUBROUTINE WOBBLE (TJD,XP,YP,POS1,POS2)
*
*     THIS SUBROUTINE CORRECTS A VECTOR IN THE ITRS (A ROTATING EARTH-
*     FIXED SYSTEM) FOR POLAR MOTION, AND ALSO CORRECTS THE LONGITUDE
*     ORIGIN (BY A TINY AMOUNT) TO THE TERRESTRIAL INTERMEDIATE ORIGIN
*     (TIO).  THE ITRS VECTOR IS THEREBY TRANSFORMED TO THE TERRESTRIAL
*     INTERMEDIATE SYSTEM, BASED ON THE TRUE (ROTATIONAL) EQUATOR AND
*     THE TERRESTRIAL INTERMEDIATE ORIGIN (TIO).  SINCE THE TRUE EQUATOR
*     IS THE PLANE ORTHOGONAL TO THE DIRECTION OF THE CELESTIAL
*     INTERMEDIATE POLE (CIP), THE COMPONENTS OF THE OUTPUT VECTOR ARE
*     REFERRED TO Z AND X AXES TOWARD THE CIP AND TIO, RESPECTIVELY.
*
*          TJD    = TT OR UT1 JULIAN DATE (IN)
*          XP     = CONVENTIONALLY-DEFINED X COORDINATE OF CELESTIAL
*                   INTERMEDIATE POLE WITH RESPECT TO ITRS POLE, IN
*                   ARCSECONDS (IN)
*          YP     = CONVENTIONALLY-DEFINED Y COORDINATE OF CELESTIAL
*                   INTERMEDIATE POLE WITH RESPECT TO ITRS POLE, IN
*                   ARCSECONDS (IN)
*          POS1   = POSITION VECTOR, GEOCENTRIC EQUATORIAL RECTANGULAR
*                   COORDINATES, REFERRED TO ITRS AXES (IN)
*          POS2   = POSITION VECTOR, GEOCENTRIC EQUATORIAL RECTANGULAR
*                   COORDINATES, REFERRED TO TRUE EQUATOR AND TIO (OUT)
*
*     NOTE 1:  IF TJD IS NEGATIVE, THE INVERSE TRANSFORMATION (TERRESTRIAL 
*     INTERMEDIATE SYSTEM TO ITRS) IS APPLIED.
*
*     NOTE 2: INPUT PARAMETERS XP, YP WERE X, Y IN NOVAS F3.0.
*     THE NAMES WERE CHANGED FOR CONSISTANCY THROUGHOUT NOVAS AND WITH
*     IERS CONVENTIONS.
*
      DOUBLE PRECISION TJD,XP,YP,POS1,POS2,PI,SECCON,T0,T,XPOLE,YPOLE,
     .     SPRIME,TIOLON,SINX,COSX,SINY,COSY,SINL,COSL,
     .     XX,YX,ZX,XY,YY,ZY,XZ,YZ,ZZ,DABS,DSIN,DCOS
      DIMENSION POS1(3), POS2(3)

      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( SECCON = 180.D0 * 3600.D0 / PI )

*     T0 = TT JULIAN DATE OF EPOCH J2000.0 (TT)
      DATA T0 / 2451545.0D0 /

      XPOLE = XP / SECCON
      YPOLE = YP / SECCON

      T = ( DABS(TJD) - T0 ) / 36525.D0

*     COMPUTE APPROXIMATE LONGITUDE OF TIO, USING EQ. (10) OF
*     LAMBERT & BIZOUARD (2002), ASTRONOMY AND ASTROPHYSICS 394,
*     317-321
      SPRIME = -47.0D-6 * T
      TIOLON = -SPRIME / SECCON
*     NOTE THAT TIOLON, THE LONGITUDE CORRECTION, IS NEGLIGIBLE FOR
*     MOST ASTRONOMICAL PURPOSES

*     COMPUTE ELEMENTS OF ROTATION MATRIX
*     EQUIVALENT TO R3(-S')R2(X)R1(Y) AS PER IERS CONVENTIONS (2003)
      SINX = DSIN ( XPOLE )
      COSX = DCOS ( XPOLE )
      SINY = DSIN ( YPOLE )
      COSY = DCOS ( YPOLE )
      SINL = DSIN ( TIOLON )
      COSL = DCOS ( TIOLON )
      XX =  COSX * COSL
      YX =  SINX * SINY * COSL + COSY * SINL
      ZX = -SINX * COSY * COSL + SINY * SINL
      XY = -COSX * SINL
      YY =  SINX * SINY * SINL + COSY * COSL
      ZY =  SINX * COSY * SINL + SINY * COSL
      XZ =  SINX
      YZ = -COSX * SINY
      ZZ =  COSX * COSY
   10 IF ( TJD .LT. 0.D0 ) GO TO 30

*     PERFORM ROTATION FROM ITRS TO TERRESTRIAL INTERMEDIATE SYSTEM
   20 POS2(1) = XX * POS1(1) + YX * POS1(2) + ZX * POS1(3)
      POS2(2) = XY * POS1(1) + YY * POS1(2) + ZY * POS1(3)
      POS2(3) = XZ * POS1(1) + YZ * POS1(2) + ZZ * POS1(3)
      GO TO 50

*     PERFORM ROTATION FROM TERRESTRIAL INTERMEDIATE SYSTEM TO ITRS
   30 POS2(1) = XX * POS1(1) + XY * POS1(2) + XZ * POS1(3)
      POS2(2) = YX * POS1(1) + YY * POS1(2) + YZ * POS1(3)
      POS2(3) = ZX * POS1(1) + ZY * POS1(2) + ZZ * POS1(3)

   50 RETURN

      END



      SUBROUTINE FRAME (POS1,K,POS2)
*
*     THIS SUBROUTINE TRANSFORMS A VECTOR FROM THE DYNAMICAL REFERENCE
*     SYSTEM TO THE INTERNATIONAL CELESTIAL REFERENCE SYSTEM (ICRS),
*     OR VICE VERSA.  THE DYNAMICAL REFERENCE SYSTEM IS BASED ON THE
*     DYNAMICAL MEAN EQUATOR AND EQUINOX OF J2000.0.  THE ICRS IS
*     BASED ON THE SPACE-FIXED ICRS AXES DEFINED BY THE RADIO CATALOG
*     POSITIONS OF SEVERAL HUNDRED EXTRAGALACTIC OBJECTS.  THE ROTATION
*     MATRIX USED HERE IS EQUIVALENT TO THAT GIVEN BY HILTON AND
*     HOHENKERK (2004), ASTRONOMY AND ASTROPHYSICS 413, 765-770,
*     EQ. (6) AND (8).
*
*          POS1   = POSITION VECTOR, EQUATORIAL RECTANGULAR
*                   COORDINATES (IN)
*          K      = DIRECTION OF ROTATION (IN)
*                   SET K < 0 FOR DYNAMICAL TO ICRS
*                   SET K > 0 FOR ICRS TO DYNAMICAL
*          POS2   = POSITION VECTOR, EQUATORIAL RECTANGULAR
*                   COORDINATES (OUT)
*
*     NOTE:  FOR GEOCENTRIC COORDINATES, THE SAME TRANSFORMATION IS
*     USED BETWEEN THE DYNAMICAL REFERENCE SYSTEM AND THE GCRS.
*
*
      DOUBLE PRECISION POS1,POS2,PI,SECCON,XI0,ETA0,DA0,
     .     XX,YX,ZX,XY,YY,ZY,XZ,YZ,ZZ
      DIMENSION POS1(3), POS2(3)
      SAVE

      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( SECCON = 180.D0 * 3600.D0 / PI )

*     XI0, ETA0, AND DA0 ARE ICRS FRAME BIASES IN ARCSECONDS TAKEN
*     FROM IERS CONVENTIONS (2003), CHAPTER 5
      DATA XI0, ETA0, DA0 / -0.0166170D0, -0.0068192D0, -0.01460D0 /

      DATA NTIMES / 0 /

      NTIMES = NTIMES + 1

*     COMPUTE ELEMENTS OF ROTATION MATRIX (TO FIRST ORDER)
      IF ( NTIMES .GT. 1 ) GO TO 20
      XX =  1.D0
      YX = -DA0  / SECCON
      ZX =  XI0  / SECCON
      XY =  DA0  / SECCON
      YY =  1.D0
      ZY =  ETA0 / SECCON
      XZ = -XI0  / SECCON
      YZ = -ETA0 / SECCON
      ZZ =  1.D0
*     INCLUDE SECOND-ORDER CORRECTIONS TO DIAGONAL ELEMENTS
      XX = 1.D0 - 0.5D0 * ( YX**2 + ZX**2 )
      YY = 1.D0 - 0.5D0 * ( YX**2 + ZY**2 )
      ZZ = 1.D0 - 0.5D0 * ( ZY**2 + ZX**2 )
   20 IF ( K .GE. 0 ) GO TO 30

*     PERFORM ROTATION FROM DYNAMICAL SYSTEM TO ICRS
      POS2(1) = XX * POS1(1) + YX * POS1(2) + ZX * POS1(3)
      POS2(2) = XY * POS1(1) + YY * POS1(2) + ZY * POS1(3)
      POS2(3) = XZ * POS1(1) + YZ * POS1(2) + ZZ * POS1(3)
      GO TO 50

*     PERFORM ROTATION FROM ICRS TO DYNAMICAL SYSTEM
   30 POS2(1) = XX * POS1(1) + XY * POS1(2) + XZ * POS1(3)
      POS2(2) = YX * POS1(1) + YY * POS1(2) + YZ * POS1(3)
      POS2(3) = ZX * POS1(1) + ZY * POS1(2) + ZZ * POS1(3)

   50 RETURN

      END



      SUBROUTINE TERRA (GLON,GLAT,HT,ST,POS,VEL)
*
*     THIS SUBROUTINE COMPUTES THE POSITION AND VELOCITY VECTORS OF
*     A TERRESTRIAL OBSERVER WITH RESPECT TO THE GEOCENTER.
*
*          GLON   = LONGITUDE OF OBSERVER WITH RESPECT TO REFERENCE
*                   MERIDIAN (EAST +) IN DEGREES (IN)
*          GLAT   = GEODETIC LATITUDE (NORTH +) OF OBSERVER
*                   IN DEGREES (IN)
*          HT     = HEIGHT OF OBSERVER IN METERS (IN)
*          ST     = LOCAL APPARENT SIDEREAL TIME AT REFERENCE MERIDIAN
*                   IN HOURS (IN)
*          POS    = POSITION VECTOR OF OBSERVER WITH RESPECT TO
*                   GEOCENTER, EQUATORIAL RECTANGULAR COORDINATES,
*                   REFERRED TO TRUE EQUATOR AND EQUINOX OF DATE,
*                   COMPONENTS IN AU (OUT)
*          VEL    = VELOCITY VECTOR OF OBSERVER WITH RESPECT TO
*                   GEOCENTER, EQUATORIAL RECTANGULAR COORDINATES,
*                   REFERRED TO TRUE EQUATOR AND EQUINOX OF DATE,
*                   COMPONENTS IN AU/DAY (OUT)
*
*     NOTE 1:  IF REFERENCE MERIDIAN IS GREENWICH AND ST=0.D0, POS
*     IS EFFECTIVELY REFERRED TO EQUATOR AND GREENWICH.
*
*     NOTE 2:  THIS SUBROUTINE IGNORES POLAR MOTION, UNLESS THE
*     OBSERVER'S LONGITUDE AND LATITUDE HAVE BEEN CORRECTED FOR IT,
*     AND VARIATION IN THE LENGTH OF DAY (ANGULAR VELOCITY OF EARTH).
*     NEGLECT OF POLAR MOTION MAY YIELD 15 METERS ERROR IN POSITION
*     AND OF ORDER 1 MILLIMETER/SEC ERROR IN VELOCITY.  NEGLECT OF
*     VARIATIONS IN LENGTH OF DAY RESULTS IN EVEN SMALLER VELOCITY
*     ERRORS.
*
*     NOTE 3:  THE TRUE EQUATOR AND EQUINOX OF DATE DO NOT FORM AN
*     INERTIAL SYSTEM.  THEREFORE, WITH RESPECT TO AN INERTIAL SYSTEM,
*     THE SMALL VELOCITY COMPONENT, OF ORDER 0.1 MILLIMETER/SEC,
*     DUE TO THE PRECESSION AND NUTATION OF THE EARTH'S AXIS, IS NOT
*     ACCOUNTED FOR HERE.
*
*
      DOUBLE PRECISION GLON,GLAT,HT,ST,POS,VEL,PI,SECCON,ERAD,F,OMEGA,
     .     AUKM,DF2,PHI,SINPHI,COSPHI,C,S,ACH,ASH,STLOCL,SINST,COSST,
     .     DSQRT,DCOS,DSIN
      DIMENSION POS(3), VEL(3)
      SAVE

      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( SECCON = 180.D0 * 3600.D0 / PI )

      DATA NTIMES / 0 /

      NTIMES = NTIMES + 1
      IF (NTIMES.EQ.1) THEN
*         GET ERAD, THE EQUATORIAL RADIUS OF EARTH IN KILOMETERS
          CALL ASTCON ('ERAD',1.D-3,ERAD)
*         GET F, THE FLATTENING FACTOR OF EARTH ELLIPSOID
          CALL ASTCON ('F',1.D0,F)
*         GET OMEGA, THE NOMINAL MEAN ROTATIONAL ANGULAR VELOCITY OF
*         EARTH IN RADIANS/SECOND
          CALL ASTCON ('ANGVEL',1.D0,OMEGA)
*         GET AUKM, THE LENGTH OF THE ASTRONOMICAL UNIT IN KILOMETERS
          CALL ASTCON ('AU',1.D-3,AUKM)
      END IF

*     COMPUTE PARAMETERS RELATING TO GEODETIC TO GEOCENTRIC CONVERSION
      DF2 = (1.D0 - F)**2
      PHI = GLAT * 3600.D0 / SECCON
      SINPHI = DSIN(PHI)
      COSPHI = DCOS(PHI)
      C = 1.D0 / DSQRT ( COSPHI**2 + DF2 * SINPHI**2 )
      S = DF2 * C
      ACH = ERAD * C + HT/1000.D0
      ASH = ERAD * S + HT/1000.D0

*     COMPUTE LOCAL SIDEREAL TIME FACTORS
      STLOCL = (ST * 54000.D0 + GLON * 3600.D0) / SECCON
      SINST = DSIN(STLOCL)
      COSST = DCOS(STLOCL)

*     COMPUTE POSITION VECTOR COMPONENTS IN KM
      POS(1) = ACH * COSPHI * COSST
      POS(2) = ACH * COSPHI * SINST
      POS(3) = ASH * SINPHI

*     COMPUTE VELOCITY VECTOR COMPONENTS IN KM/SEC
      VEL(1) = -OMEGA * ACH * COSPHI * SINST
      VEL(2) =  OMEGA * ACH * COSPHI * COSST
      VEL(3) =  0.D0

*     CONVERT POSITION AND VELOCITY COMPONENTS TO AU AND AU/DAY
      DO 20 J=1,3
      POS(J) = POS(J) / AUKM
      VEL(J) = VEL(J) / AUKM * 86400.D0
   20 CONTINUE

      RETURN

      END
      
      
      
      SUBROUTINE TIMES (TDBJD,TTJD,SECDIF)
*
*     THIS SUBROUTINE COMPUTES THE TERRESTRIAL TIME (TT) JULIAN DATE
*     CORRESPONDING TO A BARYCENTRIC DYNAMICAL TIME (TDB) JULIAN DATE.
*     THE EXPRESSION USED IN THIS VERSION IS A TRUNCATED FORM OF A
*     LONGER AND MORE PRECISE SERIES GIVEN BY FAIRHEAD & BRETAGNON
*     (1990) A&A 229, 240.  THE RESULT IS GOOD TO ABOUT 10 MICROSECONDS.
*
*          TDBJD  = TDB JULIAN DATE (IN)
*          TTJD   = TT JULIAN DATE (OUT)
*          SECDIF = DIFFERENCE TDBJD-TTJD, IN SECONDS (OUT)
*
*
      DOUBLE PRECISION TDBJD,TTJD,SECDIF,T,T0,DSIN

*     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
      DATA T0 / 2451545.00000000D0 /

      T = ( TDBJD - T0 ) / 36525.D0

*     EXPRESSION GIVEN IN USNO CIRCULAR 179, EQ. 2.6 
      SECDIF = 0.001657D0 * DSIN (  628.3076D0 * T + 6.2401D0)                     
     .       + 0.000022D0 * DSIN (  575.3385D0 * T + 4.2970D0)   
     .       + 0.000014D0 * DSIN ( 1256.6152D0 * T + 6.1969D0)            
     .       + 0.000005D0 * DSIN (  606.9777D0 * T + 4.0212D0)                                     
     .       + 0.000005D0 * DSIN (   52.9691D0 * T + 0.4444D0)     
     .       + 0.000002D0 * DSIN (   21.3299D0 * T + 5.5431D0) 
     .       + 0.000010D0 * T * DSIN ( 628.3076D0 * T + 4.2490D0)

      TTJD = TDBJD - SECDIF / 86400.D0

      RETURN

      END



      SUBROUTINE ETILT (TJD,OBLM,OBLT,EQEQ,DPSI,DEPS)
*
*     THIS SUBROUTINE COMPUTES QUANTITIES RELATED TO THE ORIENTATION
*     OF THE EARTH'S ROTATION AXIS AT JULIAN DATE TJD.
*
*          TJD    = TDB JULIAN DATE FOR ORIENTATION PARAMETERS (IN)
*          OBLM   = MEAN OBLIQUITY OF THE ECLIPTIC IN DEGREES AT
*                   DATE TJD (OUT)
*          OBLT   = TRUE OBLIQUITY OF THE ECLIPTIC IN DEGREES AT
*                   DATE TJD (OUT)
*          EQEQ   = EQUATION OF THE EQUINOXES IN TIME SECONDS AT
*                   DATE TJD (OUT)
*          DPSI   = NUTATION IN LONGITUDE IN ARCSECONDS AT
*                   DATE TJD (OUT)
*          DEPS   = NUTATION IN OBLIQUITY IN ARCSECONDS AT
*                   DATE TJD (OUT)
*
*     NOTE:  THE EQUATION OF THE EQUINOXES INCLUDES THE COMPLEMENTARY
*     TERMS.
*
*
      DOUBLE PRECISION TJD,OBLM,OBLT,EQEQ,DPSI,DEPS,PI,SECCON,
     .     T0,TLAST,T,PSI,EPS,PSICOR,EPSCOR,CTERMS,DELPSI,DELEPS,
     .     EL,ELP,F,D,OMEGA,OBM,OBT,EE,
     .     DPOLE1,DPOLE2,DX,DY,DZ,SINE,X,DP1,DP2,DP3,
     .     OBLIQ,DABS,DSIN,DCOS,EECT2000
      INTEGER ACCDIF
      DIMENSION DP1(3), DP2(3), DP3(3)
      SAVE

      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( SECCON = 180.D0 * 3600.D0 / PI )

*     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
      DATA T0 / 2451545.00000000D0 /
      DATA TLAST / 0.D0 /,   MLAST / 0 /
      DATA DELPSI, DELEPS, CTERMS, PSICOR, EPSCOR / 5 * 0.D0 /

*     FUNCTION TO COMPUTE MEAN OBLIQUITY OF THE ECLIPTIC IN ARCSECONDS
*     CAPITAINE ET AL. (2003), ASTRONOMY AND ASTROPHYSICS 412, 567-586,
*     EXPRESSION FROM EQ. (39) WITH OBLIQUITY AT J2000.0 TAKEN FROM
*     EQ. (37) OR TABLE 8
      OBLIQ(T) = ( ( ( ( -  0.0000000434D0   * T
     .                   -  0.000000576D0  ) * T
     .                   +  0.00200340D0   ) * T
     .                   -  0.0001831D0    ) * T
     .                   - 46.836769D0     ) * T + 84381.406D0

*     GET METHOD/ACCURACY MODE
      CALL GETMOD ( MODE )

*     CHECK FOR DIFFERENCE IN ACCURACY MODE FROM LAST CALL
      ACCDIF = MOD ( MODE, 2 ) - MOD ( MLAST, 2 )

      T = ( TJD - T0 ) / 36525.D0

      IF ( DABS ( TJD - TLAST ) .GT. 1.D-8 .OR. ACCDIF .NE. 0 ) THEN

*         OBTAIN NUTATION PARAMETERS IN ARCSECONDS
          CALL NOD ( T,   PSI, EPS )

*         OBTAIN COMPLEMENTARY TERMS FOR EQUATION OF THE EQUINOXES
*         IN ARCSECONDS
          IF ( MOD ( MODE, 2 ) .EQ. 0 ) THEN
*             HIGH-ACCURACY MODE
              CTERMS = EECT2000 ( TJD, 0.D0 ) * SECCON
          ELSE
*             LOW-ACCURACY MODE
              CALL FUNARG ( T,   EL, ELP, F, D, OMEGA )
*             SERIES FROM IERS CONVENTIONS (2003), CHAPTER 5,
*             TABLE 5.2C, WITH SOME ADJUSTMENTS TO COEFFICIENT VALUES
*             COPIED FROM IERS FUNCTION EECT2000, WHICH HAS A MORE
*             COMPLETE SERIES
              CTERMS =
     .          2640.96D-6 * DSIN ( OMEGA )
     .        +   63.52D-6 * DSIN ( 2.D0 * OMEGA )
     .        +   11.75D-6 * DSIN ( 2.D0 * F - 2.D0 * D + 3.D0 * OMEGA )
     .        +   11.21D-6 * DSIN ( 2.D0 * F - 2.D0 * D +        OMEGA )
     .        -    4.55D-6 * DSIN ( 2.D0 * F - 2.D0 * D + 2.D0 * OMEGA )
     .        +    2.02D-6 * DSIN ( 2.D0 * F            + 3.D0 * OMEGA )
     .        +    1.98D-6 * DSIN ( 2.D0 * F            +        OMEGA )
     .        -    1.72D-6 * DSIN ( 3.D0 * OMEGA )
     .        -    0.87D-6 * T * DSIN ( OMEGA )
*             (TERMS SMALLER THAN 2 MICROARCSECONDS OMITTED)
          END IF
          TLAST = TJD
          MLAST = MODE

      END IF

      DELPSI = PSI + PSICOR
      DELEPS = EPS + EPSCOR

*     COMPUTE MEAN OBLIQUITY OF THE ECLIPTIC IN ARCSECONDS
      OBM = OBLIQ(T)

*     COMPUTE TRUE OBLIQUITY OF THE ECLIPTIC IN ARCSECONDS
      OBT = OBM + DELEPS

*     COMPUTE EQUATION OF THE EQUINOXES IN ARCSECONDS
      EE = DELPSI * DCOS ( OBM / SECCON ) + CTERMS

*     CONVERT TO OUTPUT UNITS
      OBLM = OBM / 3600.D0
      OBLT = OBT / 3600.D0
      EQEQ = EE  / 15.D0
      DPSI = DELPSI
      DEPS = DELEPS

      RETURN


      ENTRY CELPOL (TJD,ITYPE,DPOLE1,DPOLE2)
*
*     THIS ENTRY ALLOWS FOR THE SPECIFICATION OF CELESTIAL POLE
*     OFFSETS FOR HIGH-PRECISION APPLICATIONS.  EACH SET OF OFFSETS IS
*     A CORRECTION TO THE MODELED POSITION OF THE POLE FOR A SPECIFIC
*     DATE, DERIVED FROM OBSERVATIONS AND PUBLISHED BY THE IERS.
*     THIS ENTRY, IF USED, SHOULD BE CALLED BEFORE ANY OTHER ROUTINES
*     FOR A GIVEN DATE.  VALUES OF THE POLE OFFSETS SPECIFIED VIA A CALL
*     TO THIS ENTRY WILL BE USED UNTIL EXPLICITLY CHANGED.
*
*          TJD    = TDB OR TT JULIAN DATE FOR POLE OFFSETS (IN)
*          ITYPE  = TYPE OF POLE OFFSET (IN)
*                   SET ITYPE=1 FOR CORRECTIONS TO ANGULAR COORDINATES
*                               OF MODELED POLE REFERRED TO MEAN
*                               ECLIPTIC OF DATE, THAT IS,
*                               DELTA-DELTA-PSI AND DELTA-DELTA-EPSILON
*                   SET ITYPE=2 FOR CORRECTIONS TO COMPONENTS OF
*                               MODELED POLE UNIT VECTOR WITH REFERRED
*                               TO GCRS AXES, THAT IS, DX AND DY
*          DPOLE1 = VALUE OF CELESTIAL POLE OFFSET IN FIRST COORDINATE,
*                   (DELTA-DELTA-PSI OR DX) IN MILLIARCSECONDS (IN)
*          DPOLE2 = VALUE OF CELESTIAL POLE OFFSET IN SECOND COORDINATE,
*                   (DELTA-DELTA-EPSILON OR DY) IN MILLIARCSECONDS (IN)
*
*     NOTE 1:  TJD IS USED ONLY FOR ITYPE=2, TO TRANSFORM DX AND DY TO
*     THE EQUIVALENT DELTA-DELTA-PSI AND DELTA-DELTA-EPSILON VALUES.
*
*     NOTE 2:  FOR ITYPE=2, DX AND DY ARE UNIT VECTOR COMPONENT
*     CORRECTIONS, BUT ARE EXPRESSED IN MILLIARCSECONDS SIMPLY BY
*     MULTIPLYING BY 206264806, THE NUMBER OF MILLIARCSECONDS IN ONE
*     RADIAN.
*
*
      IF ( ITYPE .EQ. 1 ) THEN

          PSICOR = DPOLE1 * 1.D-3
          EPSCOR = DPOLE2 * 1.D-3

      ELSE

          DX = DPOLE1
          DY = DPOLE2

          T = ( TJD - T0 ) / 36525.D0
*         COMPUTE SINE OF MEAN OBLIQUITY OF DATE
          SINE = DSIN ( OBLIQ(T) / SECCON )

*         THE FOLLOWING ALGORITHM, TO TRANSFORM DX AND DY TO
*         DELTA-DELTA-PSI AND DELTA-DELTA-EPSILON, IS FROM G. KAPLAN
*         (2003), USNO/AA TECHNICAL NOTE 2003-03, EQS. (7)-(9).

*         TRIVIAL MODEL OF POLE TRAJECTORY IN GCRS ALLOWS COMPUTATION
*         OF DZ
          X = ( 2004.19D0 * T ) / SECCON
          DZ = - ( X + 0.5D0 * X**3 ) * DX

*         FORM POLE OFFSET VECTOR (OBSERVED - MODELED) IN GCRS
          DP1(1) = DX * 1.D-3 / SECCON
          DP1(2) = DY * 1.D-3 / SECCON
          DP1(3) = DZ * 1.D-3 / SECCON

*         PRECESS POLE OFFSET VECTOR TO MEAN EQUATOR AND EQUINOX OF DATE
          CALL FRAME ( DP1, 1, DP2 )
          CALL PRECES ( T0, DP2, TJD, DP3 )

*         COMPUTE DELTA-DELTA-PSI AND DELTA-DELTA-EPSILON IN ARCSECONDS
          PSICOR = ( DP3(1) / SINE ) * SECCON
          EPSCOR = ( DP3(2)        ) * SECCON

      END IF

      RETURN

      END



      SUBROUTINE FUNARG (T,EL,ELP,F,D,OMEGA)
*
*     THIS SUBROUTINE COMPUTES FUNDAMENTAL ARGUMENTS (MEAN ELEMENTS)
*     OF THE SUN AND MOON.  SEE SIMON ET AL. (1994) ASTRONOMY AND
*     ASTROPHYSICS 282, 663-683, ESPECIALLY SECTIONS 3.4-3.5.
*
*          T      = TDB TIME IN JULIAN CENTURIES SINCE J2000.0 (IN)
*          EL     = MEAN ANOMALY OF THE MOON IN RADIANS
*                   AT DATE TJD (OUT)
*          ELP    = MEAN ANOMALY OF THE SUN IN RADIANS
*                   AT DATE TJD (OUT)
*          F      = MEAN LONGITUDE OF THE MOON MINUS MEAN LONGITUDE
*                   OF THE MOON'S ASCENDING NODE IN RADIANS
*                   AT DATE TJD (OUT)
*          D      = MEAN ELONGATION OF THE MOON FROM THE SUN IN
*                   RADIANS AT DATE TJD (OUT)
*          OMEGA  = MEAN LONGITUDE OF THE MOON'S ASCENDING NODE
*                   IN RADIANS AT DATE TJD (OUT)
*
*
      DOUBLE PRECISION T,EL,ELP,F,D,OMEGA,PI,SECCON,REV,DMOD

      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( SECCON = 180.D0 * 3600.D0 / PI )
      PARAMETER ( REV    = 360.D0 * 3600.D0      )

*     FUNDAMENTAL (DELAUNAY) ARGUMENTS FROM SIMON ET AL. (1994)

*     MEAN ANOMALY OF THE MOON
      EL    = DMOD (         485868.249036D0 +
     .               T*( 1717915923.2178D0 +
     .               T*(         31.8792D0 +
     .               T*(          0.051635D0 +
     .               T*(        - 0.00024470D0 )))), REV ) / SECCON

*     MEAN ANOMALY OF THE SUN
      ELP   = DMOD (        1287104.79305D0 +
     .               T*(  129596581.0481D0 +
     .               T*(        - 0.5532D0 +
     .               T*(          0.000136D0 +
     .               T*(        - 0.00001149D0 )))), REV ) / SECCON

*     MEAN ARGUMENT OF THE LATITUDE OF THE MOON
      F     = DMOD (         335779.526232D0 +
     .               T*( 1739527262.8478D0 +
     .               T*(       - 12.7512D0 +
     .               T*(       -  0.001037D0 +
     .               T*(          0.00000417D0 )))), REV ) / SECCON

*     MEAN ELONGATION OF THE MOON FROM THE SUN
      D     = DMOD (        1072260.70369D0 +
     .               T*( 1602961601.2090D0 +
     .               T*(        - 6.3706D0 +
     .               T*(          0.006593D0 +
     .               T*(        - 0.00003169D0 )))), REV ) / SECCON

*     MEAN LONGITUDE OF THE ASCENDING NODE OF THE MOON (FROM SIMON
*     SECTION 3.4(b.3), PRECESSION=5028.8200 ARCSEC/CY)
      OMEGA = DMOD (         450160.398036D0 +
     .               T*(  - 6962890.5431D0 +
     .               T*(          7.4722D0 +
     .               T*(          0.007702D0 +
     .               T*(        - 0.00005939D0 )))), REV ) / SECCON

      RETURN

      END



      SUBROUTINE REFRAC (HEIGHT,ZDOBS,REFR)
*
*     THIS SUBROUTINE COMPUTES ATMOSPHERIC REFRACTION IN ZENITH
*     DISTANCE.  THIS VERSION COMPUTES APPROXIMATE REFRACTION FOR
*     OPTICAL WAVELENGTHS.  IT CAN BE USED FOR PLANNING OBSERVATIONS
*     OR TELESCOPE POINTING, BUT SHOULD NOT BE USED FOR THE REDUCTION
*     OF PRECISE OBSERVATIONS.  BASIC ALGORITHM IS DESCRIBED IN THE
*     EXPLANATORY SUPPLEMENT TO THE ASTRONOMICAL ALMANAC, P. 144,
*     AND IS AN ADAPTATION OF A FORMULA IN BENNETT (1982), JOURNAL
*     OF NAVIGATION (ROYAL INSTITUTE) 35, 255-259.
*
*          HEIGHT = HEIGHT OF OBSERVER IN METERS (IN)
*          ZDOBS  = OBSERVED ZENITH DISTANCE IN DEGREES (IN)
*          REFR   = ATMOSPHERIC REFRACTION IN DEGREES (OUT)
*
*     NOTE:  HEIGHT IS NOT USED IF ENTRY REFDAT HAS BEEN CALLED
*     TO SPECIFY ATMOSPHERIC PRESSURE.
*
*
      DOUBLE PRECISION HEIGHT,ZDOBS,REFR,PI,DEGRAD,S,
     .     POBS,TOBS,DOBS,WLOBS,OBSP,OBST,OBSD,OBSWL,P,T,D,WL,H,R,
     .     DEXP,DTAN
      SAVE

      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( DEGRAD = PI / 180.D0           )

*     S IS APPROXIMATE SCALE HEIGHT OF ATMOSPHERE IN METERS
      DATA S / 9.1D3 /
      DATA POBS, TOBS, DOBS, WLOBS / 4 * -999.D0 /

*     COMPUTE REFRACTION ONLY FOR ZENITH DISTANCES
*     BETWEEN 0.1 AND 91 DEGREES
      IF ( ZDOBS .LT. 0.1D0 .OR. ZDOBS .GT. 91.D0 ) THEN
          REFR = 0.D0
          GO TO 77
      END IF

*     IF OBSERVED WEATHER DATA ARE AVAILABLE, USE THEM
*     OTHERWISE, USE CRUDE ESTIMATES OF AVERAGE CONDITIONS
      IF ( POBS .GE. 1.D0 .AND. TOBS .GT. -100.D0 ) THEN
          P  = POBS
          T  = TOBS
          D  = DOBS
          WL = WLOBS
      ELSE
          P  = 1010.D0 * DEXP ( -HEIGHT / S )
          T  = 10.D0
          D  =  0.D0
          WL =  0.5D0
      END IF
*     D AND WL NOT USED IN THIS VERSION

      H = 90.D0 - ZDOBS
      R = 0.016667D0 / DTAN ( ( H +  7.31D0 / ( H + 4.4D0 ) ) * DEGRAD )
      REFR = R * ( 0.28D0 * P / ( T + 273.D0 ) )

  77  RETURN


      ENTRY REFDAT (OBSP,OBST,OBSD,OBSWL)
*
*     THIS ENTRY ALLOWS FOR THE SPECIFICATION OF WEATHER OBSERVATIONS
*     AND OTHER DATA TO BE USED IN THE ATMOSPHERIC REFRACTION
*     CALCULATION.  THIS ENTRY, IF USED, SHOULD BE CALLED BEFORE
*     SUBROUTINE REFRAC OR ZDAZ FOR A GIVEN DATE/TIME.  DATA SPECIFIED
*     VIA A CALL TO THIS ENTRY WILL BE USED UNTIL EXPLICITLY CHANGED.
*
*          OBSP   = OBSERVED ATMOSPHERIC PRESSURE IN MILLIBARS (IN)
*          OBST   = OBSERVED TEMPERATURE IN DEGREES CELSIUS (IN)
*          OBSD   = OBSERVED DEW POINT IN DEGREES CELSIUS (IN)
*          OBSWL  = OBSERVING WAVELENGTH IN MICRONS (IN)
*
*     NOTE:  OBSD AND OBSWL ARE NOT USED IN THIS VERSION'S REFRACTION
*     ALGORITHM, AND CAN BE SET TO ANY VALUE.
*
*
      POBS  = OBSP
      TOBS  = OBST
      DOBS  = OBSD
      WLOBS = OBSWL
      RETURN

      END



      SUBROUTINE LIMANG (POS1,POSO,ALIMB,AFRAC)
*
*     THIS FUNCTION DETERMINES THE ANGLE OF AN OBJECT ABOVE OR BELOW
*     THE EARTH'S LIMB (HORIZON).  THE GEOMETRIC LIMB IS COMPUTED,
*     ASSUMING THE EARTH TO BE AN AIRLESS SPHERE (NO REFRACTION OR
*     OBLATENESS IS INCLUDED).  THE OBSERVER CAN BE ON OR ABOVE THE
*     EARTH.  FOR AN OBSERVER ON THE SURFACE OF THE EARTH, THIS
*     SUBROUTINE RETURNS THE APPROXIMATE UNREFRACTED ALTITUDE.
*
*          POS1   = POSITION VECTOR OF OBSERVED OBJECT, WITH RESPECT TO
*                   ORIGIN AT GEOCENTER, COMPONENTS IN AU (IN)
*          POSO   = POSITION VECTOR OF OBSERVER, WITH RESPECT TO ORIGIN
*                   AT GEOCENTER, COMPONENTS IN AU (IN)
*          ALIMB  = ANGLE OF OBSERVED OBJECT ABOVE (+) OR BELOW (-) LIMB
*                   IN DEGREES (OUT)
*          AFRAC  = NADIR ANGLE OF OBSERVED OBJECT AS A FRACTION OF
*                   APPARENT RADIUS OF LIMB (OUT)
*                   AFRAC<1.D0 MEANS BELOW THE LIMB
*                   AFRAC=1.D0 MEANS ON THE LIMB
*                   AFRAC>1.D0 MEANS ABOVE THE LIMB
*
*
      DOUBLE PRECISION POS1,POSO,ALIMB,AFRAC,PI,HALFPI,DEGCON,
     .     ERAD,AU,RADE,DISOBJ,DISOBS,APRAD,ZDLIM,COSZD,ZDOBJ,
     .     DSQRT,DACOS
      DIMENSION POS1(3), POSO(3)
      SAVE

      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( HALFPI = 0.5D0 * PI            )
      PARAMETER ( DEGCON = 180.D0 / PI           )

      DATA NTIMES / 0 /

      NTIMES = NTIMES + 1
      IF ( NTIMES .EQ. 1 ) THEN
           CALL ASTCON ( 'ERAD', 1.D0,   ERAD )
           CALL ASTCON ( 'AU', 1.D0,   AU )
           RADE = ERAD / AU
      END IF

      DISOBJ = DSQRT ( POS1(1)**2 + POS1(2)**2 + POS1(3)**2 )
      DISOBS = DSQRT ( POSO(1)**2 + POSO(2)**2 + POSO(3)**2 )

*     COMPUTE APPARENT ANGULAR RADIUS OF EARTH'S LIMB
      IF ( DISOBS .GT. RADE ) THEN
          APRAD = DASIN ( RADE / DISOBS )
      ELSE
          APRAD = HALFPI
      END IF

*     COMPUTE ZENITH DISTANCE OF EARTH'S LIMB
      ZDLIM = PI - APRAD

*     COMPUTE ZENITH DISTANCE OF OBSERVED OBJECT
      COSZD = ( POS1(1)*POSO(1) + POS1(2)*POSO(2) + POS1(3)*POSO(3) )
     .      / ( DISOBJ * DISOBS )
      IF ( COSZD .LE. -1.D0 ) THEN
         ZDOBJ = PI
      ELSE IF ( COSZD .GE. 1.D0 ) THEN
         ZDOBJ = 0.D0
      ELSE
         ZDOBJ = DACOS ( COSZD )
      END IF

*     ANGLE OF OBJECT WRT LIMB IS DIFFERENCE IN ZENITH DISTANCES
      ALIMB = ( ZDLIM - ZDOBJ ) * DEGCON

*     NADIR ANGLE OF OBJECT AS A FRACTION OF ANGULAR RADIUS OF LIMB
      AFRAC = ( PI - ZDOBJ ) / APRAD

      RETURN

      END
      
      
      
      SUBROUTINE CIOLOC ( TJD,   RACIO, K )
*
*     THIS SUBROUTINE RETURNS THE LOCATION OF THE CELESTIAL
*     INTERMEDIATE ORIGIN (CIO) FOR A GIVEN JULIAN DATE, AS A
*     RIGHT ASCENSION WITH RESPECT TO EITHER THE GCRS (GEOCENTRIC ICRS)
*     ORIGIN OR THE TRUE EQUINOX OF DATE.  THE CIO IS ALWAYS LOCATED ON
*     THE TRUE EQUATOR (=INTERMEDIATE EQUATOR) OF DATE.
*
*          TJD    = TDB JULIAN DATE (IN)
*          RACIO  = RIGHT ASCENSION OF THE CIO, IN HOURS (OUT)
*          K      = REFERENCE SYSTEM IN WHICH RIGHT ASCENSION IS
*                   GIVEN (OUT)
*                   K=1 MEANS GCRS
*                   K=2 MEANS TRUE EQUATOR AND EQUINOX OF DATE
*
*     NOTE:  IF AN EXTERNAL FILE OF CIO RIGHT ASCENSIONS IS AVAILABLE,
*     IT WILL BE USED AND K WILL BE SET TO 1.  OTHERWISE AN INTERNAL
*     COMPUTATION WILL BE USED AND K WILL BE SET TO 2.
*
*
      DOUBLE PRECISION TJD,RACIO,A,TLAST,RLAST,JD,RA,P,EQOR,DABS
      LOGICAL USEFIL
      DIMENSION JD(8), RA(8), A(1)
      SAVE

*     NUMBER OF VALUES IN ARRAYS FOR LAGRANGIAN INTERPOLATION 
      DATA M / 6 /
      
      DATA TLAST, RLAST, KLAST / 0.D0, 0.D0, 0 /
      
   3  FORMAT ( ' CIOLOC ERROR: CANNOT RETURN CIO RA VALUE FOR JD ',
     .     F10.1 )

*     CHECK IF EXTERNAL FILE EXISTS
      CALL CIORD ( 0.D0, 1,   A, A, IERR )
      USEFIL = IERR .EQ. 0

*     CHECK IF PREVIOUSLY COMPUTED RA VALUE CAN BE USED
      IF ( DABS ( TJD - TLAST ) .LE. 1.D-8 ) THEN
          RACIO = RLAST
          K = KLAST
          GO TO 77
      END IF

* --- IF EXTERNAL FILE EXISTS, INTERPOLATE RA VALUES FROM FILE --------

      IF ( USEFIL ) THEN
      
          K = 1

*         GET ARRAYS OF VALUES TO INTERPOLATE
          CALL CIORD ( TJD, M,   JD, RA, IERR )
          IF ( IERR .NE. 0 ) THEN
              WRITE ( *, 3 ) TJD
              RACIO = 99.D0
              GO TO 77
          END IF

*         PERFORM LAGRANGIAN INTERPOLATION FOR RA AT TJD
          RACIO = 0.D0
          DO 40 J = 1, M
              P = 1.D0
              DO 30 I = 1, M
                  IF ( I .EQ. J ) GO TO 30
                  P = P * ( TJD - JD(I) ) / ( JD(J) - JD(I) )
  30          CONTINUE
              RACIO = RACIO + P * RA(J)
  40      CONTINUE

          RACIO = RACIO / 54000.D0
 
* --- OTHERWISE, USE INTERNAL COMPUTATION ----------------------------

      ELSE

          K = 2

*         GET EQUATION OF THE ORIGINS
          CALL EQXRA ( TJD, 1,   EQOR )
          
          RACIO = -EQOR
 
      END IF

* ---------------------------------------------------------------------

      TLAST = TJD
      RLAST = RACIO
      KLAST = K

  77  RETURN

      END



      SUBROUTINE CIORD ( TJD, NVALS,   TLIST, RALIST, IERR )
*
*     GIVEN AN INPUT TDB JULIAN DATE AND THE NUMBER OF DATA POINTS
*     DESIRED, THIS SUBROUTINE RETURNS A SET OF JULIAN DATES AND
*     CORRESPONDING VALUES OF THE GCRS RIGHT ASCENSION OF THE CELESTIAL
*     INTERMEDIATE ORIGIN (CIO).  THE RANGE OF DATES IS CENTERED (AT LEAST
*     APPROXIMATELY) ON THE REQUESTED DATE.  THE SUBROUTINE OBTAINS THE
*     DATA FROM AN EXTERNAL DATA FILE.
*
*         TJD    = TDB JULIAN DATE (IN)
*         NVALS  = NUMBER OF JULIAN DATES AND RIGHT ASCENSION VALUES
*                  REQUESTED (NOT LESS THAN 2 OR MORE THAN 20) (IN)
*         TLIST  = ARRAY OF TDB JULIAN DATES (OUT)
*         RALIST = ARRAY OF GCRS RIGHT ASCENSIONS OF THE CIO, FOR THE
*                  JULIAN DATES IN TLIST, IN ARCSECONDS (OUT)
*         IERR   = ERROR INDICATOR (OUT)
*                  IERR=0 MEANS EVERYTHING OK
*                  IERR=1 MEANS TJD BEFORE FIRST USABLE DATE IN FILE
*                  IERR=2 MEANS TJD AFTER LAST USABLE DATE IN FILE
*                  IERR=3 MEANS BAD VALUE OF NVALS
*                  IERR=4 MEANS EXTERNAL FILE CANNOT BE FOUND
*
*     NOTE:  TJD=0.D0 WITH NVALS=1 INDICATES A SPECIAL CALL JUST TO
*     DETERMINE IF EXTERNAL FILE EXISTS.
*
*
      DOUBLE PRECISION TJD,TLIST,RALIST,T,T1,R,TBEG,TEND,TINT,DIF
      CHARACTER FILNAM*40, FILEID*(*)
      LOGICAL FILEOK
      DIMENSION TLIST(NVALS), RALIST(NVALS), T(20), R(20)
      SAVE

*     LOGICAL UNIT NUMBER AND FILE ID OF TIME SERIES OF CIO RA VALUES
      DATA LU, ITYP / 24, 1 /
      DATA FILNAM / 'CIO_RA.TXT                              ' /

      DATA NTIMES, TBEG, TEND, FILEOK / 0, 0.D0, 1.D10, .FALSE. /

   1  FORMAT ( A )
   2  FORMAT ( ' CIORD ERROR: CANNOT FIND FILE ', A )
   3  FORMAT ( ' CIORD ERROR: REQUESTED JD ', F10.1, 1X, A,
     .      ' ALLOWED JD ', F10.1 )
   4  FORMAT ( F16.6, F24.14 )

*     SPECIAL CALL JUST TO DETERMINE IF FILE EXITS
*     (NO PRINTED ERROR MESSAGE IF NOT)
      IF ( TJD .EQ. 0.D0 .AND. NVALS .EQ. 1 ) THEN
          IERR = 4
          IF ( NTIMES .EQ. 0 ) INQUIRE ( FILE=FILNAM, EXIST=FILEOK )
          IF ( FILEOK ) IERR = 0
          GO TO 77
      END IF

*     IF EXTERNAL FILE IS ALREADY KNOWN NOT TO EXIST, IMMEDIATELY
*     EXIT WITH IERR=4 
      IF ( NTIMES .GT. 0 .AND. .NOT. FILEOK ) THEN  
          WRITE ( *, 2 ) FILNAM
          IERR = 4
          GO TO 77
      END IF 

*     CHECK FOR VALID VALUE OF NVALS
      IF ( NVALS .LT. 2 .OR. NVALS .GT. 20 ) THEN
          IERR = 3
          GO TO 77
      END IF

      MIDDL = NVALS / 2
   
*     CHECK THAT REQUESTED JULIAN DATE IS WITHIN RANGE OF FILE
  10  IF ( TJD .LT. TBEG ) THEN
          WRITE ( *, 3 )  TJD, 'BEFORE FIRST', TBEG
          IERR = 1
          GO TO 77
      ELSE IF ( TJD .GT. TEND ) THEN
          WRITE ( *, 3 ) TJD, 'AFTER LAST', TEND
          IERR = 2
          GO TO 77
      END IF

      IF ( ITYP .EQ. 1 ) THEN

*         -------------------------------------------------------------
*         READ JULIAN DATES AND CIO RA VALUES FROM FORMATTED
*         SEQUENTIAL INPUT FILE

*         EACH RECORD OF THE FILE MUST CONTAIN A TDB JULIAN DATE
*         AND A CORRESPONDING CIO RA VALUE (WRT GCRS) IN ARCSECONDS

*         THE JULIAN DATES MUST INCREASE BY A FIXED INTERVAL
*         -------------------------------------------------------------

*         IF FIRST TIME, OPEN FILE AND READ INITIAL VALUES
          NTIMES = NTIMES + 1
          IF ( NTIMES .EQ. 1 ) THEN
              INQUIRE ( FILE=FILNAM, EXIST=FILEOK )
              IF ( .NOT. FILEOK ) THEN
                  WRITE ( *, 2 ) FILNAM
                  IERR = 4
                  GO TO 77
              END IF
              OPEN ( UNIT=LU, FILE=FILNAM, FORM='FORMATTED',
     .               ACCESS='SEQUENTIAL', STATUS='OLD' )
              READ ( UNIT=LU, FMT=1 )
              DO 19 I = 1, NVALS
                  READ ( UNIT=LU, FMT=4, END=50 ) T(I), R(I)
  19          CONTINUE
              TINT = NINT ( ( T(2) - T(1) ) * 1000.D0 ) / 1000.D0
              TBEG = T(MIDDL)
              IF ( TJD .LT. TBEG ) GO TO 10
          END IF

*         -------------------------------------------------------------

*         FILE READ SEQUENCE

  20      DIF = ( TJD - T(MIDDL) ) / TINT
          NDIF = DIF

*         BASIC DECISION ON HOW TO READ FILE
          IF ( DIF .GE. -0.1D0 .AND. DIF .LE. 1.1D0 ) THEN
*             NO FILE READ NECESSARY, DATA PREVIOUSLY READ CAN BE USED
              GO TO 70
          ELSE IF ( DIF .LT. 0.D0 ) THEN
*             REQUESTED JULIAN DATE BEFORE DATA PREVIOUSLY READ
              IREC = ( T(NVALS) - TBEG ) / TINT
              NBACK = 3 * NVALS
              IF ( -DIF .LE. 2 * NVALS .AND. IREC .GT. NBACK ) GO TO 34
              GO TO 25
          ELSE IF ( NDIF .GT. ( NVALS + 1 ) ) THEN
*             REQUESTED JULIAN DATE FAR AHEAD OF DATA PREVIOUSLY READ
              NSKIP = NDIF - NVALS - 1
              GO TO 30
          ELSE
*             REQUESTED JULIAN DATE A BIT AHEAD OF DATA PREVIOUSLY READ
              GO TO 40
          END IF

*         REPOSITION FILE AT BEGINNING
  25      REWIND ( UNIT=LU )
          READ ( UNIT=LU, FMT=1 )
          GO TO 36

*         FAST SKIP FORWARD
  30      DO 32 I = 1, NSKIP
              READ ( UNIT=LU, FMT=1, END=50 )
  32      CONTINUE
          GO TO 36

*         BACKSPACE FILE
  34      DO 35 I = 1, NBACK
              BACKSPACE ( UNIT=LU )
  35      CONTINUE

*         FILL UP ARRAYS
  36      DO 38 I = 1, NVALS
              READ ( UNIT=LU, FMT=4, END=50 ) T(I), R(I)
  38      CONTINUE
          GO TO 20

*         ADVANCE ARRAY DATA AND READ ONE MORE RECORD
  40      DO 44 I = 1, NVALS - 1
              T(I) = T(I+1)
              R(I) = R(I+1)
  44      CONTINUE
          READ ( UNIT=LU, FMT=4, END=50 ) T(NVALS), R(NVALS)

          GO TO 20

*         -------------------------------------------------------------

*         END OF FILE ENCOUNTERED
  50      BACKSPACE ( UNIT=LU )
          BACKSPACE ( UNIT=LU )
          READ ( UNIT=LU, FMT=4 ) TEND
          TEND = TEND - ( NVALS - MIDDL - 1 ) * TINT
          WRITE ( *, 3 ) TJD, 'AFTER LAST', TEND
          T(MIDDL) = TEND + TINT
          IERR = 2
          GO TO 77

      ELSE IF ( ITYP .EQ. 2 ) THEN

*         -------------------------------------------------------------
*         READ JULIAN DATES AND CIO RA VALUES FROM UNFORMATTED
*         DIRECT ACCESS INPUT FILE

*         EACH RECORD OF THE FILE (EXCEPT THE FIRST) MUST CONTAIN A
*         TDB JULIAN DATE AND A CORRESPONDING CIO RA VALUE (WRT GCRS)
*         IN ARCSECONDS

*         THE JULIAN DATES MUST INCREASE BY A FIXED INTERVAL

*         THE FIRST RECORD OF THE FILE IS SPECIAL AND MUST CONTAIN THE
*         TOTAL NUMBER OF RECORDS IN THE FILE
*         -------------------------------------------------------------

*         IF FIRST TIME, OPEN FILE AND READ INITIAL VALUES
          NTIMES = NTIMES + 1
          IF ( NTIMES .EQ. 1 ) THEN
              INQUIRE ( FILE=FILNAM, EXIST=FILEOK )
              IF ( .NOT. FILEOK ) THEN
                  WRITE ( *, 2 ) FILNAM
                  IERR = 4
                  GO TO 77
              END IF
              OPEN ( UNIT=LU, FILE=FILNAM, FORM='UNFORMATTED',
     .               ACCESS='DIRECT', RECL=16, STATUS='OLD' )
              READ ( UNIT=LU, REC=1 ) NRECS
              DO 59 I = 1, NVALS
                  READ ( UNIT=LU, REC=I+1 ) T(I), R(I)
  59          CONTINUE
              TINT = NINT ( ( T(2) - T(1) ) * 1000.D0 ) / 1000.D0
              TBEG = T(MIDDL)
              TEND = T(1) + ( NRECS - NVALS + MIDDL ) * TINT
              T1 = T(1)
              LREC = 1
              MAXREC = NRECS - NVALS + 1
              IF ( TJD .LT. TBEG .OR. TJD .GT. TEND ) GO TO 10
          END IF

*         -------------------------------------------------------------

*         FILE READ SEQUENCE

  60      DIF = ( TJD - T(MIDDL) ) / TINT
*         IREC IS THE DATA RECORD NUMBER (PHYSICAL RECORD NUMBER - 1)
*         OF THE FIRST RECORD IN THE SEQUENCE OF NVALS RECORDS WITH
*         THE RELEVANT DATA TO BE RETURNED
          IREC = ( ( TJD - T1 ) / TINT ) - MIDDL + 1.9D0
          IF ( IREC .LT. 1      ) IREC = 1
          IF ( IREC .GT. MAXREC ) IREC = MAXREC

*         BASIC DECISION ON HOW TO READ FILE
          IF ( DIF .GE. -0.1D0 .AND. DIF .LE. 1.1D0 ) THEN
*             NO FILE READ NECESSARY, DATA PREVIOUSLY READ CAN BE USED
              GO TO 70
          ELSE IF ( IREC .GT. LREC .AND. IREC - LREC .LE. MIDDL ) THEN
*             REQUESTED JULIAN DATE JUST AHEAD OF DATA PREVIOUSLY READ
              GO TO 62
          ELSE
*             REQUESTED JULIAN DATE IN DIFFERENT PART OF FILE
              GO TO 66
          END IF

*         ADVANCE ARRAY DATA AND READ ONE MORE RECORD
  62      DO 64 I = 1, NVALS - 1
              T(I) = T(I+1)
              R(I) = R(I+1)
  64      CONTINUE
          READ ( UNIT=LU, REC=LREC+NVALS+1 ) T(NVALS), R(NVALS)
          LREC = LREC + 1
          GO TO 60

*         GO DIRECTLY TO PROPER RECORD RANGE AND FILL UP ARRAYS
  66      DO 68 I = 1, NVALS
              READ ( UNIT=LU, REC=IREC+I ) T(I), R(I)
  68      CONTINUE
          LREC = IREC

*         -------------------------------------------------------------

      END IF

*     GOT DATA, SO FILL OUTPUT ARRAYS
  70  DO 75 I = 1, NVALS
          TLIST(I) = T(I)
          RALIST(I) = R(I)
  75  CONTINUE
      IERR = 0

  77  RETURN


      ENTRY CIOFIL ( LUNIT, FILEID, ITYPE )
*
*     THIS ENTRY ALLOWS SPECIFICATION OF THE LOGICAL UNIT NUMBER AND
*     FILE IDENTIFIER OF THE INPUT FILE CONTAINING A TIME SERIES OF CIO
*     RA VALUES.
*
*          LUNIT  = LOGIAL UNIT NUMBER TO BE USED FOR FILE (IN)
*          FILEID = FILE ID (IN)
*          ITYPE  = TYPE OF FILE (IN)
*                   SET ITYPE=1 FOR FORMATTED SEQUENTIAL FILE
*                   SET ITYPE=2 FOR UNFORMATTED BINARY FILE
*                   SET ITYPE=0 OR ANYTHING ELSE TO CLOSE THE CURRENT
*                               FILE (LUNIT AND FILEID IGNORED)
*
*     NOTE:  AFTER A CALL TO CIOFIL WITH ITYPE=0, THE ORIGINAL OR A
*     DIFFERENT FILE OF CIO RA VALUES CAN BE ACCESSED BY SUBSEQUENT
*     CALLS TO CIORD, BUT ONLY AFTER ANOTHER CALL TO CIOFIL WITH
*     ITYPE=1 OR 2.
*
*
      IF ( ITYPE .EQ. 1 .OR. ITYPE .EQ. 2 ) THEN
          LU = LUNIT
          FILNAM = FILEID
          ITYP = ITYPE
      ELSE
          CLOSE ( UNIT=LU )
          NTIMES = 0
          TBEG = 0.D0
          TEND = 1.D10
          FILEOK = .FALSE.  
      END IF

      RETURN

      END
  


      SUBROUTINE CIOBAS ( TJD, RACIO, K,   X, Y, Z )
*
*     THIS SUBROUTINE RETURNS THE ORTHONORMAL BASIS VECTORS, WITH
*     RESPECT TO THE GCRS (GEOCENTRIC ICRS), OF THE CELESTIAL
*     INTERMEDIATE SYSTEM DEFINED BY THE CELESTIAL INTERMEDIATE POLE
*     (CIP) (IN THE Z DIRECTION) AND THE CELESTIAL INTERMEDIATE ORIGIN
*     (CIO) (IN THE X DIRECTION).  A TDB JULIAN DATE AND THE RIGHT
*     ASCENSION OF THE CIO AT THAT DATE IS REQUIRED AS INPUT.  THE
*     RIGHT ASCENSION OF THE CIO CAN BE WITH RESPECT TO EITHER THE
*     GCRS ORIGIN OR THE TRUE EQUINOX OF DATE -- DIFFERENT ALGORITHMS
*     ARE USED IN THE TWO CASES. 
*
*          TJD    = TDB JULIAN DATE (IN)
*          RACIO  = RIGHT ASCENSION OF THE CIO, IN HOURS (IN)
*          K      = REFERENCE SYSTEM IN WHICH RIGHT ASCENSION IS
*                   EXPRESSED (IN)
*                   SET K=1 FOR GCRS
*                   SET K=2 FOR TRUE EQUATOR AND EQUINOX OF DATE
*          X      = UNIT VECTOR TOWARD THE CIO, EQUATORIAL RECTANGULAR
*                   COORDINATES, REFERRED TO THE GCRS (OUT)
*          Y      = UNIT VECTOR TOWARD THE Y-DIRECTION, EQUATORIAL
*                   RECTANGULAR COORDINATES, REFERRED TO THE GCRS (OUT)
*          Z      = UNIT VECTOR TOWARD NORTH CELESTIAL POLE (CIP),
*                   EQUATORIAL RECTANGULAR COORDINATES, REFERRED TO
*                   THE GCRS (OUT)
*
*
      DOUBLE PRECISION TJD,RACIO,X,Y,Z,XX,YY,ZZ,W0,W1,W2,Z0,
     .     PI,RADCON,T0,TLAST,SINRA,COSRA,XMAG,DABS,DSIN,DCOS,DSQRT
      DIMENSION X(3), Y(3), Z(3), XX(3), YY(3), ZZ(3), Z0(3),
     .     W0(3), W1(3), W2(3)
      SAVE

      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( RADCON = PI / 180.D0           )

*     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
      DATA T0 / 2451545.00000000D0 /
      DATA Z0 / 0.D0, 0.D0, 1.D0 /,   TLAST / 0.D0 /,   KLAST / 0 /
      
   3  FORMAT ( ' CIOBAS ERROR: INVALID VALUE FOR K FOR JD ',
     .     F10.1 )      

*     USE LAST-COMPUTED BASIS VECTORS IF POSSIBLE
      IF ( DABS ( TJD - TLAST ) .LE. 1.D-8 .AND. K .EQ. KLAST )
     .   GO TO 60

*     COMPUTE UNIT VECTOR Z TOWARD CELESTIAL POLE (CIP)      
      CALL NUTATE ( -TJD, Z0,   W1 )
      CALL PRECES ( TJD, W1, T0,   W2 )
      CALL FRAME ( W2, -1,   ZZ )
      
* --- RA OF CIO EXPRESSED IN GCRS -------------------------------------      
      
      IF ( K .EQ. 1 ) THEN

*         COMPUTE VECTOR X TOWARD CIO IN GCRS
          SINRA = DSIN ( RACIO * 15.D0 * RADCON )
          COSRA = DCOS ( RACIO * 15.D0 * RADCON )
          XX(1) =  ZZ(3) * COSRA
          XX(2) =  ZZ(3) * SINRA
          XX(3) = -ZZ(1) * COSRA - ZZ(2) * SINRA

*         NORMALIZE VECTOR X 
          XMAG = DSQRT ( XX(1)**2 + XX(2)**2 + XX(3)**2 )
          XX(1) = XX(1) / XMAG
          XX(2) = XX(2) / XMAG
          XX(3) = XX(3) / XMAG

*         COMPUTE UNIT VECTOR Y ORTHOGONAL TO X AND Z (Y = Z CROSS X)
          YY(1) = ZZ(2) * XX(3) - ZZ(3) * XX(2)
          YY(2) = ZZ(3) * XX(1) - ZZ(1) * XX(3)
          YY(3) = ZZ(1) * XX(2) - ZZ(2) * XX(1)
          
* --- RA OF CIO EXPRESSED IN EQUATOR-AND-EQUINOX OF DATE SYSTEM -------          
      
      ELSE IF ( K .EQ. 2 ) THEN
      
*         CONSTRUCT UNIT VECTOR TOWARD CIO
*         IN EQUATOR-AND-EQUINOX-OF-DATE SYSTEM
          W0(1) = DCOS ( RACIO * 15.D0 * RADCON )
          W0(2) = DSIN ( RACIO * 15.D0 * RADCON )
          W0(3) = 0.D0
       
*         ROTATE THE VECTOR INTO THE GCRS TO FORM UNIT VECTOR X      
          CALL NUTATE ( -TJD, W0,   W1 )
          CALL PRECES ( TJD, W1, T0,   W2 )
          CALL FRAME ( W2, -1,   XX )
          
*         COMPUTE UNIT VECTOR Y ORTHOGONAL TO X AND Z (Y = Z CROSS X)
          YY(1) = ZZ(2) * XX(3) - ZZ(3) * XX(2)
          YY(2) = ZZ(3) * XX(1) - ZZ(1) * XX(3)
          YY(3) = ZZ(1) * XX(2) - ZZ(2) * XX(1)
      
* ---------------------------------------------------------------------      
      
      ELSE
          
          WRITE ( *, 3 ) TJD
          GO TO 77
      
      END IF   

* ---------------------------------------------------------------------      

      TLAST = TJD
      KLAST = K
      
  60  DO 66 J = 1, 3
          X(J) = XX(J)
          Y(J) = YY(J)
          Z(J) = ZZ(J)
  66  CONTINUE    

  77  RETURN

      END



      SUBROUTINE EROT (DATE1,DATE2,THETA)
*
*     THIS SUBROUTINE RETURNS THE VALUE OF THE EARTH ROTATION ANGLE
*     (THETA) FOR A GIVEN UT1 JULIAN DATE.  THE EXPRESSION USED IS
*     TAKEN FROM THE NOTE TO IAU RESOLUTION B1.8 OF 2000.
*
*         DATE1  = HIGH-ORDER PART OF UT1 JULIAN DATE (IN)
*         DATE2  = LOW-ORDER PART OF UT1 JULIAN DATE (IN)
*         THETA  = EARTH ROTATION ANGLE IN DEGREES (OUT)
*
*
      DOUBLE PRECISION DATE1, DATE2, THETA, T0, THET1, THET2, THET3,
     .     DMOD

      DATA T0 / 2451545.0D0 /

*     THE ALGORITHM USED BELOW IS EQUIVALENT TO THE CANNONICAL
*     THETA = 0.7790572732640D0 + 1.00273781191135448D0 * T,
*     WHERE T IS THE TIME IN UT1 DAYS FROM 2000.0 (T=DATE1+DATE2-T0),
*     BUT IT AVOIDS MANY TWO-PI 'WRAPS' THAT DECREASE PRECISION
*     (ADOPTED FROM SOFA ROUTINE IAU_ERA00 BY PAT WALLACE; SEE ALSO
*     EXPRESSION AT TOP OF PAGE 35 OF IERS CONVENTIONS (1996))

      THET1 = 0.7790572732640D0 + 0.00273781191135448D0 * ( DATE1 - T0 )
      THET2 =                     0.00273781191135448D0 *   DATE2
      THET3 = DMOD ( DATE1, 1.D0 ) + DMOD ( DATE2, 1.D0 )
      THETA = DMOD ( THET1 + THET2 + THET3, 1.D0 ) * 360.D0
      IF ( THETA .LT. 0.D0 ) THETA = THETA + 360.D0

      RETURN

      END



      SUBROUTINE EQXRA ( TJD, K,    RAEQ )
*
*     THIS SUBROUTINE COMPUTES THE INTERMEDIATE RIGHT ASCENSION
*     OF THE EQUINOX AT JULIAN DATE TJD, USING AN ANALYTICAL EXPRESSION
*     FOR THE ACCUMULATED PRECESSION IN RIGHT ASCENSION.  FOR THE
*     TRUE EQUINOX THE RESULT IS THE EQUATION OF THE ORIGINS. 
*
*          TJD    = TDB JULIAN DATE (IN)
*          K      = EQUINOX SELECTION CODE (IN)
*                   SET K=0 FOR MEAN EQUINOX
*                   SET K=1 FOR TRUE EQUINOX (EQUATION OF THE ORIGINS)
*          RADIF  = INTERMEDIATE RIGHT ASCENSION OF THE EQUINOX,
*                   IN HOURS (+ OR -) (OUT)
*
*
      DOUBLE PRECISION TJD,RAEQ,T0,TLAST,EE,EQEQ,T,A,PRECRA,DABS
      SAVE 

*     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
      DATA T0 / 2451545.00000000D0 /
      DATA TLAST / 0.D0 /,   EE / 0.D0 / 
      
      T = ( TJD - T0 ) / 36525.D0	

*     FOR THE TRUE EQUINOX, OBTAIN THE EQUATION OF THE EQUINOXES IN
*     TIME SECONDS, WHICH INCLUDES THE 'COMPLIMENTARY TERMS'
      IF ( K .EQ. 1 ) THEN
          IF ( DABS ( TJD - TLAST ) .GT. 1.D-8 ) THEN
              CALL ETILT ( TJD,   A, A, EE, A, A )
              TLAST = TJD
          END IF
          EQEQ = EE 
      ELSE
          EQEQ = 0.D0
      END IF
      
*     PRECESSION IN RA IN ARCSECONDS TAKEN FROM CAPITAINE ET AL. (2003),
*     ASTRONOMY AND ASTROPHYSICS 412, 567-586, EQ. (42)
      PRECRA = 0.014506D0 +
     .         ( ( ( ( -    0.0000000368D0   * T
     .                 -    0.000029956D0  ) * T
     .                 -    0.00000044D0   ) * T
     .                 +    1.3915817D0    ) * T
     .                 + 4612.156534D0     ) * T

      RAEQ = - ( PRECRA / 15.D0 + EQEQ ) / 3600.D0

      RETURN
      
      END



      SUBROUTINE SETDT ( DELT )
*
*     THIS SUBROUTINE ALLOWS FOR THE SPECIFICATION OF THE DELTA-T VALUE
*     (DELTA-T = TT - UT1) TO BE USED IN THE CALCULATION OF SIDEREAL
*     TIME AND THE TERRESTRIAL-TO-CELESTIAL TRANSFORMATION.  IT ALLOWS
*     THESE CALCULATIONS TO BE PERFORMED, CORRECTLY, USING UT1 AS THE
*     TIME ARGUMENT FOR THE EARTH ROTATION ANGLE AND TDB AS THE TIME
*     ARGUMENT FOR THE PRECESSION AND NUTATION COMPONENTS.  THIS
*     SUBROUTINE, IF USED, SHOULD BE CALLED BEFORE ANY SUBROUTINE
*     RELATED TO EARTH ROTATION (E.G., SIDTIM OR TERCEL) FOR A GIVEN
*     DATE.  THE VALUE OF DELTA-T SPECIFIED HERE WILL BE USED UNTIL
*     EXPLICITLY CHANGED.
*
*          DELT   = VALUE OF DELTA-T (TT-UT1) IN SECONDS (IN)
*
*     NOTE 1:  THE COMPUTED VALUE OF SIDEREAL TIME, AND THE EQUIVALENT
*     EARTH ORIENTATION ANGLES, ARE RELATIVELY INSENSITIVE TO THE VALUE
*     OF DELTA-T: UP TO ONLY ~3 MICROARCSECONDS PER SECOND OF DELTA-T.
*     THEREFORE, FOR MANY APPLICATIONS, THIS SUBROUTINE EITHER NEED NOT
*     BE CALLED AT ALL, OR CAN BE CALLED JUST ONCE FOR A WIDE RANGE OF
*     DATES (E.G., A YEAR).  IF THIS CALL IS NOT USED, A DEFAULT
*     VALUE OF DELTA-T OF 64 SECONDS IS USED, WHICH IS APPROPRIATE TO
*     2000.0.
*
*     NOTE 2:  THE INPUT TIME ARGUMENTS TO SIDTIM AND TERCEL (TJDH AND
*     TJDL) ARE EXPRESSED IN UT1 REGARDLESS OF WHETHER THIS CALL IS
*     USED.
*
*
      DOUBLE PRECISION DELTAT, DT, DELT
      SAVE DT

*     DEFAULT VALUE OF DELTA-T IN DAYS, EQUIVALENT TO 64 SECONDS,
*     THE APPROXIMATE VALUE AT 2000.0
      DATA DT / 0.00074074D0 /

      DT = DELT / 86400.D0

      RETURN


      ENTRY GETDT ( DELTAT )

*     THIS ENTRY RETURNS THE CURRENT VALUE OF DELTA-T
*     (DELTA-T = TT - UT1), PREVIOUSLY SET BY THE USER.  THE VALUE
*     RETURNED IS TO BE USED IN THE CALCULATION OF SIDEREAL TIME AND
*     THE TERRESTRIAL-TO-CELESTIAL TRANSFORMATION.  IT ALLOWS THESE
*     CALCULATIONS TO BE PERFORMED, CORRECTLY, USING UT1 AS THE TIME
*     ARGUMENT FOR THE EARTH ROTATION ANGLE AND TDB AS THE TIME ARGUMENT
*     FOR THE PRECESSION AND NUTATION COMPONENTS.
*
*          DELTAT = VALUE OF DELTA-T (TT-UT1) IN DAYS (OUT)


      DELTAT = DT

      RETURN

      END



      SUBROUTINE SETMOD ( MODE )
*
*     THIS SUBROUTINE ALLOWS THE USER TO SPECIFY THE 'MODE' VALUE,
*     WHICH DETERMINES THE METHOD USED FOR THE COMPUTATION OF SIDEREAL
*     TIME AND THE TERRESTRIAL-TO-CELESTIAL TRANSFORMATION, AND THE
*     ACCURACY OF NUTATION AND RELATED CALCULATIONS.
*
*          MODE   = SELECTION FOR METHOD AND ACCURACY (IN)
*                   SET MODE=0 FOR CIO-BASED METHOD, FULL ACCURACY
*                   SET MODE=1 FOR CIO-BASED METHOD, REDUCED ACCURACY
*                   SET MODE=2 FOR EQUINOX-BASED METHOD, FULL ACCURACY
*                   SET MODE=3 FOR EQUINOX-BASED METHOD, REDUCED
*                                  ACCURACY
*
*     NOTE: OTHER ENTRY POINTS ARE PROVIDED TO ALLOW THE METHOD AND
*     ACCURACY TO BE SPECIFIED IN A MORE OBVIOUS WAY:
*     MODE=0 CAN BE SET BY CALL CIOTIO AND CALL HIACC
*     MODE=1 CAN BE SET BY CALL CIOTIO AND CALL LOACC
*     MODE=2 CAN BE SET BY CALL EQINOX AND CALL HIACC
*     MODE=3 CAN BE SET BY CALL EQINOX AND CALL LOACC
*
*
      SAVE IMODE, LMODE

      DATA IMODE, LMODE / 2, 2 /

      LMODE = IMODE
      IMODE = MODE

      RETURN


      ENTRY CIOTIO
      LMODE = IMODE
      IF ( IMODE .GE. 2 ) IMODE = IMODE - 2
      RETURN


      ENTRY EQINOX
      LMODE = IMODE
      IF ( IMODE .LE. 1 ) IMODE = IMODE + 2
      RETURN


      ENTRY LOACC
      LMODE = IMODE
      IF ( MOD ( IMODE, 2 ) .EQ. 0 ) IMODE = IMODE + 1
      RETURN


      ENTRY HIACC
      LMODE = IMODE
      IF ( MOD ( IMODE, 2 ) .EQ. 1 ) IMODE = IMODE - 1
      RETURN


      ENTRY RESUME
      IMODE = LMODE
      RETURN


      ENTRY GETMOD ( MODE )
*
*     THIS SUBROUTINE RETURNS THE 'MODE' VALUE, WHICH
*     DETERMINES THE METHOD USED FOR THE COMPUTATION OF SIDEREAL
*     TIME AND THE TERRESTRIAL-TO-CELESTIAL TRANSFORMATION, AND THE
*     ACCURACY OF NUTATION AND RELATED CALCULATIONS.
*
*          MODE   = SELECTION FOR METHOD AND ACCURACY (OUT)
*                   MODE=0 MEANS CIO-BASED METHOD, FULL ACCURACY
*                   MODE=1 MEANS CIO-BASED METHOD, REDUCED ACCURACY
*                   MODE=2 MEANS EQUINOX-BASED METHOD, FULL ACCURACY
*                   MODE=3 MEANS EQUINOX-BASED METHOD, REDUCED ACCURACY
*
*
      MODE = IMODE

      RETURN

      END



      SUBROUTINE GETVEC ( UNITV )
*
*     THIS SUBROUTINE ALLOWS THE USER TO RETRIEVE THE LAST COMPUTED
*     POSITION ON THE SKY AS A UNIT VECTOR.
*
*          UNITV  = UNIT VECTOR TOWARD LAST COMPUTED POSITION ON THE
*                   SKY, IN THE COORDINATE SYSTEM USED FOR THAT
*                   POSITION (OUT)
*
*
      DOUBLE PRECISION UNITV, P, POS, R, DSQRT
      DIMENSION UNITV(3), P(3), POS(3)
      SAVE P

      R = DSQRT ( P(1)**2 + P(2)**2 + P(3)**2 )

      DO 20 J = 1, 3
          UNITV(J) = P(J) / R
  20  CONTINUE

      RETURN


      ENTRY SETVEC ( POS )
*
*     THIS ENTRY STORES THE LAST COMPUTED POSITION ON THE SKY.
*
*          POS    = VECTOR TOWARD LAST COMPUTED POSITION ON THE
*                   SKY, IN THE COORDINATE SYSTEM USED FOR THAT
*                   POSITION (IN)
*
*
      DO 30 J = 1, 3
          P(J) = POS(J)
  30  CONTINUE

      RETURN

      END



      SUBROUTINE JULDAT (I,M,K,H,TJD)
*
*     THIS SUBROUTINE COMPUTES JULIAN DATE, GIVEN CALENDAR DATE AND
*     TIME.  INPUT CALENDAR DATE MUST BE GREGORIAN.  INPUT TIME VALUE
*     CAN BE IN ANY UT-LIKE TIME SCALE (UTC, UT1, TT, ETC.) - OUTPUT
*     JULIAN DATE WILL HAVE SAME BASIS.  ALGORITHM BY FLIEGEL AND
*     VAN FLANDERN.
*
*          I      = YEAR (IN)
*          M      = MONTH NUMBER (IN)
*          K      = DAY OF MONTH (IN)
*          H      = UT HOURS (IN)
*          TJD    = JULIAN DATE (OUT)
*
*
      DOUBLE PRECISION H,TJD

*     JD=JULIAN DAY NO FOR DAY BEGINNING AT GREENWICH NOON ON GIVEN DATE
      JD = K-32075+1461*(I+4800+(M-14)/12)/4+367*(M-2-(M-14)/12*12)/12
     .     -3*((I+4900+(M-14)/12)/100)/4
      TJD = JD - 0.5D0 + H/24.D0

      RETURN
      END



      SUBROUTINE CALDAT (TJD,I,M,K,H)
*
*     THIS SUBROUTINE COMPUTES CALENDAR DATE AND TIME, GIVEN JULIAN
*     DATE.  INPUT JULIAN DATE CAN BE BASED ON ANY UT-LIKE TIME SCALE
*     (UTC, UT1, TT, ETC.) - OUTPUT TIME VALUE WILL HAVE SAME BASIS.
*     OUTPUT CALENDAR DATE WILL BE GREGORIAN.  ALGORITHM BY FLIEGEL AND
*     VAN FLANDERN.
*
*          TJD    = JULIAN DATE (IN)
*          I      = YEAR (OUT)
*          M      = MONTH NUMBER (OUT)
*          K      = DAY OF MONTH (OUT)
*          H      = UT HOURS (OUT)
*
*
      DOUBLE PRECISION TJD,H,DJD,DMOD

      DJD = TJD + 0.5D0
      JD = DJD
      H = DMOD (DJD,1.D0) * 24.D0
*     JD=JULIAN DAY NO FOR DAY BEGINNING AT GREENWICH NOON ON GIVEN DATE
      L = JD + 68569
      N = 4*L/146097
      L = L - (146097*N+3)/4
*     I=YEAR, M=MONTH, K=DAY
      I = 4000*(L+1)/1461001
      L = L - 1461*I/4 + 31
      M = 80*L/2447
      K = L - 2447*M/80
      L = M / 11
      M = M + 2 - 12*L
      I = 100*(N-49) + I + L

      RETURN
      END



      SUBROUTINE ASTCON (NAME,FACTOR,CONST)
*
*     THIS SUBROUTINE SUPPLIES THE VALUES OF ASTRONOMICAL CONSTANTS.
*
*         NAME   = NAME OF CONSTANT WHOSE VALUE IS DESIRED (IN)
*                  'C'         SPEED OF LIGHT IN METERS/SECOND
*                  'C(AU/DAY)' SPEED OF LIGHT IN AU/DAY
*                  'AU'        LENGTH OF ASTRONOMICAL UNIT IN METERS
*                  'AU(SEC)'   LENGTH OF ASTRONOMICAL UNIT IN SECONDS
*                  'GS'        HELIOCENTRIC GRAVITATIONAL CONSTANT
*                                 IN METERS**3/SECOND**2
*                  'GE'        GEOCENTRIC GRAVITATIONAL CONSTANT
*                                 IN METERS**3/SECOND**2
*                  'ERAD'      EQUATORIAL RADIUS OF EARTH IN METERS
*                  'F'         FLATTENING FACTOR OF EARTH
*                  'ANGVEL'    NOMINAL MEAN ROTATIONAL ANGULAR VELOCITY
*                                 OF EARTH IN RADIANS/SECOND
*                  'MASS_SUN'  RECIPROCAL MASS OF THE SUN
*                  'MASS_EAR'  RECIPROCAL MASS OF THE EARTH
*                  'MASS_MOO'  RECIPROCAL MASS OF THE MOON
*                  'MASS_MER'  RECIPROCAL MASS OF MERCURY
*                      :             :      :        :
*                  'MASS_PLU'  RECIPROCAL MASS OF PLUTO
*         FACTOR = FACTOR BY WHICH CONSTANT VALUE IS TO BE MULTIPLIED
*                  (IN)
*         CONST  = CONSTANT VALUE AFTER MULTIPLICATION BY FACTOR (OUT)
*
*
      DOUBLE PRECISION FACTOR,CONST,C,AUSEC
      CHARACTER NAME*(*)

*     NOTE:  THESE CONSTANT VALUES ARE BASED ON THE TDB SECOND WHERE
*     APPLICABLE.

*     SPEED OF LIGHT IN METERS/SECOND IS A DEFINING PHYSICAL CONSTANT
      DATA C / 299792458.D0 /

*     LIGHT-TIME FOR ONE ASTRONOMICAL UNIT IN SECONDS, FROM DE-405
      DATA AUSEC / 499.0047838061D0 /

*     SPEED OF LIGHT IN METERS/SECOND
      IF ( NAME .EQ. 'C' ) THEN
          CONST = C

*     SPEED OF LIGHT IN AU/DAY
      ELSE IF ( NAME .EQ. 'C(AU/DAY)' ) THEN
          CONST = 86400.D0 / AUSEC

*     LENGTH OF ASTRONOMICAL UNIT IN METERS
      ELSE IF ( NAME .EQ. 'AU' ) THEN
          CONST = AUSEC * C

*     LENGTH OF ASTRONOMICAL UNIT IN SECONDS
      ELSE IF ( NAME .EQ. 'AU(SEC)' ) THEN
          CONST = AUSEC

*     HELIOCENTRIC GRAVITATIONAL CONSTANT IN METERS**3/SECOND**2, FROM
*     DE-405
      ELSE IF ( NAME .EQ. 'GS' ) THEN
          CONST = 1.32712440017987D20

*     GEOCENTRIC GRAVITATIONAL CONSTANT IN METERS**3/SECOND**2, FROM
*     DE-405
      ELSE IF ( NAME .EQ. 'GE' ) THEN
          CONST = 3.98600433D14

*     EQUATORIAL RADIUS OF EARTH IN METERS, FROM IERS CONVENTIONS (2003)
      ELSE IF ( NAME .EQ. 'ERAD' ) THEN
          CONST = 6378136.6D0

*     FLATTENING FACTOR OF EARTH, FROM IERS CONVENTIONS (2003)
      ELSE IF ( NAME .EQ. 'F' ) THEN
          CONST = 1.D0 / 298.25642D0

*     NOMINAL MEAN ROTATIONAL ANGULAR VELOCITY OF EARTH
*     IN RADIANS/SECOND, FROM IERS CONVENTIONS (2003)
      ELSE IF ( NAME .EQ. 'ANGVEL' ) THEN
          CONST = 7.2921150D-5

*     RECIPROCAL MASSES OF SOLAR SYSTEM BODIES, FROM DE-405
*     (SUN MASS / BODY MASS)
      ELSE IF ( NAME(1:4) .EQ. 'MASS' ) THEN

          CONST = 1.D0
          IF ( NAME(6:8) .EQ. 'SUN' ) CONST =         1.D0
          IF ( NAME(6:8) .EQ. 'MOO' ) CONST =  27068700.387534D0
          IF ( NAME(6:8) .EQ. 'MER' ) CONST =   6023600.D0
          IF ( NAME(6:8) .EQ. 'VEN' ) CONST =    408523.71D0
          IF ( NAME(6:8) .EQ. 'EAR' ) CONST =    332946.050895D0
          IF ( NAME(6:8) .EQ. 'MAR' ) CONST =   3098708.D0
          IF ( NAME(6:8) .EQ. 'JUP' ) CONST =      1047.3486D0
          IF ( NAME(6:8) .EQ. 'SAT' ) CONST =      3497.898D0
          IF ( NAME(6:8) .EQ. 'URA' ) CONST =     22902.98D0
          IF ( NAME(6:8) .EQ. 'NEP' ) CONST =     19412.24D0
          IF ( NAME(6:8) .EQ. 'PLU' ) CONST = 135200000.D0
          IF ( NAME(6:8) .EQ. 'EMB' ) CONST =    328900.561400D0

      END IF

      CONST = CONST * FACTOR

      RETURN

      END



      SUBROUTINE NOD (T,DPSI,DEPS)
*
*     THIS SUBROUTINE RETURNS THE VALUES FOR NUTATION IN LONGITUDE AND
*     NUTATION IN OBLIQUITY FOR A GIVEN TDB JULIAN DATE.
*
*          T     = TDB TIME IN JULIAN CENTURIES SINCE J2000.0 (IN)
*          DPSI  = NUTATION IN LONGITUDE IN ARCSECONDS (OUT)
*          DEPS  = NUTATION IN OBLIQUITY IN ARCSECONDS (OUT)
*
*
      DOUBLE PRECISION T,DPSI,DEPS,PI,SECCON,T0,T1,DP,DE
      SAVE

      PARAMETER ( PI     = 3.14159265358979324D0 )
      PARAMETER ( SECCON = 180.D0 * 3600.D0 / PI )

*     T0 = TDB JULIAN DATE OF EPOCH J2000.0 (TT)
      DATA T0 / 2451545.00000000D0 /

*     GET METHOD/ACCURACY MODE
      CALL GETMOD ( MODE )

      T1 = T * 36525.D0

*     =================================================================
*     EVALUATE NUTATION SERIES
*     RESULTING NUTATION IN LONGITUDE AND OBLIQUITY IN ARCSECONDS

*     CALL SUBROUTINE TO EVALUATE NUTATION SERIES
      IF ( MOD ( MODE, 2 ) .EQ. 0 ) THEN
*         HIGH ACCURACY MODE -- IERS SUBROUTINE
          CALL NU2000A ( T0, T1, DP, DE )
      ELSE
*         LOW ACCURACY MODE -- MODIFICATION OF IERS SUBROUTINE
          CALL NU2000K ( T0, T1, DP, DE )
      END IF
      DPSI = DP * SECCON
      DEPS = DE * SECCON

*     =================================================================

      RETURN

      END



      SUBROUTINE NU2000A ( DATE1, DATE2, DPSI, DEPS )

*+
*  - - - - - - - -
*   N U 2 0 0 0 A
*  - - - - - - - -
*
*  Nutation, IAU 2000A model (MHB_2000 without FCN).
*
*  Annexe to IERS Conventions 2000, Chapter 5
*
*  Given:
*     DATE1,DATE2    d   TT date (JD = DATE1+DATE2)
*
*  Returned:
*     DPSI,DEPS      d   nutation (luni-solar + planetary, radians)
*
*  This revision:  2002 November 25
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION DATE1, DATE2, DPSI, DEPS

*  Arcseconds to radians
      DOUBLE PRECISION DAS2R
      PARAMETER ( DAS2R = 4.848136811095359935899141D-6 )

*  Milliarcseconds to radians
      DOUBLE PRECISION DMAS2R
      PARAMETER ( DMAS2R = DAS2R / 1D3 )

*  Arc seconds in a full circle
      DOUBLE PRECISION TURNAS
      PARAMETER ( TURNAS = 1296000D0 )

*  2Pi
      DOUBLE PRECISION D2PI
      PARAMETER ( D2PI = 6.283185307179586476925287D0 )

*  Units of 0.1 microarcsecond to radians
      DOUBLE PRECISION U2R
      PARAMETER ( U2R = DAS2R/1D7 )

*  Reference epoch (J2000), JD
      DOUBLE PRECISION DJ0
      PARAMETER ( DJ0 = 2451545D0 )

*  Days per Julian century
      DOUBLE PRECISION DJC
      PARAMETER ( DJC = 36525D0 )

*  Miscellaneous
      DOUBLE PRECISION T, EL, ELP, F, D, OM, ARG, DP, DE, SARG, CARG,
     :                 DPSILS, DEPSLS,
     :                 AL, ALSU, AF, AD, AOM, ALME, ALVE, ALEA, ALMA,
     :                 ALJU, ALSA, ALUR, ALNE, APA, DPSIPL, DEPSPL
      INTEGER I, J

*  -------------------------
*  Luni-Solar nutation model
*  -------------------------

*  Number of terms in the luni-solar nutation model
      INTEGER NLS
      PARAMETER ( NLS = 678 )

*  Coefficients for fundamental arguments
      INTEGER NALS(5,NLS)

*  Longitude and obliquity coefficients
      DOUBLE PRECISION CLS(6,NLS)

*  ---------------
*  Planetary terms
*  ---------------

*  Number of terms in the planetary nutation model
      INTEGER NPL
      PARAMETER ( NPL = 687 )

*  Coefficients for fundamental arguments
      INTEGER NAPL(14,NPL)

*  Longitude and obliquity coefficients
      INTEGER ICPL(4,NPL)

*  ----------------------------------------
*  Tables of argument and term coefficients
*  ----------------------------------------

*
*  Luni-Solar argument multipliers:
*               L     L'    F     D     Om

      DATA ( ( NALS(I,J), I=1,5 ), J=  1, 10 ) /
     :          0,    0,    0,    0,    1,
     :          0,    0,    2,   -2,    2,
     :          0,    0,    2,    0,    2,
     :          0,    0,    0,    0,    2,
     :          0,    1,    0,    0,    0,
     :          0,    1,    2,   -2,    2,
     :          1,    0,    0,    0,    0,
     :          0,    0,    2,    0,    1,
     :          1,    0,    2,    0,    2,
     :          0,   -1,    2,   -2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J= 11, 20 ) /
     :          0,    0,    2,   -2,    1,
     :         -1,    0,    2,    0,    2,
     :         -1,    0,    0,    2,    0,
     :          1,    0,    0,    0,    1,
     :         -1,    0,    0,    0,    1,
     :         -1,    0,    2,    2,    2,
     :          1,    0,    2,    0,    1,
     :         -2,    0,    2,    0,    1,
     :          0,    0,    0,    2,    0,
     :          0,    0,    2,    2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J= 21, 30 ) /
     :          0,   -2,    2,   -2,    2,
     :         -2,    0,    0,    2,    0,
     :          2,    0,    2,    0,    2,
     :          1,    0,    2,   -2,    2,
     :         -1,    0,    2,    0,    1,
     :          2,    0,    0,    0,    0,
     :          0,    0,    2,    0,    0,
     :          0,    1,    0,    0,    1,
     :         -1,    0,    0,    2,    1,
     :          0,    2,    2,   -2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J= 31, 40 ) /
     :          0,    0,   -2,    2,    0,
     :          1,    0,    0,   -2,    1,
     :          0,   -1,    0,    0,    1,
     :         -1,    0,    2,    2,    1,
     :          0,    2,    0,    0,    0,
     :          1,    0,    2,    2,    2,
     :         -2,    0,    2,    0,    0,
     :          0,    1,    2,    0,    2,
     :          0,    0,    2,    2,    1,
     :          0,   -1,    2,    0,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J= 41, 50 ) /
     :          0,    0,    0,    2,    1,
     :          1,    0,    2,   -2,    1,
     :          2,    0,    2,   -2,    2,
     :         -2,    0,    0,    2,    1,
     :          2,    0,    2,    0,    1,
     :          0,   -1,    2,   -2,    1,
     :          0,    0,    0,   -2,    1,
     :         -1,   -1,    0,    2,    0,
     :          2,    0,    0,   -2,    1,
     :          1,    0,    0,    2,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J= 51, 60 ) /
     :          0,    1,    2,   -2,    1,
     :          1,   -1,    0,    0,    0,
     :         -2,    0,    2,    0,    2,
     :          3,    0,    2,    0,    2,
     :          0,   -1,    0,    2,    0,
     :          1,   -1,    2,    0,    2,
     :          0,    0,    0,    1,    0,
     :         -1,   -1,    2,    2,    2,
     :         -1,    0,    2,    0,    0,
     :          0,   -1,    2,    2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J= 61, 70 ) /
     :         -2,    0,    0,    0,    1,
     :          1,    1,    2,    0,    2,
     :          2,    0,    0,    0,    1,
     :         -1,    1,    0,    1,    0,
     :          1,    1,    0,    0,    0,
     :          1,    0,    2,    0,    0,
     :         -1,    0,    2,   -2,    1,
     :          1,    0,    0,    0,    2,
     :         -1,    0,    0,    1,    0,
     :          0,    0,    2,    1,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J= 71, 80 ) /
     :         -1,    0,    2,    4,    2,
     :         -1,    1,    0,    1,    1,
     :          0,   -2,    2,   -2,    1,
     :          1,    0,    2,    2,    1,
     :         -2,    0,    2,    2,    2,
     :         -1,    0,    0,    0,    2,
     :          1,    1,    2,   -2,    2,
     :         -2,    0,    2,    4,    2,
     :         -1,    0,    4,    0,    2,
     :          2,    0,    2,   -2,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J= 81, 90 ) /
     :          2,    0,    2,    2,    2,
     :          1,    0,    0,    2,    1,
     :          3,    0,    0,    0,    0,
     :          3,    0,    2,   -2,    2,
     :          0,    0,    4,   -2,    2,
     :          0,    1,    2,    0,    1,
     :          0,    0,   -2,    2,    1,
     :          0,    0,    2,   -2,    3,
     :         -1,    0,    0,    4,    0,
     :          2,    0,   -2,    0,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J= 91,100 ) /
     :         -2,    0,    0,    4,    0,
     :         -1,   -1,    0,    2,    1,
     :         -1,    0,    0,    1,    1,
     :          0,    1,    0,    0,    2,
     :          0,    0,   -2,    0,    1,
     :          0,   -1,    2,    0,    1,
     :          0,    0,    2,   -1,    2,
     :          0,    0,    2,    4,    2,
     :         -2,   -1,    0,    2,    0,
     :          1,    1,    0,   -2,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=101,110 ) /
     :         -1,    1,    0,    2,    0,
     :         -1,    1,    0,    1,    2,
     :          1,   -1,    0,    0,    1,
     :          1,   -1,    2,    2,    2,
     :         -1,    1,    2,    2,    2,
     :          3,    0,    2,    0,    1,
     :          0,    1,   -2,    2,    0,
     :         -1,    0,    0,   -2,    1,
     :          0,    1,    2,    2,    2,
     :         -1,   -1,    2,    2,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=111,120 ) /
     :          0,   -1,    0,    0,    2,
     :          1,    0,    2,   -4,    1,
     :         -1,    0,   -2,    2,    0,
     :          0,   -1,    2,    2,    1,
     :          2,   -1,    2,    0,    2,
     :          0,    0,    0,    2,    2,
     :          1,   -1,    2,    0,    1,
     :         -1,    1,    2,    0,    2,
     :          0,    1,    0,    2,    0,
     :          0,   -1,   -2,    2,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=121,130 ) /
     :          0,    3,    2,   -2,    2,
     :          0,    0,    0,    1,    1,
     :         -1,    0,    2,    2,    0,
     :          2,    1,    2,    0,    2,
     :          1,    1,    0,    0,    1,
     :          1,    1,    2,    0,    1,
     :          2,    0,    0,    2,    0,
     :          1,    0,   -2,    2,    0,
     :         -1,    0,    0,    2,    2,
     :          0,    1,    0,    1,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=131,140 ) /
     :          0,    1,    0,   -2,    1,
     :         -1,    0,    2,   -2,    2,
     :          0,    0,    0,   -1,    1,
     :         -1,    1,    0,    0,    1,
     :          1,    0,    2,   -1,    2,
     :          1,   -1,    0,    2,    0,
     :          0,    0,    0,    4,    0,
     :          1,    0,    2,    1,    2,
     :          0,    0,    2,    1,    1,
     :          1,    0,    0,   -2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=141,150 ) /
     :         -1,    0,    2,    4,    1,
     :          1,    0,   -2,    0,    1,
     :          1,    1,    2,   -2,    1,
     :          0,    0,    2,    2,    0,
     :         -1,    0,    2,   -1,    1,
     :         -2,    0,    2,    2,    1,
     :          4,    0,    2,    0,    2,
     :          2,   -1,    0,    0,    0,
     :          2,    1,    2,   -2,    2,
     :          0,    1,    2,    1,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=151,160 ) /
     :          1,    0,    4,   -2,    2,
     :         -1,   -1,    0,    0,    1,
     :          0,    1,    0,    2,    1,
     :         -2,    0,    2,    4,    1,
     :          2,    0,    2,    0,    0,
     :          1,    0,    0,    1,    0,
     :         -1,    0,    0,    4,    1,
     :         -1,    0,    4,    0,    1,
     :          2,    0,    2,    2,    1,
     :          0,    0,    2,   -3,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=161,170 ) /
     :         -1,   -2,    0,    2,    0,
     :          2,    1,    0,    0,    0,
     :          0,    0,    4,    0,    2,
     :          0,    0,    0,    0,    3,
     :          0,    3,    0,    0,    0,
     :          0,    0,    2,   -4,    1,
     :          0,   -1,    0,    2,    1,
     :          0,    0,    0,    4,    1,
     :         -1,   -1,    2,    4,    2,
     :          1,    0,    2,    4,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=171,180 ) /
     :         -2,    2,    0,    2,    0,
     :         -2,   -1,    2,    0,    1,
     :         -2,    0,    0,    2,    2,
     :         -1,   -1,    2,    0,    2,
     :          0,    0,    4,   -2,    1,
     :          3,    0,    2,   -2,    1,
     :         -2,   -1,    0,    2,    1,
     :          1,    0,    0,   -1,    1,
     :          0,   -2,    0,    2,    0,
     :         -2,    0,    0,    4,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=181,190 ) /
     :         -3,    0,    0,    0,    1,
     :          1,    1,    2,    2,    2,
     :          0,    0,    2,    4,    1,
     :          3,    0,    2,    2,    2,
     :         -1,    1,    2,   -2,    1,
     :          2,    0,    0,   -4,    1,
     :          0,    0,    0,   -2,    2,
     :          2,    0,    2,   -4,    1,
     :         -1,    1,    0,    2,    1,
     :          0,    0,    2,   -1,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=191,200 ) /
     :          0,   -2,    2,    2,    2,
     :          2,    0,    0,    2,    1,
     :          4,    0,    2,   -2,    2,
     :          2,    0,    0,   -2,    2,
     :          0,    2,    0,    0,    1,
     :          1,    0,    0,   -4,    1,
     :          0,    2,    2,   -2,    1,
     :         -3,    0,    0,    4,    0,
     :         -1,    1,    2,    0,    1,
     :         -1,   -1,    0,    4,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=201,210 ) /
     :         -1,   -2,    2,    2,    2,
     :         -2,   -1,    2,    4,    2,
     :          1,   -1,    2,    2,    1,
     :         -2,    1,    0,    2,    0,
     :         -2,    1,    2,    0,    1,
     :          2,    1,    0,   -2,    1,
     :         -3,    0,    2,    0,    1,
     :         -2,    0,    2,   -2,    1,
     :         -1,    1,    0,    2,    2,
     :          0,   -1,    2,   -1,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=211,220 ) /
     :         -1,    0,    4,   -2,    2,
     :          0,   -2,    2,    0,    2,
     :         -1,    0,    2,    1,    2,
     :          2,    0,    0,    0,    2,
     :          0,    0,    2,    0,    3,
     :         -2,    0,    4,    0,    2,
     :         -1,    0,   -2,    0,    1,
     :         -1,    1,    2,    2,    1,
     :          3,    0,    0,    0,    1,
     :         -1,    0,    2,    3,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=221,230 ) /
     :          2,   -1,    2,    0,    1,
     :          0,    1,    2,    2,    1,
     :          0,   -1,    2,    4,    2,
     :          2,   -1,    2,    2,    2,
     :          0,    2,   -2,    2,    0,
     :         -1,   -1,    2,   -1,    1,
     :          0,   -2,    0,    0,    1,
     :          1,    0,    2,   -4,    2,
     :          1,   -1,    0,   -2,    1,
     :         -1,   -1,    2,    0,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=231,240 ) /
     :          1,   -1,    2,   -2,    2,
     :         -2,   -1,    0,    4,    0,
     :         -1,    0,    0,    3,    0,
     :         -2,   -1,    2,    2,    2,
     :          0,    2,    2,    0,    2,
     :          1,    1,    0,    2,    0,
     :          2,    0,    2,   -1,    2,
     :          1,    0,    2,    1,    1,
     :          4,    0,    0,    0,    0,
     :          2,    1,    2,    0,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=241,250 ) /
     :          3,   -1,    2,    0,    2,
     :         -2,    2,    0,    2,    1,
     :          1,    0,    2,   -3,    1,
     :          1,    1,    2,   -4,    1,
     :         -1,   -1,    2,   -2,    1,
     :          0,   -1,    0,   -1,    1,
     :          0,   -1,    0,   -2,    1,
     :         -2,    0,    0,    0,    2,
     :         -2,    0,   -2,    2,    0,
     :         -1,    0,   -2,    4,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=251,260 ) /
     :          1,   -2,    0,    0,    0,
     :          0,    1,    0,    1,    1,
     :         -1,    2,    0,    2,    0,
     :          1,   -1,    2,   -2,    1,
     :          1,    2,    2,   -2,    2,
     :          2,   -1,    2,   -2,    2,
     :          1,    0,    2,   -1,    1,
     :          2,    1,    2,   -2,    1,
     :         -2,    0,    0,   -2,    1,
     :          1,   -2,    2,    0,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=261,270 ) /
     :          0,    1,    2,    1,    1,
     :          1,    0,    4,   -2,    1,
     :         -2,    0,    4,    2,    2,
     :          1,    1,    2,    1,    2,
     :          1,    0,    0,    4,    0,
     :          1,    0,    2,    2,    0,
     :          2,    0,    2,    1,    2,
     :          3,    1,    2,    0,    2,
     :          4,    0,    2,    0,    1,
     :         -2,   -1,    2,    0,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=271,280 ) /
     :          0,    1,   -2,    2,    1,
     :          1,    0,   -2,    1,    0,
     :          0,   -1,   -2,    2,    1,
     :          2,   -1,    0,   -2,    1,
     :         -1,    0,    2,   -1,    2,
     :          1,    0,    2,   -3,    2,
     :          0,    1,    2,   -2,    3,
     :          0,    0,    2,   -3,    1,
     :         -1,    0,   -2,    2,    1,
     :          0,    0,    2,   -4,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=281,290 ) /
     :         -2,    1,    0,    0,    1,
     :         -1,    0,    0,   -1,    1,
     :          2,    0,    2,   -4,    2,
     :          0,    0,    4,   -4,    4,
     :          0,    0,    4,   -4,    2,
     :         -1,   -2,    0,    2,    1,
     :         -2,    0,    0,    3,    0,
     :          1,    0,   -2,    2,    1,
     :         -3,    0,    2,    2,    2,
     :         -3,    0,    2,    2,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=291,300 ) /
     :         -2,    0,    2,    2,    0,
     :          2,   -1,    0,    0,    1,
     :         -2,    1,    2,    2,    2,
     :          1,    1,    0,    1,    0,
     :          0,    1,    4,   -2,    2,
     :         -1,    1,    0,   -2,    1,
     :          0,    0,    0,   -4,    1,
     :          1,   -1,    0,    2,    1,
     :          1,    1,    0,    2,    1,
     :         -1,    2,    2,    2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=301,310 ) /
     :          3,    1,    2,   -2,    2,
     :          0,   -1,    0,    4,    0,
     :          2,   -1,    0,    2,    0,
     :          0,    0,    4,    0,    1,
     :          2,    0,    4,   -2,    2,
     :         -1,   -1,    2,    4,    1,
     :          1,    0,    0,    4,    1,
     :          1,   -2,    2,    2,    2,
     :          0,    0,    2,    3,    2,
     :         -1,    1,    2,    4,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=311,320 ) /
     :          3,    0,    0,    2,    0,
     :         -1,    0,    4,    2,    2,
     :          1,    1,    2,    2,    1,
     :         -2,    0,    2,    6,    2,
     :          2,    1,    2,    2,    2,
     :         -1,    0,    2,    6,    2,
     :          1,    0,    2,    4,    1,
     :          2,    0,    2,    4,    2,
     :          1,    1,   -2,    1,    0,
     :         -3,    1,    2,    1,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=321,330 ) /
     :          2,    0,   -2,    0,    2,
     :         -1,    0,    0,    1,    2,
     :         -4,    0,    2,    2,    1,
     :         -1,   -1,    0,    1,    0,
     :          0,    0,   -2,    2,    2,
     :          1,    0,    0,   -1,    2,
     :          0,   -1,    2,   -2,    3,
     :         -2,    1,    2,    0,    0,
     :          0,    0,    2,   -2,    4,
     :         -2,   -2,    0,    2,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=331,340 ) /
     :         -2,    0,   -2,    4,    0,
     :          0,   -2,   -2,    2,    0,
     :          1,    2,    0,   -2,    1,
     :          3,    0,    0,   -4,    1,
     :         -1,    1,    2,   -2,    2,
     :          1,   -1,    2,   -4,    1,
     :          1,    1,    0,   -2,    2,
     :         -3,    0,    2,    0,    0,
     :         -3,    0,    2,    0,    2,
     :         -2,    0,    0,    1,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=341,350 ) /
     :          0,    0,   -2,    1,    0,
     :         -3,    0,    0,    2,    1,
     :         -1,   -1,   -2,    2,    0,
     :          0,    1,    2,   -4,    1,
     :          2,    1,    0,   -4,    1,
     :          0,    2,    0,   -2,    1,
     :          1,    0,    0,   -3,    1,
     :         -2,    0,    2,   -2,    2,
     :         -2,   -1,    0,    0,    1,
     :         -4,    0,    0,    2,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=351,360 ) /
     :          1,    1,    0,   -4,    1,
     :         -1,    0,    2,   -4,    1,
     :          0,    0,    4,   -4,    1,
     :          0,    3,    2,   -2,    2,
     :         -3,   -1,    0,    4,    0,
     :         -3,    0,    0,    4,    1,
     :          1,   -1,   -2,    2,    0,
     :         -1,   -1,    0,    2,    2,
     :          1,   -2,    0,    0,    1,
     :          1,   -1,    0,    0,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=361,370 ) /
     :          0,    0,    0,    1,    2,
     :         -1,   -1,    2,    0,    0,
     :          1,   -2,    2,   -2,    2,
     :          0,   -1,    2,   -1,    1,
     :         -1,    0,    2,    0,    3,
     :          1,    1,    0,    0,    2,
     :         -1,    1,    2,    0,    0,
     :          1,    2,    0,    0,    0,
     :         -1,    2,    2,    0,    2,
     :         -1,    0,    4,   -2,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=371,380 ) /
     :          3,    0,    2,   -4,    2,
     :          1,    2,    2,   -2,    1,
     :          1,    0,    4,   -4,    2,
     :         -2,   -1,    0,    4,    1,
     :          0,   -1,    0,    2,    2,
     :         -2,    1,    0,    4,    0,
     :         -2,   -1,    2,    2,    1,
     :          2,    0,   -2,    2,    0,
     :          1,    0,    0,    1,    1,
     :          0,    1,    0,    2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=381,390 ) /
     :          1,   -1,    2,   -1,    2,
     :         -2,    0,    4,    0,    1,
     :          2,    1,    0,    0,    1,
     :          0,    1,    2,    0,    0,
     :          0,   -1,    4,   -2,    2,
     :          0,    0,    4,   -2,    4,
     :          0,    2,    2,    0,    1,
     :         -3,    0,    0,    6,    0,
     :         -1,   -1,    0,    4,    1,
     :          1,   -2,    0,    2,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=391,400 ) /
     :         -1,    0,    0,    4,    2,
     :         -1,   -2,    2,    2,    1,
     :         -1,    0,    0,   -2,    2,
     :          1,    0,   -2,   -2,    1,
     :          0,    0,   -2,   -2,    1,
     :         -2,    0,   -2,    0,    1,
     :          0,    0,    0,    3,    1,
     :          0,    0,    0,    3,    0,
     :         -1,    1,    0,    4,    0,
     :         -1,   -1,    2,    2,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=401,410 ) /
     :         -2,    0,    2,    3,    2,
     :          1,    0,    0,    2,    2,
     :          0,   -1,    2,    1,    2,
     :          3,   -1,    0,    0,    0,
     :          2,    0,    0,    1,    0,
     :          1,   -1,    2,    0,    0,
     :          0,    0,    2,    1,    0,
     :          1,    0,    2,    0,    3,
     :          3,    1,    0,    0,    0,
     :          3,   -1,    2,   -2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=411,420 ) /
     :          2,    0,    2,   -1,    1,
     :          1,    1,    2,    0,    0,
     :          0,    0,    4,   -1,    2,
     :          1,    2,    2,    0,    2,
     :         -2,    0,    0,    6,    0,
     :          0,   -1,    0,    4,    1,
     :         -2,   -1,    2,    4,    1,
     :          0,   -2,    2,    2,    1,
     :          0,   -1,    2,    2,    0,
     :         -1,    0,    2,    3,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=421,430 ) /
     :         -2,    1,    2,    4,    2,
     :          2,    0,    0,    2,    2,
     :          2,   -2,    2,    0,    2,
     :         -1,    1,    2,    3,    2,
     :          3,    0,    2,   -1,    2,
     :          4,    0,    2,   -2,    1,
     :         -1,    0,    0,    6,    0,
     :         -1,   -2,    2,    4,    2,
     :         -3,    0,    2,    6,    2,
     :         -1,    0,    2,    4,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=431,440 ) /
     :          3,    0,    0,    2,    1,
     :          3,   -1,    2,    0,    1,
     :          3,    0,    2,    0,    0,
     :          1,    0,    4,    0,    2,
     :          5,    0,    2,   -2,    2,
     :          0,   -1,    2,    4,    1,
     :          2,   -1,    2,    2,    1,
     :          0,    1,    2,    4,    2,
     :          1,   -1,    2,    4,    2,
     :          3,   -1,    2,    2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=441,450 ) /
     :          3,    0,    2,    2,    1,
     :          5,    0,    2,    0,    2,
     :          0,    0,    2,    6,    2,
     :          4,    0,    2,    2,    2,
     :          0,   -1,    1,   -1,    1,
     :         -1,    0,    1,    0,    3,
     :          0,   -2,    2,   -2,    3,
     :          1,    0,   -1,    0,    1,
     :          2,   -2,    0,   -2,    1,
     :         -1,    0,    1,    0,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=451,460 ) /
     :         -1,    0,    1,    0,    1,
     :         -1,   -1,    2,   -1,    2,
     :         -2,    2,    0,    2,    2,
     :         -1,    0,    1,    0,    0,
     :         -4,    1,    2,    2,    2,
     :         -3,    0,    2,    1,    1,
     :         -2,   -1,    2,    0,    2,
     :          1,    0,   -2,    1,    1,
     :          2,   -1,   -2,    0,    1,
     :         -4,    0,    2,    2,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=461,470 ) /
     :         -3,    1,    0,    3,    0,
     :         -1,    0,   -1,    2,    0,
     :          0,   -2,    0,    0,    2,
     :          0,   -2,    0,    0,    2,
     :         -3,    0,    0,    3,    0,
     :         -2,   -1,    0,    2,    2,
     :         -1,    0,   -2,    3,    0,
     :         -4,    0,    0,    4,    0,
     :          2,    1,   -2,    0,    1,
     :          2,   -1,    0,   -2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=471,480 ) /
     :          0,    0,    1,   -1,    0,
     :         -1,    2,    0,    1,    0,
     :         -2,    1,    2,    0,    2,
     :          1,    1,    0,   -1,    1,
     :          1,    0,    1,   -2,    1,
     :          0,    2,    0,    0,    2,
     :          1,   -1,    2,   -3,    1,
     :         -1,    1,    2,   -1,    1,
     :         -2,    0,    4,   -2,    2,
     :         -2,    0,    4,   -2,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=481,490 ) /
     :         -2,   -2,    0,    2,    1,
     :         -2,    0,   -2,    4,    0,
     :          1,    2,    2,   -4,    1,
     :          1,    1,    2,   -4,    2,
     :         -1,    2,    2,   -2,    1,
     :          2,    0,    0,   -3,    1,
     :         -1,    2,    0,    0,    1,
     :          0,    0,    0,   -2,    0,
     :         -1,   -1,    2,   -2,    2,
     :         -1,    1,    0,    0,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=491,500 ) /
     :          0,    0,    0,   -1,    2,
     :         -2,    1,    0,    1,    0,
     :          1,   -2,    0,   -2,    1,
     :          1,    0,   -2,    0,    2,
     :         -3,    1,    0,    2,    0,
     :         -1,    1,   -2,    2,    0,
     :         -1,   -1,    0,    0,    2,
     :         -3,    0,    0,    2,    0,
     :         -3,   -1,    0,    2,    0,
     :          2,    0,    2,   -6,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=501,510 ) /
     :          0,    1,    2,   -4,    2,
     :          2,    0,    0,   -4,    2,
     :         -2,    1,    2,   -2,    1,
     :          0,   -1,    2,   -4,    1,
     :          0,    1,    0,   -2,    2,
     :         -1,    0,    0,   -2,    0,
     :          2,    0,   -2,   -2,    1,
     :         -4,    0,    2,    0,    1,
     :         -1,   -1,    0,   -1,    1,
     :          0,    0,   -2,    0,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=511,520 ) /
     :         -3,    0,    0,    1,    0,
     :         -1,    0,   -2,    1,    0,
     :         -2,    0,   -2,    2,    1,
     :          0,    0,   -4,    2,    0,
     :         -2,   -1,   -2,    2,    0,
     :          1,    0,    2,   -6,    1,
     :         -1,    0,    2,   -4,    2,
     :          1,    0,    0,   -4,    2,
     :          2,    1,    2,   -4,    2,
     :          2,    1,    2,   -4,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=521,530 ) /
     :          0,    1,    4,   -4,    4,
     :          0,    1,    4,   -4,    2,
     :         -1,   -1,   -2,    4,    0,
     :         -1,   -3,    0,    2,    0,
     :         -1,    0,   -2,    4,    1,
     :         -2,   -1,    0,    3,    0,
     :          0,    0,   -2,    3,    0,
     :         -2,    0,    0,    3,    1,
     :          0,   -1,    0,    1,    0,
     :         -3,    0,    2,    2,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=531,540 ) /
     :          1,    1,   -2,    2,    0,
     :         -1,    1,    0,    2,    2,
     :          1,   -2,    2,   -2,    1,
     :          0,    0,    1,    0,    2,
     :          0,    0,    1,    0,    1,
     :          0,    0,    1,    0,    0,
     :         -1,    2,    0,    2,    1,
     :          0,    0,    2,    0,    2,
     :         -2,    0,    2,    0,    2,
     :          2,    0,    0,   -1,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=541,550 ) /
     :          3,    0,    0,   -2,    1,
     :          1,    0,    2,   -2,    3,
     :          1,    2,    0,    0,    1,
     :          2,    0,    2,   -3,    2,
     :         -1,    1,    4,   -2,    2,
     :         -2,   -2,    0,    4,    0,
     :          0,   -3,    0,    2,    0,
     :          0,    0,   -2,    4,    0,
     :         -1,   -1,    0,    3,    0,
     :         -2,    0,    0,    4,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=551,560 ) /
     :         -1,    0,    0,    3,    1,
     :          2,   -2,    0,    0,    0,
     :          1,   -1,    0,    1,    0,
     :         -1,    0,    0,    2,    0,
     :          0,   -2,    2,    0,    1,
     :         -1,    0,    1,    2,    1,
     :         -1,    1,    0,    3,    0,
     :         -1,   -1,    2,    1,    2,
     :          0,   -1,    2,    0,    0,
     :         -2,    1,    2,    2,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=561,570 ) /
     :          2,   -2,    2,   -2,    2,
     :          1,    1,    0,    1,    1,
     :          1,    0,    1,    0,    1,
     :          1,    0,    1,    0,    0,
     :          0,    2,    0,    2,    0,
     :          2,   -1,    2,   -2,    1,
     :          0,   -1,    4,   -2,    1,
     :          0,    0,    4,   -2,    3,
     :          0,    1,    4,   -2,    1,
     :          4,    0,    2,   -4,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=571,580 ) /
     :          2,    2,    2,   -2,    2,
     :          2,    0,    4,   -4,    2,
     :         -1,   -2,    0,    4,    0,
     :         -1,   -3,    2,    2,    2,
     :         -3,    0,    2,    4,    2,
     :         -3,    0,    2,   -2,    1,
     :         -1,   -1,    0,   -2,    1,
     :         -3,    0,    0,    0,    2,
     :         -3,    0,   -2,    2,    0,
     :          0,    1,    0,   -4,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=581,590 ) /
     :         -2,    1,    0,   -2,    1,
     :         -4,    0,    0,    0,    1,
     :         -1,    0,    0,   -4,    1,
     :         -3,    0,    0,   -2,    1,
     :          0,    0,    0,    3,    2,
     :         -1,    1,    0,    4,    1,
     :          1,   -2,    2,    0,    1,
     :          0,    1,    0,    3,    0,
     :         -1,    0,    2,    2,    3,
     :          0,    0,    2,    2,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=591,600 ) /
     :         -2,    0,    2,    2,    2,
     :         -1,    1,    2,    2,    0,
     :          3,    0,    0,    0,    2,
     :          2,    1,    0,    1,    0,
     :          2,   -1,    2,   -1,    2,
     :          0,    0,    2,    0,    1,
     :          0,    0,    3,    0,    3,
     :          0,    0,    3,    0,    2,
     :         -1,    2,    2,    2,    1,
     :         -1,    0,    4,    0,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=601,610 ) /
     :          1,    2,    2,    0,    1,
     :          3,    1,    2,   -2,    1,
     :          1,    1,    4,   -2,    2,
     :         -2,   -1,    0,    6,    0,
     :          0,   -2,    0,    4,    0,
     :         -2,    0,    0,    6,    1,
     :         -2,   -2,    2,    4,    2,
     :          0,   -3,    2,    2,    2,
     :          0,    0,    0,    4,    2,
     :         -1,   -1,    2,    3,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=611,620 ) /
     :         -2,    0,    2,    4,    0,
     :          2,   -1,    0,    2,    1,
     :          1,    0,    0,    3,    0,
     :          0,    1,    0,    4,    1,
     :          0,    1,    0,    4,    0,
     :          1,   -1,    2,    1,    2,
     :          0,    0,    2,    2,    3,
     :          1,    0,    2,    2,    2,
     :         -1,    0,    2,    2,    2,
     :         -2,    0,    4,    2,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=621,630 ) /
     :          2,    1,    0,    2,    1,
     :          2,    1,    0,    2,    0,
     :          2,   -1,    2,    0,    0,
     :          1,    0,    2,    1,    0,
     :          0,    1,    2,    2,    0,
     :          2,    0,    2,    0,    3,
     :          3,    0,    2,    0,    2,
     :          1,    0,    2,    0,    2,
     :          1,    0,    3,    0,    3,
     :          1,    1,    2,    1,    1 /
      DATA ( ( NALS(I,J), I=1,5 ), J=631,640 ) /
     :          0,    2,    2,    2,    2,
     :          2,    1,    2,    0,    0,
     :          2,    0,    4,   -2,    1,
     :          4,    1,    2,   -2,    2,
     :         -1,   -1,    0,    6,    0,
     :         -3,   -1,    2,    6,    2,
     :         -1,    0,    0,    6,    1,
     :         -3,    0,    2,    6,    1,
     :          1,   -1,    0,    4,    1,
     :          1,   -1,    0,    4,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=641,650 ) /
     :         -2,    0,    2,    5,    2,
     :          1,   -2,    2,    2,    1,
     :          3,   -1,    0,    2,    0,
     :          1,   -1,    2,    2,    0,
     :          0,    0,    2,    3,    1,
     :         -1,    1,    2,    4,    1,
     :          0,    1,    2,    3,    2,
     :         -1,    0,    4,    2,    1,
     :          2,    0,    2,    1,    1,
     :          5,    0,    0,    0,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=651,660 ) /
     :          2,    1,    2,    1,    2,
     :          1,    0,    4,    0,    1,
     :          3,    1,    2,    0,    1,
     :          3,    0,    4,   -2,    2,
     :         -2,   -1,    2,    6,    2,
     :          0,    0,    0,    6,    0,
     :          0,   -2,    2,    4,    2,
     :         -2,    0,    2,    6,    1,
     :          2,    0,    0,    4,    1,
     :          2,    0,    0,    4,    0 /
      DATA ( ( NALS(I,J), I=1,5 ), J=661,670 ) /
     :          2,   -2,    2,    2,    2,
     :          0,    0,    2,    4,    0,
     :          1,    0,    2,    3,    2,
     :          4,    0,    0,    2,    0,
     :          2,    0,    2,    2,    0,
     :          0,    0,    4,    2,    2,
     :          4,   -1,    2,    0,    2,
     :          3,    0,    2,    1,    2,
     :          2,    1,    2,    2,    1,
     :          4,    1,    2,    0,    2 /
      DATA ( ( NALS(I,J), I=1,5 ), J=671,678 ) /
     :         -1,   -1,    2,    6,    2,
     :         -1,    0,    2,    6,    1,
     :          1,   -1,    2,    4,    1,
     :          1,    1,    2,    4,    2,
     :          3,    1,    2,    2,    2,
     :          5,    0,    2,    0,    1,
     :          2,   -1,    2,    4,    2,
     :          2,    0,    2,    4,    1 /

*
*  Luni-Solar nutation coefficients, unit 1e-7 arcsec:
*  longitude (sin, t*sin, cos), obliquity (cos, t*cos, sin)
*

      DATA ( ( CLS(I,J), I=1,6 ), J=  1, 10 ) /
     : -172064161D0, -174666D0,  33386D0, 92052331D0,  9086D0, 15377D0,
     :  -13170906D0,   -1675D0, -13696D0,  5730336D0, -3015D0, -4587D0,
     :   -2276413D0,    -234D0,   2796D0,   978459D0,  -485D0,  1374D0,
     :    2074554D0,     207D0,   -698D0,  -897492D0,   470D0,  -291D0,
     :    1475877D0,   -3633D0,  11817D0,    73871D0,  -184D0, -1924D0,
     :    -516821D0,    1226D0,   -524D0,   224386D0,  -677D0,  -174D0,
     :     711159D0,      73D0,   -872D0,    -6750D0,     0D0,   358D0,
     :    -387298D0,    -367D0,    380D0,   200728D0,    18D0,   318D0,
     :    -301461D0,     -36D0,    816D0,   129025D0,   -63D0,   367D0,
     :     215829D0,    -494D0,    111D0,   -95929D0,   299D0,   132D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J= 11, 20 ) /
     :     128227D0,     137D0,    181D0,   -68982D0,    -9D0,    39D0,
     :     123457D0,      11D0,     19D0,   -53311D0,    32D0,    -4D0,
     :     156994D0,      10D0,   -168D0,    -1235D0,     0D0,    82D0,
     :      63110D0,      63D0,     27D0,   -33228D0,     0D0,    -9D0,
     :     -57976D0,     -63D0,   -189D0,    31429D0,     0D0,   -75D0,
     :     -59641D0,     -11D0,    149D0,    25543D0,   -11D0,    66D0,
     :     -51613D0,     -42D0,    129D0,    26366D0,     0D0,    78D0,
     :      45893D0,      50D0,     31D0,   -24236D0,   -10D0,    20D0,
     :      63384D0,      11D0,   -150D0,    -1220D0,     0D0,    29D0,
     :     -38571D0,      -1D0,    158D0,    16452D0,   -11D0,    68D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J= 21, 30 ) /
     :      32481D0,       0D0,      0D0,   -13870D0,     0D0,     0D0,
     :     -47722D0,       0D0,    -18D0,      477D0,     0D0,   -25D0,
     :     -31046D0,      -1D0,    131D0,    13238D0,   -11D0,    59D0,
     :      28593D0,       0D0,     -1D0,   -12338D0,    10D0,    -3D0,
     :      20441D0,      21D0,     10D0,   -10758D0,     0D0,    -3D0,
     :      29243D0,       0D0,    -74D0,     -609D0,     0D0,    13D0,
     :      25887D0,       0D0,    -66D0,     -550D0,     0D0,    11D0,
     :     -14053D0,     -25D0,     79D0,     8551D0,    -2D0,   -45D0,
     :      15164D0,      10D0,     11D0,    -8001D0,     0D0,    -1D0,
     :     -15794D0,      72D0,    -16D0,     6850D0,   -42D0,    -5D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J= 31, 40 ) /
     :      21783D0,       0D0,     13D0,     -167D0,     0D0,    13D0,
     :     -12873D0,     -10D0,    -37D0,     6953D0,     0D0,   -14D0,
     :     -12654D0,      11D0,     63D0,     6415D0,     0D0,    26D0,
     :     -10204D0,       0D0,     25D0,     5222D0,     0D0,    15D0,
     :      16707D0,     -85D0,    -10D0,      168D0,    -1D0,    10D0,
     :      -7691D0,       0D0,     44D0,     3268D0,     0D0,    19D0,
     :     -11024D0,       0D0,    -14D0,      104D0,     0D0,     2D0,
     :       7566D0,     -21D0,    -11D0,    -3250D0,     0D0,    -5D0,
     :      -6637D0,     -11D0,     25D0,     3353D0,     0D0,    14D0,
     :      -7141D0,      21D0,      8D0,     3070D0,     0D0,     4D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J= 41, 50 ) /
     :      -6302D0,     -11D0,      2D0,     3272D0,     0D0,     4D0,
     :       5800D0,      10D0,      2D0,    -3045D0,     0D0,    -1D0,
     :       6443D0,       0D0,     -7D0,    -2768D0,     0D0,    -4D0,
     :      -5774D0,     -11D0,    -15D0,     3041D0,     0D0,    -5D0,
     :      -5350D0,       0D0,     21D0,     2695D0,     0D0,    12D0,
     :      -4752D0,     -11D0,     -3D0,     2719D0,     0D0,    -3D0,
     :      -4940D0,     -11D0,    -21D0,     2720D0,     0D0,    -9D0,
     :       7350D0,       0D0,     -8D0,      -51D0,     0D0,     4D0,
     :       4065D0,       0D0,      6D0,    -2206D0,     0D0,     1D0,
     :       6579D0,       0D0,    -24D0,     -199D0,     0D0,     2D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J= 51, 60 ) /
     :       3579D0,       0D0,      5D0,    -1900D0,     0D0,     1D0,
     :       4725D0,       0D0,     -6D0,      -41D0,     0D0,     3D0,
     :      -3075D0,       0D0,     -2D0,     1313D0,     0D0,    -1D0,
     :      -2904D0,       0D0,     15D0,     1233D0,     0D0,     7D0,
     :       4348D0,       0D0,    -10D0,      -81D0,     0D0,     2D0,
     :      -2878D0,       0D0,      8D0,     1232D0,     0D0,     4D0,
     :      -4230D0,       0D0,      5D0,      -20D0,     0D0,    -2D0,
     :      -2819D0,       0D0,      7D0,     1207D0,     0D0,     3D0,
     :      -4056D0,       0D0,      5D0,       40D0,     0D0,    -2D0,
     :      -2647D0,       0D0,     11D0,     1129D0,     0D0,     5D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J= 61, 70 ) /
     :      -2294D0,       0D0,    -10D0,     1266D0,     0D0,    -4D0,
     :       2481D0,       0D0,     -7D0,    -1062D0,     0D0,    -3D0,
     :       2179D0,       0D0,     -2D0,    -1129D0,     0D0,    -2D0,
     :       3276D0,       0D0,      1D0,       -9D0,     0D0,     0D0,
     :      -3389D0,       0D0,      5D0,       35D0,     0D0,    -2D0,
     :       3339D0,       0D0,    -13D0,     -107D0,     0D0,     1D0,
     :      -1987D0,       0D0,     -6D0,     1073D0,     0D0,    -2D0,
     :      -1981D0,       0D0,      0D0,      854D0,     0D0,     0D0,
     :       4026D0,       0D0,   -353D0,     -553D0,     0D0,  -139D0,
     :       1660D0,       0D0,     -5D0,     -710D0,     0D0,    -2D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J= 71, 80 ) /
     :      -1521D0,       0D0,      9D0,      647D0,     0D0,     4D0,
     :       1314D0,       0D0,      0D0,     -700D0,     0D0,     0D0,
     :      -1283D0,       0D0,      0D0,      672D0,     0D0,     0D0,
     :      -1331D0,       0D0,      8D0,      663D0,     0D0,     4D0,
     :       1383D0,       0D0,     -2D0,     -594D0,     0D0,    -2D0,
     :       1405D0,       0D0,      4D0,     -610D0,     0D0,     2D0,
     :       1290D0,       0D0,      0D0,     -556D0,     0D0,     0D0,
     :      -1214D0,       0D0,      5D0,      518D0,     0D0,     2D0,
     :       1146D0,       0D0,     -3D0,     -490D0,     0D0,    -1D0,
     :       1019D0,       0D0,     -1D0,     -527D0,     0D0,    -1D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J= 81, 90 ) /
     :      -1100D0,       0D0,      9D0,      465D0,     0D0,     4D0,
     :       -970D0,       0D0,      2D0,      496D0,     0D0,     1D0,
     :       1575D0,       0D0,     -6D0,      -50D0,     0D0,     0D0,
     :        934D0,       0D0,     -3D0,     -399D0,     0D0,    -1D0,
     :        922D0,       0D0,     -1D0,     -395D0,     0D0,    -1D0,
     :        815D0,       0D0,     -1D0,     -422D0,     0D0,    -1D0,
     :        834D0,       0D0,      2D0,     -440D0,     0D0,     1D0,
     :       1248D0,       0D0,      0D0,     -170D0,     0D0,     1D0,
     :       1338D0,       0D0,     -5D0,      -39D0,     0D0,     0D0,
     :        716D0,       0D0,     -2D0,     -389D0,     0D0,    -1D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J= 91,100 ) /
     :       1282D0,       0D0,     -3D0,      -23D0,     0D0,     1D0,
     :        742D0,       0D0,      1D0,     -391D0,     0D0,     0D0,
     :       1020D0,       0D0,    -25D0,     -495D0,     0D0,   -10D0,
     :        715D0,       0D0,     -4D0,     -326D0,     0D0,     2D0,
     :       -666D0,       0D0,     -3D0,      369D0,     0D0,    -1D0,
     :       -667D0,       0D0,      1D0,      346D0,     0D0,     1D0,
     :       -704D0,       0D0,      0D0,      304D0,     0D0,     0D0,
     :       -694D0,       0D0,      5D0,      294D0,     0D0,     2D0,
     :      -1014D0,       0D0,     -1D0,        4D0,     0D0,    -1D0,
     :       -585D0,       0D0,     -2D0,      316D0,     0D0,    -1D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=101,110 ) /
     :       -949D0,       0D0,      1D0,        8D0,     0D0,    -1D0,
     :       -595D0,       0D0,      0D0,      258D0,     0D0,     0D0,
     :        528D0,       0D0,      0D0,     -279D0,     0D0,     0D0,
     :       -590D0,       0D0,      4D0,      252D0,     0D0,     2D0,
     :        570D0,       0D0,     -2D0,     -244D0,     0D0,    -1D0,
     :       -502D0,       0D0,      3D0,      250D0,     0D0,     2D0,
     :       -875D0,       0D0,      1D0,       29D0,     0D0,     0D0,
     :       -492D0,       0D0,     -3D0,      275D0,     0D0,    -1D0,
     :        535D0,       0D0,     -2D0,     -228D0,     0D0,    -1D0,
     :       -467D0,       0D0,      1D0,      240D0,     0D0,     1D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=111,120 ) /
     :        591D0,       0D0,      0D0,     -253D0,     0D0,     0D0,
     :       -453D0,       0D0,     -1D0,      244D0,     0D0,    -1D0,
     :        766D0,       0D0,      1D0,        9D0,     0D0,     0D0,
     :       -446D0,       0D0,      2D0,      225D0,     0D0,     1D0,
     :       -488D0,       0D0,      2D0,      207D0,     0D0,     1D0,
     :       -468D0,       0D0,      0D0,      201D0,     0D0,     0D0,
     :       -421D0,       0D0,      1D0,      216D0,     0D0,     1D0,
     :        463D0,       0D0,      0D0,     -200D0,     0D0,     0D0,
     :       -673D0,       0D0,      2D0,       14D0,     0D0,     0D0,
     :        658D0,       0D0,      0D0,       -2D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=121,130 ) /
     :       -438D0,       0D0,      0D0,      188D0,     0D0,     0D0,
     :       -390D0,       0D0,      0D0,      205D0,     0D0,     0D0,
     :        639D0,     -11D0,     -2D0,      -19D0,     0D0,     0D0,
     :        412D0,       0D0,     -2D0,     -176D0,     0D0,    -1D0,
     :       -361D0,       0D0,      0D0,      189D0,     0D0,     0D0,
     :        360D0,       0D0,     -1D0,     -185D0,     0D0,    -1D0,
     :        588D0,       0D0,     -3D0,      -24D0,     0D0,     0D0,
     :       -578D0,       0D0,      1D0,        5D0,     0D0,     0D0,
     :       -396D0,       0D0,      0D0,      171D0,     0D0,     0D0,
     :        565D0,       0D0,     -1D0,       -6D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=131,140 ) /
     :       -335D0,       0D0,     -1D0,      184D0,     0D0,    -1D0,
     :        357D0,       0D0,      1D0,     -154D0,     0D0,     0D0,
     :        321D0,       0D0,      1D0,     -174D0,     0D0,     0D0,
     :       -301D0,       0D0,     -1D0,      162D0,     0D0,     0D0,
     :       -334D0,       0D0,      0D0,      144D0,     0D0,     0D0,
     :        493D0,       0D0,     -2D0,      -15D0,     0D0,     0D0,
     :        494D0,       0D0,     -2D0,      -19D0,     0D0,     0D0,
     :        337D0,       0D0,     -1D0,     -143D0,     0D0,    -1D0,
     :        280D0,       0D0,     -1D0,     -144D0,     0D0,     0D0,
     :        309D0,       0D0,      1D0,     -134D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=141,150 ) /
     :       -263D0,       0D0,      2D0,      131D0,     0D0,     1D0,
     :        253D0,       0D0,      1D0,     -138D0,     0D0,     0D0,
     :        245D0,       0D0,      0D0,     -128D0,     0D0,     0D0,
     :        416D0,       0D0,     -2D0,      -17D0,     0D0,     0D0,
     :       -229D0,       0D0,      0D0,      128D0,     0D0,     0D0,
     :        231D0,       0D0,      0D0,     -120D0,     0D0,     0D0,
     :       -259D0,       0D0,      2D0,      109D0,     0D0,     1D0,
     :        375D0,       0D0,     -1D0,       -8D0,     0D0,     0D0,
     :        252D0,       0D0,      0D0,     -108D0,     0D0,     0D0,
     :       -245D0,       0D0,      1D0,      104D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=151,160 ) /
     :        243D0,       0D0,     -1D0,     -104D0,     0D0,     0D0,
     :        208D0,       0D0,      1D0,     -112D0,     0D0,     0D0,
     :        199D0,       0D0,      0D0,     -102D0,     0D0,     0D0,
     :       -208D0,       0D0,      1D0,      105D0,     0D0,     0D0,
     :        335D0,       0D0,     -2D0,      -14D0,     0D0,     0D0,
     :       -325D0,       0D0,      1D0,        7D0,     0D0,     0D0,
     :       -187D0,       0D0,      0D0,       96D0,     0D0,     0D0,
     :        197D0,       0D0,     -1D0,     -100D0,     0D0,     0D0,
     :       -192D0,       0D0,      2D0,       94D0,     0D0,     1D0,
     :       -188D0,       0D0,      0D0,       83D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=161,170 ) /
     :        276D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :       -286D0,       0D0,      1D0,        6D0,     0D0,     0D0,
     :        186D0,       0D0,     -1D0,      -79D0,     0D0,     0D0,
     :       -219D0,       0D0,      0D0,       43D0,     0D0,     0D0,
     :        276D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :       -153D0,       0D0,     -1D0,       84D0,     0D0,     0D0,
     :       -156D0,       0D0,      0D0,       81D0,     0D0,     0D0,
     :       -154D0,       0D0,      1D0,       78D0,     0D0,     0D0,
     :       -174D0,       0D0,      1D0,       75D0,     0D0,     0D0,
     :       -163D0,       0D0,      2D0,       69D0,     0D0,     1D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=171,180 ) /
     :       -228D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :         91D0,       0D0,     -4D0,      -54D0,     0D0,    -2D0,
     :        175D0,       0D0,      0D0,      -75D0,     0D0,     0D0,
     :       -159D0,       0D0,      0D0,       69D0,     0D0,     0D0,
     :        141D0,       0D0,      0D0,      -72D0,     0D0,     0D0,
     :        147D0,       0D0,      0D0,      -75D0,     0D0,     0D0,
     :       -132D0,       0D0,      0D0,       69D0,     0D0,     0D0,
     :        159D0,       0D0,    -28D0,      -54D0,     0D0,    11D0,
     :        213D0,       0D0,      0D0,       -4D0,     0D0,     0D0,
     :        123D0,       0D0,      0D0,      -64D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=181,190 ) /
     :       -118D0,       0D0,     -1D0,       66D0,     0D0,     0D0,
     :        144D0,       0D0,     -1D0,      -61D0,     0D0,     0D0,
     :       -121D0,       0D0,      1D0,       60D0,     0D0,     0D0,
     :       -134D0,       0D0,      1D0,       56D0,     0D0,     1D0,
     :       -105D0,       0D0,      0D0,       57D0,     0D0,     0D0,
     :       -102D0,       0D0,      0D0,       56D0,     0D0,     0D0,
     :        120D0,       0D0,      0D0,      -52D0,     0D0,     0D0,
     :        101D0,       0D0,      0D0,      -54D0,     0D0,     0D0,
     :       -113D0,       0D0,      0D0,       59D0,     0D0,     0D0,
     :       -106D0,       0D0,      0D0,       61D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=191,200 ) /
     :       -129D0,       0D0,      1D0,       55D0,     0D0,     0D0,
     :       -114D0,       0D0,      0D0,       57D0,     0D0,     0D0,
     :        113D0,       0D0,     -1D0,      -49D0,     0D0,     0D0,
     :       -102D0,       0D0,      0D0,       44D0,     0D0,     0D0,
     :        -94D0,       0D0,      0D0,       51D0,     0D0,     0D0,
     :       -100D0,       0D0,     -1D0,       56D0,     0D0,     0D0,
     :         87D0,       0D0,      0D0,      -47D0,     0D0,     0D0,
     :        161D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         96D0,       0D0,      0D0,      -50D0,     0D0,     0D0,
     :        151D0,       0D0,     -1D0,       -5D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=201,210 ) /
     :       -104D0,       0D0,      0D0,       44D0,     0D0,     0D0,
     :       -110D0,       0D0,      0D0,       48D0,     0D0,     0D0,
     :       -100D0,       0D0,      1D0,       50D0,     0D0,     0D0,
     :         92D0,       0D0,     -5D0,       12D0,     0D0,    -2D0,
     :         82D0,       0D0,      0D0,      -45D0,     0D0,     0D0,
     :         82D0,       0D0,      0D0,      -45D0,     0D0,     0D0,
     :        -78D0,       0D0,      0D0,       41D0,     0D0,     0D0,
     :        -77D0,       0D0,      0D0,       43D0,     0D0,     0D0,
     :          2D0,       0D0,      0D0,       54D0,     0D0,     0D0,
     :         94D0,       0D0,      0D0,      -40D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=211,220 ) /
     :        -93D0,       0D0,      0D0,       40D0,     0D0,     0D0,
     :        -83D0,       0D0,     10D0,       40D0,     0D0,    -2D0,
     :         83D0,       0D0,      0D0,      -36D0,     0D0,     0D0,
     :        -91D0,       0D0,      0D0,       39D0,     0D0,     0D0,
     :        128D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :        -79D0,       0D0,      0D0,       34D0,     0D0,     0D0,
     :        -83D0,       0D0,      0D0,       47D0,     0D0,     0D0,
     :         84D0,       0D0,      0D0,      -44D0,     0D0,     0D0,
     :         83D0,       0D0,      0D0,      -43D0,     0D0,     0D0,
     :         91D0,       0D0,      0D0,      -39D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=221,230 ) /
     :        -77D0,       0D0,      0D0,       39D0,     0D0,     0D0,
     :         84D0,       0D0,      0D0,      -43D0,     0D0,     0D0,
     :        -92D0,       0D0,      1D0,       39D0,     0D0,     0D0,
     :        -92D0,       0D0,      1D0,       39D0,     0D0,     0D0,
     :        -94D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         68D0,       0D0,      0D0,      -36D0,     0D0,     0D0,
     :        -61D0,       0D0,      0D0,       32D0,     0D0,     0D0,
     :         71D0,       0D0,      0D0,      -31D0,     0D0,     0D0,
     :         62D0,       0D0,      0D0,      -34D0,     0D0,     0D0,
     :        -63D0,       0D0,      0D0,       33D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=231,240 ) /
     :        -73D0,       0D0,      0D0,       32D0,     0D0,     0D0,
     :        115D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :       -103D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         63D0,       0D0,      0D0,      -28D0,     0D0,     0D0,
     :         74D0,       0D0,      0D0,      -32D0,     0D0,     0D0,
     :       -103D0,       0D0,     -3D0,        3D0,     0D0,    -1D0,
     :        -69D0,       0D0,      0D0,       30D0,     0D0,     0D0,
     :         57D0,       0D0,      0D0,      -29D0,     0D0,     0D0,
     :         94D0,       0D0,      0D0,       -4D0,     0D0,     0D0,
     :         64D0,       0D0,      0D0,      -33D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=241,250 ) /
     :        -63D0,       0D0,      0D0,       26D0,     0D0,     0D0,
     :        -38D0,       0D0,      0D0,       20D0,     0D0,     0D0,
     :        -43D0,       0D0,      0D0,       24D0,     0D0,     0D0,
     :        -45D0,       0D0,      0D0,       23D0,     0D0,     0D0,
     :         47D0,       0D0,      0D0,      -24D0,     0D0,     0D0,
     :        -48D0,       0D0,      0D0,       25D0,     0D0,     0D0,
     :         45D0,       0D0,      0D0,      -26D0,     0D0,     0D0,
     :         56D0,       0D0,      0D0,      -25D0,     0D0,     0D0,
     :         88D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :        -75D0,       0D0,      0D0,        0D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=251,260 ) /
     :         85D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         49D0,       0D0,      0D0,      -26D0,     0D0,     0D0,
     :        -74D0,       0D0,     -3D0,       -1D0,     0D0,    -1D0,
     :        -39D0,       0D0,      0D0,       21D0,     0D0,     0D0,
     :         45D0,       0D0,      0D0,      -20D0,     0D0,     0D0,
     :         51D0,       0D0,      0D0,      -22D0,     0D0,     0D0,
     :        -40D0,       0D0,      0D0,       21D0,     0D0,     0D0,
     :         41D0,       0D0,      0D0,      -21D0,     0D0,     0D0,
     :        -42D0,       0D0,      0D0,       24D0,     0D0,     0D0,
     :        -51D0,       0D0,      0D0,       22D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=261,270 ) /
     :        -42D0,       0D0,      0D0,       22D0,     0D0,     0D0,
     :         39D0,       0D0,      0D0,      -21D0,     0D0,     0D0,
     :         46D0,       0D0,      0D0,      -18D0,     0D0,     0D0,
     :        -53D0,       0D0,      0D0,       22D0,     0D0,     0D0,
     :         82D0,       0D0,      0D0,       -4D0,     0D0,     0D0,
     :         81D0,       0D0,     -1D0,       -4D0,     0D0,     0D0,
     :         47D0,       0D0,      0D0,      -19D0,     0D0,     0D0,
     :         53D0,       0D0,      0D0,      -23D0,     0D0,     0D0,
     :        -45D0,       0D0,      0D0,       22D0,     0D0,     0D0,
     :        -44D0,       0D0,      0D0,       -2D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=271,280 ) /
     :        -33D0,       0D0,      0D0,       16D0,     0D0,     0D0,
     :        -61D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :         28D0,       0D0,      0D0,      -15D0,     0D0,     0D0,
     :        -38D0,       0D0,      0D0,       19D0,     0D0,     0D0,
     :        -33D0,       0D0,      0D0,       21D0,     0D0,     0D0,
     :        -60D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         48D0,       0D0,      0D0,      -10D0,     0D0,     0D0,
     :         27D0,       0D0,      0D0,      -14D0,     0D0,     0D0,
     :         38D0,       0D0,      0D0,      -20D0,     0D0,     0D0,
     :         31D0,       0D0,      0D0,      -13D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=281,290 ) /
     :        -29D0,       0D0,      0D0,       15D0,     0D0,     0D0,
     :         28D0,       0D0,      0D0,      -15D0,     0D0,     0D0,
     :        -32D0,       0D0,      0D0,       15D0,     0D0,     0D0,
     :         45D0,       0D0,      0D0,       -8D0,     0D0,     0D0,
     :        -44D0,       0D0,      0D0,       19D0,     0D0,     0D0,
     :         28D0,       0D0,      0D0,      -15D0,     0D0,     0D0,
     :        -51D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -36D0,       0D0,      0D0,       20D0,     0D0,     0D0,
     :         44D0,       0D0,      0D0,      -19D0,     0D0,     0D0,
     :         26D0,       0D0,      0D0,      -14D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=291,300 ) /
     :        -60D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         35D0,       0D0,      0D0,      -18D0,     0D0,     0D0,
     :        -27D0,       0D0,      0D0,       11D0,     0D0,     0D0,
     :         47D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         36D0,       0D0,      0D0,      -15D0,     0D0,     0D0,
     :        -36D0,       0D0,      0D0,       20D0,     0D0,     0D0,
     :        -35D0,       0D0,      0D0,       19D0,     0D0,     0D0,
     :        -37D0,       0D0,      0D0,       19D0,     0D0,     0D0,
     :         32D0,       0D0,      0D0,      -16D0,     0D0,     0D0,
     :         35D0,       0D0,      0D0,      -14D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=301,310 ) /
     :         32D0,       0D0,      0D0,      -13D0,     0D0,     0D0,
     :         65D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :         47D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         32D0,       0D0,      0D0,      -16D0,     0D0,     0D0,
     :         37D0,       0D0,      0D0,      -16D0,     0D0,     0D0,
     :        -30D0,       0D0,      0D0,       15D0,     0D0,     0D0,
     :        -32D0,       0D0,      0D0,       16D0,     0D0,     0D0,
     :        -31D0,       0D0,      0D0,       13D0,     0D0,     0D0,
     :         37D0,       0D0,      0D0,      -16D0,     0D0,     0D0,
     :         31D0,       0D0,      0D0,      -13D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=311,320 ) /
     :         49D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :         32D0,       0D0,      0D0,      -13D0,     0D0,     0D0,
     :         23D0,       0D0,      0D0,      -12D0,     0D0,     0D0,
     :        -43D0,       0D0,      0D0,       18D0,     0D0,     0D0,
     :         26D0,       0D0,      0D0,      -11D0,     0D0,     0D0,
     :        -32D0,       0D0,      0D0,       14D0,     0D0,     0D0,
     :        -29D0,       0D0,      0D0,       14D0,     0D0,     0D0,
     :        -27D0,       0D0,      0D0,       12D0,     0D0,     0D0,
     :         30D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -11D0,       0D0,      0D0,        5D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=321,330 ) /
     :        -21D0,       0D0,      0D0,       10D0,     0D0,     0D0,
     :        -34D0,       0D0,      0D0,       15D0,     0D0,     0D0,
     :        -10D0,       0D0,      0D0,        6D0,     0D0,     0D0,
     :        -36D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -9D0,       0D0,      0D0,        4D0,     0D0,     0D0,
     :        -12D0,       0D0,      0D0,        5D0,     0D0,     0D0,
     :        -21D0,       0D0,      0D0,        5D0,     0D0,     0D0,
     :        -29D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :        -15D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :        -20D0,       0D0,      0D0,        0D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=331,340 ) /
     :         28D0,       0D0,      0D0,        0D0,     0D0,    -2D0,
     :         17D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -22D0,       0D0,      0D0,       12D0,     0D0,     0D0,
     :        -14D0,       0D0,      0D0,        7D0,     0D0,     0D0,
     :         24D0,       0D0,      0D0,      -11D0,     0D0,     0D0,
     :         11D0,       0D0,      0D0,       -6D0,     0D0,     0D0,
     :         14D0,       0D0,      0D0,       -6D0,     0D0,     0D0,
     :         24D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         18D0,       0D0,      0D0,       -8D0,     0D0,     0D0,
     :        -38D0,       0D0,      0D0,        0D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=341,350 ) /
     :        -31D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -16D0,       0D0,      0D0,        8D0,     0D0,     0D0,
     :         29D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -18D0,       0D0,      0D0,       10D0,     0D0,     0D0,
     :        -10D0,       0D0,      0D0,        5D0,     0D0,     0D0,
     :        -17D0,       0D0,      0D0,       10D0,     0D0,     0D0,
     :          9D0,       0D0,      0D0,       -4D0,     0D0,     0D0,
     :         16D0,       0D0,      0D0,       -6D0,     0D0,     0D0,
     :         22D0,       0D0,      0D0,      -12D0,     0D0,     0D0,
     :         20D0,       0D0,      0D0,        0D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=351,360 ) /
     :        -13D0,       0D0,      0D0,        6D0,     0D0,     0D0,
     :        -17D0,       0D0,      0D0,        9D0,     0D0,     0D0,
     :        -14D0,       0D0,      0D0,        8D0,     0D0,     0D0,
     :          0D0,       0D0,      0D0,       -7D0,     0D0,     0D0,
     :         14D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         19D0,       0D0,      0D0,      -10D0,     0D0,     0D0,
     :        -34D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -20D0,       0D0,      0D0,        8D0,     0D0,     0D0,
     :          9D0,       0D0,      0D0,       -5D0,     0D0,     0D0,
     :        -18D0,       0D0,      0D0,        7D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=361,370 ) /
     :         13D0,       0D0,      0D0,       -6D0,     0D0,     0D0,
     :         17D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -12D0,       0D0,      0D0,        5D0,     0D0,     0D0,
     :         15D0,       0D0,      0D0,       -8D0,     0D0,     0D0,
     :        -11D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :         13D0,       0D0,      0D0,       -5D0,     0D0,     0D0,
     :        -18D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -35D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          9D0,       0D0,      0D0,       -4D0,     0D0,     0D0,
     :        -19D0,       0D0,      0D0,       10D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=371,380 ) /
     :        -26D0,       0D0,      0D0,       11D0,     0D0,     0D0,
     :          8D0,       0D0,      0D0,       -4D0,     0D0,     0D0,
     :        -10D0,       0D0,      0D0,        4D0,     0D0,     0D0,
     :         10D0,       0D0,      0D0,       -6D0,     0D0,     0D0,
     :        -21D0,       0D0,      0D0,        9D0,     0D0,     0D0,
     :        -15D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          9D0,       0D0,      0D0,       -5D0,     0D0,     0D0,
     :        -29D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -19D0,       0D0,      0D0,       10D0,     0D0,     0D0,
     :         12D0,       0D0,      0D0,       -5D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=381,390 ) /
     :         22D0,       0D0,      0D0,       -9D0,     0D0,     0D0,
     :        -10D0,       0D0,      0D0,        5D0,     0D0,     0D0,
     :        -20D0,       0D0,      0D0,       11D0,     0D0,     0D0,
     :        -20D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -17D0,       0D0,      0D0,        7D0,     0D0,     0D0,
     :         15D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :          8D0,       0D0,      0D0,       -4D0,     0D0,     0D0,
     :         14D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -12D0,       0D0,      0D0,        6D0,     0D0,     0D0,
     :         25D0,       0D0,      0D0,        0D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=391,400 ) /
     :        -13D0,       0D0,      0D0,        6D0,     0D0,     0D0,
     :        -14D0,       0D0,      0D0,        8D0,     0D0,     0D0,
     :         13D0,       0D0,      0D0,       -5D0,     0D0,     0D0,
     :        -17D0,       0D0,      0D0,        9D0,     0D0,     0D0,
     :        -12D0,       0D0,      0D0,        6D0,     0D0,     0D0,
     :        -10D0,       0D0,      0D0,        5D0,     0D0,     0D0,
     :         10D0,       0D0,      0D0,       -6D0,     0D0,     0D0,
     :        -15D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -22D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         28D0,       0D0,      0D0,       -1D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=401,410 ) /
     :         15D0,       0D0,      0D0,       -7D0,     0D0,     0D0,
     :         23D0,       0D0,      0D0,      -10D0,     0D0,     0D0,
     :         12D0,       0D0,      0D0,       -5D0,     0D0,     0D0,
     :         29D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :        -25D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :         22D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -18D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         15D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :        -23D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         12D0,       0D0,      0D0,       -5D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=411,420 ) /
     :         -8D0,       0D0,      0D0,        4D0,     0D0,     0D0,
     :        -19D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -10D0,       0D0,      0D0,        4D0,     0D0,     0D0,
     :         21D0,       0D0,      0D0,       -9D0,     0D0,     0D0,
     :         23D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :        -16D0,       0D0,      0D0,        8D0,     0D0,     0D0,
     :        -19D0,       0D0,      0D0,        9D0,     0D0,     0D0,
     :        -22D0,       0D0,      0D0,       10D0,     0D0,     0D0,
     :         27D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         16D0,       0D0,      0D0,       -8D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=421,430 ) /
     :         19D0,       0D0,      0D0,       -8D0,     0D0,     0D0,
     :          9D0,       0D0,      0D0,       -4D0,     0D0,     0D0,
     :         -9D0,       0D0,      0D0,        4D0,     0D0,     0D0,
     :         -9D0,       0D0,      0D0,        4D0,     0D0,     0D0,
     :         -8D0,       0D0,      0D0,        4D0,     0D0,     0D0,
     :         18D0,       0D0,      0D0,       -9D0,     0D0,     0D0,
     :         16D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :        -10D0,       0D0,      0D0,        4D0,     0D0,     0D0,
     :        -23D0,       0D0,      0D0,        9D0,     0D0,     0D0,
     :         16D0,       0D0,      0D0,       -1D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=431,440 ) /
     :        -12D0,       0D0,      0D0,        6D0,     0D0,     0D0,
     :         -8D0,       0D0,      0D0,        4D0,     0D0,     0D0,
     :         30D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :         24D0,       0D0,      0D0,      -10D0,     0D0,     0D0,
     :         10D0,       0D0,      0D0,       -4D0,     0D0,     0D0,
     :        -16D0,       0D0,      0D0,        7D0,     0D0,     0D0,
     :        -16D0,       0D0,      0D0,        7D0,     0D0,     0D0,
     :         17D0,       0D0,      0D0,       -7D0,     0D0,     0D0,
     :        -24D0,       0D0,      0D0,       10D0,     0D0,     0D0,
     :        -12D0,       0D0,      0D0,        5D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=441,450 ) /
     :        -24D0,       0D0,      0D0,       11D0,     0D0,     0D0,
     :        -23D0,       0D0,      0D0,        9D0,     0D0,     0D0,
     :        -13D0,       0D0,      0D0,        5D0,     0D0,     0D0,
     :        -15D0,       0D0,      0D0,        7D0,     0D0,     0D0,
     :          0D0,       0D0,  -1988D0,        0D0,     0D0, -1679D0,
     :          0D0,       0D0,    -63D0,        0D0,     0D0,   -27D0,
     :         -4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          0D0,       0D0,      5D0,        0D0,     0D0,     4D0,
     :          5D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :          0D0,       0D0,    364D0,        0D0,     0D0,   176D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=451,460 ) /
     :          0D0,       0D0,  -1044D0,        0D0,     0D0,  -891D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          0D0,       0D0,    330D0,        0D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :         -5D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,        0D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=461,470 ) /
     :          3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          0D0,       0D0,      5D0,        0D0,     0D0,     0D0,
     :          0D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          6D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :         -7D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :        -12D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=471,480 ) /
     :         -5D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -7D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :          7D0,       0D0,      0D0,       -4D0,     0D0,     0D0,
     :          0D0,       0D0,    -12D0,        0D0,     0D0,   -10D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -7D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=481,490 ) /
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :          0D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :          7D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :         -5D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -5D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -2D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=491,500 ) /
     :         -8D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :          9D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          6D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :         -5D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -7D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        2D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=501,510 ) /
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         -5D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          9D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          9D0,       0D0,      0D0,       -3D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=511,520 ) /
     :         -4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          8D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :          6D0,       0D0,      0D0,       -3D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=521,530 ) /
     :          3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :         -7D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          9D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -5D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :        -13D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -7D0,       0D0,      0D0,        0D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=531,540 ) /
     :         10D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         10D0,       0D0,     13D0,        6D0,     0D0,    -5D0,
     :          0D0,       0D0,     30D0,        0D0,     0D0,    14D0,
     :          0D0,       0D0,   -162D0,        0D0,     0D0,  -138D0,
     :          0D0,       0D0,     75D0,        0D0,     0D0,     0D0,
     :         -7D0,       0D0,      0D0,        4D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -2D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=541,550 ) /
     :          5D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -5D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          6D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          9D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -7D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=551,560 ) /
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          7D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -6D0,       0D0,     -3D0,        3D0,     0D0,     1D0,
     :          0D0,       0D0,     -3D0,        0D0,     0D0,    -2D0,
     :         11D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         11D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        2D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=561,570 ) /
     :         -1D0,       0D0,      3D0,        3D0,     0D0,    -1D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          0D0,       0D0,    -13D0,        0D0,     0D0,   -11D0,
     :          3D0,       0D0,      6D0,        0D0,     0D0,     0D0,
     :         -7D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :         -7D0,       0D0,      0D0,        3D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=571,580 ) /
     :          8D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         11D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          8D0,       0D0,      0D0,       -4D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         11D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -6D0,       0D0,      0D0,        3D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=581,590 ) /
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -8D0,       0D0,      0D0,        4D0,     0D0,     0D0,
     :         -7D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :          6D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :         -6D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :          6D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          6D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -2D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=591,600 ) /
     :         -5D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          6D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          0D0,       0D0,    -26D0,        0D0,     0D0,   -11D0,
     :          0D0,       0D0,    -10D0,        0D0,     0D0,    -5D0,
     :          5D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :        -13D0,       0D0,      0D0,        0D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=601,610 ) /
     :          3D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          7D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -6D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -5D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -7D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -2D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=611,620 ) /
     :         13D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :        -11D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          6D0,       0D0,      0D0,       -3D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=621,630 ) /
     :          3D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :        -12D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :          0D0,       0D0,     -5D0,        0D0,     0D0,    -2D0,
     :         -7D0,       0D0,      0D0,        4D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=631,640 ) /
     :          6D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :         -5D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         12D0,       0D0,      0D0,        0D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=641,650 ) /
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          6D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :         -6D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          6D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :          6D0,       0D0,      0D0,        0D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=651,660 ) /
     :         -6D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          7D0,       0D0,      0D0,       -4D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :         -5D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         -6D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :         -6D0,       0D0,      0D0,        3D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         10D0,       0D0,      0D0,        0D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=661,670 ) /
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          7D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          7D0,       0D0,      0D0,       -3D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :         11D0,       0D0,      0D0,        0D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :         -6D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          5D0,       0D0,      0D0,       -2D0,     0D0,     0D0 /
      DATA ( ( CLS(I,J), I=1,6 ), J=671,678 ) /
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -4D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        2D0,     0D0,     0D0,
     :          4D0,       0D0,      0D0,       -2D0,     0D0,     0D0,
     :          3D0,       0D0,      0D0,       -1D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        1D0,     0D0,     0D0,
     :         -3D0,       0D0,      0D0,        2D0,     0D0,     0D0 /

*
*  Planetary argument multipliers:
*              L   L'  F   D   Om  Me  Ve  E  Ma  Ju  Sa  Ur  Ne  pre

      DATA ( ( NAPL(I,J), I=1,14 ), J=  1, 10 ) /
     :         0,  0,  0,  0,  0,  0,  0,  8,-16,  4,  5,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -8, 16, -4, -5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  8,-16,  4,  5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  2,  2,
     :         0,  0,  0,  0,  0,  0,  0, -4,  8, -1, -5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  4, -8,  3,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0,  3, -8,  3,  0,  0,  0,  0,
     :        -1,  0,  0,  0,  0,  0, 10, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0, -2,  6, -3,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J= 11, 20 ) /
     :         0,  0,  1, -1,  1,  0,  0, -5,  8, -3,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -4,  8, -3,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  4, -8,  1,  5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -5,  6,  4,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -5,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  2, -5,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -5,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0, -2,  5,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0, -2,  5,  0,  0,  1 /
      DATA ( ( NAPL(I,J), I=1,14 ), J= 21, 30 ) /
     :         0,  0,  0,  0,  0,  0,  0,  0,  0, -2,  5,  0,  0,  2,
     :         2,  0, -1, -1,  0,  0,  0,  3, -7,  0,  0,  0,  0,  0,
     :         1,  0,  0, -2,  0,  0, 19,-21,  3,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  2, -4,  0, -3,  0,  0,  0,  0,
     :         1,  0,  0, -1,  1,  0,  0, -1,  0,  2,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0, -4, 10,  0,  0,  0,
     :        -2,  0,  0,  2,  1,  0,  0,  2,  0,  0, -5,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  3, -7,  4,  0,  0,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  0,  1,  0,  1, -1,  0,  0,  0,
     :        -2,  0,  0,  2,  1,  0,  0,  2,  0, -2,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J= 31, 40 ) /
     :        -1,  0,  0,  0,  0,  0, 18,-16,  0,  0,  0,  0,  0,  0,
     :        -2,  0,  1,  1,  2,  0,  0,  1,  0, -2,  0,  0,  0,  0,
     :        -1,  0,  1, -1,  1,  0, 18,-17,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  0,  1,  1,  0,  0,  2, -2,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -8, 13,  0,  0,  0,  0,  0,  2,
     :         0,  0,  2, -2,  2,  0, -8, 11,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -8, 13,  0,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0, -8, 12,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  8,-13,  0,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  8,-14,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J= 41, 50 ) /
     :         0,  0,  0,  0,  0,  0,  8,-13,  0,  0,  0,  0,  0,  1,
     :        -2,  0,  0,  2,  1,  0,  0,  2,  0, -4,  5,  0,  0,  0,
     :        -2,  0,  0,  2,  2,  0,  3, -3,  0,  0,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  2,  0, -3,  1,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  3, -5,  0,  2,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  2,  0, -4,  3,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0, -1,  2,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  2,  0,  0, -2,  2,  0,  0,  0,  0,  0,
     :        -1,  0,  1,  0,  1,  0,  3, -5,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J= 51, 60 ) /
     :        -1,  0,  0,  1,  0,  0,  3, -4,  0,  0,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  2,  0, -2, -2,  0,  0,  0,
     :        -2,  0,  2,  0,  2,  0,  0, -5,  9,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  0, -1,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  0,  0,  2,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  1,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  2,
     :        -1,  0,  0,  1,  0,  0,  0,  3, -4,  0,  0,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  0,  1,  0,  0,  2,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J= 61, 70 ) /
     :         0,  0,  1, -1,  2,  0,  0, -1,  0,  0,  2,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0, -9, 17,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  2,  0, -3,  5,  0,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0, -1,  2,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,
     :         1,  0,  0, -2,  0,  0, 17,-16,  0, -2,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  1, -3,  0,  0,  0,
     :        -2,  0,  0,  2,  1,  0,  0,  5, -6,  0,  0,  0,  0,  0,
     :         0,  0, -2,  2,  0,  0,  0,  9,-13,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  2,  0,  0, -1,  0,  0,  1,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J= 71, 80 ) /
     :         0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  1,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  0,  1,  0,  0,  1,  0,  0,  0,
     :         0,  0, -2,  2,  0,  0,  5, -6,  0,  0,  0,  0,  0,  0,
     :         0,  0, -1,  1,  1,  0,  5, -7,  0,  0,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  6, -8,  0,  0,  0,  0,  0,  0,
     :         2,  0,  1, -3,  1,  0, -6,  7,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  2,  0,  0,  0,  0,  1,  0,  0,  0,  0,
     :         0,  0, -1,  1,  1,  0,  0,  1,  0,  1,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  0,  2,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  1 /
      DATA ( ( NAPL(I,J), I=1,14 ), J= 81, 90 ) /
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -8, 15,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -8, 15,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0, -9, 15,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  8,-15,  0,  0,  0,  0,  0,
     :         1,  0, -1, -1,  0,  0,  0,  8,-15,  0,  0,  0,  0,  0,
     :         2,  0,  0, -2,  0,  0,  2, -5,  0,  0,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  2,  0, -5,  5,  0,  0,  0,
     :         2,  0,  0, -2,  1,  0,  0, -6,  8,  0,  0,  0,  0,  0,
     :         2,  0,  0, -2,  1,  0,  0, -2,  0,  3,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J= 91,100 ) /
     :        -2,  0,  1,  1,  0,  0,  0,  1,  0, -3,  0,  0,  0,  0,
     :        -2,  0,  1,  1,  1,  0,  0,  1,  0, -3,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  2,  0, -3,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  6, -8,  0,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  2,  0, -1, -5,  0,  0,  0,
     :        -1,  0,  0,  1,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0,
     :        -1,  0,  1,  1,  1,  0,-20, 20,  0,  0,  0,  0,  0,  0,
     :         1,  0,  0, -2,  0,  0, 20,-21,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  8,-15,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  0,-10, 15,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=101,110 ) /
     :         0,  0, -1,  1,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  0,  0,  1,  0,  0,  0,  0,
     :         0,  0,  1, -1,  2,  0,  0, -1,  0,  1,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0, -2,  4,  0,  0,  0,
     :         2,  0,  0, -2,  1,  0, -6,  8,  0,  0,  0,  0,  0,  0,
     :         0,  0, -2,  2,  1,  0,  5, -6,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  0, -1,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  1,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=111,120 ) /
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  2,
     :         0,  0,  2, -2,  1,  0,  0, -9, 13,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  7,-13,  0,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  5, -6,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  9,-17,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -9, 17,  0,  0,  0,  0,  2,
     :         1,  0,  0, -1,  1,  0,  0, -3,  4,  0,  0,  0,  0,  0,
     :         1,  0,  0, -1,  1,  0, -3,  4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  2,  0,  0, -1,  2,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=121,130 ) /
     :         0,  0, -1,  1,  1,  0,  0,  0,  2,  0,  0,  0,  0,  0,
     :         0,  0, -2,  2,  0,  1,  0, -2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  3, -5,  0,  2,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  1,  0,  0,  2,  0, -3,  1,  0,  0,  0,
     :        -2,  0,  0,  2,  1,  0,  3, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  8,-13,  0,  0,  0,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  8,-12,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0, -8, 11,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  0,  1,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0,
     :        -1,  0,  0,  0,  1,  0, 18,-16,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=131,140 ) /
     :         0,  0,  1, -1,  1,  0,  0, -1,  0, -1,  1,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  3, -7,  4,  0,  0,  0,  0,  0,
     :        -2,  0,  1,  1,  1,  0,  0, -3,  7,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  2,  0,  0, -1,  0, -2,  5,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  0,  0, -2,  5,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0, -4,  8, -3,  0,  0,  0,  0,
     :         1,  0,  0,  0,  1,  0,-10,  3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  0, -2,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  0,  0,  1,  0, 10, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  4, -8,  3,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=141,150 ) /
     :         0,  0,  0,  0,  1,  0,  0,  0,  0,  2, -5,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  0,  1,  0,  2, -5,  0,  0,  0,
     :         2,  0, -1, -1,  1,  0,  0,  3, -7,  0,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  2,  0,  0, -5,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0, -3,  7, -4,  0,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,
     :         1,  0,  0,  0,  1,  0,-18, 16,  0,  0,  0,  0,  0,  0,
     :        -2,  0,  1,  1,  1,  0,  0,  1,  0, -2,  0,  0,  0,  0,
     :         0,  0,  1, -1,  2,  0, -8, 12,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0, -8, 13,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=151,160 ) /
     :         0,  0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0,  0, -2,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -2,  2,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -1,  2,  0,  0,  0,  0,  1,
     :        -1,  0,  0,  1,  1,  0,  3, -4,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  0,  1,  1,  0,  0,  3, -4,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  0, -2,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  2,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  1 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=161,170 ) /
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  2,
     :         0,  0,  1, -1,  0,  0,  3, -6,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0, -3,  5,  0,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  2,  0, -3,  4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0, -2,  4,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0, -5,  6,  0,  0,  0,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  5, -7,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  5, -8,  0,  0,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  1,  0,  6, -8,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0, -8, 15,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=171,180 ) /
     :        -2,  0,  0,  2,  1,  0,  0,  2,  0, -3,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  1,  0,  0,  6, -8,  0,  0,  0,  0,  0,
     :         1,  0,  0, -1,  1,  0,  0, -1,  0,  1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0, -1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  1 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=181,190 ) /
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  2,
     :         0,  0,  1, -1,  2,  0,  0, -1,  0,  0, -1,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  0,  0,  0, -1,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  0,  1,  0,  0, -1,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -7, 13,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  7,-13,  0,  0,  0,  0,  0,
     :         2,  0,  0, -2,  1,  0,  0, -5,  6,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  0, -8, 11,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1, -1,  0,  2,  0,  0,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  4, -4,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=191,200 ) /
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  3,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  2,
     :        -2,  0,  0,  2,  0,  0,  3, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  2,  0,  0, -4,  8, -3,  0,  0,  0,  0,
     :         0,  0,  0,  0,  2,  0,  0,  4, -8,  3,  0,  0,  0,  0,
     :         2,  0,  0, -2,  1,  0,  0, -2,  0,  2,  0,  0,  0,  0,
     :         0,  0,  1, -1,  2,  0,  0, -1,  0,  2,  0,  0,  0,  0,
     :         0,  0,  1, -1,  2,  0,  0,  0, -2,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=201,210 ) /
     :         0,  0,  0,  0,  1,  0,  0,  1, -2,  0,  0,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  0,  1,  0,  0, -2,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  0, -2,  0,  0,  2,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  3, -6,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0, -3,  4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -3,  5,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0, -3,  5,  0,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=211,220 ) /
     :         0,  0,  2, -2,  2,  0, -3,  3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -3,  5,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2, -4,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0,  1, -4,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  2, -4,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -2,  4,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0, -3,  4,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -2,  4,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0, -2,  4,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -5,  8,  0,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=221,230 ) /
     :         0,  0,  2, -2,  2,  0, -5,  6,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -5,  8,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -5,  8,  0,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0, -5,  7,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -5,  8,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  5, -8,  0,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  2,  0,  0, -1,  0, -1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  0,  0, -1,  0,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  0, -2,  0,  1,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=231,240 ) /
     :         0,  0,  0,  0,  0,  0,  0, -6, 11,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  6,-11,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0, -1,  0,  4,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  1,  0, -4,  0,  0,  0,  0,  0,  0,
     :         2,  0,  0, -2,  1,  0, -3,  3,  0,  0,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  2,  0,  0, -2,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  0, -7,  9,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  4, -5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  1 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=241,250 ) /
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  2,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  2,
     :         0,  0,  2, -2,  2,  0,  0, -2,  0,  2,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  5,  0,  0,  2,
     :         0,  0,  0,  0,  1,  0,  3, -5,  0,  0,  0,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  3, -4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0, -3,  3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  2, -4,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  0, -4,  4,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=251,260 ) /
     :         0,  0,  1, -1,  2,  0, -5,  7,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  3, -6,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -3,  6,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0, -4,  6,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -3,  6,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0, -3,  6,  0,  0,  0,  0,  2,
     :         0,  0, -1,  1,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  2, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -5,  9,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -5,  9,  0,  0,  0,  0,  1 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=261,270 ) /
     :         0,  0,  0,  0,  0,  0,  0,  5, -9,  0,  0,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  0,  1,  0, -2,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  0, -2,  0,  2,  0,  0,  0,  0,
     :        -2,  0,  1,  1,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :         0,  0, -2,  2,  0,  0,  3, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -6, 10,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0, -6, 10,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -2,  3,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -2,  3,  0,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0, -2,  2,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=271,280 ) /
     :         0,  0,  0,  0,  0,  0,  2, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  2, -3,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  3,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  4, -8,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -4,  8,  0,  0,  0,  0,  2,
     :         0,  0, -2,  2,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -4,  7,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=281,290 ) /
     :         0,  0,  0,  0,  0,  0,  0, -4,  7,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  4, -7,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0, -2,  3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  0, -2,  0,  3,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -5, 10,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  1,  0, -1,  2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  4,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -3,  5,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -3,  5,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=291,300 ) /
     :         0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  1, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -1,  2,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0, -1,  2,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -7, 11,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -7, 11,  0,  0,  0,  0,  0,  1,
     :         0,  0, -2,  2,  0,  0,  4, -4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  2, -3,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0, -4,  4,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=301,310 ) /
     :         0,  0, -1,  1,  0,  0,  4, -5,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -4,  7,  0,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0, -4,  6,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -4,  7,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -4,  6,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -4,  6,  0,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0, -4,  5,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -4,  6,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  4, -6,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=311,320 ) /
     :        -2,  0,  0,  2,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  1, -1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -1,  0,  5,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  1, -3,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -1,  3,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -7, 12,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -1,  1,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -1,  1,  0,  0,  0,  0,  0,  1 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=321,330 ) /
     :         0,  0,  1, -1,  1,  0, -1,  0,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  1, -2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -2,  5,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -1,  0,  4,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0, -4,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0, -1,  1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -6, 10,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -6, 10,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=331,340 ) /
     :         0,  0,  2, -2,  1,  0,  0, -3,  0,  3,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -3,  7,  0,  0,  0,  0,  2,
     :        -2,  0,  0,  2,  0,  0,  4, -4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -5,  8,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  5, -8,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -1,  0,  3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -1,  0,  3,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0, -3,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  2, -4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -2,  4,  0,  0,  0,  0,  0,  1 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=341,350 ) /
     :         0,  0,  1, -1,  1,  0, -2,  3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -2,  4,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -6,  9,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -6,  9,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  6, -9,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  1,  0, -2,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0, -2,  2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -4,  6,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  4, -6,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  3, -4,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=351,360 ) /
     :         0,  0,  0,  0,  0,  0,  0, -1,  0,  2,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0, -2,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  1,  0, -1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -5,  9,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  3, -4,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -3,  4,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -3,  4,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  3, -4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  3, -4,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  1,  0,  0,  2, -2,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=361,370 ) /
     :         0,  0,  0,  0,  1,  0,  0, -1,  0,  2,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  0, -3,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  1, -5,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -1,  0,  1,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0, -1,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0, -3,  5,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0, -3,  4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  0, -2,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=371,380 ) /
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  0, -1,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0, -1,  0,  1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0, -2,  2,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -8, 14,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  2, -5,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  5, -8,  3,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  5, -8,  3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  3, -8,  3,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=381,390 ) /
     :         0,  0,  0,  0,  0,  0,  0, -3,  8, -3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0, -2,  5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -8, 12,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -8, 12,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  1, -2,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  1,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  2,  0,  0,  2,
     :         0,  0,  2, -2,  1,  0, -5,  5,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=391,400 ) /
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  3, -6,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -3,  6,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0, -3,  6,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -1,  4,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -5,  7,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -5,  7,  0,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0, -5,  6,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=401,410 ) /
     :         0,  0,  0,  0,  0,  0,  5, -7,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  0, -1,  0,  1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -1,  0,  1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0, -1,  0,  3,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  2,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -2,  6,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  1,  0,  2, -2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -6,  9,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  6, -9,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -2,  2,  0,  0,  0,  0,  0,  1 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=411,420 ) /
     :         0,  0,  1, -1,  1,  0, -2,  1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -5,  7,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  5, -7,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0, -2,  2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  4, -5,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  1, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -1,  3,  0,  0,  0,  0,  0,  1 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=421,430 ) /
     :         0,  0,  1, -1,  1,  0, -1,  2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -1,  3,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -7, 10,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -7, 10,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -4,  8,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -4,  5,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -4,  5,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  4, -5,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=431,440 ) /
     :         0,  0,  0,  0,  0,  0,  0, -2,  0,  5,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -9, 13,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -1,  5,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -2,  0,  4,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0, -4,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -2,  7,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0, -3,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=441,450 ) /
     :         0,  0,  0,  0,  0,  0, -2,  5,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0, -2,  5,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -6,  8,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -6,  8,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  6, -8,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  2,  0, -2,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -3,  9,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  5, -6,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  5, -6,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=451,460 ) /
     :         0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -5, 10,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -3,  3,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  0, -3,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=461,470 ) /
     :         0,  0,  0,  0,  0,  0,  0, -5, 13,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0, -1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0, -1,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  0, -2,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  0, -2,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  3, -2,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  3, -2,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  0, -1,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -6, 15,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -8, 15,  0,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=471,480 ) /
     :         0,  0,  0,  0,  0,  0, -3,  9, -4,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  2, -5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -2,  8, -1, -5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  6, -8,  3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=481,490 ) /
     :         0,  0,  0,  0,  0,  0,  0, -6, 16, -4, -5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -2,  8, -3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -2,  8, -3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  6, -8,  1,  5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  3, -5,  4,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -8, 11,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -8, 11,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0, -8, 11,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, 11,  0,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=491,500 ) /
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  1,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  3, -3,  0,  2,  0,  0,  0,  2,
     :         0,  0,  2, -2,  1,  0,  0,  4, -8,  3,  0,  0,  0,  0,
     :         0,  0,  1, -1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  0, -4,  8, -3,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1,  2,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  1,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -3,  7,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  0,  4,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -5,  6,  0,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=501,510 ) /
     :         0,  0,  0,  0,  0,  0, -5,  6,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  5, -6,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  5, -6,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  2,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -1,  6,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  7, -9,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  2, -1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  2, -1,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  6, -7,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  5, -5,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=511,520 ) /
     :         0,  0,  0,  0,  0,  0, -1,  4,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0, -1,  4,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -7,  9,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -7,  9,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  4, -3,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  3, -1,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -4,  4,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=521,530 ) /
     :         0,  0,  0,  0,  0,  0,  0,  2,  1,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -3,  0,  5,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -9, 12,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  3,  0, -4,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  1, -1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  7, -8,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  3,  0, -3,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=531,540 ) /
     :         0,  0,  0,  0,  0,  0,  0,  3,  0, -3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -2,  6,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -6,  7,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  6, -7,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  6, -6,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  3,  0, -2,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  3,  0, -2,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  5, -4,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  3, -2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  3, -2,  0,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=541,550 ) /
     :         0,  0,  0,  0,  0,  0,  0,  3,  0, -1,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  3,  0, -1,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  3,  0,  0, -2,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  4, -2,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  3,  0,  0, -1,  0,  0,  2,
     :         0,  0,  2, -2,  1,  0,  0,  1,  0, -1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -8, 16,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  3,  0,  2, -5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  7, -8,  3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -5, 16, -4, -5,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=551,560 ) /
     :         0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -1,  8, -3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -8, 10,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -8, 10,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0, -8, 10,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  2,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  3,  0,  1,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -3,  8,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -5,  5,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  5, -5,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=561,570 ) /
     :         0,  0,  0,  0,  0,  0,  5, -5,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  5, -5,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  7, -7,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  7, -7,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  6, -5,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  7, -8,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  5, -3,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=571,580 ) /
     :         0,  0,  0,  0,  0,  0,  4, -3,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  1,  2,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -9, 11,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -9, 11,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  4,  0, -4,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  4,  0, -3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -6,  6,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  6, -6,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  6, -6,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  4,  0, -2,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=581,590 ) /
     :         0,  0,  0,  0,  0,  0,  0,  6, -4,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  3, -1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  3, -1,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  3, -1,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  4,  0, -1,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  4,  0,  0, -2,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  5, -2,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  8, -9,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  5, -4,  0,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=591,600 ) /
     :         0,  0,  0,  0,  0,  0,  2,  1,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  2,  1,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  2,  1,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0, -7,  7,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  7, -7,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  4, -2,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  4, -2,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  4, -2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  4, -2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  5,  0, -4,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=601,610 ) /
     :         0,  0,  0,  0,  0,  0,  0,  5,  0, -3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  5,  0, -2,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -8,  8,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  8, -8,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  5, -3,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  5, -3,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -9,  9,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0, -9,  9,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0, -9,  9,  0,  0,  0,  0,  0,  1 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=611,620 ) /
     :         0,  0,  0,  0,  0,  0,  9, -9,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  6, -4,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  6,  0,  0,  0,  0,  0,  2 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=621,630 ) /
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,
     :         1,  0,  0, -2,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,
     :         1,  0,  0, -2,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0,
     :         1,  0,  0, -2,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0,
     :         1,  0,  0, -2,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,
     :        -1,  0,  0,  2,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0,
     :         1,  0,  0, -2,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=631,640 ) /
     :        -1,  0,  0,  0,  0,  0,  0,  2,  0, -3,  0,  0,  0,  0,
     :        -1,  0,  0,  0,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0,
     :        -1,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  0,  2,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0,
     :         1,  0, -1,  1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  0,  2,  0,  0,  0,  2,  0, -3,  0,  0,  0,  0,
     :        -2,  0,  0,  0,  0,  0,  0,  2,  0, -3,  0,  0,  0,  0,
     :         1,  0,  0,  0,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0,
     :        -1,  0,  1, -1,  1,  0,  0, -1,  0,  0,  0,  0,  0,  0,
     :         1,  0,  1, -1,  1,  0,  0, -1,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=641,650 ) /
     :        -1,  0,  0,  0,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0,
     :        -1,  0,  0,  2,  1,  0,  0,  2,  0, -2,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,
     :        -1,  0,  0,  2,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,
     :        -1,  0,  0,  2,  0,  0,  3, -3,  0,  0,  0,  0,  0,  0,
     :         1,  0,  0, -2,  1,  0,  0, -2,  0,  2,  0,  0,  0,  0,
     :         1,  0,  2, -2,  2,  0, -3,  3,  0,  0,  0,  0,  0,  0,
     :         1,  0,  2, -2,  2,  0,  0, -2,  0,  2,  0,  0,  0,  0,
     :         1,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0,
     :         1,  0,  0,  0,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=651,660 ) /
     :         0,  0,  0, -2,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0, -2,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0,
     :         0,  0,  2,  0,  2,  0, -2,  2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2,  0,  2,  0,  0, -1,  0,  1,  0,  0,  0,  0,
     :         0,  0,  2,  0,  2,  0, -1,  1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2,  0,  2,  0, -2,  3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  2,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,
     :         0,  0,  1,  1,  2,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :         1,  0,  2,  0,  2,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  2,  0,  2,  0, 10, -3,  0,  0,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=661,670 ) /
     :         0,  0,  1,  1,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :         1,  0,  2,  0,  2,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2,  0,  2,  0,  0,  4, -8,  3,  0,  0,  0,  0,
     :         0,  0,  2,  0,  2,  0,  0, -4,  8, -3,  0,  0,  0,  0,
     :        -1,  0,  2,  0,  2,  0,  0, -4,  8, -3,  0,  0,  0,  0,
     :         2,  0,  2, -2,  2,  0,  0, -2,  0,  3,  0,  0,  0,  0,
     :         1,  0,  2,  0,  1,  0,  0, -2,  0,  3,  0,  0,  0,  0,
     :         0,  0,  1,  1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  2,  0,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :        -2,  0,  2,  2,  2,  0,  0,  2,  0, -2,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=671,680 ) /
     :         0,  0,  2,  0,  2,  0,  2, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2,  0,  2,  0,  1, -1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2,  0,  2,  0,  0,  1,  0, -1,  0,  0,  0,  0,
     :         0,  0,  2,  0,  2,  0,  2, -2,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  2,  2,  2,  0,  0, -1,  0,  1,  0,  0,  0,  0,
     :         1,  0,  2,  0,  2,  0, -1,  1,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  2,  2,  2,  0,  0,  2,  0, -3,  0,  0,  0,  0,
     :         2,  0,  2,  0,  2,  0,  0,  2,  0, -3,  0,  0,  0,  0,
     :         1,  0,  2,  0,  2,  0,  0, -4,  8, -3,  0,  0,  0,  0,
     :         1,  0,  2,  0,  2,  0,  0,  4, -8,  3,  0,  0,  0,  0 /
      DATA ( ( NAPL(I,J), I=1,14 ), J=681,687 ) /
     :         1,  0,  1,  1,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2,  0,  2,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :         2,  0,  2,  0,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  2,  2,  2,  0,  0,  2,  0, -2,  0,  0,  0,  0,
     :        -1,  0,  2,  2,  2,  0,  3, -3,  0,  0,  0,  0,  0,  0,
     :         1,  0,  2,  0,  2,  0,  1, -1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2,  2,  2,  0,  0,  2,  0, -2,  0,  0,  0,  0 /

*
*  Planetary nutation coefficients, unit 1e-7 arcsec:
*  longitude (sin, cos), obliquity (sin, cos)
*

      DATA ( ( ICPL(I,J), I=1,4 ), J=  1, 10 ) /
     :       1440,          0,          0,          0,
     :         56,       -117,        -42,        -40,
     :        125,        -43,          0,        -54,
     :          0,          5,          0,          0,
     :          3,         -7,         -3,          0,
     :          3,          0,          0,         -2,
     :       -114,          0,          0,         61,
     :       -219,         89,          0,          0,
     :         -3,          0,          0,          0,
     :       -462,       1604,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J= 11, 20 ) /
     :         99,          0,          0,        -53,
     :         -3,          0,          0,          2,
     :          0,          6,          2,          0,
     :          3,          0,          0,          0,
     :        -12,          0,          0,          0,
     :         14,       -218,        117,          8,
     :         31,       -481,       -257,        -17,
     :       -491,        128,          0,          0,
     :      -3084,       5123,       2735,       1647,
     :      -1444,       2409,      -1286,       -771 /
      DATA ( ( ICPL(I,J), I=1,4 ), J= 21, 30 ) /
     :         11,        -24,        -11,         -9,
     :         26,         -9,          0,          0,
     :        103,        -60,          0,          0,
     :          0,        -13,         -7,          0,
     :        -26,        -29,        -16,         14,
     :          9,        -27,        -14,         -5,
     :         12,          0,          0,         -6,
     :         -7,          0,          0,          0,
     :          0,         24,          0,          0,
     :        284,          0,          0,       -151 /
      DATA ( ( ICPL(I,J), I=1,4 ), J= 31, 40 ) /
     :        226,        101,          0,          0,
     :          0,         -8,         -2,          0,
     :          0,         -6,         -3,          0,
     :          5,          0,          0,         -3,
     :        -41,        175,         76,         17,
     :          0,         15,          6,          0,
     :        425,        212,       -133,        269,
     :       1200,        598,        319,       -641,
     :        235,        334,          0,          0,
     :         11,        -12,         -7,         -6 /
      DATA ( ( ICPL(I,J), I=1,4 ), J= 41, 50 ) /
     :          5,         -6,          3,          3,
     :         -5,          0,          0,          3,
     :          6,          0,          0,         -3,
     :         15,          0,          0,          0,
     :         13,          0,          0,         -7,
     :         -6,         -9,          0,          0,
     :        266,        -78,          0,          0,
     :       -460,       -435,       -232,        246,
     :          0,         15,          7,          0,
     :         -3,          0,          0,          2 /
      DATA ( ( ICPL(I,J), I=1,4 ), J= 51, 60 ) /
     :          0,        131,          0,          0,
     :          4,          0,          0,          0,
     :          0,          3,          0,          0,
     :          0,          4,          2,          0,
     :          0,          3,          0,          0,
     :        -17,        -19,        -10,          9,
     :         -9,        -11,          6,         -5,
     :         -6,          0,          0,          3,
     :        -16,          8,          0,          0,
     :          0,          3,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J= 61, 70 ) /
     :         11,         24,         11,         -5,
     :         -3,         -4,         -2,          1,
     :          3,          0,          0,         -1,
     :          0,         -8,         -4,          0,
     :          0,          3,          0,          0,
     :          0,          5,          0,          0,
     :          0,          3,          2,          0,
     :         -6,          4,          2,          3,
     :         -3,         -5,          0,          0,
     :         -5,          0,          0,          2 /
      DATA ( ( ICPL(I,J), I=1,4 ), J= 71, 80 ) /
     :          4,         24,         13,         -2,
     :        -42,         20,          0,          0,
     :        -10,        233,          0,          0,
     :         -3,          0,          0,          1,
     :         78,        -18,          0,          0,
     :          0,          3,          1,          0,
     :          0,         -3,         -1,          0,
     :          0,         -4,         -2,          1,
     :          0,         -8,         -4,         -1,
     :          0,         -5,          3,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J= 81, 90 ) /
     :         -7,          0,          0,          3,
     :        -14,          8,          3,          6,
     :          0,          8,         -4,          0,
     :          0,         19,         10,          0,
     :         45,        -22,          0,          0,
     :         -3,          0,          0,          0,
     :          0,         -3,          0,          0,
     :          0,          3,          0,          0,
     :          3,          5,          3,         -2,
     :         89,        -16,         -9,        -48 /
      DATA ( ( ICPL(I,J), I=1,4 ), J= 91,100 ) /
     :          0,          3,          0,          0,
     :         -3,          7,          4,          2,
     :       -349,        -62,          0,          0,
     :        -15,         22,          0,          0,
     :         -3,          0,          0,          0,
     :        -53,          0,          0,          0,
     :          5,          0,          0,         -3,
     :          0,         -8,          0,          0,
     :         15,         -7,         -4,         -8,
     :         -3,          0,          0,          1 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=101,110 ) /
     :        -21,        -78,          0,          0,
     :         20,        -70,        -37,        -11,
     :          0,          6,          3,          0,
     :          5,          3,          2,         -2,
     :        -17,         -4,         -2,          9,
     :          0,          6,          3,          0,
     :         32,         15,         -8,         17,
     :        174,         84,         45,        -93,
     :         11,         56,          0,          0,
     :        -66,        -12,         -6,         35 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=111,120 ) /
     :         47,          8,          4,        -25,
     :          0,          8,          4,          0,
     :         10,        -22,        -12,         -5,
     :         -3,          0,          0,          2,
     :        -24,         12,          0,          0,
     :          5,         -6,          0,          0,
     :          3,          0,          0,         -2,
     :          4,          3,          1,         -2,
     :          0,         29,         15,          0,
     :         -5,         -4,         -2,          2 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=121,130 ) /
     :          8,         -3,         -1,         -5,
     :          0,         -3,          0,          0,
     :         10,          0,          0,          0,
     :          3,          0,          0,         -2,
     :         -5,          0,          0,          3,
     :         46,         66,         35,        -25,
     :        -14,          7,          0,          0,
     :          0,          3,          2,          0,
     :         -5,          0,          0,          0,
     :        -68,        -34,        -18,         36 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=131,140 ) /
     :          0,         14,          7,          0,
     :         10,         -6,         -3,         -5,
     :         -5,         -4,         -2,          3,
     :         -3,          5,          2,          1,
     :         76,         17,          9,        -41,
     :         84,        298,        159,        -45,
     :          3,          0,          0,         -1,
     :         -3,          0,          0,          2,
     :         -3,          0,          0,          1,
     :        -82,        292,        156,         44 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=141,150 ) /
     :        -73,         17,          9,         39,
     :         -9,        -16,          0,          0,
     :          3,          0,         -1,         -2,
     :         -3,          0,          0,          0,
     :         -9,         -5,         -3,          5,
     :       -439,          0,          0,          0,
     :         57,        -28,        -15,        -30,
     :          0,         -6,         -3,          0,
     :         -4,          0,          0,          2,
     :        -40,         57,         30,         21 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=151,160 ) /
     :         23,          7,          3,        -13,
     :        273,         80,         43,       -146,
     :       -449,        430,          0,          0,
     :         -8,        -47,        -25,          4,
     :          6,         47,         25,         -3,
     :          0,         23,         13,          0,
     :         -3,          0,          0,          2,
     :          3,         -4,         -2,         -2,
     :        -48,       -110,        -59,         26,
     :         51,        114,         61,        -27 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=161,170 ) /
     :       -133,          0,          0,         57,
     :          0,          4,          0,          0,
     :        -21,         -6,         -3,         11,
     :          0,         -3,         -1,          0,
     :        -11,        -21,        -11,          6,
     :        -18,       -436,       -233,          9,
     :         35,         -7,          0,          0,
     :          0,          5,          3,          0,
     :         11,         -3,         -1,         -6,
     :         -5,         -3,         -1,          3 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=171,180 ) /
     :        -53,         -9,         -5,         28,
     :          0,          3,          2,          1,
     :          4,          0,          0,         -2,
     :          0,         -4,          0,          0,
     :        -50,        194,        103,         27,
     :        -13,         52,         28,          7,
     :        -91,        248,          0,          0,
     :          6,         49,         26,         -3,
     :         -6,        -47,        -25,          3,
     :          0,          5,          3,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=181,190 ) /
     :         52,         23,         10,        -23,
     :         -3,          0,          0,          1,
     :          0,          5,          3,          0,
     :         -4,          0,          0,          0,
     :         -4,          8,          3,          2,
     :         10,          0,          0,          0,
     :          3,          0,          0,         -2,
     :          0,          8,          4,          0,
     :          0,          8,          4,          1,
     :         -4,          0,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=191,200 ) /
     :         -4,          0,          0,          0,
     :         -8,          4,          2,          4,
     :          8,         -4,         -2,         -4,
     :          0,         15,          7,          0,
     :       -138,          0,          0,          0,
     :          0,         -7,         -3,          0,
     :          0,         -7,         -3,          0,
     :         54,          0,          0,        -29,
     :          0,         10,          4,          0,
     :         -7,          0,          0,          3 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=201,210 ) /
     :        -37,         35,         19,         20,
     :          0,          4,          0,          0,
     :         -4,          9,          0,          0,
     :          8,          0,          0,         -4,
     :         -9,        -14,         -8,          5,
     :         -3,         -9,         -5,          3,
     :       -145,         47,          0,          0,
     :        -10,         40,         21,          5,
     :         11,        -49,        -26,         -7,
     :      -2150,          0,          0,        932 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=211,220 ) /
     :        -12,          0,          0,          5,
     :         85,          0,          0,        -37,
     :          4,          0,          0,         -2,
     :          3,          0,          0,         -2,
     :        -86,        153,          0,          0,
     :         -6,          9,          5,          3,
     :          9,        -13,         -7,         -5,
     :         -8,         12,          6,          4,
     :        -51,          0,          0,         22,
     :        -11,       -268,       -116,          5 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=221,230 ) /
     :          0,         12,          5,          0,
     :          0,          7,          3,          0,
     :         31,          6,          3,        -17,
     :        140,         27,         14,        -75,
     :         57,         11,          6,        -30,
     :        -14,        -39,          0,          0,
     :          0,         -6,         -2,          0,
     :          4,         15,          8,         -2,
     :          0,          4,          0,          0,
     :         -3,          0,          0,          1 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=231,240 ) /
     :          0,         11,          5,          0,
     :          9,          6,          0,          0,
     :         -4,         10,          4,          2,
     :          5,          3,          0,          0,
     :         16,          0,          0,         -9,
     :         -3,          0,          0,          0,
     :          0,          3,          2,         -1,
     :          7,          0,          0,         -3,
     :        -25,         22,          0,          0,
     :         42,        223,        119,        -22 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=241,250 ) /
     :        -27,       -143,        -77,         14,
     :          9,         49,         26,         -5,
     :      -1166,          0,          0,        505,
     :         -5,          0,          0,          2,
     :         -6,          0,          0,          3,
     :         -8,          0,          1,          4,
     :          0,         -4,          0,          0,
     :        117,          0,          0,        -63,
     :         -4,          8,          4,          2,
     :          3,          0,          0,         -2 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=251,260 ) /
     :         -5,          0,          0,          2,
     :          0,         31,          0,          0,
     :         -5,          0,          1,          3,
     :          4,          0,          0,         -2,
     :         -4,          0,          0,          2,
     :        -24,        -13,         -6,         10,
     :          3,          0,          0,          0,
     :          0,        -32,        -17,          0,
     :          8,         12,          5,         -3,
     :          3,          0,          0,         -1 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=261,270 ) /
     :          7,         13,          0,          0,
     :         -3,         16,          0,          0,
     :         50,          0,          0,        -27,
     :          0,         -5,         -3,          0,
     :         13,          0,          0,          0,
     :          0,          5,          3,          1,
     :         24,          5,          2,        -11,
     :          5,        -11,         -5,         -2,
     :         30,         -3,         -2,        -16,
     :         18,          0,          0,         -9 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=271,280 ) /
     :          8,        614,          0,          0,
     :          3,         -3,         -1,         -2,
     :          6,         17,          9,         -3,
     :         -3,         -9,         -5,          2,
     :          0,          6,          3,         -1,
     :       -127,         21,          9,         55,
     :          3,          5,          0,          0,
     :         -6,        -10,         -4,          3,
     :          5,          0,          0,          0,
     :         16,          9,          4,         -7 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=281,290 ) /
     :          3,          0,          0,         -2,
     :          0,         22,          0,          0,
     :          0,         19,         10,          0,
     :          7,          0,          0,         -4,
     :          0,         -5,         -2,          0,
     :          0,          3,          1,          0,
     :         -9,          3,          1,          4,
     :         17,          0,          0,         -7,
     :          0,         -3,         -2,         -1,
     :        -20,         34,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=291,300 ) /
     :        -10,          0,          1,          5,
     :         -4,          0,          0,          2,
     :         22,        -87,          0,          0,
     :         -4,          0,          0,          2,
     :         -3,         -6,         -2,          1,
     :        -16,         -3,         -1,          7,
     :          0,         -3,         -2,          0,
     :          4,          0,          0,          0,
     :        -68,         39,          0,          0,
     :         27,          0,          0,        -14 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=301,310 ) /
     :          0,         -4,          0,          0,
     :        -25,          0,          0,          0,
     :        -12,         -3,         -2,          6,
     :          3,          0,          0,         -1,
     :          3,         66,         29,         -1,
     :        490,          0,          0,       -213,
     :        -22,         93,         49,         12,
     :         -7,         28,         15,          4,
     :         -3,         13,          7,          2,
     :        -46,         14,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=311,320 ) /
     :         -5,          0,          0,          0,
     :          2,          1,          0,          0,
     :          0,         -3,          0,          0,
     :        -28,          0,          0,         15,
     :          5,          0,          0,         -2,
     :          0,          3,          0,          0,
     :        -11,          0,          0,          5,
     :          0,          3,          1,          0,
     :         -3,          0,          0,          1,
     :         25,        106,         57,        -13 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=321,330 ) /
     :          5,         21,         11,         -3,
     :       1485,          0,          0,          0,
     :         -7,        -32,        -17,          4,
     :          0,          5,          3,          0,
     :         -6,         -3,         -2,          3,
     :         30,         -6,         -2,        -13,
     :         -4,          4,          0,          0,
     :        -19,          0,          0,         10,
     :          0,          4,          2,         -1,
     :          0,          3,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=331,340 ) /
     :          4,          0,          0,         -2,
     :          0,         -3,         -1,          0,
     :         -3,          0,          0,          0,
     :          5,          3,          1,         -2,
     :          0,         11,          0,          0,
     :        118,          0,          0,        -52,
     :          0,         -5,         -3,          0,
     :        -28,         36,          0,          0,
     :          5,         -5,          0,          0,
     :         14,        -59,        -31,         -8 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=341,350 ) /
     :          0,          9,          5,          1,
     :       -458,          0,          0,        198,
     :          0,        -45,        -20,          0,
     :          9,          0,          0,         -5,
     :          0,         -3,          0,          0,
     :          0,         -4,         -2,         -1,
     :         11,          0,          0,         -6,
     :          6,          0,          0,         -2,
     :        -16,         23,          0,          0,
     :          0,         -4,         -2,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=351,360 ) /
     :         -5,          0,          0,          2,
     :       -166,        269,          0,          0,
     :         15,          0,          0,         -8,
     :         10,          0,          0,         -4,
     :        -78,         45,          0,          0,
     :          0,         -5,         -2,          0,
     :          7,          0,          0,         -4,
     :         -5,        328,          0,          0,
     :          3,          0,          0,         -2,
     :          5,          0,          0,         -2 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=361,370 ) /
     :          0,          3,          1,          0,
     :         -3,          0,          0,          0,
     :         -3,          0,          0,          0,
     :          0,         -4,         -2,          0,
     :      -1223,        -26,          0,          0,
     :          0,          7,          3,          0,
     :          3,          0,          0,          0,
     :          0,          3,          2,          0,
     :         -6,         20,          0,          0,
     :       -368,          0,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=371,380 ) /
     :        -75,          0,          0,          0,
     :         11,          0,          0,         -6,
     :          3,          0,          0,         -2,
     :         -3,          0,          0,          1,
     :        -13,        -30,          0,          0,
     :         21,          3,          0,          0,
     :         -3,          0,          0,          1,
     :         -4,          0,          0,          2,
     :          8,        -27,          0,          0,
     :        -19,        -11,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=381,390 ) /
     :         -4,          0,          0,          2,
     :          0,          5,          2,          0,
     :         -6,          0,          0,          2,
     :         -8,          0,          0,          0,
     :         -1,          0,          0,          0,
     :        -14,          0,          0,          6,
     :          6,          0,          0,          0,
     :        -74,          0,          0,         32,
     :          0,         -3,         -1,          0,
     :          4,          0,          0,         -2 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=391,400 ) /
     :          8,         11,          0,          0,
     :          0,          3,          2,          0,
     :       -262,          0,          0,        114,
     :          0,         -4,          0,          0,
     :         -7,          0,          0,          4,
     :          0,        -27,        -12,          0,
     :        -19,         -8,         -4,          8,
     :        202,          0,          0,        -87,
     :         -8,         35,         19,          5,
     :          0,          4,          2,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=401,410 ) /
     :         16,         -5,          0,          0,
     :          5,          0,          0,         -3,
     :          0,         -3,          0,          0,
     :          1,          0,          0,          0,
     :        -35,        -48,        -21,         15,
     :         -3,         -5,         -2,          1,
     :          6,          0,          0,         -3,
     :          3,          0,          0,         -1,
     :          0,         -5,          0,          0,
     :         12,         55,         29,         -6 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=411,420 ) /
     :          0,          5,          3,          0,
     :       -598,          0,          0,          0,
     :         -3,        -13,         -7,          1,
     :         -5,         -7,         -3,          2,
     :          3,          0,          0,         -1,
     :          5,         -7,          0,          0,
     :          4,          0,          0,         -2,
     :         16,         -6,          0,          0,
     :          8,         -3,          0,          0,
     :          8,        -31,        -16,         -4 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=421,430 ) /
     :          0,          3,          1,          0,
     :        113,          0,          0,        -49,
     :          0,        -24,        -10,          0,
     :          4,          0,          0,         -2,
     :         27,          0,          0,          0,
     :         -3,          0,          0,          1,
     :          0,         -4,         -2,          0,
     :          5,          0,          0,         -2,
     :          0,         -3,          0,          0,
     :        -13,          0,          0,          6 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=431,440 ) /
     :          5,          0,          0,         -2,
     :        -18,        -10,         -4,          8,
     :         -4,        -28,          0,          0,
     :         -5,          6,          3,          2,
     :         -3,          0,          0,          1,
     :         -5,         -9,         -4,          2,
     :         17,          0,          0,         -7,
     :         11,          4,          0,          0,
     :          0,         -6,         -2,          0,
     :         83,         15,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=441,450 ) /
     :         -4,          0,          0,          2,
     :          0,       -114,        -49,          0,
     :        117,          0,          0,        -51,
     :         -5,         19,         10,          2,
     :         -3,          0,          0,          0,
     :         -3,          0,          0,          2,
     :          0,         -3,         -1,          0,
     :          3,          0,          0,          0,
     :          0,         -6,         -2,          0,
     :        393,          3,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=451,460 ) /
     :         -4,         21,         11,          2,
     :         -6,          0,         -1,          3,
     :         -3,          8,          4,          1,
     :          8,          0,          0,          0,
     :         18,        -29,        -13,         -8,
     :          8,         34,         18,         -4,
     :         89,          0,          0,          0,
     :          3,         12,          6,         -1,
     :         54,        -15,         -7,        -24,
     :          0,          3,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=461,470 ) /
     :          3,          0,          0,         -1,
     :          0,         35,          0,          0,
     :       -154,        -30,        -13,         67,
     :         15,          0,          0,          0,
     :          0,          4,          2,          0,
     :          0,          9,          0,          0,
     :         80,        -71,        -31,        -35,
     :          0,        -20,         -9,          0,
     :         11,          5,          2,         -5,
     :         61,        -96,        -42,        -27 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=471,480 ) /
     :         14,          9,          4,         -6,
     :        -11,         -6,         -3,          5,
     :          0,         -3,         -1,          0,
     :        123,       -415,       -180,        -53,
     :          0,          0,          0,        -35,
     :         -5,          0,          0,          0,
     :          7,        -32,        -17,         -4,
     :          0,         -9,         -5,          0,
     :          0,         -4,          2,          0,
     :        -89,          0,          0,         38 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=481,490 ) /
     :          0,        -86,        -19,         -6,
     :          0,          0,        -19,          6,
     :       -123,       -416,       -180,         53,
     :          0,         -3,         -1,          0,
     :         12,         -6,         -3,         -5,
     :        -13,          9,          4,          6,
     :          0,        -15,         -7,          0,
     :          3,          0,          0,         -1,
     :        -62,        -97,        -42,         27,
     :        -11,          5,          2,          5 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=491,500 ) /
     :          0,        -19,         -8,          0,
     :         -3,          0,          0,          1,
     :          0,          4,          2,          0,
     :          0,          3,          0,          0,
     :          0,          4,          2,          0,
     :        -85,        -70,        -31,         37,
     :        163,        -12,         -5,        -72,
     :        -63,        -16,         -7,         28,
     :        -21,        -32,        -14,          9,
     :          0,         -3,         -1,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=501,510 ) /
     :          3,          0,          0,         -2,
     :          0,          8,          0,          0,
     :          3,         10,          4,         -1,
     :          3,          0,          0,         -1,
     :          0,         -7,         -3,          0,
     :          0,         -4,         -2,          0,
     :          6,         19,          0,          0,
     :          5,       -173,        -75,         -2,
     :          0,         -7,         -3,          0,
     :          7,        -12,         -5,         -3 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=511,520 ) /
     :         -3,          0,          0,          2,
     :          3,         -4,         -2,         -1,
     :         74,          0,          0,        -32,
     :         -3,         12,          6,          2,
     :         26,        -14,         -6,        -11,
     :         19,          0,          0,         -8,
     :          6,         24,         13,         -3,
     :         83,          0,          0,          0,
     :          0,        -10,         -5,          0,
     :         11,         -3,         -1,         -5 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=521,530 ) /
     :          3,          0,          1,         -1,
     :          3,          0,          0,         -1,
     :         -4,          0,          0,          0,
     :          5,        -23,        -12,         -3,
     :       -339,          0,          0,        147,
     :          0,        -10,         -5,          0,
     :          5,          0,          0,          0,
     :          3,          0,          0,         -1,
     :          0,         -4,         -2,          0,
     :         18,         -3,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=531,540 ) /
     :          9,        -11,         -5,         -4,
     :         -8,          0,          0,          4,
     :          3,          0,          0,         -1,
     :          0,          9,          0,          0,
     :          6,         -9,         -4,         -2,
     :         -4,        -12,          0,          0,
     :         67,        -91,        -39,        -29,
     :         30,        -18,         -8,        -13,
     :          0,          0,          0,          0,
     :          0,       -114,        -50,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=541,550 ) /
     :          0,          0,          0,         23,
     :        517,         16,          7,       -224,
     :          0,         -7,         -3,          0,
     :        143,         -3,         -1,        -62,
     :         29,          0,          0,        -13,
     :         -4,          0,          0,          2,
     :         -6,          0,          0,          3,
     :          5,         12,          5,         -2,
     :        -25,          0,          0,         11,
     :         -3,          0,          0,          1 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=551,560 ) /
     :          0,          4,          2,          0,
     :        -22,         12,          5,         10,
     :         50,          0,          0,        -22,
     :          0,          7,          4,          0,
     :          0,          3,          1,          0,
     :         -4,          4,          2,          2,
     :         -5,        -11,         -5,          2,
     :          0,          4,          2,          0,
     :          4,         17,          9,         -2,
     :         59,          0,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=561,570 ) /
     :          0,         -4,         -2,          0,
     :         -8,          0,          0,          4,
     :         -3,          0,          0,          0,
     :          4,        -15,         -8,         -2,
     :        370,         -8,          0,       -160,
     :          0,          0,         -3,          0,
     :          0,          3,          1,          0,
     :         -6,          3,          1,          3,
     :          0,          6,          0,          0,
     :        -10,          0,          0,          4 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=571,580 ) /
     :          0,          9,          4,          0,
     :          4,         17,          7,         -2,
     :         34,          0,          0,        -15,
     :          0,          5,          3,          0,
     :         -5,          0,          0,          2,
     :        -37,         -7,         -3,         16,
     :          3,         13,          7,         -2,
     :         40,          0,          0,          0,
     :          0,         -3,         -2,          0,
     :       -184,         -3,         -1,         80 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=581,590 ) /
     :         -3,          0,          0,          1,
     :         -3,          0,          0,          0,
     :          0,        -10,         -6,         -1,
     :         31,         -6,          0,        -13,
     :         -3,        -32,        -14,          1,
     :         -7,          0,          0,          3,
     :          0,         -8,         -4,          0,
     :          3,         -4,          0,          0,
     :          0,          4,          0,          0,
     :          0,          3,          1,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=591,600 ) /
     :         19,        -23,        -10,          2,
     :          0,          0,          0,        -10,
     :          0,          3,          2,          0,
     :          0,          9,          5,         -1,
     :         28,          0,          0,          0,
     :          0,         -7,         -4,          0,
     :          8,         -4,          0,         -4,
     :          0,          0,         -2,          0,
     :          0,          3,          0,          0,
     :         -3,          0,          0,          1 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=601,610 ) /
     :         -9,          0,          1,          4,
     :          3,         12,          5,         -1,
     :         17,         -3,         -1,          0,
     :          0,          7,          4,          0,
     :         19,          0,          0,          0,
     :          0,         -5,         -3,          0,
     :         14,         -3,          0,         -1,
     :          0,          0,         -1,          0,
     :          0,          0,          0,         -5,
     :          0,          5,          3,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=611,620 ) /
     :         13,          0,          0,          0,
     :          0,         -3,         -2,          0,
     :          2,          9,          4,          3,
     :          0,          0,          0,         -4,
     :          8,          0,          0,          0,
     :          0,          4,          2,          0,
     :          6,          0,          0,         -3,
     :          6,          0,          0,          0,
     :          0,          3,          1,          0,
     :          5,          0,          0,         -2 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=621,630 ) /
     :          3,          0,          0,         -1,
     :         -3,          0,          0,          0,
     :          6,          0,          0,          0,
     :          7,          0,          0,          0,
     :         -4,          0,          0,          0,
     :          4,          0,          0,          0,
     :          6,          0,          0,          0,
     :          0,         -4,          0,          0,
     :          0,         -4,          0,          0,
     :          5,          0,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=631,640 ) /
     :         -3,          0,          0,          0,
     :          4,          0,          0,          0,
     :         -5,          0,          0,          0,
     :          4,          0,          0,          0,
     :          0,          3,          0,          0,
     :         13,          0,          0,          0,
     :         21,         11,          0,          0,
     :          0,         -5,          0,          0,
     :          0,         -5,         -2,          0,
     :          0,          5,          3,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=641,650 ) /
     :          0,         -5,          0,          0,
     :         -3,          0,          0,          2,
     :         20,         10,          0,          0,
     :        -34,          0,          0,          0,
     :        -19,          0,          0,          0,
     :          3,          0,          0,         -2,
     :         -3,          0,          0,          1,
     :         -6,          0,          0,          3,
     :         -4,          0,          0,          0,
     :          3,          0,          0,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=651,660 ) /
     :          3,          0,          0,          0,
     :          4,          0,          0,          0,
     :          3,          0,          0,         -1,
     :          6,          0,          0,         -3,
     :         -8,          0,          0,          3,
     :          0,          3,          1,          0,
     :         -3,          0,          0,          0,
     :          0,         -3,         -2,          0,
     :        126,        -63,        -27,        -55,
     :         -5,          0,          1,          2 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=661,670 ) /
     :         -3,         28,         15,          2,
     :          5,          0,          1,         -2,
     :          0,          9,          4,          1,
     :          0,          9,          4,         -1,
     :       -126,        -63,        -27,         55,
     :          3,          0,          0,         -1,
     :         21,        -11,         -6,        -11,
     :          0,         -4,          0,          0,
     :        -21,        -11,         -6,         11,
     :         -3,          0,          0,          1 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=671,680 ) /
     :          0,          3,          1,          0,
     :          8,          0,          0,         -4,
     :         -6,          0,          0,          3,
     :         -3,          0,          0,          1,
     :          3,          0,          0,         -1,
     :         -3,          0,          0,          1,
     :         -5,          0,          0,          2,
     :         24,        -12,         -5,        -11,
     :          0,          3,          1,          0,
     :          0,          3,          1,          0 /
      DATA ( ( ICPL(I,J), I=1,4 ), J=681,687 ) /
     :          0,          3,          2,          0,
     :        -24,        -12,         -5,         10,
     :          4,          0,         -1,         -2,
     :         13,          0,          0,         -6,
     :          7,          0,          0,         -3,
     :          3,          0,          0,         -1,
     :          3,          0,          0,         -1 /

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Interval between fundamental epoch J2000.0 and given date (JC).
      T = ( ( DATE1-DJ0 ) + DATE2 ) / DJC

*  -------------------
*  LUNI-SOLAR NUTATION
*  -------------------

*
*  Fundamental (Delaunay) arguments from Simon et al. (1994)
*

*  Mean anomaly of the Moon.
      EL  = MOD (         485868.249036D0 +
     :            T*( 1717915923.2178D0 +
     :            T*(         31.8792D0 +
     :            T*(          0.051635D0 +
     :            T*(        - 0.00024470D0 )))), TURNAS ) * DAS2R

*  Mean anomaly of the Sun.
      ELP = MOD (        1287104.79305D0 +
     :            T*(  129596581.0481D0 +
     :            T*(        - 0.5532D0 +
     :            T*(          0.000136D0 +
     :            T*(        - 0.00001149D0 )))), TURNAS ) * DAS2R

*  Mean argument of the latitude of the Moon.
      F   = MOD (         335779.526232D0 +
     :            T*( 1739527262.8478D0 +
     :            T*(       - 12.7512D0 +
     :            T*(       -  0.001037D0 +
     :            T*(          0.00000417D0 )))), TURNAS ) * DAS2R

*  Mean elongation of the Moon from the Sun.
      D   = MOD (        1072260.70369D0 +
     :            T*( 1602961601.2090D0 +
     :            T*(        - 6.3706D0 +
     :            T*(          0.006593D0 +
     :            T*(        - 0.00003169D0 )))), TURNAS ) * DAS2R

*  Mean longitude of the ascending node of the Moon.
      OM  = MOD (         450160.398036D0 +
     :            T*(  - 6962890.5431D0 +
     :            T*(          7.4722D0 +
     :            T*(          0.007702D0 +
     :            T*(        - 0.00005939D0 )))), TURNAS ) * DAS2R

*  Initialize the nutation values.
      DP = 0D0
      DE = 0D0

*  Summation of luni-solar nutation series (in reverse order).
      DO 100 I = NLS, 1, -1

*     Argument and functions.
         ARG = MOD ( DBLE ( NALS(1,I) ) * EL  +
     :               DBLE ( NALS(2,I) ) * ELP +
     :               DBLE ( NALS(3,I) ) * F   +
     :               DBLE ( NALS(4,I) ) * D   +
     :               DBLE ( NALS(5,I) ) * OM, D2PI )
         SARG = SIN(ARG)
         CARG = COS(ARG)

*     Term.
         DP = DP + ( CLS(1,I) + CLS(2,I) * T ) * SARG
     :           +   CLS(3,I)                  * CARG
         DE = DE + ( CLS(4,I) + CLS(5,I) * T ) * CARG
     :           +   CLS(6,I)                  * SARG

 100  CONTINUE

*  Convert from 0.1 microarcsec units to radians.
      DPSILS = DP * U2R
      DEPSLS = DE * U2R

*  ------------------
*  PLANETARY NUTATION
*  ------------------

*  Mean anomaly of the Moon.
      AL   = MOD ( 2.35555598D0 + 8328.6914269554D0 * T, D2PI )

*  Mean anomaly of the Sun.
      ALSU = MOD ( 6.24006013D0 + 628.301955D0 * T, D2PI )

*  Mean argument of the latitude of the Moon.
      AF   = MOD ( 1.627905234D0 + 8433.466158131D0 * T, D2PI )

*  Mean elongation of the Moon from the Sun.
      AD   = MOD ( 5.198466741D0 + 7771.3771468121D0 * T, D2PI )

*  Mean longitude of the ascending node of the Moon.
      AOM  = MOD ( 2.18243920D0 - 33.757045D0 * T, D2PI )

*  General accumulated precession in longitude.
      APA  = ( 0.02438175D0 + 0.00000538691D0 * T ) * T

*  Planetary longitudes, Mercury through Neptune (Souchay et al. 1999).
      ALME = MOD ( 4.402608842D0 + 2608.7903141574D0 * T, D2PI )
      ALVE = MOD ( 3.176146697D0 + 1021.3285546211D0 * T, D2PI )
      ALEA = MOD ( 1.753470314D0 +  628.3075849991D0 * T, D2PI )
      ALMA = MOD ( 6.203480913D0 +  334.0612426700D0 * T, D2PI )
      ALJU = MOD ( 0.599546497D0 +   52.9690962641D0 * T, D2PI )
      ALSA = MOD ( 0.874016757D0 +   21.3299104960D0 * T, D2PI )
      ALUR = MOD ( 5.481293871D0 +    7.4781598567D0 * T, D2PI )
      ALNE = MOD ( 5.321159000D0 +    3.8127774000D0 * T, D2PI )

*  Initialize the nutation values.
      DP = 0D0
      DE = 0D0

*  Summation of planetary nutation series (in reverse order).
      DO 200 I = NPL, 1, -1

*     Argument and functions.
         ARG = MOD ( DBLE ( NAPL( 1,I) ) * AL   +
     :               DBLE ( NAPL( 2,I) ) * ALSU +
     :               DBLE ( NAPL( 3,I) ) * AF   +
     :               DBLE ( NAPL( 4,I) ) * AD   +
     :               DBLE ( NAPL( 5,I) ) * AOM  +
     :               DBLE ( NAPL( 6,I) ) * ALME +
     :               DBLE ( NAPL( 7,I) ) * ALVE +
     :               DBLE ( NAPL( 8,I) ) * ALEA +
     :               DBLE ( NAPL( 9,I) ) * ALMA +
     :               DBLE ( NAPL(10,I) ) * ALJU +
     :               DBLE ( NAPL(11,I) ) * ALSA +
     :               DBLE ( NAPL(12,I) ) * ALUR +
     :               DBLE ( NAPL(13,I) ) * ALNE +
     :               DBLE ( NAPL(14,I) ) * APA, D2PI )
         SARG = SIN(ARG)
         CARG = COS(ARG)

*     Term.
         DP = DP + DBLE( ICPL(1,I)) * SARG + DBLE( ICPL(2,I)) * CARG
         DE = DE + DBLE( ICPL(3,I)) * SARG + DBLE( ICPL(4,I)) * CARG

 200  CONTINUE

*  Convert from 0.1 microarcsec units to radians.
      DPSIPL = DP * U2R
      DEPSPL = DE * U2R

*  -----
*  TOTAL
*  -----

*  Add planetary and luni-solar components.
      DPSI = DPSIPL + DPSILS
      DEPS = DEPSPL + DEPSLS

      END



      SUBROUTINE NU2000K ( DATE1, DATE2, DPSI, DEPS )

*+
*  - - - - - - - -
*   N U 2 0 0 0 K
*  - - - - - - - -
*
*  Nutation, IAU 2000A model (MHB_2000 without FCN) MODIFIED.  Series
*     truncated for speed of execution, and using Simon et al. (1994)
*     fundamental arguments throughout.  Accuracy, compared to
*     IAU 2000 A series, is 0.1 mas in delta psi and 0.04 mas in
*     delta epsilon and delta psi sin(epsilon) over 6 centuries
*     centered at year 2000 (99% of errors less than these values).
*
*  Modified form of NU2000A, by Pat Wallace, given in subroutine annex
*  to Chapter 5 of IERS Conventions (2003).
*
*  Given:
*     DATE1,DATE2    d   TT date (JD = DATE1+DATE2)
*
*  Returned:
*     DPSI,DEPS      d   nutation (luni-solar + planetary, radians)
*
*  This revision:  2002 November 25
*                  2004 March 1     (by G. Kaplan)
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION DATE1, DATE2, DPSI, DEPS

*  Arcseconds to radians
      DOUBLE PRECISION DAS2R
      PARAMETER ( DAS2R = 4.848136811095359935899141D-6 )

*  Milliarcseconds to radians
      DOUBLE PRECISION DMAS2R
      PARAMETER ( DMAS2R = DAS2R / 1D3 )

*  Arc seconds in a full circle
      DOUBLE PRECISION TURNAS
      PARAMETER ( TURNAS = 1296000D0 )

*  2Pi
      DOUBLE PRECISION D2PI
      PARAMETER ( D2PI = 6.283185307179586476925287D0 )

*  Units of 0.1 microarcsecond to radians
      DOUBLE PRECISION U2R
      PARAMETER ( U2R = DAS2R/1D7 )

*  Reference epoch (J2000), JD
      DOUBLE PRECISION DJ0
      PARAMETER ( DJ0 = 2451545D0 )

*  Days per Julian century
      DOUBLE PRECISION DJC
      PARAMETER ( DJC = 36525D0 )

*  Miscellaneous
      DOUBLE PRECISION T, EL, ELP, F, D, OM, ARG, DP, DE, SARG, CARG,
     :                 DPSILS, DEPSLS,
     :                 ALME, ALVE, ALEA, ALMA, ALJU, ALSA, ALUR, ALNE,
     :                 APA,
     :                 DPSIPL, DEPSPL
      INTEGER I, J

*  -------------------------
*  Luni-Solar nutation model
*  -------------------------

*  Number of terms in the luni-solar nutation model
      INTEGER NLS
      PARAMETER ( NLS = 323 )

*  Coefficients for fundamental arguments
      INTEGER NALS(5,NLS)

*  Longitude and obliquity coefficients
      DOUBLE PRECISION CLS(6,NLS)

*  ---------------
*  Planetary terms
*  ---------------

*  Number of terms in the planetary nutation model
      INTEGER NPL
      PARAMETER ( NPL = 165 )

*  Coefficients for fundamental arguments
      INTEGER NAPL(14,NPL)

*  Longitude and obliquity coefficients
      DOUBLE PRECISION CPL(4,NPL)

*  ----------------------------------------
*  Tables of argument and term coefficients
*  ----------------------------------------

*
*  Luni-Solar argument multipliers:
*               L     L'    F     D     Om
*
      DATA ( ( NALS(I,J), I=1,5 ), J=  1, 20 ) /
     :          0,    0,    0,    0,    1,
     :          0,    0,    2,   -2,    2,
     :          0,    0,    2,    0,    2,
     :          0,    0,    0,    0,    2,
     :          0,    1,    0,    0,    0,
     :          0,    1,    2,   -2,    2,
     :          1,    0,    0,    0,    0,
     :          0,    0,    2,    0,    1,
     :          1,    0,    2,    0,    2,
     :          0,   -1,    2,   -2,    2,
     :          0,    0,    2,   -2,    1,
     :         -1,    0,    2,    0,    2,
     :         -1,    0,    0,    2,    0,
     :          1,    0,    0,    0,    1,
     :         -1,    0,    0,    0,    1,
     :         -1,    0,    2,    2,    2,
     :          1,    0,    2,    0,    1,
     :         -2,    0,    2,    0,    1,
     :          0,    0,    0,    2,    0,
     :          0,    0,    2,    2,    2/
      DATA ( ( NALS(I,J), I=1,5 ), J= 21, 40 ) /
     :          0,   -2,    2,   -2,    2,
     :         -2,    0,    0,    2,    0,
     :          2,    0,    2,    0,    2,
     :          1,    0,    2,   -2,    2,
     :         -1,    0,    2,    0,    1,
     :          2,    0,    0,    0,    0,
     :          0,    0,    2,    0,    0,
     :          0,    1,    0,    0,    1,
     :         -1,    0,    0,    2,    1,
     :          0,    2,    2,   -2,    2,
     :          0,    0,   -2,    2,    0,
     :          1,    0,    0,   -2,    1,
     :          0,   -1,    0,    0,    1,
     :         -1,    0,    2,    2,    1,
     :          0,    2,    0,    0,    0,
     :          1,    0,    2,    2,    2,
     :         -2,    0,    2,    0,    0,
     :          0,    1,    2,    0,    2,
     :          0,    0,    2,    2,    1,
     :          0,   -1,    2,    0,    2/
      DATA ( ( NALS(I,J), I=1,5 ), J= 41, 60 ) /
     :          0,    0,    0,    2,    1,
     :          1,    0,    2,   -2,    1,
     :          2,    0,    2,   -2,    2,
     :         -2,    0,    0,    2,    1,
     :          2,    0,    2,    0,    1,
     :          0,   -1,    2,   -2,    1,
     :          0,    0,    0,   -2,    1,
     :         -1,   -1,    0,    2,    0,
     :          2,    0,    0,   -2,    1,
     :          1,    0,    0,    2,    0,
     :          0,    1,    2,   -2,    1,
     :          1,   -1,    0,    0,    0,
     :         -2,    0,    2,    0,    2,
     :          3,    0,    2,    0,    2,
     :          0,   -1,    0,    2,    0,
     :          1,   -1,    2,    0,    2,
     :          0,    0,    0,    1,    0,
     :         -1,   -1,    2,    2,    2,
     :         -1,    0,    2,    0,    0,
     :          0,   -1,    2,    2,    2/
      DATA ( ( NALS(I,J), I=1,5 ), J= 61, 80 ) /
     :         -2,    0,    0,    0,    1,
     :          1,    1,    2,    0,    2,
     :          2,    0,    0,    0,    1,
     :         -1,    1,    0,    1,    0,
     :          1,    1,    0,    0,    0,
     :          1,    0,    2,    0,    0,
     :         -1,    0,    2,   -2,    1,
     :          1,    0,    0,    0,    2,
     :         -1,    0,    0,    1,    0,
     :          0,    0,    2,    1,    2,
     :         -1,    0,    2,    4,    2,
     :         -1,    1,    0,    1,    1,
     :          0,   -2,    2,   -2,    1,
     :          1,    0,    2,    2,    1,
     :         -2,    0,    2,    2,    2,
     :         -1,    0,    0,    0,    2,
     :          1,    1,    2,   -2,    2,
     :         -2,    0,    2,    4,    2,
     :         -1,    0,    4,    0,    2,
     :          2,    0,    2,   -2,    1/
      DATA ( ( NALS(I,J), I=1,5 ), J= 81,100 ) /
     :          2,    0,    2,    2,    2,
     :          1,    0,    0,    2,    1,
     :          3,    0,    0,    0,    0,
     :          3,    0,    2,   -2,    2,
     :          0,    0,    4,   -2,    2,
     :          0,    1,    2,    0,    1,
     :          0,    0,   -2,    2,    1,
     :          0,    0,    2,   -2,    3,
     :         -1,    0,    0,    4,    0,
     :          2,    0,   -2,    0,    1,
     :         -2,    0,    0,    4,    0,
     :         -1,   -1,    0,    2,    1,
     :         -1,    0,    0,    1,    1,
     :          0,    1,    0,    0,    2,
     :          0,    0,   -2,    0,    1,
     :          0,   -1,    2,    0,    1,
     :          0,    0,    2,   -1,    2,
     :          0,    0,    2,    4,    2,
     :         -2,   -1,    0,    2,    0,
     :          1,    1,    0,   -2,    1/
      DATA ( ( NALS(I,J), I=1,5 ), J=101,120 ) /
     :         -1,    1,    0,    2,    0,
     :         -1,    1,    0,    1,    2,
     :          1,   -1,    0,    0,    1,
     :          1,   -1,    2,    2,    2,
     :         -1,    1,    2,    2,    2,
     :          3,    0,    2,    0,    1,
     :          0,    1,   -2,    2,    0,
     :         -1,    0,    0,   -2,    1,
     :          0,    1,    2,    2,    2,
     :         -1,   -1,    2,    2,    1,
     :          0,   -1,    0,    0,    2,
     :          1,    0,    2,   -4,    1,
     :         -1,    0,   -2,    2,    0,
     :          0,   -1,    2,    2,    1,
     :          2,   -1,    2,    0,    2,
     :          0,    0,    0,    2,    2,
     :          1,   -1,    2,    0,    1,
     :         -1,    1,    2,    0,    2,
     :          0,    1,    0,    2,    0,
     :          0,   -1,   -2,    2,    0/
      DATA ( ( NALS(I,J), I=1,5 ), J=121,140 ) /
     :          0,    3,    2,   -2,    2,
     :          0,    0,    0,    1,    1,
     :         -1,    0,    2,    2,    0,
     :          2,    1,    2,    0,    2,
     :          1,    1,    0,    0,    1,
     :          1,    1,    2,    0,    1,
     :          2,    0,    0,    2,    0,
     :          1,    0,   -2,    2,    0,
     :         -1,    0,    0,    2,    2,
     :          0,    1,    0,    1,    0,
     :          0,    1,    0,   -2,    1,
     :         -1,    0,    2,   -2,    2,
     :          0,    0,    0,   -1,    1,
     :         -1,    1,    0,    0,    1,
     :          1,    0,    2,   -1,    2,
     :          1,   -1,    0,    2,    0,
     :          0,    0,    0,    4,    0,
     :          1,    0,    2,    1,    2,
     :          0,    0,    2,    1,    1,
     :          1,    0,    0,   -2,    2/
      DATA ( ( NALS(I,J), I=1,5 ), J=141,160 ) /
     :         -1,    0,    2,    4,    1,
     :          1,    0,   -2,    0,    1,
     :          1,    1,    2,   -2,    1,
     :          0,    0,    2,    2,    0,
     :         -1,    0,    2,   -1,    1,
     :         -2,    0,    2,    2,    1,
     :          4,    0,    2,    0,    2,
     :          2,   -1,    0,    0,    0,
     :          2,    1,    2,   -2,    2,
     :          0,    1,    2,    1,    2,
     :          1,    0,    4,   -2,    2,
     :         -1,   -1,    0,    0,    1,
     :          0,    1,    0,    2,    1,
     :         -2,    0,    2,    4,    1,
     :          2,    0,    2,    0,    0,
     :          1,    0,    0,    1,    0,
     :         -1,    0,    0,    4,    1,
     :         -1,    0,    4,    0,    1,
     :          2,    0,    2,    2,    1,
     :          0,    0,    2,   -3,    2/
      DATA ( ( NALS(I,J), I=1,5 ), J=161,180 ) /
     :         -1,   -2,    0,    2,    0,
     :          2,    1,    0,    0,    0,
     :          0,    0,    4,    0,    2,
     :          0,    0,    0,    0,    3,
     :          0,    3,    0,    0,    0,
     :          0,    0,    2,   -4,    1,
     :          0,   -1,    0,    2,    1,
     :          0,    0,    0,    4,    1,
     :         -1,   -1,    2,    4,    2,
     :          1,    0,    2,    4,    2,
     :         -2,    2,    0,    2,    0,
     :         -2,   -1,    2,    0,    1,
     :         -2,    0,    0,    2,    2,
     :         -1,   -1,    2,    0,    2,
     :          0,    0,    4,   -2,    1,
     :          3,    0,    2,   -2,    1,
     :         -2,   -1,    0,    2,    1,
     :          1,    0,    0,   -1,    1,
     :          0,   -2,    0,    2,    0,
     :         -2,    0,    0,    4,    1/
      DATA ( ( NALS(I,J), I=1,5 ), J=181,200 ) /
     :         -3,    0,    0,    0,    1,
     :          1,    1,    2,    2,    2,
     :          0,    0,    2,    4,    1,
     :          3,    0,    2,    2,    2,
     :         -1,    1,    2,   -2,    1,
     :          2,    0,    0,   -4,    1,
     :          0,    0,    0,   -2,    2,
     :          2,    0,    2,   -4,    1,
     :         -1,    1,    0,    2,    1,
     :          0,    0,    2,   -1,    1,
     :          0,   -2,    2,    2,    2,
     :          2,    0,    0,    2,    1,
     :          4,    0,    2,   -2,    2,
     :          2,    0,    0,   -2,    2,
     :          0,    2,    0,    0,    1,
     :          1,    0,    0,   -4,    1,
     :          0,    2,    2,   -2,    1,
     :         -3,    0,    0,    4,    0,
     :         -1,    1,    2,    0,    1,
     :         -1,   -1,    0,    4,    0/
      DATA ( ( NALS(I,J), I=1,5 ), J=201,220 ) /
     :         -1,   -2,    2,    2,    2,
     :         -2,   -1,    2,    4,    2,
     :          1,   -1,    2,    2,    1,
     :         -2,    1,    0,    2,    0,
     :         -2,    1,    2,    0,    1,
     :          2,    1,    0,   -2,    1,
     :         -3,    0,    2,    0,    1,
     :         -2,    0,    2,   -2,    1,
     :         -1,    1,    0,    2,    2,
     :          0,   -1,    2,   -1,    2,
     :         -1,    0,    4,   -2,    2,
     :          0,   -2,    2,    0,    2,
     :         -1,    0,    2,    1,    2,
     :          2,    0,    0,    0,    2,
     :          0,    0,    2,    0,    3,
     :         -2,    0,    4,    0,    2,
     :         -1,    0,   -2,    0,    1,
     :         -1,    1,    2,    2,    1,
     :          3,    0,    0,    0,    1,
     :         -1,    0,    2,    3,    2/
      DATA ( ( NALS(I,J), I=1,5 ), J=221,240 ) /
     :          2,   -1,    2,    0,    1,
     :          0,    1,    2,    2,    1,
     :          0,   -1,    2,    4,    2,
     :          2,   -1,    2,    2,    2,
     :          0,    2,   -2,    2,    0,
     :         -1,   -1,    2,   -1,    1,
     :          0,   -2,    0,    0,    1,
     :          1,    0,    2,   -4,    2,
     :          1,   -1,    0,   -2,    1,
     :         -1,   -1,    2,    0,    1,
     :          1,   -1,    2,   -2,    2,
     :         -2,   -1,    0,    4,    0,
     :         -1,    0,    0,    3,    0,
     :         -2,   -1,    2,    2,    2,
     :          0,    2,    2,    0,    2,
     :          1,    1,    0,    2,    0,
     :          2,    0,    2,   -1,    2,
     :          1,    0,    2,    1,    1,
     :          4,    0,    0,    0,    0,
     :          2,    1,    2,    0,    1/
      DATA ( ( NALS(I,J), I=1,5 ), J=241,260 ) /
     :          3,   -1,    2,    0,    2,
     :         -2,    2,    0,    2,    1,
     :          1,    0,    2,   -3,    1,
     :          1,    1,    2,   -4,    1,
     :         -1,   -1,    2,   -2,    1,
     :          0,   -1,    0,   -1,    1,
     :          0,   -1,    0,   -2,    1,
     :         -2,    0,    0,    0,    2,
     :         -2,    0,   -2,    2,    0,
     :         -1,    0,   -2,    4,    0,
     :          1,   -2,    0,    0,    0,
     :          0,    1,    0,    1,    1,
     :         -1,    2,    0,    2,    0,
     :          1,   -1,    2,   -2,    1,
     :          1,    2,    2,   -2,    2,
     :          2,   -1,    2,   -2,    2,
     :          1,    0,    2,   -1,    1,
     :          2,    1,    2,   -2,    1,
     :         -2,    0,    0,   -2,    1,
     :          1,   -2,    2,    0,    2/
      DATA ( ( NALS(I,J), I=1,5 ), J=261,280 ) /
     :          0,    1,    2,    1,    1,
     :          1,    0,    4,   -2,    1,
     :         -2,    0,    4,    2,    2,
     :          1,    1,    2,    1,    2,
     :          1,    0,    0,    4,    0,
     :          1,    0,    2,    2,    0,
     :          2,    0,    2,    1,    2,
     :          3,    1,    2,    0,    2,
     :          4,    0,    2,    0,    1,
     :         -2,   -1,    2,    0,    0,
     :          0,    1,   -2,    2,    1,
     :          1,    0,   -2,    1,    0,
     :          2,   -1,    0,   -2,    1,
     :         -1,    0,    2,   -1,    2,
     :          1,    0,    2,   -3,    2,
     :          0,    1,    2,   -2,    3,
     :         -1,    0,   -2,    2,    1,
     :          0,    0,    2,   -4,    2,
     :          2,    0,    2,   -4,    2,
     :          0,    0,    4,   -4,    4/
      DATA ( ( NALS(I,J), I=1,5 ), J=281,300 ) /
     :          0,    0,    4,   -4,    2,
     :         -2,    0,    0,    3,    0,
     :          1,    0,   -2,    2,    1,
     :         -3,    0,    2,    2,    2,
     :         -2,    0,    2,    2,    0,
     :          2,   -1,    0,    0,    1,
     :          1,    1,    0,    1,    0,
     :          0,    1,    4,   -2,    2,
     :         -1,    1,    0,   -2,    1,
     :          0,    0,    0,   -4,    1,
     :          1,   -1,    0,    2,    1,
     :          1,    1,    0,    2,    1,
     :         -1,    2,    2,    2,    2,
     :          3,    1,    2,   -2,    2,
     :          0,   -1,    0,    4,    0,
     :          2,   -1,    0,    2,    0,
     :          0,    0,    4,    0,    1,
     :          2,    0,    4,   -2,    2,
     :         -1,   -1,    2,    4,    1,
     :          1,    0,    0,    4,    1/
      DATA ( ( NALS(I,J), I=1,5 ), J=301,320 ) /
     :          1,   -2,    2,    2,    2,
     :          0,    0,    2,    3,    2,
     :         -1,    1,    2,    4,    2,
     :          3,    0,    0,    2,    0,
     :         -1,    0,    4,    2,    2,
     :         -2,    0,    2,    6,    2,
     :         -1,    0,    2,    6,    2,
     :          1,    1,   -2,    1,    0,
     :         -1,    0,    0,    1,    2,
     :         -1,   -1,    0,    1,    0,
     :         -2,    0,    0,    1,    0,
     :          0,    0,   -2,    1,    0,
     :          1,   -1,   -2,    2,    0,
     :          1,    2,    0,    0,    0,
     :          3,    0,    2,    0,    0,
     :          0,   -1,    1,   -1,    1,
     :         -1,    0,    1,    0,    3,
     :         -1,    0,    1,    0,    2,
     :         -1,    0,    1,    0,    1,
     :         -1,    0,    1,    0,    0/
      DATA ( ( NALS(I,J), I=1,5 ), J=321,323 ) /
     :          0,    0,    1,    0,    2,
     :          0,    0,    1,    0,    1,
     :          0,    0,    1,    0,    0/

*
*  Luni-Solar nutation coefficients, unit 1e-7 arcsec:
*  longitude (sin, t*sin, cos), obliquity (cos, t*cos, sin)
*
      DATA ( ( CLS(I,J), I=1,6 ), J=  1, 20 ) /
     :-172064161.D0,-174666.D0, 33386.D0, 92052331.D0, 9086.D0,15377.D0,
     : -13170906.D0,  -1675.D0,-13696.D0,  5730336.D0,-3015.D0,-4587.D0,
     :  -2276413.D0,   -234.D0,  2796.D0,   978459.D0, -485.D0, 1374.D0,
     :   2074554.D0,    207.D0,  -698.D0,  -897492.D0,  470.D0, -291.D0,
     :   1475877.D0,  -3633.D0, 11817.D0,    73871.D0, -184.D0,-1924.D0,
     :   -516821.D0,   1226.D0,  -524.D0,   224386.D0, -677.D0, -174.D0,
     :    711159.D0,     73.D0,  -872.D0,    -6750.D0,    0.D0,  358.D0,
     :   -387298.D0,   -367.D0,   380.D0,   200728.D0,   18.D0,  318.D0,
     :   -301461.D0,    -36.D0,   816.D0,   129025.D0,  -63.D0,  367.D0,
     :    215829.D0,   -494.D0,   111.D0,   -95929.D0,  299.D0,  132.D0,
     :    128227.D0,    137.D0,   181.D0,   -68982.D0,   -9.D0,   39.D0,
     :    123457.D0,     11.D0,    19.D0,   -53311.D0,   32.D0,   -4.D0,
     :    156994.D0,     10.D0,  -168.D0,    -1235.D0,    0.D0,   82.D0,
     :     63110.D0,     63.D0,    27.D0,   -33228.D0,    0.D0,   -9.D0,
     :    -57976.D0,    -63.D0,  -189.D0,    31429.D0,    0.D0,  -75.D0,
     :    -59641.D0,    -11.D0,   149.D0,    25543.D0,  -11.D0,   66.D0,
     :    -51613.D0,    -42.D0,   129.D0,    26366.D0,    0.D0,   78.D0,
     :     45893.D0,     50.D0,    31.D0,   -24236.D0,  -10.D0,   20.D0,
     :     63384.D0,     11.D0,  -150.D0,    -1220.D0,    0.D0,   29.D0,
     :    -38571.D0,     -1.D0,   158.D0,    16452.D0,  -11.D0,   68.D0/
      DATA ( ( CLS(I,J), I=1,6 ), J= 21, 40 ) /
     :     32481.D0,      0.D0,     0.D0,   -13870.D0,    0.D0,    0.D0,
     :    -47722.D0,      0.D0,   -18.D0,      477.D0,    0.D0,  -25.D0,
     :    -31046.D0,     -1.D0,   131.D0,    13238.D0,  -11.D0,   59.D0,
     :     28593.D0,      0.D0,    -1.D0,   -12338.D0,   10.D0,   -3.D0,
     :     20441.D0,     21.D0,    10.D0,   -10758.D0,    0.D0,   -3.D0,
     :     29243.D0,      0.D0,   -74.D0,     -609.D0,    0.D0,   13.D0,
     :     25887.D0,      0.D0,   -66.D0,     -550.D0,    0.D0,   11.D0,
     :    -14053.D0,    -25.D0,    79.D0,     8551.D0,   -2.D0,  -45.D0,
     :     15164.D0,     10.D0,    11.D0,    -8001.D0,    0.D0,   -1.D0,
     :    -15794.D0,     72.D0,   -16.D0,     6850.D0,  -42.D0,   -5.D0,
     :     21783.D0,      0.D0,    13.D0,     -167.D0,    0.D0,   13.D0,
     :    -12873.D0,    -10.D0,   -37.D0,     6953.D0,    0.D0,  -14.D0,
     :    -12654.D0,     11.D0,    63.D0,     6415.D0,    0.D0,   26.D0,
     :    -10204.D0,      0.D0,    25.D0,     5222.D0,    0.D0,   15.D0,
     :     16707.D0,    -85.D0,   -10.D0,      168.D0,   -1.D0,   10.D0,
     :     -7691.D0,      0.D0,    44.D0,     3268.D0,    0.D0,   19.D0,
     :    -11024.D0,      0.D0,   -14.D0,      104.D0,    0.D0,    2.D0,
     :      7566.D0,    -21.D0,   -11.D0,    -3250.D0,    0.D0,   -5.D0,
     :     -6637.D0,    -11.D0,    25.D0,     3353.D0,    0.D0,   14.D0,
     :     -7141.D0,     21.D0,     8.D0,     3070.D0,    0.D0,    4.D0/
      DATA ( ( CLS(I,J), I=1,6 ), J= 41, 60 ) /
     :     -6302.D0,    -11.D0,     2.D0,     3272.D0,    0.D0,    4.D0,
     :      5800.D0,     10.D0,     2.D0,    -3045.D0,    0.D0,   -1.D0,
     :      6443.D0,      0.D0,    -7.D0,    -2768.D0,    0.D0,   -4.D0,
     :     -5774.D0,    -11.D0,   -15.D0,     3041.D0,    0.D0,   -5.D0,
     :     -5350.D0,      0.D0,    21.D0,     2695.D0,    0.D0,   12.D0,
     :     -4752.D0,    -11.D0,    -3.D0,     2719.D0,    0.D0,   -3.D0,
     :     -4940.D0,    -11.D0,   -21.D0,     2720.D0,    0.D0,   -9.D0,
     :      7350.D0,      0.D0,    -8.D0,      -51.D0,    0.D0,    4.D0,
     :      4065.D0,      0.D0,     6.D0,    -2206.D0,    0.D0,    1.D0,
     :      6579.D0,      0.D0,   -24.D0,     -199.D0,    0.D0,    2.D0,
     :      3579.D0,      0.D0,     5.D0,    -1900.D0,    0.D0,    1.D0,
     :      4725.D0,      0.D0,    -6.D0,      -41.D0,    0.D0,    3.D0,
     :     -3075.D0,      0.D0,    -2.D0,     1313.D0,    0.D0,   -1.D0,
     :     -2904.D0,      0.D0,    15.D0,     1233.D0,    0.D0,    7.D0,
     :      4348.D0,      0.D0,   -10.D0,      -81.D0,    0.D0,    2.D0,
     :     -2878.D0,      0.D0,     8.D0,     1232.D0,    0.D0,    4.D0,
     :     -4230.D0,      0.D0,     5.D0,      -20.D0,    0.D0,   -2.D0,
     :     -2819.D0,      0.D0,     7.D0,     1207.D0,    0.D0,    3.D0,
     :     -4056.D0,      0.D0,     5.D0,       40.D0,    0.D0,   -2.D0,
     :     -2647.D0,      0.D0,    11.D0,     1129.D0,    0.D0,    5.D0/
      DATA ( ( CLS(I,J), I=1,6 ), J= 61, 80 ) /
     :     -2294.D0,      0.D0,   -10.D0,     1266.D0,    0.D0,   -4.D0,
     :      2481.D0,      0.D0,    -7.D0,    -1062.D0,    0.D0,   -3.D0,
     :      2179.D0,      0.D0,    -2.D0,    -1129.D0,    0.D0,   -2.D0,
     :      3276.D0,      0.D0,     1.D0,       -9.D0,    0.D0,    0.D0,
     :     -3389.D0,      0.D0,     5.D0,       35.D0,    0.D0,   -2.D0,
     :      3339.D0,      0.D0,   -13.D0,     -107.D0,    0.D0,    1.D0,
     :     -1987.D0,      0.D0,    -6.D0,     1073.D0,    0.D0,   -2.D0,
     :     -1981.D0,      0.D0,     0.D0,      854.D0,    0.D0,    0.D0,
     :      4026.D0,      0.D0,  -353.D0,     -553.D0,    0.D0, -139.D0,
     :      1660.D0,      0.D0,    -5.D0,     -710.D0,    0.D0,   -2.D0,
     :     -1521.D0,      0.D0,     9.D0,      647.D0,    0.D0,    4.D0,
     :      1314.D0,      0.D0,     0.D0,     -700.D0,    0.D0,    0.D0,
     :     -1283.D0,      0.D0,     0.D0,      672.D0,    0.D0,    0.D0,
     :     -1331.D0,      0.D0,     8.D0,      663.D0,    0.D0,    4.D0,
     :      1383.D0,      0.D0,    -2.D0,     -594.D0,    0.D0,   -2.D0,
     :      1405.D0,      0.D0,     4.D0,     -610.D0,    0.D0,    2.D0,
     :      1290.D0,      0.D0,     0.D0,     -556.D0,    0.D0,    0.D0,
     :     -1214.D0,      0.D0,     5.D0,      518.D0,    0.D0,    2.D0,
     :      1146.D0,      0.D0,    -3.D0,     -490.D0,    0.D0,   -1.D0,
     :      1019.D0,      0.D0,    -1.D0,     -527.D0,    0.D0,   -1.D0/
      DATA ( ( CLS(I,J), I=1,6 ), J= 81,100 ) /
     :     -1100.D0,      0.D0,     9.D0,      465.D0,    0.D0,    4.D0,
     :      -970.D0,      0.D0,     2.D0,      496.D0,    0.D0,    1.D0,
     :      1575.D0,      0.D0,    -6.D0,      -50.D0,    0.D0,    0.D0,
     :       934.D0,      0.D0,    -3.D0,     -399.D0,    0.D0,   -1.D0,
     :       922.D0,      0.D0,    -1.D0,     -395.D0,    0.D0,   -1.D0,
     :       815.D0,      0.D0,    -1.D0,     -422.D0,    0.D0,   -1.D0,
     :       834.D0,      0.D0,     2.D0,     -440.D0,    0.D0,    1.D0,
     :      1248.D0,      0.D0,     0.D0,     -170.D0,    0.D0,    1.D0,
     :      1338.D0,      0.D0,    -5.D0,      -39.D0,    0.D0,    0.D0,
     :       716.D0,      0.D0,    -2.D0,     -389.D0,    0.D0,   -1.D0,
     :      1282.D0,      0.D0,    -3.D0,      -23.D0,    0.D0,    1.D0,
     :       742.D0,      0.D0,     1.D0,     -391.D0,    0.D0,    0.D0,
     :      1020.D0,      0.D0,   -25.D0,     -495.D0,    0.D0,  -10.D0,
     :       715.D0,      0.D0,    -4.D0,     -326.D0,    0.D0,    2.D0,
     :      -666.D0,      0.D0,    -3.D0,      369.D0,    0.D0,   -1.D0,
     :      -667.D0,      0.D0,     1.D0,      346.D0,    0.D0,    1.D0,
     :      -704.D0,      0.D0,     0.D0,      304.D0,    0.D0,    0.D0,
     :      -694.D0,      0.D0,     5.D0,      294.D0,    0.D0,    2.D0,
     :     -1014.D0,      0.D0,    -1.D0,        4.D0,    0.D0,   -1.D0,
     :      -585.D0,      0.D0,    -2.D0,      316.D0,    0.D0,   -1.D0/
      DATA ( ( CLS(I,J), I=1,6 ), J=101,120 ) /
     :      -949.D0,      0.D0,     1.D0,        8.D0,    0.D0,   -1.D0,
     :      -595.D0,      0.D0,     0.D0,      258.D0,    0.D0,    0.D0,
     :       528.D0,      0.D0,     0.D0,     -279.D0,    0.D0,    0.D0,
     :      -590.D0,      0.D0,     4.D0,      252.D0,    0.D0,    2.D0,
     :       570.D0,      0.D0,    -2.D0,     -244.D0,    0.D0,   -1.D0,
     :      -502.D0,      0.D0,     3.D0,      250.D0,    0.D0,    2.D0,
     :      -875.D0,      0.D0,     1.D0,       29.D0,    0.D0,    0.D0,
     :      -492.D0,      0.D0,    -3.D0,      275.D0,    0.D0,   -1.D0,
     :       535.D0,      0.D0,    -2.D0,     -228.D0,    0.D0,   -1.D0,
     :      -467.D0,      0.D0,     1.D0,      240.D0,    0.D0,    1.D0,
     :       591.D0,      0.D0,     0.D0,     -253.D0,    0.D0,    0.D0,
     :      -453.D0,      0.D0,    -1.D0,      244.D0,    0.D0,   -1.D0,
     :       766.D0,      0.D0,     1.D0,        9.D0,    0.D0,    0.D0,
     :      -446.D0,      0.D0,     2.D0,      225.D0,    0.D0,    1.D0,
     :      -488.D0,      0.D0,     2.D0,      207.D0,    0.D0,    1.D0,
     :      -468.D0,      0.D0,     0.D0,      201.D0,    0.D0,    0.D0,
     :      -421.D0,      0.D0,     1.D0,      216.D0,    0.D0,    1.D0,
     :       463.D0,      0.D0,     0.D0,     -200.D0,    0.D0,    0.D0,
     :      -673.D0,      0.D0,     2.D0,       14.D0,    0.D0,    0.D0,
     :       658.D0,      0.D0,     0.D0,       -2.D0,    0.D0,    0.D0/
      DATA ( ( CLS(I,J), I=1,6 ), J=121,140 ) /
     :      -438.D0,      0.D0,     0.D0,      188.D0,    0.D0,    0.D0,
     :      -390.D0,      0.D0,     0.D0,      205.D0,    0.D0,    0.D0,
     :       639.D0,    -11.D0,    -2.D0,      -19.D0,    0.D0,    0.D0,
     :       412.D0,      0.D0,    -2.D0,     -176.D0,    0.D0,   -1.D0,
     :      -361.D0,      0.D0,     0.D0,      189.D0,    0.D0,    0.D0,
     :       360.D0,      0.D0,    -1.D0,     -185.D0,    0.D0,   -1.D0,
     :       588.D0,      0.D0,    -3.D0,      -24.D0,    0.D0,    0.D0,
     :      -578.D0,      0.D0,     1.D0,        5.D0,    0.D0,    0.D0,
     :      -396.D0,      0.D0,     0.D0,      171.D0,    0.D0,    0.D0,
     :       565.D0,      0.D0,    -1.D0,       -6.D0,    0.D0,    0.D0,
     :      -335.D0,      0.D0,    -1.D0,      184.D0,    0.D0,   -1.D0,
     :       357.D0,      0.D0,     1.D0,     -154.D0,    0.D0,    0.D0,
     :       321.D0,      0.D0,     1.D0,     -174.D0,    0.D0,    0.D0,
     :      -301.D0,      0.D0,    -1.D0,      162.D0,    0.D0,    0.D0,
     :      -334.D0,      0.D0,     0.D0,      144.D0,    0.D0,    0.D0,
     :       493.D0,      0.D0,    -2.D0,      -15.D0,    0.D0,    0.D0,
     :       494.D0,      0.D0,    -2.D0,      -19.D0,    0.D0,    0.D0,
     :       337.D0,      0.D0,    -1.D0,     -143.D0,    0.D0,   -1.D0,
     :       280.D0,      0.D0,    -1.D0,     -144.D0,    0.D0,    0.D0,
     :       309.D0,      0.D0,     1.D0,     -134.D0,    0.D0,    0.D0/
      DATA ( ( CLS(I,J), I=1,6 ), J=141,160 ) /
     :      -263.D0,      0.D0,     2.D0,      131.D0,    0.D0,    1.D0,
     :       253.D0,      0.D0,     1.D0,     -138.D0,    0.D0,    0.D0,
     :       245.D0,      0.D0,     0.D0,     -128.D0,    0.D0,    0.D0,
     :       416.D0,      0.D0,    -2.D0,      -17.D0,    0.D0,    0.D0,
     :      -229.D0,      0.D0,     0.D0,      128.D0,    0.D0,    0.D0,
     :       231.D0,      0.D0,     0.D0,     -120.D0,    0.D0,    0.D0,
     :      -259.D0,      0.D0,     2.D0,      109.D0,    0.D0,    1.D0,
     :       375.D0,      0.D0,    -1.D0,       -8.D0,    0.D0,    0.D0,
     :       252.D0,      0.D0,     0.D0,     -108.D0,    0.D0,    0.D0,
     :      -245.D0,      0.D0,     1.D0,      104.D0,    0.D0,    0.D0,
     :       243.D0,      0.D0,    -1.D0,     -104.D0,    0.D0,    0.D0,
     :       208.D0,      0.D0,     1.D0,     -112.D0,    0.D0,    0.D0,
     :       199.D0,      0.D0,     0.D0,     -102.D0,    0.D0,    0.D0,
     :      -208.D0,      0.D0,     1.D0,      105.D0,    0.D0,    0.D0,
     :       335.D0,      0.D0,    -2.D0,      -14.D0,    0.D0,    0.D0,
     :      -325.D0,      0.D0,     1.D0,        7.D0,    0.D0,    0.D0,
     :      -187.D0,      0.D0,     0.D0,       96.D0,    0.D0,    0.D0,
     :       197.D0,      0.D0,    -1.D0,     -100.D0,    0.D0,    0.D0,
     :      -192.D0,      0.D0,     2.D0,       94.D0,    0.D0,    1.D0,
     :      -188.D0,      0.D0,     0.D0,       83.D0,    0.D0,    0.D0/
      DATA ( ( CLS(I,J), I=1,6 ), J=161,180 ) /
     :       276.D0,      0.D0,     0.D0,       -2.D0,    0.D0,    0.D0,
     :      -286.D0,      0.D0,     1.D0,        6.D0,    0.D0,    0.D0,
     :       186.D0,      0.D0,    -1.D0,      -79.D0,    0.D0,    0.D0,
     :      -219.D0,      0.D0,     0.D0,       43.D0,    0.D0,    0.D0,
     :       276.D0,      0.D0,     0.D0,        2.D0,    0.D0,    0.D0,
     :      -153.D0,      0.D0,    -1.D0,       84.D0,    0.D0,    0.D0,
     :      -156.D0,      0.D0,     0.D0,       81.D0,    0.D0,    0.D0,
     :      -154.D0,      0.D0,     1.D0,       78.D0,    0.D0,    0.D0,
     :      -174.D0,      0.D0,     1.D0,       75.D0,    0.D0,    0.D0,
     :      -163.D0,      0.D0,     2.D0,       69.D0,    0.D0,    1.D0,
     :      -228.D0,      0.D0,     0.D0,        1.D0,    0.D0,    0.D0,
     :        91.D0,      0.D0,    -4.D0,      -54.D0,    0.D0,   -2.D0,
     :       175.D0,      0.D0,     0.D0,      -75.D0,    0.D0,    0.D0,
     :      -159.D0,      0.D0,     0.D0,       69.D0,    0.D0,    0.D0,
     :       141.D0,      0.D0,     0.D0,      -72.D0,    0.D0,    0.D0,
     :       147.D0,      0.D0,     0.D0,      -75.D0,    0.D0,    0.D0,
     :      -132.D0,      0.D0,     0.D0,       69.D0,    0.D0,    0.D0,
     :       159.D0,      0.D0,   -28.D0,      -54.D0,    0.D0,   11.D0,
     :       213.D0,      0.D0,     0.D0,       -4.D0,    0.D0,    0.D0,
     :       123.D0,      0.D0,     0.D0,      -64.D0,    0.D0,    0.D0/
      DATA ( ( CLS(I,J), I=1,6 ), J=181,200 ) /
     :      -118.D0,      0.D0,    -1.D0,       66.D0,    0.D0,    0.D0,
     :       144.D0,      0.D0,    -1.D0,      -61.D0,    0.D0,    0.D0,
     :      -121.D0,      0.D0,     1.D0,       60.D0,    0.D0,    0.D0,
     :      -134.D0,      0.D0,     1.D0,       56.D0,    0.D0,    1.D0,
     :      -105.D0,      0.D0,     0.D0,       57.D0,    0.D0,    0.D0,
     :      -102.D0,      0.D0,     0.D0,       56.D0,    0.D0,    0.D0,
     :       120.D0,      0.D0,     0.D0,      -52.D0,    0.D0,    0.D0,
     :       101.D0,      0.D0,     0.D0,      -54.D0,    0.D0,    0.D0,
     :      -113.D0,      0.D0,     0.D0,       59.D0,    0.D0,    0.D0,
     :      -106.D0,      0.D0,     0.D0,       61.D0,    0.D0,    0.D0,
     :      -129.D0,      0.D0,     1.D0,       55.D0,    0.D0,    0.D0,
     :      -114.D0,      0.D0,     0.D0,       57.D0,    0.D0,    0.D0,
     :       113.D0,      0.D0,    -1.D0,      -49.D0,    0.D0,    0.D0,
     :      -102.D0,      0.D0,     0.D0,       44.D0,    0.D0,    0.D0,
     :       -94.D0,      0.D0,     0.D0,       51.D0,    0.D0,    0.D0,
     :      -100.D0,      0.D0,    -1.D0,       56.D0,    0.D0,    0.D0,
     :        87.D0,      0.D0,     0.D0,      -47.D0,    0.D0,    0.D0,
     :       161.D0,      0.D0,     0.D0,       -1.D0,    0.D0,    0.D0,
     :        96.D0,      0.D0,     0.D0,      -50.D0,    0.D0,    0.D0,
     :       151.D0,      0.D0,    -1.D0,       -5.D0,    0.D0,    0.D0/
      DATA ( ( CLS(I,J), I=1,6 ), J=201,220 ) /
     :      -104.D0,      0.D0,     0.D0,       44.D0,    0.D0,    0.D0,
     :      -110.D0,      0.D0,     0.D0,       48.D0,    0.D0,    0.D0,
     :      -100.D0,      0.D0,     1.D0,       50.D0,    0.D0,    0.D0,
     :        92.D0,      0.D0,    -5.D0,       12.D0,    0.D0,   -2.D0,
     :        82.D0,      0.D0,     0.D0,      -45.D0,    0.D0,    0.D0,
     :        82.D0,      0.D0,     0.D0,      -45.D0,    0.D0,    0.D0,
     :       -78.D0,      0.D0,     0.D0,       41.D0,    0.D0,    0.D0,
     :       -77.D0,      0.D0,     0.D0,       43.D0,    0.D0,    0.D0,
     :         2.D0,      0.D0,     0.D0,       54.D0,    0.D0,    0.D0,
     :        94.D0,      0.D0,     0.D0,      -40.D0,    0.D0,    0.D0,
     :       -93.D0,      0.D0,     0.D0,       40.D0,    0.D0,    0.D0,
     :       -83.D0,      0.D0,    10.D0,       40.D0,    0.D0,   -2.D0,
     :        83.D0,      0.D0,     0.D0,      -36.D0,    0.D0,    0.D0,
     :       -91.D0,      0.D0,     0.D0,       39.D0,    0.D0,    0.D0,
     :       128.D0,      0.D0,     0.D0,       -1.D0,    0.D0,    0.D0,
     :       -79.D0,      0.D0,     0.D0,       34.D0,    0.D0,    0.D0,
     :       -83.D0,      0.D0,     0.D0,       47.D0,    0.D0,    0.D0,
     :        84.D0,      0.D0,     0.D0,      -44.D0,    0.D0,    0.D0,
     :        83.D0,      0.D0,     0.D0,      -43.D0,    0.D0,    0.D0,
     :        91.D0,      0.D0,     0.D0,      -39.D0,    0.D0,    0.D0/
      DATA ( ( CLS(I,J), I=1,6 ), J=221,240 ) /
     :       -77.D0,      0.D0,     0.D0,       39.D0,    0.D0,    0.D0,
     :        84.D0,      0.D0,     0.D0,      -43.D0,    0.D0,    0.D0,
     :       -92.D0,      0.D0,     1.D0,       39.D0,    0.D0,    0.D0,
     :       -92.D0,      0.D0,     1.D0,       39.D0,    0.D0,    0.D0,
     :       -94.D0,      0.D0,     0.D0,        0.D0,    0.D0,    0.D0,
     :        68.D0,      0.D0,     0.D0,      -36.D0,    0.D0,    0.D0,
     :       -61.D0,      0.D0,     0.D0,       32.D0,    0.D0,    0.D0,
     :        71.D0,      0.D0,     0.D0,      -31.D0,    0.D0,    0.D0,
     :        62.D0,      0.D0,     0.D0,      -34.D0,    0.D0,    0.D0,
     :       -63.D0,      0.D0,     0.D0,       33.D0,    0.D0,    0.D0,
     :       -73.D0,      0.D0,     0.D0,       32.D0,    0.D0,    0.D0,
     :       115.D0,      0.D0,     0.D0,       -2.D0,    0.D0,    0.D0,
     :      -103.D0,      0.D0,     0.D0,        2.D0,    0.D0,    0.D0,
     :        63.D0,      0.D0,     0.D0,      -28.D0,    0.D0,    0.D0,
     :        74.D0,      0.D0,     0.D0,      -32.D0,    0.D0,    0.D0,
     :      -103.D0,      0.D0,    -3.D0,        3.D0,    0.D0,   -1.D0,
     :       -69.D0,      0.D0,     0.D0,       30.D0,    0.D0,    0.D0,
     :        57.D0,      0.D0,     0.D0,      -29.D0,    0.D0,    0.D0,
     :        94.D0,      0.D0,     0.D0,       -4.D0,    0.D0,    0.D0,
     :        64.D0,      0.D0,     0.D0,      -33.D0,    0.D0,    0.D0/
      DATA ( ( CLS(I,J), I=1,6 ), J=241,260 ) /
     :       -63.D0,      0.D0,     0.D0,       26.D0,    0.D0,    0.D0,
     :       -38.D0,      0.D0,     0.D0,       20.D0,    0.D0,    0.D0,
     :       -43.D0,      0.D0,     0.D0,       24.D0,    0.D0,    0.D0,
     :       -45.D0,      0.D0,     0.D0,       23.D0,    0.D0,    0.D0,
     :        47.D0,      0.D0,     0.D0,      -24.D0,    0.D0,    0.D0,
     :       -48.D0,      0.D0,     0.D0,       25.D0,    0.D0,    0.D0,
     :        45.D0,      0.D0,     0.D0,      -26.D0,    0.D0,    0.D0,
     :        56.D0,      0.D0,     0.D0,      -25.D0,    0.D0,    0.D0,
     :        88.D0,      0.D0,     0.D0,        2.D0,    0.D0,    0.D0,
     :       -75.D0,      0.D0,     0.D0,        0.D0,    0.D0,    0.D0,
     :        85.D0,      0.D0,     0.D0,        0.D0,    0.D0,    0.D0,
     :        49.D0,      0.D0,     0.D0,      -26.D0,    0.D0,    0.D0,
     :       -74.D0,      0.D0,    -3.D0,       -1.D0,    0.D0,   -1.D0,
     :       -39.D0,      0.D0,     0.D0,       21.D0,    0.D0,    0.D0,
     :        45.D0,      0.D0,     0.D0,      -20.D0,    0.D0,    0.D0,
     :        51.D0,      0.D0,     0.D0,      -22.D0,    0.D0,    0.D0,
     :       -40.D0,      0.D0,     0.D0,       21.D0,    0.D0,    0.D0,
     :        41.D0,      0.D0,     0.D0,      -21.D0,    0.D0,    0.D0,
     :       -42.D0,      0.D0,     0.D0,       24.D0,    0.D0,    0.D0,
     :       -51.D0,      0.D0,     0.D0,       22.D0,    0.D0,    0.D0/
      DATA ( ( CLS(I,J), I=1,6 ), J=261,280 ) /
     :       -42.D0,      0.D0,     0.D0,       22.D0,    0.D0,    0.D0,
     :        39.D0,      0.D0,     0.D0,      -21.D0,    0.D0,    0.D0,
     :        46.D0,      0.D0,     0.D0,      -18.D0,    0.D0,    0.D0,
     :       -53.D0,      0.D0,     0.D0,       22.D0,    0.D0,    0.D0,
     :        82.D0,      0.D0,     0.D0,       -4.D0,    0.D0,    0.D0,
     :        81.D0,      0.D0,    -1.D0,       -4.D0,    0.D0,    0.D0,
     :        47.D0,      0.D0,     0.D0,      -19.D0,    0.D0,    0.D0,
     :        53.D0,      0.D0,     0.D0,      -23.D0,    0.D0,    0.D0,
     :       -45.D0,      0.D0,     0.D0,       22.D0,    0.D0,    0.D0,
     :       -44.D0,      0.D0,     0.D0,       -2.D0,    0.D0,    0.D0,
     :       -33.D0,      0.D0,     0.D0,       16.D0,    0.D0,    0.D0,
     :       -61.D0,      0.D0,     0.D0,        1.D0,    0.D0,    0.D0,
     :       -38.D0,      0.D0,     0.D0,       19.D0,    0.D0,    0.D0,
     :       -33.D0,      0.D0,     0.D0,       21.D0,    0.D0,    0.D0,
     :       -60.D0,      0.D0,     0.D0,        0.D0,    0.D0,    0.D0,
     :        48.D0,      0.D0,     0.D0,      -10.D0,    0.D0,    0.D0,
     :        38.D0,      0.D0,     0.D0,      -20.D0,    0.D0,    0.D0,
     :        31.D0,      0.D0,     0.D0,      -13.D0,    0.D0,    0.D0,
     :       -32.D0,      0.D0,     0.D0,       15.D0,    0.D0,    0.D0,
     :        45.D0,      0.D0,     0.D0,       -8.D0,    0.D0,    0.D0/
      DATA ( ( CLS(I,J), I=1,6 ), J=281,300 ) /
     :       -44.D0,      0.D0,     0.D0,       19.D0,    0.D0,    0.D0,
     :       -51.D0,      0.D0,     0.D0,        0.D0,    0.D0,    0.D0,
     :       -36.D0,      0.D0,     0.D0,       20.D0,    0.D0,    0.D0,
     :        44.D0,      0.D0,     0.D0,      -19.D0,    0.D0,    0.D0,
     :       -60.D0,      0.D0,     0.D0,        2.D0,    0.D0,    0.D0,
     :        35.D0,      0.D0,     0.D0,      -18.D0,    0.D0,    0.D0,
     :        47.D0,      0.D0,     0.D0,       -1.D0,    0.D0,    0.D0,
     :        36.D0,      0.D0,     0.D0,      -15.D0,    0.D0,    0.D0,
     :       -36.D0,      0.D0,     0.D0,       20.D0,    0.D0,    0.D0,
     :       -35.D0,      0.D0,     0.D0,       19.D0,    0.D0,    0.D0,
     :       -37.D0,      0.D0,     0.D0,       19.D0,    0.D0,    0.D0,
     :        32.D0,      0.D0,     0.D0,      -16.D0,    0.D0,    0.D0,
     :        35.D0,      0.D0,     0.D0,      -14.D0,    0.D0,    0.D0,
     :        32.D0,      0.D0,     0.D0,      -13.D0,    0.D0,    0.D0,
     :        65.D0,      0.D0,     0.D0,       -2.D0,    0.D0,    0.D0,
     :        47.D0,      0.D0,     0.D0,       -1.D0,    0.D0,    0.D0,
     :        32.D0,      0.D0,     0.D0,      -16.D0,    0.D0,    0.D0,
     :        37.D0,      0.D0,     0.D0,      -16.D0,    0.D0,    0.D0,
     :       -30.D0,      0.D0,     0.D0,       15.D0,    0.D0,    0.D0,
     :       -32.D0,      0.D0,     0.D0,       16.D0,    0.D0,    0.D0/
      DATA ( ( CLS(I,J), I=1,6 ), J=301,320 ) /
     :       -31.D0,      0.D0,     0.D0,       13.D0,    0.D0,    0.D0,
     :        37.D0,      0.D0,     0.D0,      -16.D0,    0.D0,    0.D0,
     :        31.D0,      0.D0,     0.D0,      -13.D0,    0.D0,    0.D0,
     :        49.D0,      0.D0,     0.D0,       -2.D0,    0.D0,    0.D0,
     :        32.D0,      0.D0,     0.D0,      -13.D0,    0.D0,    0.D0,
     :       -43.D0,      0.D0,     0.D0,       18.D0,    0.D0,    0.D0,
     :       -32.D0,      0.D0,     0.D0,       14.D0,    0.D0,    0.D0,
     :        30.D0,      0.D0,     0.D0,        0.D0,    0.D0,    0.D0,
     :       -34.D0,      0.D0,     0.D0,       15.D0,    0.D0,    0.D0,
     :       -36.D0,      0.D0,     0.D0,        0.D0,    0.D0,    0.D0,
     :       -38.D0,      0.D0,     0.D0,        0.D0,    0.D0,    0.D0,
     :       -31.D0,      0.D0,     0.D0,        0.D0,    0.D0,    0.D0,
     :       -34.D0,      0.D0,     0.D0,        0.D0,    0.D0,    0.D0,
     :       -35.D0,      0.D0,     0.D0,        0.D0,    0.D0,    0.D0,
     :        30.D0,      0.D0,     0.D0,       -2.D0,    0.D0,    0.D0,
     :         0.D0,      0.D0, -1988.D0,        0.D0,    0.D0,-1679.D0,
     :         0.D0,      0.D0,   -63.D0,        0.D0,    0.D0,  -27.D0,
     :         0.D0,      0.D0,   364.D0,        0.D0,    0.D0,  176.D0,
     :         0.D0,      0.D0, -1044.D0,        0.D0,    0.D0, -891.D0,
     :         0.D0,      0.D0,   330.D0,        0.D0,    0.D0,    0.D0/
      DATA ( ( CLS(I,J), I=1,6 ), J=321,323 ) /
     :         0.D0,      0.D0,    30.D0,        0.D0,    0.D0,   14.D0,
     :         0.D0,      0.D0,  -162.D0,        0.D0,    0.D0, -138.D0,
     :         0.D0,      0.D0,    75.D0,        0.D0,    0.D0,    0.D0/

*
*  Planetary argument multipliers:
*              L   L'  F   D   Om  Me  Ve  E  Ma  Ju  Sa  Ur  Ne  pre
*
      DATA ( ( NAPL(I,J), I=1,14), J=  1, 20 ) /
     :         0,  0,  0,  0,  0,  0,  0,  8,-16,  4,  5,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -8, 16, -4, -5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  8,-16,  4,  5,  0,  0,  2,
     :         0,  0,  1, -1,  1,  0,  0,  3, -8,  3,  0,  0,  0,  0,
     :        -1,  0,  0,  0,  0,  0, 10, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  4, -8,  3,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -5,  8, -3,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -5,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  2, -5,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  2, -5,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0, -2,  5,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0, -2,  5,  0,  0,  1,
     :         1,  0,  0, -2,  0,  0, 19,-21,  3,  0,  0,  0,  0,  0,
     :         1,  0,  0, -1,  1,  0,  0, -1,  0,  2,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  1,  0,  0,  2,  0, -2,  0,  0,  0,  0,
     :        -1,  0,  0,  0,  0,  0, 18,-16,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -8, 13,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -8, 13,  0,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0, -8, 12,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  8,-13,  0,  0,  0,  0,  0,  0/
      DATA ( ( NAPL(I,J), I=1,14), J= 21, 40 ) /
     :         0,  0, -1,  1,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0, -1,  2,  0,  0,  0,  0,  0,
     :        -1,  0,  0,  1,  0,  0,  3, -4,  0,  0,  0,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  0,  1,  0,  0,  1,  0,  0,  0,
     :         0,  0, -2,  2,  0,  0,  5, -6,  0,  0,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  6, -8,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  8,-15,  0,  0,  0,  0,  0,
     :         2,  0,  0, -2,  1,  0,  0, -2,  0,  3,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  2,  0, -3,  0,  0,  0,  0,
     :        -1,  0,  0,  1,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  0,  1,  0,  1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  0,  0,  1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  0, -1,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  1,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  1,
     :         0,  0,  0,  0,  1,  0,  8,-13,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  0,  0,  1,  0, 18,-16,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  0,  0, -2,  5,  0,  0,  0/
      DATA ( ( NAPL(I,J), I=1,14), J= 41, 60 ) /
     :         0,  0,  0,  0,  1,  0,  0, -4,  8, -3,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  4, -8,  3,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  0,  0,  2, -5,  0,  0,  0,
     :        -2,  0,  0,  2,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,
     :         1,  0,  0,  0,  1,  0,-18, 16,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0, -8, 13,  0,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0,  0, -2,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -2,  2,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -1,  2,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  0,  2,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  2,
     :         0,  0,  2, -2,  1,  0, -5,  6,  0,  0,  0,  0,  0,  0,
     :         0,  0, -1,  1,  0,  0,  5, -7,  0,  0,  0,  0,  0,  0,
     :        -2,  0,  0,  2,  1,  0,  0,  2,  0, -3,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0, -1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0, -1,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  1/
      DATA ( ( NAPL(I,J), I=1,14), J= 61, 80 ) /
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  0,  2,
     :        -2,  0,  0,  2,  0,  0,  3, -3,  0,  0,  0,  0,  0,  0,
     :         2,  0,  0, -2,  1,  0,  0, -2,  0,  2,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  0,  1, -2,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0,  0,  0,  0,
     :         0,  0,  1, -1,  1,  0, -3,  4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -3,  5,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0, -3,  5,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -3,  5,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2, -4,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0, -2,  4,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -5,  8,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -5,  8,  0,  0,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0, -5,  7,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -5,  8,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  5, -8,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  1,
     :         0,  0,  1, -1,  1,  0,  0, -1,  0,  2,  0,  0,  0,  0/
      DATA ( ( NAPL(I,J), I=1,14), J= 81,100 ) /
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  2,
     :         0,  0,  2, -2,  1,  0, -3,  3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  3, -6,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  1,  0,  2, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  2, -2,  1,  0,  0, -2,  0,  2,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -2,  3,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  2, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  0,  3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  3, -5,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  1, -2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  2, -3,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -4,  7,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -4,  6,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -4,  6,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  4, -6,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -1,  1,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  1, -1,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0, -1,  0,  4,  0,  0,  0,  2/
      DATA ( ( NAPL(I,J), I=1,14), J=101,120 ) /
     :         0,  0,  0,  0,  0,  0,  0, -1,  0,  3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0, -3,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -2,  4,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0, -2,  4,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -6,  9,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0, -2,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  3, -4,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  3, -4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0, -1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  0, -1,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  2, -5,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  1,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -5,  7,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -5,  7,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  1,  0,  2,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -2,  2,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  2, -2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -1,  3,  0,  0,  0,  0,  0,  1/
      DATA ( ( NAPL(I,J), I=1,14), J=121,140 ) /
     :         0,  0,  0,  0,  0,  0, -1,  3,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0, -3,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0, -2,  5,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -6,  8,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -3,  3,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  3, -3,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0, -1,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0, -1,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  3, -2,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -8, 15,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  6, -8,  3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  1,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -6, 16, -4, -5,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0, -2,  8, -3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -8, 11,  0,  0,  0,  0,  0,  2/
      DATA ( ( NAPL(I,J), I=1,14), J=141,160 ) /
     :         0,  0,  0,  0,  0,  0,  0,  1,  2,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  2,  0,  1,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -3,  7,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  0,  4,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  2, -1,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -7,  9,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  4, -4,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  1,  1,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  3,  0, -2,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  5, -4,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  3, -2,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  3,  0, -1,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  4, -2,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -8, 10,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  5, -5,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0, -9, 11,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  4,  0, -3,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  6, -6,  0,  0,  0,  0,  0,  0,
     :         0,  0,  0,  0,  0,  0,  0,  4,  0, -2,  0,  0,  0,  2/
      DATA ( ( NAPL(I,J), I=1,14), J=161,165 ) /
     :         0,  0,  0,  0,  0,  0,  3, -1,  0,  0,  0,  0,  0,  2,
     :         0,  0,  0,  0,  0,  0,  0,  4,  0, -1,  0,  0,  0,  2,
     :        -1,  0,  0,  2,  0,  0,  0,  2,  0, -2,  0,  0,  0,  0,
     :         1,  0,  2,  0,  2,  0,  0,  1,  0,  0,  0,  0,  0,  0,
     :        -1,  0,  2,  0,  2,  0,  0, -4,  8, -3,  0,  0,  0,  0/

*
*  Planetary nutation coefficients, unit 1e-7 arcsec:
*  longitude (sin, cos), obliquity (sin, cos)
*
      DATA ( ( CPL(I,J), I=1,4 ), J=  1, 20 ) /
     :    1440.D0,       0.D0,       0.D0,       0.D0,
     :      56.D0,    -117.D0,     -42.D0,     -40.D0,
     :     125.D0,     -43.D0,       0.D0,     -54.D0,
     :    -114.D0,       0.D0,       0.D0,      61.D0,
     :    -219.D0,      89.D0,       0.D0,       0.D0,
     :    -462.D0,    1604.D0,       0.D0,       0.D0,
     :      99.D0,       0.D0,       0.D0,     -53.D0,
     :      14.D0,    -218.D0,     117.D0,       8.D0,
     :      31.D0,    -481.D0,    -257.D0,     -17.D0,
     :    -491.D0,     128.D0,       0.D0,       0.D0,
     :   -3084.D0,    5123.D0,    2735.D0,    1647.D0,
     :   -1444.D0,    2409.D0,   -1286.D0,    -771.D0,
     :     103.D0,     -60.D0,       0.D0,       0.D0,
     :     -26.D0,     -29.D0,     -16.D0,      14.D0,
     :     284.D0,       0.D0,       0.D0,    -151.D0,
     :     226.D0,     101.D0,       0.D0,       0.D0,
     :     -41.D0,     175.D0,      76.D0,      17.D0,
     :     425.D0,     212.D0,    -133.D0,     269.D0,
     :    1200.D0,     598.D0,     319.D0,    -641.D0,
     :     235.D0,     334.D0,       0.D0,       0.D0/
      DATA ( ( CPL(I,J), I=1,4 ), J= 21, 40 ) /
     :     266.D0,     -78.D0,       0.D0,       0.D0,
     :    -460.D0,    -435.D0,    -232.D0,     246.D0,
     :       0.D0,     131.D0,       0.D0,       0.D0,
     :     -42.D0,      20.D0,       0.D0,       0.D0,
     :     -10.D0,     233.D0,       0.D0,       0.D0,
     :      78.D0,     -18.D0,       0.D0,       0.D0,
     :      45.D0,     -22.D0,       0.D0,       0.D0,
     :      89.D0,     -16.D0,      -9.D0,     -48.D0,
     :    -349.D0,     -62.D0,       0.D0,       0.D0,
     :     -53.D0,       0.D0,       0.D0,       0.D0,
     :     -21.D0,     -78.D0,       0.D0,       0.D0,
     :      20.D0,     -70.D0,     -37.D0,     -11.D0,
     :      32.D0,      15.D0,      -8.D0,      17.D0,
     :     174.D0,      84.D0,      45.D0,     -93.D0,
     :      11.D0,      56.D0,       0.D0,       0.D0,
     :     -66.D0,     -12.D0,      -6.D0,      35.D0,
     :      47.D0,       8.D0,       4.D0,     -25.D0,
     :      46.D0,      66.D0,      35.D0,     -25.D0,
     :     -68.D0,     -34.D0,     -18.D0,      36.D0,
     :      76.D0,      17.D0,       9.D0,     -41.D0/
      DATA ( ( CPL(I,J), I=1,4 ), J= 41, 60 ) /
     :      84.D0,     298.D0,     159.D0,     -45.D0,
     :     -82.D0,     292.D0,     156.D0,      44.D0,
     :     -73.D0,      17.D0,       9.D0,      39.D0,
     :    -439.D0,       0.D0,       0.D0,       0.D0,
     :      57.D0,     -28.D0,     -15.D0,     -30.D0,
     :     -40.D0,      57.D0,      30.D0,      21.D0,
     :     273.D0,      80.D0,      43.D0,    -146.D0,
     :    -449.D0,     430.D0,       0.D0,       0.D0,
     :      -8.D0,     -47.D0,     -25.D0,       4.D0,
     :       6.D0,      47.D0,      25.D0,      -3.D0,
     :     -48.D0,    -110.D0,     -59.D0,      26.D0,
     :      51.D0,     114.D0,      61.D0,     -27.D0,
     :    -133.D0,       0.D0,       0.D0,      57.D0,
     :     -18.D0,    -436.D0,    -233.D0,       9.D0,
     :      35.D0,      -7.D0,       0.D0,       0.D0,
     :     -53.D0,      -9.D0,      -5.D0,      28.D0,
     :     -50.D0,     194.D0,     103.D0,      27.D0,
     :     -13.D0,      52.D0,      28.D0,       7.D0,
     :     -91.D0,     248.D0,       0.D0,       0.D0,
     :       6.D0,      49.D0,      26.D0,      -3.D0/
      DATA ( ( CPL(I,J), I=1,4 ), J= 61, 80 ) /
     :      -6.D0,     -47.D0,     -25.D0,       3.D0,
     :      52.D0,      23.D0,      10.D0,     -23.D0,
     :    -138.D0,       0.D0,       0.D0,       0.D0,
     :      54.D0,       0.D0,       0.D0,     -29.D0,
     :     -37.D0,      35.D0,      19.D0,      20.D0,
     :    -145.D0,      47.D0,       0.D0,       0.D0,
     :     -10.D0,      40.D0,      21.D0,       5.D0,
     :      11.D0,     -49.D0,     -26.D0,      -7.D0,
     :   -2150.D0,       0.D0,       0.D0,     932.D0,
     :      85.D0,       0.D0,       0.D0,     -37.D0,
     :     -86.D0,     153.D0,       0.D0,       0.D0,
     :     -51.D0,       0.D0,       0.D0,      22.D0,
     :     -11.D0,    -268.D0,    -116.D0,       5.D0,
     :      31.D0,       6.D0,       3.D0,     -17.D0,
     :     140.D0,      27.D0,      14.D0,     -75.D0,
     :      57.D0,      11.D0,       6.D0,     -30.D0,
     :     -14.D0,     -39.D0,       0.D0,       0.D0,
     :     -25.D0,      22.D0,       0.D0,       0.D0,
     :      42.D0,     223.D0,     119.D0,     -22.D0,
     :     -27.D0,    -143.D0,     -77.D0,      14.D0/
      DATA ( ( CPL(I,J), I=1,4 ), J= 81,100 ) /
     :       9.D0,      49.D0,      26.D0,      -5.D0,
     :   -1166.D0,       0.D0,       0.D0,     505.D0,
     :     117.D0,       0.D0,       0.D0,     -63.D0,
     :       0.D0,      31.D0,       0.D0,       0.D0,
     :       0.D0,     -32.D0,     -17.D0,       0.D0,
     :      50.D0,       0.D0,       0.D0,     -27.D0,
     :      30.D0,      -3.D0,      -2.D0,     -16.D0,
     :       8.D0,     614.D0,       0.D0,       0.D0,
     :    -127.D0,      21.D0,       9.D0,      55.D0,
     :     -20.D0,      34.D0,       0.D0,       0.D0,
     :      22.D0,     -87.D0,       0.D0,       0.D0,
     :     -68.D0,      39.D0,       0.D0,       0.D0,
     :       3.D0,      66.D0,      29.D0,      -1.D0,
     :     490.D0,       0.D0,       0.D0,    -213.D0,
     :     -22.D0,      93.D0,      49.D0,      12.D0,
     :     -46.D0,      14.D0,       0.D0,       0.D0,
     :      25.D0,     106.D0,      57.D0,     -13.D0,
     :    1485.D0,       0.D0,       0.D0,       0.D0,
     :      -7.D0,     -32.D0,     -17.D0,       4.D0,
     :      30.D0,      -6.D0,      -2.D0,     -13.D0/
      DATA ( ( CPL(I,J), I=1,4 ), J=101,120 ) /
     :     118.D0,       0.D0,       0.D0,     -52.D0,
     :     -28.D0,      36.D0,       0.D0,       0.D0,
     :      14.D0,     -59.D0,     -31.D0,      -8.D0,
     :    -458.D0,       0.D0,       0.D0,     198.D0,
     :       0.D0,     -45.D0,     -20.D0,       0.D0,
     :    -166.D0,     269.D0,       0.D0,       0.D0,
     :     -78.D0,      45.D0,       0.D0,       0.D0,
     :      -5.D0,     328.D0,       0.D0,       0.D0,
     :   -1223.D0,     -26.D0,       0.D0,       0.D0,
     :    -368.D0,       0.D0,       0.D0,       0.D0,
     :     -75.D0,       0.D0,       0.D0,       0.D0,
     :     -13.D0,     -30.D0,       0.D0,       0.D0,
     :     -74.D0,       0.D0,       0.D0,      32.D0,
     :    -262.D0,       0.D0,       0.D0,     114.D0,
     :     202.D0,       0.D0,       0.D0,     -87.D0,
     :      -8.D0,      35.D0,      19.D0,       5.D0,
     :     -35.D0,     -48.D0,     -21.D0,      15.D0,
     :      12.D0,      55.D0,      29.D0,      -6.D0,
     :    -598.D0,       0.D0,       0.D0,       0.D0,
     :       8.D0,     -31.D0,     -16.D0,      -4.D0/
      DATA ( ( CPL(I,J), I=1,4 ), J=121,140 ) /
     :     113.D0,       0.D0,       0.D0,     -49.D0,
     :      83.D0,      15.D0,       0.D0,       0.D0,
     :       0.D0,    -114.D0,     -49.D0,       0.D0,
     :     117.D0,       0.D0,       0.D0,     -51.D0,
     :     393.D0,       3.D0,       0.D0,       0.D0,
     :      18.D0,     -29.D0,     -13.D0,      -8.D0,
     :       8.D0,      34.D0,      18.D0,      -4.D0,
     :      89.D0,       0.D0,       0.D0,       0.D0,
     :      54.D0,     -15.D0,      -7.D0,     -24.D0,
     :       0.D0,      35.D0,       0.D0,       0.D0,
     :    -154.D0,     -30.D0,     -13.D0,      67.D0,
     :      80.D0,     -71.D0,     -31.D0,     -35.D0,
     :      61.D0,     -96.D0,     -42.D0,     -27.D0,
     :     123.D0,    -415.D0,    -180.D0,     -53.D0,
     :       0.D0,       0.D0,       0.D0,     -35.D0,
     :       7.D0,     -32.D0,     -17.D0,      -4.D0,
     :     -89.D0,       0.D0,       0.D0,      38.D0,
     :       0.D0,     -86.D0,     -19.D0,      -6.D0,
     :    -123.D0,    -416.D0,    -180.D0,      53.D0,
     :     -62.D0,     -97.D0,     -42.D0,      27.D0/
      DATA ( ( CPL(I,J), I=1,4 ), J=141,160 ) /
     :     -85.D0,     -70.D0,     -31.D0,      37.D0,
     :     163.D0,     -12.D0,      -5.D0,     -72.D0,
     :     -63.D0,     -16.D0,      -7.D0,      28.D0,
     :     -21.D0,     -32.D0,     -14.D0,       9.D0,
     :       5.D0,    -173.D0,     -75.D0,      -2.D0,
     :      74.D0,       0.D0,       0.D0,     -32.D0,
     :      83.D0,       0.D0,       0.D0,       0.D0,
     :    -339.D0,       0.D0,       0.D0,     147.D0,
     :      67.D0,     -91.D0,     -39.D0,     -29.D0,
     :      30.D0,     -18.D0,      -8.D0,     -13.D0,
     :       0.D0,    -114.D0,     -50.D0,       0.D0,
     :     517.D0,      16.D0,       7.D0,    -224.D0,
     :     143.D0,      -3.D0,      -1.D0,     -62.D0,
     :      50.D0,       0.D0,       0.D0,     -22.D0,
     :      59.D0,       0.D0,       0.D0,       0.D0,
     :     370.D0,      -8.D0,       0.D0,    -160.D0,
     :      34.D0,       0.D0,       0.D0,     -15.D0,
     :     -37.D0,      -7.D0,      -3.D0,      16.D0,
     :      40.D0,       0.D0,       0.D0,       0.D0,
     :    -184.D0,      -3.D0,      -1.D0,      80.D0/
      DATA ( ( CPL(I,J), I=1,4 ), J=161,165 ) /
     :      31.D0,      -6.D0,       0.D0,     -13.D0,
     :      -3.D0,     -32.D0,     -14.D0,       1.D0,
     :     -34.D0,       0.D0,       0.D0,       0.D0,
     :     126.D0,     -63.D0,     -27.D0,     -55.D0,
     :    -126.D0,     -63.D0,     -27.D0,      55.D0/

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Interval between fundamental epoch J2000.0 and given date (JC).
      T = ( ( DATE1-DJ0 ) + DATE2 ) / DJC

*  -------------------
*  LUNI-SOLAR NUTATION
*  -------------------

*
*  Fundamental (Delaunay) arguments from Simon et al. (1994)
*
      CALL FUNARG ( T,   EL, ELP, F, D, OM )

*  Initialize the nutation values.
      DP = 0.D0
      DE = 0.D0

*  Summation of luni-solar nutation series (in reverse order).
      DO 100 I = NLS, 1, -1

*     Argument and functions.
         ARG = MOD ( DBLE ( NALS(1,I) ) * EL  +
     :               DBLE ( NALS(2,I) ) * ELP +
     :               DBLE ( NALS(3,I) ) * F   +
     :               DBLE ( NALS(4,I) ) * D   +
     :               DBLE ( NALS(5,I) ) * OM, D2PI )
         SARG = SIN(ARG)
         CARG = COS(ARG)

*     Term.
         DP = DP + ( CLS(1,I) + CLS(2,I) * T ) * SARG
     :           +   CLS(3,I)                  * CARG
         DE = DE + ( CLS(4,I) + CLS(5,I) * T ) * CARG
     :           +   CLS(6,I)                  * SARG

 100  CONTINUE

*  Convert from 0.1 microarcsec units to radians.
      DPSILS = DP * U2R
      DEPSLS = DE * U2R

*  ------------------
*  PLANETARY NUTATION
*  ------------------

*  Planetary longitudes, Mercury through Neptune, wrt mean dynamical
*  ecliptic and equinox of J2000, with high order terms omitted
*  (Simon et al. 1994, 5.8.1-5.8.8).
      ALME = MOD ( 4.402608842461D0 + 2608.790314157421D0 * T, D2PI )
      ALVE = MOD ( 3.176146696956D0 + 1021.328554621099D0 * T, D2PI )
      ALEA = MOD ( 1.753470459496D0 +  628.307584999142D0 * T, D2PI )
      ALMA = MOD ( 6.203476112911D0 +  334.061242669982D0 * T, D2PI )
      ALJU = MOD ( 0.599547105074D0 +   52.969096264064D0 * T, D2PI )
      ALSA = MOD ( 0.874016284019D0 +   21.329910496032D0 * T, D2PI )
      ALUR = MOD ( 5.481293871537D0 +    7.478159856729D0 * T, D2PI )
      ALNE = MOD ( 5.311886286677D0 +    3.813303563778D0 * T, D2PI )

*  General precession in longitude (Simon et al. 1994), equivalent
*  to 5028.8200 arcsec/cy at J2000.
      APA = ( 0.024380407358D0 + 0.000005391235D0 * T ) * T

*  Initialize the nutation values.
      DP = 0.D0
      DE = 0.D0

*  Summation of planetary nutation series (in reverse order).
      DO 200 I = NPL, 1, -1

*     Argument and functions.
         ARG = MOD ( DBLE ( NAPL( 1,I) ) * EL   +
     :               DBLE ( NAPL( 2,I) ) * ELP  +
     :               DBLE ( NAPL( 3,I) ) * F    +
     :               DBLE ( NAPL( 4,I) ) * D    +
     :               DBLE ( NAPL( 5,I) ) * OM   +
     :               DBLE ( NAPL( 6,I) ) * ALME +
     :               DBLE ( NAPL( 7,I) ) * ALVE +
     :               DBLE ( NAPL( 8,I) ) * ALEA +
     :               DBLE ( NAPL( 9,I) ) * ALMA +
     :               DBLE ( NAPL(10,I) ) * ALJU +
     :               DBLE ( NAPL(11,I) ) * ALSA +
     :               DBLE ( NAPL(12,I) ) * ALUR +
     :               DBLE ( NAPL(13,I) ) * ALNE +
     :               DBLE ( NAPL(14,I) ) * APA, D2PI )
         SARG = SIN(ARG)
         CARG = COS(ARG)

*     Term.
         DP = DP + CPL(1,I) * SARG + CPL(2,I) * CARG
         DE = DE + CPL(3,I) * SARG + CPL(4,I) * CARG

 200  CONTINUE

*  Convert from 0.1 microarcsec units to radians.
      DPSIPL = DP * U2R
      DEPSPL = DE * U2R

*  -----
*  TOTAL
*  -----

*  Add planetary and luni-solar components.
      DPSI = DPSIPL + DPSILS
      DEPS = DEPSPL + DEPSLS

      RETURN

      END



      DOUBLE PRECISION FUNCTION EECT2000 ( DATE1, DATE2 )
*+
*  - - - - - - - - -
*   E E C T 2 0 0 0
*  - - - - - - - - -
*
*  Equation of the equinoxes complementary terms, consistent with
*  IAU 2000 resolutions.
*
*  Annexe to IERS Conventions 2000, Chapter 5
*
*  Capitaine, N., Wallace, P.T., & McCarthy, D.D. (2003). Astron. &
*    Astrophys. 406, pp. 1135-1149, Table 3.
*  IERS Conventions (2010), Chapter 5, p. 60, Table 5.2e.
*    (Table 5.2e presented in the printed publication is a truncated
*    series. The full series, which is used in NOVAS, is available on
*    the IERS Conventions Center website in file tab5.2e.txt.)
*    ftp://tai.bipm.org/iers/conv2010/chapter5/
*
*  Given:
*     DATE1,DATE2   d    TT date (JD = DATE1+DATE2)
*
*  Returned:
*     EECT00        d    complementary terms (radians)
*
*  This revision:  2002 November 13
*                  References updated 2010 November 26
*
*-----------------------------------------------------------------------

      IMPLICIT NONE

      DOUBLE PRECISION DATE1, DATE2

*  2Pi
      DOUBLE PRECISION D2PI
      PARAMETER ( D2PI = 6.283185307179586476925287D0 )

*  Arcseconds to radians
      DOUBLE PRECISION DAS2R
      PARAMETER ( DAS2R = 4.848136811095359935899141D-6 )

*  Reference epoch (J2000), JD
      DOUBLE PRECISION DJ0
      PARAMETER ( DJ0 = 2451545D0 )

*  Days per Julian century
      DOUBLE PRECISION DJC
      PARAMETER ( DJC = 36525D0 )

*  Time since J2000, in Julian centuries
      DOUBLE PRECISION T

*  Miscellaneous
      INTEGER I, J
      DOUBLE PRECISION A, S0, S1
      DOUBLE PRECISION ANMP

*  Fundamental arguments
      DOUBLE PRECISION FA(14)

*  -----------------------------------------
*  The series for the EE complementary terms
*  -----------------------------------------

*  Number of terms in the series
      INTEGER NE0, NE1
      PARAMETER ( NE0=  33, NE1=  1 )

*  Coefficients of l,l',F,D,Om,LMe,LVe,LE,LMa,LJu,LSa,LU,LN,pA
      INTEGER KE0 ( 14, NE0 ),
     :        KE1 ( 14, NE1 )

*  Sine and cosine coefficients
      DOUBLE PRECISION SE0 ( 2, NE0 ),
     :                 SE1 ( 2, NE1 )

*  Argument coefficients for t^0
      DATA ( ( KE0(I,J), I=1,14), J =    1,   10 ) /
     :  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  0,  2, -2,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  0,  2, -2,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  0,  2, -2,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  0,  2,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  0,  2,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  0,  0,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  1,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  1,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0 /
      DATA ( ( KE0(I,J), I=1,14), J =   11,   20 ) /
     :  1,  0,  0,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  1,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  1,  2, -2,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  1,  2, -2,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  0,  4, -4,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  0,  1, -1,  1,  0, -8, 12,  0,  0,  0,  0,  0,  0,
     :  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  0,  2,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  1,  0,  2,  0,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  1,  0,  2,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0 /
      DATA ( ( KE0(I,J), I=1,14), J =   21,   30 ) /
     :  0,  0,  2, -2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  1, -2,  2, -3,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  1, -2,  2, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  0,  0,  0,  0,  0,  8,-13,  0,  0,  0,  0,  0, -1,
     :  0,  0,  0,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  2,  0, -2,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  1,  0,  0, -2,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  1,  2, -2,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  1,  0,  0, -2, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  0,  0,  4, -2,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0 /
      DATA ( ( KE0(I,J), I=1,14), J =   31,  NE0 ) /
     :  0,  0,  2, -2,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  1,  0, -2,  0, -3,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     :  1,  0, -2,  0, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0 /

*  Argument coefficients for t^1
      DATA ( ( KE1(I,J), I=1,14), J =    1,  NE1 ) /
     :  0,  0,  0,  0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0 /

*  Sine and cosine coefficients for t^0
      DATA ( ( SE0(I,J), I=1,2), J =    1,   10 ) /
     :            +2640.96D-6,          -0.39D-6,
     :              +63.52D-6,          -0.02D-6,
     :              +11.75D-6,          +0.01D-6,
     :              +11.21D-6,          +0.01D-6,
     :               -4.55D-6,          +0.00D-6,
     :               +2.02D-6,          +0.00D-6,
     :               +1.98D-6,          +0.00D-6,
     :               -1.72D-6,          +0.00D-6,
     :               -1.41D-6,          -0.01D-6,
     :               -1.26D-6,          -0.01D-6 /
      DATA ( ( SE0(I,J), I=1,2), J =   11,   20 ) /
     :               -0.63D-6,          +0.00D-6,
     :               -0.63D-6,          +0.00D-6,
     :               +0.46D-6,          +0.00D-6,
     :               +0.45D-6,          +0.00D-6,
     :               +0.36D-6,          +0.00D-6,
     :               -0.24D-6,          -0.12D-6,
     :               +0.32D-6,          +0.00D-6,
     :               +0.28D-6,          +0.00D-6,
     :               +0.27D-6,          +0.00D-6,
     :               +0.26D-6,          +0.00D-6 /
      DATA ( ( SE0(I,J), I=1,2), J =   21,   30 ) /
     :               -0.21D-6,          +0.00D-6,
     :               +0.19D-6,          +0.00D-6,
     :               +0.18D-6,          +0.00D-6,
     :               -0.10D-6,          +0.05D-6,
     :               +0.15D-6,          +0.00D-6,
     :               -0.14D-6,          +0.00D-6,
     :               +0.14D-6,          +0.00D-6,
     :               -0.14D-6,          +0.00D-6,
     :               +0.14D-6,          +0.00D-6,
     :               +0.13D-6,          +0.00D-6 /
      DATA ( ( SE0(I,J), I=1,2), J =   31,  NE0 ) /
     :               -0.11D-6,          +0.00D-6,
     :               +0.11D-6,          +0.00D-6,
     :               +0.11D-6,          +0.00D-6 /

*  Sine and cosine coefficients for t^1
      DATA ( ( SE1(I,J), I=1,2), J =    1,  NE1 ) /
     :               -0.87D-6,          +0.00D-6 /

* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

*  Interval between fundamental epoch J2000.0 and current date (JC).
      T = ( ( DATE1-DJ0 ) + DATE2 ) / DJC

*  Fundamental Arguments (from IERS Conventions 2000)

*  Mean Anomaly of the Moon.
      FA(1) = ANMP ( ( 485868.249036D0 +
     :               ( 715923.2178D0 +
     :               (     31.8792D0 +
     :               (      0.051635D0 +
     :               (     -0.00024470D0 )
     :               * T ) * T ) * T ) * T ) * DAS2R
     :               + MOD ( 1325D0*T, 1D0 ) * D2PI )

*  Mean Anomaly of the Sun.
      FA(2) = ANMP ( ( 1287104.793048D0 +
     :               ( 1292581.0481D0 +
     :               (      -0.5532D0 +
     :               (      +0.000136D0 +
     :               (      -0.00001149D0 )
     :               * T ) * T ) * T ) * T ) * DAS2R
     :               + MOD ( 99D0*T, 1D0 ) * D2PI )

*  Mean Longitude of the Moon minus Mean Longitude of the Ascending
*  Node of the Moon.
      FA(3) = ANMP ( (  335779.526232D0 +
     :               (  295262.8478D0 +
     :               (     -12.7512D0 +
     :               (      -0.001037D0 +
     :               (       0.00000417D0 )
     :               * T ) * T ) * T ) * T ) * DAS2R
     :               + MOD ( 1342D0*T, 1D0 ) * D2PI )

*  Mean Elongation of the Moon from the Sun.
      FA(4) = ANMP ( ( 1072260.703692D0 +
     :               ( 1105601.2090D0 +
     :               (      -6.3706D0 +
     :               (       0.006593D0 +
     :               (      -0.00003169D0 )
     :               * T ) * T ) * T ) * T ) * DAS2R
     :               + MOD ( 1236D0*T, 1D0 ) * D2PI )

*  Mean Longitude of the Ascending Node of the Moon.
      FA(5) = ANMP ( (  450160.398036D0 +
     :               ( -482890.5431D0 +
     :               (       7.4722D0 +
     :               (       0.007702D0 +
     :               (      -0.00005939D0 )
     :               * T ) * T ) * T ) * T ) * DAS2R
     :               + MOD ( -5D0*T, 1D0 ) * D2PI )

      FA( 6) = ANMP ( 4.402608842D0 + 2608.7903141574D0 * T )
      FA( 7) = ANMP ( 3.176146697D0 + 1021.3285546211D0 * T )
      FA( 8) = ANMP ( 1.753470314D0 +  628.3075849991D0 * T )
      FA( 9) = ANMP ( 6.203480913D0 +  334.0612426700D0 * T )
      FA(10) = ANMP ( 0.599546497D0 +   52.9690962641D0 * T )
      FA(11) = ANMP ( 0.874016757D0 +   21.3299104960D0 * T )
      FA(12) = ANMP ( 5.481293872D0 +    7.4781598567D0 * T )
      FA(13) = ANMP ( 5.311886287D0 +    3.8133035638D0 * T )
      FA(14) =      ( 0.024381750D0 +    0.00000538691D0 * T ) * T

*  Evaluate the EE complementary terms.
      S0 = 0D0
      S1 = 0D0

      DO I = NE0,1,-1
         A = 0D0
         DO J=1,14
            A = A + DBLE(KE0(J,I))*FA(J)
         END DO
         S0 = S0 + ( SE0(1,I)*SIN(A) + SE0(2,I)*COS(A) )
      END DO
      DO I = NE1,1,-1
         A = 0D0
         DO J=1,14
            A = A + DBLE(KE1(J,I))*FA(J)
         END DO
         S1 = S1 + ( SE1(1,I)*SIN(A) + SE1(2,I)*COS(A) )
      END DO
      EECT2000 = ( S0 + S1 * T ) * DAS2R

*  Finished.

      END



      DOUBLE PRECISION FUNCTION ANMP ( A )

*  Normalize angle into the range -pi <= A < +pi.

      IMPLICIT NONE

      DOUBLE PRECISION A

      DOUBLE PRECISION DPI, D2PI
      PARAMETER ( DPI = 3.141592653589793238462643D0,
     :            D2PI = 6.283185307179586476925287D0 )

      DOUBLE PRECISION W

      W = MOD(A,D2PI)
      IF ( ABS(W) .GE. DPI ) W = W - SIGN(D2PI,A)
      ANMP = W

      END
