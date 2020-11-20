      program fit_parallax_multi

c----------------------------
c     Fits parallax and proper motion to VLBA data

c     Input files:
c        fit_parallax_files.inp     contains names of input data files to fit
c        fit_parallax_control.inp   various controls and parameters

c----------------------------
c     Version 4d: allows flagging of data lines with "!"
c     Version 4c: fixed bug in tied proper motions
c                 earlier versions effectively only used data
c                 from last tied source for proper motion.
c     Version 4b: changes error floor to quadrature combination

c     Version 4: can tie fringe rates for different data sets 
c         (presumably same spectral channel against different QSOs)
c         together.
c     Version 3l fixes "ioef=0" typo
c     Version 3k adds separate sigma_pdf calculations for RA and Dec
c     Version 3j adds separate error floors for RA and Dec data (Jan Forbrich)
c     Version 3 accounts for the eccentricity of the Earth's orbit

c     Fixed bug in y-weighting

      implicit real*8 (a-h,o-z)

      character*48  data_file, out_file, data_files(100)
      real*8     dates(1000), x(1000), xerr(1000), y(1000), yerr(1000)
      real*8     times(1000), data(1000), model(1000), resids(1000), 
     +           res_err(1000), partls(1000,401)
      real*8     params(401), new_params(401), param_sigmas(401),
     +           scaled_sigmas(401)
      real*8     model_slope_rem
      character*16  parnames(401)
      integer    paramids(401)
      logical    print_cor, print_resids

      integer    nd_by_source(100)

      integer    pos_type(1000), source_number(1000)
      common /model_info/ t_mid,source_number,pos_type,ra_rad,dec_rad

c     Set some maximum control values, switches, and logical unit numbers...
      idebug       = 0                       ! used by least-sq routine
      itermax      = 9
      max_num_data = 1000
      max_num_files= 100                     ! sources or channels

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
c     Get input data file names...
      call get_data_file_names ( lu_control, lu_out, max_num_files,
     +                  data_files, n_files )

c     Get input control parameters...
      call control_parameters ( lu_control, lu_out, n_files,
     +                  ra, dec, t_0, err_min_ra, err_min_dec,
     +                  num_params, num_solved_for,
     +                  params, paramids, parnames )

      t_mid = t_0
c     Convert RA, Dec to radians
      ra_rad  = ra * pi/12.d0                      ! radians
      dec_rad = dec * pi/180.d0

c     =============================================================
c     Read all data files...

      i_d = 0
      do n_f = 1, n_files
c        Input data for a source (or a spectral channel) ...  

         data_file = data_files(n_f)
         call read_data( lu_data, data_file, lu_out,
     +                n_pts, dates, x, xerr, y, yerr )

         num_src = n_f            ! source / chan index
         nd_by_source(n_f) = n_pts

c        Transfer 2-D data to 1-D working arrays...
         do n = 1, n_pts

              i_d = i_d + 1
              data(i_d)    = x(n)
              times(i_d)   = dates(n)
              res_err(i_d) = xerr(n)       
c4b              if ( err_min_ra .gt. xerr(n) ) res_err(i_d)=err_min_ra
              res_err(i_d)=sqrt(err_min_ra**2 + xerr(n)**2)
              pos_type(i_d)= 1               ! indicates East offset
              source_number(i_d) = num_src   ! source number

              i_d = i_d + 1
              data(i_d)    = y(n)
              times(i_d)   = dates(n)
              res_err(i_d) = yerr(n)   
C             V2: Fixed bug following line...(used to be .gt. xerr) 
c4b              if ( err_min_dec .gt. yerr(n) ) res_err(i_d)=err_min_dec
              res_err(i_d)=sqrt(err_min_dec**2 + yerr(n)**2)
              pos_type(i_d)= 2               ! indicates North offset
              source_number(i_d) = num_src   ! source number
                                             
         enddo

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
 
                call calc_sigsq ( num_data, num_solved_for,
     +                  resids, res_err, 
     +                  sigsq_old, sigsq, sigma_pdf)

                iterpass = iterpass + 1

        enddo

        write (lu_print,2500) iterpass, sigma_pdf, 
     +                                  num_deg_freedom
2500    format (/' Iteration #',I3,': sigma_pdf =',
     +                  f10.3,' for',i4,' degrees of freedom')

c       Print out scaled uncertainties...
        write (lu_print,2398)
 2398   format(/'   Uncertainties scaled by "sigma_pdf":')
        do n_p=1, num_params
           scaled_sigmas(n_p) = sigma_pdf * param_sigmas(n_p)
           write (lu_print,2399) parnames(n_p), scaled_sigmas(n_p)
 2399      format(3x,a16,f12.6)
        enddo


C       ------------------------------------------------------------
C                               PRINTOUT
C       Printout residuals (all source/chanenls in 1 column) ?
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
c       Loop through all source/channels and make separate plotting files...
      do n_src = 1, n_files

        i0 = 1 + ( n_src - 1 )*4  ! params index before first motion parameter

c       Make ascii file containing the velocity, data, model, and resids 
c       over the velocity range used.
        do i_pos = 1, 2        ! separate files for RA an DEC models

         if ( i_pos .eq. 1 ) then
            write (out_file,4050) n_src
 4050       format('par_fit_results_ra.dat_',i3.3)
         endif
         if ( i_pos .eq. 2 ) then
            write (out_file,4060) n_src
 4060       format('par_fit_results_dec.dat_',i3.3)
         endif

         open (unit=lu_out, file=out_file, form='formatted')

         write (lu_out,4100) 
 4100    format( '!    Times     Data      Model     Resids  ScaledErr')
         do n = 1, num_data, 2
            nn = (n-1) + i_pos
            if ( n_src .eq. source_number(nn) ) then
               scaled_err = res_err(nn) * sigma_pdf
               write (lu_out,4200) times(nn),data(nn),model(nn),
     +                             resids(nn),scaled_err
 4200          format(1x,f10.3,4f10.4)
            endif
         enddo

         close ( unit=lu_out )

        enddo

c       Make another output file of the model evaluated on a finer
c       grid of velocities for smoother plotting...
        do i_pos = 1, 2        ! separate files for RA an DEC models  

         if ( i_pos .eq. 1 ) then
            write (out_file,4052) n_src
 4052       format('par_fit_model_ra.dat_',i3.3)
         endif
         if ( i_pos .eq. 2 ) then
            write (out_file,4062) n_src
 4062       format('par_fit_model_dec.dat_',i3.3)
         endif

         open (unit=lu_out_2, file=out_file, form='formatted')

         n_mod_pts = 101
         t_span = times(num_data) - times(1)
         dt     = t_span / 100.d0
         do i = 1, n_mod_pts
           data_time = times(1) + (i-1)*dt
           call calc_model ( params, data_time, i_pos, n_src,
     +                       position ) 
           write (lu_out_2,4200) data_time, position
	 enddo

         close ( unit=lu_out_2 )

        enddo

c       -------------------
c       Make even more output files of data and model with slope removed...
        do i_pos = 1, 2         ! separate files for RA an DEC models

         if ( i_pos .eq. 1 ) then
            write (out_file,4054) n_src
 4054       format('par_fit_results_desloped_ra.dat_',i3.3)
         endif
         if ( i_pos .eq. 2 ) then
            write (out_file,4064) n_src
 4064       format('par_fit_results_desloped_dec.dat_',i3.3)
         endif

         open (unit=lu_out_3, file=out_file, form='formatted')

         write (lu_out_3,4101) 
 4101    format( '!   Velocity   Data      Model     Resids',
     +          '  ScaledErr   with baseline removed.')

         do n = 1, num_data, 2
            nn = (n-1) + i_pos
            if ( n_src .eq. source_number(nn) ) then
               del_t = times(nn) - t_mid
c              slope contribution...
c              (NB: slope pivots about t_mid)
               if ( i_pos .eq. 1 ) then
c                 RA data point...
                  x_off   = params( i0 + 1 )
                  x_slope = params( i0 + 2 )
                  slope_pos = x_off + x_slope * del_t
                 else
c                 Dec data point...
                  y_off   = params( i0 + 3 )
                  y_slope = params( i0 + 4 )
                  slope_pos = y_off + y_slope * del_t
               endif
               data_slope_rem = data(nn) - slope_pos
               model_slope_rem= model(nn)- slope_pos
               scaled_err = res_err(nn) * sigma_pdf
               write (lu_out_3,4200) times(nn), data_slope_rem, 
     +                            model_slope_rem, resids(nn),
     +                            scaled_err
            endif
         enddo

         close ( unit=lu_out_3 )

        enddo

c       ---------------------------
        do i_pos = 1, 2           ! separate files for RA an DEC models
         if ( i_pos .eq. 1 ) then
            write (out_file,4056) n_src
 4056       format('par_fit_model_desloped_ra.dat_',i3.3)
         endif
         if ( i_pos .eq. 2 ) then
            write (out_file,4066) n_src
 4066       format('par_fit_model_desloped_dec.dat_',i3.3)
         endif

         open (unit=lu_out_4, file=out_file, form='formatted')
         n_mod_pts = 101
         do i = 1, n_mod_pts
           data_time = times(1) + (i-1)*dt
           del_t     = data_time - t_mid
           call calc_model ( params, data_time, i_pos, n_src,
     +                       position ) 
           if ( i_pos .eq. 1 ) then
c                 RA data point...
                  x_off   = params( i0 + 1 )
                  x_slope = params( i0 + 2 )
                  slope_pos = x_off + x_slope * del_t
              else
c                 Dec data point...
                  y_off   = params( i0 + 3 )
                  y_slope = params( i0 + 4 )
                  slope_pos = y_off + y_slope * del_t
           endif
           position_slope_removed = position - slope_pos
           write (lu_out_4,4200) data_time, position_slope_removed
	 enddo
         close ( unit=lu_out_4 )
        enddo

      enddo      !  end cycling through all sources/channels

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
     +                       n_pts, dates, x, xerr, y, yerr )

      implicit real*8 (a-h,o-z)

      character*48  data_file
      character*1   c_1
      
      real*8   dates(1000), x(1000), xerr(1000), y(1000), yerr(1000)

      call open_ascii_file ( lu_data, data_file, lu_out )

      i    = 0
      ieof = 0

C     Now read in parallax data...
      do while (ieof.eq.0)

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
 2000       format(1x,5f10.3)

         endif

      enddo

      n_pts = i
      print*,' Number of un-flagged data epochs read in: ', n_pts

      close ( unit=lu_data )

      return
      end
c
c===================================================================
        subroutine control_parameters ( lu_control, lu_out, n_files,
     +                  ra, dec, t_0, err_min_ra, err_min_dec,
     +                  num_params, num_solved_for,
     +                  params, paramids, parnames )


        implicit real*8 (a-h,o-z)

        real*8          params(401), x_paramids(401)
        character*16    parnames(401)
        integer         paramids(401)
        character*48    control_file

        control_file = 'fit_parallax_control.inp'   ! hardwired input control file name

        call open_ascii_file ( lu_control, control_file, lu_out )

        read (lu_control,8010) ra                  ! R.A. (decimal hours)
        read (lu_control,8010) dec                 ! Dec  (decimal degrees)
        read (lu_control,8010) t_0                 ! mid-time (years)
        read (lu_control,8010) err_min_ra          ! data error floor for Right Ascension (mas)
	read (lu_control,8010) err_min_dec         ! data error floor for Declination (mas)

        read (lu_control,8010) params(1), x_paramids(1) ! parallax

        do n_src = 1, n_files
           i0 = 1 + ( n_src - 1 )*4  ! starting index for motion parameters
           read (lu_control,8010) params(i0+1), x_paramids(i0+1) ! ra offset 
           read (lu_control,8010) params(i0+2), x_paramids(i0+2) ! ra slope
           read (lu_control,8010) params(i0+3), x_paramids(i0+3) ! dec offset
           read (lu_control,8010) params(i0+4), x_paramids(i0+4) ! dec slope
        enddo

8010    format(2f10.3)

        write (6,8010) ra                       ! RA (hr)
        write (6,8010) dec                      ! Dec (deg)
        write (6,8010) t_0                      ! t-mid (years)
        write (6,8010) err_min_ra               ! data error floor (mas)
	write (6,8010) err_min_dec              ! data error floor (mas)

        write (6,8010) params(1), x_paramids(1) ! parallax

        do n_src = 1, n_files
           i0 = 1 + ( n_src - 1 )*4  ! starting index for motion parameters
           write(6,8010) params(i0+1), x_paramids(i0+1) ! ra offset 
           write(6,8010) params(i0+2), x_paramids(i0+2) ! ra slope
           write(6,8010) params(i0+3), x_paramids(i0+3) ! dec offset
           write(6,8010) params(i0+4), x_paramids(i0+4) ! dec slope
        enddo

C       Set solve-for codes...
        num_params = 1 + n_files*4               ! number of parms
        num_solved_for = 0
        do n = 1, num_params
           paramids(n) = int ( x_paramids(n) )
           if (  paramids(n) .eq. 1 ) then
                 num_solved_for = num_solved_for + 1
           endif
        enddo

c       Tie parameter values together (if ID=2), 
        call uncollapse ( num_params, paramids, params )

C       Enter parameter names...
        parnames(1) = 'Parallax (mas)  '

        do n_src = 1, n_files
           i0 = 1 + ( n_src - 1 )*4  ! starting index for motion parameters
           write (parnames(i0+1),8021) n_src
 8021      format(i3.3,' X off (mas) ')
           write (parnames(i0+2),8022) n_src
 8022      format(i3.3,' X slp (m/y) ')
           write (parnames(i0+3),8023) n_src
 8023      format(i3.3,' Y off (mas) ')
           write (parnames(i0+4),8024) n_src
 8024      format(i3.3,' Y slp (m/y) ')
        enddo
        
        close(unit=lu_control)

        return
        end
c
c===================================================================
        subroutine uncollapse ( num_params, paramids, params )

c       Resets proper motion "slope" parameters if tied together
c       because same maser spot with different QSOs

c       params array contains
c          params(1)     parallax
c          params(2)     x-offset input file 1
c          params(3)     x-slope  
c          params(4)     y-offset
c          params(5)     y-slope
c          params(6)     x-offset input file 2
c          params(7)     x-slope  
c          params(8)     y-offset
c          params(9)     y-slope
c           ...

	implicit real*8 (a-h, o-z)

	real*8		params(401)
        integer         paramids(401)

c       Check if at least 2 input files (9 parameters)
        if ( num_params .ge. 9 ) then

           do j = num_params, 7, -1

c             Only odd numbered parameters (slopes) tied together
              jeven = 2 * (j / 2)
              if ( j .gt. jeven ) then

c                Check if parameter ID >= 2, if so uncollapse
                 jj = j
                 do while ( paramids(jj).ge.2 .and. jj.ge.7 )
                     jj = jj - 4
                 enddo
                 
                 params(j) = params(jj)

             endif

           enddo

        endif

        return
        end
c
c===================================================================
      subroutine collapse ( num_data, num_params, paramids, 
     +                      partls )

c     Rearranges partials if some parameters tied together

      implicit real*8 (a-h,o-z)

      real*8     partls(1000,401)
      integer    paramids(401)

      integer    pos_type(1000), source_number(1000)
      common /model_info/ t_mid,source_number,pos_type,ra_rad,dec_rad

      do i = 1, num_data

        if ( num_params .ge. 9 ) then

           do j = num_params, 7, -1

c             Only odd numbered parameters (slopes) tied together
              jeven = 2 * (j / 2)
              if ( j .gt. jeven ) then

c                Check if parameter ID >= 2, if so collapse
                 jj = j
                 do while ( paramids(jj).ge.2 .and. jj.ge.7 )
                     jj = jj - 4
                 enddo

c                Check that jj.ne.j and only move partial if source
c                group (indicated by parameter number j) matches
c                source number of data point (from i-index)
                 n_src_j = int( (j+2)/4 )
                 n_src_i = source_number(i)

                 if ( jj.ne.j .and. n_src_j.eq.n_src_i ) then
                    partls(i,jj) = partls(i,j)
                    partls(i,j)  = 0.d0
                 endif

             endif

           enddo

        endif

      enddo

      return
      end
c
c===================================================================
        subroutine get_data_file_names ( lu_control, lu_out, 
     +                  max_num_files,
     +                  data_files, n_files )


        implicit real*8 (a-h,o-z)

        character*48 control_file, file_name, blank48, data_files(100)

        blank48 = '                                                '

        control_file = 'fit_parallax_files.inp'    ! input name of file containing data files 

        call open_ascii_file ( lu_control, control_file, lu_out )

        n_files = 0
        ieof    = 0
        do while ( ieof.ge.0 .and. n_files.lt.max_num_files )

         read(lu_control,*,iostat=ieof) file_name

c        Check for end of file...
         if ( ieof.ge.0 .and. file_name.ne.blank48 ) then

            n_files = n_files + 1
            data_files( n_files ) = file_name

         endif

        enddo
        
        print *,' Number of data file names =', n_files
        
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

        integer    pos_type(1000), source_number(1000)
        common /model_info/ t_mid,source_number,pos_type,ra_rad,dec_rad

        do n = 1, num_data

           data_time = times(n)
           i_pos     = pos_type(n)
           i_src     = source_number(n)
           call calc_model ( params, data_time, i_pos, i_src,
     +                       position ) 

	   model(n) = position
	   resids(n) = data(n) - model(n)

	enddo
	
	return
	end
c
	subroutine calc_model ( params, data_time, i_pos, i_src,
     +                          position ) 

c       Calculates a position (either RA or Dec given the data_time
c       and the model parameters of the parallax and 2 sloping lines
c       (one each for the eastward and northward motions). 

	implicit real*8 ( a-h, o-z )

        real*8     params(401)

        integer    pos_type(1000), source_number(1000)
        common /model_info/ t_mid,source_number,pos_type,ra_rad,dec_rad


c       Slope contribution...
c       (NB: slope pivots about t_mid)
        del_t = data_time - t_mid

c       Get proper motion parameters for the correct source/chan
        i0 = 1 + ( i_src - 1 )*4     ! params index before 1'st motion param

        if ( i_pos .eq. 1 ) then
c             RA data... 
              x_off   = params( i0 + 1 )
              x_slope = params( i0 + 2 )
              slope_pos = x_off + x_slope*del_t
           else
c             Dec data...
              y_off   = params( i0 + 3 )
              y_slope = params( i0 + 4 )
              slope_pos = y_off + y_slope*del_t
        endif

c       Parallax offsets...
        call calc_parallax ( params, data_time, i_pos,
     +                       parallax_offset )

c       Sum contributions...
        position = slope_pos + parallax_offset

	return 
	end
c
	subroutine calc_parallax ( params, t_yr, i_pos,
     +                       parallax_offset )

c	Calculates position offset expected from the parallax
c       of a source

        implicit real*8 (a-h, o-z)

        real*8     params(401)

        integer    pos_type(1000), source_number(1000)
        common /model_info/ t_mid,source_number,pos_type,ra_rad,dec_rad

        pi    = 4.d0*atan(1.d0)
        twopi = 2.d0 * pi
        deg_to_rad = pi / 180.d0

c       longitude of sun is 0 at vernal equinox (~0.22 yrs)
        sun_long = twopi * ( t_yr - 0.22d0 )    ! radians
 
        cos_sun = cos( sun_long )
        sin_sun = sin( sun_long )
 
        obliquity = 23.4d0 * deg_to_rad         ! radians
 
        cos_obl = cos( obliquity )
        sin_obl = sin( obliquity )
 
c       Sun's position in cartesian coordinates in AU units
c       for a circular orbit
        X = cos_sun
        Y = sin_sun * cos_obl
        Z = sin_sun * sin_obl

        cos_ra = cos(ra_rad)
        sin_ra = sin(ra_rad)
        cos_dec= cos(dec_rad)
        sin_dec= sin(dec_rad)
 
        parallax = params(1)                    ! mas

c       Correct for eccentricity of the Earth's orbit
        earth_long = twopi * ( t_yr - 0.257d0 )    ! radians
        factor = 1.0d0 + 0.0167d0 * sin( earth_long )

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
	subroutine calc_partials ( num_params, params, paramids,
     +                 num_data, times, 
     +                 partls )

	implicit real*8 ( a-h,o-z )

	real*8     params(401), params_wiggled(401)
	real*8	   partls(1000,401), times(1000)

	integer	   paramids(401)

        integer    pos_type(1000), source_number(1000)
        common /model_info/ t_mid,source_number,pos_type,ra_rad,dec_rad

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
            i_src     = source_number(i_d)
            call calc_model ( params, data_time, i_pos, i_src,
     +                        position )

            do i_p = 1, num_params

c               Calc partials for ID=1 or 2...
		if ( paramids(i_p) .ge. 1 ) then

c		   Change parameters by 1 part in 10^4; but never
c			by less than a value of 10^-6 (Jy or kms)...

                   del_param = abs( params(i_p) ) * 1.d-04
                   if ( del_param .lt. 1.d-06 ) del_param = 1.d-06

                   params_wiggled(i_p) = params(i_p) + del_param
			
                   call calc_model ( params_wiggled, data_time, 
     +                                  i_pos, i_src,
     +                                  position_wiggled )

                   partls(i_d,i_p)  = (position_wiggled - position) /
     +                                 del_param

                   params_wiggled(i_p) = params(i_p) 

                endif

             enddo

	enddo

c       Rearrange partials if tied together
        call collapse ( num_data, num_params, paramids, 
     +                  partls )

	return
	end
c
	subroutine update_params(num_params,paramids,params,new_params)

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

c       Reset parameters if some are tied together (ID=2)
        call uncollapse ( num_params, paramids, params )

	return
	end
c
	subroutine calc_sigsq (num_data, num_solved_for,
     +			resids, res_err, 
     +			sigsq_old, sigsq, sigma_pdf)

C	Calculate new "sigsq"...

	implicit real*8	(a-h,o-z)

	real*8		resids(1000), res_err(1000)


c       First only RA data
	sigsq_old = sigsq
	sigsq     = 0.d0

	do i_d = 1, num_data, 2
		sigsq = sigsq + ( resids(i_d)/res_err(i_d) )**2
	enddo

	num_deg_freedom = num_data - num_solved_for

	if ( num_deg_freedom .gt. 1 ) then
              sigsq_pdf = sigsq / float(num_deg_freedom/2) 
           else
              sigsq_pdf = 0.d0
        endif

	sigma_pdf = sqrt( sigsq_pdf )

        print *,' sigma_pdf_RA = ',sigma_pdf

c       Next only Dec data
	sigsq_old = sigsq
	sigsq     = 0.d0

	do i_d = 2, num_data, 2
		sigsq = sigsq + ( resids(i_d)/res_err(i_d) )**2
	enddo

	num_deg_freedom = num_data - num_solved_for

	if ( num_deg_freedom .gt. 1 ) then
              sigsq_pdf = sigsq / float(num_deg_freedom/2) 
           else
              sigsq_pdf = 0.d0
        endif

	sigma_pdf = sqrt( sigsq_pdf )

        print *,' sigma_pdf_Dec= ',sigma_pdf

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
      DIMENSION S(1),E(1),X(1)

      character*16      VNAME(1)

      REAL*8    B(401,401),BSTORE(401,401),XIDENT(401,401),
     +          COR(401,401),ERROR(401),STDEV(401),PR(401),
     +          XHAT(401),SNEW(1000),PRODCT(1000),xhat2(401)

      integer   ID(401)

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
c       Reset parameters if some are tied together (ID=2)
        call uncollapse ( numxes, id, xhat2 )
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
 2060         FORMAT (///10X,'THE MULTIPLE CORRELATION COEFFICIENTS'/)

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
