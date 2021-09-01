        program fit_geoblocks_tropos

c       v4: skipped
c       v5: free format reading of delay-rate file, hopefully to
c           be more immune to print format changes in AIPS versions
c       v5a: added LMT coordinates (very approximate; get better ones later)
c       v5b: added AuScope antennas
c       v5c: upped dimension to allow for more than 100 sources in SU table
c       v6:  Added the effects of clock accelerations to ATMOS.FITS file
c            This will replace the normal ATMOS.FITS file with one that
c            has entries more often so as to track clock acceleration.
c       v6a: Found a test for 100 sources in read_geodetic_file;
c            changed to 300
c       v6b: Increased from 1000 to 3000 station-oriented data points
c       v7:  Expanded to allow up to 10 geoblocks
c       v7b: skipped
c       v7c: made changes as for v6c that fixed a bug in the 
c            subroutines make_delzn_2c and make_atmos_fits_w_accel
c            regarding station numbers when some ants are missing
c       v7d: checks for enough data to fit for clock and zenith delays
c       v8:  flexible number of geoblocks, up to 10
c       v8a: added WA station
c            hardwired to read max_blocks (=10) time ranges from the
c            control_file.inp   file

c       This version is for tropospheric data only (eg, dispersive
c       delay corrected X-band or high (>20 GHz) frequency data).

c       Uses Niell (1999) mapping wet function instead of approximation
c       in TMS used in fit_geoblocks_v2c.f

C	Input files:
C		1) "control_file.inp"   provides parameter controls, 
C                                                date & freq
C               2) "station_file.inp"   provides station names/numbers
C               3) "calibrator_file.inp provides 2 file names:
c                                          delay/rate data
c                                          su_table positions

	implicit real*8 ( a-h, o-z )

	real*8      new_params(156),param_sigmas(156),
     +		    params(156), param_incr(156), param_limits(156)
        character*8 parnames(156)
	integer     paramids(156)

     	real*8	    data(20000), model(20000), resids(20000),
     +		    res_err(20000), wtres(20000), gst(20000)

	real*8	    partls(20000,156)

        real*8      ut_min(10), ut_max(10), ut_blk(10)

	integer	    n_stat_a(20000), n_stat_b(20000),
     +              n_source(20000), n_block(20000)

        integer     num_per_stat_geo(12), num_per_stat_blk(12,10), 
     +              list_stations(12)

	character*32 infile, outfile, geodetic_data_file

	logical	    print_cor, print_solution, used_clock_accel

	common /constants/ c,pi,twopi,secday,sidereal_rate,
     +			deg_rad,asec_rad,frequency,wavelength

	real*8	x(12), y(12), z(12), cal_ra(300), cal_dec(300)
	common/geometry/ x,y,z,gast0,cal_ra,cal_dec,num_stations

	luprint = 6
	write (luprint,1100)
1100	format(/30x,'Program fit_geoblocks_tropos_v8a: 2020 Aug 18')

c       Dimension Limits:
	max_num_data	=20000		! max number data equations
	num_params	=  156		! max number of parameters
	num_stations    =   12          ! max number of antennas
        max_blocks      =   10          ! max number of blocks

c	------------------------------------------------------------
c		Constants, Array Geometry, and Control Parameters

	call set_constants 

c       Read "stations_file.inp"...
	call input_stations ( list_stations, n_stats_used )

c       Read "calibrator_file.inp"...
	call input_calibrators ( geodetic_data_file )

c	Read  "control_file.inp" (to get control parameters)...
	call input_controls ( max_blocks,
     +          num_sources,itermax,gain,debug,
     +		data_format, nth_print, 
     +          iyr, mon, iday,
     +          num_blocks, ut_min, ut_max, stm_mid,
     +          delay_err, rate_err, geo_reweight,
     +		num_params, params, paramids, parnames, 
     +		param_incr, param_limits) 


C       --------------------------------------------------------------------
C                    Geodetic-like data input

	infile  = geodetic_data_file           ! rates/delays for many calibrators
c       Read in rate-delay data pairs for many calibrators,
c             and add to end of "data" array
	call read_geodetic_file ( debug, infile, 
     +               max_num_data, num_data,  
     +               max_blocks, num_blocks, 
     +               ut_min, ut_max, ut_blk,
     +               gst, gst_mid,
     +               n_stat_a, n_stat_b, n_source, n_block,
     +               num_per_stat_geo, num_per_stat_blk, 
     +               num_geo_pts, num_geos,
     +               data )

	if ( stm_mid .eq. 0.d0 ) stm_mid = gst_mid

C       -------------------------------------------------------------
C                    Don't solve for parms if too little data 
	i_0 = 0                          ! start of clock offsets
	i_1 = 1 * num_stations           ! start of clock rates
	i_2 = 2 * num_stations           ! start of clock accelerations

        write (luprint,1110)
 1110   format(/'Checking antennas/blocks for enough data...')
	do i_ns = 1, num_stations
	   if ( num_per_stat_geo(i_ns) .lt. 8 ) then
              write (luprint,1120) i_ns, num_per_stat_geo(i_ns)
 1120         format(' Ant #',i2,' has',i2,' delays; >7 required')
	      paramids(i_0 + i_ns) = 0   ! clock offsets
	      paramids(i_1 + i_ns) = 0   ! clock rates
	      paramids(i_2 + i_ns) = 0   ! clock accelerations
	   endif
	enddo

c       Check each station/block for too little data...
 	do i_ns = 1, num_stations
           n_good_blocks = num_blocks
           do i_b = 1, num_blocks
              if ( num_per_stat_blk(i_ns,i_b) .lt. 4 ) then
                 write (luprint,1130) i_ns, i_b,
     +                                num_per_stat_blk(i_ns,i_b)
 1130            format(' Ant #',i2,': block',i2,' has',i2,
     +                  ' delays; >3 required')
                 n_good_blocks = n_good_blocks - 1
                 i_0 = (2 + i_b) * num_stations
                 paramids(i_0 + i_ns) = 0 ! zenith delay for station/block
	      endif
           enddo
c          Check that have enough blocks for various clock parameter fits...
           if ( n_good_blocks .lt. 3 ) then
              if ( paramids(i_2 + i_ns) .gt. 0 ) then
                write (luprint,1140) i_ns 
 1140           format(' Ant #',i2,' has <3 blocks; will not fit accel')
              endif
              paramids(i_2 + i_ns) = 0 ! clock accelerations
           endif
           if ( n_good_blocks .lt. 2 ) then
              if ( paramids(i_1 + i_ns) .gt. 0 ) then
                write (luprint,1150) i_ns 
 1150           format(' Ant #',i2,' has <2 blocks; will not fit rate')
              endif
	      paramids(i_1 + i_ns) = 0   ! clock rates
           endif
           if ( n_good_blocks .lt. 1 ) then
              if ( paramids(i_0 + i_ns) .gt. 0 ) then
                write (luprint,1160) i_ns 
 1160           format(' Ant #',i2,' has <1 block; will not fit clock')
              endif
	      paramids(i_0 + i_ns) = 0   ! clock offsets
           endif
	enddo

C       -------------------------------------------------------------
C                    Count up number of solved-for parameters
	num_params_solved_for = 0
	do i_p = 1, num_params
	   if ( paramids(i_p) .gt. 0 ) then
	      num_params_solved_for  = num_params_solved_for + 1
	   endif
	enddo
	write (6,1200) num_params_solved_for
1200	format(/' Solving for',i3,' parameters',/1x)

C       ------------------------------------------------------------
C                     Print out reference time being used...
	ut_mid_rads = (stm_mid - gast0)/sidereal_rate   ! UT (radians)
	call radians_to_hhmmss (ut_mid_rads, ut_mid_hhmmss)
	write (6,1226) ut_mid_hhmmss
 1226	format(/' Reference UT time for fits (hhmmss.s):',f12.1)

C	------------------------------------------------------------
C			Iterate Least-Squares Fitting

	idebug    = debug + 0.5
	print_cor = .false.
	iterpass  = 0          
	sigsq     = 9999.d0

	do while ( iterpass.le.itermax )

	   if ( iterpass .ge. 1 ) then

	      call calc_partials ( debug,
     +	           num_data,
     +             gst, stm_mid, 
     +             n_stat_a, n_stat_b, n_source, n_block,
     +		   num_params, params, paramids, 
     +		   partls )

c             Before calling l-s-fit, decide what it will print out...
	      if ( iterpass.eq.itermax ) then
c                   Always print out last solution and correlations
	            print_cor = .true.
		    print_solution = .true.
		 else
		    if ( mod(iterpass-1,nth_print).eq.0 .and.
     +                   iterpass .ne. 1 ) then
		         print_solution = .true.
		      else
			 print_solution = .false.
		    endif
	      endif

c             Pass error correction factor to L-S routine
              err_factor = sqrt( float(num_geos)/float(num_geo_pts) )
	      call least_squares_fit( idebug,print_solution,
     +                       print_cor,max_num_data,
     +		             num_params,num_data, err_factor,
     +			     parnames,params,paramids,resids,
     +			     partls,res_err,new_params,param_sigmas )

	      call update_params ( debug, num_params, param_incr, 
     +		             param_limits, gain, params, new_params )

	   endif

	   call calc_residuals ( debug, params, 
     +		              num_data, 
     +		              gst, stm_mid, n_source, 
     +                        n_stat_a, n_stat_b, n_block,  
     +			      data, model, resids )

	   call weight_data (num_data,
     +                       delay_err, rate_err, geo_reweight,
     +                       resids, res_err)

	   call calc_sigsq (num_data, num_params_solved_for,  
     +                  geo_reweight,
     +			resids, res_err, 
     +			sigsq_old, sigsq, sigma_pdf) 

	   if ( iterpass .ge. 1 ) then
	       write (luprint,2500) iterpass
 2500	       format (' End iteration #',i2,/)
	   endif

	   iterpass = iterpass + 1

	enddo

c	==============================================================
c				OUTPUTS
c       --------------------------------------------------------------
c       Make a copy of the "control_file.inp" with all the new 
c       parameter values...call it "control_file.out"
	call new_control_file (parnames,params,paramids,
     +              num_params, num_sources, 
     +              data_format, itermax, gain, debug,
     +              iyr,mon,iday,
     +              max_blocks, ut_min, ut_max, stm_mid,
     +              delay_err, rate_err, geo_reweight )

c       --------------------------------------------------------------
c       Print data/model/residuals...

	if ( debug .ge. 1.d0 ) then
	   call print_all_resids ( data, model, resids,
     +                     res_err, wtres, num_data, lu_print ) 
	endif

c       --------------------------------------------------------------
c       Make ascii files for plotting...
	call geo_plot_files ( num_data, 
     +                 list_stations, n_stats_used,
     +                 n_stat_a, n_stat_b, gst, 
     +                 data, model, resids, res_err )

c       --------------------------------------------------------------
c       Make file of atmospheric corrections in CLCOR/ATMO format
	ut_mid_data_rads = (gst_mid - gast0)/sidereal_rate   ! UT (radians)

        call make_delzn_2c( params, param_sigmas, delay_err,
     +                      ut_mid_rads, ut_blk, 
     +                      num_sources, num_blocks, 
     +                      list_stations, n_stats_used )

c       Check if a clock acceleration was used...if so need frequent
c       entries in ATMOS.FITS to track the changes.
        used_clock_accel = .false.
        do i_s = 1, n_stats_used
           n_s = 2*num_stations + i_s 
           if ( params(n_s) .ne. 0.d0 ) used_clock_accel = .true. 
        enddo

        if ( used_clock_accel ) then
c          Will read back the ATMOS.FILE, and make more frequent entries
           call make_atmos_fits_w_accel( params, ut_mid_rads )
        endif


	stop
	end
c
	subroutine set_constants 

	implicit real*8 (a-h,o-z)

	common /constants/ c,pi,twopi,secday,sidereal_rate,
     +			deg_rad,asec_rad,frequency,wavelength

	c 		= 2.9979248d10			! speed light (cm/sec)
	pi 		= 4.d0*atan(1.d0)		
	twopi		= 2.d0*pi 
	secday		= 86400.d0			! seconds per solar day
	sidereal_rate	= 1.002737923d0			! solar to sidereal conv

	deg_rad  	= pi/180.d0			! convert deg's to rad's
	asec_rad	= deg_rad/3600.d0		! arcsec to radians

	return
	end
c
	subroutine input_stations ( list_stations, n_stats_used ) 

C       Reads "station_file.inp" for station codes; the puts station coords
C            in common/geometry/

	implicit real*8 (a-h,o-z)

	character*32    station_file 

	character*2     stat2

	logical         used

	integer         list_stations(12)

	real*8	x(12), y(12), z(12), cal_ra(300), cal_dec(300)
	common/geometry/ x,y,z,gast0,cal_ra,cal_dec,num_stations


c       --------------------------------------------------------------------------
c                         Get station location info organized
	lu_stat = 9
	station_file = 'station_file.inp'          ! Hardwired station input file
	open (lu_stat, file=station_file, status='old')

	write (6, 3999) station_file
 3999	format(/' "',a15,'" info...')

c       Initialize station coordinates to zero so that we can test when station not used
	do i_st = 1, num_stations
	   x(i_st) = 0.d0
	   y(i_st) = 0.d0
	   z(i_st) = 0.d0
	   list_stations(i_st) = -1
	enddo

	ieof    = 0
	n_stats = 0
	do while ( ieof .ge. 0 )

	   read  (lu_stat, 1000, iostat=ieof)  stat2, i_stat_num
 1000	   format(a2,1x,i2)

	   if ( ieof .ge. 0 ) then

	      used = .false.

c             Check if station number is between 1 and num_stations...
	      if ( i_stat_num .lt. 1   .or.
     +             i_stat_num .gt. num_stations ) then
		 write (6,2000) stat2, i_stat_num
 2000		 format(/' Sub input_stations: Bad stat #: ',a2,i3)
		 stop
	      endif

C	      Set station coordinates in standard *left-handed* geodetic coordinate system:
C             X toward Greenwich longitude; Y toward 90 West longitude; Z toward N pole

	      if ( stat2 .eq. 'FD' ) then
		 x(i_stat_num) = -1324009.0026d0	! FD
		 y(i_stat_num) =  5332182.0834d0	! (meters)
		 z(i_stat_num) =  3231962.4355d0
		 used          = .true.
	      endif

	      if ( stat2 .eq. 'KP' ) then
		 x(i_stat_num) = -1995678.4969d0	! KP
		 y(i_stat_num) =  5037317.8209d0
		 z(i_stat_num) =  3357328.0825d0
		 used          = .true.
	      endif

	      if ( stat2 .eq. 'LA' ) then
		 x(i_stat_num) = -1449752.2303d0	! LA
		 y(i_stat_num) =  4975298.7034d0
		 z(i_stat_num) =  3709123.8860d0
		 used          = .true.
	      endif

	      if ( stat2 .eq. 'OV' ) then
		 x(i_stat_num) = -2409149.9782d0        ! OV
		 y(i_stat_num) =  4478573.3221d0
		 z(i_stat_num) =  3838617.3390d0
		 used          = .true.
	      endif

	      if ( stat2 .eq. 'PT' ) then
		 x(i_stat_num) = -1640953.5776d0	! PT
		 y(i_stat_num) =  5014816.1165d0
		 z(i_stat_num) =  3575411.8292d0
		 used          = .true.
	      endif

	      if ( stat2 .eq. 'BR' ) then
		 x(i_stat_num) = -2112064.8515d0	! BR
		 y(i_stat_num) =  3705356.6129d0
		 z(i_stat_num) =  4726813.7587d0
		 used          = .true.
	      endif

	      if ( stat2 .eq. 'NL' ) then
		 x(i_stat_num) =  -130872.1216d0	! NL
		 y(i_stat_num) =  4762317.2264d0
		 z(i_stat_num) =  4226850.9983d0
		 used          = .true.
	      endif

	      if ( stat2 .eq. 'HN' ) then
		 x(i_stat_num) = +1446375.2506d0	! HN
		 y(i_stat_num) =  4447939.7520d0
		 z(i_stat_num) =  4322306.0648d0
		 used          = .true.
	      endif

	      if ( stat2 .eq. 'SC' ) then
		 x(i_stat_num) = +2607848.6047d0	! SC 
		 y(i_stat_num) =  5488069.8283d0
		 z(i_stat_num) =  1932739.4616d0
		 used          = .true.
	      endif

	      if ( stat2 .eq. 'MK' ) then
		 x(i_stat_num) = -5464074.8410d0	! MK 
		 y(i_stat_num) =  2495249.3104d0
		 z(i_stat_num) =  2148296.7183d0
		 used          = .true.
	      endif

              if ( stat2.eq.'VL' .or. stat2.eq.'Y ' ) then
                 x(i_stat_num) = -1601185.3039d0        ! VLA
                 y(i_stat_num) =  5041977.1810d0
                 z(i_stat_num) =  3554875.6407d0
                 used          = .true.
              endif

              if ( stat2.eq.'EB' .or. stat2.eq.'EF' ) then
                 x(i_stat_num) =  4033947.46164d0       ! Effelsberg
                 y(i_stat_num) =  -486990.51504d0
                 z(i_stat_num) =  4900430.80011d0
                 used          = .true.
              endif

              if ( stat2 .eq. 'GB' ) then
                 x(i_stat_num) =   882589.64360d0       ! Green Bank
                 y(i_stat_num) =  4924872.32087d0
                 z(i_stat_num) =  3943729.36253d0
                 used          = .true.
              endif

              if ( stat2.eq.'JV' .or. stat2.eq.'JB' ) then
                 x(i_stat_num) =  3822626.4970d0        ! Jodrell Bank
                 y(i_stat_num) =   154105.5889d0
                 z(i_stat_num) =  5086486.2618d0
                 used          = .true.
              endif

              if ( stat2 .eq. 'MC' ) then
                 x(i_stat_num) =  4461369.9883d0        ! Medicina
                 y(i_stat_num) =  -919596.8324d0
                 z(i_stat_num) =  4449559.1894d0
                 used          = .true.
              endif

              if ( stat2 .eq. 'NT' ) then
                 x(i_stat_num) =  4934563.1230d0        ! Noto
                 y(i_stat_num) = -1321201.2698d0
                 z(i_stat_num) =  3806484.4778d0
                 used          = .true.
              endif

              if ( stat2 .eq. 'ON' ) then
                 x(i_stat_num) =  3370968.1810d0        ! Onsala 85
                 y(i_stat_num) =  -711464.9170d0
                 z(i_stat_num) =  5349664.1130d0
                 used          = .true.
              endif

              if ( stat2 .eq. 'TR' ) then
                 x(i_stat_num) =  3638558.0000d0        ! Torun
                 y(i_stat_num) = -1221967.0000d0
                 z(i_stat_num) =  5077041.0000d0
                 used          = .true.
              endif

              if ( stat2 .eq. 'WB' ) then
                 x(i_stat_num) =  3828440.6400d0        ! Westerbork
                 y(i_stat_num) =  -445226.0300d0
                 z(i_stat_num) =  5064923.0800d0
                 used          = .true.
              endif

              if ( stat2 .eq. 'HH' ) then
                 x(i_stat_num) =  5085442.7805d0        ! Hartebeestok
                 y(i_stat_num) = -2668263.4908d0
                 z(i_stat_num) = -2768697.0345d0
                 used          = .true.
              endif

              if ( stat2 .eq. 'UR' ) then
                 x(i_stat_num) =   228310.726d0         ! Urumuqui
                 y(i_stat_num) = -4631922.805d0
                 z(i_stat_num) =  4367063.964d0
                 used          = .true.
              endif

              if ( stat2 .eq. 'CM' ) then
                 x(i_stat_num) =  3920354.8000d0        ! Cambridge 32m
                 y(i_stat_num) =    -2545.7000d0
                 z(i_stat_num) =  5014285.0000d0
                 used          = .true.
              endif

              if ( stat2 .eq. 'MH' ) then
                 x(i_stat_num) =  2892579.9681d0        ! Metsahovi
                 y(i_stat_num) = -1311719.0699d0
                 z(i_stat_num) =  5512640.6897d0
                 used          = .true.
              endif

c             Very approximate values for VERA stations...
              if ( stat2 .eq. 'VM' ) then
                 x(i_stat_num) =  -3852030.d0           ! Mizusawa
                 y(i_stat_num) =  -3119313.d0
                 z(i_stat_num) =   4013805.d0
                 used          = .true.
              endif

              if ( stat2 .eq. 'VR' ) then
                 x(i_stat_num) =  -3512768.d0           ! Iriki (Kagoshima)
                 y(i_stat_num) =  -4112922.d0
                 z(i_stat_num) =   3379825.d0
                 used          = .true.
              endif

              if ( stat2 .eq. 'VO' ) then
                 x(i_stat_num) =  -4508500.d0           ! Ogasawara (Chichijima)
                 y(i_stat_num) =  -3459494.d0
                 z(i_stat_num) =   2895551.d0
                 used          = .true.
              endif

              if ( stat2 .eq. 'VS' ) then
                 x(i_stat_num) =  -3215912.d0           ! Ishigaki
                 y(i_stat_num) =  -4858713.d0
                 z(i_stat_num) =   2594166.d0
                 used          = .true.
              endif

              if ( stat2 .eq. 'HH' ) then
                 x(i_stat_num) =  5085442.7845d0        ! Hart
                 y(i_stat_num) = -2668263.4862d0
                 z(i_stat_num) = -2768697.0160d0
                 used          = .true.
              endif

              if ( stat2 .eq. 'AT' ) then
                 x(i_stat_num) = -4752447.5000d0        ! ATCA
                 y(i_stat_num) = -2790326.6000d0
                 z(i_stat_num) = -3200491.2900d0
                 used          = .true.
              endif

              if ( stat2 .eq. 'CD' ) then
                 x(i_stat_num) = -3753440.7000d0        ! Ceduna
                 y(i_stat_num) = -3912708.3000d0
                 z(i_stat_num) = -3348066.9000d0
                 used          = .true.
              endif

              if ( stat2 .eq. 'HO' ) then
                 x(i_stat_num) = -3950236.7341d0        ! Hobart
                 y(i_stat_num) = -2522347.5530d0
                 z(i_stat_num) = -4311562.5434d0
                 used          = .true.
              endif

              if ( stat2 .eq. 'MP' ) then
                 x(i_stat_num) = -4682768.6300d0        ! Mopra 
                 y(i_stat_num) = -2802619.0600d0
                 z(i_stat_num) = -3291759.9000d0
                 used          = .true.
              endif

              if ( stat2 .eq. 'PA' ) then     
                 x(i_stat_num) = -4554232.0016d0        ! Parkes
                 y(i_stat_num) = -2816758.9592d0
                 z(i_stat_num) = -3454035.8457d0
                 used          = .true.
              endif

              if ( stat2 .eq. 'HB' ) then
                 x(i_stat_num) = -3949990.687d0         ! Hobart-12m
                 y(i_stat_num) = -2522421.191d0
                 z(i_stat_num) = -4311708.164d0
                 used          = .true.
              endif

              if ( stat2 .eq. 'KE' ) then
                 x(i_stat_num) = -4147354.453d0         ! Katherine-12m
                 y(i_stat_num) = -4581542.384d0
                 z(i_stat_num) = -1573303.310d0
                 used          = .true.
              endif

              if ( stat2 .eq. 'YG' ) then
                 x(i_stat_num) = -2388896.040d0        ! Yarragedee-12m
                 y(i_stat_num) = -5043349.973d0
                 z(i_stat_num) = -3078590.989d0
                 used          = .true.
              endif

	      if ( stat2 .eq. 'WA' ) then
		 x(i_stat_num) = -5115425.60d0  	! Warkworth 30-m
		 y(i_stat_num) =  -477880.31d0
		 z(i_stat_num) = -3767042.81d0
		 used          = .true.
	      endif

	      if ( stat2 .eq. 'WW' ) then
		 x(i_stat_num) = -5115425.60d0  	! Temporarily using Warkworth 30-m
		 y(i_stat_num) =  -477880.31d0          ! coordinates for their 12-m
		 z(i_stat_num) = -3767042.81d0
		 used          = .true.
	      endif

              if ( stat2 .eq. 'YS' ) then               ! Yebes
                 x(i_stat_num) =  4848761.964d0
                 y(i_stat_num) =  -261484.496d0
                 z(i_stat_num) =  4123084.827d0
                 used          = .true.
              endif

              if ( stat2 .eq. 'YM' ) then
                 x(i_stat_num) = -3502544.259d0         ! Yamaguchi 32
                 y(i_stat_num) = -3950966.397d0
                 z(i_stat_num) =  3566381.165d0
                 used          = .true.
              endif

              if ( stat2 .eq. 'SH' ) then
                 x(i_stat_num) = -2831686.912d0         ! Shanghai
                 y(i_stat_num) =  4675733.665d0 
                 z(i_stat_num) =  3275327.685d0
                 used          = .true.
              endif

             if ( stat2 .eq. 'DA' ) then
                 x(i_stat_num) =  3829087.7497d0        ! DArnhall
                 y(i_stat_num) =   169568.7360d0
                 z(i_stat_num) =  5081082.4658d0
                 used          = .true.
              endif

             if ( stat2 .eq. 'LM' ) then
c                FOLLOWING ARE VERY APPROXIMATE (eg, +/-10,000 m)
                 x(i_stat_num) =  -767029.d0           ! LMT Mexico
                 y(i_stat_num) =  5975413.d0
                 z(i_stat_num) =  2072618.d0
                 used          = .true.
              endif

c             Check if the station code matched one of the known station names...
	      if ( used ) then

		  n_stats = n_stats + 1
		  list_stations(n_stats) = i_stat_num

		  write (6,2008) stat2, i_stat_num,
     +                x(i_stat_num), y(i_stat_num), z(i_stat_num)
 2008		  format(1x,a2,i5,3f15.4)
		else
		  write (6,2010) stat2, i_stat_num
 2010		  format(/' Sub input_stations: Station not a',
     +                   ' recognized 2-char code: ',a2,i3)
	      endif

	   endif

	enddo

        n_stats_used = n_stats

	close (unit=lu_stat)

	return
	end
c
	subroutine input_calibrators ( geodetic_data_file )

c       Reads "calibrator_file.inp" for calibrator positions.  
c       File is in correlator job (job____.fx) format 
c       Info is put in common/geometry/

	implicit real*8 (a-h,o-z)

	character*32    calibrator_file, geodetic_data_file, su_file

	common /constants/ c,pi,twopi,secday,sidereal_rate,
     +			deg_rad,asec_rad,frequency,wavelength

	real*8	x(12), y(12), z(12), cal_ra(300), cal_dec(300)
	common/geometry/ x,y,z,gast0,cal_ra,cal_dec,num_stations

	character*8  cal_name
	character*8  epoch

        max_num_sources = 300     ! dimension limit for cal_ra, cal_dec arrays

c       -----------------------------------------------------------------------
c       Input calibrator data file name and "SU" positions for geodetic like data

	lu_cal          = 9
	calibrator_file = 'calibrator_file.inp'          ! Hardwired source input file
	open (lu_cal, file=calibrator_file, status='old')
	write (6, 3999) calibrator_file
 3999	format(/' "',a15,'" info...')

c       First, read in file name for geodetic-like (rates, multi-band delays) data
	read  (lu_cal, 4020) geodetic_data_file
 4020	format(a32)
	write (6,4020 ) geodetic_data_file

c       Next, read in file name for source table
	read  (lu_cal, 4020) su_file
	write (6,4020 ) su_file

	close (unit=lu_cal)


c       ------------------------------------------------------------------
c       Next read source coordinates from su-file


	open (lu_cal, file=su_file, status='old')

	n_found = 0
	n_cal   = 0
	ieof    = 0
	do while ( ieof .ge. 0 )

c    Required format for SU_file:  in PRTAN use the following
c    inext      'su'
c    box         1 2 11 12 13 14 15
c    Will use free-format reads to avoid changes in AIPS PRTAN formatting
 
	   read  (lu_cal,*,iostat=ieof) irow,n_cal,cal_name,
     +                                  ra_deg,dec_deg,xepoch

	   if ( ieof .ge. 0 ) then

c             Select only lines with source coordinates (selected by xepoch=2000)
	      if ( abs(xepoch-2000.d0) .lt. 0.1d0 ) then

c                This line contains source coordinates...
		 backspace (unit=lu_cal)                  ! and re-read line
		 read  (lu_cal,*,iostat=ieof) irow,n_cal,cal_name,
     +                                        ra_deg,dec_deg,xepoch

		 if ( n_cal.ge.1 .and. n_cal.le.max_num_sources ) then
		       cal_ra (n_cal) = ra_deg  * deg_rad            ! radians
		       cal_dec(n_cal) = dec_deg * deg_rad            ! radians
		    else
		       print *,' INVALID SU-TABLE SOURCE NUMBER (1-300 allowed):', n_cal
		       stop
		 endif

c                Convert format to print out...
		 call radians_to_hhmmss (cal_ra(n_cal),  ra_hhmmss )
		 call radians_to_ddmmss (cal_dec(n_cal), dec_ddmmss)

		 write (6,3720) n_cal, cal_name, ra_hhmmss, dec_ddmmss
 3720		 format(' Geodetic Calibrator:',i3,1x,a8,f17.5,f16.4)

		 n_found = n_found + 1

	      endif

	   endif
	   
	enddo

	print *,'Done reading calibrator (SU) file'

	close (unit=lu_cal)

	if ( n_found .lt. 1 ) then
	   print *,' Found no calibrators: check SU-PRTAB format'
	   STOP
	endif

	return
	end
c
        subroutine input_controls( max_blocks,
     +          num_sources,num_iter,gain,debug,
     +		data_format, nth_print,
     +          iyr, mon, iday,
     +          num_blocks, ut_min, ut_max, stm_mid,
     +          delay_err, rate_err, geo_reweight,
     +		num_params, params, paramids, parnames, 
     +		param_incr, param_limits )

        implicit real*8 (a-h,o-z)

	real*8		params(156), param_limits(156), param_incr(156)
	character*8     parnames(156)
	integer		paramids(156)

	real*8          ut_min(10), ut_max(10)

	real*8          del_tau_b1(12), tau_id_b1(12), 
     +                  del_tau_b2(12), tau_id_b2(12),
     +                  del_tau_b3(12), tau_id_b3(12),
     +                  del_tau_b4(12), tau_id_b4(12),
     +                  del_tau_b5(12), tau_id_b5(12),
     +                  del_tau_b6(12), tau_id_b6(12), 
     +                  del_tau_b7(12), tau_id_b7(12),
     +                  del_tau_b8(12), tau_id_b8(12),
     +                  del_tau_b9(12), tau_id_b9(12),
     +                  del_tau_b10(12),tau_id_b10(12),
     +                  delay_offset(12), delay_id(12),
     +                  delay_rate(12), delay_rate_id(12),
     +                  delay_accel(12), delay_accel_id(12)

        character*32    control_file

	common /constants/ c,pi,twopi,secday,sidereal_rate,
     +			deg_rad,asec_rad,frequency,wavelength

	real*8	x(12), y(12), z(12), cal_ra(300), cal_dec(300)
	common/geometry/ x,y,z,gast0,cal_ra,cal_dec,num_stations

        write (6,5010)
 5010   format(//'"control_file" info...')

        control_file = 'control_file.inp'
        open (unit=8, file=control_file, status='old')

        read (8,*) data_format

        read (8,*) xnum_iter
	num_iter = xnum_iter + 0.5

        read (8,*) gain

        read (8,*) debug

C       Date and time info...
	read  (8,*) iyr, mon, iday
c       Check for reasonable values...
        if ( iyr.lt.2000  .or. iyr.gt.2050 ) then
           print *,' Unreasonable year: ',iyr
           stop
        endif
        if ( mon.lt.1  .or. mon.gt.12 ) then
           print *,' Unreasonable month: ',mon
           stop
        endif
        if ( iday.lt.1  .or. iday.gt.31 ) then
           print *,' Unreasonable day: ',iday
           stop
        endif

        call JULDA (iyr, mon, iday, 
     +              date, gmst_hr )
c       Ignore ~1 sec differences between GAST and GMST
	gast0 = gmst_hr * pi / 12.d0               ! [rad] to put in commons
	call radians_to_hhmmss ( gast0, gast0_hhmmss )

	write (6,5000) iyr, mon, iday, gast0_hhmmss
 5000	format(' Date & GAST0:',i5,2i3,5x,f13.4)

	read (8, *) frequency                ! Hz
	if ( frequency .lt. 1.d9  .or.
     +       frequency .gt. 1.d12      ) then
	   write (6,5001) frequency
 5001	   format(/' Invalid frequency: ',1pd12.3,' Hz')
	   stop
	endif
	wavelength = c/frequency		      ! cm
	freq_GHz   = frequency * 1.d-09
	write (6,5002) freq_GHz, wavelength
 5002	format(/' Observing frequency =',f7.3,' GHz;',
     +          ' wavelength =',f7.3,' cm')

c       Data uncertainties and re-weighting flag
	read (8,*) delay_err, rate_err, geo_reweight
        write (6,5003) delay_err, rate_err, geo_reweight
 5003   format(/' Delay unc =',f5.1,' cm; Rate unc =',f6.3,' Hz;',
     +          ' Automatic re-weighting flag =',f3.0)

c	read (8,8010) xnth_print
c	if ( xnth_print.gt.0.99 ) then
c	      nth_print = xnth_print + 0.5
c	   else
c	      nth_print = 1
c	endif
        nth_print = num_iter           ! hardwire to number of iterations

c       Enter referenced time for clocks 
cV8a    but NOT the number of blocks to read in
	read (8,8011) ut_hhmmss
 8011   format(f10.3,i10)

c       Calculate middle GST (if entered mid UT time .ne. 0.0)
	stm_mid = 0.d0
	if ( ut_hhmmss .ne. 0.d0 ) then
	   call hmsrad( ut_hhmmss, ut_mid_rad)
	   stm_mid = gast0 + ut_mid_rad*sidereal_rate ! radians
	endif
        write (6,5004) ut_hhmmss
 5004   format(/' Input UT center time (hhmmss):',f10.1)

c       Set number of allowed blocks 
c       If not specified, defaults to 10 blocks
c       V8a: hardwired it to read max_blocks from control_file.inp
        n_blocks = max_blocks
c        if ( n_blocks .eq. 0 ) n_blocks = 10
        if ( n_blocks .gt. max_blocks) then
           print*,' Request for ',n_blocks,' exceeds maximum of ',
     +            max_blocks
        endif

        write (6,5111) n_blocks
 5111   format(/' Will read in time ranges and zenith delay info for',
     +          i3,' geoblocks')

c       Enter UT time ranges for n_blocks...
        num_blocks = 0
        do i_b = 1, n_blocks
           read (8,8010) ut_min(i_b), ut_max(i_b)
           write (6,5005)i_b, ut_min(i_b), ut_max(i_b)
 5005      format('   Block ',i2,': ut_min and max (days):',2f7.3)
           if ( ut_min(i_b) .gt. ut_max(i_b) ) then
              print *,' *** Invalid UT time range for block!'
              stop
           endif
           if ( ut_max(i_b) .gt. 0.d0 ) num_blocks = num_blocks + 1
        enddo
        write (6,5006) num_blocks
 5006   format(' Entered time ranges for',i3,' geoblocks')

c       FOR VLBA data, SNR not currently passed to input files.
c       ignore this parameter...
c        read (8,8010) snr_limit
c        print *, 'SNR lower limit: ', snr_limit
c        snr_limit = 0.d0       

c       Station clocks...
	do i=1, num_stations                     ! Station m-b delay offsets
	        read (8,8010) delay_offset(i), delay_id(i)
	enddo
	do i=1, num_stations                     ! Station m-b delay rates
	        read (8,8010) delay_rate(i), delay_rate_id(i)
	enddo
	do i=1, num_stations                     ! Station m-b delay accels
	        read (8,8010) delay_accel(i), delay_accel_id(i)
	enddo

c       Station zenith delays
c       First zero values in case not used...
        do i = 1, num_stations
           del_tau_b1(i) = 0.d0
           tau_id_b1(i)  = 0.d0
           del_tau_b2(i) = 0.d0
           tau_id_b2(i)  = 0.d0
           del_tau_b3(i) = 0.d0
           tau_id_b3(i)  = 0.d0
           del_tau_b4(i) = 0.d0
           tau_id_b4(i)  = 0.d0
           del_tau_b5(i) = 0.d0
           tau_id_b5(i)  = 0.d0
           del_tau_b6(i) = 0.d0
           tau_id_b6(i)  = 0.d0
           del_tau_b7(i) = 0.d0
           tau_id_b7(i)  = 0.d0
           del_tau_b8(i) = 0.d0
           tau_id_b8(i)  = 0.d0
           del_tau_b9(i) = 0.d0
           tau_id_b9(i)  = 0.d0
           del_tau_b10(i) = 0.d0
           tau_id_b10(i)  = 0.d0
        enddo

c       Now read from control file...
        if ( n_blocks .ge. 1 ) then
c          Time block 1
           do i=1, num_stations ! Station zenith delay offsets
        	read (8,8010) del_tau_b1(i), tau_id_b1(i)
           enddo
        endif

        if ( n_blocks .ge. 2 ) then
c          Time block 2
           do i=1, num_stations ! Station zenith delay offsets
        	read (8,8010) del_tau_b2(i), tau_id_b2(i)
           enddo
        endif

        if ( n_blocks .ge. 3 ) then
c          Time block 3
           do i=1, num_stations ! Station zenith delay offsets
        	read (8,8010) del_tau_b3(i), tau_id_b3(i)
           enddo
        endif

        if ( n_blocks .ge. 4 ) then
c          Time block 4
           do i=1, num_stations ! Station zenith delay offsets
        	read (8,8010) del_tau_b4(i), tau_id_b4(i)
           enddo
        endif

        if ( n_blocks .ge. 5 ) then
c          Time block 5
           do i=1, num_stations ! Station zenith delay offsets
        	read (8,8010) del_tau_b5(i), tau_id_b5(i)
           enddo
        endif

        if ( n_blocks .ge. 6 ) then
c          Time block 6
           do i=1, num_stations ! Station zenith delay offsets
        	read (8,8010) del_tau_b6(i), tau_id_b6(i)
           enddo
        endif

        if ( n_blocks .ge. 7 ) then
c          Time block 7
           do i=1, num_stations ! Station zenith delay offsets
        	read (8,8010) del_tau_b7(i), tau_id_b7(i)
           enddo
        endif

        if ( n_blocks .ge. 8 ) then
c          Time block 8
           do i=1, num_stations ! Station zenith delay offsets
        	read (8,8010) del_tau_b8(i), tau_id_b8(i)
           enddo
        endif

        if ( n_blocks .ge. 9 ) then
c          Time block 9
           do i=1, num_stations ! Station zenith delay offsets
        	read (8,8010) del_tau_b9(i), tau_id_b9(i)
           enddo
        endif

        if ( n_blocks .ge. 10 ) then
c          Time block 10
           do i=1, num_stations ! Station zenith delay offsets
        	read (8,8010) del_tau_b10(i), tau_id_b10(i)
           enddo
        endif

 8010   format(8f10.3)

c       Done reading...
	close (unit=8)

c       ====================================================================
c                    Transfer values to params array...

c       ----------------
c       Clock parameters...
	i_0  = 0
c	Set parameters 1 through 12 from clock errors for 12 stations
	do i = 1, num_stations
		params(i_0 + i)       = delay_offset(i)	! cm
		param_limits(i_0 + i) = 9.9d9		! cm
		param_incr(i_0 + i)   = 9.9d9		! cm
	enddo
 
	i_0  = i_0 + num_stations
c	Set parameters 13 through 24 from clock rates for 12 stations
	do i = 1, num_stations
		params(i_0 + i)       = delay_rate(i)	! cm/hr
		param_limits(i_0 + i) = 9.9d9		! cm/hr
		param_incr(i_0 + i)   = 9.9d9		! cm/hr
	enddo

	i_0  = i_0 + num_stations
c	Set parameters 25 through 36 from clock accels for 12 stations
	do i = 1, num_stations
		params(i_0 + i)       = delay_accel(i)	! cm/hr^2
		param_limits(i_0 + i) = 9.9d9		! cm/hr^2
		param_incr(i_0 + i)   = 9.9d9		! cm/hr^2
	enddo

c       -------------------
c       Zenith delay errors for up to 10 blocks of data...
c	Set parameters 37 through 48 from block 1 zenith delays for 12 stations
        i_0  = i_0 + num_stations
        do i = 1, num_stations
		params(i_0 + i)       = del_tau_b1(i)	! cm
		param_limits(i_0 + i) = 250.d0		! cm
		param_incr(i_0 + i)   = 100.d0		! cm
        enddo

c	Set parameters 49 through 60 from block 2 zenith delays for 12 stations
        i_0  = i_0 + num_stations
        do i = 1, num_stations
		params(i_0 + i)       = del_tau_b2(i)	! cm
		param_limits(i_0 + i) = 250.d0		! cm
		param_incr(i_0 + i)   = 100.d0		! cm
        enddo

c	Set parameters 61 through 72 from block 3 zenith delays for 12 stations
        i_0  = i_0 + num_stations
        do i = 1, num_stations
		params(i_0 + i)       = del_tau_b3(i)	! cm
		param_limits(i_0 + i) = 250.d0		! cm
		param_incr(i_0 + i)   = 100.d0		! cm
        enddo

c	Set parameters 73 through 84 from block 4 zenith delays for 12 stations
        i_0  = i_0 + num_stations
        do i = 1, num_stations
		params(i_0 + i)       = del_tau_b4(i)	! cm
		param_limits(i_0 + i) = 250.d0		! cm
		param_incr(i_0 + i)   = 100.d0		! cm
        enddo

c	Set parameters 85 through 96 from block 5 zenith delays for 12 stations
        i_0  = i_0 + num_stations
        do i = 1, num_stations
		params(i_0 + i)       = del_tau_b5(i)	! cm
		param_limits(i_0 + i) = 250.d0		! cm
		param_incr(i_0 + i)   = 100.d0		! cm
        enddo

c	Set parameters 97 through 108 from block 6 zenith delays for 12 stations
        i_0  = i_0 + num_stations
        do i = 1, num_stations
		params(i_0 + i)       = del_tau_b6(i)	! cm
		param_limits(i_0 + i) = 250.d0		! cm
		param_incr(i_0 + i)   = 100.d0		! cm
        enddo

c	Set parameters 109 through 120 from block 7 zenith delays for 12 stations
        i_0  = i_0 + num_stations
        do i = 1, num_stations
		params(i_0 + i)       = del_tau_b7(i)	! cm
		param_limits(i_0 + i) = 250.d0		! cm
		param_incr(i_0 + i)   = 100.d0		! cm
        enddo

c	Set parameters 121 through 132 from block 8 zenith delays for 12 stations
        i_0  = i_0 + num_stations
        do i = 1, num_stations
		params(i_0 + i)       = del_tau_b8(i)	! cm
		param_limits(i_0 + i) = 250.d0		! cm
		param_incr(i_0 + i)   = 100.d0		! cm
        enddo

c	Set parameters 133 through 144 from block 9 zenith delays for 12 stations
        i_0  = i_0 + num_stations
        do i = 1, num_stations
		params(i_0 + i)       = del_tau_b9(i)	! cm
		param_limits(i_0 + i) = 250.d0		! cm
		param_incr(i_0 + i)   = 100.d0		! cm
        enddo

c	Set parameters 144 through 156 from block 10 zenith delays for 12 stations
        i_0  = i_0 + num_stations
        do i = 1, num_stations
		params(i_0 + i)       = del_tau_b10(i)	! cm
		param_limits(i_0 + i) = 250.d0		! cm
		param_incr(i_0 + i)   = 100.d0		! cm
        enddo

c       --------------------------------------------------------------
c                  Set parameter "IDs"...

        do i = 1, num_params
                paramids(i) = 0              ! first set all to "not solved for"
        enddo

	i_0  = 0
	i_1  = i_0 + num_stations
	i_2  = i_1 + num_stations
	i_3  = i_2 + num_stations
	i_4  = i_3 + num_stations
	i_5  = i_4 + num_stations
	i_6  = i_5 + num_stations
	i_7  = i_6 + num_stations
	i_8  = i_7 + num_stations
	i_9  = i_8 + num_stations
	i_10 = i_9 + num_stations
	i_11 = i_10+ num_stations
	i_12 = i_11+ num_stations

	do i = 1, num_stations	                     ! station parameters
	   if ( delay_id(i)      .ne. 0.0 ) paramids(i_0+i) = 1 ! solve for
	   if ( delay_rate_id(i) .ne. 0.0 ) paramids(i_1+i) = 1
	   if ( delay_accel_id(i).ne. 0.0 ) paramids(i_2+i) = 1
	   if ( tau_id_b1(i)     .ne. 0.0 ) paramids(i_3+i) = 1 
	   if ( tau_id_b2(i)     .ne. 0.0 ) paramids(i_4+i) = 1 
	   if ( tau_id_b3(i)     .ne. 0.0 ) paramids(i_5+i) = 1 
	   if ( tau_id_b4(i)     .ne. 0.0 ) paramids(i_6+i) = 1 
	   if ( tau_id_b5(i)     .ne. 0.0 ) paramids(i_7+i) = 1 
	   if ( tau_id_b6(i)     .ne. 0.0 ) paramids(i_8+i) = 1 
	   if ( tau_id_b7(i)     .ne. 0.0 ) paramids(i_9+i) = 1 
	   if ( tau_id_b8(i)     .ne. 0.0 ) paramids(i_10+i)= 1 
	   if ( tau_id_b9(i)     .ne. 0.0 ) paramids(i_11+i)= 1 
	   if ( tau_id_b10(i)    .ne. 0.0 ) paramids(i_12+i)= 1 
	enddo

c       ------------------------------------------------------------
c                       Set parameter labels

	parnames(1) = 'A1 C(cm)'
	parnames(2) = 'A2 C(cm)'
	parnames(3) = 'A3 C(cm)'
	parnames(4) = 'A4 C(cm)'
	parnames(5) = 'A5 C(cm)'
	parnames(6) = 'A6 C(cm)'
	parnames(7) = 'A7 C(cm)'
	parnames(8) = 'A8 C(cm)'
	parnames(9) = 'A9 C(cm)'
	parnames(10)= 'A10C(cm)'
	parnames(11)= 'A11C(cm)'
	parnames(12)= 'A12C(cm)'

	parnames(13)= 'A1 Rcm/h'
	parnames(14)= 'A2 Rcm/h'
	parnames(15)= 'A3 Rcm/h'
	parnames(16)= 'A4 Rcm/h'
	parnames(17)= 'A5 Rcm/h'
	parnames(18)= 'A6 Rcm/h'
	parnames(19)= 'A7 Rcm/h'
	parnames(20)= 'A8 Rcm/h'
	parnames(21)= 'A9 Rcm/h'
	parnames(22)= 'A10Rcm/h'
	parnames(23)= 'A11Rcm/h'
	parnames(24)= 'A12Rcm/h'

	parnames(25)= 'A1 Ac/h2'
	parnames(26)= 'A2 Ac/h2'
	parnames(27)= 'A3 Ac/h2'
	parnames(28)= 'A4 Ac/h2'
	parnames(29)= 'A5 Ac/h2'
	parnames(30)= 'A6 Ac/h2'
	parnames(31)= 'A7 Ac/h2'
	parnames(32)= 'A8 Ac/h2'
	parnames(33)= 'A9 Ac/h2'
	parnames(34)= 'A10Ac/h2'
	parnames(35)= 'A11Ac/h2'
	parnames(36)= 'A12Ac/h2'

	parnames(37)= 'A1 B1 cm'
	parnames(38)= 'A2 B1 cm'
	parnames(39)= 'A3 B1 cm'
	parnames(40)= 'A4 B1 cm'
	parnames(41)= 'A5 B1 cm'
	parnames(42)= 'A6 B1 cm'
	parnames(43)= 'A7 B1 cm'
	parnames(44)= 'A8 B1 cm'
	parnames(45)= 'A9 B1 cm'
	parnames(46)= 'A10B1 cm'
	parnames(47)= 'A11B1 cm'
	parnames(48)= 'A12B1 cm'

	parnames(49)= 'A1 B2 cm'
	parnames(50)= 'A2 B2 cm'
	parnames(51)= 'A3 B2 cm'
	parnames(52)= 'A4 B2 cm'
	parnames(53)= 'A5 B2 cm'
	parnames(54)= 'A6 B2 cm'
	parnames(55)= 'A7 B2 cm'
	parnames(56)= 'A8 B2 cm'
	parnames(57)= 'A9 B2 cm'
	parnames(58)= 'A10B2 cm'
	parnames(59)= 'A11B2 cm'
	parnames(60)= 'A12B2 cm'

	parnames(61)= 'A1 B3 cm'
	parnames(62)= 'A2 B3 cm'
	parnames(63)= 'A3 B3 cm'
	parnames(64)= 'A4 B3 cm'
	parnames(65)= 'A5 B3 cm'
	parnames(66)= 'A6 B3 cm'
	parnames(67)= 'A7 B3 cm'
	parnames(68)= 'A8 B3 cm'
	parnames(69)= 'A9 B3 cm'
	parnames(70)= 'A10B3 cm'
	parnames(71)= 'A11B3 cm'
	parnames(72)= 'A12B3 cm'

	parnames(73)= 'A1 B4 cm'
	parnames(74)= 'A2 B4 cm'
	parnames(75)= 'A3 B4 cm'
	parnames(76)= 'A4 B4 cm'
	parnames(77)= 'A5 B4 cm'
	parnames(78)= 'A6 B4 cm'
	parnames(79)= 'A7 B4 cm'
	parnames(80)= 'A8 B4 cm'
	parnames(81)= 'A9 B4 cm'
	parnames(82)= 'A10B4 cm'
	parnames(83)= 'A11B4 cm'
	parnames(84)= 'A12B4 cm'

	parnames(85)= 'A1 B5 cm'
	parnames(86)= 'A2 B5 cm'
	parnames(87)= 'A3 B5 cm'
	parnames(88)= 'A4 B5 cm'
	parnames(89)= 'A5 B5 cm'
	parnames(90)= 'A6 B5 cm'
	parnames(91)= 'A7 B5 cm'
	parnames(92)= 'A8 B5 cm'
	parnames(93)= 'A9 B5 cm'
	parnames(94)= 'A10B5 cm'
	parnames(95)= 'A11B5 cm'
	parnames(96)= 'A12B5 cm'

	parnames(97)= 'A1 B6 cm'
	parnames(98)= 'A2 B6 cm'
	parnames(99)= 'A3 B6 cm'
	parnames(100)= 'A4 B6 cm'
	parnames(101)= 'A5 B6 cm'
	parnames(102)= 'A6 B6 cm'
	parnames(103)= 'A7 B6 cm'
	parnames(104)= 'A8 B6 cm'
	parnames(105)= 'A9 B6 cm'
	parnames(106)= 'A10B6 cm'
	parnames(107)= 'A11B6 cm'
	parnames(108)= 'A12B6 cm'

	parnames(109)= 'A1 B7 cm'
	parnames(110)= 'A2 B7 cm'
	parnames(111)= 'A3 B7 cm'
	parnames(112)= 'A4 B7 cm'
	parnames(113)= 'A5 B7 cm'
	parnames(114)= 'A6 B7 cm'
	parnames(115)= 'A7 B7 cm'
	parnames(116)= 'A8 B7 cm'
	parnames(117)= 'A9 B7 cm'
	parnames(118)= 'A10B7 cm'
	parnames(119)= 'A11B7 cm'
	parnames(120)= 'A12B7 cm'

	parnames(121)= 'A1 B8 cm'
	parnames(122)= 'A2 B8 cm'
	parnames(123)= 'A3 B8 cm'
	parnames(124)= 'A4 B8 cm'
	parnames(125)= 'A5 B8 cm'
	parnames(126)= 'A6 B8 cm'
	parnames(127)= 'A7 B8 cm'
	parnames(128)= 'A8 B8 cm'
	parnames(129)= 'A9 B8 cm'
	parnames(130)= 'A10B8 cm'
	parnames(131)= 'A11B8 cm'
	parnames(132)= 'A12B8 cm'

	parnames(133)= 'A1 B9 cm'
	parnames(134)= 'A2 B9 cm'
	parnames(135)= 'A3 B9 cm'
	parnames(136)= 'A4 B9 cm'
	parnames(137)= 'A5 B9 cm'
	parnames(138)= 'A6 B9 cm'
	parnames(139)= 'A7 B9 cm'
	parnames(140)= 'A8 B9 cm'
	parnames(141)= 'A9 B9 cm'
	parnames(142)= 'A10B9 cm'
	parnames(143)= 'A11B9 cm'
	parnames(144)= 'A12B9 cm'

	parnames(145)= 'A1 B10cm'
	parnames(146)= 'A2 B10cm'
	parnames(147)= 'A3 B10cm'
	parnames(148)= 'A4 B10cm'
	parnames(149)= 'A5 B10cm'
	parnames(150)= 'A6 B10cm'
	parnames(151)= 'A7 B10cm'
	parnames(152)= 'A8 B10cm'
	parnames(153)= 'A9 B10cm'
	parnames(154)= 'A10B10cm'
	parnames(155)= 'A11B10cm'
	parnames(156)= 'A12B10cm'

        return
        end
c
	subroutine read_geodetic_file ( debug, infile, 
     +                  max_num_data, num_data,  
     +                  max_blocks, num_blocks, 
     +                  ut_min, ut_max, ut_blk,
     +	                gst, gst_mid, 
     +                  n_stat_a, n_stat_b, n_source, n_block,
     +			num_per_stat_geo, num_per_stat_blk, 
     +                  num_geo_pts, num_geos,
     +                  data )

C	Read in geodetic-like data for many calibrators...

	implicit real*8 (a-h,o-z)

	real*8       gst(20000), data(20000)
	integer	     n_stat_a(20000), n_stat_b(20000), 
     +               n_source(20000), n_block(20000)

        real*8       ut_min(10), ut_max(10), ut_blk(10)
        integer      n_in_blk(10)

	character*32 infile
	character*80 comment
	character*4  del_check,frate_check
	character*1  flag,blank

	common /constants/ c,pi,twopi,secday,sidereal_rate,
     +			deg_rad,asec_rad,frequency,wavelength

	integer      num_per_stat_geo(12), num_per_stat_blk(12,10)
	real*8       gst_geo(3000),delay_mb(3000),rate(3000)
	integer      n_stat_geo(3000), n_src_geo(3000), n_blk_geo(3000)
	data         max_n_geo/3000/         ! Subroutine internal dim limit

	logical      used_stat

	real*8	x(12), y(12), z(12), cal_ra(300), cal_dec(300)
	common/geometry/ x,y,z,gast0,cal_ra,cal_dec,num_stations

	lu_in	= 8
	blank	= ' '

c       -----------------------------------------------------------------
c                    Read in station oriented geodetic-like data

c	Open geodetic-like data file...
	write (6,1010) infile
1010	format(/' Openning geodetic data file "',a32,'"')
	open (unit=lu_in, file=infile, status='old')

	read  (lu_in,1001)                     ! Skip 1st line 
1001	format(a80)
	do i=1, 5                              ! Some comment lines
		read  (lu_in,1001) comment
		write (6,1001) comment
	enddo

	ieof  = 0
	n_geo = 0
        i_row      = 0
        i_row_next = 0

        do i_b = 1, max_blocks
           ut_blk(i_b)   = 0.
           n_in_blk(i_b) = 0
        enddo

	do while ( ieof.ge.0 )

	   i_ant = 0             ! to enforce skipping "page-break" lines,
	   i_src = 0             !   which seem to skip the reading of these 

c          Read STATION-ORIENTED multi-band delays and fringe rates...

c          Check for a next row number
           do while ( i_row_next.ne.i_row+1 .and. ieof.ge.0 )
              read (lu_in,*,iostat=ieof) i_row_next
           enddo
           if ( ieof.ge.0 ) then

c            Check for 'INDE' values (which can read as 0's)...
             backspace (unit=lu_in)
	     read (lu_in,*) i_row,ut_days,i_src,i_ant,
     +                      del_check,frate_check

             if ( del_check.ne.'INDE' .and. frate_check.ne.'INDE' ) then 

c               Have a potentially good line of data...backspace and re-read
                backspace (unit=lu_in)
	        read (lu_in,*) i_row,ut_days,
     +                         i_src,i_ant,del_mb,frate

c               Check for other ills before using...
	        if (  i_row.le.3000                          .and.
     +               (i_ant.ge.1 .and. i_ant.le.num_stations).and.
     +                x(i_ant).ne.0.d0                       .and.
     +               (i_src.ge.1 .and. i_src.le.300)         .and.
     +                n_geo.lt.max_n_geo )   then

 		   if ( debug .ge. 1.0 ) then
		      write (6,1002) i_row,ut_days,
     +                               i_src,i_ant,del_mb,frate
 1002                 format(' Delay-rate info:',i4,f8.5,2i5,
     +                       1pd12.3,d12.3)
		   endif

c                  Now, check if within a desired block UT time range...
                   n_b = 0
                   do i_b = 1, num_blocks
                      if ( ut_days.ge.ut_min(i_b) .and. 
     +                     ut_days.le.ut_max(i_b)       ) then
c                        Accept data and store block number
                         n_b = i_b
                      endif
                   enddo

                   if ( n_b .gt. 0 ) then

		     n_geo = n_geo + 1

  		     n_stat_geo(n_geo) = i_ant     ! Station number
		     n_src_geo(n_geo)  = i_src     ! Cal source number
                     n_blk_geo(n_geo)  = n_b       ! Geodetic block number

		     ut_rads = ut_days * twopi
		     gst_geo(n_geo) = gast0 + ut_rads*sidereal_rate ! GST (rads)

                     ut_blk(n_b)    = ut_blk(n_b)  + ut_rads
                     n_in_blk(n_b)  = n_in_blk(n_b)+ 1

		     delay_mb(n_geo)= del_mb * c                    ! MBdelay (cm)
		     rate(n_geo)    = frate  * frequency            ! fringe rate (Hz)

		   endif            ! data accepted (within a block)
 
	         endif		! "values reasonable" check

	     endif        ! 'INDE' check

           endif      ! eof check

	enddo

	num_geo_pts = n_geo
	write (6,2000) num_geo_pts
 2000	format(' Read in',i7,' station-oriented delay-rate pairs',
     +         ' within specified block UT ranges')

	close (unit=lu_in)

c       ----------------------------------------------------------------
c                  Calculate center times for each block

        do i_b = 1, max_blocks
           if ( n_in_blk(i_b) .ge. 1 ) then
              ut_blk(i_b) = ut_blk(i_b) / float( n_in_blk(i_b) )
           endif
        enddo

c       ----------------------------------------------------------------
c                  Make baseline data from station oriented data

	do n_s = 1, num_stations
	   num_per_stat_geo(n_s) = 0
           do n_b = 1, max_blocks
              num_per_stat_blk(n_s,n_b) = 0
           enddo
	enddo
	gst_mid = 0.d0

	t_diff_max = (10.d0/86400.d0)*twopi              ! 10 sec in radians

        num_data = 0
	n_1 = num_geo_pts - 1
	do n_geo = 1, n_1

	   n_src = n_src_geo(n_geo)
	   gst_n = gst_geo(n_geo)
           n_blk = n_blk_geo(n_geo)

	   m_0 = n_geo + 1
	   do m_geo = m_0, num_geo_pts

	      m_src = n_src_geo(m_geo)
	      gst_m = gst_geo(m_geo)

c             Must match source and time to make baseline data 
	      if ( n_src .eq. m_src   .and. 
     +             abs(gst_n - gst_m) .lt. t_diff_max ) then

c                Put baseline "n_geo,m_geo" geodetic like data in large "data" array
c                First the multi-band delay...

		 num_data = num_data + 1

		 n_source(num_data) = n_src
                 n_block(num_data)  = n_blk

		 i_a =  n_stat_geo(n_geo)
		 i_b =  n_stat_geo(m_geo)
		 n_stat_a(num_data) = i_a
		 n_stat_b(num_data) = i_b

c                Count number of delay-rate PAIRS here
		 num_per_stat_geo(i_a) = num_per_stat_geo(i_a) + 1
		 num_per_stat_geo(i_b) = num_per_stat_geo(i_b) + 1
		 num_per_stat_blk(i_a,n_blk)=num_per_stat_blk(i_a,n_blk)+1
		 num_per_stat_blk(i_b,n_blk)=num_per_stat_blk(i_b,n_blk)+1

		 gst(num_data)      = gst_n
		 gst_mid            = gst_mid + gst_n

c                Form delay data from station tables (sense "B minus A")
		 data(num_data)     = delay_mb(m_geo) - delay_mb(n_geo)   ! cm

c                Next the fringe rate...

		 num_data = num_data + 1

		 n_source(num_data) = n_src
                 n_block(num_data)  = n_blk

		 n_stat_a(num_data) = i_a
		 n_stat_b(num_data) = i_b

		 gst(num_data)      = gst_n

c                Form fringe rate data from station tables (sense "B minus A")
		 data(num_data)     = rate(m_geo) - rate(n_geo)           ! Hz

		 if ( debug .ge. 2 ) then
		  write (6,9000) n_src,n_stat_geo(n_geo),n_stat_geo(m_geo),
     +                          gst_n,data(num_data-1),data(num_data)
 9000		  format(' Geo-data',i3,5x,i2,'-',i2,f10.5,1pd12.3,d12.3)
		 endif

	      endif

	   enddo        ! inner loop

	enddo       ! outer loop

	num_geos = num_data 
	num_geo_pairs = num_geos / 2
	if ( num_geo_pairs .ge. 1 ) gst_mid = gst_mid/num_geo_pairs

	write (6,2010) num_geos
 2010	format(' Generated',i5,' baseline geodetic-like delay/rate pairs')

	return
	end
c
	subroutine calc_refsrc_elev ( debug, i_a, i_b, ut_days, 
     +                  elev_a, elev_b )

	implicit real*8 (a-h,o-z)

	common /constants/ c,pi,twopi,secday,sidereal_rate,
     +			deg_rad,asec_rad,frequency,wavelength

	real*8	x(12), y(12), z(12), cal_ra(300), cal_dec(300)
	common/geometry/ x,y,z,gast0,cal_ra,cal_dec,num_stations


	ut_rads = ut_days * twopi
	stm     = gast0 + ut_rads * sidereal_rate       ! GST (rads)

	ra      = ref_ra_fid				! radians
	dec     = ref_dec_fid				! (xfer to new variables)

	x_a = x(i_a)					! A-station [meters]
	y_a = y(i_a)
	z_a = z(i_a)

	x_b = x(i_b)                                    ! B-station [meters]
	y_b = y(i_b)
	z_b = z(i_b)

c       Check that station coordinates have been entered.
c       If not, do not use this datum (set elevation=-999 so that it will be
c       discarded by the read_file subroutine.
	if ( x(i_a).ne.0.d0 .and. x(i_b).ne.0.d0 ) then

c          A-station elevation
	   stat_eq = sqrt( x_a**2 + y_a**2 )
	   stat_lat = atan2( z_a, stat_eq ) ! station latitude (rad)
	   stat_w_long = atan2( y_a, x_a )  ! station west longitude (rad)

c	   local sidereal time of observation...
	   sha  = stm - ra - stat_w_long    ! source hour angle (rad)

c	   calculate source zenith angle...
	   cos_za = sin(stat_lat)*sin(dec) + 	
     +			cos(stat_lat)*cos(dec)*cos(sha)

	   zenith_angle	= acos( cos_za )/deg_rad
	   elev_a	= 90.d0 - zenith_angle

c          B-station elevation...
	   stat_eq = sqrt( x_b**2 + y_b**2 )
	   stat_lat = atan2( z_b, stat_eq ) ! station latitude (rad)
	   stat_w_long = atan2( y_b, x_b )  ! station west longitude (rad)

c	   local sidereal time of observation...
	   sha  = stm - ra - stat_w_long    ! source hour angle (rad)

c	   calculate source zenith angle...
	   cos_za = sin(stat_lat)*sin(dec) + 	
     +			cos(stat_lat)*cos(dec)*cos(sha)

	   zenith_angle	= acos( cos_za )/deg_rad
	   elev_b	= 90.d0 - zenith_angle

	   if ( debug .gt. 3.d0 ) then
	      write (6,1234) i_a, i_b, ut_days, elev_a, elev_b
 1234	      format(' Sub calc_refsrc_elev:',2i5,f10.5,2f10.1)
	   endif

	 else
c           Flag data as bad...
	    elev_a = -999.d0
	    elev_b = -999.d0
	endif

	return
	end
c
	subroutine calc_residuals ( debug, params, 
     +		              num_data, 
     +		              gst, stm_mid, n_source, 
     +                        n_stat_a, n_stat_b, n_block,  
     +			      data, model, resids )

 
C	Calculates models and resduals for entire data set.
C       This now includes phase data and geodetic-like data.

	implicit real*8 ( a-h,o-z )

	real*8	params(156)
	real*8	gst(20000),data(20000),model(20000),resids(20000)

	integer	n_stat_a(20000), n_stat_b(20000), 
     +          n_block(20000),  n_source(20000)

	common /constants/ c,pi,twopi,secday,sidereal_rate,
     +			deg_rad,asec_rad,frequency,wavelength

	real*8	x(12), y(12), z(12), cal_ra(300), cal_dec(300)
	common/geometry/ x,y,z,gast0,cal_ra,cal_dec,num_stations

       

	n_geos = num_data 
	if ( n_geos .gt. 0 ) then

c          Go through geodetic data in pairs (delay & rate)...
           n_d = 0
	   do n = 1, n_geos, 2

	      n_d = n_d + 1

	      n_geo_src = n_source(n_d)

	      call set_geo_parameters (  debug, params, 
     +                  n_d, num_sources, n_geo_src, n_block,
     +			n_stat_a, n_stat_b, gst,
     +			stm, ra, dec, 
     +			x_a,y_a,z_a, x_b,y_b,z_b,
     +			del_tau_z_a, del_tau_z_b,
     +                  clock_a, clock_b,
     +                  clock_rate_a, clock_rate_b,
     +                  clock_accel_a, clock_accel_b )

	      call calc_geo_model ( debug, stm, stm_mid, 
     +                  ra, dec, 
     +			x_a,y_a,z_a, x_b,y_b,z_b, 
     +			del_tau_z_a, del_tau_z_b,
     +                  clock_a, clock_b,
     +                  clock_rate_a, clock_rate_b,
     +                  clock_accel_a, clock_accel_b,
     +			tau_tot_cm, rate_Hz )

c             Calculate both delay & rate model and residual.
c             First, a multi-band delay datum...
	      model(n_d)  = tau_tot_cm                    ! cm
	      resids(n_d) = data(n_d) - model(n_d)

c             Next, a fringe rate datum...
	      n_d = n_d + 1
	      model(n_d)  = rate_Hz                       ! Hz
	      resids(n_d) = data(n_d) - model(n_d)

	   enddo

	endif

	return
	end
c
	subroutine set_geo_parameters ( debug, params, 
     +                  n_d, num_sources, n_geo_src, n_block,
     +			n_stat_a, n_stat_b, gst,
     +			stm, ra, dec, 
     +			x_a,y_a,z_a, x_b,y_b,z_b,
     +			del_tau_z_a, del_tau_z_b,
     +                  clock_a, clock_b,
     +                  clock_rate_a, clock_rate_b,
     +                  clock_accel_a, clock_accel_b )


c	Given data number (n_d), cal source number (n_geo_src), 
c       total number of phase (num_sources=4), and "params" array, 
c       setup parameters needed for the geodetic model calculations
	
	implicit real*8 ( a-h,o-z )

	real*8		params(156)

	real*8		gst(20000)
	integer		n_stat_a(20000), n_stat_b(20000), 
     +                  n_block(20000)

	common /constants/ c,pi,twopi,secday,sidereal_rate,
     +			deg_rad,asec_rad,frequency,wavelength

	real*8	x(12), y(12), z(12), cal_ra(300), cal_dec(300)
	common/geometry/ x,y,z,gast0,cal_ra,cal_dec,num_stations

	stm	 = gst(n_d)

C       Transfer source positions to pass back to calling routine 
	ra  = cal_ra(n_geo_src)		  ! calibrator position (rad)
	dec = cal_dec(n_geo_src)

	i_stat_a = n_stat_a(n_d)
	i_stat_b = n_stat_b(n_d)

	x_a = x(i_stat_a)					 ! meters
	y_a = y(i_stat_a)
	z_a = z(i_stat_a)

	x_b = x(i_stat_b)
	y_b = y(i_stat_b)
	z_b = z(i_stat_b)

c       Get clock parameters...
	n0  = 0                            ! start of multi-band clock offsets
	n1  = n0 + num_stations            ! start of multi-band clock rates
	n2  = n1 + num_stations            ! start of multi-band clock accelerations

	clock_a         = params(n0+i_stat_a)                    ! cm (delay/clight)
	clock_b         = params(n0+i_stat_b)                    ! cm

	clock_rate_a    = params(n1+i_stat_a)                    ! cm/hr
	clock_rate_b    = params(n1+i_stat_b)                    ! cm/hr

	clock_accel_a   = params(n2+i_stat_a)                    ! cm/hr^2
	clock_accel_b   = params(n2+i_stat_b)                    ! cm/hr^2

c       Get appropriate zenith delay values
        n_blk = n_block(n_d)

	n3  = (2 + n_blk) * num_stations     ! start of block "n_blk" zenith delays
	del_tau_z_a	= params(n3+i_stat_a)	                 ! cm of delay
	del_tau_z_b	= params(n3+i_stat_b)

	return
	end
c
	subroutine calc_geo_model ( debug, stm, stm_mid, 
     +                          ra, dec, 
     +				x_a,y_a,z_a, x_b,y_b,z_b, 
     +				del_tau_z_a, del_tau_z_b,
     +                          clock_a, clock_b,
     +                          clock_rate_a, clock_rate_b, 
     +                          clock_accel_a, clock_accel_b, 
     +				tau_tot_cm, rate_Hz )

	implicit real*8 ( a-h, o-z )

	common /constants/ c,pi,twopi,secday,sidereal_rate,
     +			deg_rad,asec_rad,frequency,wavelength


c       -----------------------------------------------------------
c	Calculate atmospheric delay error for geodetic cal data:

c	For station A...	
        del_tau_dot_a = 0.d0
	call atmosphere_shift (debug, stm,stm_mid, 
     +                  ra,dec, x_a,y_a,z_a,		        
     +		        del_tau_z_a,del_tau_dot_a,
     +			del_tau_a)                                   
	
c	For station B...	
        del_tau_dot_b = 0.d0
	call atmosphere_shift (debug, stm,stm_mid, 
     +                  ra,dec, x_b,y_b,z_b, 
     +			del_tau_z_b,del_tau_dot_b,
     +			del_tau_b)                                  
	
c	Difference for baseline delay...	
	tau_atm = del_tau_b - del_tau_a                   ! radians	
	tau_atm_cm = tau_atm * wavelength/twopi           ! cm

c       Calculate clock offset...
	del_t_ut    = (stm - stm_mid)/sidereal_rate       ! radians of UT
	del_t_hr    = del_t_ut * 12.d0/pi                 ! hours
	tau_clock_a = clock_a + clock_rate_a * del_t_hr  +
     +                0.5d0 * clock_accel_a * del_t_hr**2 ! cm
	tau_clock_b = clock_b + clock_rate_b * del_t_hr  +
     +                0.5d0 * clock_accel_b * del_t_hr**2 ! cm
	tau_clock_cm= tau_clock_b - tau_clock_a

c       Now add the atmospheric and clock delays...
	tau_tot_cm = tau_atm_cm + tau_clock_cm            ! cm

c       -----------------------------------------------------------
c	Now calculate time derivative of atmospheric delay error 
c       for fringe rate data:

	d_stm = twopi/86400.d0    ! +1 sec of time (in radians)
        stm_wiggled = stm + d_stm

c	For station A...	
	call atmosphere_shift (debug, stm_wiggled, stm_mid, 
     +                  ra,dec, x_a,y_a,z_a,		        
     +		        del_tau_z_a,del_tau_dot_a,
     +			del_tau_a_wiggled)                                   
	
c	For station B...	
	call atmosphere_shift (debug, stm_wiggled, stm_mid,
     +                  ra, dec, x_b,y_b,z_b, 
     +			del_tau_z_b, del_tau_dot_b,
     +			del_tau_b_wiggled)                                  
	
c	Difference for baseline delay...	
	tau_atm_wiggled = del_tau_b_wiggled - del_tau_a_wiggled ! rad of phs

c       Numerically evaluate time differential of delay
c          Units are  "radians of turns/1 sec of time"
	rate_radsec = tau_atm_wiggled - tau_atm 

c       Convert to cycles from radians of phase
	atm_rate_Hz = rate_radsec / twopi                     ! Hz

c       The clock drift rate also adds a term to the fringe rate:
	del_clock_rate     = (clock_rate_b + clock_accel_b*del_t_hr) -
     +                       (clock_rate_a + clock_accel_a*del_t_hr)  ! cm/hr
	del_clock_rate_sps = del_clock_rate / c / 3600.d0     ! sec/sec
	clock_rate_Hz      = del_clock_rate_sps * frequency   ! Hz

c       Add the atmospheric and clock rate terms...
	rate_Hz = atm_rate_Hz + clock_rate_Hz

	return 
	end
c
	subroutine atmosphere_shift(debug,stm,stm_mid,ra,dec,	
     +				xm,ym,zm,del_tau_z,del_tau_dot,
     +				tau_atm)

c	Returns excess atmospheric delay in radians of phase shift
c       for given station owing to a vertical path error in cm.  

c       This version assumes the excesss path follows a 
c       Niell (1999. JGR, 101, No. B2, p3227-3246) mapping function
c       assuming the residual delays are from the wet atmosphere.
c       Adopted coefficients of the model are for 30 deg station latitude.


	implicit real*8 (a-h, o-z)

	common /constants/ c,pi,twopi,secday,sidereal_rate,
     +			deg_rad,asec_rad,frequency,wavelength

	real*8	x(12), y(12), z(12), cal_ra(300), cal_dec(300)
	common/geometry/ x,y,z,gast0,cal_ra,cal_dec,num_stations

	rad_hr	= 12.d0/pi 

	stat_eq = sqrt( xm**2 + ym**2 )
	stat_lat = atan2( zm, stat_eq )		! station latitude (rad)
	stat_w_long = atan2( ym, xm )		! station west longitude (rad)

c		local sidereal time of observation...
	sha  = stm - ra - stat_w_long		! source hour angle (rad)

c		calculate source zenith angle...
	cos_za = sin(stat_lat)*sin(dec) + 	
     +			cos(stat_lat)*cos(dec)*cos(sha)

	sec_za = 1.d0 / cos_za                     
	zenith_angle	= acos( cos_za )/deg_rad        ! deg

	elevation	= 90.d0 - zenith_angle          ! deg
        el_rad          = elevation * deg_rad           ! rad

	del_tau_z_total = del_tau_z + del_tau_dot*(stm-stm_mid)*rad_hr    ! cm

        call niell ( el_rad, wmf )

	tau_cm	  = del_tau_z_total / wmf	        ! cm
	tau_waves = tau_cm / wavelength			! wavelengths 
	tau_atm	  = tau_waves * twopi			! radians of phase shift

	if ( debug .gt. 3.d0 ) then
		ut	= (stm-gast0)/(deg_rad*15.d0)	! hrs
		i_ut_hr	= ut
		ut	= (ut - i_ut_hr)*60.d0		! min
		i_ut_min= ut
		ut	= (ut - i_ut_min)*60.d0		! sec
		i_ut_sec= ut + 0.5
		write (6,1200) i_ut_hr,i_ut_min,i_ut_sec,
     +				elevation, sec_za, del_tau_z, tau_waves
1200		format(' atmos_delay: ',
     +                        3i3,f7.1,f8.2,2f9.3)
	endif

	return
	end
c
      subroutine niell ( e, wmf )

c     Niell (1999. JGR, 101, No. B2, p3227-3246) mapping function
c     assuming the residual delays are from the wet atmosphere.
c     Adopted coefficients of the model are for 30 deg station latitude.

c     Input:  e    source elevation (radians)
c     Output: wmf  wet mapping function (dimensionless) 

      implicit real*8 (a-h,o-z)

      a = 5.6795d-04         ! 30 deg latitute values from Niell '99 
      b = 1.5139d-03
      c = 4.6730d-02

      top = 1.d0 / (1.d0 + a /(1.d0 + b/(1.d0 + c)))

      sine = sin(e)
      bot = 1.d0 / (sine + a /(sine + b/(sine + c)))

      wmf = top / bot

      return
      end
c
	subroutine calc_partials ( debug,
     +	           num_data,
     +             gst, stm_mid, 
     +             n_stat_a, n_stat_b, n_source, n_block,
     +		   num_params, params, paramids, 
     +		   partls )

	implicit real*8 ( a-h,o-z )

	real*8	params(156), params_wiggled(156)
	real*8	partls(20000,156), gst(20000)

	integer	n_stat_a(20000), n_stat_b(20000), 
     +          n_source(20000), n_block(20000)
	integer	paramids(156)


C	Set secondary parameter array and zero partials array...
	do i_p = 1, num_params

		params_wiggled(i_p) = params(i_p)

		do i_d = 1, num_data
			partls(i_d,i_p) = 0.d0
		enddo

	enddo

c       ----------------------------------------------------------------
c       Now run through geodetic-like data...
	n_geos = num_data 
        n_d = 0
	if ( n_geos .gt. 0 ) then

c          Will be doing data **pairs**;set up to increment counter by 2's
	   n_d = n_d - 1           

c          Go through pairs (dalay & rate) of data...
	   do n = 1, n_geos, 2

c             Increment data record number by 2's...
	      n_d = n_d + 2

	      n_geo_src = n_source(n_d)

	      call set_geo_parameters (  debug, params, 
     +                  n_d, num_sources, n_geo_src, n_block,
     +			n_stat_a, n_stat_b, gst,
     +			stm, ra, dec, 
     +			x_a,y_a,z_a, x_b,y_b,z_b,
     +			del_tau_z_a, del_tau_z_b,
     +                  clock_a, clock_b,
     +                  clock_rate_a, clock_rate_b,
     +                  clock_accel_a, clock_accel_b )

	      call calc_geo_model ( debug, stm, stm_mid, 
     +                  ra, dec, 
     +		        x_a,y_a,z_a, x_b,y_b,z_b, 
     +			del_tau_z_a, del_tau_z_b,
     +                  clock_a, clock_b,
     +                  clock_rate_a, clock_rate_b,
     +                  clock_accel_a, clock_accel_b,
     +			tau, rate )

	      do i_p = 1, num_params

		 if ( paramids(i_p) .eq. 1 ) then

C                   Change parameters by 1 part in 10^4; but also never
C			by less than a value of 10^-4 (arcsec or cm)...
		    del_param = abs( params(i_p) ) * 1.d-04
		    if ( del_param .lt. 1.d-04 ) del_param = 1.d-04

		       params_wiggled(i_p) = params(i_p) + del_param
			
		       call set_geo_parameters (  debug, 
     +                          params_wiggled, 
     +                          n_d, num_sources, n_geo_src, n_block,
     +                          n_stat_a, n_stat_b, gst,
     +			        stm, ra, dec, 
     +			        x_a,y_a,z_a, x_b,y_b,z_b,
     +			        del_tau_z_a, del_tau_z_b,
     +                          clock_a, clock_b,
     +                          clock_rate_a, clock_rate_b,
     +                          clock_accel_a, clock_accel_b )

		       call calc_geo_model ( debug, stm, stm_mid, 
     +                          ra, dec, 
     +				x_a,y_a,z_a, x_b,y_b,z_b, 
     +				del_tau_z_a, del_tau_z_b,
     +                          clock_a, clock_b,
     +                          clock_rate_a, clock_rate_b,
     +                          clock_accel_a, clock_accel_b,
     +				tau_wiggled, rate_wiggled )

c                       First, calculate multi-band delay partial
			partls(n_d,i_p) = (tau_wiggled-tau)/
     +                                     del_param

c                       Next, calculate rate partial...
			partls(n_d+1,i_p) = (rate_wiggled-rate)/
     +                                       del_param

c                       Go back to original parameter value...
			params_wiggled(i_p)= params(i_p) 

		   endif

		enddo

	     enddo

	endif

	return
	end
c
	subroutine update_params (debug, num_params, param_incr,
     +					param_limits, gain,
     +					params, new_params )

	implicit real*8 (a-h, o-z)

	real*8		params(156), new_params(156)
	real*8		param_incr(156), param_limits(156)
	real*8		param_old(156)

	do j = 1, num_params

		param_old(j) = params(j)

		delta_param = new_params(j) - params(j) 

C		Limit magnitude of parameter change to "param_incr(j)"...
		sign = 1.d0
		if ( delta_param .lt. 0.0 ) sign = -1.d0
		if ( abs(delta_param) .gt. param_incr(j) ) 
     +				delta_param = sign * param_incr(j)

		trial_param = params(j) + gain * delta_param

		if ( trial_param .gt. 0.d0 ) then

			if ( trial_param .gt. param_limits(j) ) then
				params(j) = param_limits(j)
			   else
				params(j) = trial_param
			endif

		   else
			
			if ( abs(trial_param) .gt. param_limits(j) ) then
				params(j) = -param_limits(j)
			   else
				params(j) = trial_param
			endif

		endif		

	enddo

	return
	end
c
	subroutine calc_sigsq (num_data, num_params_solved_for,  
     +                  geo_reweight,
     +			resids, res_err, 
     +			sigsq_old, sigsq, sigma_pdf) 

C	Calculate new "sigsq"...

	implicit real*8	(a-h,o-z)

	real*8	 resids(20000), res_err(20000)

c       -----------------------------------------------------------------
c       For geodetic data...
	n_geos = num_data 	
	if ( n_geos .gt. 1 ) then

	   sigsq     = 0.d0
	   do i_d = 1, n_geos
		sigsq = sigsq + ( resids(i_d)/res_err(i_d) )**2
	   enddo

	   num_deg_freedom = n_geos - num_params_solved_for  ! not rigorous

	   if ( num_deg_freedom .gt. 1 )
     +		sigsq_pdf = sigsq / float(num_deg_freedom) 
	   sigma_pdf = sqrt( sigsq_pdf )

	   write (6,1500) sigma_pdf, num_deg_freedom, geo_reweight
 1500	   format(/' Sigma_pdf =',f10.3,' for geo data (',i6,
     +                   ' degrees of freedom); geo_reweight =',f3.0)

	endif

	return
	end
c
	subroutine weight_data (num_data, 
     +                 delay_err, rate_err, geo_reweight, 
     +                 resids, res_err)

C       Assigns data errors and, if "auto_reweight" flags set,
C       down-weights data whose residuals that are very large

	implicit real*8 (a-h,o-z)

	real*8   resids(1), res_err(1)

	common /constants/ c,pi,twopi,secday,sidereal_rate,
     +			deg_rad,asec_rad,frequency,wavelength


	n_geos = num_data 
	if ( n_geos .gt. 0 ) then

	   n_rewt_delay = 0
	   n_rewt_rate  = 0
           i_d = 0
	   do n = 1, n_geos, 2

c             Decide if using fixed weights or automatically down-weighting 
c             bad data...
	      if ( geo_reweight .ne. 0.d0 ) then

c                Delay datum re-weighted...
		 i_d = i_d + 1
		 wt_err = resids(i_d) / delay_err
		 if ( abs(wt_err) .lt. 8.d0 ) then
		    re_wt = exp((wt_err/4.d0)**2)
		   else
		    re_wt = exp(4.d0)
		 endif
		 res_err(i_d) = delay_err * re_wt                ! cm

c		 if ( re_wt .gt. 3.d0 ) print *,' Re-weighting delay datum ',
c     +                                          i_d,' by factor ',re_wt
		 if ( re_wt .gt. 3.d0 ) n_rewt_delay = n_rewt_delay + 1


c                Rate datum re-weighted
		 i_d = i_d + 1
		 wt_err = resids(i_d) / rate_err
		 if ( abs(wt_err) .lt. 8.d0 ) then
		    re_wt = exp((wt_err/4.d0)**2)
		   else
		    re_wt = exp(4.d0)
		 endif
		 res_err(i_d) = rate_err  * re_wt                ! Hz

c		 if ( re_wt .gt. 3.d0 ) print *,' Re-weighting rate datum  ',
c     +                                          i_d,' by factor ',re_wt
		 if ( re_wt .gt. 3.d0 ) n_rewt_rate = n_rewt_rate + 1

		else

		 i_d = i_d + 1
		 res_err(i_d) = delay_err                        ! cm

		 i_d = i_d + 1
		 res_err(i_d) = rate_err                         ! Hz
	
	      endif

	   enddo

	   print *,' Number delays re-weighted: ',n_rewt_delay
	   print *,' Number rates  re-weighted: ',n_rewt_rate

	endif

	return
	end
c
	subroutine new_control_file (parnames,params,paramids,
     +              num_params, num_sources, 
     +              data_format, itermax, gain, debug,
     +              iyr,mon,iday,
     +              max_blocks, ut_min, ut_max, stm_mid,
     +              delay_err, rate_err, geo_reweight )

c       Makes a new version of the "control_file.inp" called
c       "control_file.out".  This can be renamed to ".inp"
c       and used to start off where the last run ended.

	implicit real*8 ( a-h, o-z )

	real*8		params(156)
        character*8     parnames(156)
	integer		paramids(156)

	real*8          ut_min(10), ut_max(10)

	character*32    control_file, par_description

	common /constants/ c,pi,twopi,secday,sidereal_rate,
     +			deg_rad,asec_rad,frequency,wavelength

	real*8	x(12), y(12), z(12), cal_ra(300), cal_dec(300)
	common/geometry/ x,y,z,gast0,cal_ra,cal_dec,num_stations

        control_file = 'control_file.out'
        open (unit=7, file=control_file, status='unknown')

c       Write out various control parameters...
	par_description = 'data format (0=old; 1=new)'
        write (7,8010) data_format, par_description

	xitermax = itermax
	par_description = 'number interations'
        write (7,8010) xitermax, par_description

	par_description = 'loop gain'
	write (7,8010) gain, par_description

	par_description = 'debug print out level'
        write (7,8010) debug, par_description

	par_description = 'yyyy mm dd'
        write (7,8011) iyr, mon, iday, par_description

	par_description = 'observing freq (Hz)'
        write (7,8012) frequency, par_description

	par_description = 'del err(cm); rate err(Hz); re-wt'
	write (7,8040) delay_err, rate_err, geo_reweight,
     +                 par_description

	par_description = 'Reference UT_mid(hhmmss)'
        ut_mid_rad = ( stm_mid - gast0 )/sidereal_rate    ! radians
        call radians_to_hhmmss ( ut_mid_rad, ut_hhmmss )
        write (7,8013) ut_hhmmss, par_description

        do i_b = 1, max_blocks
           write (par_description,8044) i_b
 8044      format('Block',i2,' UT range (days)')
           write (7,8060) ut_min(i_b), ut_max(i_b), par_description
        enddo

c       ----------------------------------------
c       Write out 3 (clocks) + 10 (geoblocks) groups of station parameters...
	n_p = 0
	do n_param_groups = 1, 13
	   do n=1, num_stations	
	      n_p = n_p + 1
 	      write (7,8030) params(n_p), paramids(n_p),
     +                       parnames(n_p)
	   enddo
	enddo

 8010   format( f10.4,30x,'! ',a32)
 8011   format(i4,2i3,29x,'! ',a32)
 8012   format(1pd10.4,30x,'! ',a32)
 8013   format( f10.3,30x,'! ',a32)
 8020   format(2f10.6,f10.1,10x,'! ',a8,', ',a8)
 8030   format( f10.4,8x,i2,20x,'! ',a8)
 8040   format(3f10.4,10x,'! ',a32)
 8042   format(2f10.4,f10.1,10x,'! ',a32)
 8050   format(4f10.4,    '! ',a32)
 8060   format(2f10.4,20x,'! ',a32)

c       Done writing to new file...
	close (unit=7)

	return
	end
c
	subroutine make_delzn_2c(params,param_sigmas,delay_err,
     +                    ut_mid_rads, ut_blk,
     +                    num_sources, num_blocks, 
     +                    list_stations, n_stats_used )

	implicit real*8 ( a-h, o-z )

	real*8		params(156), param_sigmas(156)

        real*8          ut_blk(10)

        integer         list_stations(12)

	character*32    control_file

        integer         flag(10), sequence(10), use_block_num(12,10)

	common /constants/ c,pi,twopi,secday,sidereal_rate,
     +			deg_rad,asec_rad,frequency,wavelength

	real*8	x(12), y(12), z(12), cal_ra(300), cal_dec(300)
	common/geometry/ x,y,z,gast0,cal_ra,cal_dec,num_stations

	pi = 4.d0 * atan(1.d0)
        rad_to_hr = 12.d0/ pi

        control_file = 'ATMOS.FITS'
        open (unit=7, file=control_file, status='unknown')

c       get maximum station number (important if a station is missing)
        max_stn_num = 0
        do i_s = 1, n_stats_used
           if ( list_stations(i_s) .gt. max_stn_num ) then
              max_stn_num = list_stations(i_s)
           endif
        enddo

c       zero counter for actual number of lines to be in ATMOS.FITS file
c       (v2c drops lines if no geoblock fit and surrounding blocks OK)
	num_lines = 0        

c       set up to deal with missing blocks, or blocks with
c       very poor zenith delay uncertainties (> 1.5*delay_err)...
        do i_s = 1, max_stn_num

           do i_b = 1, num_blocks

c             Check if we have a zenith delay for this block
	      i_p = (2 + i_b)*num_stations + i_s  
	      zen_del = params(i_p)             ! cm
              zen_err = param_sigmas(i_p)       ! cm
              if ( zen_del.eq.0.d0 .or.
     +             zen_err.gt.1.5d0*delay_err ) then
c                   Missing or bad-fitting block for this station
                    sequence(i_b) = 0
                 else
                    sequence(i_b) = 1
              endif

           enddo

c          Examine the "sequence" values and set a "flag" array
c          that says which zenith delay value to use if missing values.
c          Flag values: -1 for no output line in ATMOS.FITS
c                      N>0 to use the N'th block value for this station
           call fill_in_out ( sequence, num_blocks, flag )

c          Count number of lines that will be in ATMOS.FITS file and
c          transfer flag values to 2-d array
           do i_b = 1, num_blocks
              if ( flag(i_b) .gt. 0 ) num_lines = num_lines + 1
              use_block_num(i_s,i_b) = flag(i_b)
           enddo

        enddo

	write (7,1000) num_lines
 1000	format(i5)

c       Cycle through by geodetic blocks
        do i_b = 1, num_blocks

c          Prepare data mid-time for output...
	   ut_hr = ut_blk(i_b) * rad_to_hr
	   id    = ut_hr/24.d0
	   ut_hr = ut_hr - id*24
	   ihr   = ut_hr
	   imin  = (ut_hr - ihr*1.d0) * 60.d0
	   sec   = (ut_hr - ihr*1.0d0 - imin/60.0d0) * 3600.d0

c          Delta-time between middle of data block and reference time 
c             for fit parameters
           del_hrs = ( ut_blk(i_b) - ut_mid_rads )*rad_to_hr   ! hrs

c          Station by station output for CLCOR/ATMO
           do i_s = 1, max_stn_num

c             May not write out station zenith delay if missing "interior" block
c             or extrapolate zenith delay if missing "exterior" block (v2c)
              i_b_use = use_block_num(i_s,i_b)
              if ( i_b_use .gt. 0 ) then

c                Clock values come from "equation"
c                Get clock parameter values ...
                 i_p = i_s
                 clk_del = params(i_p) ! cm

                 i_p = i_p + num_stations
                 clk_del_rate= params(i_p) ! cm/hr

                 clk_del_mid_block = clk_del + clk_del_rate*del_hrs ! cm

c                CLCOR doesn't use clock acceleration, so stating with V6 will update
c                clock and drift rate parameters to implement clock acceleration
c                in a step-wise fashion at geoblock times...
                 i_p = i_p + num_stations
                 clk_acc = params(i_p)                              ! cm/hr^2
                 clk_del_mid_block = clk_del_mid_block +            
     +                               0.5d0*clk_acc*del_hrs**2       ! cm
                 clk_del_rate = clk_del_rate + clk_acc*del_hrs      ! cm/hr

c                Next, get atmospheric parameter values for this block, 
c                which may come from an "interior" block if this one was missing...
                 i_p = (2 + i_b_use)*num_stations + i_s  
                 zen_del = params(i_p)                              ! cm

c                Change units on rates to AIPS CLCOR desired units...
                 clk_del_rate_sps14 = clk_del_rate * 1.d14 / (c*3600.d0) ! (sec/sec)*10^14

c                zen_del_rate_sps14 = zen_del_rate * 1.d14 / (c*3600.d0) ! (sec/sec)*10^14
                 zen_del_rate_sps14 = 0.d0 ! handled in CLCOR by interpolating delays

c                Output parameter values
                 write (7,1100) i_s, id, ihr, imin, sec, 
     +                          zen_del, clk_del_mid_block,   
     +                          zen_del_rate_sps14, clk_del_rate_sps14

 1100            format(i3,i4,i3,i3,f5.1,f9.3,f11.3,2f12.5)

              endif

           enddo    ! stations loop

	enddo     ! blocks loop

	close (unit=7)

	return
	end
c

	subroutine make_atmos_fits_w_accel( params, ut_mid_rads )

	implicit real*8 ( a-h, o-z )

	real*8		params(156)

c       Current dimension limits are up to 12 antennas and 10 blocks !!!!!!!!!!!!
        real*8          z_delays(12,10), times(12,10)

        integer         n_per_ant(12)

        character*24    atmos_file
        
	common /constants/ c,pi,twopi,secday,sidereal_rate,
     +			deg_rad,asec_rad,frequency,wavelength

	real*8	x(12), y(12), z(12), cal_ra(300), cal_dec(300)
	common/geometry/ x,y,z,gast0,cal_ra,cal_dec,num_stations

        n_blocks_max = 10                !!!! hardwired !!!!!!!!!!!!!!!!!!!!!!!!!

        rad_to_hr = 12.d0/ pi
        
c       zero arrays
        do n_s = 1, num_stations
           n_per_ant(n_s) = 0
           do n_t = 1, n_blocks_max
              times(n_s,n_t)    = 0.
              z_delays(n_s,n_t) = 0.
           enddo
        enddo

c       Open ATMOS.FITS file with entries at geoblock times;
c       will make more frequent entries to track antenna accelerations
        atmos_file = 'ATMOS.FITS'
        open (unit=7, file=atmos_file, status='unknown')

        read (7,*) num_lines
        if ( num_lines .lt. 1 ) then
           print*,' Sub: make_atmos_fits_w_accel: ',
     +            ' STOP, no entries in ATMOS.FITS file;',
     +            ' cannot expand for clock accelerations.'
           STOP
        endif

c       Read in and store zenith delays for each antenna at each block time
        t_min = 9.d9
        t_max = 0.d0
        do n = 1, num_lines
           read (7,*) i_ant,i_day,i_hr,i_min,sec,zen_del

           n_per_ant(i_ant) = n_per_ant(i_ant) + 1
           t_hr  = i_day*24.d0 + i_hr + i_min/60.d0 + sec/3600.d0       ! hr
           times( i_ant,n_per_ant(i_ant) ) = t_hr                       ! hr
           z_delays( i_ant,n_per_ant(i_ant) ) = zen_del                 ! cm

           if ( t_hr .lt. t_min ) t_min = t_hr
           if ( t_hr .gt. t_max ) t_max = t_hr
        enddo

c       Now rewind the ATMOS.FITS file so that it can be over-written
        rewind (unit=7)

c       ------------------------------------------------------------------
c       Find largest clock acceleration value and use this to get the
c       time increment for writing new lines.   This is needed since
c       we need to pre-calculate the number of line in the new ATMOS.FITS file
        accel_max = 0.d0
        do n_s = 1, num_stations
           abs_accel = abs( params(2*num_stations + n_s) )              ! cm/hr^2
           if ( abs_accel .gt. accel_max ) accel_max = abs_accel
        enddo

        dt_used_hr = 5.d0/60.d0                                         ! hr
        if ( accel_max .gt. 0.d0 ) then
c          Assume we have a 24 hour observing block, so 12 hr max from middle;
c          Assume we don't want greater than ~0.5 cm change between entries
           dt_hr = 12.d0*( sqrt( 1.d0+1/(accel_max*12.d0) ) -1.d0 )     ! hr    
           dt_min = dt_hr * 60.d0                                       ! min
c          Round to nearest 5 min value
           i5_dt_min = 5*int( (dt_min + 2.5d0)/5.d0 )
c          Limit range to 5 to 60 minutes... 
           if ( i5_dt_min .lt. 5  ) i5_dt_min = 5
           if ( i5_dt_min .gt. 60 ) i5_dt_min = 60
           dt_used_hr = i5_dt_min/60.d0                                 ! hr
        endif

c       Establish time range for new ATMOS.FITS file extended by +/-2 hours
c       and use integer values to avoid rounding issues...
        it_min = int( t_min )
        it_max = int( t_max ) + 1
        if ( it_min .ge. 2 ) then
              it_min = it_min - 2
           else
              it_min = 0
        endif
        it_max = it_max + 2

        it_span = it_max - it_min                                       ! hr
        n_time_entries = it_span*(60.d0/i5_dt_min) + 1

c       Count number of antennas in the geoblock
        n_ants = 0
        do n_s = 1, num_stations
           if ( n_per_ant(n_s) .gt. 0 ) n_ants = n_ants + 1
        enddo
        
c       Write first line of new ATMOS.FITS...
        n_entries = n_ants * n_time_entries
        write (7,1000) n_entries
 1000   format(i10)

c       Go through station by station and generate frequent entires for all
c       stations, even if some have zero acceleration...
        do n_s = 1, num_stations

           t_entry = it_min - dt_used_hr                                ! hr

c          Check that there is at least 1 geoblock entry for this antenna...
           if ( n_per_ant(n_s) .gt. 0 ) then

              do i_t = 1, n_time_entries

                 t_entry  = t_entry + dt_used_hr                        ! hr

                 ut_hr   = t_entry
                 id      = int( ut_hr/24.d0 )
                 ut_hr   = ut_hr - id*24
                 ihr     = int( ut_hr )
                 imin    = int( (ut_hr - ihr*1.d0) * 60.d0 )
                 sec     = (ut_hr - ihr*1.0d0 - imin/60.0d0) * 3600.d0
                 
c                Linearly interpolate zenith delays from geoblock values
                 call interp_delay ( n_s, n_per_ant, t_entry, 
     +                               times, z_delays,
     +                               zen_del )

c                Calculate clock terms from full equation
c                Delta-time between entry and reference time 
                 del_hrs  = t_entry - ut_mid_rads*rad_to_hr             ! hrs

c                Clock values come from "equation"
c                Get clock parameter values ...
                 i_p = n_s
                 clk_del = params(i_p)                                  ! cm

                 i_p = i_p + num_stations
                 clk_del_rate= params(i_p)                              ! cm/hr
                 clk_del     = clk_del + clk_del_rate*del_hrs           ! cm

c                CLCOR doesn't use clock acceleration, so stating with V6 will update
c                clock and drift rate parameters to implement clock acceleration
c                in a step-wise fashion at geoblock times...
                 i_p = i_p + num_stations
                 clk_acc = params(i_p)                                  ! cm/hr^2
                 clk_del = clk_del + 0.5d0*clk_acc*del_hrs**2           ! cm
                 clk_del_rate = clk_del_rate + clk_acc*del_hrs          ! cm/hr

c                Change units on rates to AIPS CLCOR desired units...
                 clk_del_rate_sps14 = clk_del_rate * 1.d14 / (c*3600.d0) ! (sec/sec)*10^14

c                zen_del_rate_sps14 = zen_del_rate * 1.d14 / (c*3600.d0) ! (sec/sec)*10^14
                 zen_del_rate_sps14 = 0.d0 ! handled in CLCOR by interpolating delays

c                Output interpolated values
                 write (7,1100) n_s, id, ihr, imin, sec, 
     +                          zen_del, clk_del,   
     +                          zen_del_rate_sps14, clk_del_rate_sps14
 1100            format(i3,i4,i3,i3,f5.1,f9.3,f11.3,2f12.5)

              enddo
           endif

        enddo      ! station by station loop

	close (unit=7)

	return
	end
c
        subroutine interp_delay ( n_ant, n_per_ant, t_entry,
     +                            times, z_delays,  
     +                            zen_del )

c       Interpolates (or extrapolates) zenith delays from geoblock values
c       to different times for a given antenna 

	implicit real*8 ( a-h, o-z )

        real*8          z_delays(12,10), times(12,10)
        integer         n_per_ant(12)
        logical         done

        n_entries = n_per_ant(n_ant)
        if ( n_entries .le. 0 ) then
           print*,' STOP: no zenith delays in ATMOS.FITS for antenna',
     +             n_ant
           STOP
        endif

        done = .false.

c       Cases where requested time is outside geoblock entries... 
        if ( t_entry .le. times(n_ant,1) ) then
           zen_del = z_delays(n_ant,1)
           done    = .true.
        endif
        if ( t_entry .ge. times(n_ant,n_entries) ) then 
           zen_del = z_delays(n_ant,n_entries)
           done    = .true.
        endif

        if ( .not.done ) then
c          Request time is within geoblock entries...
           do n_e = 1, n_entries-1
              if ( t_entry.ge.times(n_ant,n_e)   .and.
     +             t_entry.le.times(n_ant,n_e+1) .and.
     +             .not.done                           ) then
                 n_1 = n_e
                 n_2 = n_e + 1
                 done = .true.
              endif
           enddo
c          Have geoblock entries bracketing requested for time=t_entry...
           d_time = t_entry - times(n_ant,n_1) ! hr
           t_sep  = times(n_ant,n_2) - times(n_ant,n_1) ! hr

           z_1    = z_delays(n_ant,n_1)
           z_2    = z_delays(n_ant,n_2)
           z_sep  = z_2 - z_1   ! cm
           zen_del= z_1 + z_sep*(d_time/t_sep) ! cm
        endif

        if ( .not.done ) then
           print*,'Sub: interp_delay: something wrong:',n_ant,t_entry 
           STOP
        endif


        return
        end

c
      subroutine fill_in_out ( sequence, max, flag )

      implicit real*8 (a-h,o-z)

      integer sequence(10), flag(10)

      integer first, next, last

      if ( max .ge. 1 ) then

c        set all flags good and count number of good data...
         n_good = 0
         do i = 1, max
            flag(i) = i
            if ( sequence(i) .gt. 0 ) n_good = n_good + 1 
         enddo

c        first, check for no good data at all...
         if ( n_good .eq. 0 ) then
            do i = 1, max
               flag(i) = -1
            enddo
c           All done...
            return
         endif

c        also check for all good data...
         if ( n_good .eq. max ) then
c           All done...
            return
         endif

c        ==================================================
c        OK, have some work to do...
c        check for bracketed missing points and set flag=-1
         first = 0
         next  = 0
         i     = 0
         do while ( i .lt. max )
            i = i + 1
            if ( first.eq.0 .and. sequence(i).gt.0 ) first = i
            if ( first.gt.0 .and. i.gt.first ) then
               if ( sequence(i).gt.0 ) next  = i
               if ((next-first).eq.1 ) first = next
            endif
            if ( (next-first) .gt. 1 ) then
               do j = first+1, next-1
                  flag(j) = -1
               enddo
               first = next
               next  = 0
            endif
         enddo

c        check for missing end points and fill out (flag=interior good one)
c        two passes...forward and backward
c        forward search; find last good point...
         last = 0
         do i = 1, max
            if ( sequence(i) .gt. 0 ) last = i
         enddo
         if ( last .lt. max ) then
            do i = last+1, max
               flag(i) = last
            enddo
         endif
c        backward search; find first good point...
         first = 0
         do i = max, 1, -1
            if ( sequence(i) .gt. 0 ) first = i
         enddo
         if ( first .gt. 1 ) then
            do i = 1, first-1
               flag(i) = first
            enddo
         endif

      endif

      return
      end
c
	subroutine hmsrad ( hms, hms_rad )

C	Converts hhmmss.sss format to radians

	implicit real*8 (a-h,o-z)

	common /constants/ c,pi,twopi,secday,sidereal_rate,
     +			deg_rad,asec_rad,frequency,wavelength

	xhms = abs(hms)

	ih = xhms/10000.d0
	im = dmod(xhms,10000.d0)/100.d0
	s  = dmod(xhms,100.d0)

	hms_rad = dfloat(ih) + dfloat(im)/60.d0 + s/3600.d0
	hms_rad = hms_rad*pi/12.d0
	if ( hms .lt. 0.d0 )  hms_rad = -hms_rad

	return
	end
c       ========================================================
	subroutine dmsrad ( dms, dms_rad )

c	Converts hhmmss.sss format to radians

	implicit real*8 (a-h,o-z)

	common /constants/ c,pi,twopi,secday,sidereal_rate,
     +			deg_rad,asec_rad,frequency,wavelength

	xdms = abs(dms)

	id = xdms/10000.d0
	im = dmod(xdms,10000.d0)/100.d0
	s  = dmod(xdms,100.d0)

	dms_rad = dfloat(id) + dfloat(im)/60.d0 + s/3600.d0
	dms_rad = dms_rad*deg_rad
	if ( dms .lt. 0.d0 )  dms_rad = -dms_rad

	return
	end
c       =========================================================
        subroutine radians_to_hhmmss ( ra_rad, ra_hhmmss)

c       Input :  ra_rad       in radians
c       Output:  ra_hhmmss    as hhmmss.sssss

        implicit real*8 (a-h, o-z)

	common /constants/ c,pi,twopi,secday,sidereal_rate,
     +			deg_rad,asec_rad,frequency,wavelength

        rad_to_hr = 12.d0/ pi

	ra_hr = ra_rad * rad_to_hr

        ihr  = ra_hr
        imin = (ra_hr - ihr*1.d0) * 60.d0
        sec  = (ra_hr - ihr*1.0d0 - imin/60.0d0) * 3600.d0

        ra_hhmmss  = ihr*10000.d0 + imin*100.d0 + sec

        return
        end
c       =========================================================
        subroutine radians_to_ddmmss ( dec_rad, dec_ddmmss)

c       Input :  dec_rad       in radians
c       Output:  dec_ddmmss    as hhmmss.sssss

        implicit real*8 (a-h, o-z)

	common /constants/ c,pi,twopi,secday,sidereal_rate,
     +			deg_rad,asec_rad,frequency,wavelength

        rad_to_deg = 180.d0/ pi

	if ( dec_rad .ge. 0.d0 ) then
	     dec_sign =  1.d0
	  else
	     dec_sign = -1.d0
	endif

	dec_deg = abs( dec_rad ) / deg_rad

        ideg = dec_deg
        imin = (dec_deg - ideg*1.d0) * 60.d0
        sec  = (dec_deg - ideg*1.0d0 - imin/60.0d0) * 3600.d0

        dec_ddmmss  = ideg*10000.d0 + imin*100.d0 + sec
	dec_ddmmss  = dec_ddmmss * dec_sign

        return
        end
c
	subroutine print_all_resids ( data, model, resids,
     +                     res_err, wtres, num_data, lu_print ) 

	implicit real*8 (a-h,o-z)

	real*8	data(20000), model(20000), resids(20000),
     +		res_err(20000), wtres(20000)

c	Print out Data and Rsidual...
	lu_print = 6
	write (lu_print,3000) 
 3000	format (//'         Data       Model    Residual      ',
     +                                ' Error     Res/Err'/)
	do i = 1, num_data
	   wtres(i) = resids(i) / res_err(i)
	   write (lu_print,3100) data(i), model(i), resids(i), 
     +				res_err(i), wtres(i)
 3100	   format (1x,5f12.5)
	enddo		

	return
	end
c
	subroutine geo_plot_files ( num_data, 
     +                          list_stations, n_stats_used, 
     +                          n_stat_a, n_stat_b, gst, 
     +                          data, model, resids, res_err )

	implicit real*8 ( a-h, o-z )

	real*8	data(20000), model(20000), resids(20000),
     +		res_err(20000), gst(20000)

	integer	n_stat_a(20000), n_stat_b(20000)

        integer list_stations(12)

	character*32  outfile


	common /constants/ c,pi,twopi,secday,sidereal_rate,
     +			deg_rad,asec_rad,frequency,wavelength

	real*8	x(12), y(12), z(12), cal_ra(300), cal_dec(300)
	common/geometry/ x,y,z,gast0,cal_ra,cal_dec,num_stations

	sid_hr = 12.d0/pi/sidereal_rate   ! convert sid(rad) -> hr
	cm_ns  = 1.d09/c                  ! convert cm -> nsec

c	Make ascii files of data, model, resids (1 file per baseline 
c       per source) for geodetic-like data (multi-band delays and rates).

c       -------------------------------------------------------------
c       Open a file for every possible baseline (over-writing and erasing 
c       any old files if necessary).  This way a missing baseline still 
c       has a file, which makes plotting scripts happier.

cv1	n1 = n_stats_used - 1 
	n1 = num_stations - 1 
	do n_a = 1, n1

	   n2 = n_a + 1
cv1	   do n_b = n2, n_stats_used
	   do n_b = n2, num_stations

	      write (outfile,4000)  n_a, n_b
 4000	      format ('fit_geodetic_bl',2i2.2,'.dat')
	      open (unit=7, status='unknown', file=outfile,
     +		       form='formatted')
	      write (7,4100)  n_a, n_b
 4100	      format( '!  Baseline ',i2,'-',i2,/'!     UT(hrs)   ',
     +              '     Delay (ns) D-Model  D-Resid ',
     +               '    Rate (mHz) R-Model  R-Resid')
	      close (unit=7)

	   enddo

	enddo

c       -------------------------------------------------------------
c	Now check if there are any geodetic data 
	if ( num_data .ge. 1 ) then

c          Cycle through all baselines and collect geodetic data for
c          each baseline....
	   n1 = n_stats_used - 1 
	   do n_a = 1, n1

	      n2 = n_a + 1
	      do n_b = n2, n_stats_used

c                Get actual VLBA antenna number from parameter number
                 num_ant_a = list_stations( n_a )
                 num_ant_b = list_stations( n_b )

c                Re-open the file for baseline "num_ant_a - num_ant_b"
		 write (outfile,4000)  num_ant_a, num_ant_b
		 open (unit=7, status='old', file=outfile,
     +		       form='formatted', access='append')

c                Now search for geo data for this baseline 
c                by sifting through all geodetic data pairs ...
		 do n_d = 1, num_data, 2

		    if ( n_stat_a(n_d) .eq. num_ant_a  .and.
     +                   n_stat_b(n_d) .eq. num_ant_b      ) then

		       ut_hrs  = (gst(n_d) - gast0)*sid_hr
		       n_r = n_d + 1                       ! rate data index

c                      Convert to more standard units (simplify GRAPHIC file)
		       del_d_ns = data(n_d)  * cm_ns
		       del_m_ns = model(n_d) * cm_ns
		       del_r_ns = resids(n_d)* cm_ns
		       rate_d_mhz = data(n_r)  * 1000.d0
		       rate_m_mhz = model(n_r) * 1000.d0
		       rate_r_mhz = resids(n_r)* 1000.d0
		       write (7,4200) ut_hrs, 
     +		            del_d_ns, del_m_ns, del_r_ns,
     +                      rate_d_mhz, rate_m_mhz, rate_r_mhz
 4200		       format(1x,f12.4,4x,3f10.4,3x,3f10.2)

		    endif

		enddo

		close ( unit=7 )

	     enddo              ! b-station baseline forming loop

	  enddo             ! a-station baseline forming loop

	endif		! making sure have some geodetic data

	return
	end
c
C*JULDA -- find Julian Date and GMST at midnight for a given day
C+
      SUBROUTINE JULDA (YEAR, MONTH, DAY, DATE, gmst_hr)
      INTEGER YEAR
      INTEGER MONTH
      INTEGER DAY
      real*8  DATE
      real*8  GMST
C
C Returns DATE = Julian Date and GMST = Greenwich Mean sidereal time
C at 0 hrs U.T. on day DAY/MONTH/YEAR. Accuracy not tested; GMST is
C correct at least to nearest second
C
C History:
C  1977 Aug  5 - TJP.
C  1991 May 18 - TJP.
C-----------------------------------------------------------------------
      INTEGER MOFF(12), IC, NYRM1, JD, NDAYS
      DOUBLE PRECISION DU,SMD,T, gmst_hr
      DATA MOFF/0,31,59,90,120,151,181,212,243,273,304,334/
C
C JD number at 12 hrs UT on Jan 0 of year YEAR (Gregorian Calendar).
C
      NYRM1 = YEAR-1
      IC = NYRM1/100
      JD = 1721425 + 365*NYRM1 + NYRM1/4 - IC + IC/4
C
C Number of days from Standard Epoch 1900 Jan 0.5 (JD 2415020.0) to
C Jan 0.0 of YEAR.
C
      DU = DBLE(JD-2415020) - 0.5D0
C
C Day number; is it a leap year?
C
      NDAYS = MOFF(MONTH) + DAY
      IF (MONTH.GT.2 .AND. ( 
     1     ( MOD(YEAR,4).EQ.0 .AND. MOD(YEAR,100).NE.0 ) .OR.
     2     MOD(YEAR,400).EQ.0 ) ) NDAYS = NDAYS+1
C
C UT from Epoch to today (days), (centuries).
C
      SMD = DU+DBLE(NDAYS)
      T = SMD/36525.D0
C
C Greenwich Mean Sidereal Time.
C
      GMST = DBLE(6*3600 + 38*60) +45.836D0
     1     + 8 640 184.542D0*T  +  0.0929D0*T**2
      GMST = MOD(GMST,86400D0)

      gmst_hr = GMST / 3600.d0
C
C Julian Date.
C
      DATE = 2415020.D0+SMD
C
      END
c
      SUBROUTINE LEAST_SQUARES_FIT ( IDEBUG, print_solution,
     +                          PRINT_COR, MAXDAT,
     +                          NUMXES, N, err_factor,
     +				VNAME, X, ID, S,
     +                          C, E, XHAT2, ERROR )

C     SUBROUTINE FOR LINEAR LEAST SQUARE FITTING...
C     IT REQUIRES INPUT OF THE FOLLOWING:
C        IDEBUG    = DEBUGGING PARAMETER
C                    0 = SUPRESS PRINT OUT BEYOND SOLUTION AND CORRELATI
C                    3 = PRINT OUT ALL DATA AND MATRICES
c        print_solution = logical flag: print out least-squares solution
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

      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION C(MAXDAT,1)
      DIMENSION S(1),E(1),X(1),ID(1),XHAT2(1)

      character*8  VNAME(1),QUEST(2)

      REAL*8       B(156,156),BSTORE(156,156),XIDENT(156,156),
     +             COR(156,156),ERROR(156),STDEV(156),PR(156),
     +             XHAT(156),SNEW(20000),PRODCT(20000)

      DIMENSION LL(156),MM(156)
      DATA      MSIZE /156/


      LOGICAL PRINT_COR, print_solution

      DATA QUEST/'YES     ','NO      '/

 9999 FORMAT (1X,10(1PD13.5))
 9998 FORMAT (1X,20F6.2)
 9417 FORMAT (/1X)

      IF (IDEBUG.ge.3) then

	 WRITE (6,5522)
 5522	 FORMAT (////50X,' C MATRIX')
	 do I=1,N
	    WRITE (6,9999) (C(I,J),J=1,NUMXES)
	 enddo
	 WRITE (6,5523)
 5523	 FORMAT (////40X,'        S(I)           E(I)')
	 do I=1,N
	    WRITE (6,9999) S(I),E(I)
	 enddo

      endif

C     Count number of solve-for parameters (M)...
      M=0
      DO I=1,NUMXES
	 M=M+ID(I)
      enddo

c     Check number of data .ge. number of parameters...
      IF ( N.lt.M )  then
	 WRITE (6,2094)  N,M
 2094	 FORMAT (////' LEAST_SQUARES_FIT: # DATA POINTS (',I4,
     +                        ') < # PARAMETERS (',I2,').  STOP' )
	 STOP
      endif

c     Compress partial derivative matrix to retain only those
c     associated with solve-for parameters...
      JNEW=0
      do J=1,NUMXES
	  if ( ID(J) .ne. 0 ) then
	     JNEW=JNEW+1
	     do I=1,N
   		 C(I,JNEW)=C(I,J)
	     enddo
	  endif
      enddo

C     Weight equations by dividing each by their error: E(I)
      do I=1,N
	  SNEW(I)=S(I)/E(I)
	  do J=1,M
	     C(I,J)=C(I,J)/E(I)
	  enddo
      enddo

C     Debug printout...
      IF (IDEBUG.ge.3) then
	 WRITE (6,2006)
 2006	 FORMAT ('1'/////50X,' CSTAR MATRIX')
	 do I=1,N
 	    WRITE (6,9999) (C(I,J),J=1,M)
	 enddo
	 WRITE (6,2009)
 2009	 FORMAT (/////50X,'     SSTAR')
	 do I=1,N
	    WRITE (6,9999) SNEW(I)
	 enddo
      endif

C     Square up partial deravitive matrix...
      ITEST2=0
      BMIN10=0.D0
      BMAX10=0.D0
      DO I=1,M
	DO J=1,M
	   B(I,J)=0.D0
	enddo
      enddo

      do I=1,M
	 do J=1,M
	      do L=1,N
		 B(I,J)=B(I,J) + C(L,I)*C(L,J)
              enddo
c             Test for big range of exponents...
	      if (B(I,J).ne.0.D0) then
		 B1=DLOG10( DABS( B(I,J) ) )
		 IF (BMIN10.GT.B1) BMIN10=B1
		 IF (BMAX10.LT.B1) BMAX10=B1
              endif
 	      BSTORE(I,J)=B(I,J)
	  enddo
      enddo

c     More debug printout...
      if (IDEBUG.ge.3) then
	 WRITE (6,2010)
 2010	 FORMAT ('1'/////50X,' C*TRANSPOSE C')
	 do I=1,M
	    WRITE (6,9999) (B(I,J),J=1,M)
	 enddo
      endif


c     Warn user if the matrix has a large range of exponents...
      IF (DABS(BMAX10-BMIN10).ge.16.D0) then
	 WRITE (6,9622) BMIN10,BMAX10
 9622	 FORMAT (///1X,'********   BMIN10 = ',F6.1,
     +          'BMAX10 = ',F6.1,'   *************',///1X)
      endif

c     "Center" exponents about zero...
      BADJST=10.D0**( (BMAX10+BMIN10)/2.D0 )
      do I=1,M
	 do J=1,M
 	    B(I,J)=B(I,J)/BADJST
	 enddo
      enddo

c     About to invert the matrix...
C        THE SUBROUTINE 'MINV' IS A DOUBLE PRECISION MATRIX INVERSION
C        IT INVERTS THE MATRIX 'BB' AND STORES IT AS 'BB' (I.E. THE ORIGIN
C        MATIRX IS DESTROYED

      CALL MINV(B,M,DETERM,LL,MM,MSIZE)

c     Re-scale inverted matrix for "Centering"
      do I=1,M
	 do J=1,M
	    B(I,J)=B(I,J)/BADJST
	 enddo
      enddo

c     More debug printout...
      IF (IDEBUG.ge.3) then
	 WRITE (6,2011) DETERM
 2011	 FORMAT ('1'////'  THE DETERMINANT IS',1PD13.5)
	 WRITE (6,2022)
 2022	 FORMAT (////45X,' (C*TRANSPOSE C*) INVERSE')
	 do I=1,M
 	    WRITE (6,9999) (B(I,J),J=1,M)
	 enddo
      endif

c     Prepare to check matrix inversion by multiplying inverse matrix
c     by the original martix to obtain (hopefully) the indentity matrix...
      do I=1,M
	 do J=1,M
	    XIDENT(I,J)=0.D0
	 enddo
      enddo

      do I=1,M
	 do J=1,M
	    do L=1,M
	       XIDENT(I,J)=XIDENT(I,J) + B(I,L)*BSTORE(L,J)
	    enddo
	    if (I.ne.J) then
c               Off-diagonal elements should be near zero...
		IF (DABS(XIDENT(I,J)).GT.1.D-06) ITEST2=1
	      else
c               Diagonal elements should be near unity...
		IF (DABS(XIDENT(I,I)-1.D0).GT.1.D-06) ITEST2=-1
	    endif
	 enddo
      enddo

      if (ITEST2.eq.0) then

c        Matrix seemed to invert.
c        Calculate corrections (XHAT's) to parameters...
	 do I=1,M
	    XHAT(I)=0.D0
            do J=1,N
	       PRODCT(J)=0.D0
	       do L=1,M
   		  PRODCT(J)=PRODCT(J) + B(I,L)*C(J,L)
	       enddo
	       XHAT(I)=XHAT(I)+PRODCT(J)*SNEW(J)
	    enddo
         enddo

C        XHAT'S ARE (IF ADDED TO THE X'S) THE UNCONSTRAINED LEAST SQUARE
C        SOLUTION OF THE 'SOLVE FOR' PARAMETERS ONLY.
C        XHAT2(J) = THE LEAST SQUARES SOLTIONS (INCLUDING NON SOLVE FOR
C        PARAMETERS)
	 IN=0
	 do J=1,NUMXES
	    if (ID(J).ne.0) then
	          IN=IN+1
	          XHAT2(J)=XHAT(IN) + X(J)
	       else
    		  XHAT2(J)=X(J)
	    endif
	 enddo

c        Calculate correlation coefficients...
	 do I=1,M
	    do J=1,M
	       COR(I,J)=-BSTORE(I,J)/DSQRT(BSTORE(I,I)*BSTORE(J,J))
	    enddo
	 enddo

c       Check for unphysical negative variances...
	INOTE=0
	do I=1,M
	      IF (B(I,I).le.0.D0) then
		 B(I,I)=-B(I,I)
		 INOTE = INOTE + 1
	      endif
	enddo
	if (INOTE.gt.0) then
	   WRITE (6,2071) INOTE
 2071	   FORMAT (///' ***** THERE ARE ',I2,' NEGATIVE VARIANCES; ',
     +            'SOLUTION IS UNPHYSICAL')
	endif

c       Convert variances to standard deviations...
	do J=1,M
	   STDEV(J)=DSQRT(B(J,J))
	enddo


C       REARRANGE CALCULATED 1 SIGMA ERRORS TO APPLY TO THE SOLVED FOR
C       PARAMETERS
	IN=0
	do J=1,NUMXES
	   if (ID(J).ne.0) then
	         IN=IN+1
		 ERROR(J)=STDEV(IN)
	      else
		 ERROR(J)=0.D0
	   endif
	enddo

        if ( print_solution ) then
C          OUTPUT FOR THE UNCONSTRAINED LEAST SQUARES SOLUTION
           WRITE (6,2040)
 2040	   FORMAT (/1X,'PARAMETER    ORIGINAL VALUE   LEAST SQUARE',
     +               ' VALUE  1 SIGMA ERRORS   SOLVED FOR?')

	   do J=1,NUMXES
	      L=1
	      IF (ID(J).EQ.0) L=2
c             multiply formal errors by "err_factor" since for this problem
c             the baseline data are correlated...
              ERROR(J) = ERROR(J) * err_factor
 	      WRITE (6,2041) VNAME(J), X(J), XHAT2(J), ERROR(J), QUEST(L)
 2041	      FORMAT (2X,A8,5X,3(F13.6,5X),4X,A8)
	   enddo
        endif

c       Print correlations?
	if ( PRINT_COR )  then
C	      CALCULATE THE CORRELATION COEFFICIENTS (COR)... 
	      do I=1,M
	         do J=1,M
		    COR(I,J)=B(I,J)/DSQRT(B(I,I)*B(J,J))
		 enddo
	      enddo

	      WRITE (6,2056)
 2056	      FORMAT(/10X,' THE CORRELATION COEFFICIENTS')
	      do I=1,M
   		  WRITE (6,9998) (COR(I,J),J=1,M)
	      enddo

C  	      THE MULTIPLE CORRELATION COEFFICIENTS (PR) ARE CALCULATED
	      do I=1,M
		 PR(I)= 1.D0 - (1.d0/(BSTORE(I,I)*B(I,I)))
	      enddo
	      WRITE (6,2060)
 2060         FORMAT (///10X,'THE MULTIPLE CORRELATION COEFFICIENTS'/)
              I = 0
              do J=1,NUMXES
              	if ( ID(J).ne.0 ) then
              		I = I + 1
C                       Not sure if the true PR is the sqrt(PR)??
                        WRITE (6,2061) VNAME(J),PR(I)
 2061                   FORMAT (10X,A8,2X,F10.5)
		endif
	      enddo

         endif

       else

c         Print out warning...
	  WRITE (6,2095)  ITEST2
 2095	  FORMAT (////'  ***** ITEST2 =',I2,' ***** '
     +                   /'       MATRIX INVERSION FAILED.')
	  STOP

      endif

      RETURN
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
