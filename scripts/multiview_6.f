      program multiview

c----------------------------
c     Combines phases of for N>=2 QSOs (Target-QSO position data) to estimate
c     the expected phase offset at the origin (Xo,Yo)=(Target position),
c     as well as the slope parameters.  For inverse-MultiView the slope
c     parameters are not used, but may be good diagnostics

c     This version started with "planar_offsets.f" 

c     It solves for 3 parameters for phase via:
c           model = slope_x*offset_x + slope_y*offset_y + intercept
c     Slopes are in degrees_of_phase per degrees_of_offset
c     Intercept is in degrees_of_phase      

c     Input files:
c        multiview_control.inp   ! control file with multiview information
c        target.TBOUT            ! input file from TBOUT of SN-table for target
c        multiview.TBOUT         ! input  "          "        "       "  QSOs

c     Output file:
c        target.TBIN             ! output file in TBIN form 
  
c     MJR: started 2020 Sept 2;
c     V4a: 2021 Aug 14: scaled fitted phases by sigma_pdf (if >=4 data)
c     V4b: 2021 Aug 16: added "editing" of output phases when fitting
c          error exceeds, say, 30deg      
c     V4x  2021 Aug 21: made GREG friendly outputs for diagnostic plotting
c     V5:  Dealing with phase-wraps
c     V6:  Adding 2pt MultiView capabilities
      
c----------------------------

      implicit real*8 (a-h,o-z)

      character*48  target_infile, multiview_infile, target_outfile
      character*48  control_file

      real*8  target_times(3000), target_ants(3000)
      real*8  target_reals(3000), target_imags(3000)

      real*8  QSO_srcs(8000), QSO_times(8000), QSO_ants(8000)
      real*8  QSO_reals(8000),QSO_imags(8000)

      real*8  QSO_offset_x(8), QSO_offset_y(8)
      real*8  QSO_used_x(8), QSO_used_y(8)

      integer source_number(8), n_hits_per_QSO(8), MV_src(80)

      real*8  MV_real(80), MV_imag(80)

      real*8  MV_amp(8), MV_phase(8), MV_phase_error(8)
      real*8  MV_phase_used(8), MV_error_used(8)

      real*8  time_store(3000), phase_store(3000) 
      integer iant_store(3000)

      real*8   early_times(200), late_times(200)
      integer  wrap_ants(200), wrap_QSOs(200), wrap_turns(200)

      
c     Set some maximum control values, switches, and logical unit numbers...
      max_num_data = 3000                 ! integrations
      max_QSOs     = 8                    ! multiview QSOs
      max_wraps    = 200                  ! phase wrap times/ants

c     Lower limit for amplitude loss when vector averaging
c     complex visibilities for a QSO.  If below limit, don't use it.
      amp_limit = 0.25d0

c     Upper limit for inverse-MV phase error
      phase_err_limit = 30.d0             ! degrees

c     Arbitrarily set QSO phase errors to 10 deg, so all QSOs get
c     roughly equal weight when fitting the phase offset plane...
      do i = 1, max_QSOs
         MV_phase_error(i) = 10.d0        ! degrees
      enddo

      lu_control   = 9
      lu_print     = 6
      lu_out       = 7
      lu_data      = 8
      lu_summary   =10

      pi     = 4.d0*atan(1.d0)
      degrad = pi/180.d0

c     =============================================================
c     Get input control parameters...
      control_file = 'multiview_control.inp'
      call control_parameters ( lu_control, control_file, 
     +                          lu_print, max_QSOs, max_wraps,
     +                          n_Ifs, time_window, 
     +                          n_QSOs, source_number, 
     +                          QSO_offset_x, QSO_offset_y,
     +                          n_wraps, wrap_ants, wrap_QSOs,
     +                          wrap_turns, early_times, late_times )


      write (6,1000) n_QSOs, n_IFs, time_window
 1000 format(/' Assumes SN TBOUT file has',i2,' QSOs',
     +        ' and ',i2,' IFs',
     +       /' Multiview will use time window of +/-',f5.1,' min')

      time_window_days = time_window/ (60.d0 * 24.d0)

c     Read target data file into memory...
      target_infile  = 'target.TBOUT'
      call read_target_file ( lu_control, lu_print, lu_out,  
     +                 target_infile, target_outfile,
     +                 max_num_data, n_IFs,
     +                 max_ant_num, npts_target, 
     +                 target_times, target_ants,
     +                 target_reals, target_imags )

c     Read MultiView data file into memory...
      multiview_infile = 'multiview.TBOUT'
      call read_multiview_file ( lu_control, lu_print, 
     +                           multiview_infile, max_num_data, n_IFs,
     +                           n_QSOs, source_number,
     +                           max_ant_num, npts_QSO, 
     +                           QSO_srcs, QSO_times, 
     +                           QSO_ants, QSO_reals, QSO_imags )
      
c     =============================================================
c     At each antenna and target time, collect multiview phases  
c     within the input time window, and do multiview fitting...

c     First, document output files...
      write (lu_summary,1100)
 1100 format('!Ant UT(days)  phase(deg) +/-      ',
     +     'X-slope (deg/deg) +/-   Y-slope (deg/deg) +/-',
     +     '   for MultiView fits')

      n_store = 0
      do k = 1, max_ant_num     ! antennas

         iant   = k
         lu_ant = 10 + k
         lu_ant_phases = 100 + k
         lu_for_plots = 200 + iant
         
         write (lu_print,2801) k
 2801    format(///'! Processing antenna #',i2,' ********************')
         write (lu_ant,2800) k
         write (lu_ant_phases,2800) k
 2800    format('! Processing antenna #',i2)

         write (lu_ant,1200)
 1200    format('!Ant        UT(days)    Amp   phs        Amp   phs   ',
     +          '     Amp   phs        Amp   phs    for QSOs')
         write (lu_ant_phases,1300)
 1300    format('!Ant        UT(days)  phs      phs',
     +          '      phs      phs    for QSOs')
         write (lu_for_plots,1400)
 1400    format('!Ant UT(days)    phs        X-slope (deg/deg)',
     +          ' Y-slope   from MultiView fits')
         
         
         do i = 1, npts_target  !  times

            if ( target_ants(i) .eq. k ) then

c              get target time for center of MV window
               t_center = target_times(i)

               n = 0
               do j = 1, npts_QSO                     ! MV QSOs

                  dt = abs( QSO_times(j) - t_center ) ! days
                  if ( dt.le.time_window_days .and.
     +                 QSO_ants(j).eq.k             ) then

                     n = n + 1

                     MV_src(n)  = QSO_srcs(j)
                     MV_real(n) = QSO_reals(j)
                     MV_imag(n) = QSO_imags(j)

                  endif         ! time window check

               enddo            ! MV QSOs data

               n_MV = n

c              Check that we have at least 2 MV QSOs...
c              If >=3 will fit a plane; if 2 a straight line
               call count_data( n_MV, n_QSOs, MV_src, 
     +                          n_hits_per_QSO, num_QSOs ) 
               if ( num_QSOs .ge. 2 ) then

                  write (lu_print,2900) t_center, n_MV, num_QSOs
 2900             format(' Target time',f9.5,' has a total of',i3,
     +                      ' MV data involving',i2,' QSOs.')

                  call average_by_QSO ( n_MV, MV_src, 
     +                                     MV_real, MV_imag,
     +                                     n_QSOs,
     +                                     MV_amp, MV_phase )

c                 Fix up phase wraps..
                  call fix_wraps ( n_wraps, wrap_ants,
     +                                wrap_QSOs, wrap_turns,
     +                                early_times, late_times,
     +                                iant, t_center, MV_phase )
                     
                  write (lu_ant,2950) k, t_center, 
     +                 (MV_amp(iQ), MV_phase(iQ), iQ=1,n_QSOs)
 2950             format(i3,7x,f10.5,8(f8.3,f6.0,3x))
                  write (lu_ant_phases,2952) k, t_center, 
     +                     (MV_phase(iQ), iQ=1,n_QSOs)
 2952             format(i3,7x,f10.5,8(f6.0,3x))

c                 --------------------------------------------
c                 Now fit plane to phases vs (X,Y) QSO-offsets...

c                 First, toss bad/missing data...
                  call compress( n_QSOs, amp_limit,
     +                     MV_amp, MV_phase, MV_phase_error, 
     +                     QSO_offset_x, QSO_offset_y,
     +                     n_good_QSOs, 
     +                     MV_phase_used, MV_error_used,
     +                     QSO_used_x, QSO_used_y )

c                 "compress" checks for low amplitudes after 
c                 vector averaging, so can lose a QSO. 
c                 Need to check (again) for at least 2 good ones...

c                 -------------------------------------------- 
c                 If >= 3, fit plane to phases vs (X,Y) QSO-offsets...
                  if ( n_good_QSOs .ge. 3 ) then

c                    Fit phases to a plane
                     call fit_plane(n_good_QSOs, 
     +                        MV_phase_used, MV_phase_error, 
     +                        QSO_used_x, QSO_used_y,
     +                        plane_slope_x, plane_slope_y, phs_mv,
     +                        err_slope_x, err_slope_y, phs_mv_err)

                     if ( phs_mv_err .gt. phase_err_limit ) then
                        write (lu_print,2955) iant, t_center,
     +                                              phs_mv_err
 2955                   format(' Rejected MV fit: ant#',i2,
     +                         ' at time',f8.5,': phase error=',f7.1)
                     endif
                   
                  endif
                  
c                 --------------------------------------------- 
c                 If =2, fit straight line 
                  if ( num_QSOs .eq. 2 ) then

                     x1 = QSO_used_x(1)
                     y1 = QSO_used_y(1)
                     
                     x2 = QSO_used_x(2)
                     y2 = QSO_used_y(2)

                     phs1 = MV_phase_used(1)
                     phs2 = MV_phase_used(2)
                     
                     call two_point_phase ( x1,y1,phs1, x2,y2,phs2,
     +                    phs_mv, plane_slope_x, plane_slope_y )

c                    Enter "dummy" values for some parameter errors...
                     phs_mv_err = 10.d0
                     err_slope_x = 1.d0
                     err_slope_y = 1.d0

                  endif

c                 ---------------------------------------------
c                 If <2, print warning                 
                  if ( num_QSOs .lt. 2 ) then
                     write (lu_print,3000) num_QSOs, t_center
 3000                format(' Not enough QSOs (',
     +                    i1,') at UT',f10.5,' to fit MV plane.')
c                    Flag for no outputs...
                     phs_mv_err = 999.d0
                  endif

               endif
                  
               if ( phs_mv_err .le. phase_err_limit ) then
c                 Store for later output
                  n_store = n_store + 1
                  time_store(n_store) = t_center
                  iant_store(n_store) = iant
                  phase_store(n_store)= phs_mv
                           
                  write (lu_summary,2957) iant, t_center,
     +                 phs_mv, phs_mv_err, 
     +                 plane_slope_x, err_slope_x,
     +                 plane_slope_y, err_slope_y
 2957             format(i3,f10.5,2f8.0,5x,2f10.3,5x,2f10.3)

                  write (lu_for_plots,2958) iant, t_center,
     +                 phs_mv, plane_slope_x, plane_slope_y
 2958             format(i3,f10.5,f8.0,5x,f10.3,5x,f10.3)

                  write (lu_print,2960)  
     +                 phs_mv,   plane_slope_x, plane_slope_y, 
     +                 phs_mv_err, err_slope_x,   err_slope_y
 2960             format(' MV fit: ',f10.0, 2f10.3,
     +                             /'      +/-',f10.0, 2f10.3)
               endif

            endif   ! target antenna number
               
         enddo   ! target times

      enddo   ! antennas

c     Output to target.TBIN file...
      target_outfile = 'target.TBIN'
      open ( file=target_outfile, unit=lu_out )

      call write_output( lu_print, lu_data, lu_out, target_infile,
     +                   n_store, time_store, iant_store, phase_store)



      end

c===================================================================
      subroutine read_target_file ( lu_data, lu_print, lu_out, 
     +                              data_file, out_file,
     +                              max_num_data, n_IFs,
     +                              max_ant_num, n_pts, 
     +                              times, antennas,
     +                              reals, imags )
  
      implicit real*8 (a-h,o-z)

      character*48  data_file, out_file

      character*512 c384
      character*8   c8
      equivalence (c384,c8)
      
      
      real*8  times(3000), antennas(3000), reals(3000), imags(3000)

      pi = 4.d0*atan(1.d0)
      rad_deg = 180.d0/pi

c     Open input SN file created with TBOUT
      call open_ascii_file ( lu_data, data_file, lu_print )

c     Skip the introductory matierial...
      ieof = 0
      c8   = '        '
      do while ( ieof.ge.0 .and. c8.ne.'***BEGIN' )
         read (lu_data,1000,iostat=ieof) c384
 1000    format(a384)
      enddo

c     Ready to read SN table data...
      print*,' Rec#  UT(hrs)  Ant#  Src#    Phase(deg)'

      n    = 0
      nrec_old = 0
      max_ant_num = 0
      do while ( ieof.ge.0 .and. c8.ne.'***END*P') 

c        Test if new record...if just next IF, skip
         read (lu_data,*,iostat=ieof) nrec
         if ( ieof.ge.0 .and. nrec.ne.nrec_old ) then

            nrec_old = nrec
            backspace ( unit=lu_data )

            read (lu_data,*) nrec, t, dt, nsrc, nant, nsub,
     +            nfreq, fr, node, delmb, disp, ddisp, vreal, vimag,
     +            delay, rate, wt, irefant

            if ( n.lt.max_num_data ) then 
            
               if ( vreal.gt.-999.0 .and. vimag.gt.-999.0 ) then

                  n = n + 1

                  times(n) = t                    ! hours
                  antennas(n) = nant
                  reals(n) = vreal
                  imags(n) = vimag

                  phase = atan2(vimag,vreal)*rad_deg
                  write (6,2000) nrec,t,nant,nsrc,phase
 2000             format(i5,f10.6,i5,i5,5x,f7.1)

                  if ( nant .gt. max_ant_num ) max_ant_num = nant

               endif

            endif

         endif

      enddo

      close ( unit=lu_data )

      n_pts = n
      print*,' Number of target records accepted: ',n_pts
   
      if ( n_pts .ge. max_num_data ) then
         print*,' STOP: hit dimension limit (max_num_data)'
      endif

      return
      end

c===================================================================
      subroutine write_output( lu_print, lu_data, lu_out, data_file,
     +           n_store, time_store, iant_store, phase_store )
  
      implicit real*8 (a-h,o-z)

      character*48  data_file

      real*8  time_store(3000), phase_store(3000) 
      integer iant_store(3000)

      character*512 c384
      character*8   c8
      equivalence (c384,c8)

      logical found
      
      pi = 4.d0*atan(1.d0)
      deg_rad = pi/180.d0

c     Maximum allowed time difference between this multiview phase
c     and that found in the input file
      t_diff_allowed = 30.d0/86400.d0         ! days

c     Re-open input SN file created with TBOUT
      call open_ascii_file ( lu_data, data_file, lu_print )

c     Copy the introductory matierial...
      ieof = 0
      c8   = '        '
      do while ( ieof.ge.0 .and. c8.ne.'***BEGIN' )
         read (lu_data,1000,iostat=ieof) c384
         write(lu_out,1000) c384
 1000    format(a384)
      enddo

c     Read SN entry, copy most of it and update phases to output file...
c     Also, sequentially number the output lines
      n_out =0
      do while ( ieof.ge.0 .and. c8.ne.'***END*P') 

            read (lu_data,*,iostat=ieof) nrec, t, dt, nsrc, nant, nsub,
     +            nfreq, fr, node, delmb, disp, ddisp, vreal, vimag,
     +            delay, rate, wt, irefant

            if ( ieof.ge.0 .and. c8.ne.'***END*P' ) then 

               found = .false.
               n = 0
               do while ( .not.found .and. n.lt.n_store )

                  n = n + 1

                  t_diff = abs( t - time_store(n)  )        ! days    
                  if ( t_diff.le.t_diff_allowed .and. 
     +                 nant.eq.iant_store(n) ) found = .true.

               enddo

               if ( .not.found ) then 
                  write (lu_print,1500) t, nant 
 1500             format('No multiview phase at ',f10.5,'for ant ',i2)
               endif

               if ( found ) then

c                 Have matching time and antenna, and so
c                 replace reals and imags with multiview version.
c                 NB: assumes both polarizations identical
                  phsrad = phase_store(n) * deg_rad
                  vreal  = cos( phsrad ) 
                  vimag  = sin( phsrad )

c                 zero delays and rates
                  delay = 0.d0
                  rate  = 0.d0

c                 Write output record...
                  n_out = n_out + 1
                  write (lu_out,2000) n_out, t, dt, nsrc, nant, nsub,
     +               nfreq, fr, node, 
     +               delmb, disp, ddisp, vreal, vimag,
     +               delay, rate, wt, irefant,
     +               delmb, disp, ddisp, vreal, vimag,
     +               delay, rate, wt, irefant
 2000             format(i8,e24.15,e15.6,4i11,e15.6,i11,8e15.6,i11,
     +                   8e15.6,i11)

               endif

         endif

      enddo

      write (lu_out,3000) 
 3000 format('***END*PASS***')

      close ( unit=lu_data )
      close ( unit=lu_out  )

      return
      end

c===================================================================
      subroutine read_multiview_file ( lu_data, lu_print, 
     +           data_file, max_num_data, n_IFs,
     +           n_QSOs, source_number, 
     +           max_ant_num, n_pts, 
     +           QSO_srcs, QSO_times, 
     +           QSO_ants, QSO_reals, QSO_imags )
  
      implicit real*8 (a-h,o-z)

      character*48  data_file
      character*8   c8
      
      real*8  QSO_srcs(8000), QSO_times(8000), QSO_ants(8000)
      real*8  QSO_reals(8000), QSO_imags(8000)

      integer source_number(8)

      pi = 4.d0*atan(1.d0)
      rad_deg = 180.d0/pi

      call open_ascii_file ( lu_data, data_file, lu_print )

c     Skip the introductory matierial...
      ieof = 0
      c8   = '        '
      do while ( ieof.ge.0 .and. c8.ne.'***BEGIN' )
         read (lu_data,*,iostat=ieof) c8
      enddo

c     Ready to read SN table data...
      print*,' Rec#  UT(hrs)  Ant# Src# MV#    Phase(deg)'

      n    = 0
      nrec_old = 0
      do while ( ieof.ge.0 .and. c8.ne.'***END*P') 

c        Test if new record...if just next IF, skip
         read (lu_data,*,iostat=ieof) nrec
         if ( ieof.ge.0 .and. nrec.ne.nrec_old ) then

            nrec_old = nrec
            backspace ( unit=lu_data )

            read (lu_data,*) nrec, t, dt, nsrc, nant, nsub,
     +            nfreq, fr, node, delmb, disp, ddisp, vreal, vimag,
     +            delay, rate, wt, irefant

            if ( n.lt.max_num_data ) then 
            
               do j = 1, n_QSOs

                  if ( nsrc .eq. source_number(j) ) then

                     if ( vreal.gt.-999.0 .and. vimag.gt.-999.0 ) then
                     
                        n = n + 1

                        QSO_srcs(n)  = j         ! MV QSO index
                        QSO_ants(n)  = nant      ! antenna
                        QSO_times(n) = t         ! time (hours)
                        QSO_reals(n) = vreal     ! real
                        QSO_imags(n) = vimag     ! imag

                        phase = atan2(vimag,vreal)*rad_deg
                        write (6,1000) nrec,t,nant,nsrc,j,phase
 1000                   format(i5,f10.6,i5,i5,'-->',i1,5x,f7.1)

                     endif

                  endif

               enddo

c              Skip other IFs...
               do m =  1, n_IFs-1
                  read(lu_data,*) nrec
               enddo

            endif

         endif

      enddo

      close ( unit=lu_data )

      n_pts = n
      print*,' Number of multiview records accepted: ', n_pts
   
      if ( n_pts .ge. max_num_data ) then
         print*,' STOP: hit dimension limit (max_num_data)'
      endif

      return
      end

c==========================================================
      subroutine count_data( n_MV, n_QSOs, MV_src, 
     +                       n_hits_per_QSO, num_QSOs ) 

      implicit real*8 (a-h,o-z)

      integer  MV_src(80), n_hits_per_QSO(8)

c     Zero counter...
      do n = 1, n_QSOs
         n_hits_per_QSO(n) = 0
      enddo

c     Count number for each QSO
      do i = 1, n_MV
         i_QSO = MV_src(i)
         n_hits_per_QSO(i_QSO) = n_hits_per_QSO(i_QSO) + 1
      enddo

c     Count number of QSOs with "hits"
      num_QSOs = 0
      do n = 1, n_QSOs
         if ( n_hits_per_QSO(n) .ge. 1 ) num_QSOs =  num_QSOs + 1
      enddo

      write (6,1000) num_QSOs, (n_hits_per_QSO(i), i=1,n_QSOs)
 1000 format(' Number useful QSOs =',i2,'; individual QSO "hits" =',8i5)

      return
      end

c===================================================================
      subroutine average_by_QSO ( n_MV, MV_src, 
     +                            MV_real, MV_imag,
     +                            n_QSOs,
     +                            MV_amp, MV_phase )

c     Takes MV QSO visibilities for each antenna within a time window
c     of a target record and averages to a single complex visibility
c     in order to prepare to fit a multiview plane 

      implicit real*8 (a-h,o-z)
      
      integer MV_src(80)
      real*8  MV_real(80), MV_imag(80)

      real*8  avg_real(8), avg_imag(8)
      real*8  MV_amp(8), MV_phase(8)

      pi = 4.d0*atan(1.d0)
      rad_deg = 180.d0/pi

      do i = 1, n_QSOs

         n = 0
         temp_real = 0.d0
         temp_imag = 0.d0

         do j = 1, n_MV

c           Sift data for each QSO...and vector average
            if ( MV_src(j) .eq. i ) then

               n = n + 1
               temp_real = temp_real + MV_real(j)
               temp_imag = temp_imag + MV_imag(j)

            endif

         enddo

         if ( n .ge. 1 ) then
               avg_real(i) = temp_real / n
               avg_imag(i) = temp_imag / n
               amp   = sqrt( avg_real(i)**2 + avg_imag(i)**2 )
               phase = atan2( avg_imag(i),avg_real(i) )*rad_deg
            else
c              Flag no data for this QSO
               avg_real(i) = -999.d0
               avg_imag(i) = -999.d0
               amp   = 0.d0
               phase = -999.d0
         endif

         write (6,1000) i, n, amp, phase
 1000    format('DEBUG: QSO#',i2,' has',i2,' visibs averaged to give',
     +          ' amp & phase =',f7.3,f7.0)

         MV_amp(i)   = amp
         MV_phase(i) = phase                   ! deg

      enddo

      return
      end

c==================================================================
      subroutine fix_wraps ( n_wraps, wrap_ants,
     +                       wrap_QSOs, wrap_turns,
     +                       early_times, late_times,
     +                       iant, t_center, MV_phase )

c     If phase wraps are specified in multiview_control.inp
c     for an antenna and QSO within a time range, fixes them
c     by adding specified turns.
      
      implicit real*8 (a-h,o-z)
      
      real*8   early_times(200), late_times(200)
      integer  wrap_ants(200), wrap_QSOs(200), wrap_turns(200)

      real*8   MV_phase(8)

      if ( n_wraps .ge. 1 ) then
         
c        Cycle through all specified wraps, checking for antenna/time
         do m = 1, n_wraps

            if ( wrap_ants(m) .eq. iant ) then

cv6            t_center_hrs = t_center * 24.d0             ! hours
cv6            Use decimal days for time range               
               if ( t_center.ge.early_times(m) .and.
     +              t_center.le.late_times(m)        ) then

                 iQSO = wrap_QSOs(m)
                 MV_phase(iQSO) = MV_phase(iQSO) + wrap_turns(m)*360.d0

               endif   ! time check

            endif   ! antenna check

         enddo   ! input wraps

      endif   ! at least 1 wrap to consider
      
      return
      end
      
c==================================================================

      subroutine compress( n_QSOs, amp_limit,
     +                     MV_amp, MV_phase, MV_phase_error, 
     +                     QSO_offset_x, QSO_offset_y,
     +                     n_good, 
     +                     MV_phase_used, MV_error_used,
     +                     QSO_used_x, QSO_used_y )

c     Removes bad/missing data and returns "compressed" arrays

      implicit real*8 (a-h,o-z)

      real*8  QSO_offset_x(8), QSO_offset_y(8)
      real*8  QSO_used_x(8), QSO_used_y(8)

      real*8  MV_amp(8), MV_phase(8), MV_phase_error(8)
      real*8  MV_phase_used(8), MV_error_used(8)

      n_good = 0
      do n = 1, n_QSOs

         if ( MV_amp(n) .ge. amp_limit ) then

            n_good = n_good + 1

            MV_phase_used(n_good) = MV_phase(n)
            MV_error_used(n_good) = MV_phase_error(n)

            QSO_used_x(n_good)    = QSO_offset_x(n)
            QSO_used_y(n_good)    = QSO_offset_y(n)

         endif

      enddo

      return
      end

c==================================================================
      subroutine fit_plane ( num_data, data, res_err, 
     +                       offsets_x, offsets_y,
     +                       slope_x,     slope_y,     const,
     +                       err_slope_x, err_slope_y, err_const )

c     Fits a plane to data (eg, z-height) versus (x,y)-offset positions
c               model = slope_x*offset_x + slope_y*offset_y + intercept
c     with parameters     ^^^                ^^^                ^^^
 
c     Currently dimensioned for up to 8 measurements (ie, up to 8 QSOs)
c     and 3 parameters

      implicit real*8 (a-h,o-z)

      real*8       data(8), res_err(8), model(8), resids(8)
      real*8       offsets_x(8), offsets_y(8)

      real*8       params(3), partls(8,3)
      real*8       new_params(3), param_sigmas(3)
      integer      paramids(3)

      character*16 parnames(3)

      logical      print_cor

      num_params   = 3                ! hardwired parameters
      idebug       = 0
      print_cor    = .false.
      max_num_data = 8

      num_solved_for = num_params

      parnames(1) = 'X-slope[deg/deg]'
      parnames(2) = 'Y-slope[deg/deg]'
      parnames(3) = 'Offset [deg]    '

c     set initial parameter values and solve-for codes
      do n = 1, num_params
         params(n)   = 0.d0
         paramids(n) = 1
      enddo

c     Now do linear least-squares fitting...
      call calc_partials ( num_data, offsets_x, offsets_y, 
     +                     partls )

      call calc_residuals ( params, 
     +                      num_data, data, offsets_x, offsets_y,
     +                      model, resids ) 
 
      call least_squares_fit( idebug, print_cor, max_num_data,
     +                   num_params, num_data,
     +                   parnames, params, paramids, resids,
     +                   partls, res_err, new_params, param_sigmas )

      call update_params ( num_params, params, new_params )

      call calc_residuals ( params, 
     +                      num_data, data, offsets_x, offsets_y,
     +                      model, resids ) 
 
      call calc_sigsq ( num_data, num_solved_for,
     +                  resids, res_err, 
     +                  sigsq, sigma_pdf)

c     Transfer fitted parameters from array to individual variables...
      slope_x = params(1)
      slope_y = params(2)
      const   = params(3)

c     Scale parameter uncertainties by sigma_pdf (if >=4 data used)
c     otherwise double uncertainties based on assumed errors of 10deg      
      scale = 2.d0                 
      if ( sigma_pdf .gt. 0.01d0 ) scale = max(sigma_pdf,0.2d0)
      err_slope_x = param_sigmas(1)*scale
      err_slope_y = param_sigmas(2)*scale
      err_const   = param_sigmas(3)*scale

c     Printout residuals?
      if ( abs(const) .gt. 80.d0 ) then 
         print*,' '
         print*,'    Data      Model     Resids    Errors'
         do n = 1, num_data
            write (6,1000) data(n), model(n), resids(n), res_err(n)
 1000       format(4f10.3)
         enddo
         write (6,1100) sigma_pdf
 1100    format(' SQRT(reduced_chi-squared) =',f7.3)
      endif

      return
      end

c======================================================================

      subroutine calc_partials ( n_pts, offsets_x, offsets_y, 
     +                           partls )

      implicit real*8 (a-h,o-z)

      real*8       offsets_x(8), offsets_y(8)
      real*8       params(3), partls(8,3) 

      do n = 1, n_pts

            partls(n,1) = offsets_x(n) 
            partls(n,2) = offsets_y(n) 
            partls(n,3) = 1.d0 

      enddo

      return
      end

c======================================================================

      subroutine calc_residuals ( params, 
     +                num_data, data, offsets_x, offsets_y,
     +                model, resids ) 

      implicit real*8 (a-h,o-z)

      real*8       data(8), model(8), resids(8)
      real*8       offsets_x(8), offsets_y(8)
      real*8       params(3) 

      do n = 1, num_data

         model(n) = params(1)*offsets_x(n) +
     +              params(2)*offsets_y(n) +
     +              params(3)

         resids(n)= data(n) - model(n)

      enddo

      return
      end

c======================================================================
	subroutine update_params(num_params, params, new_params)

c       Takes parameter values before (params) and after (new_params)
c       least squares fit; adjusts the changes by the "gain";
c       and replaces the adjusted new parameters in the "params" array.

	implicit real*8 (a-h, o-z)

	real*8		params(3), new_params(3)
 
        gain = 1.0                     ! hardwired gain parameter

	do j = 1, num_params

		delta_param = new_params(j) - params(j) 
		trial_param = params(j) + gain * delta_param
		params(j) = trial_param

	enddo

	return
	end

c======================================================================
	subroutine calc_sigsq (num_data, num_solved_for,
     +			resids, res_err, 
     +			sigsq, sigma_pdf)

C	Calculate new "sigsq"...

	implicit real*8	(a-h,o-z)

	real*8		resids(8), res_err(8)


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

c======================================================================
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

      REAL*8    B(3,3),BSTORE(3,3),XIDENT(3,3),
     +          COR(3,3),ERROR(3),STDEV(3),PR(3),
     +          XHAT(3),SNEW(8),PRODCT(8),xhat2(3)

      integer   ID(3)

      DIMENSION LL(3),MM(3)
      DATA      MSIZE /3/

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
ccccc        call uncollapse ( numxes, id, xhat2 )
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

c=======================================================================
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

c===================================================================
      subroutine control_parameters ( lu_control, control_file, 
     +                   lu_print, max_QSOs, max_wraps,
     +                   n_IFs, time_window, 
     +                   n_QSOs, source_number, 
     +                   QSO_offset_x, QSO_offset_y, 
     +                   n_wraps, wrap_ants, wrap_QSOs,
     +                   wrap_turns, early_times, late_times )
      
      implicit real*8 (a-h,o-z)

      real*8  QSO_offset_x(8), QSO_offset_y(8)

      character*48    control_file

      integer         source_number(8)

      character*13    name, QSO_name(8)

      character*1     c1

      real*8          early_times(200), late_times(200)
      integer         wrap_ants(200), wrap_QSOs(200), wrap_turns(200)

      
      call open_ascii_file ( lu_control, control_file, lu_print )

c     Read through comments...
      ieof = 0
      c1   = '!'
      do while ( ieof.ge.0 .and. c1.eq.'!' )
         read (lu_control,*,iostat=ieof) c1 
      enddo

c     Read in # of IFs (used to skip reads)
      backspace ( unit=lu_control )
      read (lu_control,*) n_IFs

c     Read in time window for multiview...
      c1   = '!'
      do while ( ieof.ge.0 .and. c1.eq.'!' )
         read (lu_control,*,iostat=ieof) c1 
      enddo
      backspace ( unit=lu_control )
      read (lu_control,*) time_window

c     Read QSO offsets...
      c1   = '!'
      do while ( ieof.ge.0 .and. c1.eq.'!' )
         read (lu_control,*,iostat=ieof) c1 
      enddo
      backspace ( unit=lu_control )
      
      c1   = ' '
      n    = 0
      do while ( ieof.ge.0 .and. c1.ne.'!' )
         read (lu_control,*,iostat=ieof) c1 
         if ( c1 .ne. '!' ) then
            backspace ( unit=lu_control )
            read (lu_control,*,iostat=ieof) number, name, dx, dy 
            if ( ieof.ge.0 .and. n.lt.max_QSOs) then
               n = n + 1
               source_number(n)= number
               QSO_name(n)     = name
               QSO_offset_x(n) = dx ! deg
               QSO_offset_y(n) = dy ! deg
            endif
         endif
      enddo
      n_QSOs = n

c     Read in phase-wrap fixes
      c1   = '!'
      do while ( ieof.ge.0 .and. c1.eq.'!' )
         read (lu_control,*,iostat=ieof) c1 
      enddo
      if ( ieof .ge. 0 )  backspace ( unit=lu_control )

      m   = 0
      do while ( ieof .ge. 0 )
         read (lu_control,*,iostat=ieof) iant, iQSO,
     +                                   t_early, t_late, iturns
         if ( ieof.ge.0 .and. m.lt.max_wraps ) then
            m = m + 1
            wrap_ants(m)   = iant           ! antenna number
            wrap_QSOs(m)   = iQSO           ! QSO number
            early_times(m) = t_early        ! early time (decimal days) 
            late_times(m)  = t_late         ! late time
            wrap_turns(m)  = iturns         ! turns of phase
         endif
      enddo
      n_wraps = m

      close(unit=lu_control)

c     Document...
      if ( n_QSOs .ge. 2 ) then
            do n = 1, n_QSOs
               write (lu_print,1000) QSO_name(n), source_number(n),  
     +                               QSO_offset_x(n), QSO_offset_y(n)
 1000          format(' MultiView QSO ',a13,' #:',i3,
     +                '; (X,Y) offset (deg)=',2f6.2)
            enddo
         else
            print*,' STOP: Need >=2 Multiview QSOs. Found ',n_QSOs
            stop
      endif

      if ( n_wraps .ge. 1 ) then
            do m = 1, n_wraps
               write (lu_print,2000) wrap_ants(m), wrap_QSOs(m),
     +              early_times(m), late_times(m), wrap_turns(m)
 2000          format(/' Phase wrap fixes: ','Ant #',i2,'; QSO',i3,
     +            '; between times (days)',2f6.2,' add',i3,' turns')
            enddo
      endif

      
      return
      end

c======================================================================
	subroutine open_ascii_file ( lu_in, ascii_file, lu_print )

C	Opens ascii file and reads and prints comments (if first character is "!") 
C       Leaves file open, ready to read data.

	character*48       ascii_file

	character*80       comment
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
 1020	   format(a)

	       if ( c_1 .eq. '!' ) then
		   write (6,1030) comment
 1030		   format(1x,a)
		 else
		   backspace (unit=lu_in)
	       endif

	enddo

	write (6,1040)                ! print a blank line at end
 1040	format(' End of Comment lines',/1x)

	return
	end


c================================================================= 
      subroutine two_point_phase ( x1,y1,phs1, x2,y2,phs2,
     +                             phs, phs_slope_x, phs_slope_y )

c     inputs: phases at two (x,y) points   
c     output: interpolated phase at (x,y) nearest (0,0) 

      implicit real*8 (a-h,o-z)

      dx = x2 - x1                          ! deg 
      dy = y2 - y1

      if ( abs(dx) .lt. 1.d-06 ) dx = 1.d-06
      if ( abs(dy) .lt. 1.d-06 ) dy = 1.d-06

      slope = dy / dx

      x = 1.d0/(slope + 1.d0/slope) * (-y1 + slope*x1)
      y = -(1.d0/slope)*x

      d1 = sqrt( (x1-x)**2 + (y1-y)**2 )
      d2 = sqrt( (x2-x)**2 + (y2-y)**2 )
      d  = d1 + d2                          ! deg 

      del_p = phs2 - phs1

      phs = phs1 + del_p*(d1/d)             ! deg

c     Only have one phase slope...arbitrarily put in "x" slot      
      phs_slope   = del_p/d                 ! deg/deg 
      phs_slope_x = phs_slope               ! deg/deg  
      phs_slope_y = 0.d0                    ! deg/deg 

      return
      end


