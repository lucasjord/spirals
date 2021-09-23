      program pconv_test
      
c     Fixes VLBI data with Jones matrices; specifically 
c     mixed polarizations to a circular basis

c     V3a: started with v3; changed input from R,I to A,phs
c     V3b: input a single baseline to output; this allows
c          a long list of data to be examined in time order
c     V4:  rotate linears (X,Y) -> (H,V)
c     V5:  various simplifations to do only the conversion
c          of raw data to a circular basis
c     ================================================
c     V6:  adapted sim_obs_5.f to apply to real data
c     V7:  added phase shifts after PCONV
c     V8:  more changes to phase shifts after PCONV to deal
c          with "non-standard" feed polarizations
c          (eg, (L,R) instead of (R,L) or (V,H) insted of (H,V)
c     V9:  skipped
c     V10: using algorithms from sim_obs_10.f
c     V11:   "      "        "   sim_obs_11.f, which handles   
c          linear*linear correctly.
c     V12a: add grid-search for electronic phases for linear*linear
c     V12b: 12a had a bug which appeared as PAs changing by simply
c           re-running the program.  Caused by month vs mmmm.
c     V13:  changed format from real,imag,wt to amp,iphs,wt
c           now reading in the 0730 data in "amp_phas.uvprt"
c           and experimenting with KE's polarization orientation.

      implicit real*8 (a-h,o-z)

c     Measured visibilities
      real*8   Em_real(2,2), Em_imag(2,2)

      real*8   Em_amp(2,2), Em_phs(2,2)
      integer iEm_phs(2,2)

      real*8   temp_real(2,2), temp_imag(2,2)

c     Real visibilities
      real*8  E_amp(2,2), E_phs(2,2), E_real(2,2), E_imag(2,2)

c     Jones (X,Y) to (H,V) matrices
      real*8  g_x(4), g_y(4), psi_a(4), psi_b(4), PA(4)

c     Feed type flags (antenna dependent)
      logical linear_flag(4), RCP_LCP(4), bad_correlator

c     Modified manual phase cals (baseline dependent)                   
      real*8  phs_cal_RR(6),phs_cal_LL(6),phs_cal_RL(6),phs_cal_LR(6)

      character*24  descriptor
 
      character*48  file_name
 
      character*1   slash
      
      character*12  time
      
      integer year, month, day
      common /obs_data/ ra, dec, year, month, day


      pi = 4.d0*atan(1.d0)
      degrad = pi/180.d0

      lu_data = 8
      
c     -------------------------------------------------------
c     Raw correlator data from S001e: 

      num_ants = 4              ! important to get this to match the data
      num_blines = num_ants*(num_ants-1)/2

      n_refant = 1              ! must have circular feeds              

c     Day 0 of observation:
      year  = 2020
      month =   04
      day   =   24

c     Source coordinates
      ra  =  1.1257963530D+02 * degrad                ! rad
      dec = -1.1686833503D+01 * degrad                ! rad

c     Use first scans for a "modified manual phase-cal"
c     Get modified manual phase cal (mmpc) data...                      
      file_name = 'amp_phs.uvprt'
      open ( unit=lu_data, file=file_name )

c     Get modified manual phase cal (mmpc) data...                      
      call get_mmpc_data ( lu_data, num_blines, n_refant,
     +              linear_flag,
     +              phs_cal_RR, phs_cal_LL, phs_cal_RL, phs_cal_LR )
      
c     -------------------------------------------------------
c     Read in data
      i_d  = 0
      ieof = 0
      do while ( ieof .ge. 0 )

c        test read a line, looking for a slash in column 3
         read (lu_data,1000,iostat=ieof) slash
 1000    format(2x,a1)
         if ( ieof.ge.0 .and. slash.eq.'/' ) then

c           have a data line...
            backspace ( unit=lu_data )

c     Example input line for real/imag uvprt output...
c12345678901234567890123456789012345678901234567890123456789012345678901234567890
c0/03:45:03.2  1- 2  -26434.93   17514.16  -20470.73   2.037   4 0.0014   2.224  -1 0.0020   1.917 176 0.0016   1.929  15 0.0018


            read (lu_data,2000) time, nant1, nant2, 
     +                          Em_amp(1,1), iEm_phs(1,1),
     +                          Em_amp(2,2), iEm_phs(2,2),
     +                          Em_amp(1,2), iEm_phs(1,2),
     +                          Em_amp(2,1), iEm_phs(2,1)
 2000       format(1x,a12,1x,i2,1x,i2,33x,4(f8.3,i4,7x))

c           Convert to real/imag (to preserve historical code)
            Em_phs(1,1) = iEm_phs(1,1)
            Em_phs(2,2) = iEm_phs(2,2)
            Em_phs(1,2) = iEm_phs(1,2)
            Em_phs(2,1) = iEm_phs(2,1)
            call ap2ri ( Em_amp,  Em_phs,
     +                   Em_real, Em_imag )
            
            call check_correlators ( Em_real, Em_imag,
     +                               bad_correlator )
            
            if ( .not.bad_correlator ) then

               i_d = i_d + 1
               
c              Print out input data
               call ri2ap ( Em_real, Em_imag,
     +                      Em_amp,  Em_phs )
               write (6,3000) time, nant1, nant2,
     +                Em_amp(1,1),Em_phs(1,1), Em_amp(2,2),Em_phs(2,2),
     +                Em_amp(1,2),Em_phs(1,2), Em_amp(2,1),Em_phs(2,1)
 3000          format(1x,a12,1x,i2,1x,i2,18x,4(f8.3,f8.0,5x))

c              Get antenna calibration information, and put in arrays
               call ant_info ( nant1, time, 
     +                         linear_flag, RCP_LCP,
     +                         psi_a, psi_b,
     +                         g_x, g_y, PA )
               call ant_info ( nant2, time,
     +                         linear_flag, RCP_LCP,
     +                         psi_a, psi_b,
     +                         g_x, g_y, PA )

c          write(6,1239) time,nant1,nant2,PA     
c 1239     format('DEBUG: ',a10,3x,2i3,5x,4f10.0)
          
c              Remove manual phase-cal from all 4 correlator phases
c              and add in appropriate 90deg phase shifts
               call modified_manphscal ( nant1, nant2, linear_flag,
     +                                   psi_a, psi_b,
     +                                   phs_cal_RR, phs_cal_LL,
     +                                   phs_cal_RL, phs_cal_LR,
     +                                   Em_phs)

               call ap2ri ( Em_amp,  Em_phs,
     +                      Em_real, Em_imag )

c              Convert to circular basis...
               call pconv ( Em_real, Em_imag, 
     +                      nant1, nant2, linear_flag, g_x, g_y,
     +                      E_real, E_imag )

               call ri2ap ( E_real, E_imag,
     +                      E_amp,  E_phs )

c              Correct for linear feed order (H,V) or (V,H)
               call fix_linear_90deg ( nant1, nant2,
     +                                 linear_flag, RCP_LCP, psi_a,
     +                                 E_phs )

c              Remove parallactic angle effect on phase in 
c              converted circular basis
               PA1 = PA(nant1)
               PA2 = PA(nant2)
               call fix_PA ( PA1, PA2,
     +                       E_phs )

c              Unwrap phases (place between -180 and +180 deg)
               call unwrap ( E_phs )

c              print out converted data              
               write (6,4000) time, nant1, nant2, PA1, PA2,
     +                  E_amp(1,1),E_phs(1,1), E_amp(2,2),E_phs(2,2),
     +                  E_amp(1,2),E_phs(1,2), E_amp(2,1),E_phs(2,1)
 4000          format(1x,a12,1x,i2,1x,i2,2x,2f7.0,2x,4(f8.3,f8.0,5x))
               print*,' '

            endif   ! good data check

         endif   ! eof check

      enddo   ! reading all data

      close ( unit=lu_data )

      end

c==================================================================
      subroutine get_mmpc_data ( lu_data, num_blines, n_refant,
     +           linear_flag,
     +           phs_cal_RR, phs_cal_LL, phs_cal_RL, phs_cal_LR )

c     Use first group of scans for modified manual phase-cals;
c     these are baseline dependent and they must be in the correct
c     (ie, 1-2, 1-3, 1-4, 2-3, 2-4, 3-4). 
c     Also, for now using a poor short-cut algorithm with refant=1
c     so refant < any linear-feed antenna number

      implicit real*8 (a-h,o-z)

      real*8  Em_amp(2,2), Em_phs(2,2)
      integer iEm_phs(2,2)

c     Following are antenna dependent arrays      
      real*8  g_x(4), g_y(4), psi_a(4), psi_b(4), PA(4)
      logical linear_flag(4), RCP_LCP(4)

c     Following are baseline dependent arrays
      real*8 amp_obs_RR(6), amp_obs_LL(6), amp_obs_RL(6), amp_obs_LR(6)
      real*8 phs_obs_RR(6), phs_obs_LL(6), phs_obs_RL(6), phs_obs_LR(6)

      real*8 pre_cal_RR(6), pre_cal_LL(6), pre_cal_RL(6), pre_cal_LR(6)
      real*8 phs_cal_RR(6), phs_cal_LL(6), phs_cal_RL(6), phs_cal_LR(6)

      integer na_1(6), na_2(6)

      character*1  slash
      character*12  time
      integer year, month, day
      common /obs_data/ ra, dec, year, month, day
      
      ieof = 0
      i_d  = 0
      do while ( ieof.ge.0 .and. i_d.lt.num_blines )

c        test read a line, looking for a slash in column 3
         read (lu_data,1000,iostat=ieof) slash
 1000    format(2x,a1)
         if ( ieof.ge.0 .and. slash.eq.'/' ) then

c           have a data line...
            backspace ( unit=lu_data )

c           Input a UVPRT record
            read (lu_data,2000) time, nant1, nant2, 
     +                          Em_amp(1,1), iEm_phs(1,1),
     +                          Em_amp(2,2), iEm_phs(2,2),
     +                          Em_amp(1,2), iEm_phs(1,2),
     +                          Em_amp(2,1), iEm_phs(2,1)
 2000       format(1x,a12,1x,i2,1x,i2,33x,4(f8.3,i4,7x))

c           Convert phases from integer to real*8
            Em_phs(1,1) = iEm_phs(1,1)
            Em_phs(2,2) = iEm_phs(2,2)
            Em_phs(1,2) = iEm_phs(1,2)
            Em_phs(2,1) = iEm_phs(2,1)

            i_d = i_d + 1

c           Store antenna numbers for mmpc scan.
            na_1(i_d)    = nant1
            na_2(i_d)    = nant2

c           Get antenna calibration information and calculate
c           parallactic angles
            call ant_info ( nant1, time,
     +                      linear_flag, RCP_LCP,
     +                      psi_a, psi_b,
     +                      g_x, g_y, PA )
            call ant_info ( nant2, time,
     +                      linear_flag, RCP_LCP,
     +                      psi_a, psi_b,
     +                      g_x, g_y, PA )

c            write(6,4356) i_d,time,nant1,nant2,PA
c 4356       format('DEBUG: ',i2,3x,a10,3x,2i2,3x,4f7.1)
            
c           Store observed amp and phases in baseline arrays
            amp_obs_RR(i_d) = Em_amp(1,1)
            amp_obs_LL(i_d) = Em_amp(2,2)
            amp_obs_RL(i_d) = Em_amp(1,2)
            amp_obs_LR(i_d) = Em_amp(2,1)

            phs_obs_RR(i_d) = Em_phs(1,1)
            phs_obs_LL(i_d) = Em_phs(2,2)
            phs_obs_RL(i_d) = Em_phs(1,2)
            phs_obs_LR(i_d) = Em_phs(2,1)

         endif

      enddo

      rewind ( unit=lu_data )

c     Check that we got all 6 baselines from 4 antennas... 
      if ( i_d .ne. num_blines ) then
         print*,' STOP: number of mmpc recs (',i_d,
     +        ') < number of baselines (',num_blines
         STOP
      endif

c     ----------------------------
c     Pre_calibrate: as if a normal manual phase cal already done
      call get_pre_calibrate (
     +         phs_obs_RR, phs_obs_LL, phs_obs_RL, phs_obs_LR,
     +         pre_cal_RR, pre_cal_LL, pre_cal_RL, pre_cal_LR )

c     Except zero phases for a linear*linear baseline
      do n_b = 1, num_blines
         n1 = na_1(n_b)
         n2 = na_2(n_b)
         if ( linear_flag(n1) .and. linear_flag(n2) ) then
            pre_cal_RR(n_b)=0.d0
            pre_cal_LL(n_b)=0.d0
            pre_cal_RL(n_b)=0.d0
            pre_cal_LR(n_b)=0.d0
         endif
      enddo

c     Now actually pre-calibrate the observed phases
      do n_b = 1, num_blines
         call pre_calibrate ( n_b,
     +         pre_cal_RR, pre_cal_LL, pre_cal_RL, pre_cal_LR,
     +         phs_obs_RR, phs_obs_LL, phs_obs_RL, phs_obs_LR )
      enddo
      
c     ----------------------------
c     Go through and generate mmpc phases... 
      do n = 1, num_blines

         n1 = na_1(n)
         n2 = na_2(n)
         PAdiff = PA(n2) - PA(n1)             ! deg

c        ------------------------------         
c        For circular*circular, get electronic phase by taking 
c        observed phase and correcting for PAdiff ... 
         if (.not.linear_flag(n1) .and. .not.linear_flag(n2)) then
            phs_cal_RR(n) = phs_obs_RR(n) - PAdiff
            phs_cal_LL(n) = phs_obs_LL(n) + PAdiff
            phs_cal_RL(n) = phs_obs_RL(n) + PAdiff
            phs_cal_LR(n) = phs_obs_LR(n) - PAdiff
         endif

c        ------------------------------         
c        For mixed polarization data...
         if ( .not.linear_flag(n1) .and. linear_flag(n2) ) then
            phs_cal_RR(n) = phs_obs_RR(n) - PAdiff
            phs_cal_LL(n) = phs_obs_LL(n) + PAdiff
            phs_cal_RL(n) = phs_obs_RL(n) - PAdiff
            phs_cal_LR(n) = phs_obs_LR(n) + PAdiff
         endif
         if ( linear_flag(n1) .and. .not.linear_flag(n2) ) then
            phs_cal_RR(n) = phs_obs_RR(n) - PAdiff
            phs_cal_LL(n) = phs_obs_LL(n) + PAdiff
            phs_cal_RL(n) = phs_obs_RL(n) + PAdiff
            phs_cal_LR(n) = phs_obs_LR(n) - PAdiff
         endif

c        ------------------------------         
c        But for linear*linear want only feed electronic phases,
c        which cannot come from the observed phases,
c        because they are either 0 or 180deg
         if ( linear_flag(n1) .and. linear_flag(n2) ) then
            n_bline = n
            call grid_search ( n1, n2, n_bline, 
     +           PA, g_x, g_y, linear_flag, RCP_LCP, psi_a, psi_b,
     +           amp_obs_RR, amp_obs_LL, amp_obs_RL, amp_obs_LR,
     +           phs_obs_RR, phs_obs_LL, phs_obs_RL, phs_obs_LR,
     +           phs_cal_RR, phs_cal_LL, phs_cal_RL, phs_cal_LR )
         endif

         call unwrap_1D ( n, phs_cal_RR )
         call unwrap_1D ( n, phs_cal_LL )
         call unwrap_1D ( n, phs_cal_RL )
         call unwrap_1D ( n, phs_cal_LR )
         write (6,1100) n1, n2,
     +                  phs_cal_RR(n), phs_cal_LL(n),
     +                  phs_cal_RL(n), phs_cal_LR(n),
     +                  psi_a(n1), psi_a(n2)
 1100    format(' Modified manual phase-cals for baseline ',
     +          i1,'-',i1,': ',4f7.0,'; psi_a values: ',2f7.0)

      enddo

      return
      end

c==================================================================
      subroutine get_pre_calibrate ( phs_RR, phs_LL, phs_RL, phs_LR,
     +                               cal_RR, cal_LL, cal_RL, cal_LR )

c     Sets up for a "standard" (antenna-based) manual phase-cal
c     Assumes baselines are in standard order
c     Assumes refant is #1
c     Input 6 baseline observed phases
c     Output 6 baseline calibration phases

      implicit real*8 (a-h,o-z)

      real*8 phs_RR(6), phs_LL(6), phs_RL(6), phs_LR(6)
      real*8 cal_RR(6), cal_LL(6), cal_RL(6), cal_LR(6)

      cal_RR(1) = phs_RR(1)                      ! bline 1-2
      cal_RR(2) = phs_RR(2)                      ! 1-3 
      cal_RR(3) = phs_RR(3)                      ! 1-4
      cal_RR(4) = phs_RR(2)-phs_RR(1)            ! 2-3 
      cal_RR(5) = phs_RR(3)-phs_RR(1)            ! 2-4
      cal_RR(6) = phs_RR(3)-phs_RR(2)            ! 3-4

      cal_LL(1) = phs_LL(1)                      ! bline 1-2
      cal_LL(2) = phs_LL(2)                      ! 1-3
      cal_LL(3) = phs_LL(3)                      ! 1-4 
      cal_LL(4) = phs_LL(2)-phs_LL(1)            ! 2-3 
      cal_LL(5) = phs_LL(3)-phs_LL(1)            ! 2-4 
      cal_LL(6) = phs_LL(3)-phs_LL(2)            ! 3-4 

      cal_RL(1) = phs_RL(1)                      ! bline 1-2
      cal_RL(2) = phs_RL(2)                      ! 1-3
      cal_RL(3) = phs_RL(3)                      ! 1-4 
      cal_RL(4) = phs_RL(2)-phs_RL(1)            ! 2-3 
      cal_RL(5) = phs_RL(3)-phs_RL(1)            ! 2-4
      cal_RL(6) = phs_RL(3)-phs_RL(2)            ! 3-4 

      cal_LR(1) = phs_LR(1)                      ! bline 1-2
      cal_LR(2) = phs_LR(2)                      ! 1-3 
      cal_LR(3) = phs_LR(3)                      ! 1-4
      cal_LR(4) = phs_LR(2)-phs_LR(1)            ! 2-3 
      cal_LR(5) = phs_LR(3)-phs_LR(1)            ! 2-4
      cal_LR(6) = phs_LR(3)-phs_LR(2)            ! 3-4 

      return
      end

c==================================================================
      subroutine pre_calibrate ( n_bline,
     +                           cal_RR, cal_LL, cal_RL, cal_LR,
     +                           phs_RR, phs_LL, phs_RL, phs_LR )

      implicit real*8 (a-h,o-z)

      real*8 phs_RR(6), phs_LL(6), phs_RL(6), phs_LR(6)
      real*8 cal_RR(6), cal_LL(6), cal_RL(6), cal_LR(6)

      n = n_bline

c      write (6,1000) n, phs_RR(n), phs_LL(n), phs_RL(n), phs_LR(n)
c      write (6,1000) n, cal_RR(n), cal_LL(n), cal_RL(n), cal_LR(n)

      phs_RR(n_bline) = phs_RR(n_bline) - cal_RR(n_bline)
      phs_LL(n_bline) = phs_LL(n_bline) - cal_LL(n_bline)
      phs_RL(n_bline) = phs_RL(n_bline) - cal_RL(n_bline)
      phs_LR(n_bline) = phs_LR(n_bline) - cal_LR(n_bline)

c      write (6,1000) n, phs_RR(n), phs_LL(n), phs_RL(n), phs_LR(n)
c 1000 format(i5,5x,4f7.0)
c      print*,' '

      return
      end
      
c================================================================== 
      subroutine grid_search ( nant1, nant2, n_bline, 
     +           PA, g_x, g_y, linear_flag, RCP_LCP, psi_a, psi_b,
     +           amp_obs_RR, amp_obs_LL, amp_obs_RL, amp_obs_LR,
     +           phs_obs_RR, phs_obs_LL, phs_obs_RL, phs_obs_LR,
     +           phs_cal_RR, phs_cal_LL, phs_cal_RL, phs_cal_LR )

c     Using a grid of phases (1 for each correlator), to calibrate and
c     convert one scan on a linear*linear baseline to a circular basis.
c     It keeps track best phases based on maximizing converted RR and
c     LL amplituded.

c     Since only the phase **differences** between the feeds at each
c     antenna are important, can arbitrarily set the first ("R")
c     feed phases to zero

      implicit real*8 (a-h,o-z)
      
c     Jones (X,Y) to (H,V) matrices
      real*8  g_x(4), g_y(4), psi_a(4), psi_b(4), PA(4)

c     Feed type flags
      logical linear_flag(4), RCP_LCP(4)

c     Observed amp and phases on baselines for one scan
      real*8  amp_obs_RR(6),amp_obs_LL(6),amp_obs_RL(6),amp_obs_LR(6)
      real*8  phs_obs_RR(6),phs_obs_LL(6),phs_obs_RL(6),phs_obs_LR(6)

c     Modified manual phase cals for baseline
      real*8  phs_try_RR(6),phs_try_LL(6),phs_try_RL(6),phs_try_LR(6)
      real*8  phs_cal_RR(6),phs_cal_LL(6),phs_cal_RL(6),phs_cal_LR(6)

c     Matrix formatted visibilities
      real*8   Em_amp(2,2), Em_phs(2,2), Em_real(2,2), Em_imag(2,2)
      real*8   E_amp(2,2), E_phs(2,2), E_real(2,2), E_imag(2,2)


c     Transfer amplitudes to matrix format            
      Em_amp(1,1) = amp_obs_RR(n_bline)
      Em_amp(2,2) = amp_obs_LL(n_bline)
      Em_amp(1,2) = amp_obs_RL(n_bline)
      Em_amp(2,1) = amp_obs_LR(n_bline)
      
c     Since only the phase **differences** between the feeds 
c     at each antenna are important, can arbitrarily set the first (R)
c     feed phases to zero.     
      phs_R1 = 0.d0
      phs_R2 = 0.d0

c     Only need to search over 180deg range
      ratio_max = 0.d0
      do ii = 1, 181
         phs_L1 = float(ii-91)

         do jj = 1, 181
            phs_L2 = float(jj-91)
            
            phs_try_RR(n_bline) = phs_R2 - phs_R1
            phs_try_LL(n_bline) = phs_L2 - phs_L1
            phs_try_RL(n_bline) = phs_L2 - phs_R1
            phs_try_LR(n_bline) = phs_R2 - phs_L1

c           Transfer/restore phases (which will be changed in
c           subroutine modified_manphscal)        
            Em_phs(1,1) = phs_obs_RR(n_bline)
            Em_phs(2,2) = phs_obs_LL(n_bline)
            Em_phs(1,2) = phs_obs_RL(n_bline)
            Em_phs(2,1) = phs_obs_LR(n_bline)

c           Start conversion to circular basis...
c           Remove trial phase-cal from all 4 correlator phases
            call modified_manphscal ( nant1, nant2, 
     +           linear_flag, psi_a, psi_b,
     +           phs_try_RR, phs_try_LL, phs_try_RL, phs_try_LR,
     +           Em_phs )

            call ap2ri ( Em_amp,  Em_phs,
     +                   Em_real, Em_imag )

c           Convert to circular basis...
            call pconv ( Em_real, Em_imag,
     +                   nant1, nant2, linear_flag, g_x, g_y, 
     +                   E_real, E_imag )

            call ri2ap ( E_real, E_imag,
     +                   E_amp,  E_phs )

c           Want phases that give maximum converted RR and LL amplitudes
            amp   = E_amp(1,1) + E_amp(2,2)
            ratio = amp / ( E_amp(1,2)+E_amp(2,1) )
            if ( ratio .gt. ratio_max ) then
               ratio_max = ratio
               
c               print*,ii,jj,ratio,E_amp(1,1),E_amp(2,2)
               
c              Store trial phase to keep the best one
               phs_cal_RR(n_bline) = phs_try_RR(n_bline)
               phs_cal_LL(n_bline) = phs_try_LL(n_bline)
               phs_cal_RL(n_bline) = phs_try_RL(n_bline)
               phs_cal_LR(n_bline) = phs_try_LR(n_bline)
            endif

         enddo

      enddo

      write (6,1000) nant1, nant2, n_bline,
     +               phs_cal_RR(n_bline),phs_cal_LL(n_bline),
     +               phs_cal_RL(n_bline),phs_cal_LR(n_bline)
 1000 format(' Sub grid_search: baseline',i2,'-',i2,' = index',i2,
     +       ': mmpc phases =',4f7.0)
      
      return
      end
      
c==================================================================     
      subroutine get_bline_number ( na_1, na_2,
     +                              n_refant, n_linear, num_blines,
     +                              n_bline )

c     TEMPORARY ROUTINE...will need to generalize later                 

      integer na_1(6), na_2(6)

      do n = 1, num_blines

         if ( na_1(n).eq.n_refant .and. na_2(n).eq.n_linear ) then
            n_bline = n
         endif

      enddo

      return
      end
      
c==================================================================
      subroutine ant_info ( nant, time,
     +                      linear_flag, RCP_LCP,
     +                      psi_a, psi_b,
     +                      g_x, g_y, PA )

c     Uses antenna number to assign antenna characteristics,
c     including feed orientations (eg, psi=90deg for Horizontal feed)
c     and calculates parallactic angle at time of observation.
     .

c     ******* Experiment specific antenna numbering and information.
c     Will need changing ************

      implicit real*8 (a-h,o-z)
      
      real*8  julian_day

c     Jones (X,Y) to (H,V) matrices
      real*8  g_x(4), g_y(4), psi_a(4), psi_b(4), PA(4)

c     Feed type flags
      logical linear_flag(4), RCP_LCP(4)

      integer day_of_obs
      
      character*12  time

      integer year, month, day
      common /obs_data/ ra, dec, year, month, day

      pi = 4.d0*atan(1.d0)
      hr_rad = pi/12.d0
      sidereal_rate = 1.002737923d0	

c     decode observing UT from AIPS printout...
      read (time,1000) iday, ihr, imin, sec
 1000 format(i1,1x,i2,1x,i2,1x,f4.1)

      day_of_obs = day + iday

c      write(6,1234) year,month,day,day_of_obs
c 1234 format('DEBUG: ',3i5,5x,i5)
      
      call get_gst ( year, month, day_of_obs, 
     +               julian_day, gmst0_hr )

      gmst0_rad  = gmst0_hr * hr_rad                      ! rad

      ut_hr   = ihr + imin/60.d0 + sec/3600.d0            ! hr
      ut_rad  = ut_hr * hr_rad                            ! rad
      stm     = gmst0_rad + ut_rad * sidereal_rate        ! GST (rads)

      pang = -999.d0
      
      if ( nant .eq. 1 ) then
c        Antenna 1=CD: circular 
         linear_flag(nant) = .false.
         RCP_LCP(nant) = .true.           ! flag for first receiver
         psi_a(nant) = 0.d0
         psi_b(nant) = 0.d0
         g_x(nant) = 1.00d0
         g_y(nant) = 0.98d0
         xm = -3753443.6952d0   ! Ceduna 30-m
         ym = -3912709.8136d0
         zm = -3348066.4774d0
	 call calc_parang( xm, ym, zm, ra, dec, stm, 
     +                     azimuth, elevation, pang)
      endif
      
      if ( nant .eq. 2 ) then
c        Antenna 2=HB: linear 
         linear_flag(nant) = .true.
         RCP_LCP(nant) = .true.           ! not used, but set anyway
         psi_a(nant) = 90.d0           ! Horizontal
         psi_b(nant) =  0.d0           ! Vertical
         g_x(nant) = 0.385d0
         g_y(nant) = 0.339d0
         xm = -3949991.0621d0   ! Hobart 12-m
         ym = -2522421.2644d0
         zm = -4311707.7651d0
	 call calc_parang( xm, ym, zm, ra, dec, stm, 
     +                     azimuth, elevation, pang)
      endif
            
      if ( nant .eq. 3 ) then
c        Antenna 4=KE: linear 
         linear_flag(nant) = .true.
         RCP_LCP(nant) = .true.           ! not used, but set anyway
         psi_a(nant) =    0.d0           ! Vertical !!!! 
         psi_b(nant) =   90.d0           ! Horizontal

         g_x(nant) = 0.356d0
         g_y(nant) = 0.256d0
         xm = -4147354.9058d0   ! Katherine 12-m
         ym = -4581542.3113d0
         zm = -1573302.8326d0
	 call calc_parang( xm, ym, zm, ra, dec, stm, 
     +                     azimuth, elevation, pang)
      endif
      
      if ( nant .eq. 4 ) then
c        Antenna 5=WA: circular 
         linear_flag(nant) = .false.
         RCP_LCP(nant) = .true.           ! flag for first receiver
         psi_a(nant) = 0.d0
         psi_b(nant) = 0.d0
         g_x(nant) = 1.21d0
         g_y(nant) = 0.94d0
         xm = -5115425.8614d0   ! Warkworth 30-m
         ym =  -477880.2417d0
         zm = -3767041.9812d0
	 call calc_parang( xm, ym, zm, ra, dec, stm, 
     +                     azimuth, elevation, pang)
c        NB: replace pang with az+el for Warkworth
         pang = azimuth + elevation
      endif
      
      if ( pang .lt. -998.d0 ) then
         print*,'STOP: invalid antenna number ',nant
         STOP
      endif

      call unwrap_scalar ( pang )
      PA(nant) = pang
      
      return
      end

c==============================================================
      subroutine calc_parang( xm, ym, zm, ra, dec, stm, 
     +                        azimuth, elevation, parang)

c     Calculates source azimuth, elevation, and parallactic angle

      implicit real*8 ( a-h,o-z )

      pi = 4.d0*atan(1.d0)
      degrad = pi/180.d0

c     station coordinates
      stat_eq = sqrt( xm**2 + ym**2 )
      stat_lat = atan2( zm, stat_eq ) ! station latitude (rad)
      stat_w_long = atan2( ym, xm )   ! station west longitude (rad)

      sinlat  = sin(stat_lat)
      coslat  = cos(stat_lat)

c     local sidereal time of observation...
      sha  = stm - ra - stat_w_long   ! source hour angle (rad)
      cossha  = cos(sha) 
      sinsha  = sin(sha)

c     calculate source zenith angle...
      cosdec  = cos(dec)
      sindec  = sin(dec)
      cosza = sinlat*sindec + coslat*cosdec*cossha
      sinza = sin( acos(cosza) )

      zenith_angle = acos( cosza )/degrad ! degrees
      elevation    = 90.d0 - zenith_angle !    "

c     Calculate source azimuth angle...
      cosaz = ( sindec - cosza*sinlat ) / ( sinza*coslat )
      sinaz = -sinsha*cosdec / sinza
      az = atan2( sinaz, cosaz )/degrad ! degrees

c     Put between 0 and 360 degrees
      if ( az .lt. 0.d0 ) az = az + 360.d0
      azimuth = az

c     Finally get parallactic angle...
      sinpa = -sinaz * coslat / cosdec
      cospa = -cossha * cosaz - sinsha * sinaz * sinlat
      pa_rad = atan2( sinpa, cospa )    ! rad
      parang = pa_rad / degrad          ! deg

      return
      end

c======================================================================
      subroutine modified_manphscal ( nant1, nant2, linear_flag,
     +                                psi_a, psi_b,
     +                                phs_cal_RR, phs_cal_LL,
     +                                phs_cal_RL, phs_cal_LR,
     +                                Em_phs )

c     Input: feed information and modified manual phase-cals,
c     and measured phases (Em_phs).
      
c     Operations: removes modified manual phase-cals and for 
c     linear feeds adds/subtracts feed orientation (psi_a, psi_b).
c     Eg, psi=90deg for Horizontal or 0deg for Vertical.

c     NOTE: if both antennas are linear, the feed orientation is
c     not applied.
      
c     Output: Em_phs after operations

      implicit real*8 (a-h,o-z)

      real*8   Em_phs(2,2)

c     Antenna-dependent parameters
      logical  linear_flag(4)
      real*8   psi_a(4), psi_b(4)

c     Baseline-dependent values
      real*8   phs_cal_RR(6),phs_cal_LL(6),phs_cal_RL(6),phs_cal_LR(6)

      n1 = nant1            ! shortened name 
      n2 = nant2

c     get baseline index from antennas...                              
c     ********  Currently hardwired for 4 stations  ***********
      if ( n1.eq.1 .and. n2.eq.2 ) n_b = 1         ! baseline index
      if ( n1.eq.1 .and. n2.eq.3 ) n_b = 2         ! baseline index
      if ( n1.eq.1 .and. n2.eq.4 ) n_b = 3         ! baseline index
      if ( n1.eq.2 .and. n2.eq.3 ) n_b = 4         ! baseline index
      if ( n1.eq.2 .and. n2.eq.4 ) n_b = 5         ! baseline index
      if ( n1.eq.3 .and. n2.eq.4 ) n_b = 6         ! baseline index

c     Apply modified manual phase-cals ...         
      Em_phs(1,1) = Em_phs(1,1) - phs_cal_RR(n_b)
      Em_phs(2,2) = Em_phs(2,2) - phs_cal_LL(n_b)
      Em_phs(1,2) = Em_phs(1,2) - phs_cal_RL(n_b)
      Em_phs(2,1) = Em_phs(2,1) - phs_cal_LR(n_b)

c     Rotate phase for linear feed orientation, 
c     (ie, by 90 degrees for Horizontal feed),
c     if and only if one of the two antennas are linearly polarized 
c     (but not both).
      if ( .not.linear_flag(n1) .or. .not.linear_flag(n2) ) then

         if ( linear_flag(n1) ) Em_phs(1,1) = Em_phs(1,1) - psi_a(n1)
         if ( linear_flag(n2) ) Em_phs(1,1) = Em_phs(1,1) + psi_a(n2)

         if ( linear_flag(n1) ) Em_phs(2,2) = Em_phs(2,2) - psi_b(n1)
         if ( linear_flag(n2) ) Em_phs(2,2) = Em_phs(2,2) + psi_b(n2)

         if ( linear_flag(n1) ) Em_phs(1,2) = Em_phs(1,2) + psi_a(n1)
         if ( linear_flag(n2) ) Em_phs(1,2) = Em_phs(1,2) - psi_b(n2)

         if ( linear_flag(n1) ) Em_phs(2,1) = Em_phs(2,1) + psi_b(n1)
         if ( linear_flag(n2) ) Em_phs(2,1) = Em_phs(2,1) - psi_a(n2)

      endif

      return
      end

c=======================================================================
      subroutine fix_linear_90deg ( nant1, nant2,
     +                              linear_flag, RCP_LCP, psi_a,
     +                              E_phs )

c     Removes a 90deg phase-shift introduced in mixed polarization data
c     (ie, one is linear and the other circular) when converted to 
c     circular, and when the first linear feed is horizontal 
c     (position angle psi = 90deg EofN).

c     UNSURE WHAT TO DO WITH CROSS-HANDS DATA...FOR NOW NOTHING

      implicit real*8 (a-h,o-z)

      real*8  E_amp(2,2), E_phs(2,2)

      real*8  psi_a(4)
      logical linear_flag(4), RCP_LCP(4)

      n1 = nant1
      n2 = nant2

      if ( linear_flag(n1) .and. .not.linear_flag(n2) ) then
         c2 = 1.d0
         if ( .not.RCP_LCP(n2) ) c2 = -1.d0
         E_phs(1,1) = E_phs(1,1) + c2*psi_a(n1)               ! deg
         E_phs(2,2) = E_phs(2,2) - c2*psi_a(n1)
      endif

      if ( linear_flag(n2) .and. .not.linear_flag(n1) ) then
         c1 = 1.d0
         if ( .not.RCP_LCP(n1) ) c1 = -1.d0
         E_phs(1,1) = E_phs(1,1) - c1*psi_a(n2)                   ! deg
         E_phs(2,2) = E_phs(2,2) + c1*psi_a(n2)
      endif
      
c     Finally if (H,V)*(L,R) or (L,R)*(H,V) need to remove 2xPsi
      if ( linear_flag(n1) .and. .not.RCP_LCP(n2) ) then
         E_phs(1,1) = E_phs(1,1) - 2.d0*psi_a(n1)
         E_phs(2,2) = E_phs(2,2) + 2.d0*psi_a(n1)
      endif

      if ( linear_flag(n2) .and. .not.RCP_LCP(n1) ) then
         E_phs(1,1) = E_phs(1,1) - 2.d0*psi_a(n2)
         E_phs(2,2) = E_phs(2,2) + 2.d0*psi_a(n2)
      endif

      return
      end

c=======================================================================
      subroutine fix_PA ( PA1, PA2,
     +                    E_phs )

c     Subtracts (PA2-PA1) from RR data                       
c     Adds      (PA2-PA1) to   LL data                       

c     UNSURE WHAT TO DO WITH CROSS-HANDS DATA...FOR NOW NOTHING         

      implicit real*8 (a-h,o-z)

      real*8  PA(4), E_amp(2,2), E_phs(2,2)

      PAdiff = PA2 - PA1

      E_phs(1,1) = E_phs(1,1) - PAdiff   ! deg                          
      E_phs(2,2) = E_phs(2,2) + PAdiff

c     Unsure about the following!!!!                                    
      E_phs(1,2) = E_phs(1,2) - PAdiff   ! deg                          
      E_phs(2,1) = E_phs(2,1) + PAdiff

      return
      end
      
c=======================================================================
      subroutine fix_parang ( linear_flag, RCP_LCP, PA, psi_a, 
     +                        E_amp, E_phs )

c     This is currently not used...but some of the code may be
c     needed for non-standard cases, eg, (L,R) or (V,H) feeds.
c     ****** Will need to add antenna numbers and change
c     to antenna indexing of input arrays ********************
      
c     Assuming data is now (RR,LL), remove effects of parallactic 
c     angle rotating the phase.   For an original circular feed, 
c     the correction is to subtract PA value from RCP and add to LCP.
c     For an original linear feed, the correction is 2xPA.

c     UNSURE WHAT TO DO WITH CROSS-HANDS DATA...FOR NOW NOTHING

      implicit real*8 (a-h,o-z)

      real*8  PA(2), E_amp(2,2), E_phs(2,2), psi_a(2)

      logical linear_flag(2), RCP_LCP(2)

      f1 = 1.d0
      if ( linear_flag(1) ) then
c        Factor of 2 needed for originally linear feed and (R,L)
c        Factor of 0 needed for originally linear feed and (L,R)
         if ( RCP_LCP(2) )      f1 = 2.d0
         if ( .not.RCP_LCP(2) ) f1 = 0.d0
c        If (H,V)*(H,V), do nothing...          
         if ( linear_flag(2) .and. abs(psi_a(1)).lt.1.d0 ) f1 = 0.d0
      endif

      f2 = 1.d0
      if ( linear_flag(2) ) then
c        Factor of 2 needed for originally linear feed and (R,L)
c        Factor of 0 needed for originally linear feed and (L,R)
         if ( RCP_LCP(1) )      f2 = 2.d0
         if ( .not.RCP_LCP(1) ) f2 = 0.d0
c        Again, if (H,V)*(H,V), do nothing...  
         if ( linear_flag(1) .and. abs(psi_a(2)).lt.1.d0 ) f2 = 0.d0
      endif

c     Sign of phase correction depends on (RCP,LCP) vs (LCP,RCP)
      c1 = 1.d0
      if ( .not.linear_flag(1) .and. .not.RCP_LCP(1) ) c1 = -1.d0
      c2 = 1.d0
      if ( .not.linear_flag(2) .and. .not.RCP_LCP(2) ) c2 = -1.d0

      E_phs(1,1) = E_phs(1,1) - (f2*c2*PA(2) - f1*c1*PA(1))     ! deg
      E_phs(2,2) = E_phs(2,2) + (f2*c2*PA(2) - f1*c1*PA(1))

c     Unsure about the following!!!!
      E_phs(1,2) = E_phs(1,2) - (f1*c1*PA(1) - f2*c2*PA(2))     ! deg
      E_phs(2,1) = E_phs(2,1) + (f2*c2*PA(2) - f1*c1*PA(1))

      return
      end

c===============================================================
      subroutine fix_mmpc_PAs ( linear_flag, RCP_LCP, PAdif_mmpc,
     +                          E_amp, E_phs )

c     ****This is not used, but keep code around if needed later*****

c     Assuming data is now (RR,LL), remove effects of parallactic
c     angle rotating the phase that is left in the data from the
c     modified manual phase cal (mmpc).
c     Only need to do this for linear x circular (and vice versa) data

c     UNSURE WHAT TO DO WITH CROSS-HANDS DATA.

      implicit real*8 (a-h,o-z)

      real*8  PA(2), E_amp(2,2), E_phs(2,2)

      logical linear_flag(2), RCP_LCP(2)

c     Seem to need additional correction for linear x circular case:
c     remove parallactic angles from the modified manual phase cal scan
      if ( ( linear_flag(1) .and. .not.linear_flag(2) ) .or.
     +     ( linear_flag(2) .and. .not.linear_flag(1) ) )   then

         E_phs(1,1) = E_phs(1,1) + PAdif_mmpc          ! deg
         E_phs(2,2) = E_phs(2,2) - PAdif_mmpc

c        Unsure about the cross-hands!!!!

      endif

      return
      end
      
c=============================================================
      subroutine unwrap ( E_phs )

      implicit real*8 (a-h,o-z)

      real*8  E_phs(2,2)

      do i = 1, 2
         do j = 1, 2

            do while ( E_phs(i,j) .lt. -180.d0 )
               E_phs(i,j) = E_phs(i,j) + 360.d0
            enddo

            do while ( E_phs(i,j) .gt. +180.d0 )
               E_phs(i,j) = E_phs(i,j) - 360.d0
            enddo

         enddo
      enddo

      return
      end

c=============================================================
      subroutine unwrap_1D ( n, phases )

      implicit real*8 (a-h,o-z)

      real*8  phases(6)

      do while ( phases(n) .lt. -180.d0 )
         phases(n) = phases(n) + 360.d0
      enddo

      do while ( phases(n) .gt. +180.d0 )
         phases(n) = phases(n) - 360.d0
      enddo

      return
      end

c=============================================================
      subroutine unwrap_scalar ( phase )

      implicit real*8 (a-h,o-z)

      do while ( phase .lt. -180.d0 )
         phase = phase + 360.d0
      enddo

      do while ( phase .gt. +180.d0 )
         phase = phase - 360.d0
      enddo

      return
      end

c==================================================================
      subroutine pconv ( Em_real, Em_imag, 
     +                   nant1, nant2, linear_flag, g_x, g_y, 
     +                   E_real, E_imag )

      implicit real*8 (a-h,o-z)

      real*8  Em_real(2,2),Em_imag(2,2), E_real(2,2), E_imag(2,2)
      real*8  temp_real(2,2),  temp_imag(2,2)
      real*8  temp2_real(2,2), temp2_imag(2,2)
      real*8  temp3_real(2,2), temp3_imag(2,2)

c     Jones (X,Y) to (H,V) matrices 
      real*8  g_x(4), g_y(4), PA(4)
      real*8  Jpa1_real(2,2), Jpa1_imag(2,2) 
      real*8  Jpa2_real(2,2), Jpa2_imag(2,2) 

c     Marti-Vidal's "C" matrices to convert linear to circular
      logical linear_flag(4)

      real*8  CL2C1_real(2,2), CL2C1_imag(2,2) 
      real*8  CL2C2_real(2,2), CL2C2_imag(2,2) 
      real*8  CL2C2H_real(2,2),CL2C2H_imag(2,2), CL2C2T_imag(2,2) 

      real*8  J1_real(2,2), J1_imag(2,2)
      real*8  J1inv_real(2,2), J1inv_imag(2,2)

      real*8  J2_real(2,2),    J2_imag(2,2)
      real*8  J2T_real(2,2),   J2T_imag(2,2)
      real*8  J2H_real(2,2),   J2H_imag(2,2)
      real*8  J2Hinv_real(2,2),J2Hinv_imag(2,2) 

      character*24  descriptor

      pi = 4.d0*atan(1.d0)
      degrad = pi/180.d0      

c     Document measured visibilities with R & L phases shifted...
      descriptor = 'E-phs in (real & imag)  '
      call print_matrix (descriptor, Em_real, Em_imag)

c     -----------------------------------------------------
c     Get "hybrid" linear-to-circular corrections for both ants:
c     "C" matrix from Eq. 3.5 of Marti-Vidal (arXiv:1504.06579v1)
c     Note, this matrix cannot be inverted, so can't use it as
c     a standard Jones matrix
      call convert_L2C ( nant1, linear_flag, 
     +                   CL2C1_real, CL2C1_imag )
      descriptor = 'CL2C1 matrix'
      call print_matrix (descriptor, CL2C1_real, CL2C1_imag) 

      call convert_L2C ( nant2, linear_flag,
     +                   CL2C2_real, CL2C2_imag )
      descriptor = 'CL2C2 matrix'
      call print_matrix (descriptor, CL2C2_real, CL2C2_imag) 

c     -----------------------------------------------------
c     Get Jones gain corrections for both antennas
      call jones_wo_PA ( nant1, linear_flag, g_x, g_y,
     +                J1_real, J1_imag )
      descriptor = 'J1 matrix'
      call print_matrix (descriptor, J1_real, J1_imag)
 
c     Invert Jones matrix
      call invert_complex_2x2matrix ( J1_real,    J1_imag,
     +                                J1inv_real, J1inv_imag )
      descriptor = 'J1inv matrix'
      call print_matrix (descriptor, J1inv_real, J1inv_imag) 

      call jones_wo_PA ( nant2, linear_flag, g_x, g_y, 
     +                J2_real, J2_imag )
c     Get Hermitian (transpose and complex conjugate)
      call hermitian_complex_2x2 ( J2_real,  J2_imag,
     +                             J2H_real, J2H_imag )
      descriptor = 'J2H matrix'
      call print_matrix (descriptor, J2H_real, J2H_imag) 

c     Invert Jones matrix
      call invert_complex_2x2matrix ( J2H_real,    J2H_imag,
     +                                J2Hinv_real, J2Hinv_imag )
      descriptor = 'J2Hinv matrix'
      call print_matrix (descriptor, J2Hinv_real, J2Hinv_imag) 


c-----------------------------------------------------------
c                 Apply Jones matrix calibrations
c     First measured visib x antenna 2...
      call multiply_complex_2x2matrix ( Em_real,     Em_imag, 
     +                                  J2Hinv_real, J2Hinv_imag,
     +                                  temp_real,   temp_imag )

c     Then x antenna 1...
      call multiply_complex_2x2matrix ( J1inv_real, J1inv_imag, 
     +                                  temp_real,  temp_imag,
     +                                  temp2_real, temp2_imag )

c-----------------------------------------------------------
c                 Apply "hybrid" L2C matrix calibrations
c     1st antenna conversion
      call multiply_complex_2x2matrix ( CL2C1_real, CL2C1_imag, 
     +                                  temp2_real, temp2_imag,
     +                                  temp_real,  temp_imag )

c     2nd antenna conversion
c     Get Hermitian (transpose and complex conjugate)
      call hermitian_complex_2x2 (CL2C2_real,  CL2C2_imag,
     +                            CL2C2H_real, CL2C2H_imag )
      descriptor = 'CL2C2H (real & imag)'
      call print_matrix (descriptor, CL2C2H_real, CL2C2H_imag)


      call multiply_complex_2x2matrix ( temp_real,  temp_imag, 
     +                                  CL2C2H_real,CL2C2H_imag,
     +                                  E_real, E_imag )

c     Remove 90 phase shift from LL=(2,2) for unknown reason...
c      call ninety_deg_shift ( linear_flag, 
c     +                        temp3_real, temp3_imag,
c     +                        E_real, E_imag ) 

      descriptor = 'E-true (real & imag)'
      call print_matrix (descriptor, E_real, E_imag)


      return
      end

c----------------------------------------------------------------
      subroutine jones_wo_PA ( Nant, linear_flag, g_x, g_y, 
     +                         J_real, J_imag )

c     Calculates complex Jones matrix for linear feeds
c     for antenna (Nant) with real gain amplitudes (g_x and g_y).

c     TESTING if       
c     the X,Y feeds are rotated by parallactic angle (PA)
c     with respect to Horizontal,Vertical.
      
c     For circularly polarized antennas, PA is set to zero,
c     but the gains are still used.

      implicit real*8 (a-h,o-z)

      real*8  J_real(2,2), J_imag(2,2)

      real*8  g_x(4), g_y(4)
      logical linear_flag(4)

      character*24 descriptor

      pi = 4.d0*atan(1.d0)
      degrad = pi/180.d0

c      if ( linear_flag(Nant) ) then
c           parang = PA(Nant)      ! deg
c           parang = 0.d0  !!!!!!!!!!!!!!!!! TESTING ***************
c        else
c           parang = 0.d0
c      endif
      parang = 0              ! NO PARANG IN THIS VERSION
      pa_rad = parang * degrad             ! rad
      cospa  = cos(pa_rad)
      sinpa  = sin(pa_rad)

      xgain = g_x(Nant)
      ygain = g_y(Nant)

      J_real(1,1) = xgain * cospa
      J_real(1,2) = xgain * sinpa
      J_real(2,1) =-ygain * sinpa
      J_real(2,2) = ygain * cospa

      J_imag(1,1) = 0.d0
      J_imag(1,2) = 0.d0
      J_imag(2,1) = 0.d0
      J_imag(2,2) = 0.d0

      write (descriptor,1000) Nant, xgain, ygain, parang
 1000 format('Ant ',i2,':',2f5.2,f5.1)
      call print_matrix (descriptor, J_real, J_imag)

      return
      end


c----------------------------------------------------------------
      subroutine convert_L2C ( Nant, linear_flag,
     +                         CL2C_real, CL2C_imag )

c     Calculates complex matrix that is used to convert linearXcircular
c     to circularXcircular (CL2C) visibilities for antenna (Nant),
c     assuming Horizontal,Vertical orientations.  
C     (So, before applying L2C, rotate from X,Y to H,V by the 
c     parallactic angle.)

      implicit real*8 (a-h,o-z)

      real*8  CL2C_real(2,2), CL2C_imag(2,2)

      logical linear_flag(2)

      character*24 descriptor

      scale = 1.d0/sqrt(2.d0)

c     Check linear_flag...
      if ( linear_flag(Nant) ) then
c           Linear feeds...mix to circular
            CL2C_real(1,1) = scale
            CL2C_real(1,2) = 0.d0
            CL2C_real(2,1) = scale
            CL2C_real(2,2) = 0.d0

            CL2C_imag(1,1) = 0.d0
            CL2C_imag(1,2) =-scale
            CL2C_imag(2,1) = 0.d0
            CL2C_imag(2,2) = scale
         else
c           Circular feeds...use identity matrix
            CL2C_real(1,1) = 1.d0
            CL2C_real(1,2) = 0.d0
            CL2C_real(2,1) = 0.d0
            CL2C_real(2,2) = 1.d0

            CL2C_imag(1,1) = 0.d0
            CL2C_imag(1,2) = 0.d0
            CL2C_imag(2,1) = 0.d0
            CL2C_imag(2,2) = 0.d0
      endif

      write (descriptor,1000) Nant, linear_flag(Nant)
 1000 format('Antenna',i3,': LinFeed? ',L1)
      call print_matrix (descriptor, CL2C_real, CL2C_imag)
 

      return
      end

c---------------------------------------------------------------
      subroutine hermitian_complex_2x2 ( M_real,  M_imag,
     +                                   MH_real, MH_imag )

c     Input complex matrix: M_real  + i M_imag
c     Output Hermetian    : MH_real + i MH_imag

      real*8 M_real(2,2), M_imag(2,2), MH_real(2,2), MH_imag(2,2)

      real*8 MT_imag(2,2)

      call transpose_2x2 ( M_real,  MH_real )

      call transpose_2x2 ( M_imag,  MT_imag )
      call negate_2x2    ( MT_imag, MH_imag )   

      return
      end

c---------------------------------------------------------------
      subroutine multiply_complex_2x2matrix ( A, B, C, D,
     +                                        E, F )

c     Input a complex matrices (A+iB) and (C+iD),
c     where A, B, C and D are real 2x2 matrices.
c     Output is the product complex matrix (E+iF)

      real*8 A(2,2), B(2,2), C(2,2), D(2,2), E(2,2), F(2,2)

c     Real part matrix
      E(1,1) = A(1,1)*C(1,1) - B(1,1)*D(1,1) +
     +         A(1,2)*C(2,1) - B(1,2)*D(2,1)

      E(1,2) = A(1,1)*C(1,2) - B(1,1)*D(1,2) +
     +         A(1,2)*C(2,2) - B(1,2)*D(2,2)

      E(2,1) = A(2,1)*C(1,1) - B(2,1)*D(1,1) +
     +         A(2,2)*C(2,1) - B(2,2)*D(2,1)

      E(2,2) = A(2,1)*C(1,2) - B(2,1)*D(1,2) +
     +         A(2,2)*C(2,2) - B(2,2)*D(2,2)

c     Imaginary part matrix
      F(1,1) = A(1,1)*D(1,1) + B(1,1)*C(1,1) +
     +         A(1,2)*D(2,1) + B(1,2)*C(2,1)

      F(1,2) = A(1,1)*D(1,2) + B(1,1)*C(1,2) +
     +         A(1,2)*D(2,2) + B(1,2)*C(2,2)

      F(2,1) = A(2,1)*D(1,1) + B(2,1)*C(1,1) +
     +         A(2,2)*D(2,1) + B(2,2)*C(2,1)

      F(2,2) = A(2,1)*D(1,2) + B(2,1)*C(1,2) +
     +         A(2,2)*D(2,2) + B(2,2)*C(2,2)


      return
      end

c---------------------------------------------------------------
      subroutine invert_complex_2x2matrix ( A, B,
     +                                      C, D )

c     Input a complex matrix (A+iB). 
c     Returns its inverse (C+iD),
c     where A, B, C and D are real 2x2 matrices

c     If complex M = A + iB,  then its inverse is given by
c     Minv = (A + B*Ainv*B)inv - i(B + A*Binv*A)inv

      real*8 A(2,2), B(2,2), C(2,2), D(2,2)

      real*8 Ainv(2,2), Binv(2,2), temp1(2,2), temp2(2,2)

      call invert_2x2 ( A, Ainv )
      call invert_2x2 ( B, Binv )
      
c     Get reals 2x2 of inverse...
      call multiply_2x2 ( Ainv, B, temp1 )
      call multiply_2x2 ( B, temp1, temp2 )
      call add_2x2      ( A, temp2, temp1 )
      call invert_2x2   ( temp1, C )

c     Get imaginaries 2x2 of inverse...
      call multiply_2x2 ( Binv, A, temp1 )
      call multiply_2x2 ( A, temp1, temp2 )
      call add_2x2      ( B, temp2, temp1 )
      call invert_2x2   ( temp1, temp2 )
      call negate_2x2   ( temp2, D )

      return
      end

c     ----------------------------------------------------------
      subroutine add_2x2 (x,y,z)

      implicit real*8 (a-h,o-z)

c     adds real 2x2 matrices x + y = z

      real*8  x(2,2), y(2,2), z(2,2)

c     Index convention:
c       1,1  1,2  
c       2,1  2,2  


      do i = 1, 2

         do j = 1, 2

               z(i,j) = x(i,j) + y(i,j)

         enddo

      enddo

      return
      end

c     ----------------------------------------------------------
      subroutine negate_2x2 (x, z)

      implicit real*8 (a-h,o-z)

c     negates a real 2x2 matrice -x = z

      real*8  x(2,2), z(2,2)

c     Index convention:
c       1,1  1,2  
c       2,1  2,2  


      do i = 1, 2

         do j = 1, 2

               z(i,j) = -x(i,j)

         enddo

      enddo

      return
      end

c     ----------------------------------------------------------
      subroutine multiply_2x2 (x,y, z)

      implicit real*8 (a-h,o-z)

c     multiplies real 2x2 matrices x * y = z

      real*8  x(2,2), y(2,2), z(2,2)

c     Index convention:
c       1,1  1,2  
c       2,1  2,2  


      do i = 1, 2

         do j = 1, 2

            z(i,j) = 0.d0

            do k = 1, 2
               z(i,j) = z(i,j) + x(i,k) * y(k,j)
            enddo

         enddo

      enddo

      return
      end

c     ----------------------------------------------------------
      subroutine transpose_2x2 (M, MT)

      implicit real*8 (a-h,o-z)

c     transposes a real 2x2 matrix  M -> MT

      real*8  M(2,2), MT(2,2)

c     Index convention:
c       1,1  1,2  
c       2,1  2,2  

      MT(1,1) = M(1,1)
      MT(1,2) = M(2,1)
      MT(2,1) = M(1,2)
      MT(2,2) = M(2,2)

      return
      end

c     ----------------------------------------------------------
      subroutine invert_2x2 (M, Minv)

      implicit real*8 (a-h,o-z)

c     inverts a real 2x2 matrix:  M -> Minv

      real*8  M(2,2), Minv(2,2)

c     Index convention:
c       1,1  1,2  
c       2,1  2,2  

      a = M(1,1)
      b = M(1,2)
      c = M(2,1)
      d = M(2,2)

      det  = a*d - b*c
      if ( det .ne. 0.d0 ) then

            det1 = 1.d0 / det 

         else

c           Check for all zero matrix, which is OK,
c           but can't be inverted.
            if ( a.eq.0.d0 .and. b.eq.0.d0 .and.
     +           c.eq.0.d0 .and. d.eq.0.d0      ) then
c                 OK, since input matrix is all zeros
                  det1 = 0.d0
               else
c                 not OK...warn and stop
                  print*,'Determinant zero...STOP'
                  write (6,9000) a, b, c, d
 9000             format(/' Singular matrix:',2f10.3,
     +                   /'                 ',2f10.3)
                  stop
            endif

      endif

      Minv(1,1) = det1*d
      Minv(1,2) = det1*(-b)
      Minv(2,1) = det1*(-c)
      Minv(2,2) = det1*a

      return
      end

c     ----------------------------------------------------------
      subroutine print_matrix (descriptor, M_real, M_imag) 

      real*8  M_real(2,2), M_imag(2,2)

      character*24  descriptor

      lu_print = -1                  ! kill printout
      if ( lu_print .ge. 6 ) then
         write (6,1000) descriptor,
     +               M_real(1,1), M_real(1,2), 
     +                  M_imag(1,1), M_imag(1,2),
     +               M_real(2,1), M_real(2,2),
     +                  M_imag(2,1), M_imag(2,2)
 1000    format(/a24,':',2f10.3,5x,2f10.3, 
     +          /25x    ,2f10.3,5x,2f10.3)
      endif

      return
      end

c     ----------------------------------------------------------
      subroutine ri2ap ( E_real, E_imag,
     +                   E_amp,  E_phs )

c     Converts real,imag  to   amp,phs(deg)

      implicit real*8 (a-h,o-z)

      real*8  E_real(2,2), E_imag(2,2), E_amp(2,2), E_phs(2,2)

      pi = 4.d0*atan(1.d0)
      raddeg = 180.d0/pi      

      do i = 1, 2

         do j = 1, 2

            E_amp(i,j) = sqrt( E_real(i,j)**2 + E_imag(i,j)**2 )
            E_phs(i,j) = atan2(E_imag(i,j), E_real(i,j) )*raddeg

         enddo

      enddo

      return
      end

c     ----------------------------------------------------------
      subroutine ap2ri ( E_amp,  E_phs,
     +                   E_real, E_imag )

c     Converts  amp,phs(deg)  to  real,imag   

      implicit real*8 (a-h,o-z)

      real*8  E_amp(2,2), E_phs(2,2)

      real*8  E_real(2,2), E_imag(2,2)

      pi = 4.d0*atan(1.d0)
      degrad = pi/180.d0

      do i = 1, 2

         do j = 1, 2

            phs_rad = E_phs(i,j) * degrad
            E_real(i,j) = E_amp(i,j) * cos(phs_rad)
            E_imag(i,j) = E_amp(i,j) * sin(phs_rad)

         enddo

      enddo

      return
      end

c     ----------------------------------------------------------
      subroutine ninety_deg_shift ( linear_flag, 
     +                              Ein_real,  Ein_imag,
     +                              Eout_real, Eout_imag ) 

c     Adds or subtracts 90deg from LL=(2,2), if "linear_flag"
c     is set for antenna 1 or 2.
c     Note, if both linear, then nothing is done.

      implicit real*8 (a-h,o-z)

      real*8  Ein_real(2,2), Ein_imag(2,2)

      real*8  Eout_real(2,2), Eout_imag(2,2)

      logical linear_flag(2)

c     Default is no change...
      do i = 1, 2
         do j = 1, 2
            Eout_real(i,j) = Ein_real(i,j)
            Eout_imag(i,j) = Ein_imag(i,j)
         enddo
      enddo

c     Check if antenna 1 was linearly polarized
      if ( linear_flag(1) .and. .not.linear_flag(2) ) then
         Eout_real(2,2) =  Ein_imag(2,2)
         Eout_imag(2,2) = -Ein_real(2,2)
      endif

c     Check if antenna 2 was linearly polarized
      if ( linear_flag(2) .and. .not.linear_flag(1) ) then
         Eout_real(2,2) = -Ein_imag(2,2)
         Eout_imag(2,2) =  Ein_real(2,2)
      endif

      return
      end

c==============================================================
      subroutine check_correlators ( Em_real, Em_imag,
     +                               bad_correlator )

c     Checks for bad correlation (with 0.0 amp and 0 phs)
      
      real*8   Em_real(2,2), Em_imag(2,2)

      logical  bad_correlator

      bad_correlator = .false.

      do i = 1, 2
         do j = 1, 2
            if ( Em_real(i,j).eq.0.d0 .and. Em_imag(i,j).eq.0 ) then
               bad_correlator = .true.
            endif
         enddo
      enddo
      
      return
      end
      
c==============================================================
      SUBROUTINE get_gst (YEAR, MONTH, DAY, DATE, gmst_hr)

C Returns DATE = Julian Date and GMST = Greenwich Mean sidereal time
C at 0 hrs U.T. on day DAY/MONTH/YEAR. Accuracy not tested; GMST is
C correct at least to nearest second
C
C History:
C  1977 Aug  5 - TJP.
C  1991 May 18 - TJP.
c  2006 Sep 11 - MJR.
C-----------------------------------------------------------------------

      INTEGER YEAR
      INTEGER MONTH
      INTEGER DAY

      DOUBLE PRECISION DATE
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
      return
      end
