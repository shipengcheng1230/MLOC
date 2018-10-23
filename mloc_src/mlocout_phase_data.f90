!*****************************************************************************************
subroutine mlocout_phase_data (it)

! Standard printed output files:
!   Phase data (.phase_data)
!   Large absolute residuals (.xdat)
!   Large cluster residuals (.lres)
!   Depth phases (.depth_phase)
! Optional:
!   direct calibration data (.dcal_phase_data)
!   Lg data (.lg)
!   T-phase data (.tp)
!   Limited distance range (.oldr)

! Also fills the arrays to calculate the reading error for output to .rderr

! January 26, 1989 by eric bergman
! Modified 9/20/94 by eab.
! 11/17/2017: format change for v10.4.0, adding deployment code

   implicit none
        
   include 'mloc.inc'
   
   logical :: bdp
   integer :: it, iev, n_constrained, n_input, indx(ntmax0), i, ii, j, it1, io_dcal, lunit
   real :: prmax(0:itmax1), dts(0:itmax1), dtslim, deltiev(ntmax0), h_mult, deltkm
   real :: x_dep_constrained(nevmax), x_dep_input(nevmax), xmed, cluster_default_depth, lon_out
   character(len=1) :: bdp_flag
   character(len=4) :: rpres, reason
   character(len=8) :: rdp_phase, depth_free_fixed
   character(len=21) :: qtest
   character(len=60) :: depth_set1, depth_set2, depth_set
   character(len=100) :: outfil, outfil2, outfil3, outfil4, oldr_filnam
   character(len=120) :: gfmt, bfmt
   character(len=132) :: msg
   
   data rpres/'PRES'/
   
   it1 = it + 1
   if (nstep .eq. 0) it1 = 0 ! Forward modelling
   
   prmax(0) = wind2
   do i = 1,itmax1
      prmax(i) = wind2
   end do
   
   ! Initialize reading errors
   indexq = 0
   qres = 0.
!    do i = 1,nqmax
!       indexq(i) = 0
!       do j = 1,nevmax
!          qres(i,j) = 0.
!       end do
!    end do

   ! Output files
   outfil = trim(outfile)//'.phase_data'
   if (verbose_screen) then
      write (msg,'(3a,i3)') 'mlocout_phase_data: opening ', trim(outfil), ' on unit ', io_out
      call fyi (trim(msg))
   end if
   open (io_out,file=outfil,status='new')
   
   outfil2 = trim(outfile)//'.xdat' ! List of unflagged readings that fail PRES
   if (verbose_screen) then
      write (msg,'(3a,i3)') 'mlocout_phase_data: opening ', trim(outfil2), ' on unit ', io_xdat
      call fyi (trim(msg))
   end if
   open (io_xdat,file=outfil2,status='new')
   
   outfil3 = trim(outfile)//'.depth_phases'
   if (verbose_screen) then
      write (msg,'(3a,i3)') 'mlocout_phase_data: opening ', trim(outfil3), ' on unit ', io_depth_phase
      call fyi (trim(msg))
   end if
   open (io_depth_phase,file=outfil3,status='new')
   
   ! Output file for data used in direct calibration of the hypocentroid
   if (direct_cal) then
      outfil4 = trim(outfile)//'.dcal_phase_data'
      io_dcal = lunit()
      if (verbose_screen) then
         write (msg,'(3a,i3)') 'mlocout_phase_data: opening ', trim(outfil4), ' on unit ', io_dcal
         call fyi (trim(msg))
      end if
      open (io_dcal,file=outfil4,status='new')
   end if
               
   ! Output for limited distance range
   if (oldr_out) then
      oldr_filnam = trim(outfile)//'.oldr'
      if (verbose_screen) then
         write (msg,'(3a,i3)') 'mlocout_phase_data: opening ', trim(oldr_filnam), ' on unit ', io_oldr
         call fyi (trim(msg))
      end if
      open (io_oldr,file=oldr_filnam,status='new')
   end if 
   
   ! Format assignments for individual phase data
   if (it .eq. 0) then
      gfmt = '(a1,a5,1x,a8,1x,a1,a1,a8,1x,f4.2,f7.2,i5,f6.1,f7.2,f6.2,f6.1,1f7.2                ,t134,a8,1x,a3,1x,a8,1x,2i5)'
      bfmt = '(a1,a5,1x,a8,1x,a1,a1,a8,1x,f4.2,f7.2,i5,f6.1,f7.2,f6.2,f6.1,1f7.2,t112,a4        ,t134,a8,1x,a3,1x,a8,1x,2i5)'
   else if (it .eq. 1) then
      gfmt = '(a1,a5,1x,a8,1x,a1,a1,a8,1x,f4.2,f7.2,i5,f6.1,f7.2,f6.2,f6.1,2f7.2,t112,2f7.4,f7.2,t134,a8,1x,a3,1x,a8,1x,2i5)'
      bfmt = '(a1,a5,1x,a8,1x,a1,a1,a8,1x,f4.2,f7.2,i5,f6.1,f7.2,f6.2,f6.1,2f7.2,t112,a4        ,t134,a8,1x,a3,1x,a8,1x,2i5)'
   else if (it .eq. 2) then
      gfmt = '(a1,a5,1x,a8,1x,a1,a1,a8,1x,f4.2,f7.2,i5,f6.1,f7.2,f6.2,f6.1,3f7.2,t112,2f7.4,f7.2,t134,a8,1x,a3,1x,a8,1x,2i5)'
      bfmt = '(a1,a5,1x,a8,1x,a1,a1,a8,1x,f4.2,f7.2,i5,f6.1,f7.2,f6.2,f6.1,3f7.2,t112,a4        ,t134,a8,1x,a3,1x,a8,1x,2i5)'
   else if (it .eq. 3) then
      gfmt = '(a1,a5,1x,a8,1x,a1,a1,a8,1x,f4.2,f7.2,i5,f6.1,f7.2,f6.2,f6.1,4f7.2,t112,2f7.4,f7.2,t134,a8,1x,a3,1x,a8,1x,2i5)'
      bfmt = '(a1,a5,1x,a8,1x,a1,a1,a8,1x,f4.2,f7.2,i5,f6.1,f7.2,f6.2,f6.1,4f7.2,t112,a4        ,t134,a8,1x,a3,1x,a8,1x,2i5)'
   else if (it .eq. 4) then
      gfmt = '(a1,a5,1x,a8,1x,a1,a1,a8,1x,f4.2,f7.2,i5,f6.1,f7.2,f6.2,f6.1,5f7.2,t112,2f7.4,f7.2,t134,a8,1x,a3,1x,a8,1x,2i5)'
      bfmt = '(a1,a5,1x,a8,1x,a1,a1,a8,1x,f4.2,f7.2,i5,f6.1,f7.2,f6.2,f6.1,5f7.2,t112,a4        ,t134,a8,1x,a3,1x,a8,1x,2i5)'
   end if
   
   ! Some information about event depths
   n_constrained = 0
   n_input = 0
   cluster_default_depth = -99.
   do iev = 1, nev
      if (depset_pr(iev) .ne. 'c' .and. depset_pr(iev) .ne. 'u' .and. depset_pr(iev) .ne. ' ') then
         n_constrained = n_constrained + 1
         x_dep_constrained(n_constrained) = depthp(iev,it1)
      end if
      if (depset_pr(iev) .eq. 'c') then
         cluster_default_depth = depthp(iev,it1)
      end if
      x_dep_input(iev) = depth_inp(iev,1)        
   end do
   if (cluster_default_depth .gt. 0.) write (io_depth_phase,'(a,f6.1)') 'Cluster default depth: ', cluster_default_depth
   if (n_constrained .ge. 3) then
      call mdian1 (x_dep_constrained, n_constrained, xmed)
      write (io_depth_phase,'(a,f6.1)') 'Median of constrained depths: ', xmed
      median_constrained_depths = xmed
   end if
   if (nev .ge. 3) then
      call mdian1 (x_dep_input, nev, xmed)
      write (io_depth_phase,'(a,f6.1)') 'Median of input file depths: ', xmed
   end if
   
   ! Good phase data for each event.
   if (direct_cal) then
      write (io_dcal,'(a,t165,a)')&
       ' STA   NETWORK    PHASE    READ  DELTA AZIM  RAY     WGT   STA  P RESIDUALS FOR ITERATION #',' '
      write (io_dcal,'(a,t113,a)')&
       ' CODE                       ERR             PARAM         CORR  *INP*   *0*    *1*    *2*    *3*    *4*',&
       ' DTMPH  DTMPC    ECI AUTHOR   CHA PHASE0     MNF IDIF'
   end if
   if (oldr_out) then
      write (io_oldr,'(a,t165,a)')&
       ' STA   NETWORK    PHASE    READ  DELTA AZIM  RAY     WGT   STA  P RESIDUALS FOR ITERATION #',' '
      write (io_oldr,'(a,t113,a)')&
       ' CODE                       ERR             PARAM         CORR  *INP*   *0*    *1*    *2*    *3*    *4*',&
       ' DTMPH  DTMPC    ECI AUTHOR   CHA PHASE0     MNF IDIF'
   end if
   write (io_depth_phase,'(a,t165,a)')&
    ' STA   NETWORK    PHASE    READ  DELTA AZIM  RAY     WGT   STA  P RESIDUALS FOR ITERATION #',' '
   write (io_depth_phase,'(a,t113,a)')&
    ' CODE                       ERR             PARAM         CORR  *INP*   *0*    *1*    *2*    *3*    *4*',&
    ' DTMPH  DTMPC    ECI AUTHOR   CHA PHASE0     MNF IDIF'
   
   do iev = 1,nev ! Loop over events
         
      ! Information on how depths were set
      if (depset_pr(iev) .eq. 'c') then
         depth_set1 = ', from cluster default depth'
      else if (depset_pr(iev) .eq. 'd') then
         depth_set1 = ', from depth phases'
      else if (depset_pr(iev) .eq. 'e') then
         depth_set1 = ', from engineering information'
      else if (depset_pr(iev) .eq. 'f') then
         depth_set1 = ', from fault modeling (e.g., InSAR, GPS)'
      else if (depset_pr(iev) .eq. 'l') then
         depth_set1 = ', from local-distance readings'
      else if (depset_pr(iev) .eq. 'm') then
         depth_set1 = ', from mloc solution with free depth'
      else if (depset_pr(iev) .eq. 'n') then
         depth_set1 = ', from near-source readings'
      else if (depset_pr(iev) .eq. 'r') then
         depth_set1 = ', from relocation with free depth outside mloc'
      else if (depset_pr(iev) .eq. 'u') then
         depth_set1 = ', unconstrained'
      else if (depset_pr(iev) .eq. 'w') then
         depth_set1 = ', from waveform analysis'
      else
         depth_set1 = ', from unknown source with code '//depset_pr(iev)
      end if
      if (hdepthshift .lt. 0.01) then
         depth_set2 = ' '
      else
         depth_set2 = ', perturbed in command file'
      end if
      if (mindx(iev,3) .gt. 0) then
         depth_free_fixed = ', free'
         depth_set = trim(depth_free_fixed)
      else
         depth_free_fixed = ', fixed'
         depth_set = trim(depth_free_fixed)//trim(depth_set1)//trim(depth_set2)
      end if
   
      write (io_out,'(2a,t165,a)')&
       ' **************************************************************************************',&   
       '*****************************************************************************', '*'
      
      write (io_out,'(a,i3,t30,a30,5x,a9,t165,a)') ' CLUSTER EVENT ', iev, evtnam(iev), 'GOOD DATA',' '
      write (io_depth_phase,'(/a,i3,t30,a30,t65,a,t80,a,f5.1,a,t165,a)') ' CLUSTER EVENT ', iev, evtnam(iev),&
       'GOOD DATA','Depth = ', depthp(iev,it1), depth_set, ' '
      if (direct_cal) write (io_dcal,'(/a,i3,t30,a30,t65,a,t80,a,f5.1,a,t165,a)') ' CLUSTER EVENT ', iev, evtnam(iev),&
       'GOOD DATA','Depth = ', depthp(iev,it1), depth_set, ' '
      if (oldr_out) write (io_oldr,'(/a,i3,t30,a30,t65,a,t80,a,f5.1,a,t165,a)') ' CLUSTER EVENT ', iev, evtnam(iev),&
       'GOOD DATA','Depth = ', depthp(iev,it1), depth_set, ' '
      
      write (io_out,'(1x,a,t165,a)') trim(infile(iev)), ' ' ! Data file
      
      ! Start location
      lon_out = lon_inp(iev)
      call set_longitude_range (lon_out, longitude_range)         
      write (io_out,'(1x,a6,i4,4(i3),f5.1,f8.3,f9.3,f6.1,f4.1,t165,a)') 'Input ', iyre (iev), mone(iev), idye(iev),&
       hour_inp(iev), min_inp(iev), sec_inp(iev), lat_inp(iev), lon_out, depth_inp(iev,1), rmag(iev), ' '
      
      ! Final location, unshifted
      lon_out = lonp(iev,it1)
      call set_longitude_range (lon_out, longitude_range)
      write (io_out,'(1x,a6,i4,4(i3),f5.1,f8.3,f9.3,f6.1,f4.1,t165,a)') 'Final ', iyre (iev), mone(iev), idye(iev),&
       hourp(iev,it1), minp(iev,it1), secp(iev,it1), latp(iev,it1), lon_out, depthp(iev,it1), rmag(iev), ' '
      
      ! Header line
      write (io_out,'(t165,a)') ' '
      write (io_out,'(a,t165,a)')&
       ' STA   NETWORK    PHASE    READ  DELTA AZIM  RAY     WGT   STA  P RESIDUALS FOR ITERATION #',' '
      write (io_out,'(a,t113,a)')&
       ' CODE                       ERR             PARAM         CORR  *INP*   *0*    *1*    *2*    *3*    *4*',&
       ' DTMPH  DTMPC    ECI AUTHOR   CHA PHASE0     MNF IDIF'
      write (io_out,'(t165,a)') ' '
      
      ! Index table on delta for sorting
      do i=1,nst(iev)
         deltiev(i) = delt(iev,i)
      end do
      call indexx (nst(iev), deltiev, indx)
      
      h_mult = depthp(iev,it1)*3. ! Multiple of focal depth, used to select near-source readings for depth control
      
      do ii = 1,nst(iev) ! Loop over readings
      
         i = indx(ii)
         deltkm = delt(iev,i)*111. ! Used to select near-source readings for depth control 
         
         if (.not.fltrhres(iev,i) .and. .not.fltrhflag(iev,i) .and. .not.fltrcres(iev,i) .and. .not.fltrcflag(iev,i)) then ! Not filtered
            do j = 0,it
               dts(j) = dt(iev,i,j) - s(iev,i,j)
            end do
            
            if (it .eq. 0) then ! Forward problem
            
               ! All readings written to phase_data file
               write (io_out,gfmt) ' ', stname(iev,i), deployment(iev,i), fcode(iev,i), no_phid(iev,i), phase(iev,i),&
                sdread(iev,i), delt(iev,i), nint(azes(iev,i)), psd(iev,i), weight(iev,i), stacs(iev,i), resisc(iev,i), dts(0),&
                readsrc(iev,i), channel(iev,i), phase0(iev,i), mnf_line(iev,i), diff_line(iev,i)
                
               ! Near-source readings written to the depth_phases file
               if (deltkm .le. max(h_mult,100.)) then
                  write (io_depth_phase,gfmt) ' ', stname(iev,i), deployment(iev,i), fcode(iev,i), no_phid(iev,i), phase(iev,i),&
                   sdread(iev,i), delt(iev,i), nint(azes(iev,i)), psd(iev,i), weight(iev,i), stacs(iev,i), resisc(iev,i),&
                   dts(0), readsrc(iev,i), channel(iev,i), phase0(iev,i), mnf_line(iev,i), diff_line(iev,i)
               end if

               ! Depth phases written to the depth_phases file                   
               if (phase(iev,i) .eq. 'pP      ' .or. phase(iev,i) .eq. 'sP      ' .or. phase(iev,i) .eq. 'pwP     ') then
                  if (bdp(stname(iev,i))) then ! Check for stations known to report bad depth phases
                     bdp_flag = '*'
                  else
                     bdp_flag = ' '
                  end if
                  if (rdp_res(iev,i,0) .gt. -100.) then
                     rdp_phase = phase(iev,i)(1:2)//'-P    '
                     write (io_depth_phase,gfmt) bdp_flag, stname(iev,i), deployment(iev,i), fcode(iev,i), no_phid(iev,i),&
                      rdp_phase, sdread(iev,i), delt(iev,i), nint(azes(iev,i)), psd(iev,i), weight(iev,i), stacs(iev,i),&
                      resisc(iev,i), rdp_res(iev,i,0), readsrc(iev,i), channel(iev,i), phase0(iev,i), mnf_line(iev,i),&
                      diff_line(iev,i)
                  else
                     rdp_phase = phase(iev,i)
                     write (io_depth_phase,gfmt) bdp_flag, stname(iev,i), deployment(iev,i), fcode(iev,i), no_phid(iev,i),&
                      rdp_phase, sdread(iev,i), delt(iev,i), nint(azes(iev,i)), psd(iev,i), weight(iev,i), stacs(iev,i),&
                      resisc(iev,i), dts(0), readsrc(iev,i), channel(iev,i), phase0(iev,i), mnf_line(iev,i), diff_line(iev,i)
                  end if
               end if
                
               ! Direct calibration readings written to dcal_phase_data file
               if (direct_cal .and. .not.fltrh(iev,i)) write (io_dcal,gfmt) ' ', stname(iev,i), deployment(iev,i), fcode(iev,i),&
                no_phid(iev,i), phase(iev,i), sdread(iev,i), delt(iev,i), nint(azes(iev,i)), psd(iev,i), weight(iev,i),&
                stacs(iev,i), resisc(iev,i), dts(0), readsrc(iev,i), channel(iev,i), phase0(iev,i), mnf_line(iev,i),&
                diff_line(iev,i)
                
               ! Output for limited distance range file
               if (oldr_out .and. delt(iev,i) .ge. dist1 .and. delt(iev,i) .le. dist2) write (io_oldr,gfmt) ' ',&
                stname(iev,i), deployment(iev,i), fcode(iev,i), no_phid(iev,i), phase(iev,i), sdread(iev,i), delt(iev,i),&
                nint(azes(iev,i)), psd(iev,i), weight(iev,i), stacs(iev,i), resisc(iev,i), dts(0), readsrc(iev,i),&
                channel(iev,i), phase0(iev,i), mnf_line(iev,i), diff_line(iev,i)

            else ! Iterated results
            
               ! All readings written to phase_data file
               write (io_out,gfmt) ' ', stname(iev,i), deployment(iev,i), fcode(iev,i), no_phid(iev,i), phase(iev,i),&
                sdread(iev,i), delt(iev,i), nint(azes(iev,i)), psd(iev,i), weight(iev,i), stacs(iev,i), resisc(iev,i),&
                (dts(j), j=0,it), dtmph(iev,i), dtmpc(iev,i), eci(iev,i), readsrc(iev,i), channel(iev,i), phase0(iev,i),&
                mnf_line(iev,i), diff_line(iev,i)
                
               ! Near-source readings written to the depth_phases file
               if (deltkm .le. max(h_mult,100.)) then
                  write (io_depth_phase,gfmt) ' ', stname(iev,i), deployment(iev,i), fcode(iev,i), no_phid(iev,i),&
                   phase(iev,i), sdread(iev,i), delt(iev,i), nint(azes(iev,i)), psd(iev,i), weight(iev,i), stacs(iev,i),&
                   resisc(iev,i), (dts(j), j=0,it), dtmph(iev,i), dtmpc(iev,i), eci(iev,i), readsrc(iev,i), channel(iev,i),&
                   phase0(iev,i), mnf_line(iev,i), diff_line(iev,i)
               end if
                
               ! Depth phases written to the depth_phases file                   
               if (phase(iev,i) .eq. 'pP      ' .or. phase(iev,i) .eq. 'sP      ' .or. phase(iev,i) .eq. 'pwP     ') then
                  if (bdp(stname(iev,i))) then ! Check for stations known to report bad depth phases
                     bdp_flag = '*'
                  else
                     bdp_flag = ' '
                  end if
                  if (rdp_res(iev,i,0) .gt. -100.) then
                     rdp_phase = phase(iev,i)(1:2)//'-P    '
                     write (io_depth_phase,gfmt) bdp_flag, stname(iev,i), deployment(iev,i), fcode(iev,i), no_phid(iev,i),&
                      rdp_phase, sdread(iev,i), delt(iev,i), nint(azes(iev,i)), psd(iev,i), weight(iev,i), stacs(iev,i),&
                      resisc(iev,i), (rdp_res(iev,i,j), j=0,it), dtmph(iev,i), dtmpc(iev,i), eci(iev,i), readsrc(iev,i),&
                      channel(iev,i), phase0(iev,i), mnf_line(iev,i), diff_line(iev,i)
                  else
                     rdp_phase = phase(iev,i)
                     write (io_depth_phase,gfmt) bdp_flag, stname(iev,i), deployment(iev,i), fcode(iev,i), no_phid(iev,i),&
                      rdp_phase, sdread(iev,i), delt(iev,i), nint(azes(iev,i)), psd(iev,i), weight(iev,i), stacs(iev,i),&
                      resisc(iev,i), (dts(j), j=0,it), dtmph(iev,i), dtmpc(iev,i), eci(iev,i), readsrc(iev,i), channel(iev,i),&
                      phase0(iev,i), mnf_line(iev,i), diff_line(iev,i)
                  end if
               end if
                
               ! Direct calibration readings written to dcal_phase_data file
               if (direct_cal .and. .not.fltrh(iev,i)) write (io_dcal,gfmt) ' ', stname(iev,i), deployment(iev,i), fcode(iev,i),&
                no_phid(iev,i), phase(iev,i), sdread(iev,i), delt(iev,i), nint(azes(iev,i)), psd(iev,i), weight(iev,i),&
                stacs(iev,i), resisc(iev,i), (dts(j), j=0,it), dtmph(iev,i), dtmpc(iev,i), eci(iev,i), readsrc(iev,i),&
                channel(iev,i), phase0(iev,i), mnf_line(iev,i), diff_line(iev,i)
  
               ! Output for limited distance range file
               if (oldr_out .and. delt(iev,i) .ge. dist1 .and. delt(iev,i) .le. dist2) write (io_oldr,gfmt) ' ',&
                stname(iev,i), deployment(iev,i), fcode(iev,i), no_phid(iev,i), phase(iev,i), sdread(iev,i),&
                delt(iev,i), nint(azes(iev,i)), psd(iev,i), weight(iev,i), stacs(iev,i), resisc(iev,i), (dts(j), j=0,it),&
                dtmph(iev,i), dtmpc(iev,i), eci(iev,i), readsrc(iev,i), channel(iev,i), phase0(iev,i), mnf_line(iev,i),&
                diff_line(iev,i)

               ! Output of readings with large cluster residuals
               if (lresout) then
                  if (abs(eci(iev,i)) .gt. lres) then
                     write (io_lres,'(i3,1x,a20,1x,a8,1x,a8,i5,1x,f8.2,1x,a,2x,a,1x,i6)') iev, sad(iev,i),&
                      readsrc(iev,i), phase0 (iev,i), mnf_line(iev,i), eci(iev,i), infile20(iev), phase(iev,i),&
                      diff_line(iev,i)
                  end if
               end if
               
            end if
         end if
         
         ! Fill in array of residuals for calculation of reading errors
         if (connected(iev,i)) then
            qtest = stname(iev,i)//deployment(iev,i)//phase(iev,i)
            do j = 1,nqc
               if (qtest .eq. qname1(j) .and. idiff(iev,i) .eq. 0 .and. indexq(j) .lt. n_qres_max) then
                  indexq(j) = indexq(j) + 1
                  if (calibration) then ! Indirect calibration
                     qres(j,indexq(j)) = dt_cal(iev,i) - s_cal(iev,i)
!                        print *, iev, i, qtest, j, indexq(j), qres(j,indexq(j)), dt_cal(iev,i), s_cal(iev,i)
                  else ! Direct calibration or uncalibrated
                     qres(j,indexq(j)) = dts(it)
                  end if
                  iqiev(j,indexq(j)) = iev ! array of event numbers association with residuals, needed for tomo3
                  exit
               end if
            end do
         end if
     
      end do ! End of loop over readings
      write (io_out,'(t165,a)') ' '
      write (io_out,'(t165,a)') ' '
      
   end do ! End of loop over events
   
   ! Bad phase data for each event
   
   do iev = 1,nev ! Loop over events
   
      ! Information on how depths were set
      if (depset_pr(iev) .eq. 'c') then
         depth_set1 = ', from cluster default depth'
      else if (depset_pr(iev) .eq. 'd') then
         depth_set1 = ', from depth phases'
      else if (depset_pr(iev) .eq. 'e') then
         depth_set1 = ', from engineering information'
      else if (depset_pr(iev) .eq. 'f') then
         depth_set1 = ', from fault modeling (e.g., InSAR, GPS)'
      else if (depset_pr(iev) .eq. 'l') then
         depth_set1 = ', from local-distance readings'
      else if (depset_pr(iev) .eq. 'm') then
         depth_set1 = ', from mloc solution with free depth'
      else if (depset_pr(iev) .eq. 'n') then
         depth_set1 = ', from near-source readings'
      else if (depset_pr(iev) .eq. 'r') then
         depth_set1 = ', from relocation with free depth outside mloc'
      else if (depset_pr(iev) .eq. 'u') then
         depth_set1 = ', unconstrained'
      else if (depset_pr(iev) .eq. 'w') then
         depth_set1 = ', from waveform analysis'
      else
         depth_set1 = ', from unknown source with code '//depset_pr(iev)
      end if
      if (hdepthshift .lt. 0.01) then
         depth_set2 = ' '
      else
         depth_set2 = ', perturbed in command file'
      end if

      write (io_out,'(2a,t165,a)')&
       ' **************************************************************************************',&   
       '*****************************************************************************', '*'
      write (io_out,'(a,i3,t30,a30,5x,a8,t165,a)') ' CLUSTER EVENT ', iev, evtnam(iev), 'BAD DATA', ' '
      write (io_depth_phase,'(/a,i3,t30,a30,t65,a,t80,a,f5.1,2a,t165,a)') ' CLUSTER EVENT ', iev, evtnam(iev),&
       'BAD DATA','Depth = ', depthp(iev,it1), trim(depth_set1), trim(depth_set2), ' '
      if (direct_cal) write (io_dcal,'(/a,i3,t30,a30,t65,a,t80,a,f5.1,2a,t165,a)') ' CLUSTER EVENT ', iev, evtnam(iev),&
       'BAD DATA','Depth = ', depthp(iev,it1), trim(depth_set1), trim(depth_set2), ' '
      if (oldr_out) write (io_oldr,'(/a,i3,t30,a30,t65,a,t80,a,f5.1,2a,t165,a)') ' CLUSTER EVENT ', iev, evtnam(iev),&
       'BAD DATA','Depth = ', depthp(iev,it1), trim(depth_set1), trim(depth_set2), ' '
      write (io_out,'(1x,a,t165,a)') infile(iev), ' ' ! Data file
      
      ! Start location         
      lon_out = lon_inp(iev)
      call set_longitude_range (lon_out, longitude_range)
      write (io_out,'(1x,a6,i4,4(i3),f5.1,f8.3,f9.3,f6.1,f4.1,t165,a)') 'Input ', iyre (iev), mone(iev), idye(iev),&
       hour_inp(iev), min_inp(iev), sec_inp(iev), lat_inp(iev), lon_out, depth_inp(iev,1), rmag(iev), ' '
       
      ! Final location, unshifted
      lon_out = lonp(iev,it1)
      call set_longitude_range (lon_out, longitude_range)
      write (io_out,'(1x,a6,i4,4(i3),f5.1,f8.3,f9.3,f6.1,f4.1,t165,a)') 'Final ', iyre (iev), mone(iev), idye(iev),&
       hourp(iev,it1), minp(iev,it1), secp(iev,it1), latp(iev,it1), lon_out, depthp(iev,it1), rmag(iev), ' '
      
      ! Missing stations
      if (nmiss(iev) .gt. 0) then
         write (io_out,'(t165,a)') ' '
         write (io_out,'(a,t165,a)') ' Missing stations:', ' '
         do i = 1,nmiss(iev) 
            write (io_out,'(5x,a,t165,a)') missta(iev,i), ' '
         end do
      end if
      
      ! Header line 
      write (io_out,'(t165,a)') ' '
      write (io_out,'(2a,t165,a/2a,t134,a)') ' STA   NETWORK    PHASE    READ  DELTA AZIM  RAY     WGT   STA ',&
       ' P RESIDUALS FOR ITERATION #                    WHY', ' ',&
       ' CODE                       ERR             PARAM         CORR ',&
       '*INP*    *0*    *1*    *2*    *3*    *4*        BAD', 'AUTHOR   CHA PHASE0     MNF IDIF'
      write (io_out,'(t165,a)') ' '
      
      ! Index table on delta for sorting
      do i = 1,nst(iev)
         deltiev(i) = delt(iev,i)
      end do
      call indexx (nst(iev), deltiev, indx)
      
      do ii=1,nst(iev) ! Loop over readings
         i = indx(ii)
         
         if (fcode(iev,i) .eq. 'm') cycle ! Don't print readings from missing stations
         
         ! Reasons for badness
         ! Only consider if a reading has been flagged or has a large residual
         if (fltrhres(iev,i) .or. fltrhflag(iev,i) .or. fltrcres(iev,i) .or. fltrcflag(iev,i)) then
            do j = 0,it
               dts(j) = dt(iev,i,j) - s(iev,i,j)
            end do
            
            if (fltrhres(iev,i) .or. fltrcres(iev,i)) then ! Large residual
               reason = rpres
!                  if (.not.fltrhflag(iev,i) .and. .not.fltrcflag(iev,i) .and. delt(iev,i) .gt. windloclim) then ! Not already flagged
               if (.not.fltrhflag(iev,i) .and. .not.fltrcflag(iev,i)) then ! Not already flagged
                  write (io_xdat,'(i3,1x,a6,1x,a8,1x,a8,1x,i5,2(1x,a))') iev, stname(iev,i),&
                   readsrc(iev,i), phase0(iev,i), mnf_line(iev,i), 'x', infile20(iev)
               end if
            end if
            
            if (fltrhflag(iev,i) .or. fltrcflag(iev,i)) then ! Flagged for some reason
               dtslim = 2.0*ttsprd(iev,i)
               if (abs(dts(it)-ttoff(iev,i)) .le. dtslim .and. (fcode(iev,i) .eq. 'x' .or. fcode(iev,i) .eq. 's')) then
                  reason = '?   ' ! Maybe investigate why small residual is flagged
               else
                  reason = '    ' ! Leave blank for flagged phases
               end if
            end if
                           
            write (io_out,bfmt) ' ', stname(iev,i), deployment(iev,i), fcode(iev,i), no_phid(iev,i), phase(iev,i),&
             sdread(iev,i), delt(iev,i), nint(azes(iev,i)), psd(iev,i), weight(iev,i), stacs(iev,i), resisc(iev,i),&
             (dts(j), j=0,it), reason, readsrc(iev,i), channel(iev,i), phase0(iev,i), mnf_line(iev,i), diff_line(iev,i)
             
            ! Near-source readings written to the depth_phases file
            if (deltkm .le. max(h_mult,100.)) then
               write (io_depth_phase,bfmt) ' ', stname(iev,i), deployment(iev,i), fcode(iev,i), no_phid(iev,i), phase(iev,i),&
                sdread(iev,i), delt(iev,i), nint(azes(iev,i)), psd(iev,i), weight(iev,i), stacs(iev,i), resisc(iev,i),&
                (dts(j), j=0,it), reason, readsrc(iev,i), channel(iev,i), phase0(iev,i), mnf_line(iev,i), diff_line(iev,i)
            end if
             
            ! Depth phases written to the depth_phases file                   
            if (phase(iev,i) .eq. 'pP      ' .or. phase(iev,i) .eq. 'sP      ' .or. phase(iev,i) .eq. 'pwP     ') then
               if (bdp(stname(iev,i))) then ! Check for stations known to report bad depth phases
                  bdp_flag = '*'
               else
                  bdp_flag = ' '
               end if
               if (rdp_res(iev,i,0) .gt. -100.) then
                  rdp_phase = phase(iev,i)(1:2)//'-P    '
                  write (io_depth_phase,bfmt) bdp_flag, stname(iev,i), deployment(iev,i), fcode(iev,i), no_phid(iev,i),&
                   rdp_phase, sdread(iev,i), delt(iev,i), nint(azes(iev,i)), psd(iev,i), weight(iev,i), stacs(iev,i),&
                   resisc(iev,i), (rdp_res(iev,i,j), j=0,it), reason, readsrc(iev,i), channel(iev,i), phase0(iev,i),&
                   mnf_line(iev,i), diff_line(iev,i)
               else
                  rdp_phase = phase(iev,i)
                  write (io_depth_phase,bfmt) bdp_flag, stname(iev,i), deployment(iev,i), fcode(iev,i), no_phid(iev,i),&
                   rdp_phase, sdread(iev,i), delt(iev,i), nint(azes(iev,i)), psd(iev,i), weight(iev,i), stacs(iev,i),&
                   resisc(iev,i), (dts(j), j=0,it), reason, readsrc(iev,i), channel(iev,i), phase0(iev,i), mnf_line(iev,i),&
                   diff_line(iev,i)
               end if
            end if
             
            ! Direct calibration readings written to dcal_phase_data file
            if (direct_cal .and. .not.fltrhdelt(iev,i) .and. (fltrhres(iev,i) .or. fltrhflag(iev,i)))&
             write (io_dcal,bfmt) ' ', stname(iev,i), deployment(iev,i), fcode(iev,i), no_phid(iev,i),&
             phase(iev,i), sdread(iev,i), delt(iev,i), nint(azes(iev,i)), psd(iev,i), weight(iev,i), stacs(iev,i),&
             resisc(iev,i), (dts(j), j=0,it), reason, readsrc(iev,i), channel(iev,i), phase0(iev,i), mnf_line(iev,i),&
             diff_line(iev,i)
             
            ! Output for limited distance range file
            if (oldr_out .and. delt(iev,i) .ge. dist1 .and. delt(iev,i) .le. dist2) write (io_oldr,bfmt) ' ',&
             stname(iev,i), deployment(iev,i), fcode(iev,i), no_phid(iev,i), phase(iev,i), sdread(iev,i), delt(iev,i),&
             nint(azes(iev,i)), psd(iev,i), weight(iev,i), stacs(iev,i), resisc(iev,i), (dts(j), j=0,it), reason,&
             readsrc(iev,i), channel(iev,i), phase0(iev,i), mnf_line(iev,i), diff_line(iev,i)
            
         end if
      end do ! End of loop over readings
      write (io_out,'(t165,a)') ' '
      write (io_out,'(t165,a)') ' '
      
   end do ! End of loop over events
   
   close (io_out)
   close (io_xdat)
   close (io_depth_phase)
   if (lresout) close (io_lres)
   if (direct_cal) close (io_dcal)
   if (oldr_out) close (io_oldr)
   
   return
   
end subroutine mlocout_phase_data
