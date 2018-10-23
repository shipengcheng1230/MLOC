!***********************************************************************************************************************************
subroutine proc_anno (iev, command_help)

! Annotations, which will be printed at the end of lines in HDF files.
! Limited to 20 characters

   implicit none

   include 'mloc.inc'

   integer :: iev
   logical :: command_help
   character(len=132) :: msg

   if (command_help) then
      write (*,'(4(a/))')&
      '"anno" (annotation) is used to provide some extra textual information about an',&
      'event. The string will be printed at the end of the corresponding line in the',&
      'HDF files. Maximum 20 characters. This annotation over-writes an annotation given in',&
      'the event record of an MNF file, in which case a warning will be given.'
      return
   end if

   annotation(iev) = params(1:20)

   if (verbose_screen) then
      write (msg,'(a,i3,1x,a)') 'proc_anno: event ', iev, trim(annotation(iev))
      call fyi (trim(msg))
   end if

   return
   
end subroutine proc_anno


!***********************************************************************************************************************************
subroutine proc_bdps (command_help)

! File of stations that are suspected off reporting bogus depth phases

   implicit none

   include 'mloc.inc'

   logical :: ex, command_help

   if (command_help) then
      write (*,'(8(a/))')&
       '"bdps" (bogus depth phase stations) specifies the pathname of a file containing',&
       'a list of stations that are suspected of reporting bogus depth phases, i.e., depth phase',&
       'arrival times that are generated from theoretical arrivals relative to a preliminary',&
       'hypocenter (e.g., the PDE). Depth phase readings from listed stations will be plotted',&
       'at a smaller size in "tt6" plots and will have an asterisk next to their entries in',&
       'the ".depth_phases" file. The file format for the list can be anything, as long as the',&
       'station code is left-justified in columns 2:6. No additional data, e.g., coordinates or',&
       'station name, are required but they are useful for identification.'
      return
   end if

   if (params(1:1) .eq. ' ') then
      call warnings ('proc_bdps: "bdps" requires a pathname as argument')
   else
      inquire (file=trim(params),exist=ex)
      if (ex) then
         bdp_list = .true.
         bdp_filnam = trim(params)
         if (verbose_screen) call fyi ('proc_bdps: bogus depth phase station file: '//trim(bdp_filnam))
      else
         call warnings ('proc_bdps: file '//trim(params)//' does not exist')
      end if
   end if

   return
   
end subroutine proc_bdps


!***********************************************************************************************************************************
subroutine proc_bias (command_help)

! Hypocentroid bias correction (see Eq. 86 in Jordan & Sverdrup, 1981)
! Can be toggled by using the command without an argument, or set explicitly with arguments "on" or "off".

   implicit none

   include 'mloc.inc'

   logical :: command_help

   if (command_help) then
      write (*,'(7(a/))')&
       '"bias" (bias correction) determines if a correction will be made for a minor source',&
       'of bias in the hypocentroid, related to the fact that some events have more',&
       'readings than others. See Jordan & Sverdrup for details. By default, or if the',&
       'argument is blank or "on", the correction is made. If the argument is "off" the',&
       'correction is not made. Correcting the bias increases the variance of the',&
       'hypocentroid slightly. Compared to the main source of bias in the hypocentroid,',&
       'i.e., travel time model inadequacies, this is really negligible.'
      write (*,'(a,l1/)') 'Current value: ', bias_corr
      return
   end if

   if (params(1:1) .eq. ' ') then
      bias_corr = .not.bias_corr
   else if (params(1:2) .eq. 'on') then
      bias_corr = .true.
   else if (params(1:3) .eq. 'off') then
      bias_corr = .false.
   else
      call warnings ('proc_bias: illegal argument')
   end if

   if (verbose_screen) then
      if (bias_corr) then
         call fyi ('proc_bias: hypocentroid bias correction is on')
      else
         call fyi ('proc_bias: hypocentroid bias correction is off')
      end if
   end if

return

end subroutine proc_bias

!***********************************************************************************************************************************
subroutine proc_bloc (command_help)

! BAYESLOC format output file.
! by Ezgi Karasözen 3/9/2016

   implicit none

   include 'mloc.inc'

   logical :: command_help

   if (command_help) then
      write (*,'(4(a/))')&
       '"bloc" (BayesLoc formatted output file) specifies that an output file (with suffix',&
       '".bloc") will be written. The format is designed for easy import into BayesLoc.',&
       'Default is "off". Issuing the command with no argument toggles the current state, or the',&
       'arguments "on" or "off" can be used to explicitly set the state.'
      return
   end if

   if (params(1:1) .eq. ' ') then
      blocout = .not.blocout
   else if (params(1:2) .eq. 'on') then
      blocout = .true.
   else if (params(1:3) .eq. 'off') then
      blocout = .false.
   else
      call warnings ('proc_bloc: illegal argument')
   end if

   if (verbose_screen) then
      if (blocout) then
         call fyi ('proc_bloc: will create a BAYESLOC file')
      else
         call fyi ('proc_bloc: will not create a BAYESLOC file')
      end if
   end if

   return
   
end subroutine proc_bloc
!***********************************************************************************************************************************
subroutine proc_bptc (command_help)

! Bounce point topography correction for pP, sP and pwP.
! No correction is made for zero focal depth, regardless of topography.
! Can be toggled by using the command without an argument, or set explicitly with arguments "on" or "off".

   implicit none

   include 'mloc.inc'

   logical :: command_help

   if (command_help) then
      write (*,'(9(a/))')&
       '"bptc" (bounce point topography correction) makes a correction to theoretical',&
       'travel times of pP and sP for the topography at the bounce point. For areas',&
       'above sea level the times of pP and sP are increased. For oceanic areas the travel',&
       'times of pP and sP are reduced because the reflection point (seafloor) is below',&
       'sea level. If the water depth is greater than 1.5 km the travel time of the pwP',&
       'phase is calculated by adding a term for the propagation time in the water column',&
       'to the pP time. No calculation is done for swP. The command can be toggled on and',&
       'off by issuing the command without an argument, or set explicitly with arguments',&
       '"on" and "off"'
      write (*,'(a,l1/)') 'Current value: ', bptc
      return
   end if

   if (params(1:1) .eq. ' ') then
      bptc = .not.bptc
   else if (params(1:2) .eq. 'on') then
      bptc = .true.
   else if (params(1:3) .eq. 'off') then
      bptc = .false.
   else
      call warnings ('proc_bptc: illegal argument')
   end if

   if (verbose_screen) then
      if (bptc) then
         call fyi ('proc_bptc: bounce point topography correction is on')
      else
         call fyi ('proc_bptc: bounce point topography correction is off')
      end if
   end if

return

end subroutine proc_bptc


!***********************************************************************************************************************************
subroutine proc_cal (iev, comd, command_help)

! Calibration event

   implicit none

   include 'mloc.inc'

   integer :: iev
   real :: t11, t12, t22
   logical :: command_help
   character(len=4) :: comd
   character(len=132) :: msg

   if (command_help) then
      write (*,'(14(a/))')&
       '"calX" (calibration event) is used to declare a calibration location for the current',&
       'event. Epicenter is always assumed to be calibrated. The 4th character in the command',&
       'is used to specify which of the other hypocentral parameters are to be considered',&
       'calibrated:',&
       '  e : the epicenter alone is calibrated (e.g., from InSAR)',&
       '  f : focal depth is also calibrated, but not origin time',&
       '  t : origin time is constrained by local-distance data but it cannot be considered',&
       '      to be calibrated because focal depth is not calibrated',&
       '  h : epicenter, depth, and origin time are calibrated',&
       'The command takes 11 arguments: hour, minutes, and seconds of the origin time, plus',&
       'latitude, longitude and depth, in that order, followed by 3 arguments giving the',&
       '90% confidence ellipse and 2 arguments for the uncertainties in depth and OT.',&
       'The sequence for confidence ellipse is: azimuth of semi-minor axis (clockwise from north)',&
       'semi-minor axis length, semi-major axis length (both in km).'
      return
   end if

   call decode (params, value, nvar)
   
   if (nvar .eq. 11) then
      write (msg,'(a,i3,1x,2a)') 'proc_cal: event ', iev, trim(evtnam(iev)), ' is declared as a calibration event'
      call fyi (trim(msg))
      cal_par(iev) = comd(4:4)
      if (cal_par(iev) .eq. 'i') then ! Old form of the command
         call warnings ('proc_cal: the form "cali" is no longer supported; converted to "cale"')
         cal_par(iev) = 'e'
      end if
      if (cal_par(iev) .ne. 'e' .and. cal_par(iev) .ne. 'f' .and. cal_par(iev) .ne. 't' .and. cal_par(iev) .ne. 'h') then
         call warnings ('proc_cal: command "'//comd//'" not recognized')
         cal_par(iev) = ' '
      end if
      if (cal_par(iev) .eq. 'h') ot_cal = .true. ! Origin time for the cluster can be calibrated
      calibration = .true.
      cal_event(iev,1) = .true.
      ncal(1) = ncal(1) + 1
      cal_lat(iev,1) = value(4)
      cal_lon(iev,1) = value(5)
      cal_dep(iev,1) = value(6)
      cal_hr(iev,1) = int(value(1))
      cal_min(iev,1) = int(value(2))
      cal_sec(iev,1) = value(3)
      ! Check for old-style parameterization of confidence ellipse (as covariance elements)
      if (value(8) .lt. 1.0e-3) then
         write (msg,'(a,3f10.3)') 'proc_cal: Old-style parameterization of confidence ellipse? ', value(7), value(8), value(9)
         call warnings (trim(msg))
      end if
      ! Covariances
      call ell2cv (1.0, value(7), value(8), value(9), t11, t12, t22)
      rcv(iev,1,1,1) = t11 ! Latitude
      rcv(iev,1,1,2) = t12 ! Lat-Lon cross term
      rcv(iev,1,2,1) = t12 ! Lat-Lon cross term
      rcv(iev,1,2,2) = t22 ! Longitude
      rcv(iev,1,3,3) = value(10)*value(10) ! Depth
      rcv(iev,1,4,4) = value(11)*value(11) ! OT
      if (verbose_screen) then
         write (msg,'(a,3f8.3,2i3,f6.2,6f8.2)') 'proc_cali: ', cal_lat(iev,1), cal_lon(iev,1), cal_dep(iev,1), cal_hr(iev,1),&
          cal_min(iev,1), cal_sec(iev,1), rcv(iev,1,1,1), rcv(iev,1,1,2), rcv(iev,1,2,1), rcv(iev,1,2,2), rcv(iev,1,3,3),&
          rcv(iev,1,4,4)
         call fyi (trim(msg))
      end if
   else
      call warnings ('proc_cal: "cal" requires 11 arguments')
   end if

   return
   
end subroutine proc_cal


!***********************************************************************************************************************************
subroutine proc_ccat (command_help)

! COMCAT output file. This is a format designed especially for import into the USGS/NEIC "COMCAT" catalog server.

   implicit none

   include 'mloc.inc'

   logical :: command_help
   character(len=160) :: command_line

   if (command_help) then
      write (*,'(7(a/))')&
       '"ccat" (output file for COMCAT) specifies that a .comcat file will be written.',&
       'This is a single file similar to the .datf_mnf file, but with additional information',&
       'needed for uploading to the USGS/NEIC COMCAT server. The command takes a single optional',&
       'argument, the filename (assumed to be in the data directory) of a text file containing',&
       'a commentary on the nature of the cluster and the relocation procedures and results.',&
       'The commentary file should be hard-wrapped at a reasonable line-length (72 characters',&
       'is good). The commentary will be added to the .comcat file as comment records.'
      return
   end if
   
   comcatout = .true.

   if (params(1:1) .ne. ' ') then
      commentary_fname = trim(params)
   else
      call warnings ('proc_ccat: no commentary file has been specified')
      commentary_fname = 'none'
   end if
   
   ! Create a subdirectory for ComCat files
   ccat_folder = trim(datadir)//dirsym//trim(basename)//'_comcat'
   command_line = 'mkdir '//trim(ccat_folder)
   call system (trim(command_line))

   if (verbose_screen) call fyi ('proc_ccat: will create a COMCAT file using the commentary file '//trim(commentary_fname))

   return
   
end subroutine proc_ccat


!***********************************************************************************************************************************
subroutine proc_cfil (cmndfil, command_help)

! Specify a command file

   implicit none

   include 'mloc.inc'

   logical :: cmndfil, ex, command_help
   character(len=100) :: filnam

   if (command_help) then
      write (*,'(5(a/))')&
       '"cfil" takes one argument, the name of a command file. All the commands in that',&
       'file will be processed before control is returned for interactive command',&
       'processing. Multiple command files can be invoked (by separate "cfil" commands), but',&
       'they must all be found in the data directory specified when the program starts.',&
       'The "cfil" command cannot be issued from a command file, only interactively.'
      return
   end if

   if (cmndfil) then
      call warnings ('proc_cfil: "cfil" cannot be called from a command file')
   else if (params(1:1) .ne. ' ') then
      filnam = trim(datadir)//dirsym//trim(params)
      inquire (file=filnam,exist=ex)
      if (ex) then
         cmndfil = .true.
         open (io_cfil,file=filnam,status='old')
         if (verbose_screen) call fyi ('proc_cfil: file '//trim(filnam)//' opened as command file')
      else
         call warnings ('proc_cfil: '//trim(filnam)//' does not exist')
      end if
   else
      call warnings ('proc_cfil: "cfil" requires a file name as argument')
   end if

   return
   
end subroutine proc_cfil


!***********************************************************************************************************************************
subroutine proc_clim (command_help)

! Epicentral distance limits for cluster vectors.
! Up to three distance ranges can be specified.

   implicit none

   include 'mloc.inc'

   integer :: i, j
   logical :: command_help
   character(len=132) :: msg

   if (command_help) then
      write (*,'(5(a/))')&
       '"clim" (cluster vector limits) specifies a set of up to three ranges in epicentral',&
       'distance that will be used to estimate the cluster vectors. This allows skipping',&
       'readings at short distances, in the Pdiff range, and near the caustic for PKP',&
       'phases, for example. Since each range requires two values, this command requires',&
       '2, 4, or 6 arguments.'
      write (*,'(a,3(2f6.1,2x)/)') 'Current values : ', (clim(i,1),clim(i,2),i=1,nclim)
      return
   end if

   if (params(1:1) .eq. ' ') then
      nclim = 1
      clim(1,1) = 1.0
      clim(1,2) = 180.0
   else 
      call decode (params, value, nvar)
      if (nvar .eq. 2 .or. nvar .eq. 4 .or. nvar .eq. 6) then
         nclim = nvar/2
         do i = 1,nclim
            j = 2*i
            clim(i,1) = value(j-1)
            clim(i,2) = value(j)
            if (verbose_screen) then
               write (msg,'(a,i1,a,2f6.1)') 'proc_clim: distance interval ', i, ':',clim(i,1), clim(i,2)
               call fyi (trim(msg))
            end if
         end do
      else
         call warnings ('proc_clim: "clim" requires 1-3 pairs of parameters')
      end if
   end if

   return
   
end subroutine proc_clim


!***********************************************************************************************************************************
subroutine proc_comm (command_help)

! Comment line

   implicit none

   include 'mloc.inc'

   logical :: command_help

   if (command_help) then
      write (*,'(5(a/))')&
       '"comm" (comment) declares a comment line. It can be used to disable a command in',&
       'the command file or simply to add some additional information. It was formerly used',&
       'as a convenient way to take an event out of the cluster temporarily without losing',&
       'track of it completely, but that usage has been superceded by the "kill" command for',&
       'blocks of events and the "kill" argument to the "memb" command for single events.'
      return
   end if

   if (verbose_screen) call fyi ('proc_comm: '//trim(params))

   return
   
end subroutine proc_comm


!***********************************************************************************************************************************
subroutine proc_corr (command_help)

! Station elevation correction

   implicit none

   include 'mloc.inc'

   logical :: command_help

   if (command_help) then
      write (*,'(9(a/))')&
       '"corr" (station elevation correction) controls station elevation corrections.',&
       'A single integer argument is required. The allowed arguments are:',&
       '   0 - No correction',&
       '   1 - Station elevation corrections are made',&
       'Regardless of the usage of the "corr" command, ellipticity corrections are always',&
       'applied. In the case of direct calibration the choice of station elevation correction',&
       'changes the reference plane for focal depth. If station corrections are made, the',&
       'reference plane is the geoid. If not, the reference plane is, roughly, the surface',&
       '(average elevation) of the source region.'
      write (*,'(a,i1/)') 'Current value: ', tt_corr
      return
   end if

   call decode (params, value, nvar)
   
   select case (nint(value(1)))
      case (0)
         tt_corr = 0
         if (verbose_screen) call fyi ('proc_corr: no station corrections')
      case (1)
         tt_corr = 1
         if (verbose_screen) call fyi ('proc_corr: elevation corrections are made')
      case default
         call warnings ('"corr" requires one argument (0 or 1)')
   end select
   
   return
   
end subroutine proc_corr


!***********************************************************************************************************************************
subroutine proc_cptf (command_help)

! Color palette table for plotting topography in GMT

   implicit none

   include 'mloc.inc'

   logical :: command_help, ex
   character(len=132) :: filename

   if (command_help) then
      write (*,'(3(a/))')&
       '"cptf" (color palette table file) defines a color palette table to use for plotting',&
       'topography. The argument is only the name of the desired color palette table; the path',&
       'is defined in the main program.'
      return
   end if

   if (params(1:1) .eq. ' ') then
      call warnings ('"cptf" requires one argument, the name of the color palette table file')
   else
      filename = trim(cpt_path)//dirsym//trim(params)
      inquire (file=filename,exist=ex)
      if (ex) then
         cpt_file = trim(params)
         if (verbose_screen) call fyi ('proc_cptf: color palette is '//trim(filename)) 
      else
         call warnings ('proc_cptf: file '//trim(filename)//' does not exist')
      end if
   end if

   return
   
end subroutine proc_cptf


!***********************************************************************************************************************************
subroutine proc_ctyp (command_help)

! Calibration type.

   implicit none

   include 'mloc.inc'

   logical :: command_help

   if (command_help) then
      write (*,'(8(a/))')&
      '"ctyp" specifies calibration type. This applies only to indirect calibration mode,',&
      'and only to calibration events. It determines the way in which the uncertainty',&
      'for calibration events is calculated. There are three options, specified by an',&
      'integer value:',&
      '   1 = traditional: use the uncertainty of the supplied calibration data',&
      '   2 = systematic: add calibration shift uncertainty to cluster vector uncertainty',&
      '   3 = optimal: traditional or systematic, whichever has smaller semi-major axis',&
      'The default is "systematic".'
      return
   end if

   call decode (params, value, nvar)
   
   if (nvar .eq. 1) then
      if (value(1) .ge. 1 .and. value(1) .le. 3) then
         icaltype = nint(value(1))
         select case (icaltype)
            case (1)
               icaltype_pr = 'traditional '
            case (2)
               icaltype_pr = 'systematic  '
            case (3)
               icaltype_pr = 'optimal     '
         end select
      else
         call warnings ('"ctyp": illegal argument')
      end if
   else
      call warnings ('proc_ctyp: "ctyp" requires one parameter, an integer (1-3) specifying calibration type')
   end if
   
   if (verbose_screen) call fyi ('proc_ctyp: calibration type: '//trim(icaltype_pr))

   return
   
end subroutine proc_ctyp


!***********************************************************************************************************************************
subroutine proc_cvff (command_help)

! Specification of an additional uncertainty ('fudge factor') in cluster vectors

   implicit none

   include 'mloc.inc'

   logical :: command_help
   character(len=132) :: msg

   if (command_help) then
      write (*,'(5(a/))')&
       '"cvff" specifies a radius (km) for additional uncertainty that will be added to the',&
       'cluster vector covariance matrices to account for bias from non-gaussian components',&
       'of the arrival time data sets. This is completely ad hoc (the acronym stands for',&
       'Cluster Vector Fudge Factor), but it goes in the right direction. A value of about',&
       '1.0 km has been inferred from some tests, but further testing is needed.'
      return
   end if

   call decode (params, value, nvar)
   
   if (nvar .eq. 1) then
      radius_cvff = value(1)
      if (verbose_screen) then
         write (msg,'(a,f6.1)') 'proc_cvff: fudge factor ', radius_cvff
         call fyi (trim(msg))
      end if
   else
      call warnings ('proc_cvff: "cvff" requires one parameter, a radius, in km')
   end if

   return
   
end subroutine proc_cvff


!***********************************************************************************************************************************
subroutine proc_cvou (command_help)

! Covariance matrix output file

   implicit none

   include 'mloc.inc'

   logical :: command_help

   if (command_help) then
      write (*,'(3(a/))')&
       '"cvou" specifies that the full cluster vector covariance matrix and some other',&
       'information will be written to a .cv file. This is needed for certain',&
       'specialized statistical tests. There are no arguments.'
      return
   end if

   if (params(1:1) .eq. ' ') then
      cv_out = .not.cv_out
   else if (params(1:2) .eq. 'on') then
      cv_out = .true.
   else if (params(1:3) .eq. 'off') then
      cv_out = .false.
   else
      call warnings ('proc_cvou: illegal argument')
   end if

   if (cv_out) then
      if (verbose_screen) call fyi ('proc_cvou: will create a covariance matrix (.cv) file')
   else
      if (verbose_screen) call fyi ('proc_cvou: will not create a covariance matrix (.cv) file')
   end if

   return
   
end subroutine proc_cvou


!***********************************************************************************************************************************
subroutine proc_cvtt (command_help)

! Cluster vector travel time error term

   implicit none

   include 'mloc.inc'

   logical :: command_help
   character(len=132) :: msg

   if (command_help) then
      write (*,'(9(a/))')&
      '"cvtt" adds an additional variance to differential time data. When using',&
      'differential arrival times from waveform cross-correlation, a pure reading',&
      'uncertainty is provided, but the actual variance includes a component from',&
      'inadequate theoretical differential TTs. This command allows that term to be',&
      'defined. It is only relevant if differential time data are used. The correction',&
      'is never added to bulletin data because this error term is already absorbed in',&
      'the empirical reading error. The two parameters are vmr (velocity model variance),',&
      'a percent value of velocity variance in the TT model, (default value 0.05, or 5%),',&
      'and cscale, a measure of the scale of the cluster in km (default 10 km).'
      return
   end if

   call decode (params, value, nvar)
   
   if (nvar .eq. 2) then
      vmr = value(1)
      cscale = value(2)
      if (verbose_screen) then
         write (msg,'(a,f4.2,f6.1)') 'proc_cvtt:', vmr, cscale
         call fyi (trim(msg))
      end if
   else
      call warnings ('proc_cvtt: "cvtt" requires two arguments')
   end if

   return
   
end subroutine proc_cvtt

!***********************************************************************************************************************************
subroutine proc_datf (command_help)

! DATF output file

   implicit none

   include 'mloc.inc'

   logical :: command_help

   if (command_help) then
      write (*,'(6(a/))')&
       '"datf" (DATF formatted output file) specifies that a .datf output will be written.',&
       'This is a single file with all the events written in MNF format, but with phase',&
       'IDs and flags as they are at the end of HD analysis. Flags and phase IDs can be',&
       'changed during the HD analysis. Arguments "on" and "off" can be used to set the',&
       'logical state explicitly. Issuing the command without an argument toggles the',&
       'logical state.'
      return
   end if

   if (params(1:1) .eq. ' ') then
      datfout = .not.datfout
   else if (params(1:2) .eq. 'on') then
      datfout = .true.
   else if (params(1:3) .eq. 'off') then
      datfout = .false.
   else
      call warnings ('proc_datf: illegal argument')
   end if

   if (datfout) then
      if (verbose_screen) call fyi ('proc_datf: will create a DATF file')
   else
      if (verbose_screen) call fyi ('proc_datf: will not create a DATF file')
   end if

   return
   
end subroutine proc_datf


!***********************************************************************************************************************************
subroutine proc_dbug (command_help)

! Debug mode, extra output in some routines

   implicit none

   include 'mloc.inc'

   logical :: command_help

   if (command_help) then
      write (*,'(4(a/))')&
       '"dbug" (debug) specifies that extra information about the HD analysis',&
       'will be written to the output window. Issuing the command without arguments',&
       'toggles the current state. The arguments "on" and "off" may be used to set the',&
       'state explicitly. Default is "off".'
      return
   end if

   if (params(1:1) .eq. ' ') then
      debug = .not.debug
   else if (params(1:2) .eq. 'on') then
      debug = .true.
   else if (params(1:3) .eq. 'off') then
      debug = .false.
   else
      call warnings ('proc_dbug: illegal argument')
   end if

   if (debug) then
      if (verbose_screen) call fyi ('proc_dbug: debugging output is turned on')
   else
      if (verbose_screen) call fyi ('proc_dbug: debugging output is turned off')
   end if

   return
   
end subroutine proc_dbug


!***********************************************************************************************************************************
subroutine proc_dcal (command_help)

! Direct calibration mode

   implicit none

   include 'mloc.inc'

   logical :: command_help

   if (command_help) then
      write (*,'(12(a/))')&
      '"dcal" specifies direct calibration mode, i.e., the hypocentroid is being located',&
      ' with data (usually near-distance data) that has minimal travel time error and so',&
      ' the hypocentroid can be taken as bias-free. The confidence ellipse of the',&
      ' hypocentroid is added to those of the cluster vectors to obtain final estimate of',&
      ' uncertainty for each event''s location. Four new output files are created:',&
      '   .dcal_phase_data - listing the data used for the hypocentroid',&
      '   _dcal.bash - GMT script, showing stations and raypaths used for hypocentroid',&
      '   _dcal.ps - Postscript file, from the _dcal.bash script',&
      '   .hdf_dcal - hdf with calibrated locations',&
      'Confidence ellipses in the .hdf_dcal file are cumulative, from hypocentroid and',&
      'cluster vector for each event. Can be toggled by using the command without an',&
      'argument, or set explicitly with arguments "on" or "off".'
      return
   end if

   if (params(1:1) .eq. ' ') then
      direct_cal = .not.direct_cal
   else if (params(1:2) .eq. 'on') then
      direct_cal = .true.
   else if (params(1:3) .eq. 'off') then
      direct_cal = .false.
   else
      call warnings ('proc_dcal: illegal argument')
   end if

   if (direct_cal) then
      dcal_pr = 'on '
   else
      dcal_pr = 'off'
   end if
   
   if (verbose_screen) call fyi ('proc_dcal: direct calibration is '//dcal_pr)

   return
   
end subroutine proc_dcal


!***********************************************************************************************************************************
subroutine proc_dem1 (command_help)

! Plot a digital elevation model in GMT.
! This command handles a lower-rez, globally-available DEM.

! References
! ETOPO1 <http://www.ngdc.noaa.gov/mgg/global/global.html>
! GLOBE <http://www.ngdc.noaa.gov/mgg/topo/globe.html>
! GINA (no longer hosted on-line)

   implicit none

   include 'mloc.inc'

   logical :: command_help

   if (command_help) then
      write (*,'(7(a/))')&
       '"dem1" (digital elevation model) specifies that commands to plot topography will be',&
       'added to the GMT script. Three digital elevation data sets are available, which are',&
       'selected by the argument:',&
       '   etopo1 = ETOPO1 1-minute gridded topography and bathymetry',&
       '   gina = GINA global gridded topography (30 arc sec) and bathymetry (2 arc min)',&
       '   globe = GLOBE 0.5-minute (~1 km) gridded topography',&
       'If no argument is given, plotting of topography is turned off'
      return
   end if

   if (params(1:1) .eq. ' ') then
      plot_globe = .false.
      plot_etopo1 = .false.
      plot_gina = .false.
   else if (params(1:6) .eq. 'etopo1') then
      plot_etopo1 = .true.
      plot_globe = .false.
      plot_gina = .false.
   else if (params(1:4) .eq. 'gina') then
      plot_gina = .true.
      plot_etopo1 = .false.
      plot_globe = .false.
   else if (params(1:5) .eq. 'globe') then
      plot_globe = .true.
      plot_etopo1 = .false.
      plot_gina = .false.
   else
      call warnings ('proc_dem1: illegal argument')
   end if

   plot_dem1 = (plot_globe .or. plot_etopo1 .or. plot_gina)
   
   if (plot_globe) then
      if (verbose_screen) call fyi ('proc_dem1: topography plotted using GLOBE')
   else if (plot_etopo1) then
      if (verbose_screen) call fyi ('proc_dem1: topography plotted using ETOPO1')
   else if (plot_gina) then
      if (verbose_screen) call fyi ('proc_dem1: topography plotted using GINA')            
   else
      if (verbose_screen) call fyi ('proc_dem1: topography will not be plotted')
   end if

   return
   
end subroutine proc_dem1


!***********************************************************************************************************************************
subroutine proc_dem2 (command_help)

! Plot a high-resolution digital elevation model in GMT.

   implicit none

   include 'mloc.inc'

   logical :: command_help, ex

   if (command_help) then
      write (*,'(6(a/))')&
       '"dem2" (digital elevation model) defines a DEM to be used in the base plot and other',&
       'plots that benefit from high-resolution topography. It assumes the user has provided',&
       'a ".grd" file and the argument to the command gives the pathname of that file, relative',&
       'to the MLOC working directory. A good website to obtain such files is:',&
       '   <http://topex.ucsd.edu/gmtsar/demgen/>',&
       'If no argument is given, plotting of topography is turned off'
      return
   end if

   if (params(1:1) .eq. ' ') then
      plot_dem2 = .false.
      if (verbose_screen) call fyi ('proc_dem2: hi-rez topography will not be plotted')
   else
      inquire (file=trim(params),exist=ex)
      if (ex) then
         dem2_filename = trim(params)
         plot_dem2 = .true.
         if (verbose_screen) call fyi ('proc_dem2: hi-rez topography plotted from '//trim(dem2_filename)) 
      else
         plot_dem2 = .false.
         call warnings ('proc_dem2: file '//trim(params)//' does not exist')
      end if
   end if

   return
   
end subroutine proc_dem2


!***********************************************************************************************************************************
subroutine proc_dep (iev, comd, command_help)

! Starting focal depth.
! The 4th character in the command is taken as the code for how starting focal depth was set (see mlocout_hdf.f90).
! This does not imply that depth is a fixed parameter. That is determined by commands 'frec' and 'freh'.
! c = cluster default depth 
! d = depth phases
! e = engineered (man-made explosion)
! f = fault model (InSAR, GPS, etc.)
! i = input data file
! l = local distance readings (more than 2-3 focal depths)
! m = mloc solution (with free depth)
! n = near-source station readings
! r = relocation (outside mloc) with free depth
! u = unknown
! w = waveform analysis
! With a couple of exceptions, the 4th character is carried as a flag in HDF files.
! Uncertainty in assigned depth can be specified by adding a second and/or third parameter. If one value
! is given in addition to the depth, it is taken as a symmetric depth uncertainty (+/-). If two additional
! values are given, the first is taken as + uncertainty (deeper) and the second is taken as - uncertainty
! (shallower). These values are not used in the relocation, but they are passed through in the HD
! file for future reference.

   implicit none

   include 'mloc.inc'

   integer :: iev
   logical :: command_help
   character(len=4) :: comd
   character(len=132) :: msg

   if (command_help) then
      write (*,'(26(a/))')&
       '"depX" (depth) specifies a starting focal depth. The 4th character is used to',&
       'specify the kind of information used to constrain the depth. At least one',&
       'argument is required, a positive depth value in km, and it can optionally take one',&
       'or two additional arguments to assign an uncertainty to that depth.',&
       '  Uncertainty in the assigned depth can be symmetric or asymmetric. If a single',&
       'value is appended to the depth estimate it is taken as a symmetric uncertainty. If',&
       'two values are appended the first is taken as plus (deeper), the second as minus',&
       '(shallower) uncertainty, in km. Both values should be positive.',&
       '  Recognized values for the 4th character are:',&
       '   c = cluster default depth ',&
       '   d = depth phases',&
       '   e = engineered (man-made explosion)',&
       '   f = fault model (InSAR, GPS, etc.)',&
       '   l = local-distance readings (more than 2-3 focal depths)',&
       '   m = mloc solution with free depth',&
       '   n = near-source station readings',&
       '   r = relocation (outside mloc) with free depth',&
       '   u = unconstrained (just a guess)',&
       '   w = waveform analysis',&
       'The command accepts any character in the 4th position so custom values can be used.',&
       '  The "depX" command can be issued multiple times, first to set a default depth',&
       'for all events ("depc"), then to set selected events at different depths. If it is issued',&
       'before any events are declared (by "memb" commands), it applies to all events in the',&
       'cluster that do not have a constrained depth declared in the input file. If issued after',&
       'a "memb" command it applies only to the current event, and it overrides the depth read',&
       'from the input file or from an HDF file.'
      return
   end if

   call decode (params, value, nvar)
   
   if (iev .eq. 0) then ! Cluster default depth
      depth_default = value(1)
      depth_default_c = 'c'
      if (nvar .gt. 1) call warnings ('proc_dep: cluster default depth does not take uncertainties')
   else ! Individual events
      if (nvar .ge. 1 .and. nvar .le. 3) then
         depth_cf(iev,1) = value(1)
         if (nvar .eq. 2) then
            depth_cf(iev,2) = value(2)
            depth_cf(iev,3) = value(2)
         else if (nvar .eq. 3) then
            depth_cf(iev,2) = value(2)
            depth_cf(iev,3) = value(3)
         end if
         if (comd(4:4) .eq. 't') then ! Special handling for backward compatibility
            depth_cf_c(iev) = 'u'
         else
            depth_cf_c(iev) = comd(4:4)
         end if
         if (verbose_screen) then
            write (msg,'(a,i3,1x,a,3f8.2,a)') 'proc_dep: event ', iev, depth_cf_c(iev), depth_cf(iev,1), depth_cf(iev,2),&
             depth_cf(iev,3)
            call fyi (trim(msg))
         end if
      else
         call warnings ('proc_dep: "dep " requires 1-3 parameters')
      end if      
   end if

   return
   
end subroutine proc_dep


!***********************************************************************************************************************************
subroutine proc_diff (command_help)

! Differential time data file

   implicit none

   include 'mloc.inc'

   logical :: ex, command_help

   if (command_help) then
      write (*,'(3(a/))')&
       '"diff" (differential time data) specifies the pathname of a file containing',&
       'differential time data for pairs of events in the cluster. The command can be',&
       'issued more than once, but only the last instance will be processed.'
      return
   end if

   if (params(1:1) .eq. ' ') then
      call warnings ('proc_diff: "diff" requires a pathname as argument')
   else
      inquire (file=trim(params),exist=ex)
      if (ex) then
         diffdat = .true.
         diffdatfilnam = trim(params)
         if (verbose_screen) call fyi ('proc_diff: differential time data: '//trim(diffdatfilnam))
      else
         call warnings ('proc_diff: file '//trim(params)//' does not exist')
      end if
   end if

   return
   
end subroutine proc_diff


!***********************************************************************************************************************************
subroutine proc_ellp (command_help)

! Plot an ellipse.

   implicit none

   include 'mloc.inc'

   integer :: i
   logical :: command_help
   character(len=132) :: msg

   if (command_help) then
      write (*,'(8(a/))')&
       '"ellp" (ellipse) is used to specify the parameters of an ellipse that will be',&
       'plotted in the map plots by GMT. The command takes five arguments:',&
       '   latitude of the center point',&
       '   longitude of the center point',&
       '   azimuth of the semi-major axis',&
       '   full length of the major axis, in km',&
       '   full length of the minor axis, in km',&
       'The command can be issued up to 10 times.'
      return
   end if

   call decode (params, value, nvar)
   
   if (nvar .eq. 5) then
      if (n_ellipse .lt. n_ellipse_nmax) then
         ellipse_plot = .true.
         n_ellipse = n_ellipse + 1
         do i = 1,5
            ellipse_data(n_ellipse,i) = value(i)
         end do
         call set_longitude_range (ellipse_data(n_ellipse,2), longitude_range)
         if (verbose_screen) then
            write (msg,'(a,5f8.3)') 'proc_ellp: ', ellipse_data(n_ellipse,1), ellipse_data(n_ellipse,2),&
             ellipse_data(n_ellipse,3), ellipse_data(n_ellipse,4), ellipse_data(n_ellipse,5)
            call fyi (trim(msg))
         end if
      else
         call warnings ('proc_ellp: maximum number of ellipses reached')
      end if
   else
      call warnings ('proc_ellp: "ellp" requires five parameters')
   end if

   return
   
end subroutine proc_ellp


!*****************************************************************************************
subroutine proc_epap (command_help)

! Make a plot of empirical path anomalies for a specific phase.

   implicit none
   
   include 'mloc.inc'
   
   logical :: command_help
   character(len=132) :: msg
   character(len=8) :: param1, param2, param3, param4
   
   if (command_help) then
      write (*,'(9(a/))')&
       '"epap" (empirical path anomaly plot) controls creation of a map of empirical path',&
       'anomalies for a specific phase. A symbol is plotted at each station that recorded',&
       'the phase of interest at least twice. The symbol is color-coded according to the',&
       'sign of the empirical path anomaly and the size is proportional to its absolute',&
       'value. Ray paths can optionally be drawn from the hypocentroid to the symbols. The',&
       'command takes three arguments: phase name, epicentral distance limit (degrees),',&
       'and a flag (0 or 1) which controls the plotting of raypaths. The epicentral',&
       'distance is used to set the boundaries of the map. The command can be issued',&
       'up to 6 times.'
      return
   end if

   if (params(1:1) .ne. ' ') then
      call decode2 (params, nvar, param1, param2, param3, param4)
      if (nvar .ne. 3) then
         write (msg,'(a,i1)') 'proc_epap: incorrect number of arguments = ', nvar
         call warnings (trim(msg))
         return
      end if
      if (verbose_screen) then
         call fyi ('proc_epap: phase = '//param1)
         call fyi ('proc_epap: distance limit   = '//param2)
         call fyi ('proc_epap: plot raypaths  = '//param3)
      end if
      epa_plot = .true.
      n_epa_plot = n_epa_plot + 1
      if (n_epa_plot .le. n_epa_plot_max) then
         epa_plot_phase(n_epa_plot) = param1
         call decode (param2, value, nvar)
         epa_plot_distance(n_epa_plot) = value(1)
         epa_plot_raypath(n_epa_plot) = (param3(1:1) .eq. '1')
      else
         call warnings ('proc_epap: maximum number of instances reached')
         n_epa_plot = n_epa_plot_max
      end if
   else
      call warnings ('proc_epap: "epap" requires three arguments')
   end if

end subroutine proc_epap


!***********************************************************************************************************************************
subroutine proc_eplt (command_help)

! Determines if a "confidence ellipse" plot is made.
! Can be toggled by using the command without an argument, or set explicitly with arguments "on" or "off".

   implicit none

   include 'mloc.inc'

   logical :: command_help

   if (command_help) then
      write (*,'(5(a/))')&
       '"eplt" (ellipse plot) specifies whether a "confidence ellipse" plot will be made.',&
       'This plot includes confidence ellipses for relative location but not event',&
       'numbers or relocation vectors. Issuing the command without arguments',&
       'toggles the current state. The arguments "on" and "off" may be used to set the',&
       'state explicitly. Default is "off".'
      return
   end if

   if (params(1:1) .eq. ' ') then
      eplt = .not.eplt
   else if (params(1:2) .eq. 'on') then
      eplt = .true.
   else if (params(1:3) .eq. 'off') then
      eplt = .false.
   else
      call warnings ('proc_eplt: illegal argument')
   end if

   if (eplt) then
      if (verbose_screen) call fyi ('proc_eplt: a confidence ellipse plot will be made')
   else
      if (verbose_screen) call fyi ('proc_eplt: a confidence ellipse plot will not be made')
   end if

   return
   
end subroutine proc_eplt


!***********************************************************************************************************************************
subroutine proc_even (iev, command_help)

! Event name

   implicit none

   include 'mloc.inc'

   integer :: iev
   logical :: command_help

   if (command_help) then
      write (*,'(3(a/))')&
       '"even" (event) is used to specify the name of the current event. This normally has',&
       'a form like YYYYMMDD.HHMM.SS (i.e., 16 characters). Maximum 30 characters. "even"',&
       'takes one argument.'
      return
   end if

   if (params(1:1) .ne. ' ') then
      evtnam(iev) = trim(params)
      if (verbose_screen) call fyi ('proc_even: '//evtnam(iev))
   else
      call warnings ('proc_even: "even" requires one argument, an event name')
   end if

   return
   
end subroutine proc_even


!***********************************************************************************************************************************
subroutine proc_fdhp (command_help)

! Determines if a histogram of focal depths is made.
! Can be toggled by using the command without an argument, or set explicitly with arguments "on" or "off".

   implicit none

   include 'mloc.inc'

   logical :: command_help

   if (command_help) then
      write (*,'(4(a/))')&
       '"fdhp" (focal depth histogram plot) specifies whether a plot will be made',&
       'of a histogram of focal depths. Issuing the command without arguments',&
       'toggles the current state. The arguments "on" and "off" may be used to set the',&
       'state explicitly. Default is "off".'
      return
   end if

   if (params(1:1) .eq. ' ') then
      fdhp = .not.fdhp
   else if (params(1:2) .eq. 'on') then
      fdhp = .true.
   else if (params(1:3) .eq. 'off') then
      fdhp = .false.
   else
      call warnings ('proc_fdhp: illegal argument')
   end if

   if (fdhp) then
      if (verbose_screen) call fyi ('proc_fdhp: a histogram of focal depths will be made')
   else
      if (verbose_screen) call fyi ('proc_fdhp: a histogram of focal depths will not be made')
   end if

   return
   
end subroutine proc_fdhp


!***********************************************************************************************************************************
subroutine proc_flag (command_help)

! Use of data flags

   implicit none

   include 'mloc.inc'

   logical :: command_help

   if (command_help) then
      write (*,'(3(a/))')&
       '"flag" (data flags) determines if data flags encountered in input files will be',&
       'honored ("on") or ignored ("off"). Issuing the command with no argument toggles',&
       'the current state. Default is "on".'
      return
   end if

   if (params(1:1) .eq. ' ') then
      dflag = .not.dflag
   else if (params(1:2) .eq. 'on') then
      dflag = .true.
   else if (params(1:3) .eq. 'off') then
      dflag = .false.
   else
      call warnings ('proc_flag: illegal argument')
   end if

   if (dflag) then
      dflag_pr = 'on '
   else
      dflag_pr = 'off'
   end if
   if (verbose_screen) call fyi ('proc_flag: Use of data flags is '//dflag_pr)

   return
   
end subroutine proc_flag


!***********************************************************************************************************************************
subroutine proc_fmap (command_help)

! Fault map for GMT scripts.

   implicit none

   include 'mloc.inc'

   logical :: ex, command_help

   if (command_help) then
      write (*,'(6(a/))')&
       '"fmap" specifies a file of digitized faults or other linear features to be plotted',&
       'in GMT. The standard GMT format - longitude, latitude - is used, in free format.',&
       'The symbol ">" can be used to separate multiple line segments. The pathname',&
       'relative to mloc is given. The normal place for these files is the directory',&
       '"/tables/faults". The command can be issued multiple times. Formatting for specific',&
       'fault types is done with the command "fplt".'
      return
   end if

   if (params(1:1) .eq. ' ') then
      call warnings ('proc_fmap: "fmap" requires a pathname as argument')
   else
      if (n_fault_map .lt. n_fault_map_nmax) then
         inquire (file=trim(params),exist=ex)
         if (ex) then
            fault_map = .true.
            n_fault_map = n_fault_map + 1
            fault_map_filename(n_fault_map) = trim(params)
            psxy_wsf(n_fault_map) = psxy_wsf0
            if (verbose_screen) call fyi ('proc_fmap: fault map in file: '//trim(fault_map_filename(n_fault_map)))
            if (verbose_screen) call fyi ('   ...with plot options: '//psxy_wsf(n_fault_map))
         else
            call warnings ('proc_fmap: file '//trim(params)//' does not exist')
         end if
      else
         call warnings ('proc_fmap: maximum number of fault map files reached!')
      end if
   end if

   return
   
end subroutine proc_fmap


!***********************************************************************************************************************************
subroutine proc_fplt (command_help)

! Plotting options for specific types of faults in GMT (see command FMAP).

   implicit none

   include 'mloc.inc'

   logical :: command_help

   if (command_help) then
      write (*,'(6(a/))')&
       '"fplt" (fault plotting) is used to specify plotting options for specific types of',&
       'faults in GMT (see command "fmap"). The argument for the command is a text string',&
       'with fully-formed "-W -Sf" options for psxy. A blank argument resets the default',&
       'string. The current value of psxy_wsf0 is loaded into a fault-specific variable',&
       'for each fault tmap that is read. Reissue the "fplt" command before a new "fmap"',&
       'command to change the plot options for that fault.'
      return
   end if

   if (params(1:1) .eq. ' ') then
      psxy_wsf0 = '-Wthick/255/0/0'
      if (verbose_screen) call fyi ('proc_fplt: fault plotting options reset to '//trim(psxy_wsf0))
   else
      psxy_wsf0 = trim(params)
   end if

   return
   
end subroutine proc_fplt

          
!***********************************************************************************************************************************
subroutine proc_frec (iev, command_help)

! Free parameters for cluster vectors

   implicit none

   include 'mloc.inc'

   integer :: iev, i
   logical :: command_help
   character(len=132) :: msg

   if (command_help) then
      write (*,'(7(a/))')&
       '"frec" (free parameters) is used to specify which location parameters will be free',&
       'and which will be fixed for cluster vectors. It requires four arguments that are',&
       'either 0 or 1. 0 indicates a fixed parameter, 1 indicates a free parameter. The',&
       'order is latitude, longitude, depth, origin time. If "frec" is invoked before any',&
       'events are declared, the values apply to all events. If invoked after an event is',&
       'declared, the values only apply to the current event.'
      return
   end if

   call decode (params, value, nvar)
   
   if (nvar .eq. 4) then
      if (iev .eq. 0) then ! Apply to all events
         do i = 1,nevmax
            if (nint(value(1)) .eq. 0) then
               latf(i) = .false.
               fixpr(i,1) = 'fixed'
            else if (nint(value(1)) .eq. 1) then
               latf(i) = .true.
               fixpr(i,1) = 'free '
            else
               call warnings ('proc_frec: error in input')
            end if
            if (nint(value(2)) .eq. 0) then
               lonf(i) = .false.
               fixpr(i,2) = 'fixed'
            else if (nint(value(2)) .eq. 1) then
               lonf(i) = .true.
               fixpr(i,2) = 'free '
            else
               call warnings ('proc_frec: error in input')
            end if
            if (nint(value(3)) .eq. 0) then
               depthf(i) = .false.
               fixpr(i,3) = 'fixed'
            else if (nint(value(3)) .eq. 1) then
               depthf(i) = .true.
               fixpr(i,3) = 'free '
            else
               call warnings ('proc_frec: error in input')
            end if
            if (nint(value(4)) .eq. 0) then
               timef(i) = .false.
               fixpr(i,4) = 'fixed'
            else if (nint(value(4)) .eq. 1) then
               timef(i) = .true.
               fixpr(i,4) = 'free '
            else
               call warnings ('proc_frec: error in input')
            end if
         end do
         if (verbose_screen) then
            write (msg,'(2a,1x,a,1x,a,1x,a)') 'proc_frec: all events ', fixpr(1,1), fixpr(1,2), fixpr(1,3), fixpr(1,4)
            call fyi (trim(msg))
         end if
      else ! Application to individual events
         if (nint(value(1)) .eq. 0) then
            latf(iev) = .false.
            fixpr(iev,1) = 'fixed'
         else if (nint(value(1)) .eq. 1) then
            latf(iev) = .true.
            fixpr(iev,1) = 'free '
         else
            call warnings ('proc_frec: error in input')
         end if
         if (nint(value(2)) .eq. 0) then
            lonf(iev) = .false.
            fixpr(iev,2) = 'fixed'
         else if (nint(value(2)) .eq. 1) then
            lonf(iev) = .true.
            fixpr(iev,2) = 'free '
         else
            call warnings ('proc_frec: error in input')
         end if
         if (nint(value(3)) .eq. 0) then
            depthf(iev) = .false.
            fixpr(iev,3) = 'fixed'
         else if (nint(value(3)) .eq. 1) then
            depthf(iev) = .true.
            fixpr(iev,3) = 'free '
         else
            call warnings ('proc_frec: error in input')
         end if
         if (nint(value(4)) .eq. 0) then
            timef(iev) = .false.
            fixpr(iev,4) = 'fixed'
         else if (nint(value(4)) .eq. 1) then
            timef(iev) = .true.
            fixpr(iev,4) = 'free '
         else
            call warnings ('proc_frec: error in input')
         end if
         if (verbose_screen) then
            write (msg,'(a,i3,1x,a,1x,a,1x,a,1x,a)') 'proc_frec: event ', iev, fixpr(iev,1), fixpr(iev,2), fixpr(iev,3),&
             fixpr(iev,4)
            call fyi (trim(msg))
         end if
      end if
   else
      call warnings ('proc_frec: "frec" requires four parameters')
   end if

   return
   
end subroutine proc_frec


!***********************************************************************************************************************************
subroutine proc_freh (command_help)

! Free parameters for hypocentroid

   implicit none

   include 'mloc.inc'

   logical :: command_help
   character(len=132) :: msg

   if (command_help) then
      write (*,'(4(a/))')&
       '"freh" (free parameters) is used to specify which location parameters will be free',&
       'and which will be fixed for the hypocentroid. It requires four arguments that are',&
       'either 0 or 1. 0 indicates a fixed parameter, 1 indicates a free parameter. The',&
       'order is latitiude, longitude, depth, origin time.'
      return
   end if
  
   call decode (params, value, nvar)
   
   if (nvar .eq. 4) then
      if (nint(value(1)) .eq. 0) then
         latfh = .false.
         fixprh(1) = 'fixed'
      else if (nint(value(1)) .eq. 1) then
         latfh = .true.
         fixprh(1) = 'free '
      else
         call warnings ('proc_freh: error in input')
      end if
      if (nint(value(2)) .eq. 0) then
         lonfh = .false.
         fixprh(2) = 'fixed'
      else if (nint(value(2)) .eq. 1) then
         lonfh = .true.
         fixprh(2) = 'free '
      else
         call warnings ('proc_freh: error in input')
      end if
      if (nint(value(3)) .eq. 0) then
         depthfh = .false.
         fixprh(3) = 'fixed'
      else if (nint(value(3)) .eq. 1) then
         depthfh = .true.
         fixprh(3) = 'free '
      else
         call warnings ('proc_freh: error in input')
      end if
      if (nint(value(4)) .eq. 0) then
         timefh = .false.
         fixprh(4) = 'fixed'
      else if (nint(value(4)) .eq. 1) then
         timefh = .true.
         fixprh(4) = 'free '
      else
         call warnings ('proc_freh: error in input')
      end if
      if (verbose_screen) then
         write (msg,'(a,1x,a,1x,a,1x,a,1x,a)') 'proc_freh: ', fixprh(1), fixprh(2), fixprh(3), fixprh(4)
         call fyi (trim(msg))
      end if
   else
      call warnings ('proc_freh: "freh" requires four parameters')
   end if

   return
   
end subroutine proc_freh


!***********************************************************************************************************************************
subroutine proc_help (command_help)

! Command summaries and full help

   implicit none

   include 'mloc.inc'

   logical :: command_help

   if (command_help) then
      write (*,'(6(a/))')&
       '"help" is used to provide detailed information about the purpose and usage of the',&
       'commands in mloc. If it is issued without an argument a list of all commands is',&
       'returned, organized by general category, and with a short descriptions of the function',&
       'of each command. A particular command may be given as an argument, in which case a more',&
       'detailed description of that command is returned. Only a single command can be given as',&
       'argument.'
      return
   end if

   write (*,'(a)') 'Calibration'
   write (*,'(t3,a)') 'cal : calibration data for current event'
   write (*,'(t3,a)') 'ctyp: treatment of calibration events in indirect calibration'
   write (*,'(t3,a)') 'dcal: direct calibration'

   write (*,'(a)') 'Informational'
   write (*,'(t3,a)') 'anno: annotation'
   write (*,'(t3,a)') 'comm: comment line'
   write (*,'(t3,a)') 'dbug: extra logging for debugging'
   write (*,'(t3,a)') 'help: details about commands'
   write (*,'(t3,a)') 'revi: review current relocation control parameters'
   write (*,'(t3,a)') 'vlog: set verbose mode for logging'
   write (*,'(t3,a)') 'vscr: set verbose mode for screen display'

   write (*,'(a)') 'Input'
   write (*,'(t3,a)') 'cfil: specify a command file'
   write (*,'(t3,a)') 'diff: differential time data'
   write (*,'(t3,a)') 'even: event name'
   write (*,'(t3,a)') 'inpu: specify an event data file'
   write (*,'(t3,a)') 'kill: kill a block of events'
   write (*,'(t3,a)') 'memb: start a new event (or kill an event with the "kill" argument)'

   write (*,'(a)') 'Inversion'
   write (*,'(t3,a)') 'bias: hypocentroid bias correction'
   write (*,'(t3,a)') 'clim: epicentral distance limits for cluster vectors'
   write (*,'(t3,a)') 'flag: use of data flags'
   write (*,'(t3,a)') 'frec: free parameters for cluster vectors'
   write (*,'(t3,a)') 'freh: free parameters for the hypocentroid'
   write (*,'(t3,a)') 'hlim: epicentral distance limits for hypocentroid'
   write (*,'(t3,a)') 'phyp: set "use only P arrivals for hypocentroid" flag'
   write (*,'(t3,a)') 'pttt: perfect theoretical travel times for hypocentroid'
   write (*,'(t3,a)') 'run : begin the relocation process'
   write (*,'(t3,a)') 'shcl: hypocentroid convergence limits'
   write (*,'(t3,a)') 'step: number of iterations to run'
   write (*,'(t3,a)') 'stop: stop processing'
   write (*,'(t3,a)') 'tikh: Tikhonov regularization'
   write (*,'(t3,a)') 'weig: residual weighting by reading error'

   write (*,'(a)') 'Miscellaneous'
   write (*,'(t3,a)') 'lonr: set longitude range (-180 to 180 or 0 to 360)'

   write (*,'(a)') 'Output'
   write (*,'(t3,a)') 'bloc: BayesLoc output file'
   write (*,'(t3,a)') 'ccat: COMCAT output file'
   write (*,'(t3,a)') 'datf: output MNF bulletin of all events in their final form'
   write (*,'(t3,a)') 'lres: LRES file of large cluster residuals'
   write (*,'(t3,a)') 'mdou: map_dat output file for GMT'
   write (*,'(t3,a)') 'mech: focal mechanism data (from MNF event files) written to a file'
   write (*,'(t3,a)') 'oldr: output of phase readings over a limited distance range'
   write (*,'(t3,a)') 'puke: PUKE formatted output file'
   write (*,'(t3,a)') 'subc: select a subcluster based on data for direct calibration'
   write (*,'(t3,a)') 'tomo: tomography output files'
   write (*,'(t3,a)') 'ttou: Empirical TTs for a specific phase'   

   write (*,'(a)') 'Phase Identification'
   write (*,'(t3,a)') 'phid: toggle phase re-identification'
   write (*,'(t3,a)') 'ppri: phases that cannot be renamed'
   write (*,'(t3,a)') 'skip: skip readings of a given phase (also by station and author)'

   write (*,'(a)') 'Plotting'
   write (*,'(t3,a)') 'cptf: color palette table for topography'
   write (*,'(t3,a)') 'dem1: plot regional-scale topography in GMT'
   write (*,'(t3,a)') 'dem2: plot high-resolution topography in smaller-scale GMT plots'
   write (*,'(t3,a)') 'ellp: plot an ellipse'
   write (*,'(t3,a)') 'epap: plot of empirical path anomalies'
   write (*,'(t3,a)') 'eplt: make a map of locations with only confidence ellipses'
   write (*,'(t3,a)') 'fdhp: make a histogram of focal depths'
   write (*,'(t3,a)') 'fmap: digital fault map for GMT script'
   write (*,'(t3,a)') 'fplt: plotting of faults'
   write (*,'(t3,a)') 'plot: plotting of selected events'
   write (*,'(t3,a)') 'pltt: travel time plots'
   write (*,'(t3,a)') 'rdpp: relative depth phase plots for individual events'
   write (*,'(t3,a)') 'splt: make a seismicity map, same-size symbols for locations'
   write (*,'(t3,a)') 'star: plot a star to highlight a specific event'
   write (*,'(t3,a)') 'stat: plot a triangle to indicate a station location'
   write (*,'(t3,a)') 'tt5e: single-event local distance (tt5) plot'
   write (*,'(t3,a)') 'tt5s: single-station local distance (tt5) plot'
   write (*,'(t3,a)') 'vect: which event shift vectors to plot'
   write (*,'(t3,a)') 'xsec: cross-sections'

   write (*,'(a)') 'Residuals and Uncertainties'
   write (*,'(t3,a)') 'cvff: Cluster vector fudge factor (non-gaussian contribution to uncertainty)'
   write (*,'(t3,a)') 'cvou: Output full covariance matrices'
   write (*,'(t3,a)') 'cvtt: cluster vector travel time error'
   write (*,'(t3,a)') 'mare: minimum allowed reading error'
   write (*,'(t3,a)') 'rels: set reading errors for local stations'
   write (*,'(t3,a)') 'rfil: reading error file (.rderr)'
   write (*,'(t3,a)') 'tfil: travel time spread file (.ttsprd)'
   write (*,'(t3,a)') 'wind: windowing of residuals'

   write (*,'(a)') 'Starting locations'
   write (*,'(t3,a)') 'dep : set starting focal depth'
   write (*,'(t3,a)') 'lat : set starting latitude'
   write (*,'(t3,a)') 'long: set starting longitude'
   write (*,'(t3,a)') 'pert: perturbation to all starting locations'
   write (*,'(t3,a)') 'rhdf: read starting locations from an HDF file'
   write (*,'(t3,a)') 'time: set starting origin time'

   write (*,'(a)') 'Stations'
   write (*,'(t3,a)') 'bdps: list of stations suspected of reporting bogus depth phases'
   write (*,'(t3,a)') 'nsmd: search NEIC station metadata for missing station codes'
   write (*,'(t3,a)') 'radf: read agency and deployment fields to resolve station code conflicts'
   write (*,'(t3,a)') 'skip: skip readings from a given station (also by phase and author)'
   write (*,'(t3,a)') 'sstn: supplemental station file'

   write (*,'(a)') 'Travel Times'
   write (*,'(t3,a)') 'bptc: Bounce point topography correction'
   write (*,'(t3,a)') 'corr: Station elevation correction'
   write (*,'(t3,a)') 'lgtt: Lg travel time calculation'
   write (*,'(t3,a)') 'lmod: local velocity model'
   write (*,'(t3,a)') 'secv: station elevation correction velocities'
   write (*,'(t3,a)') 'taup: Global TT model, using Tau-P formulation'
   write (*,'(t3,a)') 'tptt: T-phase travel time calculation'

   write (*,'(/t3,a/)') 'For further information, follow the help command with a command name'

   return
   
end subroutine proc_help         
      

!***********************************************************************************************************************************
subroutine proc_hlim (command_help)

! Epicentral distance limits for hypocentroid

   implicit none

   include 'mloc.inc'

   integer :: i, j
   logical :: command_help
   character(len=132) :: msg

   if (command_help) then
      write (*,'(5(a/))')&
       '"hlim" (hypocentroid limits) specifies a set of up to three ranges in epicentral',&
       'distance that will be used to estimate the hypocentroid. This allows skipping',&
       'readings at short distances, in the Pdiff range, and near the caustic for PKP',&
       'phases. Since each range requires two values, this command requires 2, 4, or 6',&
       'arguments.'
      return
   end if

   if (params(1:1) .eq. ' ') then
      nhlim = 1
      hlim(1,1) = 30.0
      hlim(1,2) = 90.0
   else 
      call decode (params, value, nvar)
      if (nvar .eq. 2 .or. nvar .eq. 4 .or. nvar .eq. 6) then
         nhlim = nvar/2
         do i = 1,nhlim
            j = 2*i
            hlim(i,1) = value(j-1)
            hlim(i,2) = value(j)
            if (verbose_screen) then
               write (msg,'(a,i1,a,2f6.1)') 'proc_hlim: distance interval ', i, ':',hlim(i,1), hlim(i,2)
               call fyi (trim(msg))
            end if
         end do
      else
         call warnings ('proc_hlim: "hlim" requires 1-3 pairs of parameters')
      end if
   end if

   return
   
end subroutine proc_hlim


!***********************************************************************************************************************************
subroutine proc_inpu (iev, command_help)

! Specify an input data file for the current event.

   implicit none

   include 'mloc.inc'

   integer :: iev
   logical :: ex, command_help
   character(len=132) :: msg

   if (command_help) then
      write (*,'(4(a/))')&
       '"inpu" (input) is used to specify the name of the file containing the phase arrival',&
       'time information for an event. The command requires a single argument. The standard',&
       'form of the base filename is YYYYMMDD.HHMM.SS. The only supported format is MNF (MLOC',&
       'Native Format). The maximum length of the input filename is 20 characters.'
      return
   end if

   if (iev .eq. 0) then
      call warnings ('proc_inpu: "memb" must be invoked before "inpu"')
   else if (params(1:1) .ne. ' ') then
      if (len(trim(params)) .gt. 20) then
         call warnings ('proc_inpu: invalid filename length: '//trim(params))
      end if
      infile20(iev) = trim(params)
      if (verbose_screen) then
         write (msg,'(a,i3,1x,a)') 'proc_inpu: ', iev, infile20(iev)
         call fyi (trim(msg))
      end if
      infile(iev) = trim(datadir)//dirsym//infile20(iev)
      inquire (file=infile(iev),exist=ex)
      if (.not.ex) then
         call warnings ('proc_inpu: file '//trim(infile(iev))//' does not exist')
         infile(iev) = ' '
      end if
   else
      call warnings ('proc_inpu: "inpu" requires one argument, a filename')
   end if            

   return
   
end subroutine proc_inpu


!***********************************************************************************************************************************
subroutine proc_kill (command_help)

! Kill a block of events

   implicit none

   include 'mloc.inc'

   logical :: command_help

   if (command_help) then
      write (*,'(4(a/))')&
      '"kill" provides a convenient way to kill blocks of events in a cluster. It',&
      'requires one argument, either "on" or "off". After the "kill on" command is issued,',&
      'no other commands will be processed until the "kill off" command is processed. To kill',&
      'a single event, the "memb" command with the argument "kill" is preferred.'
      return
   end if
   
   if (params(1:2) .eq. 'on') then
      kill_all = .true.
   else if (params(1:3) .eq. 'off') then
      kill_all = .false.
   else
      call warnings ('proc_kill: illegal argument "'//trim(params)//'" (see the MEMB command to kill a single event)')
   end if

   return
   
end subroutine proc_kill


!***********************************************************************************************************************************
subroutine proc_lat (iev, command_help)

! Starting latitude

   implicit none

   include 'mloc.inc'

   integer :: iev, i
   logical :: command_help
   character(len=132) :: msg

   if (command_help) then
      write (*,'(7(a/))')&
       '"lat" (latitude) specifies a starting latitude. It requires one argument, a latitude',&
       'value. If it is issued before any events are declared (by "memb" commands), it',&
       'applies to all events in the cluster. Otherwise it applies only to the current',&
       'event. The "lat" command can be issued many times, first to set a default latitude',&
       'for all events, then to set selected events at different latitudes. If the "lat"',&
       'command is not applied to an event, the latitude read from the input data file will',&
       'be used as the starting value.'
      return
   end if
   
   call decode (params, value, nvar)
   
   if (nvar .eq. 1) then
      if (iev .eq. 0) then ! Apply to all events
         do i = 1,nevmax
            lat_cf(i) = value(1)
         end do
         if (verbose_screen) then
             write (msg,'(a,f8.3)') 'proc_lat: all events ', lat_cf(1)
             call fyi (trim(msg))
         end if
      else ! Application to individual events
         lat_cf(iev) = value(1)
         if (verbose_screen) then
            write (msg,'(a,i3,f8.3)') 'proc_lat: event ', iev, lat_cf(iev)
            call fyi (trim(msg))
         end if
      end if
   else
      call warnings ('proc_lat: "lat" requires one parameter')
   end if

   return
   
end subroutine proc_lat


!***********************************************************************************************************************************
subroutine proc_lgtt (command_help)

! Lg travel time calculations.

   implicit none

   include 'mloc.inc'

   logical :: command_help
   character(len=132) :: msg

   if (command_help) then
      write (*,'(2(a/))')&
       '"lgtt" specifies an intercept (sec), slope (sec/degree) and minimum epicentral',&
       ' distance (degrees) to be used in calculating travel tiems for Lg.'
      return
   end if

   call decode (params, value, nvar)
   
   if (nvar .eq. 3) then
      lg_a = value(1)
      lg_b = value(2)
      lg_min = value(3)
      if (verbose_screen) then
         write (msg,'(a,3f8.2)') 'proc_lgtt:', lg_a, lg_b, lg_min
         call fyi (trim(msg))
      end if
   else
      call warnings ('proc_lgtt: "lgtt" requires three arguments')
   end if

   return
   
end subroutine proc_lgtt


!***********************************************************************************************************************************
subroutine proc_lmod (command_help)

! Local velocity model

   implicit none

   include 'mloc.inc'

   integer :: ierr
   logical :: command_help, op, ex
   character(len=100) :: filnam
   character(len=132) :: msg

   if (command_help) then
      write (*,'(7(a/))')&
       '"lmod" (local model) specifies a local velocity model to be used to calculate travel',&
       'times at short epicentral distance ranges. Local velocity models are normally',&
       'stored in the "tables/crust/" subdirectory. The format of local velocity model',&
       'files is the format used by HYPOSAT. If a custom velocity model is not specified,',&
       'the global tau-p model is used for all distances. The command takes one argument,',&
       'the pathname (relative to mloc) of the model file. The distance and depth ranges for',&
       'which it will be used are specified in the file itself.'
      return
   end if

   if (params(1:1) .eq. ' ') then
      call warnings ('proc_lmod: "lmod" requires a pathname as argument')
   else 
      filnam = trim(params)
      inquire (file=filnam,exist=ex)
      if (.not.ex) then
         call warnings ('proc_lmod: '//trim(filnam)//' does not exist')
         dlimlocmod = 0.
         zlimlocmod = 0.
         locmod = .false.
         return
      end if
      open (io_locmod,file=filnam,status='old')
      inquire (unit=io_locmod,opened=op,name=locmodfname)
      locmod = op
      if (op) then
         call rd_loc_mod (io_locmod, dlimlocmod, zlimlocmod, ierr) ! Read local velocity model
         if (ierr .eq. 0) then
            if (verbose_screen) then
               write (msg,'(2a,2f8.1)') 'proc_lmod: ', trim(locmodfname), dlimlocmod, zlimlocmod
               call fyi (trim(msg))
            end if
         else
            call warnings ('proc_lmod: an error occurred reading the file '//trim(locmodfname))
            dlimlocmod = 0.
            zlimlocmod = 0.
            locmod = .false.
         end if
         close (io_locmod)
      else
         call warnings ('proc_lmod: file '//trim(locmodfname)//' was not opened')
         dlimlocmod = 0.
         zlimlocmod = 0.
         locmod = .false.
      end if
   end if

   return
   
end subroutine proc_lmod


!***********************************************************************************************************************************
subroutine proc_long (iev, command_help)

! Starting longitude

   implicit none

   include 'mloc.inc'

   integer :: iev, i
   logical :: command_help
   character(len=132) :: msg

   if (command_help) then
      write (*,'(7(a/))')&
       '"long" (longitude) specifies a starting longitude. It requires one argument, a',&
       'longitude value. If it is issued before any events are declared (by "memb"',&
       'commands), it applies to all events in the cluster. Otherwise it applies only to',&
       'the current event. The "long" command can be issued many times, first to set a',&
       'default longitude for all events, then to set selected events at different',&
       'longitudes. If the "long" command is not applied to an event, the longitude read',&
       'from the input data file will be used as the starting value.'
      return
   end if

   call decode (params, value, nvar)
   
   if (nvar .eq. 1) then
      if (iev .eq. 0) then ! Apply to all events
         do i = 1,nevmax
            lon_cf(i) = value(1)
         end do
         if (verbose_screen) then
            write (msg,'(a,f8.3)') 'proc_long: all events ', lon_cf(1)
            call fyi (trim(msg))
         end if
      else ! Application to individual events
         lon_cf(iev) = value(1)
         if (verbose_screen) then
            write (msg,'(a,i3,f8.3)') 'proc_long: event ', iev, lon_cf(iev)
            call fyi (trim(msg))
         end if
      end if
   else
      call warnings ('proc_long: "long" requires one parameter')
   end if

   return
   
end subroutine proc_long


!*****************************************************************************************
subroutine proc_lonr (command_help)

! Longitude range

   implicit none
   
   include 'mloc.inc'
   
   logical :: command_help
   character(len=132) :: msg

   if (command_help) then
      write (*,'(7(a/))')&
       '"lonr" (longitude range) specifies the range to which longitudes should be converted.',&
       'The command takes a single integer argument that gives the center of the longitude range.',&
       'Two values are accepted:',&
       '  0  : -180° ≤ longitude < 180° (default)',&
       '  180:    0° ≤ longitude < 360°',&
       'A value should be chosen that keeps longitudes for all events in the cluster in the',&
       'same range.'
      return
   end if

   call decode (params, value, nvar)
    
   if (nvar .eq. 1) then 
      if (abs(value(1)) .lt. 1.) then
         longitude_range = 0
      else if (abs(value(1)-180.) .lt. 1) then
         longitude_range = 180
      else
         write (msg,'(a,i4)') 'proc_lonr: argument not recognized: ', value(1)
         call warnings (trim(msg))
      end if
   else
      call warnings ('proc_lonr: "lonr" requires one argument')
   end if
   
   return
   
end subroutine proc_lonr


!*****************************************************************************************
subroutine proc_lres (command_help)

! Output of file with large cluster residuals

   implicit none

   include 'mloc.inc'

   logical :: command_help
   character(len=132) :: msg

   if (command_help) then
      write (*,'(3(a/))')&
       '"lres" (large residual) specifies that a .lres output file will be written,',&
       'containing all readings with cluster residuals (eci) larger than the value given',&
       'by the single argument.'
      return
   end if
    
   if (params(1:1) .eq. ' ') then
      lres = 4.0
   else
      call decode (params, value, nvar)
      if (nvar .eq. 1) then
         lres = value(1)
         lresout = .true.
         if (verbose_screen) then
            write (msg,'(a,f6.2)') 'proc_lres: ', lres
            call fyi (trim(msg))
         end if
      else
         call warnings ('proc_lres: "lres" requires one argument, the lower limit on large residuals')
      end if 
   end if

   return
   
end subroutine proc_lres


!***********************************************************************************************************************************
subroutine proc_mare (command_help)

! Minimum allowed reading errors.
! Do not confuse with command "rels" which sets the reading errors explicitly for local
! distance readings.

   implicit none

   include 'mloc.inc'

   logical :: command_help
   character(len=132) :: msg

   if (command_help) then
      write (*,'(8(a/))')&
       '"mare" (minimum allowed reading error) specifies the the minimum reading errors that',&
       'will be allowed, regardless of what value is read from a rderr file. Three arguments',&
       'are required:',&
       '   1 - the minimum value for local phases',&
       '   2 - the minimum allowed value for phases beyond local distance',&
       '   3 - the minimum allowed value for teleseismic depth phases',&
       'Local phases are defined as those with "g" or "b" as the second character of the',&
       'phase name.'
      return
   end if

   call decode (params, value, nvar)
   
   if (nvar .eq. 3) then
      rderr_min_loc = value(1)
      rderr_min = value(2)
      rderr_min_depth = value(3)
      if (rderr_min_loc .gt. rderr_min) then
         write (msg,'(a)') 'proc_mare: this may be an old "mare" command - the order of arguments has changed'
         call warnings (trim(msg))
      end if
      if (verbose_screen) then
         write (msg,'(a,3f6.2)') 'proc_mare: ', rderr_min_loc, rderr_min, rderr_min_depth
         call fyi (trim(msg))
      end if
   else
      call warnings ('proc_mare: "mare" requires three arguments')
   end if 

   return
   
end subroutine proc_mare


!***********************************************************************************************************************************
subroutine proc_mdou (command_help)

! Output file '.map_dat' for import into GMT scripts.

   implicit none

   include 'mloc.inc'

   logical :: command_help

   if (command_help) then
      write (*,'(4(a/))')&
       '"mdpou" (map_dat output file) specifies that a special file containing lat-lon data',&
       'will be written. It will have the filename suffix ".map_dat". This file is designed',&
       'for easy import by a GMT script for additional plotting. If no argument is given,',&
       'the status is toggled. The arguments "off" and "on" can be used to set the state explicitly.'
      return
   end if

   if (params(1:1) .eq. ' ') then
      md_out = .not.md_out
   else if (params(1:2) .eq. 'on') then
      md_out = .true.
   else if (params(1:3) .eq. 'off') then
      md_out = .false.
   else
      call warnings ('proc_mdpou: illegal argument')
   end if

   if (md_out) then
      if (verbose_screen) call fyi ('proc_mdou: will create a .map_dat file')
   else
      if (verbose_screen) call fyi ('proc_mdou: will not create a .map_dat file')
   end if

   return
   
end subroutine proc_mdou


!***********************************************************************************************************************************
subroutine proc_mech (iev, command_help)

! Focal mechanism data
! Data is provided in free format by the user, and passed to an output file along with hypocentral data for plotting.

   implicit none

   include 'mloc.inc'

   integer :: iev
   logical :: command_help
   character(len=132) :: msg

   if (command_help) then
      write (*,'(5(a/))')&
       '"mech" allows focal mechanism or moment tensor data to be specified and carried',&
       'through to an output file, where it is appended to the relocated hypocenter info.',&
       'The intended use is to create a data file for plotting focal mechanisms. However',&
       'the argument to the command is simply a text string, so it can be any format.',&
       'MLOC makes no use of the information contained in the text string.'
      return
   end if

   if (params(1:1) .eq. ' ') then
      call warnings ('proc_mech: "mech" requires a text string as argument')
   else
      if (.not.focal_mech) focal_mech = .true.
      focal_mech_iev(iev) = .true.
      focal_mech_line(iev) = trim(params)
      if (verbose_screen) then
         write (msg,'(a,i3,1x,a)') 'proc_mech: event ', iev, focal_mech_line(iev)
         call fyi (trim(msg))
      end if
   end if

   return
   
end subroutine proc_mech


!***********************************************************************************************************************************
subroutine proc_memb (iev, command_help)

! Start the specification of a new event

   implicit none

   include 'mloc.inc'

   integer :: iev
   logical :: command_help
   character(len=132) :: msg

   if (command_help) then
      write (*,'(7(a/))')&
       '"memb" (member) is used to define a new event in the cluster or to cause an event in',&
       'the command file to be skipped (killed). The command takes no argument if a new event',&
       'is being defined. When the command is issued, the event counter is incremented and',&
       'specification of a new event starts. To use this command to kill an event, add the',&
       'argument "kill". Any text following "kill" is ignored, but it can be used for a',&
       'comment about the reason for killing the event. This is the preferred way to kill',&
       'a single event; see the command "kill" to skip blocks of events.'
      return
   end if
   
   if (kill_one) kill_one=.false. ! Re-initialize the memb/kill function after one event.

   if (params(1:4) .eq. 'kill') then
      kill_one = .true.
      kill_reason = params
      return
   end if

   if (params(1:1) .eq. ' ') then
      iev = iev + 1
   else 
      call warnings ('proc_memb: "memb" takes no argument except "kill"')
   end if
   
   if (iev .gt. nev) nev = iev ! nev is incremented here.
   
   if (nev .gt. nevmax) then
      call warnings ('proc_memb: maximum number of events reached')
      iev = nevmax
      nev = nevmax
   end if
   
   if (verbose_screen) then
      write (msg,'(a,i3)') 'proc_memb: event ', iev
      call fyi (trim(msg))
   end if
      
   return
   
end subroutine proc_memb


!***********************************************************************************************************************************
subroutine proc_nsmd (command_help)

! Check NEIC station metadata file for missing station codes
! Can be toggled by using the command without an argument, or set explicitly with arguments "on" or "off".

   implicit none

   include 'mloc.inc'

   character(len=132) :: nsmdfile, msg
   logical :: command_help, ex

   if (command_help) then
      write (*,'(10(a/))')&
       '"nsmd" (NEIC station metadata) determines if the NEIC station metadata file will',&
       'be searched for missing station codes after the master station file and any',&
       'supplemental station files have been searched. Any station codes that are found',&
       'in the NEIC metadata file will NOT be used in the current run; a supplemental',&
       'station file must be created and referenced in a subsequent run. The body of a',&
       'supplemental station file in "NEIC" format (isstn=5) for any matching station codes',&
       'will be found in the .log file. If there are multiple instances of a station code',&
       'with different coordinates they are all listed. Default is no search. Logical state',&
       'is toggled by issuing the command without an argument, or set explicitly with "on"',&
       'or "off" as the argument.'
      write (*,'(a,l1/)') 'Current value: ', nsmd
      return
   end if

   if (params(1:1) .eq. ' ') then
      nsmd = .not.nsmd
   else if (params(1:2) .eq. 'on') then
      nsmd = .true.
   else if (params(1:3) .eq. 'off') then
      nsmd = .false.
   else
      call warnings ('proc_nsmd: illegal argument')
   end if

   if (verbose_screen) then
      if (nsmd) then
         call fyi ('proc_nsmd: NEIC station metadata will be searched for missing station codes')
      else
         call fyi ('proc_nsmd: NEIC station metadata will not be searched for missing station codes')
      end if
   end if
   
   ! Check that the NEIC station metadata file is present
   if (nsmd) then
      nsmdfile = trim(station_path)//dirsym//'neic_stn.dat'
      inquire (file=nsmdfile,exist=ex)
      if (.not.ex) then
         write (msg, '(3a)') 'proc_nsmd: NEIC station metadata file (', trim(nsmdfile), ') does not exist'
         call warnings (trim(msg))
         nsmd = .false.
      end if
   end if

   return

end subroutine proc_nsmd


!***********************************************************************************************************************************
subroutine proc_oldr (command_help)

! Output file of phase readings over a limited distance range.

   implicit none

   include 'mloc.inc'

   logical :: command_help
   character(len=132) :: msg

   if (command_help) then
      write (*,'(3(a/))')&
       '"oldr" (output limited distance range) specifies an epicentral distance range for which',&
       'an output file will be created to list all phase readings in that range. Two arguments are',&
       'required'
      return
   end if

   call decode (params, value, nvar)
   
   if (nvar .eq. 2) then
      dist1 = value(1)
      dist2 = value(2)
      if (dist1 .gt. dist2) then
         call warnings ('proc_oldr: incorrect order of arguments')
         oldr_out = .false.
      else if (dist1 .lt. 0. .or. dist2 .gt. 180.) then
         write (msg,'(a,2f10.3)') 'proc_oldr: invalid argument: ', dist1, dist2
         call warnings (trim(msg))
         oldr_out = .false.
      else
         oldr_out = .true.
      end if
      if (verbose_screen) then
         if (oldr_out) then
            write (msg,'(a,2f6.2)') 'proc_oldr: epicentral distance range: ', dist1, dist2
            call fyi (trim(msg))
         end if
      end if
   else
      call warnings ('proc_oldr: "oldr" requires two arguments')
   end if 

   return
   
end subroutine proc_oldr


!***********************************************************************************************************************************
subroutine proc_pert (command_help)

! Perturb the starting location of all events

   implicit none

   include 'mloc.inc'

   logical :: command_help
   character(len=132) :: msg

   if (command_help) then
      write (*,'(7(a/))')&
       '"pert" (perturb starting locations) is used to specify a perturbation in location',&
       'parameters that will be applied to all events in the cluster. Four arguments are',&
       'required, in the order latitude, longitude, depth, and origin time. Latitude and',&
       'longitude perturbations are in decimal degrees. This is used primarily in cases',&
       'where local data are to be used to estimate the hypocentroid and the starting',&
       'locations are biased enough to cause convergence problems, but it can also be used',&
       'to force at least one iteration when the program is converging immediately'
      return
   end if

   call decode (params, value, nvar)
   
   if (nvar .eq. 4) then
      hlatshift = value(1)
      hlonshift = value(2)
      hdepthshift = value(3)
      htimeshift = value(4)
      if (verbose_screen) then
         write (msg,'(a,4f8.3)') 'proc_pert: hypocentroid perturbed by ', hlatshift, hlonshift, hdepthshift, htimeshift
         call fyi (trim(msg))
      end if
   else
      call warnings ('proc_pert: "pert" requires four parameters')
   end if

   return
   
end subroutine proc_pert


!***********************************************************************************************************************************
subroutine proc_phid (command_help)

! Phase re-identification.

   implicit none

   include 'mloc.inc'

   logical :: command_help

   if (command_help) then
      write (*,'(3(a/))')&
       '"phid" (phase identification) is used to toggle phase re-identification, which is',&
       'done after the data are read in and once more after the first iteration. By default,',&
       'it is turned on.'
      return
   end if

   if (params(1:1) .eq. ' ') then
      phid = .not.phid
   else if (params(1:2) .eq. 'on') then
      phid = .true.
   else if (params(1:3) .eq. 'off') then
      phid = .false.
   else
      call warnings ('proc_phid: illegal argument')
   end if

   return
   
end subroutine proc_phid


!***********************************************************************************************************************************
subroutine proc_phyp (command_help)

! Only P arrivals used for hypocentroid.
! It also keeps S-P readings if direct calibration is being used.

   implicit none

   include 'mloc.inc'

   logical :: command_help

   if (command_help) then
      write (*,'(6(a/))')&
       '"phyp" (P hypocentroid) specifies that only P phases will be used to estimate the',&
       'hypocentroid. This is normally used with "hlim" to specify only arrivals between',&
       '30 and 90 degrees, which provides a consistent, minimally biased estimate of the',&
       'hypocentroid. If direct calibration is being used, S-P readings will be retained',&
       'as well. If no argument is given, the status is toggled. The arguments "off" and',&
       '"on" can be used to set the state explicitly. Default = "on".'
      return
   end if

   if (params(1:1) .eq. ' ') then
      ponly = .not.ponly
   else if (params(1:2) .eq. 'on') then
      ponly = .true.
   else if (params(1:3) .eq. 'off') then
      ponly = .false.
   else
      call warnings ('proc_phyp: illegal argument')
   end if

   if (ponly) then
      if (verbose_screen) call fyi ('proc_phyp: using only first-arriving P for hypocentroid')
   else
      if (verbose_screen) call fyi ('proc_phyp: using all phases for hypocentroid')
   end if

   return
   
end subroutine proc_phyp


!***********************************************************************************************************************************
subroutine proc_plot (iev, command_help)

! Selected events for plotting

   implicit none

   include 'mloc.inc'

   integer :: iplot, iev
   logical :: command_help
   character(len=132) :: msg

   if (command_help) then
      write (*,'(6(a/))')&
       '"plot" (selective plotting) specifies that the current event is selected for',&
       'plotting in a secondary GMT script (_selN.bash). All cluster events are always',&
       'plotted in the main GMT script (_base.bash). The command takes one argument,',&
       'which specifies which selected-event plots this event will be included in. Limit',&
       'is 9 selected-event plots, and an event can belong to more than one by using the',&
       '"plot" command several times before the following "memb" command.'
      return
   end if

   if (nev .gt. 0) then
      call decode (params, value, nvar)
      if (nvar .eq. 1) then
         iplot = nint(value(1))
         if (iplot .ge. 1 .and. iplot .le. nplot_max) then
            plot(iev,iplot) = .true.
            if (verbose_screen) then
               write (msg,'(a,i3,a,i2)') 'proc_plot: event ', iev, 'will be plotted in plot ', iplot
               call fyi (trim(msg))
            end if
         else
            call warnings ('proc_plot: invalid plot number')
         end if
      else
         call warnings ('proc_plot: "plot" requires one argument, a valid plot number')
      end if
   else
      call warnings ('proc_plot: "plot" can only be issued after a "memb" command')
   end if

   return
   
end subroutine proc_plot


!***********************************************************************************************************************************
subroutine proc_pltt (command_help)

! Create GMT script for travel time plot. There are nine plots that can be made, by giving
! the indices of the desired plots as parameters:
! 1 = tt1, Summary travel time plot, 0-180 deg (default if no arguments are given)
! 2 = tt2, Entire teleseismic P branch, reduced to the theoretical P time
! 3 = tt3, Around the PKP caustic
! 4 = tt4, Near-source data
! 5 = tt5, local data, 0-4.0 deg.
! 6 = tt6, local-regional data, 0-30 deg. This can be done with reduced velocity.
! 7 = tt7, Sg, Sb, Sn, and Lg, 0-15 deg. This can be done with reduced velocity.
! 8 = tt8, Summry of relative depth phases (pP-P and sP-P) for all events
! 9 = tt9, S-P times

   implicit none

   include 'mloc.inc'

   integer :: i
   logical :: command_help
!  character(len=1) :: answer
   character(len=132) :: msg

   if (command_help) then
      write (*,'(16(a/))')&
       '"pltt" specifies one or more TT vs distance plots to be made. There are nine plot',&
       'types available, which are specified by the indices 1-9 after the command. Up to nine',&
       'arguments can be given in a single instance of the command:',&
       '  1 = Summary TT plot, full distance range 0-180 degrees',&
       '  2 = Teleseismic P residuals, 16-120 degrees',&
       '  3 = PKP branches around the caustic',&
       '  4 = Near-source residuals',&
       '  5 = Local distance, 0-4.0 degrees',&
       '  6 = Local-regional distances, 0-30 degrees (reduced velocities)',&
       '  7 = Local and regional S phases (reduced velocities)',&
       '  8 = Summary plot of relative depth phases, pP-P and sP-P, for all events',&
       '  9 = S-P times',&
       'If no arguments are given, all plots are turned off. If a custom crustal model has',&
       'been specified for short ranges, the travel times will be calculated from it in the',&
       'applicable distance range. In all cases, travel times are calculated for the depth',&
       'of the hypocentroid.'
      return
   end if

   call decode (params, value, nvar)
   
   if (nvar .eq. 0) then
      tt1 = .false.
      tt2 = .false.
      tt3 = .false.
      tt4 = .false.
      tt5 = .false.
      tt6 = .false.
      tt7 = .false.
      tt8 = .false.
      tt9 = .false.
   else if (nvar .le. 9) then
      do i = 1,nvar
         if (nint(value(i)) .eq. 1) tt1 = .true.
         if (nint(value(i)) .eq. 2) tt2 = .true.
         if (nint(value(i)) .eq. 3) tt3 = .true.
         if (nint(value(i)) .eq. 4) tt4 = .true.
         if (nint(value(i)) .eq. 5) tt5 = .true.
         if (nint(value(i)) .eq. 6) tt6 = .true.
         if (nint(value(i)) .eq. 7) tt7 = .true.
         if (nint(value(i)) .eq. 8) tt8 = .true.
         if (nint(value(i)) .eq. 9) tt9 = .true.
      end do
!     if (tt6 .or. tt7) then
!        write (*,'(t3,a)',advance='no') 'Reduced velocities? (y or n): '
!        read (*,'(a)') answer
!        reduced = (answer .eq. 'y')
!     end if
   else
      call warnings ('proc_pltt: "pltt" takes a maximum of nine arguments')
   end if

   if (verbose_screen) then
      write (msg,'(a,9(l1,1x))') 'proc_pltt: travel time plots ', tt1, tt2, tt3, tt4, tt5, tt6, tt7, tt8, tt9
      call fyi (trim(msg))
   end if

   return
   
end subroutine proc_pltt


!***********************************************************************************************************************************
subroutine proc_ppri (command_help)

! Prevent phase re-identification for a specific phase.
! One argument, the phase name which should not be changed. This command can be issued
! multiple times.

   implicit none

   include 'mloc.inc'

   logical :: command_help

   if (command_help) then
      write (*,'(3(a/))')&
       '"ppri" specifies  phases that cannot be renamed. It takes one argument, a phase',&
       'name. The command can be issued multiple times, up to the limit specified by',&
       '"n_no_phreid_max".'
      return
   end if

   if (params(1:1) .ne. ' ') then
      n_no_phreid = n_no_phreid + 1
      if (n_no_phreid .le. n_no_phreid_max) then
         no_phreid(n_no_phreid) = trim(params)
         if (verbose_screen) call fyi ('proc_ppri: phase '//trim(no_phreid(n_no_phreid))//' will not be re-identified')
      else
         call warnings ('proc_ppri: maximum number of instances reached')
      end if
   else
      call warnings ('proc_ppri: "ppri" requires a single argument')
   end if

   return
   
end subroutine proc_ppri


!***********************************************************************************************************************************
subroutine proc_pttt (command_help)

! Perfect theoretical travel times for the hypocentroid

   implicit none

   include 'mloc.inc'

   logical :: command_help

   if (command_help) then
      write (*,'(3(a/))')&
       '"pttt" (perfect theoretical travel times) allows the variance of the travel time,',&
       'model (ttsprd) to be set to zero ("on") for the purpose of calculating the hypocentroid',&
       'and its uncertainty. Phase spread values (either default or read from a .ttsprd file)',&
       'are still used in the windowing algorithm. Issuing the command with no argument toggles',&
       'the current state. Default is "off".'
      return
   end if

   if (params(1:1) .eq. ' ') then
      pttt = .not.pttt
   else if (params(1:2) .eq. 'on') then
      pttt = .true.
   else if (params(1:3) .eq. 'off') then
      pttt = .false.
   else
      call warnings ('proc_pttt: illegal argument')
   end if

   if (pttt) then
      pttt_pr = 'on '
   else
      pttt_pr = 'off'
   end if
   if (verbose_screen) call fyi ('proc_pttt: assumption of perfect travel times is '//data_weight_pr)

   return
   
end subroutine proc_pttt


!***********************************************************************************************************************************
subroutine proc_puke (command_help)

! PUKE format output file (derived from Engdahl's PICK format).
! One of several output file options to provide the results in useful ways

   implicit none

   include 'mloc.inc'

   logical :: command_help

   if (command_help) then
      write (*,'(5(a/))')&
       '"puke" (PUKE formatted output file) specifies that a .puke output will be written.',&
       'This is a close cousin of the "PICK" format that Engdahl sometimes uses to transfer',&
       'relocation results to other researchers. If no argument is given, the status is',&
       'toggled. The arguments "off" and "on" can be used to set the state explicitly.'
      return
   end if

   if (params(1:1) .eq. ' ') then
      pukeout = .not.pukeout
   else if (params(1:2) .eq. 'on') then
      pukeout = .true.
   else if (params(1:3) .eq. 'off') then
      pukeout = .false.
   else
      call warnings ('proc_puke: illegal argument')
   end if

   if (verbose_screen) then
      if (pukeout) then
         call fyi ('proc_puke: will create a PUKE file')
      else
         call fyi ('proc_puke: will not create a PUKE file')
      end if
   end if

   return
   
end subroutine proc_puke


!***********************************************************************************************************************************
subroutine proc_radf (command_help)

! Read agency and deployment fields from station files and event data files.

   implicit none

   include 'mloc.inc'

   logical :: command_help

   if (command_help) then
      write (*,'(4(a/))')&
       '"radf" (Read Agency and Deployment Fields) specifies that agency and deployment fields',&
       'in event data files (.mnf files) and station lists will be read and used to resolve',&
       'station code conflicts. If no argument is given, the status is toggled. The arguments',&
       '"off" and "on" can be used to set the state explicitly. Default is "off".'
      return
   end if

   if (params(1:1) .eq. ' ') then
      read_ad = .not.read_ad
   else if (params(1:2) .eq. 'on') then
      read_ad = .true.
   else if (params(1:3) .eq. 'off') then
      read_ad = .false.
   else
      call warnings ('proc_radf: illegal argument')
   end if

   if (verbose_screen) then
      if (read_ad) then
         call fyi ('proc_radf: will read and use agency and deployment fields')
      else
         call fyi ('proc_radf: will not read agency and deployment fields')
      end if
   end if

   return
   
end subroutine proc_radf


!*****************************************************************************************
subroutine proc_rdpp (command_help)

! Makes a plot of the relative depth phase data for a single event. May be issued multiple
! times. Plots are put into a subdirectory.

   implicit none

   include 'mloc.inc'

   logical :: command_help
   character(len=132) :: msg

   if (command_help) then
      write (*,'(3(a/))')&
       '"rdpp" makes a plot of relative depth phase (pP-P, sP-P) data for a single event.',&
       'The command takes a single argument, the event name (as given in the "even" command.',&
       'The command can be issued multiple times, up to the limit specified by "n_rdpp_max".'
      return
   end if

   if (params(1:1) .ne. ' ') then
      if (len_trim(params) .le. 30) then
         if (verbose_screen) call fyi ('proc_rdpp: relative depth phase plot for event '//trim(params))
         rdpp = .true.
         n_rdpp = n_rdpp + 1
         if (n_rdpp .le. n_rdpp_max) then
            rdpp_evt(n_rdpp) = trim(params)
         else
            call warnings ('proc_rdpp: maximum number of instances reached')
            n_rdpp = n_rdpp_max
         end if
      else
         msg = 'proc_rdpp: event name must be 30 characters or less ('//trim(params)//')'
         call warnings (trim(msg))
      end if
   else
      call warnings ('proc_rdpp: "rdpp" requires one argument, an event name')
   end if

   return
   
end subroutine proc_rdpp


!***********************************************************************************************************************************
subroutine proc_rels (command_help)

! Reading errors for local stations
! This command sets a common value for local distance phases within a specified distance.
! Do not confuse with command "mare" which only sets a lower limit for empirical reading errors.

   implicit none

   include 'mloc.inc'

   logical :: command_help
   character(len=132) :: msg

   if (command_help) then
      write (*,'(5(a/))')&
       '"rels" (reading error local stations) specifies the reading errors that will be used',&
       'for direct-arriving phases Pg, Pb, Sg, and Sb within a specified distance. This is',&
       'mainly used for single event location of events with local data, e.g., calibration',&
       'events. Three arguments are required, the reading error for P and S phases, respectively,',&
       'and the distance range (epicentral degrees) in which these values apply.'
      write (*,'(/a,2(f5.2,a),f4.2,a)') 'Current values: ', rderr_loc_p, ' (P) /', rderr_loc_s, '(S) /', rderr_loc_delt, ' deg'
      return
   end if

   call decode (params, value, nvar)
   
   if (nvar .eq. 3) then
      rels_set = .true.
      rderr_loc_p = value(1)
      rderr_loc_s = value(2)
      rderr_loc_delt = value(3)
      if (verbose_screen) then
         write (msg,'(a,2f6.2)') 'proc_rels: ', rderr_loc_p, rderr_loc_s, rderr_loc_delt
         call fyi (trim(msg))
      end if
   else
      call warnings ('proc_rels: "rels" requires three arguments')
   end if 

   return
   
end subroutine proc_rels


!***********************************************************************************************************************************
subroutine proc_revi (command_help)

! Review of currently-selected parameters and options

   implicit none

   include 'mloc.inc'

   integer :: jev, i, j
   logical :: command_help

   if (command_help) then
      write (*,'(3(a/))')&
       '"revi" (review) is used to provide detailed information about the current events',&
       'in the cluster and the settings that influence how the relocation will be done. It',&
       'takes no arguments.'
      return
   end if

   write (*,'(t3,a,i3,a,i3,a)') 'Number of events in cluster = ', nev, '  (', ncal(1), ' calibration events)'
   if (nev .ge. 1) then
      do jev = 1,nev
         write (*,'(/t3,a,i3,2(2x,a))') 'Event ', jev, trim(evtnam(jev)), trim(infile(jev))
         write (*,'(t3,8a)') '(frec) lat/lon/depth/OT: ', fixpr(jev,1),'/',fixpr(jev,2),'/',fixpr(jev,3),'/',fixpr(jev,4)
         if (cal_event(jev,1)) then
            write (*,'(t3,a,3f10.3,2i4,f10.3)') '(cali) ', cal_lat(jev,1), cal_lon(jev,1), cal_dep(jev,1), cal_hr(jev,1),&
             cal_min(jev,1), cal_sec(jev,1)
         end if
      end do
   end if
   write (*,'(/t3,a/)') 'Parameters common to all events (alphabetic by command name):'     
   if (bdp_list) write (*,'(t3,a,t40,2a)') '(bdps) bogus depth phase stations', ': ', trim(bdp_filnam)
   write (*,'(t3,a,t40,a,l1)') '(bias) hypocentroid bias-correction', ': ', bias_corr
   write (*,'(t3,a,t40,a,l1)') '(bloc) BayesLoc file output', ': ', blocout
   write (*,'(t3,a,t40,a,l1)') '(bptc) bounce point topography correction', ': ', bptc
   write (*,'(t3,a,t40,a,l1)') '(ccat) COMCAT file output', ': ', comcatout
   write (*,'(t3,a,t40,a,3(2f6.1,2x))') '(clim) cluster distance limits', ': ', (clim(i,1),clim(i,2),i=1,nclim)
   select case (tt_corr)
      case (0)
         write (*,'(t3,a,t40,a)') '(corr) station corrections', ': none'
      case (1)
         write (*,'(t3,a,t40,a)') '(corr) station corrections', ': elevation corrections'
   end select
   if (calibration) then
      select case (icaltype)
         case (1)
            write (*,'(t3,a,t40,a)') '(ctyp) indirect calibration type', ': traditional'
         case (2)
            write (*,'(t3,a,t40,a)') '(ctyp) indirect calibration type', ': systematic'
         case (3)
            write (*,'(t3,a,t40,a)') '(ctyp) indirect calibration type', ': optimal'
      end select
   end if
   write (*,'(t3,a,t40,a,f6.2)') '(cvff) radius', ': ', radius_cvff
   write (*,'(t3,a,t40,a,l1)') '(cvou) covariance matrix output', ': ', cv_out
   write (*,'(t3,a,t40,a,l1)') '(datf) DATF file output', ': ', datfout
   write (*,'(t3,a,t40,2a)') '(dcal) direct calibration', ': ', dcal_pr
   if (plot_gina) write (*,'(t3,a,t40,a)') '(dem1) topography/bathymetry', ': GINA DEM'
   if (plot_globe) write (*,'(t3,a,t40,a)') '(dem1) topography', ': GLOBE DEM'
   if (plot_etopo1) write (*,'(t3,a,t40,a)') '(dem1) topography/bathymetry', ': ETOPO1 DEM'
   if (plot_dem2) write (*,'(t3,a,t40,2a)') '(dem2) high-rez topography', ': ', trim(dem2_filename)
   if (diffdat) write (*,'(t3,a,t40,2a)') '(diff) differential times data', ': ', trim(diffdatfilnam)
   if (epa_plot) then
      do i = 1,n_epa_plot
         write (*,'(t3,a,t40,2a,f5.1,1x,l1)') '(epap) emirical path anomaly plot', ': ',&
          epa_plot_phase(i), epa_plot_distance(i), epa_plot_raypath(i)
      end do
   end if
   write (*,'(t3,a,t40,a,l1)') '(eplt) confidence ellipse plot', ': ', eplt
   write (*,'(t3,a,t40,a,l1)') '(fdhp) focal depth histogram', ': ', fdhp
   write (*,'(t3,a,t40,2a)') '(flag) use data flags', ': ', dflag_pr
   if (fault_map) then
      do i = 1,n_fault_map
         write (*,'(t3,a,t40,2a)') '(fmap) fault map file', ': ', trim(fault_map_filename(i))
         write (*,'(t3,a,t40,2a)') '(fplt) fault plotting options',': ', psxy_wsf(i)
      end do
   end if
   write (*,'(t3,a,t40,8a)') '(freh) hypocentroid lat/lon/depth/OT', ': ', fixprh(1), '/', fixprh(2), '/', fixprh(3),&
    '/', fixprh(4)
   write (*,'(t3,a,t40,a,3(2f6.1,2x))') '(hlim) hypocentroid limits', ': ', (hlim(i,1),hlim(i,2),i=1,nhlim)
   write (*,'(t3,a,t40,a,2f6.2)') '(lgtt) Lg travel times', ': ', lg_a, lg_b
   if (locmod) write (*,'(t3,a,t40,2a)') '(lmod) local velocity model', ': ', trim(locmodfname)
   write (*,'(t3,a,t40,a,i1)') '(lonr) center of longitude range', ':', longitude_range
   if (lresout) then
      write (*,'(t3,a,t40,a,l1,a,f5.2,a)') '(lres) LRES file output',': ', lresout, ' (LRES = ', lres, ')'
   else
      write (*,'(t3,a,t40,a,l1)') '(lres) LRES file output',': ', lresout
   end if
   write (*,'(t3,a,t40,3(a,f5.2),a)') '(mare) minimum reading errors', ': ', rderr_min_loc, ' (local) /',&
    rderr_min, ' (general) /', rderr_min_depth, ' (depth phases)'
   write (*,'(t3,a,t40,a,l1)') '(mdou) map_dat output file', ': ', md_out
   write (*,'(t3,a,t40,a,l1)') '(nsmd) NEIC station metadata search', ': ', nsmd
   if (oldr_out) then
      write (*,'(t3,a,t40,a,l1,2f7.2)') '(oldr) output limited distance range', ': ', oldr_out, dist1, dist2
   else
      write (*,'(t3,a,t40,a,l1)') '(oldr) output limited distance range', ': ', oldr_out
   end if
   write (*,'(t3,a,t40,a,4f7.2)') '(pert) perturbation', ': ', hlatshift, hlonshift, hdepthshift, htimeshift
   write (*,'(t3,a,t40,a,l1)') '(phid) phase re-identification', ': ', phid
   write (*,'(t3,a,t40,a,l1)') '(phyp) only P for hypocentroid', ': ', ponly
   if (n_no_phreid .gt. 0) then
      do i = 1,n_no_phreid
         write (*,'(t3,a,t40,2a)') '(ppri) phase renaming prevented for', ': ', no_phreid(i)
      end do
   end if
   write (*,'(t3,a,t40,2a)') '(pttt) perfect theoretical TTs', ': ', pttt_pr
   write (*,'(t3,a,t40,a,l1)') '(puke) PUKE file output', ': ', pukeout
   write (*,'(t3,a,t40,a,l1)') '(radf) Read agency and deployment fields', ': ', read_ad
   if (rdpp) then
      do i = 1,n_rdpp
         write (*,'(t3,a,t40,a)') '(rdpp) relative depth phase plot for ', trim(rdpp_evt(i))
      end do
   end if
   if (rels_set) then
      write (*,'(t3,a,t40,a,3(f5.2,a))') '(rels) reading errors local phases', ': ', rderr_loc_p,&
       ' (P) / ', rderr_loc_s, ' (S) /', rderr_loc_delt, ' (deg)'
   else
      write (*,'(t3,a,t40,a)') '(rels) reading errors local phases', ': not set'
   end if
   write (*,'(t3,a,t40,2a)') '(rfil) reading error file', ': ', trim(rderrfname)
   if (read_hdf) write (*,'(t3,a,t40,2a)') '(rhdf) starting locations', ': ', trim(rhdf_filnam)
   write (*,'(t3,a,t40,2f5.2)') '(secv) station elevation correction vel', ': ', pcrvel, scrvel
   write (*,'(t3,a,t40,a,f6.3,2f4.1)') '(shcl) hypocentroid convergence',': ', cl_epi_h, cl_dep_h, cl_ot_h
   if (n_skip .gt. 0) then
      do i = 1,n_skip
         write (*,'(t3,a,t40,a,1x,a,1x,a)') '(skip) skipping station/phase/author', ': ', (skip_params(i,j),j=1,3)
      end do
   end if
   write (*,'(t3,a,t40,a,l1)') '(splt) seismicity map', ': ', splt
   if (suppstn) then
      do i = 1,n_supp_stn_file
         write (*,'(t3,a,t40,2a)') '(sstn) supplemental station file', ': ', trim(suppfilnam(i))
      end do
   end if
   if (subc_set) then
      write (*,'(t3,a,t40,a,f6.2,2(1x,i3))') '(subc) subcluster selection', ':', subc_delt, subc_nmin, subc_nconnect
   else
      write (*,'(t3,a,t40,a)') '(subc) subcluster selection', ': not set'
   end if
   write (*,'(t3,a,t40,2a)') '(taup) global TT model', ': ', trim(taup_model)
   write (*,'(t3,a,t40,a,f6.2)') '(tikh) Tikhonov regularization', ': ', tikhonov_factor
   write (*,'(t3,a,t40,a,2f6.2)') '(tptt) T-phase travel times', ': ', tphase_a, tphase_b
   write (*,'(t3,a,t40,2a)') '(tfil) TT spread file', ': ', trim(ttsprdfname)
   if (tt5e) then
      do i = 1,n_tt5e
         write (*,'(t3,a,t40,a)') '(tt5e) single event tt5 plot for ', trim(tt5e_evt(i))
      end do
   end if
   if (tt5s) then
      do i = 1,n_tt5s
         write (*,'(t3,a,t40,a)') '(tt5s) single station tt5 plot for ', trim(tt5s_sta(i))
      end do
   end if
   if (ttou) then
      do i = 1,n_ttou
         write (*,'(t3,a,t40,a)') '(ttou) empirical TTs for ', trim(ttou_phase(i))
      end do
   end if
   write (*,'(t3,a,t40,a,4l1)') '(vect) vectors plotted', ': ', vectors(1), vectors(2), vectors(3), vectors(4)
   write (*,'(t3,a,t40,a,l1)') '(vlog) verbose logging mode', ': ', verbose_log
   write (*,'(t3,a,t40,a,l1)') '(vscr) verbose screen mode', ': ', verbose_screen
   write (*,'(t3,a,t40,2a)') '(weig) weight residuals', ': ', data_weight_pr
   write (*,'(t3,a,t40,2(a,f4.1))') '(wind) TT window', ': ', wind1, '/', wind2
   write (*,'(/t3,a/)') 'Continue with commands:'

   return
   
end subroutine proc_revi


!***********************************************************************************************************************************
subroutine proc_rfil (command_help)

! Empirical reading error file. 

   implicit none

   include 'mloc.inc'

   logical :: ex, command_help

   if (command_help) then
      write (*,'(4(a/))')&
      '"rfil" specifies a .rderr file to be opened for reading the empirically-derived',&
      'reading errors of specific station-phases. The file is expected to be in the',&
      'data directory, so the argument of the command is just the filename, not the',&
      'full pathname. Revert to the default values by entering the argument "default".'
      return
   end if

   rderrfname = trim(datadir)//dirsym//trim(params)
   inquire (file=rderrfname,exist=ex)
   if (ex) then
      read_rderr = .true.
      if (verbose_screen) call fyi ('proc_rfil: using '//trim(rderrfname))
   else
      if (trim(params) .eq. 'default') then
         rderrfname = 'default'
         read_rderr = .false.
      else
         call warnings ('proc_rfil: file '//trim(rderrfname)//' does not exist')
      end if
   end if

   return
   
end subroutine proc_rfil


!***********************************************************************************************************************************
subroutine proc_rhdf (command_help)

! Starting locations from an HDF file.

   implicit none

   include 'mloc.inc'

   logical :: ex, command_help

   if (command_help) then
      write (*,'(6(a/))')&
       '"rhdf" specifies an HDF file that will be read to obtain starting locations. The',&
       'hypocentral values read from the HDF file are over-ridden by any commands "lat", "long",',&
       '"time", or "depX" issued later in the command file. If any cluster events missing',&
       'from the HDF file, starting locations for those events will revert to those of the',&
       'corresponding input file. The command takes one argument, the name of the HDF file,',&
       'which must be stored in the data directory with the input files.'
      return
   end if

   rhdf_filnam = trim(datadir)//dirsym//trim(params)
   inquire (file=rhdf_filnam,exist=ex)
   if (ex) then
      read_hdf = .true.
      if (verbose_screen) call fyi ('proc_rhdf: using: '//trim(rhdf_filnam))
   else
      call warnings ('proc_rhdf: file '//trim(rhdf_filnam)//' does not exist')
   end if

   return
   
end subroutine proc_rhdf


!***********************************************************************************************************************************
subroutine proc_run (cmndfil, command_help)

! Begin the relocation

   implicit none

   include 'mloc.inc'

   logical :: cmndfil, command_help
   character(len=132) :: msg

   if (command_help) then
      write (*,'(4(a/))')&
       '"run" ends the command processing phase and starts the HD relocation process.',&
       'There are a few additional queries that provide important controls over the',&
       'relocation process, but no further opportunities to add events to the cluster or',&
       'issue any of the standard commands.'
      return
   end if

   if (cmndfil) close (io_cfil) ! In case the run command is given in a command file

   return
   
end subroutine proc_run


!***********************************************************************************************************************************
subroutine proc_secv (command_help)

! Set velocities for station elevation corrections

   implicit none

   include 'mloc.inc'

   logical :: command_help
   character(len=132) :: msg

   if (command_help) then
      write (*,'(5(a/))')&
       '"secv" (station elevation correction velocities) specifies the velocities that will be used',&
       'for station elevation corrections. The command takes two arguments, the velocity for P and S.'
      write (*,'(/a,2(f5.2,a))') 'Current values: ', pcrvel, ' (P) /', scrvel, ' (S)'
      return
   end if

   call decode (params, value, nvar)
   
   if (nvar .eq. 2) then
      pcrvel = value(1)
      scrvel = value(2)
      if (verbose_screen) then
         write (msg,'(a,2f6.2)') 'proc_secv: ', pcrvel, scrvel
         call fyi (trim(msg))
      end if
   else
      call warnings ('proc_secv: "secv" requires two arguments')
   end if 

   return
   
end subroutine proc_secv


!***********************************************************************************************************************************
subroutine proc_shcl (command_help)

! Set new convergence limits for the hypocentroid components.

   implicit none
   
   include 'mloc.inc'

   logical :: command_help
   character(len=132) :: msg

   if (command_help) then
      write (*,'(14(a/))')&
       '"shcl" (set hypocentroid convergence limits) modifies the criteria used to decide if',&
       'convergence has been reached, i.e., the change in each parameter from one iteration',&
       'to the next is smaller than the convergence limit. This command only deals with',&
       'the limits relating to the hypocentroid (origin time, epicenter and focal depth).',&
       'Actual convergence requires all parameters for the hypocentroid and each cluster vector',&
       'to satisfy their corresponding convergence criteria, but it is often the case that',&
       'convergence is prevented by oscillations in the hypocentroid origin time, sometimes',&
       'coupled with instability in one of the epicentral parameters. The "epicenter"',&
       'limit is applied separately to latitude and longitude. It should be noted that achieving',&
       'convergence by increasing the limits may not produce a reliable solution, but there are',&
       'advantages to having a converged solution for investigating why it was necessary.',&
       ' ',&
       'The command takes three arguments, convergence limits for epicenter (degrees), focal',&
       'depth (km) and origin time (sec).'
      write (*,'(/a,f6.3,a,f4.1,a,f4.1,a)') 'Current values: ', cl_epi_h, ' (deg), ', cl_dep_h, '(km), ', cl_ot_h, ' (sec)'
      return
   end if

   call decode (params, value, nvar)
   
   if (nvar .eq. 3) then
      cl_epi_h = value(1) ! Hypocentroid lat-lon, degrees
      cl_dep_h = value(2) ! Hypocentroid depth, km
      cl_ot_h = value(3) ! Hypocentroid OT, sec
      if (verbose_screen) then
         write (msg,'(a,3f6.3)') 'proc_shcl: ', cl_epi_h, cl_dep_h, cl_ot_h
         call fyi (trim(msg))
      end if
   else
      call warnings ('proc_shcl: "shcl" requires three arguments')
   end if 


   return
   
end subroutine proc_shcl


!***********************************************************************************************************************************
subroutine proc_skip (command_help)

! Specify readings that will be skipped on the basis of station, phase or author. This version of SKIP allows specification of all
! three parameters, or wild cards ("*") can be used to specify all instances of the parameter. Readings with blank phase name can be
! skipped by specifying the phase code as "blank" (no quotes).

   implicit none

   include 'mloc.inc'

   logical :: command_help
   character(len=8) :: param1, param2, param3, param4
   character(len=132) :: msg

   if (command_help) then
      write (*,'(11(a/))')&
       '"skip" specifies a triplet of station-phase-author for which all readings are flagged',&
       '(with "s"). Skipped readings will not to be used in the relocation for cluster',&
       'vectors or hypocentroid. The command takes three arguments, a station code, a phase',&
       'code and an author code. Wildcards ("*") are supported:',&
       '  "skip GRMI * *" means skip all readings from GRMI',&
       '  "skip GRMI Pg *" means skip all Pg readings from GRMI',&
       '  "skip GRMI Pg ARGhods" means skip all Pg readings from GRMI by the author "ARGhods"',&
       '  "skip * * ARGhods" means skip all readings from the author "ARGhods"',&
       'Readings with blank phase code can be skipped by giving "blank" (without quotes) as',&
       'phase code. The command can be issued multiple times, up to the limit specified by',&
       '"n_skip_max".'
      return
   end if

   if (params(1:1) .ne. ' ') then
      call decode2 (params, nvar, param1, param2, param3, param4)
      if (nvar .ne. 3) then
         write (msg,'(a,i1)') 'proc_skip: incorrect number of arguments = ', nvar
         call warnings (trim(msg))
         return
      end if
      if (verbose_screen) then
         call fyi ('proc_skip: station skip code = '//param1)
         call fyi ('proc_skip: phase skip code   = '//param2)
         call fyi ('proc_skip: author skip code  = '//param3)
      end if
      skip = .true.
      n_skip = n_skip + 1
      if (n_skip .le. n_skip_max) then
         skip_params(n_skip,1) = param1
         skip_params(n_skip,2) = param2
         skip_params(n_skip,3) = param3
      else
         call warnings ('proc_skip: maximum number of instances reached')
         n_skip = n_skip_max
      end if
   else
      call warnings ('proc_skip: "skip" requires three arguments')
   end if

   return
   
end subroutine proc_skip


!***********************************************************************************************************************************
subroutine proc_splt (command_help)

! Determines if a "seismicity" plot is made.
! Can be toggled by using the command without an argument, or set explicitly with arguments "on" or "off".

   implicit none

   include 'mloc.inc'

   logical :: command_help

   if (command_help) then
      write (*,'(5(a/))')&
       '"splt" (seismicity plot) specifies whether a "seismicity map" will be made.',&
       'This type of plot indicates event locations by open circles of 1 km diameter',&
       'and relocation vectors, but does not carry event numbers. Issuing the command',&
       'without arguments toggles the current state. The arguments "on" and "off" may be',&
       'used to set the state explicitly. Default is "off".'
      return
   end if

   if (params(1:1) .eq. ' ') then
      splt = .not.splt
   else if (params(1:2) .eq. 'on') then
      splt = .true.
   else if (params(1:3) .eq. 'off') then
      splt = .false.
   else
      call warnings ('proc_splt: illegal argument')
   end if

   if (splt) then
      if (verbose_screen) call fyi ('proc_splt: a seismicity map will be made')
   else
      if (verbose_screen) call fyi ('proc_splt: a seismicity map will not be made')
   end if

   return
   
end subroutine proc_splt


!***********************************************************************************************************************************
subroutine proc_spou (command_help)

! Output file of S-P data.

   implicit none

   include 'mloc.inc'

   logical :: command_help

   if (command_help) then
      write (*,'(4(a/))')&
       '"spou" (S-P output file) specifies that a special file',&
       'containing S-P data will be written. It will have the filename suffix "sp".',&
       'If no argument is given, the status is toggled. The arguments "off" and "on" can',&
       'be used to set the state explicitly.'
      return
   end if

   if (params(1:1) .eq. ' ') then
      sp_out = .not.sp_out
   else if (params(1:2) .eq. 'on') then
      sp_out = .true.
   else if (params(1:3) .eq. 'off') then
      sp_out = .false.
   else
      call warnings ('proc_spou: illegal argument')
   end if

   if (sp_out) then
      if (verbose_screen) call fyi ('proc_spou: will create a .sp file')
   else
      if (verbose_screen) call fyi ('proc_spou: will not create a .sp file')
   end if

   return
   
end subroutine proc_spou


!***********************************************************************************************************************************
subroutine proc_sstn (command_help)

! Supplementary station file

   implicit none

   include 'mloc.inc'

   logical :: ex, command_help

   if (command_help) then
      write (*,'(15(a/))')&
       '"sstn" (supplemental station data) specifies the pathname of a file containing',&
       'coordinates for stations not in the main station file. These files are often',&
       'stored in the "tables/stn/" subdirectory, but can be stored elsewhere, including',&
       'the cluster data directory. Several formats are supported, specified by an integer',&
       'in the first column of the first line:',&
       '  0 - New master station file format, from April 28, 2014',&
       '  1 - ISC FFB format, (deg-min-sec*10)',&
       '  2 - SEISAN format, (deg-decimal min)',&
       '  3 - simplified "mloc" format, decimal degrees, geographic coordinates',&
       '  4 - China Seismic Bureau format',&
       '  5 - NEIC format',&
       '  6 - MSU format (used by Kevin Mackey)',&
       '  9 - Former master station file format (used index "0")',&
       'There is a support document describing the formats in detail. The command can be',&
       'called multiple times, up to the value of n_supp_stn_file_max.'
      return
   end if

   if (params(1:1) .eq. ' ') then
      call warnings ('proc_sstn: "sstn" requires a pathname as argument')
   else
      if (n_supp_stn_file .lt. n_supp_stn_file_max) then
         inquire (file=trim(params),exist=ex)
         if (ex) then
            if (index(trim(params),trim(station_master)) .gt. 0) then
               call warnings ('proc_sstn: master station file cannot be called from this command')
               return
            end if
            suppstn = .true.
            n_supp_stn_file = n_supp_stn_file + 1
            suppfilnam(n_supp_stn_file) = trim(params)
            if (verbose_screen) call fyi ('proc_sstn: supplementary station file: '//trim(suppfilnam(n_supp_stn_file)))
         else
            call warnings ('proc_sstn: file '//trim(params)//' does not exist')
         end if
      else
         call warnings ('proc_sstn: maximum number of supplemental station files reached')
      end if
   end if

   return
   
end subroutine proc_sstn


!***********************************************************************************************************************************
subroutine proc_star (iev, command_help)

! Plot a star at the location of a specific event.

   implicit none

   include 'mloc.inc'

   integer :: iev
   logical :: command_help
   character(len=132) :: msg

   if (command_help) then
      write (*,'(5(a/))')&
       '"star" (star plot) is used to specify the plotting of a solid red star at the location',&
       'of the current event in the map plots by GMT. This is normally used to highlight',&
       'events of special interest, such as a mainshock and major aftershocks. The command takes one',&
       'argument, the size of the circumscribing circle (recommended sizes are in the range 0.2-0.5).',&
       'The command can be issued up to 10 times.'
      return
   end if

   call decode (params, value, nvar)
   
   if (nvar .eq. 1) then
      if (n_star .lt. n_star_nmax) then
         star_plot = .true.
         n_star = n_star + 1
         iev_star(n_star) = iev
         star_size(n_star) = value(1)
         if (verbose_screen) then
            write (msg,'(a,i3,f8.3)') 'proc_star: ', iev_star(n_star), star_size(n_star)
            call fyi (trim(msg))
         end if
      else
         call warnings ('proc_star: maximum number of stars reached')
      end if
   else
      call warnings ('proc_star: "star" requires one parameter')
   end if

   return
   
end subroutine proc_star


!***********************************************************************************************************************************
subroutine proc_stat (command_help)

! Plot a station location as a triangle.

   implicit none

   include 'mloc.inc'

   integer :: i
   logical :: command_help
   character(len=132) :: msg

   if (command_help) then
      write (*,'(6(a/))')&
       '"stat" (station plot) is used to specify the parameters of a triangle that will be',&
       'plotted in the map plots by GMT. The command takes three arguments:',&
       '   latitude of the center point',&
       '   longitude of the center point',&
       '   size of the circumscribing circle (0.30 is a good choice)',&
       'The command can be issued up to 30 times.'
      return
   end if

   call decode (params, value, nvar)
   
   if (nvar .eq. 3) then
      if (n_stat .lt. n_stat_nmax) then
         stat_plot = .true.
         n_stat = n_stat + 1
         do i = 1,3
            stat_data(n_stat,i) = value(i)
         end do
         call set_longitude_range (stat_data(n_stat,2), longitude_range)
         if (verbose_screen) then
            write (msg,'(a,3f8.3)') 'proc_stat: ', stat_data(n_stat,1), stat_data(n_stat,2), stat_data(n_stat,3)
            call fyi (trim(msg))
         end if
      else
         call warnings ('proc_stat: maximum number of stations for plotting reached')
      end if
   else
      call warnings ('proc_stat: "stat" requires three parameters')
   end if

   return
   
end subroutine proc_stat


!***********************************************************************************************************************************
subroutine proc_step (command_help)

! Stop after a given number of iterations.
! "step 0" will not run any inversions.

   implicit none

   include 'mloc.inc'

   logical :: command_help
   character(len=132) :: msg

   if (command_help) then
      write (*,'(5(a/))')&
       '"step" specifies the number of iterations to run. It takes one argument, an integer',&
       'between 0 and the maximum number of iterations allowed (4). If the convergence',&
       'criteria are met before the specified number of iterations, the program ends as',&
       'usual. "step 0" is equivalent to the old "fwd" command. The residuals to the',&
       'starting locations will be calculated and the program will exit. No inversion or',&
       'iteration will be done'
      return
   end if

   call decode (params, value, nvar)
   
   if (nvar .eq. 1) then
      if (int(value(1)) .ge. 0 .and. int(value(1)) .le. itmax) then
         nstep = int(value(1))
         if (verbose_screen) then
            write (msg,'(a,i1,a)') 'proc_step: will execute ', nstep, ' iterations'
            call fyi (trim(msg))
         end if
      else
         write (msg,'(a,i1)') 'proc_step: number of steps must be between 0 and ', itmax
         call warnings (trim(msg))
      end if
   else
      call warnings ('proc_step: "step" takes one argument')
   end if

   return
   
end subroutine proc_step


!***********************************************************************************************************************************
subroutine proc_stop (command_help)

! Close open files and stop the run

   implicit none

   include 'mloc.inc'

   logical :: op, command_help

   if (command_help) then
      write (*,'(4(a/))')&
       '"stop" aborts the current run while still in the command processing phase. It is',&
       'normally used when some aspect of the reading of a command file (or interactive',&
       'input) has gone wrong. It is not needed in normal operation, and cannot be issued',&
       'after the "run" command has been given.'
      return
   end if

   inquire (io_in,opened=op)
   if (op) close (io_in)
   inquire (io_cfil,opened=op)
   if (op) close (io_cfil)
   inquire (io_out,opened=op)
   if (op) close (io_out)
   inquire (io_taup,opened=op)
   if (op) close (io_taup)

   stop
   
end subroutine proc_stop


!*****************************************************************************************
subroutine proc_subc (command_help)

! Subcluster.
! Used to select events that meet given requirements about the number of readings within a
! certain distance and connectivity with other events in the cluster, and creates the main
! body of the corresponding command file, which is written into the .log file. The purpose
! is to facilitate the creation of a subcluster that is especially suited for direct
! calibration. Those results can be turned around and used as calibration events for indirect
! calibration of the larger cluster.

   implicit none

   include 'mloc.inc'

   logical :: command_help
   character(len=132) :: msg

   if (command_help) then
      write (*,'(7(a/))')&
       '"subc" (subcluster) specifies the parameters to use to select events that would be',&
       'most suitable for a high-quality subcluster with direct calibration. Three parameters',&
       'are required:',&
       '   1 - the maximum epicentral distance',&
       '   2 - the minimum number of readings within that distance',&
       '   3 - the minimum number of station-phases in common with other events in the cluster',&
       'The main body of the corresponding command file is written to the .log file.'
      return
   end if

   call decode (params, value, nvar)
   
   if (nvar .eq. 3) then
      subc_set = .true.
      subc_delt = value(1)
      subc_nmin = nint(value(2))
      subc_nconnect = nint(value(3))
      if (verbose_screen) then
         write (msg,'(a,f6.2,2(1x,i3))') 'proc_subc: ', subc_delt, subc_nmin, subc_nconnect
         call fyi (trim(msg))
      end if
   else
      call warnings ('proc_subc: "subc" requires three arguments')
   end if 

   return
   
end subroutine proc_subc


!***********************************************************************************************************************************
subroutine proc_taup (command_help)

! Global travel-time model, using Tau-P formulation

   implicit none

   include 'mloc.inc'

   logical :: command_help

   if (command_help) then
      write (*,'(3(a/))')&
       '"taup" (Tau-P) specifies a global travel time model using the tau-p formulation.',&
       'Currently, ak135 is the only option. The command takes a single argument, the model',&
       'name. The associated ".tbl" and ".hed" files must be stored in the "tables/tau-p" subdirectory.'
      return
   end if
   
   if (params(1:1) .ne. ' ') then
      taup_model = trim(params)
      if (verbose_screen) call fyi ('proc_taup: using '//trim(taup_model))
   else
      call warnings ('proc_taup: "taup" requires one argument, the name of the velocity model')
   end if

   return
   
end subroutine proc_taup


!***********************************************************************************************************************************
subroutine proc_tfil (command_help)

! Travel-time spread file.

   implicit none

   include 'mloc.inc'

   logical :: ex, command_help

   if (command_help) then
      write (*,'(4(a/))')&
      '"tfil" specifies a .ttsprd file to be opened for reading the empirically-derived',&
      'spreads of different phases. The file is expected to be in the data directory, so',&
      'the argument of the command is just the filename, not the full pathname. Revert to',&
      'the default values by entering the argument "default".'
   end if

   ttsprdfname = trim(datadir)//dirsym//trim(params)
   inquire (file=ttsprdfname,exist=ex)
   if (ex) then
      read_ttsprd = .true.
      if (verbose_screen) call fyi ('proc_tfil: using '//trim(ttsprdfname))
   else
      if (trim(params) .eq. 'default') then
         ttsprdfname = 'default'
         read_ttsprd = .false.
      else
         call warnings ('proc_tfil: file '//trim(ttsprdfname)//' does not exist')
         read_ttsprd = .false.
      end if
   end if

   return
   
end subroutine proc_tfil


!***********************************************************************************************************************************
subroutine proc_tikh (command_help)

! Tikhonov regularization of the cluster vector perturbations

   implicit none

   include 'mloc.inc'

   logical :: command_help
   character(len=132) :: msg

   if (command_help) then
      write (*,'(5(a/))')&
      '"tikh" specifies the value to be used for Tikhonov regularization of the',&
      'perturbations for cluster vectors. If a value of 0 is used (default) there is no',&
      'regularization. For data sets exhibiting convergence problems it may be helpful to',&
      'set a non-zero (positive) value. Determination of the optimal value is a non-trivial',&
      'task, but values in the range 0.2-0.6 seem to be about right for this application.'
      return
   end if

   call decode (params, value, nvar)
   
   if (nvar .eq. 1) then
      tikhonov_factor = value(1)
      if (tikhonov_factor .lt. 0.) then
         call warnings ('proc_tikh: the Tikhonov regularization factor must be greater than zero')
         tikhonov_factor = 0.
      end if
      if (tikhonov_factor .gt. 1.0) then
         call fyi ('proc_tikh: a value greater than 1.0 may cause excessive damping')
      end if
      if (verbose_screen) then
         write (msg,'(a,f4.2)') 'proc_tikh: Tikhonov regularization factor = ', tikhonov_factor
         call fyi (trim(msg))
      end if
   else
      call warnings ('proc_tikh: "tikh" requires one argument')
   end if

   return
   
end subroutine proc_tikh


!***********************************************************************************************************************************
subroutine proc_time (iev, command_help)

! Starting origin time

   implicit none

   include 'mloc.inc'

   integer :: iev
   real :: hms2s
   logical :: command_help
   character(len=132) :: msg

   if (command_help) then
      write (*,'(5(a/))')&
       '"time" specifies a starting origin time. It requires three arguments (hours,',&
       'minutes, seconds). Unlike other location parameters there is never a need to set',&
       'all starting origin times to a common value, so it can only be aplied to the current',&
       'event. If the "time" command is not applied to an event, the time read from the input',&
       'data file or an HDF file will be used as the starting value.'
      return
   end if

   call decode (params, value, nvar)
   
   if (nvar .eq. 3) then
      if (iev .eq. 0) then ! This is never needed and not allowed
         call warnings ('proc_time: "time" can only be issued after an event is declared')
      else ! Application to individual events
         if (value(1) .ge. 0.0 .and. value(1) .le. 23.0 .and.&
             value(2) .ge. 0.0 .and. value(2) .le. 59.0 .and.&
             value(3) .ge. 0.0 .and. value(3) .le. 59.99) then
            time_cf(iev) = hms2s(nint(value(1)),nint(value(2)),value(3))
            if (verbose_screen) then
               write (msg,'(a,i3,3f6.2,f8.2)') 'proc_time: event ', iev, value(1), value(2), value(3), time_cf(iev)
               call fyi (trim(msg))
            end if
         else
            call warnings ('proc_time: OT is outside bounds')
         end if
      end if
   else
      call warnings ('proc_time: "time" requires three parameters')
   end if

   return
   
end subroutine proc_time


!***********************************************************************************************************************************
subroutine proc_tomo (command_help)

! Output files for tomography.

   implicit none

   include 'mloc.inc'

   integer :: lp, itomo_in, lenb
   logical :: command_help
   character(len=8) :: tomo_phase_in
   character(len=132) :: msg

   if (command_help) then
      write (*,'(7(a/))')&
       '"tomo" (tomography) specifies output files for tomography. This command requires two',&
       'arguments, a phase name and a flag for which kind of data to extract:',&
       '   1 = Extract all readings of the specified phase',&
       '   2 = Extract only readings which were used for the cluster vectors',&
       '   3 = Extract empirical path anomalies',&
       'The "tomo" command can be issued multiple times, up to the limit specified by the',&
       'variable "nitomomax".'
      return
   end if

   lp = lenb(params)
   read (params(lp:lp),'(i1)') itomo_in
   if (itomo_in .ge. 1 .and. itomo_in .le. 3) then
      tomo_phase_in = params(1:lp-2)
   else
      call warnings ('proc_tomo: error in input line')
   end if
   nitomo = nitomo + 1
   if (nitomo .le. nitomomax) then
      tomo_phase(nitomo) = tomo_phase_in
      itomo(nitomo) = itomo_in
      if (verbose_screen) then
         write (msg,'(t3,2a,1x,i1)') 'proc_tomo: tomography output: ', trim(tomo_phase(nitomo)), itomo(nitomo)
         call fyi (trim(msg))
      end if
   else
      call warnings ('proc_tomo: maximum number of instances reached')
   end if

   return
   
end subroutine proc_tomo


!***********************************************************************************************************************************
subroutine proc_tptt (command_help)

! T-phase travel time calculations.

   implicit none

   include 'mloc.inc'

   logical :: command_help
   character(len=132) :: msg

   if (command_help) then
      write (*,'(2(a/))')&
       '"tptt" specifies an intercept (sec) and slope (sec/degree) to be used in calculating',&
       'travel tiems for T-phases.'
      return
   end if

   call decode (params, value, nvar)
   
   if (nvar .eq. 2) then
      tphase_a = value(1)
      tphase_b = value(2)
      if (verbose_screen) then
         write (msg,'(a,2f8.2)') 'proc_tptt:', tphase_a, tphase_b
         call fyi (trim(msg))
      end if
   else
      call warnings ('proc_tptt: "tptt" requires two arguments')
   end if

   return
   
end subroutine proc_tptt


!*****************************************************************************************
subroutine proc_tt5e (command_help)

! Special case of the 'pltt' command that allows a type 5 TT plot (local distance, out ot 4°)
! to be made for a single event. This is useful for investigating problem events. The command
! can be issued multiple times. The plots are put into a subdirectory.

   implicit none

   include 'mloc.inc'

   logical :: command_help
   character(len=132) :: msg

   if (command_help) then
      write (*,'(3(a/))')&
       '"tt5e" allows a type 5 plot (local distance, out to 4°) to be made for a single event.',&
       'The command takes a single argument, the event name (as given in the "even" command.',&
       'The command can be issued multiple times, up to the limit specified by "n_tt5e_max".'
      return
   end if

   if (params(1:1) .ne. ' ') then
      if (len_trim(params) .le. 30) then
         if (verbose_screen) call fyi ('proc_tt5e: single station tt5 plot for event '//trim(params))
         tt5e = .true.
         n_tt5e = n_tt5e + 1
         if (n_tt5e .le. n_tt5e_max) then
            tt5e_evt(n_tt5e) = trim(params)
         else
            call warnings ('proc_tt5e: maximum number of instances reached')
            n_tt5e = n_tt5e_max
         end if
      else
         msg = 'proc_tt5e: event name must be 30 characters or less ('//trim(params)//')'
         call warnings (trim(msg))
      end if
   else
      call warnings ('proc_tt5e: "tt5e" requires one argument, an event name')
   end if

   return
   
end subroutine proc_tt5e


!*****************************************************************************************
subroutine proc_tt5s (command_help)

! Special case of the 'pltt' command that allows a type 5 TT plot (local distance, out ot 4°)
! to be made for a single station. This is useful for investigating the phase IDs near the
! cross-over distance for Pg-Pn and Sg-Sn. The command can be issued multiple times. The plots
! are put into a subdirectory, as for the relative depth phase plots (type 8).

   implicit none

   include 'mloc.inc'

   logical :: command_help
   character(len=132) :: msg

   if (command_help) then
      write (*,'(3(a/))')&
       '"tt5s" allows a type 5 plot (local distance, out to 4°) to be made for a single station.',&
       'The command takes a single argument, the station code (case-sensitive). The command',&
       'can be issued multiple times, up to the limit specified by "n_tt5s_max".'
      return
   end if

   if (params(1:1) .ne. ' ') then
      if (len_trim(params) .le. 6) then
         if (verbose_screen) call fyi ('proc_tt5s: single station tt5 plot for station '//params)
         tt5s = .true.
         n_tt5s = n_tt5s + 1
         if (n_tt5s .le. n_tt5s_max) then
            tt5s_sta(n_tt5s) = trim(params)
            if (ichar(tt5s_sta(n_tt5s)(1:1)) .ge. 97)  call warnings ('proc_tt5s: station codes are case-sensitive')
         else
            call warnings ('proc_tt5s: maximum number of instances reached')
            n_tt5s = n_tt5s_max
         end if
      else
         msg = 'proc_tt5s: station code must be 6 characters or less ('//trim(params)//')'
         call warnings (trim(msg))
      end if
   else
      call warnings ('proc_tt5s: "tt5s" requires one argument, a station code')
   end if

   return
   
end subroutine proc_tt5s


!***********************************************************************************************************************************
subroutine proc_ttou (command_help)

! Empirical travel times output file

   implicit none

   include 'mloc.inc'

   character(len=132) :: msg
   logical :: command_help

   if (command_help) then
      write (*,'(3(a/))')&
       '"ttou" (TT output) specifies that an output file of empirical travel time data will',&
       ' be written for a specific phase. The command may be issued multiple times, up to',&
       ' the limit in n_ttou_max. S-P is supported but the format of the output file is different.'
      return
   end if

   if (params(1:1) .ne. ' ') then
      if (len_trim(params) .le. 8) then
         if (verbose_screen) call fyi ('proc_ttou: output file of empirical TTs for phase '//params)
         ttou = .true.
         n_ttou = n_ttou + 1
         if (n_ttou .le. n_ttou_max) then
            ttou_phase(n_ttou) = trim(params)
         else
            call warnings ('proc_ttou: maximum number of instances reached')
            n_ttou = n_ttou_max
         end if
      else
         msg = 'proc_ttou: phase code must be 8 characters or less ('//trim(params)//')'
         call warnings (trim(msg))
      end if
   else
      call warnings ('proc_ttou: "ttou" requires one argument, a phase code')
   end if

   return
   
end subroutine proc_ttou


!***********************************************************************************************************************************
subroutine proc_vect (command_help)

! Plotting of relocation vectors.

   implicit none

   include 'mloc.inc'

   integer :: i
   logical :: command_help
   character(len=132) :: msg

   if (command_help) then
      write (*,'(6(a/))')&
       '"vect" (vectors) controls the plotting of four types of relocation vectors:',&
       '   1) From data file epicenter to final location, either calibrated or uncalibrated, in black.',&
       '   2) From starting location to final location (but not with calibration shift), in green.',&
       '   3) Calibration shift, for indirect calibration, in red.',&
       '   4) Residual calibration shift, for indirect calibration, in blue.',&
       'The command takes four arguments which must be either 0 or 1, to set the plotting of these',&
       'vectors. The default state is TRUE for all four vectors.'
      return
   end if

   call decode (params, value, nvar)
   
   if (nvar .eq. 4) then
      do i = 1,4
         if (nint(value(i)) .eq. 0) then
            vectors(i) = .false.
         else if (nint(value(i)) .eq. 1) then
            vectors(i) = .true.
         else
            call warnings ('proc_vect: illegal argument')
         end if
      end do
   else
      call warnings ('proc_vect: "vect" requires four parameters, integers (0 or 1) specifying relocation vector plotting')
   end if

   if (verbose_screen) then
      write (msg,'(a,4(l1,1x))') 'proc_vect: vector plotting: ', vectors(1), vectors(2), vectors(3), vectors(4)
      call fyi (trim(msg))
   end if

   return
   
end subroutine proc_vect


!***********************************************************************************************************************************
subroutine proc_vlog (command_help)

! Verbose_log mode, extra output to log file

   implicit none

   include 'mloc.inc'

   logical :: command_help

   if (command_help) then
      write (*,'(4(a/))')&
       '"vlog" (verbose log) specifies that extra information about the HD analysis will',&
       'be written to the log file. Issuing the command without arguments toggles the',&
       'current state. The arguments "on" and "off" may be used to set the state explicitly.',&
       'Default is "on".'
      return
   end if

   if (params(1:1) .eq. ' ') then
      verbose_log = .not.verbose_log
   else if (params(1:2) .eq. 'on') then
      verbose_log = .true.
   else if (params(1:3) .eq. 'off') then
      verbose_log = .false.
   else
      call warnings ('proc_vlog: illegal argument')
   end if

   if (verbose_log) then
      if (verbose_screen) call fyi ('proc_vlog: using verbose mode for log')
   else
      if (verbose_screen) call fyi ('proc_vlog: using quiet mode for log')
   end if

   return
   
end subroutine proc_vlog


!***********************************************************************************************************************************
subroutine proc_vscr (command_help)

! Verbose_screen mode, extra output to standard output

   implicit none

   include 'mloc.inc'
   include 'ttlim.inc'

   logical :: command_help

   if (command_help) then
      write (*,'(4(a/))')&
       '"vscr" (verbose screen) specifies that extra information about the HD analysis',&
       'will be written to the output window. Issuing the command without arguments',&
       'toggles the current state. The arguments "on" and "off" may be used to set the',&
       'state explicitly. Default is "off".'
      return
   end if

   if (params(1:1) .eq. ' ') then
      verbose_screen = .not.verbose_screen
   else if (params(1:2) .eq. 'on') then
      verbose_screen = .true.
   else if (params(1:3) .eq. 'off') then
      verbose_screen = .false.
   else
      call warnings ('proc_vscr: illegal argument')
   end if
   verbose_taup = verbose_screen

   if (verbose_screen) call fyi ('proc_vscr: using verbose mode for screen')

   return
   
end subroutine proc_vscr


!***********************************************************************************************************************************
subroutine proc_weig (command_help)

! Weighting of residuals by inverse of errors

   implicit none

   include 'mloc.inc'

   logical :: command_help

   if (command_help) then
      write (*,'(3(a/))')&
       '"weig" (weight) determines if the data (residuals) are weighted equally ("off"),',&
       'or if the data will be weighted inversely to their reading error ("on"). Issuing',&
       'the command with no argument toggles the current state. Default is "on".'
      return
   end if

   if (params(1:1) .eq. ' ') then
      data_weight = .not.data_weight
   else if (params(1:2) .eq. 'on') then
      data_weight = .true.
   else if (params(1:3) .eq. 'off') then
      data_weight = .false.
   else
      call warnings ('proc_weig: illegal argument')
   end if

   if (data_weight) then
      data_weight_pr = 'on '
   else
      data_weight_pr = 'off'
   end if
   if (verbose_screen) call fyi ('proc_weig: weighting of residuals is '//data_weight_pr)

   return
   
end subroutine proc_weig


!***********************************************************************************************************************************
subroutine proc_wind (command_help)

! Windowing of residuals

   implicit none

   include 'mloc.inc'

   logical :: command_help
   character(len=132) :: msg

   if (command_help) then
      write (*,'(14(a/))')&
       '"wind" (windowing) specifies two parameters (and an optional third one) used to define',&
       'windows for each phase that are used to cut large residuals out of the problem. The',&
       'parameters given here are multipliers for the travel time spread assigned to each',&
       'phase. The travel time spreads can be taken from default values or read from a .ttsprd',&
       'file from a previous run. Default values are 3 and 4. With these values, all residuals',&
       'within 3*tsprd(phase) will be given full weight. Residuals between 3*tsprd(phase) and',&
       '4*tsprd(phase) have weights that taper smoothly to zero with a 1-cosine function. The',&
       '.ttsprd file also carries a baseline adjustment for each phase that helps keep the',&
       'window function centered over the actual distribution of that phase, rather than',&
       'centering it on the ak135 theoretical time for that phase. The optional third parameter',&
       'is an epicentral distance below which the windows are expanded by a factor of 2. This',&
       'is used in direct calibration to help keep good readings from being lost because of a',&
       'poor starting location. The default value is 1.2 degrees. Weights assigned via the',&
       '"wind" command are written to the .phase_data output file.'
      return
   end if
   
   call decode (params, value, nvar)
   
   if (nvar .eq. 2 .or. nvar .eq. 3) then
      wind1 = value(1)
      wind2 = value(2)
      if (wind2 .lt. wind1) then
         call warnings ('proc_wind: wind2 must be ≥ wind1')
         wind2 = wind1
      end if
      if (nvar .eq. 3) windloclim = value(3)
      if (verbose_screen) then
         if (nvar .eq. 2) write (msg,'(a,2f4.1)') 'proc_wind: ', wind1, wind2
         if (nvar .eq. 3) write (msg,'(a,3f4.1)') 'proc_wind: ', wind1, wind2, windloclim
         call fyi (trim(msg))
      end if
   else
      call warnings ('proc_wind: "wind" requires two or three parameters')
   end if 

   return
   
end subroutine proc_wind


!***********************************************************************************************************************************
subroutine proc_xsec (command_help)

! Plotting of cross-sections (depth profiles)

   implicit none

   include 'mloc.inc'

   logical :: command_help
   character(len=132) :: msg

   if (command_help) then
      write (*,'(8(a/))')&
       '"xsec" controls the plotting of cross-sections. The command takes six arguments:',&
       '   Latitude of the first end-point ',&
       '   Longitude of the first end-point ',&
       '   Latitude of the second end-point ',&
       '   Longitude of the second end-point ',&
       '   Depth (km) of the cross-section (all sections start at zero depth) ',&
       '   Full width (km) of the cross-section ',&
       'Multiple cross-sections can be defined, up the limit set by the parameter "n_xsec_max"'
      return
   end if

   call decode (params, value, nvar)
   
   if (nvar .eq. 6) then
      if (n_xsec .lt. n_xsec_max) then
         n_xsec = n_xsec + 1
         xsec_lat1(n_xsec) = value(1)
         xsec_lon1(n_xsec) = value(2)
         xsec_lat2(n_xsec) = value(3)
         xsec_lon2(n_xsec) = value(4)
         xsec_depth(n_xsec) = value(5)
         xsec_width(n_xsec) = value(6)
         if (verbose_screen) then
            write (msg,'(a,6f8.3)') 'proc_xsec: ', xsec_lat1(n_xsec), xsec_lon1(n_xsec), xsec_lat2(n_xsec), xsec_lon2(n_xsec),&
             xsec_depth(n_xsec), xsec_width(n_xsec)
            call fyi (trim(msg))
         end if
      else
         call warnings ('proc_xsec: maximum number of cross sections reached')
      end if
   else
      call warnings ('proc_xsec: "xsec" requires six arguments')
   end if

   return
   
end subroutine proc_xsec


!***********************************************************************
subroutine parse (linein, comd, params)

! This subroutine breaks the string 'line' into the command 'comd' and a 
! string (params) containing all of the parameter entries.
! 7/8/2018: If present, the last instance of character "!" and any following text are removed before processing

   implicit none

!  include 'mloc.inc'

   integer :: i, j, k, nf, lenb, nn, cmd_comment
   character(len=4) :: comd, buf
   character(len=80) :: linein, line
   character(len=76) :: params

   params(1:len(params)) = ' '
   comd(1:4) = '    '
   line = ' '
   
   cmd_comment = index(linein,'!',.true.)
   if (cmd_comment .eq. 0) then
      line = linein
   else if (cmd_comment .gt. 1) then
      k = cmd_comment - 1
      line = linein(1:k)
   end if

   ! Find the command, skipping leading blanks
   do i = 1,lenb(line)
      if (line(i:i) .ne. ' ') then
         nn = i
         buf = line(i:i+3)
         comd = buf(1:lenb(buf))
         exit
      end if
   end do

   ! Find blank between command and parameters
   do j = nn,len(line)
      if (line(j:j) .eq. ' ') then
         nf = j
         exit
      end if
   end do
   
   ! Find parameter string
   do j = nf,len(line)
      if (line(j:j) .ne. ' ') then
         k = lenb(line)
         params(1:k-j+1) = line(j:k)
         exit
      end if
   end do 

   return
   
end subroutine parse


!***********************************************************************
subroutine decode (string, xvar, nvar)

!  subroutine to decode character string 'string' into 'nvar' real
!  variables in 'xvar' with blanks or commas as the delimiting 
!  characters.

!  no subroutines called

   implicit none
      
!  include 'mloc.inc'

   character(len=*) :: string
   character(len=22) :: value
   character(len=8) :: fmt
   real :: xvar(12)
   integer :: i, j, k, l, npt, nin, nvar, nd, nf, nl, iv, nw, lenb, nlen, ios
   
   fmt(1:1) = '('
   fmt(5:5) = '.'
   fmt(8:8) = ')'

   !  find first nonblank character

   nlen=len(string)
   do i=1,nlen
      if (string(i:i) .ne. ' ') then
         nf=i
         go to 10
      end if
   end do

   !  null string

   nvar=0 
   return

   !  find last nonblank character

10 continue
   nl=lenb(string)
   l=1
   k=nf

   !  decode each character in turn

20 continue

   !  find a delimiter character

   do j=k,nl+1
      if ((string(j:j) .eq. ' ') .or. (string(j:j) .eq. ',')) then
         value(1:len(value))=string(k:j-1)
         nw=j-k
         if (nw .eq. 0) then
            k=j+1
            go to 20
         end if

         !  search for e-use e format if found

         nin=index(value,'e')
         if (nin .ne. 0) then
            fmt(2:2)='e'
            write (fmt(3:4),'(i2)',iostat=ios) nw
            if (ios .gt. 0) go to 40
            npt=index(value,'.')
            if (npt .ne. 0) then
               nd=nin-npt-1
            else
               nd=0
            end if
            write (fmt(6:7),'(i2)',iostat=ios) nd
            if (ios .gt. 0) go to 40
            read (value,fmt,iostat=ios) xvar(l)
            if (ios .ne. 0) go to 40
            l=l+1
            k=j+1
            go to 20
         end if

         !  search for decimal point-use f format if found

         nin=index(value,'.')
         fmt(2:2)='f'
         if (nin .ne. 0) then
            nd=lenb(value)-nin
            write (fmt(3:4),'(i2)',iostat=ios) nw
            if (ios .gt. 0) go to 40
            write (fmt(6:7),'(i2)',iostat=ios) nd
            if (ios .gt. 0) go to 40
            read (value,fmt,iostat=ios) xvar(l)
            if (ios .ne. 0) go to 40
            l=l+1
            k=j+1
            go to 20
         else

            !  no decimal point-put one in

            iv=lenb(value)+1
            value(iv:iv)='.'
            nw=nw+1
            nd=0
            write (fmt(3:4),'(i2)',iostat=ios) nw
            if (ios .gt. 0) go to 40
            write (fmt(6:7),'(i2)',iostat=ios) nd
            if (ios .gt. 0) go to 40
            read (value,fmt,iostat=ios) xvar(l)
            if (ios .ne. 0) go to 40
            l=l+1
            k=j+1
            go to 20
         end if
      end if
   end do

   !  no delimiter found-done

   nvar=l-1 
   return

   !  error return

40 continue
   nvar=-1
   
   return

end subroutine decode


!*****************************************************************************************
subroutine decode2 (string_in, nvar, param1, param2, param3, param4)

! Breaks the input "string" into nvar substrings.
! Substrings are assumed to be separated by blanks and have a maximum length of 8 characters.
! A maximum of four substrings is supported.

   implicit none
   
   character(len=*) :: string_in
   character(len=80) :: string
   character(len=8) :: param1, param2, param3, param4
   integer :: nvar, nblank, index_blank(4), i, n
   
   param1 = ' '
   param2 = ' '
   param3 = ' '
   param4 = ' '
   index_blank = 0
   nvar = 0
   
   string = adjustl(string_in)
   n = len_trim(string)
   if (n .eq. 0) then
      call warnings ('decode2: empty parameter string')
      return
   end if

   ! Number of substrings and indices of the blanks
   nblank = 0
   do i = 1,n
      if (string(i:i) .eq. ' ') then
         nblank = nblank + 1
         index_blank(nblank) = i
         if (nblank .gt. 3) then
            call warnings ('decode2: too many separators')
            nvar = 0
            return
         end if
      end if
   end do
   nvar = nblank + 1
   index_blank(nvar) = n + 1
   
   if (nvar .ge. 1) then
      param1 = string(1:index_blank(1)-1)
   end if
   
   if (nvar .ge. 2) then
      param2 = string(index_blank(1)+1:index_blank(2)-1)
   end if
   
   if (nvar .ge. 3) then
      param3 = string(index_blank(2)+1:index_blank(3)-1)
   end if
   
   if (nvar .eq. 4) then
      param4 = string(index_blank(3)+1:index_blank(4)-1)
   end if
         
   return

end subroutine decode2
