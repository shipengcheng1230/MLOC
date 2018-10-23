      subroutine mlocout_kml (it)
      
      ! Creates a kml file of the epicenters that can be displayed in Google Earth
      ! Expects to find icons in a folder named "tables/kml" under the executable
      ! mloc file. Icons are color-coded by depth for events with depth constraint:
      !   0 -  9 km - red
      !  10 - 19 km - green
      !  20 - 29 km - skyblue
      !  30+     km - blue
      ! A yellow icon is used for events set to the cluster default depth
      
      implicit none
      
      include 'mloc.inc'
      
      integer it, it1, iev, irange
      character*100 outfil
      character*30 dep_con
      character*24 cdate
      character*20 color
      real xlat, xlon, xdep, xmagms, xmagmw, xlonmin, xlonmax, ylatmin, ylatmax
      real xlath, xlonh
      real delta, deltdg, deltkm, azeqst, azesdg, azsteq, azsedg
      character(len=132) :: msg
      
      it1 = it + 1
      color = ' '
            
      outfil = trim(outfile)//'.kml'
      if (verbose_screen) then
         write (msg,'(3a,i3)') 'mlocout_kml: opening ', trim(outfil), ' on unit ', io_out
         call fyi (trim(msg))
      end if
      open (io_out,file=outfil,status='new')
      
      ! Focus point for KML is the iterated hypocentroid. Don't worry about correction in the case of indirect calibration
      xlath = lath(it1)
      xlonh = lonh(it1)
      call set_longitude_range (xlonh, 0) ! Google Earth requires longitude -180 < xlon < 180
      
      ! Calculate range for KML "LookAt" command from the diagonal length of the cluster bounding box
      call map_boundaries (it1, 0, xlonmin, xlonmax, ylatmin, ylatmax)
      call delaz (ylatmin, xlonmin, ylatmax, xlonmax, delta, deltdg, deltkm, azeqst, azesdg, azsteq, azsedg, 0)
      irange = nint(deltkm*1.e3)
            
      call kml_prelude (io_out, xlath, xlonh, irange)
      
      do iev = 1,nev ! Loop over events

         write (cdate,'(i4,a,i2,a,i2,1x,i2,a,i2,a,f6.3)') iyre(iev), '/', mone(iev), '/', idye(iev),&
          hourp(iev,it1), ':', minp(iev,it1), ':', secp(iev,it1)

         if (calibration) then
            xlat = latp_cal(iev)
            xlon = lonp_cal(iev)
            xdep = depthp_cal(iev)
         else
            xlat = latp(iev,it1)
            xlon = lonp(iev,it1)
            xdep = depthp(iev,it1)
         end if
         
         call set_longitude_range (xlon, 0) ! Google Earth requires longitude -180 < xlon < 180
         
         if (depset_pr(iev) .eq. 'c') then
            color = 'yellow'
            dep_con = ' cluster default'
         else
            dep_con = ' constrained ('//depset_pr(iev)//')'
            if (xdep .lt. 10.0) then
               color = 'red'
            else if (xdep .lt. 20.0) then
               color = 'green'
            else if (xdep .lt. 30.0) then
               color = 'skyblue'
            else
               color = 'blue'
            end if
         end if
    
         xmagms = 0.
         xmagmw = 0.

         call kml_point (io_out, iev, cdate, xlat, xlon, xdep, rmag(iev), xmagms, xmagmw, color, dep_con)
      
      end do
  
      call kml_coda (io_out)
      
      close (io_out)
  
      return
      end
    

!*****************************************************************************************
      subroutine kml_point (iout, iev, cdate, xlat, xlon, xdep, xmagmb, xmagms, xmagmw, color, dep_con)
       
      implicit none
      
      integer iout, iev
      character*30 dep_con
      character*24 cdate
      character*21 icon
      character*20 color
      character*8 cxlat, cxlon
      character*5 cxdep
      character*6 cscale
      character*3 cxiev, cxmagmb, cxmagms, cxmagmw
      real xlon, xlat, xdep, xmag, xmagmb, xmagms, xmagmw, scale, scale0, scale_factor
      
      scale0 = 0.3
      scale_factor = 1.75

      icon = '#'//color
      
      xmag = amax1(xmagmb, xmagms, xmagmw)
      if (xmag .gt. 2.) then
         scale = scale0 + ((xmag - 2.0)/5.0)*scale_factor
      else
         scale = scale0
      end if
      write (cscale,'(f6.3)') scale
      
      write (cxiev,'(i3)') iev
      write (cxlat,'(f8.3)') xlat
      write (cxlon,'(f8.3)') xlon
      write (cxdep,'(f5.1)') xdep
      
      write (iout,'(a)')  '   <Placemark>'
      write (iout,'(3a)') '      <name>', trim(adjustl(cxiev)), '</name>'
      write (iout,'(3a)') '      <description><![CDATA[Event : <b>', trim(adjustl(cxiev)), '</b>'
      write (iout,'(3a)')   '      <br>Time : <b>', cdate, '</b>'
      if (xmagmb .lt. 0.1 .and. xmagms .lt. 0.1 .and. xmagmw .lt. 0.1) then
         write (iout,'(a)') '      <br>Magnitude : <b>Unknown</b>'
      else
         if (xmagmb .gt. 0.) then
            write (cxmagmb,'(f3.1)') xmagmb
            write (iout,'(4a)') '      <br>Magnitude : <b>', trim(cxmagmb), ' mb', '</b>'
         end if
         if (xmagms .gt. 0.) then
            write (cxmagms,'(f3.1)') xmagms
            write (iout,'(4a)') '      <br>Magnitude : <b>', trim(cxmagmb), ' Ms', '</b>'
         end if
         if (xmagmw .gt. 0.) then
            write (cxmagmw,'(f3.1)') xmagmw
            write (iout,'(4a)') '      <br>Magnitude : <b>', trim(cxmagmb), ' Mw', '</b>'
         end if
      end if
      write (iout,'(3a)') '      <br>Latitude : <b>', trim(adjustl(cxlat)), '</b>'     
      write (iout,'(3a)') '      <br>Longitude : <b>', trim(adjustl(cxlon)), '</b>'     
      write (iout,'(4a)') '      <br>Depth : <b>', trim(adjustl(cxdep)), trim(dep_con), '</b>]]></description>'     
      write (iout,'(3a)') '      <styleUrl>', trim(icon), '</styleUrl>'
      write (iout,'(a)')  '      <Style>'
      write (iout,'(a)')  '         <IconStyle>'
      write (iout,'(3a)') '            <scale>', trim(adjustl(cscale)), '</scale>'
      write (iout,'(a)')  '            <hotSpot x="0.5" y="0.5" xunits="fraction" yunits="fraction" />'
      write (iout,'(a)')  '         </IconStyle>'
      write (iout,'(a)')  '      </Style>'
      write (iout,'(a)')  '      <Point>'
      write (iout,'(5a)') '         <coordinates>', trim(adjustl(cxlon)), ',', trim(adjustl(cxlat)), ',0</coordinates>'
      write (iout,'(a)')  '      </Point>'
      write (iout,'(a)')  '   </Placemark>'
  
      return
      end
    

!*****************************************************************************************
      subroutine kml_prelude (iout, xlath, xlonh, irange)
      
      implicit none
      
      integer iout, irange
      real xlath, xlonh
      character*8 cxlath, cxlonh, cirange
      
      write (cxlath,'(f8.3)') xlath
      write (cxlonh,'(f8.3)') xlonh
      write (cirange,'(i8)') irange
      
      write (iout,'(a)') '<?xml version="1.0" encoding="UTF-8"?>'
      write (iout,'(a)') '<kml xmlns="http://www.opengis.net/kml/2.2">'
      write (iout,'(a)') '<Document>'
      
      write (iout,'(a)') '   <Style id="i_red">'
      write (iout,'(a)') '      <IconStyle>'
      write (iout,'(a)') '         <Icon>'
      write (iout,'(a)') '            <href>../tables/kml/i_red.png</href>'
      write (iout,'(a)') '         </Icon>'
      write (iout,'(a)') '      </IconStyle>'
      write (iout,'(a)') '   </Style>'
      write (iout,'(a)') '   <Style id="a_red">'
      write (iout,'(a)') '      <IconStyle>'
      write (iout,'(a)') '         <Icon>'
      write (iout,'(a)') '            <href>../tables/kml/a_red.png</href>'
      write (iout,'(a)') '         </Icon>'
      write (iout,'(a)') '      </IconStyle>'
      write (iout,'(a)') '   </Style>'
      write (iout,'(a)') '   <StyleMap id="red">'
      write (iout,'(a)') '      <Pair>'
      write (iout,'(a)') '         <key>normal</key>'
      write (iout,'(a)') '         <styleUrl>#i_red</styleUrl>'
      write (iout,'(a)') '      </Pair>'
      write (iout,'(a)') '      <Pair>'
      write (iout,'(a)') '         <key>highlight</key>'
      write (iout,'(a)') '         <styleUrl>#a_red</styleUrl>'
      write (iout,'(a)') '      </Pair>'
      write (iout,'(a)') '  </StyleMap>'
      
      write (iout,'(a)') '   <Style id="i_yellow">'
      write (iout,'(a)') '      <IconStyle>'
      write (iout,'(a)') '         <Icon>'
      write (iout,'(a)') '            <href>../tables/kml/i_yellow.png</href>'
      write (iout,'(a)') '         </Icon>'
      write (iout,'(a)') '      </IconStyle>'
      write (iout,'(a)') '   </Style>'
      write (iout,'(a)') '   <Style id="a_yellow">'
      write (iout,'(a)') '      <IconStyle>'
      write (iout,'(a)') '         <Icon>'
      write (iout,'(a)') '            <href>../tables/kml/a_yellow.png</href>'
      write (iout,'(a)') '         </Icon>'
      write (iout,'(a)') '      </IconStyle>'
      write (iout,'(a)') '   </Style>'
      write (iout,'(a)') '   <StyleMap id="yellow">'
      write (iout,'(a)') '      <Pair>'
      write (iout,'(a)') '         <key>normal</key>'
      write (iout,'(a)') '         <styleUrl>#i_yellow</styleUrl>'
      write (iout,'(a)') '      </Pair>'
      write (iout,'(a)') '      <Pair>'
      write (iout,'(a)') '         <key>highlight</key>'
      write (iout,'(a)') '         <styleUrl>#a_yellow</styleUrl>'
      write (iout,'(a)') '      </Pair>'
      write (iout,'(a)') '   </StyleMap>'
      
      write (iout,'(a)') '   <Style id="i_green">'
      write (iout,'(a)') '      <IconStyle>'
      write (iout,'(a)') '         <Icon>'
      write (iout,'(a)') '            <href>../tables/kml/i_green.png</href>'
      write (iout,'(a)') '         </Icon>'
      write (iout,'(a)') '      </IconStyle>'
      write (iout,'(a)') '   </Style>'
      write (iout,'(a)') '   <Style id="a_green">'
      write (iout,'(a)') '      <IconStyle>'
      write (iout,'(a)') '         <Icon>'
      write (iout,'(a)') '            <href>../tables/kml/a_green.png</href>'
      write (iout,'(a)') '         </Icon>'
      write (iout,'(a)') '      </IconStyle>'
      write (iout,'(a)') '   </Style>'
      write (iout,'(a)') '   <StyleMap id="green">'
      write (iout,'(a)') '      <Pair>'
      write (iout,'(a)') '         <key>normal</key>'
      write (iout,'(a)') '         <styleUrl>#i_green</styleUrl>'
      write (iout,'(a)') '      </Pair>'
      write (iout,'(a)') '      <Pair>'
      write (iout,'(a)') '         <key>highlight</key>'
      write (iout,'(a)') '         <styleUrl>#a_green</styleUrl>'
      write (iout,'(a)') '      </Pair>'
      write (iout,'(a)') '   </StyleMap>'
      
      write (iout,'(a)') '   <Style id="i_blue">'
      write (iout,'(a)') '      <IconStyle>'
      write (iout,'(a)') '         <Icon>'
      write (iout,'(a)') '            <href>../tables/kml/i_blue.png</href>'
      write (iout,'(a)') '         </Icon>'
      write (iout,'(a)') '      </IconStyle>'
      write (iout,'(a)') '   </Style>'
      write (iout,'(a)') '   <Style id="a_blue">'
      write (iout,'(a)') '      <IconStyle>'
      write (iout,'(a)') '         <Icon>'
      write (iout,'(a)') '            <href>../tables/kml/a_blue.png</href>'
      write (iout,'(a)') '         </Icon>'
      write (iout,'(a)') '      </IconStyle>'
      write (iout,'(a)') '   </Style>'
      write (iout,'(a)') '   <StyleMap id="blue">'
      write (iout,'(a)') '      <Pair>'
      write (iout,'(a)') '         <key>normal</key>'
      write (iout,'(a)') '         <styleUrl>#i_blue</styleUrl>'
      write (iout,'(a)') '      </Pair>'
      write (iout,'(a)') '      <Pair>'
      write (iout,'(a)') '         <key>highlight</key>'
      write (iout,'(a)') '         <styleUrl>#a_blue</styleUrl>'
      write (iout,'(a)') '      </Pair>'
      write (iout,'(a)') '   </StyleMap>'
      
      write (iout,'(a)') '   <Style id="i_skyblue">'
      write (iout,'(a)') '      <IconStyle>'
      write (iout,'(a)') '         <Icon>'
      write (iout,'(a)') '            <href>../tables/kml/i_skyblue.png</href>'
      write (iout,'(a)') '         </Icon>'
      write (iout,'(a)') '      </IconStyle>'
      write (iout,'(a)') '   </Style>'
      write (iout,'(a)') '   <Style id="a_skyblue">'
      write (iout,'(a)') '      <IconStyle>'
      write (iout,'(a)') '         <Icon>'
      write (iout,'(a)') '            <href>../tables/kml/a_skyblue.png</href>'
      write (iout,'(a)') '         </Icon>'
      write (iout,'(a)') '      </IconStyle>'
      write (iout,'(a)') '   </Style>'
      write (iout,'(a)') '   <StyleMap id="skyblue">'
      write (iout,'(a)') '      <Pair>'
      write (iout,'(a)') '         <key>normal</key>'
      write (iout,'(a)') '         <styleUrl>#i_skyblue</styleUrl>'
      write (iout,'(a)') '      </Pair>'
      write (iout,'(a)') '      <Pair>'
      write (iout,'(a)') '         <key>highlight</key>'
      write (iout,'(a)') '         <styleUrl>#a_skyblue</styleUrl>'
      write (iout,'(a)') '      </Pair>'
      write (iout,'(a)') '   </StyleMap>'
      
      write (iout,'(a)') '   <LookAt>'
      write (iout,'(3a)') '      <longitude>',trim(adjustl(cxlonh)),'</longitude>'
      write (iout,'(3a)') '      <latitude>',trim(adjustl(cxlath)),'</latitude>'
      write (iout,'(3a)') '      <range>',trim(adjustl(cirange)),'</range>'
      write (iout,'(a)') '      <tilt>0</tilt>'
      write (iout,'(a)') '      <heading>0</heading>'
      write (iout,'(a)') '   </LookAt>'
  
      return
      end
    

!*****************************************************************************************
      subroutine kml_coda (iout)
      
      implicit none
      
      integer iout
  
      write (iout,'(a)') '</Document>'
      write (iout,'(a)') '</kml>'
      
      return
      end
