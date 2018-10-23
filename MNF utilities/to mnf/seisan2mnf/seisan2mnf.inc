   common /eventline/ event_usage, evid, event_geog
   character(len=100) :: event_geog
   character(len=10) :: evid
   character(len=1) :: event_usage

   common /hypocenterline/ hypo_usage, hypo_year, hypo_month, hypo_day, hypo_hour, hypo_minute,&
    hypo_seconds, latitude, longitude, depth, hypo_author, hypo_orid
   character(len=10) :: hypo_orid
   character(len=8) :: hypo_author
   character(len=1) :: hypo_usage
   integer :: hypo_year, hypo_month, hypo_day, hypo_hour, hypo_minute
   real :: hypo_seconds, latitude, longitude, depth
   
   common /magnitudeline/ magnitude_usage, magnitude, magnitude_scale, magnitude_author, magnitude_orid
   character(len=10) :: magnitude_orid
   character(len=8) :: magnitude_author
   character(len=5) :: magnitude_scale
   character(len=1) :: magnitude_usage
   real :: magnitude

   common /phaseline/ phase_usage, station_flag, station_in, station_out, distance_pr, azeq_pr, phase_flag,&
    phase, phase_year, phase_month, phase_day, phase_hour, phase_minute, phase_hour_p, phase_minute_p, phase_seconds,&
    phase_seconds_p, iptim_pr, residual_pr, agency, deployment, location, channel, phase_author, arrid_pr

   character(len=10) :: arrid_pr
   character(len=8) :: phase_author, phase, deployment
   character(len=6) :: station_in, station_out, distance_pr
   character(len=5) :: residual_pr, agency
   character(len=3) :: channel, azeq_pr
   character(len=2) :: iptim_pr, location
   character(len=1) :: phase_usage, station_flag, phase_flag
   integer :: phase_year, phase_month, phase_day, phase_hour, phase_minute, phase_hour_p, phase_minute_p
   real :: phase_seconds, phase_seconds_p
