module mnf_def

   save

   character(len=100) :: event_annotation
   character(len=1) :: event_usage

   character(len=10) :: hypo_orid
   character(len=8) :: hypo_author
   character(len=5) :: smaj, smin, ot_err
   character(len=3) :: iaz
   character(len=1) :: hypo_usage, depth_code
   integer :: hypo_year, hypo_month, hypo_day, hypo_hour, hypo_minute
   real :: hypo_seconds, latitude, longitude, depth, depth_err
   
   character(len=10) :: magnitude_orid
   character(len=8) :: magnitude_author
   character(len=5) :: magnitude_scale
   character(len=1) :: magnitude_usage
   real :: magnitude

   character(len=10) :: arrid_pr
   character(len=8) :: phase_author, phase, deployment
   character(len=6) :: station, distance_pr
   character(len=5) :: residual_pr, agency
   character(len=3) :: channel, azeq_pr
   character(len=2) :: iptim_pr, location
   character(len=1) :: phase_usage, station_flag, phase_flag
   integer :: phase_year, phase_month, phase_day, phase_hour, phase_minute
   real :: phase_seconds
   
   character(len=10) :: evid
   character(len=6) :: version = '1.3.3 '
   
end module mnf_def