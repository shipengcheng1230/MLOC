   ! Program limits:

   integer, parameter :: nevmax = 200 ! max number of events in the cluster
   integer, parameter :: ntmax0 = 5000 ! max number TT residuals for a single event before any filtering
   integer, parameter :: ntmax1 = 35000 ! max number TT residuals used in inversion
   integer, parameter :: nmax1 = 23000 ! max number of stations in the station list (plus any auxilliary list)
   integer, parameter :: nqmax = 4000 ! max number of distinct station-phases used
   integer, parameter :: nsprdmax = 60 ! Maximum number of entries in the ttsprd file
   integer, parameter :: itmax = 4 ! Maximum number of iterations
   integer, parameter :: ndiffmax = 500 ! Maximum number of differential times
   integer, parameter :: nplot_max = 9 ! Maximum number of selected-event plots
   integer, parameter :: nhdfmax = 200 ! Maximum number of events in an HDF file for input as starting locations.
   integer, parameter :: n_supp_stn_file_max = 8 ! Maximum number of supplemental station files that can be read
   integer, parameter :: n_fault_map_nmax = 4 ! Maximum number of digital fault map files that can be read
   integer, parameter :: nitomomax = 10 ! Maximum number of tomography output files
   integer, parameter :: n_no_phreid_max = 6 ! Maximum number of phases that can be prevented from phase re-identification
   integer, parameter :: n_skip_max = 10 ! Maximum number of stations or phases that can be skipped
   integer, parameter :: n_qres_max = 400 ! Maximum number of samples allowed for estimating empirical reading error
   integer, parameter :: n_miss_sta_max = 400 ! Maximum number of missing stations for any one event or all events together
   integer, parameter :: n_ellipse_nmax = 20 ! Maximum number of individual ellipses that can plotted
   integer, parameter :: n_star_nmax = 10 ! Maximum number of individual stars that can plotted
   integer, parameter :: n_stat_nmax = 30 ! Maximum number of station locations that can plotted
   integer, parameter :: n_xsec_max = 4 ! Maximum number of cross sections that can plotted
   integer, parameter :: n_bdp_max = 40 ! Maximum number of stations in the 'bad depth phase' list
   integer, parameter :: n_tt5s_max = 20 ! Maximum number of single-station tt5 plots
   integer, parameter :: n_tt5e_max = 20 ! Maximum number of single-event tt5 plots
   integer, parameter :: n_rdpp_max = 30 ! Maximum number of relative depth phase plots
   integer, parameter :: n_epa_plot_max = 6 ! Maximum number of empirical path anomaly plots
   integer, parameter :: n_psdre_max = 20 ! Maximum number of phase-specific default reading error values
   integer, parameter :: n_ttou_max = 6 ! Maximum number of output files with TT data for specific phases
   integer, parameter :: bp_nlon = 4321 ! Number of longitude values of bouncepoint topography data file
   integer, parameter :: bp_nlat = 2161 ! Number of latitude values of bouncepoint topography data file
      
   ! Derived parameters

   integer, parameter :: mtmax = nevmax*4 ! max number of free parameters
   integer, parameter :: mtmax2 = nevmax*2
   integer, parameter :: itmax1 = itmax + 1
   integer, parameter :: ntqmax = nevmax*2 ! max number of instances of a particular station-phase
   
   real, parameter :: pi = 3.1415926536
   real, parameter :: radius = 6371.
   real, parameter :: rpd = 1.7453293e-2
   real, parameter :: xdpr = 57.29577951
   real, parameter :: dgkmla = xdpr/radius
   
   ! Depth-related
   real :: depth_inp, depth_hdf, depth_cf, depth_default, depthp_plus, depthp_minus, median_constrained_depths
   character(len=1) :: depth_inp_c, depth_hdf_c, depth_cf_c, depth_default_c
   common /depth/ depth_inp(nevmax,3), depth_hdf(nevmax,3), depth_cf(nevmax,3), depth_default,&
    depthp_plus(nevmax), depthp_minus(nevmax), median_constrained_depths,&
    depth_inp_c(nevmax), depth_hdf_c(nevmax),&
    depth_cf_c(nevmax), depth_default_c

   real :: alphah, xl1h, xl2h, alphac, xl1c, xl2c, ccv, alic, blic, kcritc, hcv, alphacg, xl1cg, xl2cg
   logical :: cv_out
   common /ellipse/ alphah, xl1h, xl2h, alphac(nevmax), xl1c(nevmax), xl2c(nevmax), ccv(nevmax,4,4), alic(nevmax),&
    blic(nevmax), kcritc(nevmax), hcv(4,4), alphacg(nevmax), xl1cg(nevmax), xl2cg(nevmax), cv_out
   
   integer :: n_ellipse
   real :: ellipse_data
   logical :: ellipse_plot
   common /ellipse2/ ellipse_plot, n_ellipse, ellipse_data(n_ellipse_nmax,5)

   integer :: n_star, iev_star
   real :: star_size
   logical :: star_plot
   common /star/ star_plot, n_star, iev_star(n_star_nmax), star_size(n_star_nmax)
   
   integer :: n_stat
   real :: stat_data
   logical :: stat_plot
   common /station_plot/ stat_plot, n_stat, stat_data(n_stat_nmax,3)
   
   character(len=80) :: basename, version
   character(len=8) :: mloc_author
   common /title/ basename, version, mloc_author

   integer :: jb, ntq, ntqi
   real :: vnhat, sighatj, wq, wq2
   common /wq/ wq(nqmax), wq2(nqmax), jb(ntmax1), vnhat(ntmax1), sighatj(ntmax1), ntq(ntqmax,nqmax), ntqi(nqmax)

   integer :: nev, mtiev, mindx
   character(len=10) :: evid
   character(len=30) :: evtnam
   character(len=6) :: evid_source
   character(len=76) :: event_comment
   common /event/ nev, mtiev(nevmax), mindx(nevmax,4), evtnam(nevmax), evid(nevmax),&
    evid_source(nevmax), event_comment(nevmax)

   real :: a
   common /deriv/ a(nevmax,ntmax0,4)

   integer :: ndat, ndatdl, ndatpr, ndatfl, nmiss, ndatc, n_miss_sta, n_miss_sta_total
   character(len=20) :: missta, n_miss_sta_list
   common /ndat/ ndat(nevmax,0:itmax1), ndatdl(nevmax,0:itmax1), ndatpr(nevmax,0:itmax1), ndatfl(nevmax), nmiss(nevmax),&
    missta(nevmax,n_miss_sta_max), ndatc(nevmax,0:itmax1), n_miss_sta_total, n_miss_sta(n_miss_sta_max),&
    n_miss_sta_list(n_miss_sta_max)

   integer :: kcode
   real :: tto, dt, delt, azes, azse, psd, stacs, sdstcs, elcr, s, rel_depth_phase, rdp_res
   common /station/ tto(nevmax,ntmax0), dt(nevmax,ntmax0,0:itmax1), delt(nevmax,ntmax0), azes(nevmax,ntmax0),&
    azse(nevmax,ntmax0), psd(nevmax,ntmax0), stacs(nevmax,ntmax0), sdstcs(nevmax,ntmax0), elcr(nevmax,ntmax0),&
    s(nevmax,ntmax0,0:itmax1), kcode(nevmax,ntmax0), rel_depth_phase(nevmax,ntmax0), rdp_res(nevmax,ntmax0,0:itmax1)
   
   integer :: ittabl, istcor, iwght, ibias, longitude_range
   common /flag/ ittabl, istcor, iwght, ibias, longitude_range

   integer :: hourh, hourp, minh, minp, nstep
   real :: lath, lonh, latp, lonp, depthh, sech, lat_cf, lat_hdf, lon_cf, lon_hdf, time_cf, time_hdf, depthp, secp
   real :: latpgc ! event latitude in geocentric coordinates.
   real :: lathgc ! hypocentroid latitude in geocentric coordinates.
   real :: hlatshift, hlonshift, hdepthshift, htimeshift
   logical :: latf, lonf, depthf, timef, latfh, lonfh, depthfh, timefh
   character(len=1) :: depset_pr
   character(len=5) :: fixpr, fixprh
   common /mloc0/ lath(0:itmax1), lathgc(0:itmax1), lonh(0:itmax1), depthh(0:itmax1), hourh(0:itmax1),&
    minh(0:itmax1), sech(0:itmax1), lat_cf(nevmax), lat_hdf(nevmax), lon_cf(nevmax), lon_hdf(nevmax), time_cf(nevmax),&
    time_hdf(nevmax), latp(nevmax,0:itmax1), latpgc(nevmax,0:itmax1), lonp(nevmax,0:itmax1), depthp(nevmax,0:itmax1),&
    hourp(nevmax,0:itmax1), minp(nevmax,0:itmax1), secp(nevmax,0:itmax1), latf(nevmax), lonf(nevmax), depthf(nevmax),&
    timef(nevmax), hlatshift, hlonshift, hdepthshift,&
    htimeshift, latfh, lonfh, depthfh, timefh, nstep, depset_pr(nevmax), fixprh(4), fixpr(nevmax,4)
   
   integer :: nqmth
   real :: delx0, sdxhath, ehatsqh, shath, dtmph, ehiev, dimph, bcorr
   common /hypcntr/ delx0(4,0:itmax1), sdxhath(4), ehatsqh(0:itmax1), shath(0:itmax1), dtmph(nevmax,ntmax0),&
    ehiev(nevmax,0:itmax1), dimph(12), bcorr(4), nqmth(0:itmax1)

   integer :: ntc, nqc
   real :: dxp, sdxhatc, ehatsqc, shatc, dtmpc, eciev, dimpc, eci, dimpciev, shatsqci, radius_cvff, tikhonov_factor
   double precision :: vhatc, shatsqc
   common /cluster/ dxp(nevmax,4,0:itmax1), sdxhatc(nevmax,4), ehatsqc(0:itmax1), shatc(0:itmax1), dtmpc(nevmax,ntmax0),&
    eciev(nevmax,0:itmax1), dimpc(12), eci(nevmax,ntmax0), ntc, nqc, shatsqc, shatsqci(nevmax), vhatc(mtmax,mtmax),&
    dimpciev(nevmax,12), radius_cvff, tikhonov_factor

   integer :: io_in, io_out, io_log, io_pdf, io_sout, io_gmt, io_xdat, io_locmod, io_taup, io_cal, io_cfil,&
    io_rhdf, io_depth_phase, io_ttsprd, io_plog, io_stn_log, io_bdp, io_bayes, io_tt
   character(len=20) :: infile20
   character(len=100) :: infile, outfile, logfile, datadir, mloc_path
   common /io/ io_in, infile(nevmax), io_out, outfile, io_log, logfile, io_pdf, io_sout, infile20(nevmax),&
    io_gmt, io_xdat, io_locmod, io_taup, io_cal, io_cfil, io_rhdf, io_depth_phase, io_ttsprd, io_plog,&
    io_stn_log, io_bdp, io_bayes, io_tt, datadir, mloc_path

   real :: wind1, wind2, windloclim
   common /windlimits/ wind1, wind2, windloclim

   integer :: nstat1, n_supp_stn_file, jdate_on, jdate_off, kstat, nkstat, n_supp_stn, n_failed_date_range, nstn_used
   real :: stalat, stalon, ahgtr, sca0, sd1
   logical :: suppstn, nsmd, read_ad
   character(len=8) :: sta_author, sta_deployment
   character(len=5) :: nstr1, sta_agency
   character(len=13) :: stn_dcal_used 
   character(len=40) :: station_master, nsmd_file
   character(len=100) :: suppfilnam
   common /st1/ nstr1(nmax1), stalat(nmax1), stalon(nmax1), ahgtr(nmax1), sta_author(nmax1), jdate_on(nmax1), jdate_off(nmax1),&
    sca0(nmax1), sd1(nmax1), kstat(nmax1), nstat1, n_failed_date_range, nkstat, suppstn, n_supp_stn_file,&
    n_supp_stn, suppfilnam(n_supp_stn_file_max), station_master, stn_dcal_used(nmax1), nstn_used, nsmd, nsmd_file,&
    sta_agency(nmax1), sta_deployment(nmax1), read_ad
   
   integer :: n_dupe, n_dupe_pure, n_dupe_minor, n_dupe_significant, n_dupe_conflict
   character(len=40) :: duplication(nmax1)
   common /duplicate/ n_dupe, n_dupe_pure, n_dupe_minor, n_dupe_significant, n_dupe_conflict, duplication
               
   real :: stladg, stlndg, ahgts, sd1s, sd2s, sca0s, sca1s, scb1s, sca2s, scb2s
   common /stcorri/ stladg(nevmax,ntmax0), stlndg(nevmax,ntmax0), ahgts(nevmax,ntmax0), sd1s(nevmax,ntmax0),&
    sd2s(nevmax,ntmax0), sca0s(nevmax,ntmax0), sca1s(nevmax,ntmax0), scb1s(nevmax,ntmax0), sca2s(nevmax,ntmax0),&
    scb2s(nevmax,ntmax0)

   integer :: nevt, mone, idye, iyre, hour_inp, min_inp
   real :: sec_inp, lat_inp, lon_inp, time_inp, rmag
   character(len=2) :: mmag
   character(len=8) :: hypo_author, magnitude_author
   common /iscepi/ nevt(nevmax), mone(nevmax), idye(nevmax), iyre(nevmax), hour_inp(nevmax), min_inp(nevmax),&
    sec_inp(nevmax), lat_inp(nevmax), lon_inp(nevmax), time_inp(nevmax), rmag(nevmax), mmag(nevmax), hypo_author(nevmax),&
    magnitude_author(nevmax)

   integer :: nst, ipay, ipamo, ipad, ipah, ipam, iptim, mnf_line
   real :: pas, resisc, sdread, ttsprd, ttoff, weight
   logical :: dflag ! if true, use data flags (x, p, d, k, l)
   logical :: rel_phase ! Relative phase flag
   character(len=1) :: fcode
   character(len=3) :: channel
   character(len=5) :: stname, agency
   character(len=8) :: phase, readsrc, phase0, deployment
   character(len=10) :: dist_az
   character(len=20) :: sad
   character(len=47) :: adslcaa
   common /iscphas/ nst(nevmax), mnf_line(nevmax,ntmax0), stname(nevmax,0:ntmax0),ipay(nevmax,0:ntmax0), ipamo(nevmax,0:ntmax0),&
    ipad(nevmax,ntmax0), ipah(nevmax,ntmax0), ipam(nevmax,ntmax0), pas(nevmax,ntmax0),&
    phase(nevmax,ntmax0), resisc(nevmax,ntmax0), sdread(nevmax,ntmax0), ttsprd(nevmax,ntmax0), ttoff(nevmax,ntmax0),&
    fcode(nevmax,ntmax0), channel(nevmax,ntmax0), iptim(nevmax,ntmax0), readsrc(nevmax,ntmax0),&
    dflag, phase0(nevmax,ntmax0), adslcaa(nevmax,ntmax0), agency(nevmax,ntmax0), deployment(nevmax,ntmax0),&
    sad(nevmax,ntmax0), weight(nevmax,ntmax0), rel_phase(nevmax,ntmax0), dist_az(nevmax,ntmax0)

   character(len=3) :: weigpr, latpr, lonpr, biaspr, data_weight_pr, dflag_pr, pttt_pr
   character(len=5) :: depthpr
   character(len=11) :: timepr
   character(len=20) :: tablpr
   character(len=80) :: stacpr
   common /print/ tablpr(2), stacpr(2), weigpr(0:1), latpr(nevmax), lonpr(nevmax), depthpr(nevmax), timepr(nevmax),&
    biaspr(0:1), data_weight_pr, dflag_pr, pttt_pr

   character(len=1) :: dirsym
   character(len=80) :: taup_path, ellip_path, station_path, cpt_path, dem_path
   common /directory/ taup_path, ellip_path, station_path, cpt_path, dem_path, dirsym

   ! Calibration data
   integer :: ncal, cal_min, cal_hr, icaltype
   real :: cal_lat, cal_lon, cal_dep,cal_sec, otsr, rcv, scv, lath_epa, lonh_epa, depthh_epa,&
    latp_cal, lonp_cal, depthp_cal, otsp_cal,del_cal_lat, del_cal_lon, del_cal_dep, del_cal_tim,&
    s_cal, dt_cal, delt_cal, azes_cal, azse_cal, elcr_cal, stacs_cal, sdstcs_cal, w12, w3, w4, accv, sgcv, cal_level
   logical :: calibration, cal_event, ot_cal
   character(len=12) :: icaltype_pr
   character(len=4) :: cal_code_input
   character(len=1) :: cal_par
   common /cal/ ncal(3), cal_min(nevmax,3), cal_hr(nevmax,3), cal_lat(nevmax,3), cal_lon(nevmax,3), cal_dep(nevmax,3),&
    cal_sec(nevmax,3), otsr(nevmax), rcv(nevmax,3,4,4), scv(nevmax,4,4), sgcv(nevmax,4,4), latp_cal(nevmax), lonp_cal(nevmax),&
    depthp_cal(nevmax), otsp_cal(nevmax),lath_epa, lonh_epa, depthh_epa, del_cal_lat, del_cal_lon, del_cal_dep, del_cal_tim,&
    s_cal(nevmax,ntmax0), dt_cal(nevmax,ntmax0), delt_cal(nevmax,ntmax0), azes_cal(nevmax,ntmax0), azse_cal(nevmax,ntmax0),&
    elcr_cal(nevmax,ntmax0), stacs_cal(nevmax,ntmax0), sdstcs_cal(nevmax,ntmax0), w12(nevmax), w3(nevmax), w4(nevmax),&
    accv(nevmax,4,4), cal_event(nevmax,3), calibration, icaltype, icaltype_pr, cal_par(nevmax), ot_cal, cal_code_input(nevmax),&
    cal_level
   
   integer :: io_rderr, indexq, nqname, iqiev
   real :: rdsigma, qres, rderr0, rderr_min, rderr_min_loc, rderr_min_depth, qlat, qlon, qelev,&
    rderr_loc_p, rderr_loc_s, rderr_loc_delt, epa, ere
   logical :: read_rderr, rels_set
   character(len=21) :: qname, qname1
   character(len=100) :: rderrfname
   common /rderr/ qname(nqmax), qname1(nqmax), indexq(nqmax), qres(nqmax,n_qres_max), iqiev(nqmax,n_qres_max),&
    rdsigma(nqmax), rderr0(nqmax), qlat(nqmax), qlon(nqmax), qelev(nqmax), epa(nqmax), ere(nqmax), nqname,&
    io_rderr, rderr_min, rderr_min_loc, rderr_min_depth, rderrfname, read_rderr, rels_set, rderr_loc_p,&
    rderr_loc_s, rderr_loc_delt
   
   integer :: n_no_phreid
   logical :: phid ! .true. if phase identification is to be done
   logical :: phidird ! Individual reading flag to allow change of phase ID if .true.
   character(len=1) :: no_phid ! Flag to prevent phase reidentification for particular readings
   character(len=8) :: no_phreid
   common /pdf/ phid, phidird(nevmax,ntmax0), no_phreid(n_no_phreid_max), n_no_phreid, no_phid(nevmax,ntmax0)
   
   logical :: pukeout ! .true. if PUKE file should be created.
   common /puke/ pukeout

   logical :: blocout ! .true. if BAYESLOC file should be created.
   real :: eqr
   common /bloc/ blocout, eqr(nevmax)
   
   logical :: sp_out ! .true. if an S-P file should be created.
   common /sminusp/ sp_out
   
   integer :: io_datf, io_dat0 ! logical unit for writing the .datf and .dat0 files.
   logical :: datfout ! .true. if .datf file should be created.
   common /datf/ datfout, io_dat0, io_datf
   
   integer :: io_comcat
   logical :: comcatout
   character(len=132) :: ccat_folder
   character(len=80) :: commentary_fname
   common /comcat/ comcatout, io_comcat, commentary_fname, ccat_folder
   
   integer :: io_lres ! logical unit for LRES file.
   real :: lres ! cutoff value for cluster residual.
   logical :: lresout ! .true. if LRES file should be created.
   common /lres/ lresout, lres, io_lres
   
   character(len=20) :: annotation
   character(len=180) :: hdfline
   common /hdfline/ hdfline(nevmax), annotation(nevmax)
   
   logical :: verbose_screen, verbose_log, debug
   common /verbose/ verbose_screen, verbose_log, debug
   
   real :: rdbt, alphagt, s1gt, s2gt, depgt, otgt, gtlevel
   common /gtshift/ rdbt, alphagt, s1gt, s2gt, depgt, otgt, gtlevel
   
   integer :: indph
   real :: dlimlocmod, zlimlocmod
   logical :: locmod
   character(len=100) :: locmodfname
   common /localmodel/ locmod, dlimlocmod, zlimlocmod, indph, locmodfname
   
   integer :: nclim, nhlim
   real :: clim, hlim
   common /chlim/ nclim, nhlim, clim(3,2), hlim(3,2)
   
   logical :: fltrcdelt, fltrcres, fltrcflag, fltrc, fltrhdelt, fltrhres, fltrhflag, fltrh, connected
   common /fltr/ fltrcdelt(nevmax,ntmax0), fltrcres(nevmax,ntmax0), fltrcflag(nevmax,ntmax0), fltrc(nevmax,ntmax0),&
    fltrhdelt(nevmax,ntmax0), fltrhres(nevmax,ntmax0), fltrhflag(nevmax,ntmax0), fltrh(nevmax,ntmax0),&
    connected(nevmax,ntmax0)
   
   logical :: ponly
   common /phyp/ ponly
   
   integer :: nsprd
   real :: sprd, ttoffset
   logical :: read_ttsprd
   character(len=8) :: sprdph
   character(len=100) :: ttsprdfname
   common /sprd/ nsprd, sprdph(nsprdmax), sprd(nsprdmax), ttoffset(nsprdmax), ttsprdfname, read_ttsprd
   
   integer :: ishell, gmt_version
   logical :: plot, tt1, tt2, tt3, tt4, tt5, tt6, tt7, tt8, tt9, reduced, vectors(4), splt, eplt, md_out, fdhp
   character(len=132) :: gmt_script_dir
   common /gmt/ gmt_version, plot(nevmax,0:nplot_max), tt1, tt2, tt3, tt4, tt5, tt6, tt7, tt8, tt9, reduced, vectors,&
    ishell, splt, eplt, fdhp, md_out, gmt_script_dir
   
   integer :: ndiff, idiff, idiff0, io_diffdat, diff_line
   logical :: diffdat
   character(len=80) :: diffdatfilnam
   common /diff/ ndiff, idiff(nevmax,ntmax0), diff_line(nevmax,ntmax0), idiff0(nqmax), diffdat, io_diffdat, diffdatfilnam
   
   real :: vmr, cscale
   common /cvtt/ vmr, cscale
   
   integer :: tt_corr
   logical :: bias_corr ! Bias correction of hypocentroid
   logical :: data_weight ! Weighting data by inverse of error
   logical :: pttt ! perfect theoretical travel times assumed for hypocentroid
   common /control/ bias_corr, data_weight, tt_corr, pttt
   
   character(len=20) :: taup_model
   common /taup_mod/ taup_model
   
   integer :: nrhdf, hdf_yr4, hdf_mon, hdf_day
   real :: hdf_lat, hdf_lon, hdf_dep, hdf_dep_plus, hdf_dep_minus, hdf_time
   logical :: read_hdf
   character(len=10) :: hdf_evid
   character(len=60) :: rhdf_filnam
   character(len=2) :: hdf_dep_code
   common /hdf/ read_hdf, nrhdf, hdf_yr4(nhdfmax), hdf_mon(nhdfmax), hdf_day(nhdfmax), hdf_lat(nhdfmax), hdf_lon(nhdfmax),&
    hdf_dep(nhdfmax), hdf_dep_plus(nhdfmax), hdf_dep_minus(nhdfmax), hdf_time(nhdfmax), hdf_evid(nhdfmax), rhdf_filnam,&
    hdf_dep_code(nhdfmax)
   
   ! Lg travel time calculations
   real :: lg_a, lg_b, lg_min
   logical :: lg_out
   common /lg/ lg_out, lg_a, lg_b, lg_min
   
   ! T-phase data
   real :: tphase_a, tphase_b
   logical :: tphase_out
   common /tp/ tphase_out, tphase_a, tphase_b
   
   ! Digital fault map
   integer :: n_fault_map
   logical :: fault_map
   character(len=30) :: psxy_wsf0, psxy_wsf
   character(len=100) :: fault_map_filename
   common /fault/ fault_map, fault_map_filename(n_fault_map_nmax), n_fault_map, psxy_wsf0, psxy_wsf(n_fault_map_nmax)
   
   ! Direct calibration
   real :: alphadc, xl1dc, xl2dc, ddepdc, dotdc
   logical :: direct_cal
   character(len=3) :: dcal_pr
   common /dcal/ direct_cal, alphadc(nevmax), xl1dc(nevmax), xl2dc(nevmax), ddepdc(nevmax), dotdc(nevmax), dcal_pr
   
   ! Plot topography in GMT
   logical :: plot_dem1, plot_dem2, plot_globe, plot_etopo1, plot_gina
   character(len=132) :: dem2_filename
   character(len=40) :: cpt_file
   common /topography/ plot_dem1, plot_dem2, plot_globe, plot_etopo1, plot_gina, dem2_filename, cpt_file
   
   ! Tomography output files
   integer :: itomo, nitomo
   character(len=8) :: tomo_phase
   common /tomography/ tomo_phase(nitomomax), itomo(nitomomax), nitomo
   
   ! Skip specified stations, phases and authors
   integer :: n_skip
   logical :: skip
   character(len=8) :: skip_params
   common /skipsta/ skip, skip_params(n_skip_max,3), n_skip
      
   ! Cross sections
   integer :: n_xsec
   real :: xsec_lat1, xsec_lon1, xsec_lat2, xsec_lon2, xsec_depth, xsec_width
   common /xsec/ n_xsec, xsec_lat1(n_xsec_max), xsec_lon1(n_xsec_max), xsec_lat2(n_xsec_max), xsec_lon2(n_xsec_max),&
    xsec_depth(n_xsec_max), xsec_width(n_xsec_max)
   
   ! Focal Mechanisms
   integer :: io_focal_mech
   logical :: focal_mech, focal_mech_iev
   character(len=100) :: focal_mech_line
   common /focal_mechanism/ io_focal_mech, focal_mech, focal_mech_iev(nevmax), focal_mech_line(nevmax)
   
   ! Kill events
   integer :: kill_count
   logical :: kill_one, kill_all
   character(len=76) :: kill_reason
   common /kill_event/ kill_count, kill_one, kill_all, kill_reason
   
   ! Command arguments
   integer :: nvar
   real :: value
   character(len=76) :: params
   common /command_args/ params, value(12), nvar
   
   ! Convergence limits
   real :: cl_epi_h, cl_epi_c, cl_dep_h, cl_dep_c, cl_ot_h, cl_ot_c
   integer :: convergence_test_index
   common /conv_limits/ cl_epi_h, cl_epi_c, cl_dep_h, cl_dep_c, cl_ot_h, cl_ot_c, convergence_test_index
   
   ! List of stations suspected of reporting bogus depth phases
   logical :: bdp_list
   integer :: n_bdp
   character(len=108) :: bdp_station
   character(len=100) :: bdp_filnam
   common /bogus/ bdp_list, bdp_filnam, n_bdp, bdp_station(n_bdp_max)
   
   ! Single-station tt5 plots
   logical :: tt5s
   integer :: n_tt5s
   character(len=5) :: tt5s_sta
   common /tt5s_plot/ tt5s, n_tt5s, tt5s_sta(n_tt5s_max)

   ! Single-event tt5 plots
   logical :: tt5e
   integer :: n_tt5e
   character(len=30) :: tt5e_evt
   common /tt5e_plot/ tt5e, n_tt5e, tt5e_evt(n_tt5e_max)
   
   ! Relative depth phases
   logical :: rdpp
   integer :: n_rdpp, nrdp
   character(len=30) :: rdpp_evt
   common /rdpp_plot/ rdpp, n_rdpp, rdpp_evt(n_rdpp_max), nrdp
   
   ! Empirical path anomaly plots
   logical :: epa_plot, epa_plot_raypath
   integer :: n_epa_plot
   real :: epa_plot_distance
   character(len=8) :: epa_plot_phase
   common /epa_plots/ epa_plot, n_epa_plot, epa_plot_phase(n_epa_plot_max), epa_plot_distance(n_epa_plot_max),&
    epa_plot_raypath(n_epa_plot_max)

   ! Selection of a direct calibration subcluster
   integer :: subc_nmin, subc_nconnect
   real :: subc_delt
   logical :: subc_set
   common /subc/ subc_set, subc_delt, subc_nmin, subc_nconnect
   
   ! Phase-specific default reading error values
   character(len=8) :: psdre_phase
   real :: psdre
   integer :: n_psdre
   common /readerr_default/ n_psdre, psdre_phase(n_psdre_max), psdre(n_psdre_max)
   
   ! Special output file for limited distance range (usually used for Pg/Pn cross-over range)
   real :: dist1, dist2
   logical :: oldr_out
   integer :: io_oldr
   common /oldr/ dist1, dist2, oldr_out, io_oldr
   
   ! Bounce point topography
   integer*2 :: bp_topo
   logical :: bptc
   common /bp/ bp_topo(bp_nlon,bp_nlat), bptc
   
   ! Station elevation correction
   real :: pcrvel, scrvel
   common /sec/ pcrvel, scrvel
   
   ! Empirical TT output files
   integer :: n_ttou
   character(len=8) :: ttou_phase
   logical :: ttou
   common /ttou_files/ ttou, n_ttou, ttou_phase(n_ttou_max)
