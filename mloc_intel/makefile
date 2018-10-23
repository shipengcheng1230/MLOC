OBJFILES=   cal_shift.o \
			dcal.o \
			dsvd2.o \
			hyposat_loc.o \
			libtau.o \
			mloc_commands.o \
            mloc.o \
            mlocinv.o \
            mlocio_mnf.o \
            mloclib_date_time.o \
            mloclib_geog.o \
            mloclib_gmt.o \
            mloclib_inv.o \
            mloclib_messages.o \
            mloclib_phases.o \
            mloclib_set.o \
            mloclib_stations.o \
            mloclib_statistics.o \
            mloclib_tt.o \
            mloclib.o \
            mlocout_bloc.o \
            mlocout_comcat.o \
            mlocout_cv.o \
            mlocout_gmt.o \
            mlocout_gmt5.o \
            mlocout_hdf.o \
            mlocout_kml.o \
            mlocout_phase_data.o \
            mlocout_puke.o \
            mlocout_rderr.o \
            mlocout_sp.o \
            mlocout_summary.o \
            mlocout_tomo.o \
            mlocout_ttsprd.o \
            mlocset.o

mloc_i: $(OBJFILES)
	ifort $(OBJFILES) -o mloc_i

FFLAGS=  -c  -warn all -save   -fp-stack-check
# -c (compile but don't link)
# -W (warn more verbosely than usual)
# -Wall (enable commonly used warnings)
# NOTE: libtau needs an extra compiler flag concerning assumed byte record length

# Location of source files
SRCDIR=../mloc_src

cal_shift.o: $(SRCDIR)/cal_shift.f90 $(SRCDIR)/mloc.inc $(SRCDIR)/cal_shift.inc
	ifort $(FFLAGS) $(SRCDIR)/cal_shift.f90

dcal.o: $(SRCDIR)/dcal.f90 $(SRCDIR)/mloc.inc 
	ifort $(FFLAGS) $(SRCDIR)/dcal.f90

dsvd2.o: $(SRCDIR)/dsvd2.f90 $(SRCDIR)/mloc.inc 
	ifort $(FFLAGS) $(SRCDIR)/dsvd2.f90

hyposat_loc.o: $(SRCDIR)/hyposat_loc.f90 $(SRCDIR)/hyposat_loc.inc 
	ifort $(FFLAGS) $(SRCDIR)/hyposat_loc.f90

libtau.o: $(SRCDIR)/libtau.f90 $(SRCDIR)/ttlim.inc 
	ifort $(FFLAGS) -assume byterecl $(SRCDIR)/libtau.f90

mloc_commands.o: $(SRCDIR)/mloc_commands.f90 $(SRCDIR)/mloc.inc 
	ifort $(FFLAGS) $(SRCDIR)/mloc_commands.f90

mloc.o: $(SRCDIR)/mloc.f90 $(SRCDIR)/mloc.inc $(SRCDIR)/libtau.inc $(SRCDIR)/ttlim.inc
	ifort $(FFLAGS) $(SRCDIR)/mloc.f90

mlocinv.o: $(SRCDIR)/mlocinv.f90 $(SRCDIR)/mloc.inc 
	ifort $(FFLAGS) $(SRCDIR)/mlocinv.f90

mlocio_diff.o: $(SRCDIR)/mlocio_diff.f90 $(SRCDIR)/mloc.inc 
	ifort $(FFLAGS) $(SRCDIR)/mlocio_diff.f90

mlocio_mnf.o: $(SRCDIR)/mlocio_mnf.f90 $(SRCDIR)/mloc.inc 
	ifort $(FFLAGS) $(SRCDIR)/mlocio_mnf.f90

mloclib_date_time.o: $(SRCDIR)/mloclib_date_time.f90 $(SRCDIR)/mloc.inc 
	ifort $(FFLAGS) $(SRCDIR)/mloclib_date_time.f90

mloclib_geog.o: $(SRCDIR)/mloclib_geog.f90 $(SRCDIR)/mloc.inc 
	ifort $(FFLAGS) $(SRCDIR)/mloclib_geog.f90

mloclib_gmt.o: $(SRCDIR)/mloclib_gmt.f90 $(SRCDIR)/mloc.inc 
	ifort $(FFLAGS) $(SRCDIR)/mloclib_gmt.f90

mloclib_inv.o: $(SRCDIR)/mloclib_inv.f90 $(SRCDIR)/mloc.inc 
	ifort $(FFLAGS) $(SRCDIR)/mloclib_inv.f90

mloclib_messages.o: $(SRCDIR)/mloclib_messages.f90 $(SRCDIR)/mloc.inc 
	ifort $(FFLAGS) $(SRCDIR)/mloclib_messages.f90

mloclib_phases.o: $(SRCDIR)/mloclib_phases.f90 $(SRCDIR)/mloc.inc 
	ifort $(FFLAGS) $(SRCDIR)/mloclib_phases.f90

mloclib_set.o: $(SRCDIR)/mloclib_set.f90 $(SRCDIR)/mloc.inc 
	ifort $(FFLAGS) $(SRCDIR)/mloclib_set.f90

mloclib_stations.o: $(SRCDIR)/mloclib_stations.f90 $(SRCDIR)/mloc.inc 
	ifort $(FFLAGS) $(SRCDIR)/mloclib_stations.f90

mloclib_statistics.o: $(SRCDIR)/mloclib_statistics.f90 $(SRCDIR)/mloc.inc 
	ifort $(FFLAGS) $(SRCDIR)/mloclib_statistics.f90

mloclib_tt.o: $(SRCDIR)/mloclib_tt.f90 $(SRCDIR)/mloc.inc 
	ifort $(FFLAGS) $(SRCDIR)/mloclib_tt.f90

mloclib.o: $(SRCDIR)/mloclib.f90 $(SRCDIR)/mloc.inc $(SRCDIR)/ellip.inc
	ifort $(FFLAGS) $(SRCDIR)/mloclib.f90

mlocout_bloc.o: $(SRCDIR)/mlocout_bloc.f90 $(SRCDIR)/mloc.inc 
	ifort $(FFLAGS) $(SRCDIR)/mlocout_bloc.f90

mlocout_comcat.o: $(SRCDIR)/mlocout_comcat.f90 $(SRCDIR)/mloc.inc 
	ifort $(FFLAGS) $(SRCDIR)/mlocout_comcat.f90

mlocout_cv.o: $(SRCDIR)/mlocout_cv.f90 $(SRCDIR)/mloc.inc 
	ifort $(FFLAGS) $(SRCDIR)/mlocout_cv.f90

mlocout_gmt.o: $(SRCDIR)/mlocout_gmt.f90 $(SRCDIR)/mloc.inc 
	ifort $(FFLAGS) $(SRCDIR)/mlocout_gmt.f90

mlocout_gmt5.o: $(SRCDIR)/mlocout_gmt5.f90 $(SRCDIR)/mloc.inc 
	ifort $(FFLAGS) $(SRCDIR)/mlocout_gmt5.f90

mlocout_hdf.o: $(SRCDIR)/mlocout_hdf.f90 $(SRCDIR)/mloc.inc 
	ifort $(FFLAGS) $(SRCDIR)/mlocout_hdf.f90

mlocout_kml.o: $(SRCDIR)/mlocout_kml.f90 $(SRCDIR)/mloc.inc 
	ifort $(FFLAGS) $(SRCDIR)/mlocout_kml.f90

mlocout_phase_data.o: $(SRCDIR)/mlocout_phase_data.f90 $(SRCDIR)/mloc.inc 
	ifort $(FFLAGS) $(SRCDIR)/mlocout_phase_data.f90

mlocout_puke.o: $(SRCDIR)/mlocout_puke.f90 $(SRCDIR)/mloc.inc 
	ifort $(FFLAGS) $(SRCDIR)/mlocout_puke.f90

mlocout_rderr.o: $(SRCDIR)/mlocout_rderr.f90 $(SRCDIR)/mloc.inc 
	ifort $(FFLAGS) $(SRCDIR)/mlocout_rderr.f90

mlocout_sp.o: $(SRCDIR)/mlocout_sp.f90 $(SRCDIR)/mloc.inc 
	ifort $(FFLAGS) $(SRCDIR)/mlocout_sp.f90

mlocout_summary.o: $(SRCDIR)/mlocout_summary.f90 $(SRCDIR)/mloc.inc 
	ifort $(FFLAGS) $(SRCDIR)/mlocout_summary.f90

mlocout_tomo.o: $(SRCDIR)/mlocout_tomo.f90 $(SRCDIR)/mloc.inc 
	ifort $(FFLAGS) $(SRCDIR)/mlocout_tomo.f90

mlocout_ttsprd.o: $(SRCDIR)/mlocout_ttsprd.f90 $(SRCDIR)/mloc.inc 
	ifort $(FFLAGS) $(SRCDIR)/mlocout_ttsprd.f90

mlocset.o: $(SRCDIR)/mlocset.f90 $(SRCDIR)/mloc.inc 
	ifort $(FFLAGS) $(SRCDIR)/mlocset.f90

clean:
	rm mloc_i *.o


