AM_FCFLAGS =				\
	-I $(top_builddir)/src/3d/ 	\
	-I $(top_builddir)/unit-tests	\
	-I $(top_builddir)/src/utils	\
	-I $(top_builddir)/src/fft	\
	-I $(top_builddir)/src/netcdf

AM_LDFLAGS =						\
	$(top_builddir)/src/3d/parcels/*.o 		\
	$(top_builddir)/src/3d/utils/*.o		\
	$(top_builddir)/src/3d/fields/*.o		\
	$(top_builddir)/src/3d/inversion/*.o

# 21 November 2021
# gnu.org/software/automake/manual/html_node/Libtool-Convenience-Libraries.html
noinst_LTLIBRARIES = libcombi.la
libcombi_la_SOURCES =
libcombi_la_LIBADD = 					\
	$(top_builddir)/unit-tests/libunit_test.la	\
	$(top_builddir)/src/utils/libepic_utils.la	\
	$(top_builddir)/src/fft/libepic_fft.la		\
	$(top_builddir)/src/netcdf/libepic_netcdf.la

unittestsdir = $(prefix)/unit-tests

unittests_PROGRAMS = 				\
	test_jacobi_1				\
	test_jacobi_2				\
	test_jacobi_3				\
	test_jacobi_4				\
	test_ellipsoid_split 			\
	test_ellipsoid_split_merge		\
	test_parcel_init_3d			\
	test_fft_3d				\
	test_fft_2_3d				\
	test_trilinear				\
	test_laplace_correction_3d		\
	test_gradient_correction_3d		\
	test_parcel_correction_3d		\
	test_ellipsoid_merge_1			\
	test_ellipsoid_merge_2			\
	test_grid2par				\
	test_ellipsoid_reflection		\
	test_free_slip_3d			\
	test_vor2vel				\
	test_velgradg				\
	test_vtend				\
	test_lapinv0_1				\
	test_lapinv0_2				\
	test_lapinv1				\
	test_diffz0				\
	test_diffz1

dataroot_DATA =		\
	numpy_eigh.py 	\
	V.asc	   	\
	D.asc	   	\
	B.asc

test_jacobi_1_SOURCES = test_jacobi_1.f90
test_jacobi_1_LDADD = libcombi.la

test_jacobi_2_SOURCES = test_jacobi_2.f90
test_jacobi_2_LDADD = libcombi.la

test_jacobi_3_SOURCES = test_jacobi_3.f90
test_jacobi_3_LDADD = libcombi.la

test_jacobi_4_SOURCES = test_jacobi_4.f90
test_jacobi_4_LDADD = libcombi.la

test_ellipsoid_split_SOURCES = test_ellipsoid_split.f90
test_ellipsoid_split_LDADD = libcombi.la

test_ellipsoid_split_merge_SOURCES = test_ellipsoid_split_merge.f90
test_ellipsoid_split_merge_LDADD = libcombi.la

test_parcel_init_3d_SOURCES = test_parcel_init_3d.f90
test_parcel_init_3d_LDADD = libcombi.la

test_fft_3d_SOURCES = test_fft_3d.f90
test_fft_3d_LDADD = libcombi.la

test_fft_2_3d_SOURCES = test_fft_2_3d.f90
test_fft_2_3d_LDADD = libcombi.la

test_trilinear_SOURCES = test_trilinear.f90
test_trilinear_LDADD = libcombi.la

test_laplace_correction_3d_SOURCES = test_laplace_correction_3d.f90
test_laplace_correction_3d_LDADD = libcombi.la

test_gradient_correction_3d_SOURCES = test_gradient_correction_3d.f90
test_gradient_correction_3d_LDADD = libcombi.la

test_parcel_correction_3d_SOURCES = test_parcel_correction_3d.f90
test_parcel_correction_3d_LDADD = libcombi.la

test_ellipsoid_merge_1_SOURCES = test_ellipsoid_merge_1.f90
test_ellipsoid_merge_1_LDADD = libcombi.la

test_ellipsoid_merge_2_SOURCES = test_ellipsoid_merge_2.f90
test_ellipsoid_merge_2_LDADD = libcombi.la

test_ellipsoid_reflection_SOURCES = test_ellipsoid_reflection.f90
test_ellipsoid_reflection_LDADD = libcombi.la

test_grid2par_SOURCES = test_grid2par.f90
test_grid2par_LDADD = libcombi.la

test_free_slip_3d_SOURCES = test_free_slip_3d.f90
test_free_slip_3d_LDADD = libcombi.la

test_vor2vel_SOURCES = test_vor2vel.f90
test_vor2vel_LDADD = libcombi.la

test_velgradg_SOURCES = test_velgradg.f90
test_velgradg_LDADD = libcombi.la

test_vtend_SOURCES = test_vtend.f90
test_vtend_LDADD = libcombi.la

test_lapinv0_1_SOURCES = test_lapinv0_1.f90
test_lapinv0_1_LDADD = libcombi.la

test_lapinv0_2_SOURCES = test_lapinv0_2.f90
test_lapinv0_2_LDADD = libcombi.la

test_lapinv1_SOURCES = test_lapinv1.f90
test_lapinv1_LDADD = libcombi.la

test_diffz0_SOURCES = test_diffz0.f90
test_diffz0_LDADD = libcombi.la

test_diffz1_SOURCES = test_diffz1.f90
test_diffz1_LDADD = libcombi.la
