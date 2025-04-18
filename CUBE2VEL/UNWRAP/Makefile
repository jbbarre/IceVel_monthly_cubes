
all: tre

one:
	rm *.o *.so

two:
	gcc -shared -fPIC -O3 -c unw2vel_3d_sp_wrap.c
	gfortran -shared -O3 -fPIC -c unw2vel_3d_sp.f90
	gcc -shared -lgfortran -O3 -o unw2vel_3d_sp.so unw2vel_3d_sp.o unw2vel_3d_sp_wrap.o
tow:
	gcc -shared -fPIC -O3 -c unw2vel_sp_wrap.c
	gfortran -shared -O3 -fPIC -c unw2vel_sp.f90
#	gfortran -shared -fPIC -O3 -c SORT.f
#	gfortran -shared -fPIC -O3 -c MEDIAN.f
	gcc -shared -lgfortran -O3 -o unw2vel_sp.so unw2vel_sp.o unw2vel_sp_wrap.o
tre:
	gcc -shared -fPIC -O3 -c unw2vel_feathering_sp_wrap.c
	gfortran -shared -O3 -fPIC -c unw2vel_feathering_sp.f90
	gcc -shared -lgfortran -O3 -o unw2vel_feathering_sp.so unw2vel_feathering_sp_wrap.o unw2vel_feathering_sp.o

