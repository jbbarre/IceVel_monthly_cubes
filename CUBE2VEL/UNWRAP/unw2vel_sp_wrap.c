/* unw2vel_wrapc.c */

/* This function is the wrapper routine called by IDL */

#include <stdio.h>

int unw2vel_sp(int argc, void *argv[])

{
    extern void unw2vel_sp1_(); /* Declare The Fortran Routine */


    float *v, *inc, *velox, *veloy, *stdx, *stdy;
    float *slx, *sly, *heading;
    short int *nx, *ny, *nz, *nn;

    v = (float *) argv[0]; 
    inc = (float *) argv[1];
    slx = (float *) argv[2];
    sly = (float *) argv[3];
    heading = (float *) argv[4];
    nx = (short int *) argv[5];
    ny = (short int *) argv[6];
    nz = (short int *) argv[7];
    velox = (float *) argv[8];
    veloy = (float *) argv[9];
    stdx = (float *) argv[10];
    stdy = (float *) argv[11];
    nn = (short int *) argv[12];


    unw2vel_sp1_(v,inc,slx,sly,heading,nx,ny,nz, \
            velox,veloy,stdx,stdy, \
            nn);


    return 1;

}

