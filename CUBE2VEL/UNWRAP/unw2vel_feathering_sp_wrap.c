/* unw2vel_wrapc.c */

/* This function is the wrapper routine called by IDL */

#include <stdio.h>

int unw2vel_feathering_sp(int argc, void *argv[])

{
    extern void unw2vel_feathering_sp1_(); /* Declare The Fortran Routine */


    float *v, *inc, *velox, *veloy, *stdx, *stdy, *noise;
    float *errvx, *errvy;
    float *slx, *sly, *heading, *feathering, *angle;
    short int *nx, *ny, *nz, *nn;

    v = (float *) argv[0]; 
    inc = (float *) argv[1];
    slx = (float *) argv[2];
    sly = (float *) argv[3];
    heading = (float *) argv[4];
    feathering = (float *) argv[5];
    noise = (float *) argv[6];
    nx = (short int *) argv[7];
    ny = (short int *) argv[8];
    nz = (short int *) argv[9];
    velox = (float *) argv[10];
    veloy = (float *) argv[11];
    stdx = (float *) argv[12];
    stdy = (float *) argv[13];
    errvx = (float *) argv[14];
    errvy = (float *) argv[15];
    nn = (short int *) argv[16];
    angle = (float *) argv[17];

    unw2vel_feathering_sp1_(v,inc,slx,sly,heading,feathering,noise,nx,ny,nz, \
            velox,veloy,stdx,stdy,errvx,errvy, \
            nn,angle);


    return 1;

}

