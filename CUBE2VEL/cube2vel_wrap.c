/* cube2vel_wrapc.c */

/* This function is the wrapper routine called by IDL */

#include <stdio.h>

int cube2vel(int argc, void *argv[])

{
    extern void cube2vel1_(); /* Declare The Fortran Routine */


    float *vx, *vy, *velox, *veloy, *stdx, *stdy;
    float *errx, *erry, *stdx2, *stdy2, *errfx, *errfy;
    float *w_vx, *w_vy;
    short int *nx, *ny, *nz, *nn, *nn2;

    vx = (float *) argv[0]; 
    vy = (float *) argv[1];
    nx = (short int *) argv[2];
    ny = (short int *) argv[3];
    nz = (short int *) argv[4];
    errx = (float *) argv[5];
    erry = (float *) argv[6];
    velox = (float *) argv[7];
    veloy = (float *) argv[8];
    stdx = (float *) argv[9];
    stdy = (float *) argv[10];
    stdx2 = (float *) argv[11];
    stdy2 = (float *) argv[12];
    errfx = (float *) argv[13];
    errfy = (float *) argv[14];
    w_vx = (float *) argv[15];
    w_vy = (float *) argv[16];
    nn = (short int *) argv[17];
    nn2 = (short int *) argv[18];

    cube2vel1_(vx,vy,nx,ny,nz,errx,erry, \
            velox,veloy,stdx,stdy,stdx2, \
            stdy2,errfx,errfy,w_vx,w_vy,nn,nn2);


    return 1;

}

