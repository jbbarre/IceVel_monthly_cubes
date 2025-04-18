/* cube2vel_wrapc.c */

/* This function is the wrapper routine called by IDL */

#include <stdio.h>

int cube2vel_sp(int argc, void *argv[])

{
    extern void cube2vel_sp1_(); /* Declare The Fortran Routine */


    float *vx, *vy, *velox, *veloy, *stdx, *stdy;
    float *errx, *erry, *stdx2, *stdy2, *errfx, *errfy;
    float *w_vx, *w_vy, *weightx, *weighty;
    short int *nx, *ny, *nz, *nn, *nn2;

    vx = (float *) argv[0]; 
    vy = (float *) argv[1];
    nx = (short int *) argv[2];
    ny = (short int *) argv[3];
    nz = (short int *) argv[4];
    errx = (float *) argv[5];
    erry = (float *) argv[6];
    weightx = (float *) argv[7];
    weighty = (float *) argv[8];
    velox = (float *) argv[9];
    veloy = (float *) argv[10];
    stdx = (float *) argv[11];
    stdy = (float *) argv[12];
    stdx2 = (float *) argv[13];
    stdy2 = (float *) argv[14];
    errfx = (float *) argv[15];
    errfy = (float *) argv[16];
    w_vx = (float *) argv[17];
    w_vy = (float *) argv[18];
    nn = (short int *) argv[19];
    nn2 = (short int *) argv[20];

    cube2vel_sp1_(vx,vy,nx,ny,nz,errx,erry,weightx,weighty, \
            velox,veloy,stdx,stdy,stdx2, \
            stdy2,errfx,errfy,w_vx,w_vy,nn,nn2);


    return 1;

}

