!   cube2vel.f90
!
!   This is the routine that does the work. It in principle 
!   needs to know nothing about its call from IDL.
!
!   Gray et al. 2011 for 3d vel mapping with unwrapped phases
!
SUBROUTINE unw2vel_3d_sp1(v,inc,slx,sly,heading,nx,ny,nz, &
velox,veloy,veloz,stdx,stdy,stdz,nn)

IMPLICIT NONE

INTEGER*2 nx, ny, nz
REAL*4 v(nx,ny,nz), inc(nx,ny,nz), heading(nx,ny,nz)
REAL*4 velox(nx,ny), veloy(nx,ny), veloz(nx,ny), med
REAL*4 stdx(nx,ny), stdy(nx,ny),slx(nx,ny),sly(nx,ny)
REAL*4 stdz(nx,ny)
REAL*4 Xa,Xb,Xc,dtor,vra,vrb,vrc
REAL*4 THETAa,THETAb,THETAc
REAL*4 slopex,slopey,sinalphak_square
REAL*4 invdetT,invT1
REAL*4 A,B,C,D,E
REAL*4 F,G,H,KK ! KK because K is used for loop
REAL*4 vx(10000),vy(10000),vz(10000)
INTEGER*2 nn(nx,ny)
INTEGER*2 i,j,k,l,m,k0
INTEGER*4 cnt,num

dtor = 3.1415927/180.
! COMPUTE MEDIAN on third DIMENSION
DO i = 1, nx
	
    DO j = 1, ny

        NUM=1
        vx(:)=0.
        vy(:)=0.
        vz(:)=0.

        !print *,i,j
        DO k = 1, nz-2
            IF ((v(i,j,k) .NE. 0) .AND. (NUM .LE. 10000)) THEN
                Xa = 6.2831853 - heading(i,j,k)*dtor
                THETAa = inc(i,j,k)
               
                vra = v(i,j,k)*sin(THETAa) !! ground range -> slant range

                DO l = k+1,nz-1 

                    IF ((v(i,j,l) .NE. 0) .AND. (NUM .LE. 10000)) THEN ! {{{
                        
                        Xb = 6.2831853 - heading(i,j,l)*dtor        
                        THETAb = inc(i,j,l)
                       
                        vrb = v(i,j,l)*sin(THETAb) !! ground range -> slant range
                        
                        IF (ABS(min(6.2831853 - abs(Xa - Xb),abs(Xa - Xb))) .GE. 20.*dtor) THEN ! {{{

                            DO m = l+1,nz
                                IF ((v(i,j,m) .NE. 0) .AND. (NUM .LE. 10000)) THEN

                                    Xc = 6.2831853 - heading(i,j,m)*dtor
                                        
                                    IF ((ABS(min(6.2831853 - abs(Xa - Xc),abs(Xa - Xc))) .GE. 20.*dtor) .AND. &
                                        (ABS(min(6.2831853 - abs(Xb - Xc),abs(Xb -Xc))) .GE. 20.*dtor)) THEN
                                        THETAc=inc(i,j,m)
                                        vrc=v(i,j,m)*sin(THETAc) !! ground range -> slant range

                                        invT1 = 1. / (sin(Xa-Xc)/tan(THETAb) + &
                                                sin(Xc-Xb)/tan(THETAa) + &
                                                sin(Xb-Xa)/tan(THETAc))
                                        
                                       ! vz(num) = invT1 * (sin(Xc-Xa)/sin(THETAb)*vrb + &
                                       !         sin(Xa-Xb)/sin(THETAc)*vrc + &
                                       !         sin(Xb-Xc)/sin(THETAa)*vra)

                                        invdetT = 1./(sin(THETAa)*sin(THETAb)*sin(THETAc)) * invT1

                                        A =  sin(Xb)*sin(THETAb)*cos(THETAc) - cos(THETAb)*sin(Xc)*sin(THETAc)
                                        B = -cos(THETAb)*cos(Xc)*sin(THETAc) + cos(Xb)*sin(THETAb)*cos(THETAc)
                                        C =  sin(THETAb)*sin(THETAc)*sin(Xb-Xc)
                                        D =  cos(THETAa)*sin(Xc)*sin(THETAc) - sin(Xa)*sin(THETAa)*cos(THETAc)
                                        E = -cos(Xa)*sin(THETAa)*cos(THETAc) + cos(THETAa)*cos(Xc)*sin(THETAc)
                                        F =  sin(THETAa)*sin(THETAc)*sin(Xc-Xa)
                                        G =  sin(Xa)*sin(THETAa)*cos(THETAb) - cos(THETAa)*sin(Xb)*sin(THETAb)
                                        H = -cos(THETAa)*cos(Xb)*sin(THETAb) + cos(Xa)*sin(THETAa)*cos(THETAb) 
                                        !! K => H error in Gray 2011 sup. materials
                                        KK = sin(THETAa)*sin(THETAb)*sin(Xa-Xb)
                                        !! KK -> K because k is used for loop

                                        !! vra,b,c is displacement in ground range.
                                        !! Gray 2011 uses slant range (multiply
                                        !           by sin(inc).

                                        vx(num) = invdetT * (A*vra + D*vrb + G*vrc)
                                        vy(num) = invdetT * (B*vra + E*vrb + H*vrc)
                                        vz(num) = invdetT * (C*vra + F*vrb + KK*vrc)
                                        num = num+1
                                    ENDIF
                                ENDIF
                            ENDDO
                        ENDIF ! }}}

                    ENDIF ! }}}
                ENDDO
            ENDIF
        ENDDO ! k

        cnt=COUNT((vx .ne. 0) .and. (vy .ne. 0) .and. (vy .ne. 0))
        nn(i,j) = cnt

        IF (cnt .le. 2) THEN
            velox(i,j) = sum( VX , mask =  ((vx .ne. 0) .and. &
            (vy .ne. 0) .and. (vy .ne. 0)) ) / cnt
            veloy(i,j) = sum( VY , mask =  ((vx .ne. 0) .and. &
            (vy .ne. 0) .and. (vy .ne. 0))) / cnt
            veloz(i,j) = sum( VZ , mask =  ((vx .ne. 0) .and. &
            (vy .ne. 0) .and. (vy .ne. 0))) / cnt
        ENDIF
        IF (cnt .GE. 3) THEN

            CALL MDIAN( PACK( VX , mask = ((vx .ne. 0) .and. (vy .ne. 0) .and. &
            (vy .ne. 0))),cnt,med)
            velox(i,j)=med
            stdx(i,j) = sqrt( sum( (vx - med)**2 , mask = ((vx .ne. 0) .and. &
            (vy .ne. 0) .and. (vy .ne. 0))) / (cnt-1))

            CALL MDIAN( PACK( VY , mask = ((vx .ne. 0) .and. (vy .ne. 0) .and. &
            (vy .ne. 0))),cnt,med)
            veloy(i,j)=med
            stdy(i,j) = sqrt( sum( (vy - med)**2 , mask = ((vx .ne. 0) .and. & 
            (vy .ne. 0) .and. (vy .ne. 0))) / (cnt-1))

            CALL MDIAN( PACK( VZ , mask = ((vx .ne. 0) .and. (vy .ne. 0) .and. &
            (vy .ne. 0))),cnt,med)
            veloz(i,j)=med
            stdz(i,j) = sqrt( sum( (vz - med)**2 , mask = ((vx .ne. 0) .and. &
            (vy .ne. 0) .and. (vy .ne. 0))) / (cnt-1))

        ENDIF
        



    ENDDO ! j
ENDDO ! i

END SUBROUTINE unw2vel_3d_sp1


SUBROUTINE MDIAN(X,N,XMED)
real X(N)
call hpsort(N,X)
N2=N/2
if (2*N2.eq.N) then
    XMED = 0.5*(X(N2)+X(N2+1))
else
    XMED = X(N2+1)
endif
return
END

SUBROUTINE HPSORT(N,RA)
  real RA(N)
  L=N/2+1
  IR=N
  !The index L will be decremented from its initial value during the
  !"hiring" (heap creation) phase. Once it reaches 1, the index IR 
  !will be decremented from its initial value down to 1 during the
  !"retirement-and-promotion" (heap selection) phase.
10 continue
  if(L > 1)then
    L=L-1
    RRA=RA(L)
  else
    RRA=RA(IR)
    RA(IR)=RA(1)
    IR=IR-1
    if(IR.eq.1)then
      RA(1)=RRA
      return
    end if
  end if
  I=L
  J=L+L
20 if(J.le.IR)then
  if(J < IR)then
    if(RA(J) < RA(J+1))  J=J+1
  end if
  if(RRA < RA(J))then
    RA(I)=RA(J)
    I=J; J=J+J
  else
    J=IR+1
  end if
  goto 20
  end if
  RA(I)=RRA
  goto 10
END

