!   cube2vel.f90
!
!   This is the routine that does the work. It in principle 
!   needs to know nothing about its call from IDL.
!
SUBROUTINE unw2vel_feathering_sp1(v,inc,slx,sly,heading,feathering,noise,nx,ny,nz, &
velox,veloy,stdx,stdy,errvx,errvy,nn,angle)

IMPLICIT NONE

INTEGER*2 nx, ny, nz
REAL*4 v(nx,ny,nz), inc(nx,ny,nz), heading(nx,ny,nz)
REAL*4 feathering(nx,ny,nz)
REAL*4 noise(nx,ny,nz)
REAL*4 velox(nx,ny), veloy(nx,ny), med
REAL*4 stdx(nx,ny), stdy(nx,ny),slx(nx,ny),sly(nx,ny)
REAL*4 errvx(nx,ny), errvy(nx,ny)
REAL*4 alpha,beta,dtor,vra,vrd,tetadd,tetaaa,pi
REAL*4 slopex,slopey,sinalphak_square
REAL*4 DET,angle,vtmpx,vtmpy,medx,medy
!REAL*4 dvra, dvrd
REAL*4 UNITY(2,2),A(2,2),B(2,2)
REAL*4 AB(2,2),ABC(2,2),C(2,2),D(2,2)
REAL*4 INVABC(2,2)
REAL*4 vx(10000),vy(10000),WEIGHT(10000)
REAL*4 SIGMA_VX(10000),SIGMA_VY(10000)
INTEGER*2 nn(nx,ny)
INTEGER*2 i,j,k,l,m,k0
INTEGER*4 cnt,num

UNITY(1,1)=1.
UNITY(1,2)=0.
UNITY(2,1)=0.
UNITY(2,2)=1.

VX(:)=0.
VY(:)=0.
WEIGHT(:)=0.
SIGMA_VX(:)=0.
SIGMA_VY(:)=0.

dtor = 3.1415927/180.
pi = 3.1415927
! COMPUTE MEDIAN on third DIMENSION
DO i = 1, nx
	
    DO j = 1, ny
        NUM=1
        vx(:)=0.
        vy(:)=0.

        DO k = 1, nz-1
            IF ((v(i,j,k) .NE. 0) .AND. (NUM .LE. 10000)) THEN

                DO l = k+1,nz 
                    IF ((v(i,j,l) .NE. 0) .AND. (NUM .LE. 10000)) THEN ! {{{
                        
                        alpha = heading(i,j,l)*dtor
                        beta = heading(i,j,k)*dtor

                        IF (modulo(ABS(min(360.*dtor - abs(alpha - beta),abs(alpha -beta))),pi) .GE. angle*dtor) THEN ! {{{
                            
                            alpha = alpha - beta
                            vra = v(i,j,k)
                            vrd = v(i,j,l)

                            !dvra = 1.
                            !dvrd = 1.

                            tetadd = inc(i,j,k)
                            tetaaa = inc(i,j,l)

                            slopex = slx(i,j)
                            slopey = sly(i,j)

                            A(1,1) = cos(beta)
                            A(2,1) = cos(alpha+beta)
                            A(1,2) = sin(beta)
                            A(2,2) = sin(alpha+beta)

                            sinalphak_square = max(sin(alpha)*sin(alpha),0.000001)
                            
                            B(1,1) = 1. / sinalphak_square
                            B(2,1) = -cos(alpha) / sinalphak_square
                            B(1,2) = -cos(alpha) / sinalphak_square
                            B(2,2) = 1. / sinalphak_square

                            AB = MATMUL(B,A)
                            !AB(1,1) = B(1,1)*A(1,1) + B(2,1)*A(1,2)
                            !AB(2,1) = B(1,1)*A(2,1) + B(2,1)*A(2,2)
                            !AB(1,2) = B(1,2)*A(1,1) + B(2,2)*A(1,2)
                            !AB(2,2) = B(1,2)*A(2,1) + B(2,2)*A(2,2)

                            C(1,1) = slopex/tan(tetaaa)
                            C(2,1) = slopey/tan(tetaaa)
                            C(1,2) = slopex/tan(tetadd)
                            C(2,2) = slopey/tan(tetadd)

                            ABC = MATMUL(C,AB)

                            ABC = UNITY - ABC
                            !ABC(1,1) = UNITY(1,1) - C(1,1)*AB(1,1) + C(2,1)*AB(1,2)
                            !ABC(2,1) = UNITY(2,1) - C(1,1)*AB(2,1) + C(2,1)*AB(2,2)
                            !ABC(1,2) = UNITY(1,2) - C(1,2)*AB(1,1) + C(2,2)*AB(1,2)
                            !ABC(2,2) = UNITY(2,2) - C(1,2)*AB(2,1) + C(2,2)*AB(2,2)
                            
                            DET = (ABC(1,1)*ABC(2,2)-ABC(2,1)*ABC(1,2))
                            INVABC(1,1) = ABC(2,2) / DET
                            INVABC(2,1) = -ABC(2,1) / DET
                            INVABC(1,2) = -ABC(1,2) / DET
                            INVABC(2,2) = ABC(1,1) / DET
                           
                            D = MATMUL(INVABC,AB)
                            !D(1,1) = AB(1,1)*INVABC(1,1) + AB(2,1)*INVABC(1,2)
                            !D(2,1) = AB(1,1)*INVABC(2,1) + AB(2,1)*INVABC(2,2)
                            !D(1,2) = AB(1,2)*INVABC(1,1) + AB(2,2)*INVABC(1,2)
                            !D(2,2) = AB(1,2)*INVABC(2,1) + AB(2,2)*INVABC(2,2)
                            
                            VX(num) = vra*D(1,1)+vrd*D(2,1)
                            VY(num) = vra*D(1,2)+vrd*D(2,2)
                            SIGMA_VX(num) = abs(noise(i,j,k)*D(1,1))+ &
                                abs(noise(i,j,l)*D(2,1))
                            SIGMA_VY(num) = abs(noise(i,j,k)*D(1,2))+ &
                                abs(noise(i,j,l)*D(2,2))

                            WEIGHT(num) = feathering(i,j,k)*feathering(i,j,l)

                                    !*&
                            !modulo(ABS(min(360.*dtor - abs(alpha - beta),abs(alpha -beta))),pi)*&
                            !modulo(ABS(min(360.*dtor - abs(alpha - beta),abs(alpha -beta))),pi)
                            !ABS(min(360.*dtor - abs(alpha - beta),abs(alpha -beta)))

                            num = num+1
                        ENDIF ! }}}
                    ENDIF ! }}}
                ENDDO
            ENDIF
        ENDDO ! k

        cnt=COUNT( ((vx .ne. 0) .and. (vy .ne. 0)))

        nn(i,j) = COUNT((v(i,j,:) .ne. 0))
            
        if (cnt .eq. 1) THEN
           velox(i,j) = sum(vx)
           veloy(i,j) = sum(vy)
           errvx(i,j) = sum(SIGMA_VX)
           errvy(i,j) = sum(SIGMA_VY)
        ENDIF

        if (cnt .ge. 2) THEN
        velox(i,j) = sum( VX*WEIGHT , mask =  ((vx .ne. 0) .and. (vy .ne. 0)))/&
                sum(WEIGHT , mask =  ((vx .ne. 0) .and. (vy .ne. 0)))
   
        veloy(i,j) = sum( VY*WEIGHT , mask =  ((vx .ne. 0) .and. (vy .ne. 0)))/&
                sum(WEIGHT , mask =  ((vx .ne. 0) .and. (vy .ne. 0)))

        errvx(i,j) = sqrt( sum( (SIGMA_VX*WEIGHT)**2 , mask =  ((vx .ne. 0) .and. (vy .ne.0)))/&
                sum( (WEIGHT)**2, mask =  ((vx .ne. 0) .and. (vy .ne.0))) )
        errvy(i,j) = sqrt( sum( (SIGMA_VY*WEIGHT)**2 , mask =  ((vx .ne. 0) .and. (vy .ne. 0)))/&
                        sum( (WEIGHT)**2, mask =  ((vx .ne. 0) .and. (vy .ne. 0))) )

        ENDIF

        if (cnt .ge. 3) THEN
            stdx(i,j) = sqrt( sum( (vx - velox(i,j))**2 , mask = ((vx .ne. 0) .and. (vy .ne. 0))) / (cnt-1))
            stdy(i,j) = sqrt( sum( (vy - veloy(i,j))**2 , mask = ((vx .ne. 0) .and. (vy .ne. 0))) / (cnt-1))

            IF (cnt .GT. 3) THEN 
                !medx=velox(i,j)
                !medy=veloy(i,j)
                !CALL MDIAN( PACK( vx, mask= (vx .ne. 0) .and. (vy .ne. 0)),cnt,medx)
                !CALL MDIAN( PACK( vy, mask= (vx .ne. 0) .and. (vy .ne. 0)),cnt,medy)

                !vtmpx = 0.1*abs(medx)
                !vtmpx = 3.*stdx(i,j)
                !IF (vtmpx .LT. 3) THEN
                !    vtmpx=3
                !ENDIF

                !vtmpy = 0.1*abs(medy)
                !vtmpy = 3.*stdy(i,j)
                !IF (vtmpy .LT. 3) THEN
                !    vtmpy=3
                !ENDIF

                !DO k = 1, num
                !IF (vx(k) .NE. 0) THEN
                !    IF ((ABS(vx(k)-medx) .GT. vtmpx ) .OR. &
                !    (ABS(vy(k)-medy) .GT. vtmpy )) THEN
                !        vx(k)=0
                !        vy(k)=0
                !        weight(k)=0
                !    ENDIF
                !ENDIF
                !ENDDO
            ENDIF
        ENDIF

        !cnt=COUNT( ((vx .ne. 0) .and. (vy .ne. 0)))
        !nn(i,j) = cnt

        !if (cnt .ge. 1) THEN
        !    velox(i,j) = sum( VX*WEIGHT , mask =  ((vx .ne. 0).and. (vy .ne. 0)))/&
        !        sum(WEIGHT , mask = ((vx .ne. 0) .and. (vy .ne. 0)))

        !    veloy(i,j) = sum( VY*WEIGHT, mask =  ((vx .ne. 0) .and. (vy .ne. 0)))/&
        !        sum(WEIGHT , mask =  ((vx .ne. 0) .and. (vy .ne. 0)))

        !    errvx(i,j) = sqrt( sum( (SIGMA_VX*WEIGHT)**2 , mask =  ((vx .ne. 0) .and. (vy .ne.0)))/&
        !        sum( (WEIGHT)**2, mask =  ((vx .ne. 0) .and. (vy .ne.0))) )
        !    errvy(i,j) = sqrt( sum( (SIGMA_VY*WEIGHT)**2 , mask =  ((vx .ne. 0) .and. (vy .ne. 0)))/&
        !        sum( (WEIGHT)**2, mask =  ((vx .ne. 0) .and. (vy .ne. 0))) )
        !ENDIF

        !IF (cnt .le. 2) THEN
        !    velox(i,j) = sum( VX , mask =  ((vx .ne. 0) .and. (vy .ne. 0))) / cnt
        !    veloy(i,j) = sum( VY , mask =  ((vx .ne. 0) .and. (vy .ne. 0))) / cnt
        !ENDIF

        !IF (cnt .GE. 3) THEN

        !   CALL MDIAN( PACK( VX , mask = ((vx .ne. 0) .and. (vy .ne. 0))),cnt,med)
        !    velox(i,j)=med
        !    stdx(i,j) = sqrt( sum( (vx - med)**2 , mask = ((vx .ne. 0) .and. (vy .ne. 0))) / (cnt-1))


        !    CALL MDIAN( PACK( VY , mask = ((vx .ne. 0) .and. (vy .ne. 0))),cnt,med)
        !    veloy(i,j)=med
        !    stdy(i,j) = sqrt( sum( (vy - med)**2 , mask = ((vx .ne. 0) .and. (vy .ne. 0))) / (cnt-1))

        !ENDIF

    ENDDO ! j
ENDDO ! i

END SUBROUTINE unw2vel_feathering_sp1


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

