!   cube2vel.f90
!
!   This is the routine that does the work. It in principle 
!   needs to know nothing about its call from IDL.
!
SUBROUTINE cube2vel_sp1(vx0,vy0,nx,ny,nz,errx0,erry0,weightx0,weighty0, &
velox,veloy,stdx,stdy,stdx2,stdy2,errfx,errfy,w_vx,w_vy,nn,nn2)

IMPLICIT NONE

INTEGER*2 nx, ny, nz, nz2
REAL*4 vx0(nx,ny,nz), vy0(nx,ny,nz)
REAL*4 tmpx(5*nz), tmpy(5*nz), tmp(5*nz)
REAL*4 vx(nx,ny,5*nz), vy(nx,ny,5*nz)
REAL*4 errx0(nz), erry0(nz)
REAL*4 weightx0(nz), weighty0(nz)
REAL*4 velox(nx,ny), veloy(nx,ny), med
REAL*4 stdx(nx,ny), stdy(nx,ny)
REAL*4 stdx2(nx,ny), stdy2(nx,ny)
REAL*4 w_vx(nx,ny), w_vy(nx,ny)
REAL*4 errfx(nx,ny), errfy(nx,ny)
REAL*4 errx(5*nz), erry(5*nz), weightx(5*nz)
REAL*4 weighty(5*nz)
INTEGER*2 nn(nx,ny)
INTEGER*2 nn2(nx,ny)
INTEGER*2 w(nz)
INTEGER*2 i,j,k,k0
INTEGER*4 cnt
LOGICAL mask(5*nz)
!print *,'toto'

w_vx(1:nx,1:ny)=0.
w_vy(1:nx,1:ny)=0.
velox(1:nx,1:ny)=0
veloy(1:nx,1:ny)=0
stdx(1:nx,1:ny)=0
stdy(1:nx,1:ny)=0
stdx2(1:nx,1:ny)=0
stdy2(1:nx,1:ny)=0

vx(1:nx,1:ny,1:5*nz)=0.
vy(1:nx,1:ny,1:5*nz)=0
vx(1:nx,1:ny,1:nz) = vx0
vy(1:nx,1:ny,1:nz) = vy0

!print *,nx,ny,nz

DO i = 2, nx-1
    DO j = 2, ny-1
        DO k = 1,nz 
            IF ((vx0(i-1,j,k) .NE. 0) .AND. (vx0(i+1,j,k) .NE. 0)) THEN
                vx(i,j,nz+k) = (vx0(i-1,j,k) + vx0(i+1,j,k)) / 2
                vy(i,j,nz+k) = (vy0(i-1,j,k) + vy0(i+1,j,k)) / 2
            ENDIF
        ENDDO

        DO k = 1,nz
            IF ((vx0(i,j-1,k) .NE. 0) .AND. (vx0(i,j+1,k) .NE. 0)) THEN
                vx(i,j,2*nz+k) = (vx0(i,j-1,k) + vx0(i,j+1,k)) / 2
                vy(i,j,2*nz+k) = (vy0(i,j-1,k) + vy0(i,j+1,k)) / 2
            ENDIF
        ENDDO

        DO k = 1,nz
            IF ((vx0(i-1,j-1,k) .NE. 0) .AND. (vx0(i+1,j+1,k) .NE. 0)) THEN
                vx(i,j,3*nz+k) = (vx0(i-1,j-1,k) + vx0(i+1,j+1,k)) / 2
                vy(i,j,3*nz+k) = (vy0(i-1,j-1,k) + vy0(i+1,j+1,k)) / 2
            ENDIF
        ENDDO

        DO k = 1,nz
            IF ((vx0(i+1,j-1,k) .NE. 0) .AND. (vx0(i-1,j+1,k) .NE. 0)) THEN
                vx(i,j,4*nz+k) = (vx0(i+1,j-1,k) + vx0(i-1,j+1,k)) / 2
                vy(i,j,4*nz+k) = (vy0(i+1,j-1,k) + vy0(i-1,j+1,k)) / 2
            ENDIF
        ENDDO
    ENDDO
ENDDO

errx(1:nz) = errx0
errx(nz+1:2*nz) = errx0
errx(2*nz+1:3*nz) = errx0
errx(3*nz+1:4*nz) = errx0
errx(4*nz+1:5*nz) = errx0

erry(1:nz) = erry0
erry(nz+1:2*nz) = erry0
erry(2*nz+1:3*nz) = erry0
erry(3*nz+1:4*nz) = erry0
erry(4*nz+1:5*nz) = erry0

weightx(1:nz) = weightx0
weightx(nz+1:2*nz) = weightx0/4.
weightx(2*nz+1:3*nz) = weightx0/4.
weightx(3*nz+1:4*nz) = weightx0/9.
weightx(4*nz+1:5*nz) = weightx0/9.

weighty(1:nz) = weighty0
weighty(nz+1:2*nz) = weighty0/4.
weighty(2*nz+1:3*nz) = weighty0/4.
weighty(3*nz+1:4*nz) = weighty0/9.
weighty(4*nz+1:5*nz) = weighty0/9.
!
nz2 = 5*nz
! COMPUTE MEDIAN on third DIMENSION
DO i = 1, nx
	DO j = 1, ny

    ! COMPUTE MEDIAN FOR VX
	tmpx(1:nz2) = vx(i,j,1:nz2)
    tmpy(1:nz2) = vy(i,j,1:nz2)
	cnt=COUNT( (tmpx .ne. 0) .and. (tmpy .ne. 0))
   
    ! COUNT ORIGINAL vx0 \= 0 
    nn(i,j)=COUNT(vx0(i,j,1:nz) .ne. 0)

  ! WEIGHTED STD X
    
    errfx(i,j) = sqrt( SUM( (errx0*weightx0)**2 , mask = ( vx0(i,j,1:nz) .ne. 0 ))/&
    SUM( weightx0 , mask = ( vx0(i,j,1:nz) .ne. 0 ) )**2 ) 

    errfy(i,j) = sqrt(  SUM( (erry0*weighty0)**2, mask = ( vx0(i,j,1:nz) .ne. 0 ))/&
    SUM( weighty0 , mask = ( vx0(i,j,1:nz) .ne. 0 ) )**2 )
       
    IF (cnt .GT. 1) THEN
	    CALL MDIAN( PACK( tmpx, mask= (tmpx .ne. 0) .and. (tmpy .ne. 0)),cnt,med)
        velox(i,j)=med
	ELSE IF (cnt .EQ. 1) THEN
	    velox(i,j)=SUM( tmpx, mask= (tmpx .ne. 0) .and. (tmpy .ne. 0))
	ELSE
	    velox(i,j)=0 
	ENDIF
    IF (cnt .GT. 3) THEN    
        ! COMPUTE STD FOR VX
        med = sum( tmpx , mask =  (tmpx .ne. 0) .and. (tmpy .ne. 0)) / cnt
        stdx(i,j) = sqrt( sum( (tmpx - med)**2 , mask = (tmpx .ne. 0) .and. (tmpy .ne. 0)) / (cnt-1))
    ENDIF
    
    ! COMPUTE MEDIAN FOR VY
	IF (cnt .GT. 1) THEN
        CALL MDIAN( PACK( tmpy, mask= (tmpx .ne. 0) .and. (tmpy .ne. 0)),cnt,med)
        veloy(i,j)=med
	ELSE IF (CNT .EQ. 1) THEN
	    veloy(i,j)=SUM(tmpy, mask= (tmpx .ne. 0) .and. (tmpy .ne. 0))
	ELSE
	    veloy(i,j)=0
	ENDIF
    IF (cnt .GT. 3) THEN
        ! COMPUTE STD FOR VY
        med = sum( tmpy , mask =  (tmpx .ne. 0) .and. (tmpy .ne. 0)) / cnt
        stdy(i,j) = sqrt( sum( (tmpy - med)**2 , mask = (tmpx .ne. 0) .and. (tmpy .ne. 0)) / (cnt-1))

        ! REMOVE VALUES MORE THAN 2 STD AWAY FOR MEDIAN
	    DO k = 1, nz2
	        IF (vx(i,j,k) .NE. 0) THEN
	            IF ( ((vx(i,j,k)-velox(i,j)) .GT. (2*stdx(i,j))) .OR. &
                ((vy(i,j,k)-veloy(i,j)) .GT. (2*stdy(i,j))) ) THEN
	                vx(i,j,k)=0
	                vy(i,j,k)=0
                ENDIF
	        ENDIF
	    ENDDO
    ENDIF	

    
    tmpx = vx(i,j,1:nz2)
    tmpy = vy(i,j,1:nz2)
    cnt=COUNT( (tmpx .ne. 0) .and. (tmpy .ne. 0))
    mask=( (tmpx .ne. 0) .and. (tmpy .ne. 0) .and. (weightx .ne. 0) .and. (weighty .ne. 0) )
    nn2(i,j)=cnt
    ! WEIGHTED MEAN X
    IF (cnt .NE. 0) THEN
        
        w_vx(i,j) = sum( tmpx*weightx , mask = mask ) / sum( weightx, mask = mask )

        ! WEIGHTED STD X
        !errfx(i,j) = sqrt( sum( (errx*weightx)**2, mask = ( (tmpx .ne. 0) .and. (tmpy .ne. 0))) / &
        !sum( weightx, mask = ( (tmpx .ne. 0) .and. (tmpy .ne. 0)))**2 )

        IF (cnt .GE. 2) THEN
            stdx2(i,j) = sqrt( cnt * sum( weightx*(tmpx - w_vx(i,j))**2 , mask = mask ) / &
            sum( weightx , mask = mask ) / (cnt-1) )
        ENDIF


        ! WEIGHTED MEAN Y
        w_vy(i,j) = sum( tmpy*weighty, mask = mask ) / sum( weighty, mask = mask )

        IF (cnt .GE. 2) THEN
            stdy2 (i,j) = sqrt( cnt * sum( weighty * &
                (tmpy - w_vy(i,j))**2 , mask = mask ) / &
                sum( weighty , mask = mask ) / (cnt-1) )
        ENDIF
    ENDIF
	ENDDO
ENDDO

END SUBROUTINE cube2vel_sp1


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

