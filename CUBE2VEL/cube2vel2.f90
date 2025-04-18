!   cube2vel.f90
!
!   This is the routine that does the work. It in principle 
!   needs to know nothing about its call from IDL.
!
SUBROUTINE cube2vel1(vx,vy,nx,ny,nz,errx,erry, &
velox,veloy,stdx,stdy,stdx2,stdy2,errfx,errfy,w_vx,w_vy,nn,nn2)

IMPLICIT NONE

INTEGER*2 nx, ny, nz
REAL*4 vx(nx,ny,nz), vy(nx,ny,nz), tmp(nz), tmp2(nz)
REAL*4 velox(nx,ny), veloy(nx,ny), med
REAL*4 stdx(nx,ny), stdy(nx,ny)
REAL*4 stdx2(nx,ny), stdy2(nx,ny)
REAL*4 w_vx(nx,ny), w_vy(nx,ny)
REAL*4 errfx(nx,ny), errfy(nx,ny)
REAL*4 errx(nz), erry(nz), weight(nz)
INTEGER*2 nn(nx,ny)
INTEGER*2 nn2(nx,ny)
INTEGER*2 w(nz)
INTEGER*2 i,j,k
INTEGER*4 cnt
!print *,'toto'


! COMPUTE MEDIAN on third DIMENSION
DO i = 1, nx
	DO j = 1, ny

    ! COMPUTE MEDIAN FOR VX
	tmp = vx(i,j,1:nz)
	cnt=COUNT(tmp .ne. 0)
    nn(i,j)=cnt
	IF (cnt .GT. 1) THEN
	    CALL MEDIAN( PACK( tmp, mask= tmp .ne. 0),cnt,0,med)
	    velox(i,j)=med
	ELSE IF (COUNT(tmp .ne. 0) .EQ. 1) THEN
	    velox(i,j)=SUM( tmp, mask= tmp .ne. 0)
	ELSE
	    velox(i,j)=0 
	ENDIF
	
    ! COMPUTE STD FOR VX
	med = sum( tmp , mask =  tmp .ne. 0) / cnt
	stdx(i,j) = sqrt( sum( (tmp - med)**2 , mask = tmp .ne. 0) / (cnt-1))
	
    ! COMPUTE MEDIAN FOR VY
	tmp = vy(i,j,1:nz) 
	cnt=COUNT(tmp .ne. 0)
	IF (cnt .GT. 1) THEN
	    CALL MEDIAN( PACK( tmp, mask= tmp .ne. 0),cnt,0,med)
	    veloy(i,j)=med
	ELSE IF (COUNT(tmp .ne. 0) .EQ. 1) THEN
	    veloy(i,j)=SUM(tmp, mask= tmp .ne. 0)
	ELSE
	    veloy(i,j)=0
	ENDIF
	
    ! COMPUTE STD FOR VY
	med = sum( tmp , mask =  tmp .ne. 0) / cnt
	stdy(i,j) = sqrt( sum( (tmp - med)**2 , mask = tmp .ne. 0) / (cnt-1))
	
    ! REMOVE VALUES MORE THAN 2 STD AWAY FOR MEDIAN
	IF (CNT .GE. 2) THEN
	    DO k = 1, nz
	        IF (vx(i,j,k) .NE. 0) THEN
	            IF ( ((vx(i,j,k)-velox(i,j)) .GT. (2*stdx(i,j))) .OR. &
                ((vy(i,j,k)-veloy(i,j)) .GT. (2*stdy(i,j))) ) THEN
	                vx(i,j,k)=0
	                vy(i,j,k)=0
                ENDIF
	        ENDIF
	    ENDDO
    ENDIF	


    tmp = vx(i,j,1:nz)
    cnt=COUNT(tmp .ne. 0)
    weight = 1 / errx**2
    nn2(i,j)=cnt
    ! WEIGHTED MEAN X
    w_vx(i,j) = sum( tmp*weight , mask = tmp .ne. 0) / sum( weight, mask = tmp.ne. 0)

    ! WEIGHTED STD X
    errfx(i,j) = sqrt( sum( (errx*weight)**2, mask = tmp .ne. 0) / sum( weight, mask = tmp.ne. 0)**2 )

    stdx2(i,j) = sqrt( cnt * sum( weight * &
        (tmp - w_vx(i,j))**2 , mask = tmp .ne. 0 ) / &
        sum( weight , mask = tmp .ne. 0) / (cnt-1) )

    tmp = vy(i,j,1:nz)
    cnt=COUNT(tmp .ne. 0)
    weight = 1 / erry**2

    ! WEIGHTED MEAN Y
    w_vy(i,j) = sum( tmp*weight , mask = tmp .ne. 0) / sum( weight, mask = tmp.ne. 0)
    errfy(i,j) = sqrt( sum( (erry*weight)**2, mask = tmp .ne. 0) / sum( weight, mask = tmp.ne. 0)**2 )

    stdy2(i,j) = sqrt( cnt * sum( weight * &
        (tmp - w_vy(i,j))**2 , mask = tmp .ne. 0 ) / &
        sum( weight , mask = tmp .ne. 0) / (cnt-1) )

	ENDDO
ENDDO



END


