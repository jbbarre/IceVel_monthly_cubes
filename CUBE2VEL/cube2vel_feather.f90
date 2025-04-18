!   cube2vel.f90
!
!   This is the routine that does the work. It in principle 
!   needs to know nothing about its call from IDL.
!
SUBROUTINE cube2vel_feather_sp1(vx,vy,nx,ny,nz,errx,erry,weightx,weighty, &
velox,veloy,stdx,stdy,stdx2,stdy2,errfx,errfy,w_vx,w_vy,nn,nn2)

IMPLICIT NONE

INTEGER*2 nx, ny, nz
REAL*4 vx(nx,ny,nz), vy(nx,ny,nz)
REAL*4 errx(nz), erry(nz)
REAL*4 weightx(nx,ny,nz), weighty(nx,ny,nz)
REAL*4 velox(nx,ny), veloy(nx,ny), med
REAL*4 stdx(nx,ny), stdy(nx,ny)
REAL*4 stdx2(nx,ny), stdy2(nx,ny)
REAL*4 w_vx(nx,ny), w_vy(nx,ny)
REAL*4 errfx(nx,ny), errfy(nx,ny)
INTEGER*2 nn(nx,ny)
INTEGER*2 nn2(nx,ny)
INTEGER*2 w(nz)
INTEGER*2 i,j,k,k0
INTEGER*4 cnt
LOGICAL mask(nz)
!print *,'toto'

w_vx(1:nx,1:ny)=0.
w_vy(1:nx,1:ny)=0.
velox(1:nx,1:ny)=0
veloy(1:nx,1:ny)=0
stdx(1:nx,1:ny)=0
stdy(1:nx,1:ny)=0
stdx2(1:nx,1:ny)=0
stdy2(1:nx,1:ny)=0

! COMPUTE MEDIAN on third DIMENSION
DO i = 1, nx
	DO j = 1, ny

    ! COMPUTE MEDIAN FOR VX
	cnt=COUNT( (vx(i,j,1:nz) .ne. 0) .and. (vy(i,j,1:nz) .ne. 0))
   
    ! COUNT ORIGINAL vx \= 0 
    nn(i,j)=COUNT(vx(i,j,1:nz) .ne. 0)

  ! WEIGHTED STD X
    
    errfx(i,j) = sqrt( SUM( (errx*weightx(i,j,1:nz))**2 , mask = ( vx(i,j,1:nz) .ne. 0 ))/&
    SUM( weightx(i,j,1:nz) , mask = ( vx(i,j,1:nz) .ne. 0 ) )**2 ) 

    errfy(i,j) = sqrt(  SUM( (erry*weighty(i,j,1:nz))**2, mask = ( vx(i,j,1:nz) .ne. 0 ))/&
    SUM( weighty(i,j,1:nz) , mask = ( vx(i,j,1:nz) .ne. 0 ) )**2 )
       
    IF (cnt .GT. 1) THEN
	    CALL MDIAN( PACK( vx(i,j,1:nz), mask= (vx(i,j,1:nz) .ne. 0) .and. (vy(i,j,1:nz) .ne. 0)),cnt,med)
        velox(i,j)=med
	ELSE IF (cnt .EQ. 1) THEN
	    velox(i,j)=SUM( vx(i,j,1:nz), mask= (vx(i,j,1:nz) .ne. 0) .and. (vy(i,j,1:nz) .ne. 0))
	ELSE
	    velox(i,j)=0 
	ENDIF
    IF (cnt .GT. 3) THEN    
        ! COMPUTE STD FOR VX
        med = sum( vx(i,j,1:nz) , mask =  (vx(i,j,1:nz) .ne. 0) .and. (vy(i,j,1:nz) .ne. 0)) / cnt
        stdx(i,j) = sqrt( sum( (vx(i,j,1:nz) - med)**2 , mask = (vx(i,j,1:nz) .ne. 0) .and. (vy(i,j,1:nz) .ne. 0)) / (cnt-1))
    ENDIF
    
    ! COMPUTE MEDIAN FOR VY
	IF (cnt .GT. 1) THEN
        CALL MDIAN( PACK( vy(i,j,1:nz), mask= (vx(i,j,1:nz) .ne. 0) .and. (vy(i,j,1:nz) .ne. 0)),cnt,med)
        veloy(i,j)=med
	ELSE IF (CNT .EQ. 1) THEN
	    veloy(i,j)=SUM(vy(i,j,1:nz), mask= (vx(i,j,1:nz) .ne. 0) .and. (vy(i,j,1:nz) .ne. 0))
	ELSE
	    veloy(i,j)=0
	ENDIF
    IF (cnt .GT. 3) THEN
        ! COMPUTE STD FOR VY
        med = sum( vy(i,j,1:nz) , mask =  (vx(i,j,1:nz) .ne. 0) .and. (vy(i,j,1:nz) .ne. 0)) / cnt
        stdy(i,j) = sqrt( sum( (vy(i,j,1:nz) - med)**2 , mask = (vx(i,j,1:nz) .ne. 0) .and. (vy(i,j,1:nz) .ne. 0)) / (cnt-1))

        ! REMOVE VALUES MORE THAN 2 STD AWAY FOR MEDIAN
	    DO k = 1, nz
	        IF (vx(i,j,k) .NE. 0) THEN
	            IF ( (ABS(vx(i,j,k)-velox(i,j)) .GT. 100.) .OR. &
                (ABS(vy(i,j,k)-veloy(i,j)) .GT. 100. ) ) THEN
	                vx(i,j,k)=0
	                vy(i,j,k)=0
                ENDIF
	        ENDIF
	    ENDDO
    ENDIF	

    cnt=COUNT( (vx(i,j,1:nz) .ne. 0) .and. (vy(i,j,1:nz) .ne. 0))
    mask = ( (vx(i,j,1:nz) .ne. 0) .and. (vy(i,j,1:nz) .ne. 0))
    nn2(i,j)=cnt
    ! WEIGHTED MEAN X
    IF (cnt .NE. 0) THEN
        
        w_vx(i,j) = sum( vx(i,j,1:nz)*weightx(i,j,1:nz) , mask = mask ) / sum( weightx(i,j,1:nz) , mask = mask )

        IF (cnt .GE. 2) THEN
            stdx2(i,j) = sqrt( cnt * sum( weightx(i,j,1:nz) * &
            (vx(i,j,1:nz) - w_vx(i,j))**2 , mask = mask) / &
            sum( weightx(i,j,1:nz) , mask = mask ) / (cnt-1) )
        ENDIF


        ! WEIGHTED MEAN Y
        w_vy(i,j) = sum( vy(i,j,1:nz)*weighty(i,j,1:nz) , mask = mask ) / sum( weighty(i,j,1:nz), mask = mask )

        IF (cnt .GE. 2) THEN
            stdy2 (i,j) = sqrt( cnt * sum( weighty(i,j,1:nz) * &
                (vy(i,j,1:nz) - w_vy(i,j))**2 , mask = mask ) / &
                sum( weighty(i,j,1:nz) , mask = mask) / (cnt-1) )
        ENDIF
    ENDIF
	ENDDO
ENDDO

END SUBROUTINE cube2vel_feather_sp1


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

