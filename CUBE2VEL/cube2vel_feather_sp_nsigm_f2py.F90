! FILE:
SUBROUTINE cube2vel_feather_nsigm_sp1(vx0,vy0,errx0,erry0,weightx0,weighty0,nsigma, & ! input
        velox,veloy,stdx,stdy,stdx2,stdy2,errfx,errfy,w_vx,w_vy,nn,nn2, & ! output
        nx,ny,nz) ! hidden

    implicit none;

    INTEGER*2 nx, ny, nz
    REAL*4 vx0(nx,ny,nz), vy0(nx,ny,nz)
    REAL*4 errx0(nz), erry0(nz)
    REAL*4 weightx0(nx,ny,nz), weighty0(nx,ny,nz) 
    INTEGER*2 nz2
    REAL*4 velox(nx,ny), veloy(nx,ny)
    REAL*4 med,medx,medy,nsigma
    REAL*4 vtmpx,vtmpy
    REAL*4 stdx(nx,ny), stdy(nx,ny)
    REAL*4 stdx2(nx,ny), stdy2(nx,ny)
    REAL*4 w_vx(nx,ny), w_vy(nx,ny)
    REAL*4 errfx(nx,ny), errfy(nx,ny)
    REAL*4 errx(5*nz), erry(5*nz)
    REAL*4 weightx(5*nz), weighty(5*nz)
    REAL*4 vx(5*nz), vy(5*nz)
    INTEGER*2 nn(nx,ny),nn2(nx,ny)
    INTEGER*2 w(nz)
    INTEGER*2 i,j,k,k0
    INTEGER*4 cnt
    LOGICAL mask(5*nz)
    !print *,'toto'

    !f2py intent(in) vx0,vy0,errx0,erry0,weightx0,weighty0,nsigma
    !f2py intent(hide),depend(vx0) :: nx=shape(vx0,0), ny=shape(vx0,1), nz=shape(vx0,2)
    !f2py intent(out) :: velox,veloy,stdx,stdy,stdx2,stdy2,errfx,errfy,w_vx,w_vy,nn,nn2

    w_vx(1:nx,1:ny)=0.
    w_vy(1:nx,1:ny)=0.
    velox(1:nx,1:ny)=0
    veloy(1:nx,1:ny)=0
    stdx(1:nx,1:ny)=0
    stdy(1:nx,1:ny)=0
    stdx2(1:nx,1:ny)=0
    stdy2(1:nx,1:ny)=0

    print *,nx,ny,nz
    !vx(1:nx,1:ny,1:5*nz)=0.
    !vy(1:nx,1:ny,1:5*nz)=0
    !vx(1:nx,1:ny,1:nz) = vx0
    !vy(1:nx,1:ny,1:nz) = vy0

    !print *,nx,ny,nz

    !DO i = 2, nx-1
    !DO j = 2, ny-1
    !DO k = 1,nz 
    !IF ((vx0(i-1,j,k) .NE. 0) .AND. (vx0(i+1,j,k) .NE. 0)) THEN
    !    vx(i,j,nz+k) = (vx0(i-1,j,k) + vx0(i+1,j,k)) / 2
    !    vy(i,j,nz+k) = (vy0(i-1,j,k) + vy0(i+1,j,k)) / 2
    !ENDIF
    !ENDDO

    !DO k = 1,nz
    !IF ((vx0(i,j-1,k) .NE. 0) .AND. (vx0(i,j+1,k) .NE. 0)) THEN
    !    vx(i,j,2*nz+k) = (vx0(i,j-1,k) + vx0(i,j+1,k)) / 2
    !    vy(i,j,2*nz+k) = (vy0(i,j-1,k) + vy0(i,j+1,k)) / 2
    !ENDIF
    !ENDDO

    !DO k = 1,nz
    !IF ((vx0(i-1,j-1,k) .NE. 0) .AND. (vx0(i+1,j+1,k) .NE. 0)) THEN
    !    vx(i,j,3*nz+k) = (vx0(i-1,j-1,k) + vx0(i+1,j+1,k)) / 2
    !    vy(i,j,3*nz+k) = (vy0(i-1,j-1,k) + vy0(i+1,j+1,k)) / 2
    !ENDIF
    !ENDDO

    !DO k = 1,nz
    !IF ((vx0(i+1,j-1,k) .NE. 0) .AND. (vx0(i-1,j+1,k) .NE. 0)) THEN
    !    vx(i,j,4*nz+k) = (vx0(i+1,j-1,k) + vx0(i-1,j+1,k)) / 2
    !    vy(i,j,4*nz+k) = (vy0(i+1,j-1,k) + vy0(i-1,j+1,k)) / 2
    !ENDIF
    !ENDDO
    !ENDDO
    !ENDDO

    !errx(1:nz) = errx0
    !errx(nz+1:2*nz) = errx0
    !errx(2*nz+1:3*nz) = errx0
    !errx(3*nz+1:4*nz) = errx0
    !errx(4*nz+1:5*nz) = errx0

    !erry(1:nz) = erry0
    !erry(nz+1:2*nz) = erry0
    !erry(2*nz+1:3*nz) = erry0
    !erry(3*nz+1:4*nz) = erry0
    !erry(4*nz+1:5*nz) = erry0

    !weightx(1:nx,1:ny,1:nz) = weightx0/4.
    !weightx(1:nx,1:ny,nz+1:2*nz) = weightx0/4.
    !weightx(1:nx,1:ny,2*nz+1:3*nz) = weightx0/4.
    !weightx(1:nx,1:ny,3*nz+1:4*nz) = weightx0/4.
    !weightx(1:nx,1:ny,4*nz+1:5*nz) = weightx0/4.

    !weighty(1:nx,1:ny,1:nz) = weighty0/4.
    !weighty(1:nx,1:ny,nz+1:2*nz) = weighty0/4.
    !weighty(1:nx,1:ny,2*nz+1:3*nz) = weighty0/4.
    !weighty(1:nx,1:ny,3*nz+1:4*nz) = weighty0/4.
    !weighty(1:nx,1:ny,4*nz+1:5*nz) = weighty0/4.

    nz2 = 5*nz
    ! COMPUTE MEDIAN on third DIMENSION
    DO i = 1, nx
    DO j = 1, ny

    ! re-initiliaze vx,vy,weights
    vx(:)=0.
    vy(:)=0.
    weightx(:)=0.
    weighty(:)=0.

    ! COMPUTE MEDIAN FOR VX
    vx(1:nz) = vx0(i,j,1:nz)
    vy(1:nz) = vy0(i,j,1:nz)
    weightx(1:nz) = weightx0(i,j,1:nz)
    weighty(1:nz) = weighty0(i,j,1:nz)

    DO k = 1, nz
        IF ( (i .GT. 1) .AND. (i .LT. nx) ) THEN
            IF ((vx0(i-1,j,k) .NE. 0) .AND. (vx0(i+1,j,k) .NE. 0)) THEN
                vx(nz+k) = (vx0(i-1,j,k) + vx0(i+1,j,k)) / 2
                vy(nz+k) = (vy0(i-1,j,k) + vy0(i+1,j,k)) / 2
                weightx(nz+k) = weightx0(i,j,k)/4.
                weighty(nz+k) = weighty0(i,j,k)/4.
            ENDIF
        ENDIF
        IF ( (j .GT. 1) .AND. (j .LT. ny) ) THEN
            IF ( (vx0(i,j-1,k) .NE. 0) .AND. (vx0(i,j+1,k) .NE. 0) ) THEN
                vx(2*nz+k) = (vx0(i,j-1,k) + vx0(i,j+1,k)) / 2    
                vy(2*nz+k) = (vy0(i,j-1,k) + vy0(i,j+1,k)) / 2
                weightx(2*nz+k) = weightx0(i,j,k)/4.
                weighty(2*nz+k) = weighty0(i,j,k)/4.
            ENDIF
        ENDIF

        IF ( (i .GT. 1) .AND. (j .LT. ny) ) THEN
            IF ((vx0(i-1,j-1,k) .NE. 0) .AND. (vx0(i+1,j+1,k) .NE. 0) ) THEN
                vx(3*nz+k) = (vx0(i-1,j-1,k) + vx0(i+1,j+1,k)) / 2
                vy(3*nz+k) = (vx0(i-1,j-1,k) + vx0(i+1,j+1,k)) / 2
                weightx(3*nz+k) = weightx0(i,j,k)/4.
                weighty(3*nz+k) = weighty0(i,j,k)/4.
            ENDIF
        ENDIF 
        
        IF ( (j .GT. 1) .AND. (i .LT. nx) ) THEN
            IF ((vx0(i+1,j-1,k) .NE. 0) .AND. (vx0(i-1,j+1,k) .NE. 0) ) THEN
                vx(4*nz+k) = (vx0(i+1,j-1,k) + vx0(i-1,j+1,k)) / 2
                vy(4*nz+k) = (vy0(i+1,j-1,k) + vy0(i-1,j+1,k)) / 2
                weightx(4*nz+k) = weightx0(i,j,k)/4.
                weighty(4*nz+k) = weighty0(i,j,k)/4.
            ENDIF
        ENDIF
    ENDDO
    cnt=COUNT( (vx .ne. 0) .and. (vy .ne. 0))

    ! COUNT ORIGINAL vx0 \= 0 
    nn(i,j)=COUNT(vx0(i,j,1:nz) .ne. 0)

    ! WEIGHTED STD X

    errfx(i,j) = sqrt( SUM( (errx0*weightx0(i,j,1:nz))**2 , mask = ( vx0(i,j,1:nz) .ne. 0 ))/&
        SUM( weightx0(i,j,1:nz) , mask = ( vx0(i,j,1:nz) .ne. 0 ) )**2 ) 

    errfy(i,j) = sqrt(  SUM( (erry0*weighty0(i,j,1:nz))**2, mask = ( vx0(i,j,1:nz) .ne. 0 ))/&
        SUM( weighty0(i,j,1:nz) , mask = ( vx0(i,j,1:nz) .ne. 0 ) )**2 )

    IF (cnt .GT. 1) THEN
        CALL MDIAN( PACK( vx, mask= (vx .ne. 0) .and. (vy .ne. 0)),cnt,med)
        velox(i,j) = med !sum( vx*weightx , mask = (vx .ne. 0) .and. (vy .ne. 0)) / sum( weighty, mask = (vx .ne. 0) .and. (vy .ne. 0)) !med
    ELSE IF (cnt .EQ. 1) THEN
        velox(i,j) = SUM( vx, mask= (vx .ne. 0) .and. (vy .ne. 0))
    ELSE
        velox(i,j)=0 
    ENDIF
    IF (cnt .GT. 3) THEN    
        ! COMPUTE STD FOR VX
        !med = sum( vx , mask =  (vx .ne. 0) .and. (vy .ne. 0)) / cnt
        stdx(i,j) = sqrt( sum( (vx - velox(i,j))**2 , mask = (vx .ne. 0) .and. (vy .ne. 0)) / (cnt-1))
    ENDIF
    vtmpx=stdx(i,j)!0.1*abs(med) 
    IF (vtmpx .LT. 1) THEN 
        vtmpx = 1
    ENDIF
    !IF (vtmpx .GT. 50.) THEN
    !    vtmpx = 50.
    !ENDIF
    ! COMPUTE MEDIAN FOR VY
    IF (cnt .GT. 1) THEN
        CALL MDIAN( PACK( vy, mask= (vx .ne. 0) .and. (vy .ne. 0)),cnt,med)
        veloy(i,j) = med !sum( vy*weighty , mask = (vx .ne. 0) .and. (vy .ne. 0)) / sum( weighty, mask = (vx .ne. 0) .and. (vy .ne. 0)) !med
    ELSE IF (CNT .EQ. 1) THEN
        veloy(i,j)=SUM(vy, mask= (vx .ne. 0) .and. (vy .ne. 0))
    ELSE
        veloy(i,j)=0
    ENDIF
    IF (cnt .GT. 3) THEN
        ! COMPUTE STD FOR Vy
        !med = sum( vy , mask =  (vx .ne. 0) .and. (vy .ne. 0)) / cnt
        stdy(i,j) = sqrt( sum( (vy - veloy(i,j))**2 , mask = (vx .ne. 0) .and. (vy .ne. 0)) / (cnt-1))
    ENDIF
    vtmpy = stdy(i,j) !0.1*abs(med)
    IF (vtmpy .LT. 1) THEN
        vtmpy = 1
    ENDIF
    !IF (vtmpy .GT. 50.) THEN
    !    vtmpy = 50.
    !ENDIF
    !print *,i,j 
    IF (cnt .GT. 3) THEN

        ! REMOVE VALUES MORE THAN 1 STD AWAY FOR MEDIAN
        DO k = 1, nz2
        IF (vx(k) .NE. 0) THEN
            !print *,vx(k),velox(i,j),vtmpx,vy(k),veloy(i,j),vtmpy
            IF ( (ABS(vx(k)-velox(i,j)) .GT. nsigma*vtmpx ) .OR. &
                (ABS(vy(k)-veloy(i,j)) .GT. nsigma*vtmpy ) ) THEN

            
            vx(k)=0
            vy(k)=0
            ENDIF
        ENDIF
        ENDDO
    ENDIF	
    
cnt=COUNT( (vx .ne. 0) .and. (vy .ne. 0))
mask = ( (vx .ne. 0) .and. (vy .ne. 0))
nn2(i,j)=COUNT( (vx(1:nz) .ne. 0) .and. (vy(1:nz) .ne. 0))

! WEIGHTED MEAN X
IF (cnt .NE. 0) THEN


    w_vx(i,j) = sum( vx*weightx , mask = mask ) / sum( weightx , mask = mask )

    IF (cnt .GE. 2) THEN
        stdx2(i,j) = sqrt( cnt * sum( weightx * &
            (vx - w_vx(i,j))**2 , mask = mask) / &
            sum( weightx , mask = mask ) / (cnt-1) )
    ENDIF


    ! WEIGHTED MEAN Y
    w_vy(i,j) = sum( vy*weighty , mask = mask ) / sum( weighty, mask = mask )

    IF (cnt .GE. 2) THEN
        stdy2(i,j) = sqrt( cnt * sum( weighty * &
            (vy - w_vy(i,j))**2 , mask = mask ) / &
            sum( weighty , mask = mask) / (cnt-1) )
    ENDIF
ENDIF
ENDDO
ENDDO

END SUBROUTINE cube2vel_feather_nsigm_sp1


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

