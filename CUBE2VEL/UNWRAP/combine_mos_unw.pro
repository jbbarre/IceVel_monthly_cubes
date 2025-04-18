
path='/mnt/pennell-z0/eric/ANTARCTICA/'

SENSORS = [$
    ;'ALOS_2006',$
    'ALOS_2007',$
    'ALOS_2008',$
    'ALOS_2009',$
    ;'ALOS_2010',$
    'envisat2007',$
    'envisat2008',$
    'envisat2009',$
    ;'ERS1992',$
    ;'ERS1994',$
    'ERS1996',$
    'RADARSAT1997',$
    ;'RADARSAT2000',$
    ;'RADARSAT2002',$
    ;'RADARSAT2002-2003',$
    ;'RADARSAT2000-2001',$
    ;'RADARSAT2003-2004',$
    ;'RADARSAT2004-2005',$
    ;'RADARSAT2005-2006',$
    'RADARSAT-2_2009',$
    'RADARSAT-2_2011',$
    'RADARSAT-2_2013',$
    'RADARSAT-2_2015',$
    'RADARSAT-2_2016',$
    ;'SENTINEL1',$
    ;'SENTINEL1/24d',$
    ;'TDX_2011',$
    ;'TDX_2012',$
    'TDX_2013']

geo0=load_geo_param()

a=fltarr(2,2)
b=fltarr(2,2)
c=fltarr(2,2)
unity=a & unity[*,0]=[1.,0.] & unity[*,1]=[0.,1.]
cr = string("15b) ;"
FORm="($,i5,a,a)"

vel_final = complexarr(geo0.npix,geo0.nrec)

;;size of the cubes
nrec = 250
npix = 250

;FOR i=0.,geo.npix-1,990. DO BEGIN
;    FOR j=0.,geo.nrec-1,990. DO BEGIN
FOR i=0.,geo0.npix-1,nrec-5 DO BEGIN
    FOR j=0.,geo0.nrec-1,npix-5 DO BEGIN

        xmin = float(i)*geo0.posting + geo0.xmin
        xmax = xmin + ( float(npix) - 1.) * geo0.posting
        ymax = geo0.ymax - float(j)*geo0.posting
        ymin = ymax - ( float(nrec) - 1.) * geo0.posting

        geoloc = {xmin:xmin,ymax:ymax,posting:geo0.posting,nrec:nrec,npix:npix}

        print,'==========================================================================================='
        print
        print,'unw_cubes/c_x'+string(i,FORmat='(I05)')+ $
            '_y'+string(j,FORmat='(I05)')+'_post'+string(geo0.posting,FORmat='(I04)')+'.nc'
        print

        ;; load list of cubes {{{
        cube=[]
        FOR k=0,n_elements(SENSORS)-1 DO BEGIN

            cube_file = 'unw_cubes/c_x'+string(i,FORmat='(I05)')+ $
                '_y'+string(j,FORmat='(I05)')+'_post'+string(geo0.posting,FORmat='(I04)')+'.nc'
            print,path+SENSORS[k]+'/MOSAIC/'+cube_file,file_test(path+SENSORS[k]+'/MOSAIC/'+cube_file)
            IF file_test(path+SENSORS[k]+'/MOSAIC/'+cube_file) THEN $
                cube = [cube,nc2unw_cube(path+SENSORS[k]+'/MOSAIC/'+cube_file)]
        ENDFOR

        IF n_elements(cube) ge 1 THEN BEGIN

            v_final = complexarr(npix,nrec)

            dem = geocode_im('/mnt/pennell-z0/eric/DATA_SERVER/ANTARCTICA/DEM/BEDMAP2/bedmap2_surface.dat',geoloc=geoloc,/float)
            slx = convol(dem,[-1.,0.,1.])/(2.*geoloc.posting)
            sly = convol(dem,transpose([1,0,-1]))/(2.*geoloc.posting) 

            dem = 0

            ;; loop over x ,y axis 
            FOR ii=00,npix-1 DO BEGIN

                IF (ii mod 100) eq 0 THEN print,FORm=FORm,ii/999.*100.,' %',cr

                FOR jj=0,nrec-1 DO BEGIN

                    VX = []
                    VY = []
                    ;; loop over phases
                    FOR kk=0,n_elements(cube)-2 DO BEGIN
                        IF cube[kk].v[ii,jj] ne 0 and finite(cube[kk].v[ii,jj]) ne 0 THEN BEGIN

                            FOR mm=kk+1,n_elements(cube)-1 DO BEGIN
                                IF cube[mm].v[ii,jj] ne 0 and finite(cube[mm].v[ii,jj]) ne 0 THEN BEGIN

                                    alpha = cube[mm].heading[ii,jj]*!dtor
                                    beta  = cube[kk].heading[ii,jj]*!dtor

                                    ;; angle > 30 degrees -> combine phase {{{                                        
                                    IF ABS(min( [360.*!dtor - abs(alpha - beta),abs(alpha - beta) ])) GE 30*!dtor THEN BEGIN

                                        alpha=alpha-beta

                                        vra = cube[kk].v[ii,jj]
                                        vrd = cube[mm].v[ii,jj]

                                        tetadd = cube[kk].inc[ii,jj]
                                        tetaaa = cube[mm].inc[ii,jj]

                                        slopex = slx[ii,jj]
                                        slopey = sly[ii,jj]

                                        A[*,0] = [cos(beta), cos(alpha+beta)]
                                        A[*,1] = [sin(beta), sin(beta+alpha)]

                                        sinalphak_square = sin(alpha)*sin(alpha) > 0.000001

                                        B[*,0] = [ 1. , -cos(alpha) ] / sinalphak_square
                                        B[*,1] = [ -cos(alpha) , 1. ] / sinalphak_square

                                        AB = B#A

                                        C[*,0] = [slopex/tan(tetaaa) , slopey/tan(tetaaa)]
                                        C[*,1] = [slopex/tan(tetadd) , slopey/tan(tetadd)]

                                        D=AB#invert(unity-C#AB)

                                        VX = [VX , vra*D[0,0] + vrd*D[1,0]]
                                        VY = [VY , vra*D[0,1] + vrd*D[1,1]]
                                    ENDIF
                                    ;; }}}
                                ENDIF
                            ENDFOR
                        ENDIF
                    ENDFOR

                    IF n_elements(vx) ge 3 THEN V_final[ii,jj] = complex(median(vx),median(vy))
                    IF n_elements(vx) eq 2 THEN V_final[ii,jj] = complex(mean(vx),mean(vy))
                    IF n_elements(vx) eq 1 THEN V_final[ii,jj] = complex(vx,vy)
                ENDFOR
            ENDFOR

            IF total(float(v_final)) ne 0 THEN vel_final[i,j] = v_final[0:(npix-1)<(geo0.npix-1L-i),0:(nrec-1)<(geo0.nrec-1L-j)]

        ENDIF
    ENDFOR
ENDFOR

save,filename='unw_final.sav',vel_final

stop
end
