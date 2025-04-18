function unw2vel_gre_feather_sp_fortran,icube,jcube,geo0,npix=npix,nrec=nrec,SENSORS=SENSORS,verbose=verbose

    path='/mnt/pennell-z0/eric/GREENLAND/'

    IF NOT KEYWORD_SET(SENSORS) THEN BEGIN
        SENSORS = [$
            'alos2009_greenland',$
            'envisat_greenland_2010',$
            'RADARSAT_greenland_2009',$
            'ERS1992_greenland']
    ENDIF
    IF NOT KEYWORD_SET(npix) THEN npix=250
    IF NOT KEYWORD_SET(nrec) THEN nrec=250

    xmin = float(icube)*geo0.posting + geo0.xmin
    xmax = xmin + ( float(npix) - 1.) * geo0.posting
    ymax = geo0.ymax - float(jcube)*geo0.posting
    ymin = ymax - ( float(nrec) - 1.) * geo0.posting

    geoloc = {xmin:xmin,ymax:ymax,posting:geo0.posting,nrec:nrec,npix:npix}

    IF KEYWORD_SET(verbose) THEN BEGIN
        print,'==========================================================================================='
        print
        print,'unw_cubes/c_x'+string(icube,FORmat='(I05)')+ $
            '_y'+string(jcube,FORmat='(I05)')+'_post'+string(geo0.posting,FORmat='(I04)')+'.nc'
        print
    ENDIF

    cube=[]
    FOR k=0,n_elements(SENSORS)-1 DO BEGIN

        cube_file = 'unw_cubes/c_x'+string(icube,FORmat='(I05)')+ $
            '_y'+string(jcube,FORmat='(I05)')+'_post'+string(geo0.posting,FORmat='(I04)')+'.nc'
        IF keyword_set(verbose) then print,path+SENSORS[k]+'/MOSAIC/'+cube_file,file_test(path+SENSORS[k]+'/MOSAIC/'+cube_file)
        IF file_test(path+SENSORS[k]+'/MOSAIC/'+cube_file) THEN $
            cube = [cube,nc2unw_cube(path+SENSORS[k]+'/MOSAIC/'+cube_file)]
        help,cube
    ENDFOR

    IF n_elements(cube) ge 1 THEN BEGIN

        v=cube.v
        
        feathering = v
        feathering[*]=1.
        feathering[0:3,*,*]=0.
        feathering[*,0:3,*]=0.

        w=where(finite(v) eq 0,cnt)
        if cnt ne 0 then begin
            v[w]=0.
            feathering[w]=0.
        endif

        w=where(v eq 0,cnt)
        if cnt ne 0 then feathering[w]=0.

        for i=0,n_elements(cube)-1 do $
            feathering[*,*,i] = (smooth( reform(feathering[*,*,i]), 30 ) - 0.5)*2.

        help,feathering

        heading = cube.heading
        w=where(finite(heading) eq 0,cnt)
        if cnt ne 0 then heading[w]=0.

        inc = cube.inc
        w=where(finite(inc) eq 0,cnt)
        if cnt ne 0 then inc[w]=0.

        dem = geocode_im('/mnt/pennell-z0/eric/DATA_SERVER/GREENLAND/DEM/gimpdem100.dat',geoloc=geoloc,/float)
        slx = convol(dem,[-1.,0.,1.])/(2.*geoloc.posting)
        sly = convol(dem,transpose([1,0,-1]))/(2.*geoloc.posting)

        dem = 0

        nn = size(cube.v)

        velox=fltarr(250,250)
        veloy=fltarr(250,250)
        stdx=fltarr(250,250)
        stdy=fltarr(250,250)
        nn1=intarr(250,250)
        
        ;stop

        cube=[]
        a=call_external('/home/jeremie/ST_RELEASE/MOSAIC/UNWRAP/unw2vel_feathering_sp.so','unw2vel_feathering_sp', $
            v, $
            inc, $
            slx, $
            sly, $
            heading, $
            feathering, $
            nn[1], $
            nn[2], $
            nn[3], $
            velox, $
            veloy, $
            stdx, $
            stdy,$
            nn1)
    endif else begin
        return,{ $
            vx:fltarr(npix,nrec),$
            vy:fltarr(npix,nrec),$
            cnt:intarr(npix,nrec),$
            stdx:fltarr(npix,nrec),$
            stdy:fltarr(npix,nrec)}
    endelse
    
    IF total(velox) ne 0 THEN BEGIN
        return,{ $
            vx: velox , $
            vy: veloy , $
            cnt: nn1  , $
            stdx: stdx  , $
            stdy: stdy }
    ENDIF

end
