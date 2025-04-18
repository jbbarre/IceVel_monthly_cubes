function findnoise,source
    
    source_decomp = strsplit(source,'/',/extract)
    for i=0,n_elements(source_decomp)-1 do source_decomp[i]=strmid(source_decomp[i],0,3)    
    ; L-band sensor
    w=where(source_decomp eq 'alo' or source_decomp eq 'ALO', cnt)
    if cnt ge 1 then begin
        ;; L-BAND sensor
        return, 1.0 / 4. /!dpi * 0.23622047 *365.25
    endif else begin
        return, 1.0 / 4. /!dpi * 0.056274579 *365.25
    ENDELSE
end

function unw2vel_feather_sp_fortran,icube,jcube,geo0,npix=npix,nrec=nrec,verbose=verbose,debug=debug

    path='/mnt/pennell-z0/eric/ANTARCTICA/'

    IF NOT KEYWORD_SET(npix) THEN npix=250
    IF NOT KEYWORD_SET(nrec) THEN nrec=250

    xmin = float(icube)*geo0.posting + geo0.xmin
    xmax = xmin + ( float(npix) - 1.) * geo0.posting
    ymax = geo0.ymax - float(jcube)*geo0.posting
    ymin = ymax - ( float(nrec) - 1.) * geo0.posting

    geoloc = {xmin:xmin,ymax:ymax,posting:geo0.posting,nrec:nrec,npix:npix}

    IF KEYWORD_SET(verbose) THEN BEGIN ;; {{{
        print,'==========================================================================================='
        print
        print,'unw_cubes/u_x'+string(icube,FORmat='(I05)')+ $
            '_y'+string(jcube,FORmat='(I05)')+'_post'+string(geo0.posting,FORmat='(I04)')+'.nc'
        print
    ENDIF ;; }}}

    cube=[]

    cube_file = 'unw_cubes/u_x'+string(icube,FORmat='(I05)')+ $
        '_y'+string(jcube,FORmat='(I05)')+'_post'+string(geo0.posting,FORmat='(I04)')+'.nc'
    IF NOT FILE_TEST('/u/pennell-z0/eric/ANTARCTICA/MOSAIC/450m/UNW/MOSAIC/'+cube_file) THEN GOTO,jump_skip

    if keyword_set(verbose) then print,'cube = concat_unw_cubes(',icube-245,',',jcube-245,',3,3)'
    cube = concat_unw_cubes(icube-245,jcube-245,3,3) ;nc2unw_cube('/u/pennell-z0/eric/ANTARCTICA/MOSAIC/450m/UNW/MOSAIC/'+cube_file)
    if keyword_set(debug) then  begin
        help,cube
        stop
    endif
    
    IF n_elements(cube) ge 2 THEN BEGIN

        v=cube.v
        w=where(finite(v) eq 0 or v eq 0,cnt)
        if cnt ne 0 then begin
            v[w] = 0.
        endif
        for i=0,n_elements(cube)-1 do v[*,*,i] = median(v[*,*,i],3)
        feathering = v
        feathering[*]=1.

        w=where(finite(v) eq 0 or v eq 0,cnt)
        if cnt ne 0 then begin
            v[w] = 0.
            feathering[w]=0.
        endif

        w=where(v eq 0,cnt)
        if cnt ne 0 then feathering[w]=0.
        w=0

        inc = cube.inc
        ;for i=0,n_elements(cube)-1 do inc[*,*,i] = median(inc[*,*,i],3)
        inc = inc[245:245+249,245:245+249,*]

        v = v[245:245+249,245:245+249,*]
        w=where(v eq 0,cnt)
        if cnt ne 0 then inc[w]=0.

        feather_size = 80
        noise = fltarr(250,250,n_elements(cube))
        f=feathering
        feathering = fltarr(250,250,n_elements(cube))
        for i=0,n_elements(cube)-1 do begin
            ; find sensor
            noise[*,*,i] = findnoise(cube[i].source) / ABS( cube[i].date2-cube[i].date1)/sin(inc[*,*,i])
            feathering[*,*,i] = (((MORPH_DISTANCE(reform(f[*,*,i]), NEIGHBOR_SAMPLING = 3))<feather_size)/float(feather_size))[245:245+249,245:245+249]
        ENDFoR
        f=0

        heading = cube.heading
        for i=0,n_elements(cube)-1 do heading[*,*,i] = median(heading[*,*,i],5)

        heading = heading[245:245+249,245:245+249,*]
        w=where(v eq 0,cnt)
        if cnt ne 0 then heading[w]=0.
        if cnt ne 0 then noise[w]=0.

        dem = geocode_im('/u/pennell-z0/eric/DATA_SERVER/ANTARCTICA/DEM/TanDEM-X_500m/TDX_DEM_500m.filtered.ASTER_PEN.BEDMAP2.dat',/float,geoloc=geoloc)
        
        slx = convol(dem,[-1.,0.,1.])/(2.*geoloc.posting)
        sly = convol(dem,transpose([1,0,-1]))/(2.*geoloc.posting)

        w=where(total(total(v,1),1) ne 0,cnt)
        if cnt eq 0 then goto,jump_skip ;; no good data
        if cnt ne 0 then begin
            inc = inc[*,*,w]
            v = v[*,*,w]
            heading = heading[*,*,w]
            feathering = feathering[*,*,w]
            noise = noise[*,*,w]
        endif
        dem = 0

        nn = size(v)
        if nn[0] ne 3 then goto,jump_skip 
        if nn[3] le 1 then goto,jump_skip ;; only one image - cannot reconstruct 2d speed - need 2 maps

        velox=fltarr(250,250)
        veloy=fltarr(250,250)
        stdx=fltarr(250,250)
        stdy=fltarr(250,250)
        errx=fltarr(250,250)
        erry=fltarr(250,250)
        nn1=intarr(250,250)
        cube=[]
        if keyword_set(debug) then stop
        
        if keyword_set(verbose) then help,v,inc,slx,sly,heading,feathering,noise
              
        a=call_external('/home/jeremie/ST_RELEASE/MOSAIC/FORTRAN/UNWRAP/unw2vel_feathering_sp.so','unw2vel_feathering_sp', $
            v, $
            inc, $
            slx, $
            sly, $
            heading, $
            feathering, $
            noise, $
            nn[1], $
            nn[2], $
            nn[3], $
            velox, $
            veloy, $
            stdx, $
            stdy,$
            errx, $
            erry, $
            nn1,$
            25.) ;; angle minimum between tracks
        
        if keyword_set(debug) then stop
    endif else begin
        jump_skip:
        return,{ $
            vx:fltarr(npix,nrec),$
            vy:fltarr(npix,nrec),$
            cnt:intarr(npix,nrec),$
            stdx:fltarr(npix,nrec),$
            stdy:fltarr(npix,nrec),$
            errx:fltarr(npix,nrec),$
            erry:fltarr(npix,nrec)}
    endelse

    IF total(velox) ne 0 THEN BEGIN
        return,{ $
            vx: velox , $
            vy: veloy , $
            cnt: nn1  , $
            stdx: stdx  , $
            stdy: stdy, $
            errx: errx, $
            erry: erry}
    ENDIF ELSE BEGIN
        goto,jump_skip
    ENDELSE 

end
