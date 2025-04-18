
geo0=load_geo_param()

cr = string("15b) ;"
FORm="($,i5,a,a)"

v = fltarr(geo0.npix,geo0.nrec)
inc = fltarr(geo0.npix,geo0.nrec)
heading = fltarr(geo0.npix,geo0.nrec)

vcount = intarr(geo0.npix,geo0.nrec)
inccount = intarr(geo0.npix,geo0.nrec)
headcount = intarr(geo0.npix,geo0.nrec)
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

        cube_file = 'unw_cubes/c_x'+string(i,FORmat='(I05)')+ $
            '_y'+string(j,FORmat='(I05)')+'_post'+string(geo0.posting,FORmat='(I04)')+'.nc'
        IF file_test(cube_file) THEN BEGIN
            cube = nc2unw_cube(cube_file)

            ;; loop over x ,y axis 
            FOR kk=0,n_elements(cube)-1 DO BEGIN
                
                inc0=reform(cube[kk].inc[*,*])
                w=where(finite(inc0) eq 0,cnt)
                if cnt ne 0 then inc0[w]=0.
                inc[i:i+249,j:j+249] += inc0
                inccount[i:i+249,j:j+249] += (inc0 NE 0)

                heading0=reform(cube[kk].heading[*,*])
                w=where(finite(heading0) eq 0,cnt)
                 if cnt ne 0 then heading0[w]=0.
                heading[i:i+249,j:j+249] += heading0
                headcount[i:i+249,j:j+249] += (heading0 NE 0)

                v0 = reform(cube[kk].v[*,*])
                w=where(finite(v0) eq 0,cnt)
                 if cnt ne 0 then v0[w]=0.
                v[i:i+249,j:j+249] += v0
                vcount[i:i+249,j:j+249] += (v0 NE 0)

            ENDFOR

        ENDIF

    ENDFOR
ENDFOR


stop
end
