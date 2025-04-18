function cube2vel_feather_sp_fortran,cube,weight,debug=debug
    ;cube=cube0

    common env,svn_rep
    on_error, 2
    v = cube.v
    w=where(finite(abs(v)) eq 0,cnt)
    if cnt ne 0 then v[w] = complex(0.,0.)
    for i=0,n_elements(cube)-1 do v[*,*,i] = complex( median(float(v[*,*,i]),3),  median(imaginary(v[*,*,i]),3))

    feathering = abs(v)
    feathering[*]=1.
    ;feathering[0:1,*,*]=0.
    ;feathering[*,0:1,*]=0.

    w=where(finite(abs(v)) eq 0,cnt)
    if cnt ne 0 then feathering[w]=0.

    w=where(v eq 0,cnt)
    if cnt ne 0 then feathering[w]=0.
    w=0

    nx = n_elements(cube[0].v[*,0])
    ny = n_elements(cube[0].v[0,*])
    
    feather_size = 80
    for i=0,n_elements(cube)-1 do begin
        f = fltarr(250+2*feather_size,250+2*feather_size)
        for j=0,feather_size-1 do $
            f[j,feather_size:feather_size+249]=feathering[0,*,i]
        for j=0,feather_size-1 do $
            f[feather_size:feather_size+249,j]=feathering[*,0,i]
        for j=feather_size+250,250+2*feather_size-1 do $
            f[j,feather_size:feather_size+249]=feathering[249,*,i]
        for j=feather_size+250,250+2*feather_size-1 do $
            f[feather_size:feather_size+249,j]=feathering[*,249,i]
        f[feather_size,feather_size]=reform(feathering[*,*,i]) ;; take care of image edges ..
        f[0:feather_size-1,0:feather_size-1]=feathering[0,0,i]
        f[0:feather_size-1,feather_size+250:*]=feathering[0,249,i]
        f[feather_size+250:*,feather_size+250:*]=feathering[249,249,i]
        f[feather_size+250:*,0:feather_size-1]=feathering[249,0,i]
        feathering[*,*,i] = (((MORPH_DISTANCE(f, NEIGHBOR_SAMPLING = 3))<feather_size)/float(feather_size))[feather_size:feather_size+249,feather_size:feather_size+249]
    ENDFoR
    f=0
        ;feathering[*,*,i] = ( (smooth( reform(feathering[*,*,i]), 30 , /EDGE_TRUNCATE ) - 0.5) > 0 )*2.
    
    IF KEYWORD_SET(debug) then help,weight
    IF KEYWORD_SET(debug) then help,feathering
    weight_mapx = feathering
    feathering=0
    weight_mapy=weight_mapx

    for i=0,n_elements(cube)-1 do weight_mapx[*,*,i] = reform(weight_mapx[*,*,i])*float(weight[i])
    for i=0,n_elements(cube)-1 do weight_mapy[*,*,i] = reform(weight_mapy[*,*,i])*imaginary(weight[i])
    IF KEYWORD_SET(debug) then help,weight_mapx

    IF n_elements(cube) eq 1 THEN BEGIN

        w=where(float(v) ne 0)
        nn = lonarr( n_elements(v[*,0]),n_elements(v[0,*]))
        nn[w]=1

        return,{$
            weighted_mean :v, $
            weighted_std  :complexarr( n_elements(v[*,0]),n_elements(v[0,*]) ),$
            n_filt_values :nn,$
            med           :v,$
            std           : cube[0].err*nn,$
            err           : cube[0].err*nn,$
            n_values      :nn}
    ENDIF

    nn = size(v)
    if nn[0] eq 3 and nn[1] eq nx and nn[2] eq ny and nn[3] ge 2 then begin

        velox = fltarr(nn[1],nn[2])
        veloy = fltarr(nn[1],nn[2])
        stdx = fltarr(nn[1],nn[2]);,/nozero)
        stdy = fltarr(nn[1],nn[2]);,/nozero)
        stdx2 = fltarr(nn[1],nn[2]);,/nozero)
        stdy2 = fltarr(nn[1],nn[2]);,/nozero)
        errfx = fltarr(nn[1],nn[2]);,/nozero)
        errfy = fltarr(nn[1],nn[2]);,/nozero)
        w_vx = fltarr(nn[1],nn[2]);,/nozero)
        w_vy = fltarr(nn[1],nn[2]);,/nozero)
        nn1 = intarr(nn[1],nn[2]);,/nozero)
        nn2 = intarr(nn[1],nn[2]);,/nozero)

        velox=fltarr(nx,ny)
        veloy=fltarr(nx,ny)
        stdx=fltarr(nx,ny)
        stdy=fltarr(nx,ny)
        stdx2=fltarr(nx,ny)
        stdy2=fltarr(nx,ny)
        errfx=fltarr(nx,ny)
        errfy=fltarr(nx,ny)
        w_vx=fltarr(nx,ny)
        w_vy=fltarr(nx,ny)
        nn1=intarr(nx,ny)
        nn2=intarr(nx,ny)

        IF KEYWORD_SET(DEBUG) THEN BEGIN
            print,'====='
            print,'Dimension ',nn[1],nn[2],nn[3]
            help,velox,veloy,stdx,stdy,stdy2,errfx,errfy,w_vx,w_vy,nn1,nn2
            help,v
            help,float(cube[*].err)
            help,float(imaginary(cube[*].err))
            print,'---=====---'
        ENDIF    
        a=call_external(svn_rep+'/MOSAIC/FORTRAN/CUBE2VEL/cube2vel_feather_sp.so','cube2vel_feather_sp', $
            float(v), $
            floaT(imaginary(v)), $
            nn[1], $
            nn[2], $
            nn[3], $
            float(cube[*].err), $
            float(imaginary(cube[*].err)), $
            weight_mapx, $
            weight_mapy, $
            velox, $
            veloy, $
            stdx, $
            stdy,$
            stdx2,$
            stdy2,$
            errfx,$
            errfy,$
            w_vx,$
            w_vy,$
            nn1,$
            nn2)
    endif else begin
        return,0
    endelse
    return,{weighted_mean :complex(w_vx,w_vy), $
    weighted_std  :complex(stdx2,stdy2),$
    n_filt_values :nn2,$
    med           :complex(velox,veloy),$
    std           :complex(stdx,stdy),$
    err           :complex(errfx,errfy),$
    n_values      :nn1}
end
