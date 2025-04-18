function cube2vel_sp_fortran,cube,weight,debug=debug
    ;cube=cube0

    on_error, 2
    v = cube.v
    w=where(finite(abs(v)) eq 0,cnt)
    if cnt ne 0 then v[w] = complex(0.,0.)

;    IF n_elements(cube) eq 2 THEN BEGIN
;        
;        v[w] = complex(!values.f_nan,!values.f_nan)
;        
;        nn = total(finite(float(v)),3)
;        
;        v0=mean(v,dim=3,/nan)
;        w=where(finite(abs(v0)) eq 0,cnt)
;        if cnt ne 0 then v0[w] = complex(0.,0.)
;
;        return,{$
;            weighted_mean :v0, $
;            weighted_std  :complexarr( n_elements(v0[*,0]),n_elements(v0[0,*]) ) + cube[0].err,$
;            n_filt_values :nn,$
;            med           :v0,$
;            std           :complexarr( n_elements(v0[*,0]),n_elements(v0[0,*]) ) + cube[0].err,$
;            err           :complexarr( n_elements(v0[*,0]),n_elements(v0[0,*]) ) + cube[0].err,$
;            n_values      :nn}
;
;    ENDIF
;
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
    if nn[1] eq 250 and nn[2] eq 250 and nn[3] ge 2 then begin

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

        velox=fltarr(250,250)
        veloy=fltarr(250,250)
        stdx=fltarr(250,250)
        stdy=fltarr(250,250)
        stdx2=fltarr(250,250)
        stdy2=fltarr(250,250)
        errfx=fltarr(250,250)
        errfy=fltarr(250,250)
        w_vx=fltarr(250,250)
        w_vy=fltarr(250,250)
        nn1=intarr(250,250)
        nn2=intarr(250,250)

        IF KEYWORD_SET(DEBUG) THEN BEGIN
            print,'====='
            print,'Dimension ',nn[1],nn[2],nn[3]
            help,velox,veloy,stdx,stdy,stdy2,errfx,errfy,w_vx,w_vy,nn1,nn2
            help,v
            help,float(cube[*].err)
            help,float(imaginary(cube[*].err))
            print,'---=====---'
        ENDIF    
        a=call_external('/home/jeremie/ST_RELEASE/MOSAIC/CUBE2VEL/cube2vel_sp.so','cube2vel_sp', $
            float(v), $
            floaT(imaginary(v)), $
            nn[1], $
            nn[2], $
            nn[3], $
            float(cube[*].err), $
            float(imaginary(cube[*].err)), $
            float(weight), $
            float(imaginary(weight)), $
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
