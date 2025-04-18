function cube2vel_fortran,cube
    ;cube=cube0

    IF n_elements(cube) eq 1 THEN BEGIN
        return,{$
            weighted_mean :reFORm(cube[0].v), $
            weighted_std  :complexarr(n_elements(cube[0].v[*,0]),n_elements(cube[0].v[0,*]) ) + cube[0].err,$
            n_filt_values :finite(abs(cube[0].v)),$
            med           :reform(cube[0].v),$
            std           :complexarr(n_elements(cube[0].v[*,0]),n_elements(cube[0].v[0,*]) ) + cube[0].err,$
            err           :complexarr(n_elements(cube[0].v[*,0]),n_elements(cube[0].v[0,*]) ) + cube[0].err,$
            n_values      :finite(abs(cube[0].v))}

    ENDIF

    v = cube.v
    w=where(finite(abs(v)) eq 0,cnt)
    if cnt ne 0 then v[w] = complex(0.,0.)
    nn = size(v)
    if nn[1] eq 250 or nn[2] eq 250 then begin

        velox=fltarr(nn[1],nn[2],/nozero)
        veloy=fltarr(nn[1],nn[2],/nozero)
        stdx=fltarr(nn[1],nn[2],/nozero)
        stdy=fltarr(nn[1],nn[2],/nozero)
        stdx2=fltarr(nn[1],nn[2],/nozero)
        stdy2=fltarr(nn[1],nn[2],/nozero)
        errfx=fltarr(nn[1],nn[2],/nozero)
        errfy=fltarr(nn[1],nn[2],/nozero)
        w_vx=fltarr(nn[1],nn[2],/nozero)
        w_vy=fltarr(nn[1],nn[2],/nozero)
        nn1 = intarr(nn[1],nn[2],/nozero)
        nn2=intarr(nn[1],nn[2],/nozero)

        a=call_external('/home/jeremie/ST_RELEASE/MOSAIC/CUBE2VEL/cube2vel.so','cube2vel', $
            float(v), $
            imaginary(v), $
            nn[1], $
            nn[2], $
            nn[3], $
            float(cube[*].err), $
            imaginary(cube[*].err), $
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
            nn2,/unload)
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
