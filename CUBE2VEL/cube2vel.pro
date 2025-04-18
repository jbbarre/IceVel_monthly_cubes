function cube2vel,cube0
    cube=cube0

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
    
    velo = complex(median(float(cube[*].v),dim=3,/even), median(imaginary(cube[*].v),dim=3,/even))

    nn = total(finite(abs(cube[*].v)),3)
    std = complex(stddev(float(cube[*].v),dim=3,/nan), stddev(imaginary(cube[*].v),dim=3,/nan))
    
    ;; eliminate value far from median solution ; keep IF only one or two values {{{
    tmp = float(cube[*].v)
    FOR i=0,n_elements(cube)-1 DO tmp[*,*,i] = 1. /  ( ( (abs(float(cube[i].v) - float(velo)) LE (2.*float(std) )) AND $
        (abs(imaginary(cube[i].v) - imaginary(velo)) LE (2.*imaginary(std)) ) ) OR (nn EQ 1) OR (nn EQ 0) OR (nn EQ 2) )
    
    w=where(finite(tmp) eq 0,cnt)
    IF cnt ne 0 THEN tmp[w]=!values.f_nan
    w=0
    
    cube[*].v = cube[*].v*tmp
    tmp=0
    ;; }}}

    ;; weighted mean 

    weight = cube[*].v ;; == weight 
    FOR i=0,n_elements(cube)-1 DO $
        weight[*,*,i] = finite(cube[i].v) *  complex(1./float(cube[i].err)^2,1./imaginary(cube[i].err)^2)
   
    w_mean = total(cube[*].v*weight,3,/nan)/total(weight,3) ;; weighted mean
    err0 = cube[*].v
    FOR i=0,n_elements(cube)-1 DO $
        err0[*,*,i] = finite(cube[i].v) * cube[i].err
    
    errx = sqrt( total( (float(weight)*float(err0))^2 ,3,/nan) / total( float(weight) ,3)^2 )
    erry = sqrt( total( (imaginary(weight)*imaginary(err0))^2 ,3,/nan) / total( imaginary(weight) ,3)^2 )
    err = complex(errx,erry)
    err0 = 0
    ;; weighted standard deviation
    nn2 = total(finite(abs(cube[*].v)),3)
    
    std_temp = complex( cmreplicate(float(w_mean),n_elements(cube[*].v[0,0])) , cmreplicate(imaginary(w_mean),n_elements(cube[*].v[0,0])) )
    std2 = sqrt( nn2 * total( weight*(cube[*].v - std_temp)^2,3,/nan) / total(weight,3) / (nn2-1) )
    std_temp=0
    weight=0
    return,{weighted_mean :w_mean, $
            weighted_std  :std2,$
            n_filt_values :nn2,$
            med           :velo,$
            std           :std,$
            err           :err,$
            n_values      :nn}
end
