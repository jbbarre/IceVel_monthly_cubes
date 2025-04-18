
geo=load_geo_param()

a=fltarr(2,2)
b=fltarr(2,2)
c=fltarr(2,2)
unity=a & unity[*,0]=[1.,0.] & unity[*,1]=[0.,1.]
cr = string("15b) ;"
FORm="($,i5,a,a)"

vel_final = complexarr(geo.npix,geo.nrec)
dvra = 5d-2 ;; assuming 5 centimeters precision FOR now
dvrb = 5d-2 
dvrc = 5d-2

FOR i=8910.,geo.npix-1,990. DO BEGIN
    FOR j=4950.,geo.nrec-1,990. DO BEGIN

        xmin = float(i)*geo.posting + geo.xmin
        xmax = xmin + (999.)*geo.posting
        ymax = geo.ymax - float(j)*geo.posting
        ymin = ymax - 999.*geo.posting

        geoloc = {xmin:xmin,ymax:ymax,posting:geo.posting,nrec:1000,npix:1000}

        print,'==========================================================================================='
        print
        print,'unw_cubes/unw_cube_'+strcompress(long(xmin/1d3),/r)+'-'+strcompress(long(ymax/1d3),/r)+'.sav'
        print

        ;; load list of cubes {{{
        list_of_cubes = []
        
        IF file_test('../../ALOS_2007/MOSAIC/unw_cubes/unw_cube_'+strcompress(long(xmin/1d3),/r)+'-'+strcompress(long(ymax/1d3),/r)+'.sav') THEN $
            list_of_cubes = [list_of_cubes,'../../ALOS_2007/MOSAIC/unw_cubes/unw_cube_'+strcompress(long(xmin/1d3),/r)+'-'+strcompress(long(ymax/1d3),/r)+'.sav']
        IF file_test('../../ERS1996/MOSAIC/unw_cubes/unw_cube_'+strcompress(long(xmin/1d3),/r)+'-'+strcompress(long(ymax/1d3),/r)+'.sav') THEN $
            list_of_cubes = [list_of_cubes,'../../ERS1996/MOSAIC/unw_cubes/unw_cube_'+strcompress(long(xmin/1d3),/r)+'-'+strcompress(long(ymax/1d3),/r)+'.sav']
        IF file_test('../../ALOS_2008/MOSAIC/unw_cubes/unw_cube_'+strcompress(long(xmin/1d3),/r)+'-'+strcompress(long(ymax/1d3),/r)+'.sav') THEN $
            list_of_cubes = [list_of_cubes,'../../ALOS_2008/MOSAIC/unw_cubes/unw_cube_'+strcompress(long(xmin/1d3),/r)+'-'+strcompress(long(ymax/1d3),/r)+'.sav']
        IF file_test('../../RADARSAT1997/MOSAIC/unw_cubes/unw_cube_'+strcompress(long(xmin/1d3),/r)+'-'+strcompress(long(ymax/1d3),/r)+'.sav') THEN $
            list_of_cubes = [list_of_cubes,'../../RADARSAT1997/MOSAIC/unw_cubes/unw_cube_'+strcompress(long(xmin/1d3),/r)+'-'+strcompress(long(ymax/1d3),/r)+'.sav']
        IF file_test('../../envisat2007/MOSAIC/unw_cubes/unw_cube_'+strcompress(long(xmin/1d3),/r)+'-'+strcompress(long(ymax/1d3),/r)+'.sav') THEN $
            list_of_cubes = [list_of_cubes,'../../envisat2007/MOSAIC/unw_cubes/unw_cube_'+strcompress(long(xmin/1d3),/r)+'-'+strcompress(long(ymax/1d3),/r)+'.sav']
        IF file_test('../../envisat2008/MOSAIC/unw_cubes/unw_cube_'+strcompress(long(xmin/1d3),/r)+'-'+strcompress(long(ymax/1d3),/r)+'.sav') THEN $
            list_of_cubes = [list_of_cubes,'../../envisat2008/MOSAIC/unw_cubes/unw_cube_'+strcompress(long(xmin/1d3),/r)+'-'+strcompress(long(ymax/1d3),/r)+'.sav']
        IF file_test('../../envisat2009/MOSAIC/unw_cubes/unw_cube_'+strcompress(long(xmin/1d3),/r)+'-'+strcompress(long(ymax/1d3),/r)+'.sav') THEN $
            list_of_cubes = [list_of_cubes,'../../envisat2009/MOSAIC/unw_cubes/unw_cube_'+strcompress(long(xmin/1d3),/r)+'-'+strcompress(long(ymax/1d3),/r)+'.sav']
        IF file_test('../../RADARSAT-2_2015/MOSAIC/unw_cubes/unw_cube_'+strcompress(long(xmin/1d3),/r)+'-'+strcompress(long(ymax/1d3),/r)+'.sav') THEN $            
            list_of_cubes = [list_of_cubes,'../../RADARSAT-2_2015/MOSAIC/unw_cubes/unw_cube_'+strcompress(long(xmin/1d3),/r)+'-'+strcompress(long(ymax/1d3),/r)+'.sav']
        IF file_test('../../RADARSAT-2_2013/MOSAIC/unw_cubes/unw_cube_'+strcompress(long(xmin/1d3),/r)+'-'+strcompress(long(ymax/1d3),/r)+'.sav') THEN $                    
            list_of_cubes = [list_of_cubes,'../../RADARSAT-2_2013/MOSAIC/unw_cubes/unw_cube_'+strcompress(long(xmin/1d3),/r)+'-'+strcompress(long(ymax/1d3),/r)+'.sav']
        IF file_test('../../RADARSAT-2_2011/MOSAIC/unw_cubes/unw_cube_'+strcompress(long(xmin/1d3),/r)+'-'+strcompress(long(ymax/1d3),/r)+'.sav') THEN $
            list_of_cubes = [list_of_cubes,'../../RADARSAT-2_2011/MOSAIC/unw_cubes/unw_cube_'+strcompress(long(xmin/1d3),/r)+'-'+strcompress(long(ymax/1d3),/r)+'.sav']
        IF file_test('../../RADARSAT-2_2009/MOSAIC/unw_cubes/unw_cube_'+strcompress(long(xmin/1d3),/r)+'-'+strcompress(long(ymax/1d3),/r)+'.sav') THEN $
            list_of_cubes = [list_of_cubes,'../../RADARSAT-2_2009/MOSAIC/unw_cubes/unw_cube_'+strcompress(long(xmin/1d3),/r)+'-'+strcompress(long(ymax/1d3),/r)+'.sav']
        ;; }}}
        
        IF n_elements(list_of_cubes) ne 0 THEN BEGIN ; {{{

            cube0 = []
            FOR il=0,n_elements(list_of_cubes)-1 DO BEGIN
                restore,list_of_cubes[il]
                cube0 = [cube0,cube]
            ENDFOR

            cube = cube0
            cube0=0

            IF n_elements(ll) ge 1 THEN BEGIN

                v_final = fltarr(3,1000,1000)
                v_final_med = fltarr(3,1000,1000)


                ;; loop over x ,y axis {{{ 
                FOR ii=0,999 DO BEGIN

                    progress_bar,ii/999.*100. 
                    FOR jj=0,999 DO BEGIN
                        ;print,jj
                        ;IF (jj mod 100) eq 0 THEN print,FORm=FORm,jj/999.*100.,' %',cr
                        VX = []
                        VY = []
                        VZ = []

                        DVX = []
                        DVY = []
                        DVZ = []

                        w=where(finite(cube[*].v[ii,jj]) ne 0,cnt)
                        IF cnt ne 0 THEN BEGIN
                        head = 2.*!pi - reFORm(cube[w].heading[ii,jj])*!dtor
                        inc = cube[w].inc[ii,jj]
                        v = cube[w].v[ii,jj]

                        ;; loop over phases
                        FOR kk=0,cnt-3 DO BEGIN
                            vra = v[kk]
                            Xa = head[kk]
                            THETAa = inc[kk]

                            FOR mm=kk+1,cnt-2 DO BEGIN

                                vrb = v[mm]
                                Xb = head[mm]
                                THETAb = inc[mm]
                                IF (ABS(min( [360.*!dtor - abs(Xa - Xb),abs(Xa - Xb) ])) GE 7*!dtor) THEN BEGIN

                                    ;FOR oo=mm+1,cnt-1 DO BEGIN
                                    oo=mm+1
                                    vrc = v[oo:*]
                                    Xc = head[oo:*]
                                    THETAc = inc[oo:*]

                                    woo = where( (ABS(min( [360.*!dtor - abs(Xa - Xc),abs(Xa - Xc) ])) GE 7*!dtor) AND $
                                        (ABS(min( [360.*!dtor - abs(Xb - Xc),abs(Xb - Xc) ])) GE 7*!dtor),cnt_oo)

                                    IF cnt_oo ne 0 THEN BEGIN

                                        vrc=vrc[woo]
                                        xc=xc[woo]
                                        THETAc=THETAc[woo]
                                        ;; angle > 30 degrees -> combine phase {{{                                        
                                        ;IF (ABS(min( [360.*!dtor - abs(Xa - Xc),abs(Xa - Xc) ])) GE 10*!dtor) AND $
                                        ;    (ABS(min( [360.*!dtor - abs(Xb - Xc),abs(Xb - Xc) ])) GE 10*!dtor) THEN BEGIN

                                        ;    THETAc = inc[oo]

                                        invT1 = 1. / (sin(Xa-Xc)/tan(THETAb) + sin(Xc-Xb)/tan(THETAa) + sin(Xb-Xa)/tan(THETAc))

                                        vz = [ vz , invT1 * ( sin(Xc-Xa)/sin(THETAb)*vrb + $
                                            sin(Xa-Xb)/sin(THETAc)*vrc + $
                                            sin(Xb-Xc)/sin(THETAa)*vra)]

                                        invdetT = 1./ (sin(THETAa)*sin(THETAb)*sin(THETAc)) * invT1

                                        A =  sin(Xb)*sin(THETAb)*cos(THETAc) - cos(THETAb)*sin(Xc)*sin(THETAc)
                                        B = -cos(THETAb)*cos(Xc)*sin(THETAc) + cos(Xb)*sin(THETAb)*cos(THETAc)
                                        D =  cos(THETAa)*sin(Xc)*sin(THETAc) - sin(Xa)*sin(THETAa)*cos(THETAc)
                                        E = -cos(Xa)*sin(THETAa)*cos(THETAc) + cos(THETAa)*cos(Xc)*sin(THETAc)
                                        G =  sin(Xa)*sin(THETAa)*cos(THETAb) - cos(THETAa)*sin(Xb)*sin(THETAb)
                                        H = -cos(THETAa)*cos(Xb)*sin(THETAb) + cos(Xa)*sin(THETAa)*cos(THETAb) ;; K => H error in Gray et al. 2011 sup. materials

                                        vx = [vx, invdetT * (A*vra + D*vrb + G*vrc) ]
                                        vy = [vy ,invdetT * (B*vra + E*vrb + H*vrc) ]

                                        dvx = [ dvx , abs(invdetT) * sqrt( (A*dvra)^2 + (D*dvrb)^2 + (G*dvrc)^2 ) ]
                                        dvy = [ dvy , abs(invdetT) * sqrt( (B*dvra)^2 + (E*dvrb)^2 + (H*dvrc)^2 ) ]
                                        dvz = [ dvz , abs(invT1) * sqrt( (sin(Xc-Xa)/sin(THETAb)*dvrb)^2 + $
                                            (sin(Xa-Xb)/sin(THETAc)*vrc)^2 + $
                                            (sin(Xb-Xc)/sin(THETAa)*vra)^2) ]
                                    ENDIF
                                    ;; }}}
                                    ;ENDFOR
                                ENDIF
                            ENDFOR
                        ENDFOR

                        IF n_elements(vx) ge 2 THEN $
                            ;; eliminate outliers
                            V_final[*,ii,jj] = [total(vx*dvx^2)/total(dvx^2),total(vy*dvy^2)/total(dvy^2),total(vz*dvz^2)/total(dvz^2)]
                        IF n_elements(vx) ge 3 THEN V_final_med[*,ii,jj] = [median(vx),median(vy),median(vz)]
                        ;IF n_elements(vx) eq 2 THEN V_final_med[*,ii,jj] = [mean(vx),mean(vy),mean(vz)]
                        IF n_elements(vx) eq 1 THEN V_final[*,ii,jj] = [vx,vy,vz]
                    ENDIF
                    ENDFOR
                ENDFOR ;; }}}
                stop
                print,total(float(v_final))
                IF total(float(v_final)) ne 0 THEN vel_final[i,j] = v_final[0:999<(geo.npix-1L-i),0:999<(geo.nrec-1L-j)]

            ENDIF
        ENDIF ;; }}}
        stop
    ENDFOR
ENDFOR

save,filename='unw_final.sav',vel_final

stop
end
