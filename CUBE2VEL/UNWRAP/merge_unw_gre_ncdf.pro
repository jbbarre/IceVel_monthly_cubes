pro merge_unw_gre_ncdf,display=display 

path='/mnt/pennell-z0/eric/GREENLAND/'

SENSORS = [$
    'alos2007_greenland',$
    'alos2009_greenland',$
    'envisat_greenland_2010',$
    'RADARSAT-2_2014_greenland',$
    'RADARSAT_greenland_2009',$
    'ERS1992_greenland']

geo0=load_geo_param()

cr = string("15b) ;"
FORm="($,i5,a,a)"

cd,current=pwd

;;size of the cubes
nrec = 250
npix = 250

unik=strcompress(RANDOMU(Seed, 1,/long),/r)

SHMMap, 'shm_vx_final'+unik,/float,Dimension=[geo0.npix,geo0.nrec]
SHMMap, 'shm_vy_final'+unik,/float,Dimension=[geo0.npix,geo0.nrec]
SHMMap, 'shm_stdx_final'+unik,/float,Dimension=[geo0.npix,geo0.nrec]
SHMMap, 'shm_stdy_final'+unik,/float,Dimension=[geo0.npix,geo0.nrec]
SHMMap, 'shm_cnt_final'+unik,/long,Dimension=[geo0.npix,geo0.nrec]

vx_final = SHMVAR('shm_vx_final'+unik)
vy_final = SHMVAR('shm_vy_final'+unik)
stdx_final = SHMVAR('shm_stdx_final'+unik)
stdy_final = SHMVAR('shm_stdy_final'+unik)
cnt_final = SHMVAR('shm_cnt_final'+unik)

;; initialisation mthread {{{
IF not keyword_set(nthreads) THEN nthreads=6
FOR i=0,nthreads-1 DO BEGIN
    bridge = (i eq 0)? create_struct('n'+string(i,FORmat='(I2)'), Obj_New('IDL_IDLBridge')) : $
        create_struct('n'+string(i,FORmat='(I2)'), Obj_New('IDL_IDLBridge'),bridge) ;;; add output='' FOR debug
ENDFOR

FOR ithread=0,nthreads-1 DO BEGIN
    print,'ithread init :',ithread
    ;; prepare thread {{{
    bridge.(ithread)->EXECUTE, '@' + PREF_GET('IDL_STARTUP')
    bridge.(ithread)->setvar,'pwd',pwd
    bridge.(ithread)->execute,'cd,pwd'
    ;; }}}

    struct_pass,geo0,bridge.(ithread)
    bridge.(ithread)->setvar,'npix',npix
    bridge.(ithread)->setvar,'nrec',nrec
    bridge.(ithread)->setvar,'sensors',sensors

    bridge.(ithread)->EXECUTE,"SHMMap, 'shm_vx_final"+unik+"',/float,Dimension=[geo0.npix,geo0.nrec]"
    bridge.(ithread)->EXECUTE,"SHMMap, 'shm_vy_final"+unik+"',/float,Dimension=[geo0.npix,geo0.nrec]"
    bridge.(ithread)->EXECUTE,"SHMMap, 'shm_stdx_final"+unik+"',/float,Dimension=[geo0.npix,geo0.nrec]"
    bridge.(ithread)->EXECUTE,"SHMMap, 'shm_stdy_final"+unik+"',/float,Dimension=[geo0.npix,geo0.nrec]"
    bridge.(ithread)->EXECUTE,"SHMMap, 'shm_cnt_final"+unik+"',/long,Dimension=[geo0.npix,geo0.nrec]"

    bridge.(ithread)->EXECUTE,"vx_final = SHMVAR('shm_vx_final"+unik+"')"
    bridge.(ithread)->EXECUTE,"vy_final = SHMVAR('shm_vy_final"+unik+"')"
    bridge.(ithread)->EXECUTE,"stdx_final = SHMVAR('shm_stdx_final"+unik+"')"
    bridge.(ithread)->EXECUTE,"stdy_final = SHMVAR('shm_stdy_final"+unik+"')"
    bridge.(ithread)->EXECUTE,"cnt_final = SHMVAR('shm_cnt_final"+unik+"')"


ENDFOR

ithread = 0
Bridges_RUNNING = bytarr(nthreads)

!Quiet=1
print
print,'MULTI-THREADS :',nthreads,' x=running, o=free'
out = ''
FOR i=0,nthreads-1 DO out += string(i,FORmat='(I02)')+' '
print,out
progress_mthread,Bridges_RUNNING
;; }}}

if keyword_set(display) then winjer,geo0.npix/25.,geo0.nrec/25.,title=pwd

;07595_y07840
FOR i=735.,3500,nrec-5 DO BEGIN;geo0.npix-1,nrec-5 DO BEGIN
        FOR j=2940.,3500,npix-5 DO BEGIN;geo0.nrec-1,npix-5 DO BEGIN

        percent = (float(i)*geo0.nrec+j*npix)/(float(geo0.npix)*float(geo0.nrec)) * 100.
        msg = strcompress(i,/r)+' '+strcompress(j,/r)
        ;; check free thread {{{
        IF ithread LT nthreads THEN IF Bridges_RUNNING[ithread] THEN GOTO,JUMP_FOR_FREE_THREAD
        IF ithread GE nthreads THEN BEGIN
            JUMP_FOR_FREE_THREAD:
            NOFREE_BRIDGE=1
            WHILE NOFREE_BRIDGE DO BEGIN

                FOR ith = 0, nthreads-1 DO BEGIN

                    IF bridge.(ith)->status() eq 2 or bridge.(ith)->status() eq 0 THEN BEGIN

                        bridge.(ith)->execute,'tmp = n_elements(result)'
                        tmp = bridge.(ith)->getvar('tmp')

                        IF tmp ne 0 then BEGIN
                            bridge.(ith)->execute,'vx_final[i,j] = result.vx[0:(npix-1L)<(geo0.npix-1L-i),0:(nrec-1L)<(geo0.nrec-1L-j)]'
                            bridge.(ith)->execute,'vy_final[i,j] = result.vy[0:(npix-1L)<(geo0.npix-1L-i),0:(nrec-1L)<(geo0.nrec-1L-j)]'
                            bridge.(ith)->execute,'stdx_final[i,j] = result.stdx[0:(npix-1L)<(geo0.npix-1L-i),0:(nrec-1L)<(geo0.nrec-1L-j)]'
                            bridge.(ith)->execute,'stdy_final[i,j] = result.stdy[0:(npix-1L)<(geo0.npix-1L-i),0:(nrec-1L)<(geo0.nrec-1L-j)]'
                            bridge.(ith)->execute,'cnt_final[i,j]  = result.cnt[0:(npix-1L)<(geo0.npix-1L-i),0:(nrec-1L)<(geo0.nrec-1L-j)]'
                            bridge.(ith)->execute,'result=0'
                            bridge.(ith)->execute,'delvar,result'
                            bridge.(ith)->execute,'delvar,i'
                            bridge.(ith)->execute,'delvar,j'

                        ENDIF

                        Bridges_RUNNING[ith]=0
                        ithread=min(where(Bridges_RUNNING eq 0))
                        NOFREE_BRIDGE=0

                        ;IF NOFREE_BRIDGE eq 0 THEN GOTO,JUMP_ENDWHILE
                    ENDIF
                    WAIT,0.05
                ENDFOR
            ENDWHILE

        ENDIF
        ;; }}}
        bridge.(ithread)->setvar,'i',i
        bridge.(ithread)->setvar,'j',j
        bridge.(ithread)->execute,'result=unw2vel_gre_feather_sp_fortran(i,j,geo0,sensors=sensors)',/NOWAIT
        ;result=unw2vel_gre_feather_sp_fortran(i,j,geo0,sensors=sensors,/verb)
        ;stop
        Bridges_RUNNING[ithread]=1
        ithread += 1

        progress_mthread,Bridges_RUNNING,percent=percent,msg=msg
        IF j eq 0 and keyword_set(display) THEN BEGIN
            loadct,15,/silent
            tvscl,(abs(complex(congrid( vx_final, geo0.npix/25., geo0.nrec/25.), congrid( vy_final, geo0.npix/25., geo0.nrec/25.)))<3000)^0.35 ,/nan
        ENDIF
    ENDFOR
ENDFOR

;; write NC {{{

long_name = ['Ice velocity in x direction',$
    'Ice velocity in y direction',$
    'Std velocity in x direction',$
    'Std velocity in y direction',$
    'Cnt velocity']

standard_name = ['land_ice_x_velocity',$
    'land_ice_y_velocity',$
    'std_ice_x_velocity',$
    'std_ice_y_velocity',$
    'cnt_ice_velocity']

units = ['meter/year','meter/year','meter/year','meter/year','meter/year','meter/year','']

save_ncdf_ps_south,'vel_final.nc',geo0,$
    vx=vx_final,$
    vy=vy_final,$
    stdx=stdx_final,$
    stdy=stdy_final,$
    cnt=nn_final, long_name=long_name,$
    standard_name=standard_name,units=units

 write_jpeg,'unw_final.jpg',bytscl((abs(congrid(complex(vx_final,vy_final),1244,1244))<3000)^0.35)

stop
end
