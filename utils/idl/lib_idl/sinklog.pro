; get wc of files
spawn, 'wc SinkLog', res
res=strtrim(res, 1)
nline=long(res)

time=0.0
npart=0L
openr, fp, 'SinkLog', /get_lun
readf, fp, time, npart
point_lun, fp, 0
data=dblarr(2+10*npart,nline)
readf, fp, data
free_lun, fp

time=data[0,*]
if (npart eq 1) then begin
	msink=data[2,*]
endif else begin
	msink=dblarr(npart,nline)
	for n=0,npart-1 do msink[n,*]=data[2+10*n,*]
endelse

; cut out non-unique parts
;time=make_monotonic(time, msink)

;plot, time, msink-msink[0]

;lfit=linfit(time, msink)

g=6.67d-8
cs=1.881d4
rhoinf=1.0d-25
vinf=0.0
msun=2.0d33
yrsec=365.25*24.*60.*60.

n=n_elements(msink)
mdot=(msink[1:n-1]-msink[0:n-2])/(time[1:n-1]-time[0:n-2])
thalf=0.5*(time[1:n-1]+time[0:n-2])

;mdotth=exp(1.5)*!pi*rhoinf*(g*msun)^2/(cs^3)
;lambda = exp(1.5)/4
;mdotth = 4*!pi*rhoinf*g^2*msink[0]^2* $
;	sqrt((lambda^2*cs^2+vinf^2)/(cs^2+vinf^2)^4)

;print, 'Mdot = ', lfit[1]*yrsec/msun
;print, 'Mdot theor = ', mdotth*yrsec/msun
;print, 'Mdot err = ', (lfit[1]-mdotth)/mdotth

end
