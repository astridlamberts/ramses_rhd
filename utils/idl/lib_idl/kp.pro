function kp, t

tcut=[0., 375., 575., 675., 960., 1100.]
slope=[7.0, 0.7, 0.1, 0.3, -3.0, 0.0]
base=[0.3, 7.3, 3.0, 2.8, 3.1, 0.1]

nint=n_elements(tcut)
kp=t
for n=0,nint-1 do begin
    if (n eq nint-1) then begin
        idx=where(t ge tcut[n])
        if idx[0] ne -1 then kp[idx] = base[n]
    endif else begin
        idx=where((t ge tcut[n]) and (t lt tcut[n+1]))
        if idx[0] ne -1 then $
          kp[idx] = base[n] + $
                    slope[n]*(t[idx] - tcut[n])/(tcut[n+1]-tcut[n])
    endelse
endfor
return, kp

end
