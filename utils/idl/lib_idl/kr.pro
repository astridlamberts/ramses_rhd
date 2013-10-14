function kr, t

tcut=[0., 350., 600., 700., 950., 1000.]
slope=[4.4, 0.0, 0.0, 0.0, -0.15, 0.0]
base=[0.1, 3.9, 0.7, 0.25, 0.25, 0.1]

nint=n_elements(tcut)
kr=t
for n=0,nint-1 do begin
    if (n eq nint-1) then begin
        idx=where(t ge tcut[n])
        if idx[0] ne -1 then kr[idx] = base[n]
    endif else begin
        idx=where((t ge tcut[n]) and (t lt tcut[n+1]))
        if idx[0] ne -1 then $
          kr[idx] = base[n] + $
                    slope[n]*(t[idx] - tcut[n])/(tcut[n+1]-tcut[n])
    endelse
endfor
return, kr

end
