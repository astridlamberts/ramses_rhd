function brighttemp_newt, t

common brigttemp_common, freq, i

h=6.67d-27
c=3d10

return, (planck(t, freq)-i)/(2*h*freq^3/c^2)

end


function brighttemp, i_input, nu, tguess=tguess

; return brightness temperature for intensity i observed at frequency nu

common brigttemp_common, freq, i

freq=nu
i=i_input
ni=n_elements(i)
if not keyword_set(tguess) then t=replicate(10, ni) $
else if n_elements(tguess) eq 1 then t=replicate(tguess, ni) $
else t=tguess

return, newton(t, 'brighttemp_newt')

end

