function pdfth, x, mach, int=int

sigma=sqrt(alog(1+mach^2/4.0))
mn=-sigma^2/2

if not keyword_set(int) then $
    return, 1/sqrt(2*!pi*sigma^2)*exp(-0.5*((x-mn)/sigma)^2) $
else $
    return, 1/sqrt(2*!pi*sigma^2)*exp(-0.5*((x-mn)/sigma)^2)*(x[1]-x[0])

end

