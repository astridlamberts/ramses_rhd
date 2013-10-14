function alphabondiinterp, x

; Note: physical density is alpha(x)*rhoinf, velocity is
; lambda/(x^2 alpha(x))

alphatab=[820.254, 701.882, 600.752, 514.341, 440.497, 377.381, $
323.427, 277.295, 237.845, 204.1, 175.23, 150.524, 129.377, 111.27, $
95.7613, 82.4745, 71.0869, 61.3237, 52.9498, 45.7644, 39.5963, $
34.2989, 29.7471, 25.8338, 22.4676, 19.5705, 17.0755, 14.9254, $
13.0714, 11.4717, 10.0903, 8.89675, 7.86467, 6.97159, 6.19825, $
5.52812, 4.94699, 4.44279, 4.00497, 3.6246, 3.29395, 3.00637, $
2.75612,  2.53827, 2.34854, 2.18322, 2.03912, 1.91344, 1.80378, $
1.70804, 1.62439] 


lambda = 0.25*exp(1.5)

x1 = abs(x)

itab = floor(50*alog(x1/0.01)/alog(2./0.01))
itab1 = itab * ((itab ge 0) and (itab lt 50))

wgt = (alog(x1) - alog(0.01) - itab/50.0*alog(2./0.01)) / $
		( 1.0/50.0 * alog(2./0.01) )

res = dblarr(n_elements(x))
res = res + (itab lt 0) * lambda/sqrt(2.0*x1^3)
sub=where(itab ge 50)
if (sub[0] ne -1) then res = res + (itab ge 50) * exp(1.0/x1)
res = res + ((itab ge 0) and (itab lt 50)) * $
	(alphatab[itab1]*(1.0-wgt) + alphatab[itab1+1]*wgt)

return, res

end
