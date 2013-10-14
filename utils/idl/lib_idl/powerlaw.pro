pro powerlaw, x, a, f, pder

pder=fltarr(n_elements(x),2)
pder[*,0]=x^(a[1])
pder[*,1]=a[0]*a[1]*x^(a[1]-1)

f=a[0]*x^(a[1])

end