function rhoetemp, rho, vx, vy, vz, e, gamma=gamma, mu=mu, ke=ke, saveke=saveke

if not keyword_set(mu) then mu=14./6.*1.67d-24
if not keyword_set(gamma) then gamma=5./3.
kb=1.38d-16

vx2=amr_multiply(vx, vz)
vy2=amr_multiply(vy, vy)
vz2=amr_multiply(vz, vz)
v2tmp=amr_add(vx2, vy2)
v2=amr_add(v2tmp, vz2)
amr_free, vx2
amr_free, vy2
amr_free, vz2
amr_free, v2tmp
twoke=amr_multiply(rho, v2)
ke=amr_divide(twoke, 2.0)
amr_free, twoke
ie=amr_subtract(e, ke)
if not keyword_set(saveke) then amr_free, ke

esp=amr_divide(ie, rho)
temp=amr_multiply(esp, (gamma-1)*mu/kb)

return, temp

end
