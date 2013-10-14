function planck, t, nu

h=6.67d-27
c=3d10
kb=1.38d-16

return, 2*h*nu^3/c^2/(exp(h*nu/(kb*t))-1)

end
