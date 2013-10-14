function orb_param, m1, m2, r1, r2, v1, v2

; Procedure to compute the orbital parameters of a system of two
; objects with the given mass, position, and velocity. The routine
; returns a, the semi-major axis, and e, the eccentricity. All
; units are assumed to be CGS.

; note for reference: for vectors a and b, transpose(a)#b == a . b

G = 6.67d-8

; compute the reduced mass
mu = m1*m2/(m1+m2)

; compute the separation vector and its derivative in Cartesian coordinates
rvec = r1 - r2
rvecdot = v1 - v2

; compute the separation between the particles
r = sqrt(transpose(rvec)#rvec)

; compute the rate at which the separation of the two particles is changing,
; given by r.rdot/r
rdot = transpose(rvec)#rvecdot/r

; compute the velocity in the cm frame, (r1 - r2)^2
v = sqrt(transpose(rvecdot)#rvecdot)

; compute the rotation rate in the cm frame, using
; v^2 = rdot^2 + r^2 thetadot^2
thetadot = sqrt(v^2-rdot^2)/r

; compute the force constant
k = G*m1*m2

; compute the cm frame angular momentum
l = mu*r^2*thetadot

; compute the cm frame energy
en = 0.5*mu*v^2 - k/r

; compute the semi-major axis
a = k/(2*abs(en))

; compute the eccentricity
e = sqrt(1 + 2*en*l^2/(mu*k^2))

; return
return, [a, e]

end




