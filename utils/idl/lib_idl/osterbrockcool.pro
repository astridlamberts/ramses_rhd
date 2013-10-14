function osterbrockcool, t

h=6.63d-27
c=3.0d10
kb=1.38d-16
angstrom=1d-8
micron=1d-4

; Ions included are oxygen, neon, and nitrogen, in the singly and 
; doubly ionized states. Use solar abundances, with 80% of gas assumed
; to be in singly ionized state, 20% in doubly ionized
xO=7d-4
xNe=9d-5
xN=9d-5
xSingle=0.8
xDouble=1.0-xSingle

lc=0

; We include all transitions from the ground states of these ions.
; Statistical weigth: a state is described by (2S+1) L or (2S+1) L J,
; where L = 0, 1, 2, 3... corresponds to S, P, D, F, ... The weight
; is (2S+1)(2L+1) for a state without a J number, or (2J+1) for a state
; with a J number

; singly ionized oxygen, 4S -> 2D
omega=1.34
lambda=3727.*angstrom
wgt1=4
energy=h*c/lambda
q12=8.629d-6/sqrt(t)*omega/wgt1*exp(-energy/(kb*t))
cool=xO*xSingle*q12*energy
lc=lc+cool

; singly ionized oxygen, 4S -> 2P
omega=0.4
lambda=2470.*angstrom
wgt1=4
energy=h*c/lambda
q12=8.629d-6/sqrt(t)*omega/wgt1*exp(-energy/(kb*t))
cool=xO*xSingle*q12*energy
lc=lc+cool

; doubly ionized oxygen, 3P -> 1D
omega=2.29
lambda=5000.*angstrom
wgt1=9
energy=h*c/lambda
q12=8.629d-6/sqrt(t)*omega/wgt1*exp(-energy/(kb*t))
cool=xO*xDouble*q12*energy
lc=lc+cool

; doubly ionized oxygen, 3P -> 1S
omega=0.29
lambda=2325.*angstrom
energy=h*c/lambda
wgt1=9
q12=8.629d-6/sqrt(t)*omega/wgt1*exp(-energy/(kb*t))
cool=xO*xDouble*q12*energy
lc=lc+cool

; doubly ionized oxygen, 3P0 -> 3P1
omega=0.55
lambda=88.356*micron
energy=h*c/lambda
wgt1=1
q12=8.629d-6/sqrt(t)*omega/wgt1*exp(-energy/(kb*t))
cool=xO*xDouble*q12*energy
lc=lc+cool

; doubly ionized oxygen, 3P0 -> 3P2
omega=0.27
lambda=32.661*micron
energy=h*c/lambda
wgt1=1
q12=8.629d-6/sqrt(t)*omega/wgt1*exp(-energy/(kb*t))
cool=xO*xDouble*q12*energy
lc=lc+cool

; singly ionized neon, 1P(1/2) -> 1P(3/2)
omega=0.28
lambda=12.814*micron
energy=h*c/lambda
wgt=2
q12=8.629d-6/sqrt(t)*omega/wgt1*exp(-energy/(kb*t))
cool=xNe*xSingle*q12*energy
lc=lc+cool

; doubly ionized neon, 3P -> 1D
omega=1.36
lambda=3950.*angstrom
energy=h*c/lambda
wgt1=9
q12=8.629d-6/sqrt(t)*omega/wgt1*exp(-energy/(kb*t))
cool=xNe*xDouble*q12*energy
lc=lc+cool

; doubly ionized neon, 3P -> 1S
omega=0.15
lambda=1800.*angstrom
energy=h*c/lambda
wgt1=9
q12=8.629d-6/sqrt(t)*omega/wgt1*exp(-energy/(kb*t))
cool=xNe*xDouble*q12*energy
lc=lc+cool

; doubly ionized neon, 3P0 -> 3P1
omega=0.24
lambda=36.02*micron
energy=h*c/lambda
wgt1=1
q12=8.629d-6/sqrt(t)*omega/wgt1*exp(-energy/(kb*t))
cool=xNe*xDouble*q12*energy
lc=lc+cool

; doubly ionized neon, 3P0 -> 3P2
omega=0.21
lambda=10.86*micron
energy=h*c/lambda
wgt1=1
q12=8.629d-6/sqrt(t)*omega/wgt1*exp(-energy/(kb*t))
cool=xNe*xDouble*q12*energy
lc=lc+cool

; singly ionized nitrogen, 3P -> 1D
omega=2.64
lambda=6550.*angstrom
energy=h*c/lambda
wgt1=9
q12=8.629d-6/sqrt(t)*omega/wgt1*exp(-energy/(kb*t))
cool=xN*xSingle*q12*energy
lc=lc+cool

; singly ionized nitrogen, 3P -> 1S
omega=0.29
lambda=3065.*angstrom
wgt1=9
q12=8.629d-6/sqrt(t)*omega/wgt1*exp(-energy/(kb*t))
cool=xN*xSingle*q12*energy
lc=lc+cool

; singly ionized nitrogen, 3P0 -> 3P1
omega=0.41
lambda=205.5*micron
energy=h*c/lambda
wgt1=1
q12=8.629d-6/sqrt(t)*omega/wgt1*exp(-energy/(kb*t))
cool=xN*xSingle*q12*energy
lc=lc+cool

; singly ionized nitrogen, 3P0 -> 3P2
omega=0.27
lambda=76.5*micron
energy=h*c/lambda
wgt1=1
q12=8.629d-6/sqrt(t)*omega/wgt1*exp(-energy/(kb*t))
cool=xN*xSingle*q12*energy
lc=lc+cool

; doubly ionized nitrogen, 1P(1/2)->1P(3/2)
omega=1.45
lambda=57.343*micron
energy=h*c/lambda
wgt1=2
q12=8.629d-6/sqrt(t)*omega/wgt1*exp(-energy/(kb*t))
cool=xN*xDouble*q12*energy
lc=lc+cool

; multiply by common factor
lc = lc*8.63d-6/sqrt(t)

; Also include hydrogen free-free cooling
cool=1.42d-27*sqrt(t)*1.3
lc=lc+cool

return, lc

end
