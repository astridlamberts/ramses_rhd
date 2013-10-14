function ucbplnoise,npts,alpha
;+
; NAME:
;      ucbplnoise
;
;
; PURPOSE:
;     Produce 1-d "timeline" of power law noise
;
;
; CATEGORY:
;     math, statistics, probability, fft, simulation
;
;
; CALLING SEQUENCE:
;     ucbplnoise,npoints,logrhythmic_index
;
;
; INPUTS:
;     npoints = number of sample points to be produced
;     logrhythmic_index = index of power law fit to fft of output
;
; OPTIONAL INPUTS:
; KEYWORD PARAMETERS:
; OUTPUTS:
; OPTIONAL OUTPUTS:
; RESTRICTIONS:
;
;
;
; PROCEDURE:
;    Procedure by Mark Kruholz of UCB
;    Just create an array of the appropriate number of (complex)
;points in k-space. Assign each point at positive k in the array a Gaussian
;random amplitude with a dispersion of unity (I'm sure there's an IDL
;routine to do that), then multiply by factor k^alpha, where alpha is
;chosen for whatever power law you want to have. Also assign each point at
;positive k a random phase from 0 to 2 pi. Then set all the negative k
;points equal to the complex conjugates of the corresponding positive k
;points (to force the results in real space to be purely real). Then just
;inverse Fourier transform (using the FFT routine) the array and set the
;normalization to whatever you want. The result is an array of Gaussian
;random noise with a power spectrum k^alpha.

;
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;      Coded by Bruce Grossan 17Apr03.  Copyright Bruce Grossan
;      2003-2050.  No unauthorized use permitted without donating
;      goodies to either Grossan or a UCB astro tea. All math stolen
;      from Mark.  Doesn't appear to work, though.
;
;-

; Create an complex array in K-space
transform=complexarr(npts)

; with gausian random values
transform=complex(randomn(seed,n_elements(transform)),imaginary(transform)) 

; multiply by k to the alpha-th power
transform[1:(npts+1)/2]=transform[1:(npts+1)/2] * $
        (findgen((npts+1)/2)+1)^alpha
transform[0]=0.0

; multiply by a complex random phase
phase=randomu(seed,npts)*2*!pi
stuff=complex(0,phase)
transform=exp(stuff)*transform

;adjust for IDL style fft spectra
if npts mod 2 eq 0 then begin
   transform(npts/2+1:npts-1)=$
      conj(transform[npts/2-1-indgen(npts/2-1)])
      transform(npts/2)=abs(transform(npts/2)) ; nyquist term must be real 
endif else $
   transform(npts/2+1:npts-1)=$
      conj(transform[npts/2-indgen(npts/2)])



transform(0)=complex(abs(transform(0))) ;because we want const term real

tl=fft(transform,/inv)
return,tl
end

