function corrfunction, mat1, mat2, auto=auto

mat1ft=fft(mat1, -1)
if keyword_set(auto) then mat2ft=conj(mat1ft) $
else mat2ft=conj(fft(mat2, -1))
prod=mat1ft*mat2ft

end
