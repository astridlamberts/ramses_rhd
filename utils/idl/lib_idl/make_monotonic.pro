function make_monotonic, lst, l1, l2, l3, l4, l5, l6, l7, l8

nlst=n_elements(lst)
copylist=lonarr(nlst)-1
lastval=lst[0]
listptr=0L
copyptr=0L
newblock=0
for i=1L, nlst-1 do begin
    if (lst[i] gt lastval) then begin
	if (newblock) then begin
	    listptr=i
	    newblock=0
	endif
	lastval = lst[i]
	continue
    endif
    if (newblock eq 0) then begin
	newblock=1
	copylist[copyptr] = listptr
	copylist[copyptr+1] = i-1
	copyptr = copyptr+2
    endif
endfor
if (lastval eq lst[nlst-1]) then begin
    copylist[copyptr] = listptr
    copylist[copyptr+1] = i-1
    copyptr = copyptr+2
endif

outsize=0L
for i=0L,copyptr/2-1 do begin
   outsize=outsize+copylist[2*i+1]-copylist[2*i]+1
endfor

out=lst[0:outsize-1]
if (n_elements(l1) ne 0) then l1out=l1[0:outsize-1]
if (n_elements(l2) ne 0) then l2out=l2[0:outsize-1]
if (n_elements(l3) ne 0) then l3out=l3[0:outsize-1]
if (n_elements(l4) ne 0) then l4out=l4[0:outsize-1]
if (n_elements(l5) ne 0) then l5out=l5[0:outsize-1]
if (n_elements(l6) ne 0) then l6out=l6[0:outsize-1]
if (n_elements(l7) ne 0) then l7out=l7[0:outsize-1]
if (n_elements(l8) ne 0) then l8out=l8[0:outsize-1]

copyptr1=0
for i=0L,copyptr/2-1 do begin
    blocksize=copylist[2*i+1]-copylist[2*i]+1
    out[copyptr1:copyptr1+blocksize-1] = lst[copylist[2*i]:copylist[2*i+1]]
    if (n_elements(l1) ne 0) then $
	l1out[copyptr1:copyptr1+blocksize-1]=l1[copylist[2*i]:copylist[2*i+1]]
    if (n_elements(l2) ne 0) then $
	l2out[copyptr1:copyptr1+blocksize-1]=l2[copylist[2*i]:copylist[2*i+1]]
    if (n_elements(l3) ne 0) then $
	l3out[copyptr1:copyptr1+blocksize-1]=l3[copylist[2*i]:copylist[2*i+1]]
    if (n_elements(l4) ne 0) then $
	l4out[copyptr1:copyptr1+blocksize-1]=l4[copylist[2*i]:copylist[2*i+1]]
    if (n_elements(l5) ne 0) then $
	l5out[copyptr1:copyptr1+blocksize-1]=l5[copylist[2*i]:copylist[2*i+1]]
    if (n_elements(l6) ne 0) then $
	l6out[copyptr1:copyptr1+blocksize-1]=l6[copylist[2*i]:copylist[2*i+1]]
    if (n_elements(l7) ne 0) then $
	l7out[copyptr1:copyptr1+blocksize-1]=l7[copylist[2*i]:copylist[2*i+1]]
    if (n_elements(l8) ne 0) then $
	l8out[copyptr1:copyptr1+blocksize-1]=l8[copylist[2*i]:copylist[2*i+1]]
    copyptr1=copyptr1+blocksize
endfor

if (n_elements(l1) ne 0) then l1=l1out
if (n_elements(l2) ne 0) then l2=l2out
if (n_elements(l3) ne 0) then l3=l3out
if (n_elements(l4) ne 0) then l4=l4out
if (n_elements(l5) ne 0) then l5=l5out
if (n_elements(l6) ne 0) then l6=l6out
if (n_elements(l7) ne 0) then l7=l7out
if (n_elements(l8) ne 0) then l8=l8out

return, out

end

