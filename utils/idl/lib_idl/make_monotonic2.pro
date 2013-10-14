function make_monotonic2, lst

len=n_elements(lst)
minlst=lonarr(len+1)
maxlst=lonarr(len+1)
minptr=0L
maxptr=0L
blockstart=lonarr(len)
blockend=lonarr(len)
blockptr=0L

; first find all the local maxima and minima
minlst[minptr]=0
minptr=minptr+1
for i=1L, len-1 do begin
    if (lst[i] le lst[i-1]) then begin
	minlst[minptr]=i
	minptr=minptr+1
	maxlst[maxptr]=i-1
	maxptr=maxptr+1
    endif
endfor
maxlst[maxptr]=len-1
maxptr=maxptr+1

; work backwards through the minima and maxima
for i=minptr-1,0,-1 do begin

    blockend[blockptr]=
    blockstar[blockptr]=min
