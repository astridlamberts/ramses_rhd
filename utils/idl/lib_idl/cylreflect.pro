function cylreflect, arr

sz=size(arr)
xsz=sz[1]
zsz=sz[2]
reflarr=dblarr(2*xsz,2*zsz)

reflarr[xsz:2*xsz-1,zsz:2*zsz-1]=arr
reflarr[0:xsz-1,zsz:2*zsz-1]=reverse(arr,1)
reflarr[xsz:2*xsz-1,0:zsz-1]=reverse(arr,2)
reflarr[0:xsz-1,0:zsz-1]=reverse(reflarr[xsz:2*xsz-1,0:zsz-1],1)

return, reflarr

end
