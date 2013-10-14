function floorarr, arr, cut=cut, fac=fac

if not keyword_set(cut) then cut=0.0
if not keyword_set(fac) then fac=0.1

arr1=arr
arr1[where(arr le cut)] = fac*min(arr[where(arr gt cut)])

return, arr1

end
