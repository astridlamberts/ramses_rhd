implicit none
real :: r ! autres possibilites : tableau ( real :: r(3,3) par exemple), ou double precision

!initialisation de la graine
CALL init_random_seed()        
            
!appel de la procÃ©dure qui gÃ©nÃ¨re les nombres aleatoires
CALL RANDOM_NUMBER(r)
	print *,r

end 


! exemple d'initialisation du generateur en utilisant l'heure
SUBROUTINE init_random_seed()
implicit none
INTEGER :: i, n, heure
INTEGER, DIMENSION(:), ALLOCATABLE :: graine

! demande au systeme quelle taille de graine il faut          
CALL RANDOM_SEED(size = n)
! on cree le tableau graine
ALLOCATE(graine(n))

! on regarde l'heure
CALL SYSTEM_CLOCK(COUNT=heure)
          
! on utilise l'heure pour initialiser la graine
graine = heure + 37 * (/ (i - 1, i = 1, n) /)
! (/ (i - 1, i = 1, n) /) veut dire : tableau, de cases 1 a n
! dans la case i, valeur i-1
! apres, on multiplie toutes les valeurs par 37, et on ajoute l'heure

! on met la graine dans le generateur
CALL RANDOM_SEED(PUT = graine)
          
! pour que la graine qu'on a utilisee ne prenne pas de place
DEALLOCATE(graine)

END SUBROUTINE
     

