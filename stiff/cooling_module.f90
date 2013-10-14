!====================================================================================================
!                             SUPERMEGAGIGAMODULECOOLINGQUIDEPOTE
!====================================================================================================
! Les subroutines et fonctions d'interet general sont :
!
! ROUTINES A APPELER PAR LE CODE HYDRO
!
!    subroutine set_model(...) : pour choisir le modele de cooling et ses parametres
!
!    subroutine set_table(aexp) : pour creer la table avec les parametres par defaut 
!          Plus pratique a appeler que cmp_table(nH_min,nH_max,T2_min,T2_max,nbin_n,nbin_T,aexp)
!
!    subroutine solve_cooling(...) : pour calculer le cooling
!
! ROUTINE A MODIFIER SI NECESSAIRE
!
!    function J0simple(aexp) : donne le J0 en fonction du redshift dans les modeles Teyssier ou Theuns
!
! AUTRES ROUTINES POUR LES PROGRAMMES D'ANALYSE
!
!    interpolate_table(nH,Tsurmu,n_spec,T,mu) : un utilitaire. Sachant que la table est connue,
!          les abondances n_spec, la temperature T et le poids moleculaire mu sont interpoles pour
!          des valeurs donnees de nH et Tsurmu. 
!
!====================================================================================================
module cooling_module
  use amr_parameters
  implicit none
  logical :: verbose_cooling=.false.

  real(kind=8),parameter ::smallnum_cooling= 1d-30
  real(kind=8),parameter ::twopi   = 6.2831853d0
  real(kind=8),parameter ::hplanck = 6.6262000d-27
  real(kind=8),parameter ::eV      = 1.6022000d-12
  real(kind=8),parameter ::kB      = 1.3806200d-16
  real(kind=8),parameter ::clight  = 2.9979250d+10
  real(kind=8),parameter ::Gyr     = 3.1536000d+16
  real(kind=8),parameter ::X       = 0.76
  real(kind=8),parameter ::Y       = 0.24 
  real(kind=8),parameter ::rhoc    = 1.8800000d-29
  real(kind=8),parameter ::mH      = 1.6600000d-24
  real(kind=8),parameter ::mu_mol  = 1.2195D0
  integer,parameter::HI      = 1
  integer,parameter::HEI     = 2
  integer,parameter::HEII    = 3

  ! Les parametres de la table par defaut
  integer,parameter     :: nbin_T_fix=101
  integer,parameter     :: nbin_n_fix=161
  real(kind=8),parameter:: nH_min_fix=1.d-10
  real(kind=8),parameter:: nH_max_fix=1.d+6
  real(kind=8),parameter:: T2_min_fix=1.d-5
  real(kind=8),parameter:: T2_max_fix=1.d+5
  
  type cooling_table
     integer::n1
     integer::n2
     real(kind=8),dimension(:)    ,pointer::nH
     real(kind=8),dimension(:)    ,pointer::T2
     real(kind=8),dimension(:)    ,pointer::T2eq
     real(kind=8),dimension(:,:)  ,pointer::cool
     real(kind=8),dimension(:,:)  ,pointer::heat
     real(kind=8),dimension(:,:)  ,pointer::mu
     real(kind=8),dimension(:,:,:),pointer::n_spec
  end type cooling_table

  type(cooling_table)::table,table2
  logical, parameter :: if_species_abundances=.true. ! Utilisation de table%n_spec si necessaire

  real(kind=8),parameter :: dumfac_ion_theuns=2.d0    ! Facteur correctif de Theuns et al.
  real(kind=8),parameter :: dumfac_rec_theuns=0.75D0  ! idem
  real(kind=8) :: dumfac_ion=dumfac_ion_theuns
  real(kind=8) :: dumfac_rec=dumfac_rec_theuns


  ! On DOIT AVOIR OU teyssier OU theuns OU madau OU weinberg OU weinbergint avec un OU exclusif
  logical :: teyssier=.false.         
  logical :: theuns=.false.           
  logical :: madau=.false.           
  logical :: weinberg=.false. 
  logical :: weinbergint=.true. 

  ! Si teyssier ou theuns :
  real(kind=8) :: J0in=1.d-22  ! J0 default 
  real(kind=8) :: J0min=1.d-29 ! Valeur minimale du J0 (saturation a grand redshift)
  real(kind=8) :: aexp_ref=0.0001        
  real(kind=8) :: J0min_ref=2.77168510365299962D-25 ! J0min_ref precalcule pour
                                                    ! H0=70, omegab=0.04, omega0=0.3, omegaL=0.7
  logical :: high_z_realistic_ne=.true. ! Calcul du J0min de telle sorte que le n_e soit
                                        ! realiste a grand z. J0min=J0min_ref/(aexp/aexp_ref)^2
  real(kind=8) :: alpha=1.d0   ! J(nu) \propto nu^{-alpha} 
  ! Si madau ou weinberg ou weinbergint :
  real(kind=8) :: normfacJ0=0.74627   ! Facteur de normalisation pour J0 pour un J(nu,z) de type haardt et Madau
                                      ! Ce facteur la est celui utilise par Dave et al. pour LCDM
 

  logical, parameter :: if_cooling_functions=.true. ! Sauvegarde des termes de cooling/heating dans les
                                                    ! variables en dessous
  real(kind=8) ::cb1s,cb2s,cb3s,ci1s,ci2s,ci3s,cr1s,cr2s,cr3s,cds,ce1s,ce3s,ch1s,ch2s,ch3s,cocs,cohs
  real(kind=8) ::cool_out, heat_out

  ! Les heating et photoionization rates de Dave et al. pour le J0 derniere version de HM (weinberg ou weinbergint si
  !   if_read_weinberg=.true. (voir plus bas) dans ce dernier cas)
  real(kind=8),allocatable, dimension(:,:)::table_weinberg   ! Table d'interpolation en input
  character(len=128), parameter :: table_weinberg_name='TREECOOL' ! Nom du fichier avec les donnees
  integer,parameter :: luweinberg=21                         ! unit pour lire le fichier
  integer :: Nweinberg                                       ! Nombre de bins en redshift

  ! Les coefficients d'interpolation des heating rates de Dave et al. (weinbergint)
  logical,parameter :: if_read_weinberg=.false. ! .true. pour lire le fichier table_weinberg_name
                                                ! puis interpoler par un polynome
                                                ! .false. pour utiliser les valeurs des coefficients
                                                ! precalcules listes plus bas
  integer,parameter :: Norderweinberg=7         ! Ordre+1 du polynome d'interpolation (NE PAS CHANGER)
  real(kind=8) :: coefweinberg(Norderweinberg,6)= reshape( &
 &                    (/ -0.31086729929951613D+002, 0.34803667059463761D+001,-0.15145716066316397D+001, &
 &                        0.54649951450632972D+000,-0.16395924120387340D+000, 0.25197466148524143D-001, &
 &                       -0.15352763785487806D-002, &
 &                       -0.31887274113252204D+002, 0.44178493140927095D+001,-0.20158132553082293D+001, &  
 &                        0.64080497292269134D+000,-0.15981267091909040D+000, 0.22056900050237707D-001, &
 &                       -0.12837570029562849D-002, &
 &                       -0.35693331167978656D+002, 0.20207245722165794D+001,-0.76856976101363744D-001, &
 &                       -0.75691470654320359D-001,-0.54502220282734729D-001, 0.20633345104660583D-001, & 
 &                       -0.18410307456285177D-002, &
 &                       -0.56967559787460921D+002, 0.38601174525546353D+001,-0.18318926655684415D+001, &
 &                        0.67360594266440688D+000,-0.18983466813215341D+000, 0.27768907786915147D-001, &
 &                       -0.16330066969315893D-002, &
 &                       -0.56977907250821026D+002, 0.38686249565302266D+001,-0.13330942368518774D+001, &  
 &                        0.33988839029092172D+000,-0.98997915675929332D-001, 0.16781612113050747D-001, &
 &                       -0.11514328893746039D-002, &
 &                       -0.59825233828609278D+002, 0.21898162706563347D+001,-0.42982055888598525D+000, &
 &                        0.50312144291614215D-001,-0.61550639239553132D-001, 0.18017109270959387D-001, & 
 &                       -0.15438891584271634D-002 /), (/Norderweinberg,6/) )

contains 

!=======================================================================
subroutine set_model(Nmodel,J0in_in,J0min_in,alpha_in,normfacJ0_in, &
 &                   correct_cooling,realistic_ne, &
 &                   h,omegab,omega0,omegaL,astart_sim,T2_sim)
!=======================================================================
! Nmodel(integer) =1 : Teyssier : ancien choix de l'evolution et de la forme du J(nu,z)
!                  2 : Theuns   : pareil mais avec les fonctions interpolees de Theuns (+ rapide)
!                  3 : Madau    : J(nu,z) de Theuns et al. 1998 avec les anciennes mesures de 
!                                 Haardt et Madau (HM)
!                  4 : Weinberg : J(nu,z) de Dave et al. 1999 avec les nouvelles mesure de HM 
!                                 lues dans le fichier table_weinberg_name
!                  5 : idem 4 mais interpole interpole de maniere polynomiale : RECOMMANDE
!                 -1 : defaut defini dans le module 
! J0in_in (dble) : valeur du J0 utilisee pour Teyssier et Theuns
!            Exemple : J0in_in=1.d-22
!            J0in_in <= 0 utilise le defaut defini dans le module
! J0min_in (dble) : valeur du J0min ou J0min_ref (voir option realistic_ne) 
!            utilisee dans tous les modeles a grand redshift 
!            Exemple : J0min_in=1.d-29
!            J0min_in <= 0 utilise le defaut defini dans le module
! alpha_in (dble) : valeur de l'indice spectral du J(nu) \propto nu^{-alpha}
!            Exemple : alpha=1.
!            alpha_in < 0 utilise le defaut defini dans le module
! normfacJ0_in (dble) : valeur du facteur de normalisation dans le cas des
!            spectres de Haardt et Madau. C'est un nombre de l'ordre de
!            l'unite en general plus petit que 1.
!            Exemple : normfacJ0_in=0.74627
!            normfacJ0_in prend le defaut defini dans le module
! correct_cooling (integer) : 0 : pas de correction
!                             1 : correction de Theuns et al 98
!                            -1 : defaut defini dans le module
! realistic_ne (integer) : 0 : pas de n_e realiste a grand redshift :
!                              Le J0min reste le meme quel que soit le redshift
!                              (J0min=J0min_in si celui-ci est > 0)
!                          1 : n_e realiste a grand redshift : J0min proportionnel a 1/a^2 
!                              egal initialement a J0min_ref pour a=aexp_ref=0.0001
!                              (J0min_ref=J0min_in si celui-ci est > 0)
!                          2 : RECOMMANDE : pareil que 1, mais J0min_ref est calcule de 
!                              maniere iterative pour avoir le bon n_e a z=19. 
!                              Le J0min_in n'est pas relevant dans ce cas la. 
! h (dble)          : H0/100
! omegab (dble)     : omega baryons
! omega0 (dble)     : omega matiere total
! omegaL (dble)     : omega Lambda
! astart_sim (dble) : redshift auquel on veut commencer la simulation
! T2_sim     (dble) : ce sera en output, le T/mu en K a ce redshift pour des regions de contraste
!                     de densite nul. 
!
! NOTE :
! Dans les cas madau, ou weinberg ou weinbergint, le J0 a grand redshift est calcule comme 
! dans l'option theuns :
!   madau :       pour z >= 15 ou quand le taux trouve est plus petit que celui donne par 
!                 l'option theuns=.true.
!   weinberg :    quand on sort de la table des taux
!   weinbergint : pour z >= 8.5 ou quand le taux trouve est plus petit que celui donne 
!                 par l'option theuns=.true.
!=======================================================================
  implicit none
  real(kind=8) :: J0in_in,J0min_in,alpha_in,normfacJ0_in,astart_sim,T2_sim
  real(kind=8) :: J0min_ref_calc,h,omegab,omega0,omegaL
  integer :: Nmodel,correct_cooling,realistic_ne

  real(kind=8) :: astart,aend,dasura,T2end,mu,ne

  if (Nmodel /= -1) then
     teyssier=.false.
     theuns=.false.
     madau=.false.
     weinberg=.false.
     weinbergint=.false.
     if (Nmodel==1) then
        teyssier=.true.
     elseif (Nmodel==2) then
        theuns=.true.
     elseif (Nmodel==3) then
        madau=.true.
     elseif (Nmodel==4) then
        weinberg=.true.
     elseif (Nmodel==5) then
        weinbergint=.true.
     else
        write(*,*) 'ERROR in set_model : wrong value of Nmodel'
        write(*,*) 'Nmodel =',Nmodel
        STOP
     endif
  endif
  if (J0in_in > 0.d0) J0in=J0in_in
  if (alpha_in > 0.d0) alpha=alpha_in
  if (normfacJ0_in > 0.d0) normfacJ0=normfacJ0_in
  if (correct_cooling == 0) then
     dumfac_ion=1.d0
     dumfac_rec=1.d0
  elseif (correct_cooling == 1) then
     dumfac_ion=dumfac_ion_theuns
     dumfac_rec=dumfac_rec_theuns
  elseif (correct_cooling /= -1) then
     write(*,*) 'ERROR in set_model : wrong value of correct_cooling'
     write(*,*) 'correct_cooling =',correct_cooling
     STOP
  endif
  if (realistic_ne == 0) then
     astart=5.d-4
     high_z_realistic_ne=.false.
     if (J0min_in > 0.d0) J0min=J0min_in
  elseif (realistic_ne == 1) then
     astart=aexp_ref
     high_z_realistic_ne=.true.
     if (J0min_in > 0.d0) J0min_ref=J0min_in
  elseif (realistic_ne == 2) then
     astart=aexp_ref
     high_z_realistic_ne=.true.
     call compute_J0min(h,omegab,omega0,omegaL,J0min_ref_calc)
     J0min_ref=J0min_ref_calc
  else
     write(*,*) 'ERROR in set_model : wrong value of realistic_ne'
     write(*,*) 'realistic_ne =',realistic_ne
     STOP
  endif

  if (astart_sim < astart) then
     write(*,*) 'ERROR in set_model : astart_sim is too small.'
     write(*,*) 'astart     =',astart
     write(*,*) 'astart_sim =',astart_sim
     STOP
  endif

  ! Calcul de la temperature initiale
  aend=astart_sim
  dasura=0.02d0
  call evol_single_cell(astart,aend,dasura,h,omegab,omega0,omegaL, &
 &                          -1.0d0,T2end,mu,ne,.false.)
  if (verbose_cooling) write(*,*) 'Starting temperature in K :',T2end*mu
  T2_sim=T2end 
end subroutine set_model

!=======================================================================
subroutine set_table(aexp)
!=======================================================================
  implicit none

  real(kind=8) :: aexp
  
  integer :: nbin_n,nbin_T
  real(kind=8) :: nH_min,nH_max,T2_min,T2_max

  nH_min=nH_min_fix
  nH_max=nH_max_fix
  T2_min=T2_min_fix
  T2_max=T2_max_fix
  nbin_n=nbin_n_fix
  nbin_T=nbin_T_fix

  call cmp_table(nH_min,nH_max,T2_min,T2_max,nbin_n,nbin_T,aexp)
end subroutine set_table

!=======================================================================
subroutine output_cool(filename)
!=======================================================================
  implicit none
  character(LEN=80)::filename

  open(unit=10,file=filename,form='unformatted')
  write(10)table%n1,table%n2
  write(10)table%nH
  write(10)table%T2
  write(10)table%T2eq
  write(10)table%cool
  write(10)table%heat
  write(10)table%mu
  if (if_species_abundances) write(10)table%n_spec
  close(10)

end subroutine output_cool

!=======================================================================
subroutine interpolate_table(nH,Tsurmu,n_spec,T,mu)
!=======================================================================
! Une subroutine utilitaire pour calculer les abondances (n_spec), la
! temperature (T) et le poids moleculaire (mu) pour une valeur donnee
! de nH et de Tsurmu, sachant que la table et nspec sont connus.
! Pour utiliser cette subroutine, il faut avoir if_species_abundances=.true.
!=======================================================================
  implicit none
  real(kind=8) :: nH, Tsurmu, mu
  real(kind=8) :: n_spec(1:6),T
  real(kind=8) :: facT,facH,T2eq
  integer :: i_T2,i_nH,j
  real(kind=8) :: w1T,w2T,w1H,w2H,w11,w12,w21,w22,dlog_nH,dlog_T2


  dlog_nH = dble(table%n1-1)/(table%nH(table%n1)-table%nH(1))
  dlog_T2 = dble(table%n2-1)/(table%T2(table%n2)-table%T2(1))
  facH  = log10(nH)
  i_nH=  MIN(MAX(int((facH-table%nH(1))*dlog_nH)+1,1),table%n1-1) 
  w1H   = (table%nH(i_nH+1)-facH)*dlog_nH
  w2H   = (facH-table%nH(i_nH  ))*dlog_nH
  T2eq=table%T2eq(i_nH)**w1H*table%T2eq(i_nH+1)**w2H
  facT  = log10(Tsurmu/T2eq)
  i_T2 = MIN(MAX(int((facT-table%T2(1))*dlog_T2)+1,1),table%n2-1)
  w1T   = (table%T2(i_T2+1)-facT)*dlog_T2
  w2T   = (facT-table%T2(i_T2  ))*dlog_T2
  w11=w1T*w1H
  w12=w1T*w2H
  w21=w2T*w1H
  w22=w2T*w2H
  do j=1,6
     n_spec(j)=10.d0**( log10(table%n_spec(i_nH  ,i_T2  ,j))*w11 &
          &            +log10(table%n_spec(i_nH  ,i_T2+1,j))*w21 &
          &            +log10(table%n_spec(i_nH+1,i_T2  ,j))*w12 &
          &            +log10(table%n_spec(i_nH+1,i_T2+1,j))*w22 )
  enddo
  mu=10.d0**( log10(table%mu(i_nH  ,i_T2  ))*w11 &
    &        +log10(table%mu(i_nH  ,i_T2+1))*w21 &
    &        +log10(table%mu(i_nH+1,i_T2  ))*w12 &
    &        +log10(table%mu(i_nH+1,i_T2+1))*w22 )
  T=Tsurmu*mu
end subroutine interpolate_table


!=======================================================================
subroutine evol_single_cell(astart,aend,dasura,h,omegab,omega0,omegaL, &
 &                          J0min_in,T2end,mu,ne,if_write_result)
!=======================================================================
! astart : valeur du facteur d'expansion au debut du calcul
! aend   : valeur du facteur d'expansion a la fin du calcul
! dasura : la valeur de da/a entre 2 pas de temps
! h      : la valeur de H0/100 
! omegab : la valeur de Omega baryons
! omega0 : la valeur de Omega matiere (total)
! omegaL : la valeur de Omega Lambda
! J0min_in : la valeur du J0min a injecter :
!          Si high_z_realistic_ne alors c'est J0min a a=astart qui
!          est considere
!          Sinon, c'est le J0min habituel.
!          Si J0min_in <=0, les parametres par defaut ou predefinis
!          auparavant sont pris pour le J0min.
! T2end  : Le T/mu en output
! mu     : le poids moleculaire en output
! ne     : le ne en output
! if_write_result : .true. pour ecrire l'evolution de la temperature
!          et de n_e sur l'ecran.
!=======================================================================
  implicit none
  real(kind=8)::astart,aend,T2end,h,omegab,omega0,omegaL,J0min_in,ne,dasura
  logical :: if_write_result

  real(kind=8)::aexp,daexp,dt_cool,coeff
  real(kind=8)::T2_com,T2_old,T2,T2_left,T2_right,err_T2
  real(kind=8)::nH_com,nH
  
  real(kind=8),dimension(1:3)::t_rad_spec,h_rad_spec
  real(kind=8) ::mu
  real(kind=8) ::cool_tot,heat_tot
  real(kind=8) ::diff
  integer::niter
  real(kind=8) :: n_spec(1:6)


  if (J0min_in > 0.d0) then
     if (high_z_realistic_ne) then
        J0min_ref = J0min_in
        aexp_ref = astart
     else
        J0min = J0min_in
     endif
  endif
  aexp = astart
  T2_com = 2.726d0 / aexp * aexp**2 / mu_mol
  nH_com = omegab*rhoc*h**2*X/mH
  do while (aexp < aend)
     daexp = dasura*aexp
     dt_cool=daexp/(aexp*100.*h*3.2408608e-20*HsurH0(1.0/aexp-1.,omega0,omegaL,1.-omega0-omegaL))
     
     nH = nH_com/aexp**3
     T2_old = T2_com/aexp**2

     ! Compute radiative ionization and heating rates
     call set_rates(t_rad_spec,h_rad_spec,aexp)
     
     ! Iteration to find new T2
     err_T2=1.
     T2_left=1.d-2
     T2_right=1.d8
     niter=0
     coeff = 2.*nH*X/3./kB
     do while (err_T2 > 1.d-10.and.niter <= 100)
        T2=0.5*(T2_left+T2_right)        
        call cmp_cooling(T2,nH,t_rad_spec,h_rad_spec,cool_tot,heat_tot,mu,aexp,n_spec)
        diff = coeff*(heat_tot-cool_tot) +(T2_old-T2)/dt_cool
        if(diff>0.)then 
           T2_left =0.5*(T2_left+T2_right)
           T2_right=T2_right
        else
           T2_left =T2_left
           T2_right=0.5*(T2_left+T2_right)
        end if
        err_T2=abs(T2_right-T2_left)/T2_left
        niter=niter+1
     end do
     if (niter > 100) then
        write(*,*) 'ERROR in evol_single_cell : too many iterations'
        STOP
     endif

     T2_com=T2*aexp**2

     aexp = aexp + daexp

     if (if_write_result) write(*,'(4(1pe10.3))')aexp,nH,T2_com*mu/aexp**2,n_spec(1)/nH

  end do

  T2end=T2
  ne=n_spec(1)/nH
end subroutine evol_single_cell

!=======================================================================
subroutine compute_J0min(h,omegab,omega0,omegaL,J0min_in)
!=======================================================================
  implicit none
  real(kind=8) :: omega0,omegaL,h,omegab,ne_to_find,mu
  real(kind=8) :: h0,astart,aend,J0min_in,T2end,ne
  real(kind=8) :: J0min_left,J0min_right,err_J0min,diff,xval,dasura
  integer :: niter
  logical :: if_write_result=.false.

  xval=sqrt(omega0)/(h*omegab)
  ne_to_find=1.2d-5*xval ! From the book of Peebles p. 173
  astart=aexp_ref
  aend=0.05 ! Attention a cette valeur. Il ne faut pas que la reionisation
            ! commence avant (pour le moment c'est le cas, car dans tous
            ! les modeles consideres, la reionisation commence a z < 15)
  dasura=0.05

  err_J0min=1.
  J0min_left=1d-20
  J0min_right=1d-30
  niter=0

  do while (err_J0min > 1.d-3 .and. niter <= 100)
     J0min_in=0.5*(J0min_left+J0min_right)     
     call evol_single_cell(astart,aend,dasura,h,omegab,omega0,omegaL, &
 &                          J0min_in,T2end,mu,ne,if_write_result)
     diff=ne-ne_to_find
     if (diff>0.d0) then
        J0min_left=0.5*(J0min_left+J0min_right)
        J0min_right=J0min_right
     else
        J0min_left=J0min_left
        J0min_right=0.5*(J0min_left+J0min_right)
     endif
     err_J0min=abs(J0min_right-J0min_left)/J0min_left
     niter=niter+1
  enddo
  if (niter > 100) then
     write(*,*) 'ERROR in compute_J0min : too many iterations'
     STOP
  endif

  if (verbose_cooling)  write(*,*) 'J0min found ',J0min_in
end subroutine compute_J0min

!=======================================================================
subroutine solve_cooling(nH,T2,dt,dT2dt,ncell)
!=======================================================================
  implicit none  
  integer::ncell
  real(kind=8)::dt
  real(kind=8),dimension(1:ncell)::nH,T2,dT2dt
    
  real(kind=8)::facT,facH,dlog_nH,dlog_T2,coeff,T2eq
  real(kind=8)::cool,heat,w1T,w2T,w1H,w2H,w11,w12,w21,w22
  real(kind=8),dimension(1:ncell)::err,rgt,lft,tau,tau_old,rgt_new,lft_new
  real(kind=8),dimension(1:ncell)::time,time_new,rho,rho_new,tau_old_new

  integer::i,i_T2,i_nH,iter,n,n_active
  integer,dimension(1:ncell)::ind,ind_new,iii

  logical::tau_negative
  
  
  dlog_nH = dble(table%n1-1)/(table%nH(table%n1)-table%nH(1))
  dlog_T2 = dble(table%n2-1)/(table%T2(table%n2)-table%T2(1))
  
  do i=1,ncell
     rho(i)=log10(nH(i))
  end do
  
  ! Implicit method with dichotomy
  do i=1,ncell
     facH=rho(i)
     i_nH=  MIN(MAX(int((facH-table%nH(1))*dlog_nH)+1,1),table%n1-1) 
     w1H   = (table%nH(i_nH+1)-facH)*dlog_nH
     w2H   = (facH-table%nH(i_nH  ))*dlog_nH
     T2eq=table%T2eq(i_nH)**w1H*table%T2eq(i_nH+1)**w2H
     tau(i) = T2(i)/T2eq
     coeff = 2.*nH(i)*X/3./kB/T2eq
     time(i) = dt*coeff
  end do
  
  do i=1,ncell
     if(tau(i) >=1.)then
        lft(i)=1.
        rgt(i)=tau(i)
     else
        lft(i)=tau(i)
        rgt(i)=1.
     end if
  end do
  
  tau_old=tau
  err=1.
  iter=0
  n=ncell
  do i = 1,n
     ind(i)=i
  end do
  
  ! Loop over active cells
  do while (n > 0)
     
     iter=iter+1
     if (iter > 50) then
        write(*,*) 'Too many iterations in solve_cooling'
        STOP
     endif
     
     ! Check positivity 
     tau_negative=.false.
     do i=1,n
        tau(i) = 0.5*(lft(i)+rgt(i))
        if (tau(i) <= 0.) tau_negative=.true.
     end do
     if (tau_negative) then
        write(*,*) 'ERROR in solve_cooling :'
        write(*,*) 'Temperature is negative'
        STOP
     endif
     
     do i=1,n
        facT  = log10(tau(i))
        facH  = rho(i)
        i_T2 = MIN(MAX(int((facT-table%T2(1))*dlog_T2)+1,1),table%n2-1)
        i_nH=  MIN(MAX(int((facH-table%nH(1))*dlog_nH)+1,1),table%n1-1) 
        w1T   = (table%T2(i_T2+1)-facT)*dlog_T2
        w2T   = (facT-table%T2(i_T2  ))*dlog_T2
        w1H   = (table%nH(i_nH+1)-facH)*dlog_nH
        w2H   = (facH-table%nH(i_nH  ))*dlog_nH
        w11=w1T*w1H
        w12=w1T*w2H
        w21=w2T*w1H
        w22=w2T*w2H
        heat = 10.d0**(table%heat(i_nH  ,i_T2  )*w11 &
             &        +table%heat(i_nH  ,i_T2+1)*w21 &
             &        +table%heat(i_nH+1,i_T2  )*w12 &
             &        +table%heat(i_nH+1,i_T2+1)*w22)
        cool = 10.d0**(table%cool(i_nH  ,i_T2  )*w11 &
             &        +table%cool(i_nH  ,i_T2+1)*w21 &
             &        +table%cool(i_nH+1,i_T2  )*w12 &
             &        +table%cool(i_nH+1,i_T2+1)*w22)
        err(i) = tau(i) - tau_old(i) - (heat-cool)*time(i)
     end do
     
     do i=1,n
        if(err(i) >= 0.)then
           rgt(i)=tau(i)
        else
           lft(i)=tau(i)
        end if
     end do
     
     do i=1,n
        err(i)=ABS(rgt(i)-lft(i))/tau(i)
     end do
     
     n_active=count(mask=err(1:n) > 1.d-5)
     if (n_active < n .and. n_active > 0) then
        iii(1:n_active)=pack(array=(/ (i,i=1,n) /),mask=err(1:n) > 1.d-5)
        iii(n_active+1:n)=pack(array=(/ (i,i=1,n) /),mask=err(1:n) <= 1.d-5)
        do i=1,n
           ind_new(i)=ind(iii(i))
           lft_new(i)=lft(iii(i))
           rgt_new(i)=rgt(iii(i))
           tau_old_new(i)=tau_old(iii(i))
           time_new(i)=time(iii(i))
           rho_new(i)=rho(iii(i))
        enddo
        ind(1:n)=ind_new(1:n)
        lft(1:n)=lft_new(1:n)
        rgt(1:n)=rgt_new(1:n)
        tau_old(1:n)=tau_old_new(1:n)
        time(1:n)=time_new(1:n)
        rho(1:n)=rho_new(1:n)
     endif
     n=n_active
  end do
  ! End loop over active cells

  ! Check positivity 
  tau_negative=.false.
  do i=1,ncell
     tau(i) = 0.5*(lft(i)+rgt(i))
     if (tau(i) <= 0.) tau_negative=.true.
  end do
  
  if (tau_negative) then
     write(*,*) 'ERROR in solve_cooling :'
     write(*,*) 'Temperature is negative'
     STOP
  endif
  
    ! Compute dtau/dt
  do i=1,ncell
     tau(i)=tau(i)-tau_old(i)
  end do

  ! Store results in array dT2dt
!CDIR NODEP
  do i=1,ncell
     facH=rho(i)
     i_nH=  MIN(MAX(int((facH-table%nH(1))*dlog_nH)+1,1),table%n1-1) 
     w1H   = (table%nH(i_nH+1)-facH)*dlog_nH
     w2H   = (facH-table%nH(i_nH  ))*dlog_nH
     dT2dt(ind(i))=tau(i)*(table%T2eq(i_nH)**w1H*table%T2eq(i_nH+1)**w2H)
  end do
  
end subroutine solve_cooling

!=======================================================================
function J0simple(aexp)
!=======================================================================
! Le J0 dans le cas teyssier ou theuns
!=======================================================================
  real(kind=8) :: J0simple,aexp
  if (aexp .lt. 1.d0/7.d0) then
     J0simple=0.d0
  elseif (aexp .lt. 1.d0/4.d0)then
     J0simple=4.d0*aexp
  elseif (aexp .lt. 1.d0/3.d0)then
     J0simple=1.d0
  else
     J0simple=1.d0/(3.*aexp)**3
  endif
  J0simple=max(J0simple*J0in,J0min)
  return
end function J0simple

!=======================================================================
subroutine cmp_table(nH_min,nH_max,T2_min,T2_max,nbin_n,nbin_T,aexp)
!=======================================================================
  implicit none
  include 'mpif.h'

  real(kind=8)::nH_min,nH_max,T2_min,T2_max,aexp
  integer::nbin_n,nbin_T
  integer::myid,ncpu,ierr
  integer::i_n,i_T
  real(kind=8),dimension(1:3)::t_rad_spec,h_rad_spec
  integer :: i,j,n1,n2
  logical :: first=.true.
  save first
  
  call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,ncpu,ierr)

  if(.not.first)then
     deallocate(table%cool)
     deallocate(table%heat)
     deallocate(table%mu)
     deallocate(table%T2)
     deallocate(table%nH)
     deallocate(table%T2eq)
     if (if_species_abundances) deallocate(table%n_spec)
  else
     first=.false.
  endif
  
  table%n1=nbin_n
  table%n2=nbin_T
  allocate(table%cool(nbin_n,nbin_T))
  allocate(table%heat(nbin_n,nbin_T))
  allocate(table%mu  (nbin_n,nbin_T))
  allocate(table%T2eq(nbin_n))
  allocate(table2%cool(nbin_n,nbin_T))
  allocate(table2%heat(nbin_n,nbin_T))
  allocate(table2%mu  (nbin_n,nbin_T))
  allocate(table2%T2eq(nbin_n))
  allocate(table%nH  (nbin_n))
  allocate(table%T2  (nbin_T))
  if (if_species_abundances) allocate(table%n_spec(nbin_n,nbin_T,1:6))
  if (if_species_abundances) allocate(table2%n_spec(nbin_n,nbin_T,1:6))
  
  do i_n=1,nbin_n
     table%nH(i_n)=nH_min*10.**(dble(i_n-1)/dble(nbin_n-1) &
            & *log10(nH_max/nH_min))
  end do
  
  do i_T=1,nbin_T
     table%T2(i_T)=T2_min*10.**(dble(i_T-1)/dble(nbin_T-1) &
          & *log10(T2_max/T2_min))
  end do
  
  ! Compute radiative ionization and heating rates
  call set_rates(t_rad_spec,h_rad_spec,aexp)

  ! Create the table
  table%T2eq=0.0
  table%mu=0.0
  table%cool=0.0
  table%heat=0.0
  if (if_species_abundances) table%n_spec=0.0
  do i_n = myid+1,nbin_n,ncpu
     call iterate(i_n,t_rad_spec,h_rad_spec,nbin_T,aexp)
  end do
  call MPI_ALLREDUCE(table%T2eq,table2%T2eq,nbin_n,&
       &MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(table%mu,table2%mu,nbin_n*nbin_T,&
       &MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(table%cool,table2%cool,nbin_n*nbin_T,&
       &MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call MPI_ALLREDUCE(table%heat,table2%heat,nbin_n*nbin_T,&
       &MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  if (if_species_abundances)then
     call MPI_ALLREDUCE(table%n_spec,table2%n_spec,6*nbin_n*nbin_T,&
          &MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  end if

  ! Convert to log scale
  table%nH   = log10(table%nH)
  table%T2   = log10(table%T2)
  table%cool = log10(table2%cool)
  table%heat = log10(table2%heat)
  if (if_species_abundances) table%n_spec=table2%n_spec
  table%T2eq = table2%T2eq
  table%mu   = table2%mu

  deallocate(table2%cool)
  deallocate(table2%heat)
  deallocate(table2%mu)
  deallocate(table2%T2eq)
  if (if_species_abundances) deallocate(table2%n_spec)

end subroutine cmp_table

!=======================================================================
subroutine set_rates(t_rad_spec,h_rad_spec,aexp)
!=======================================================================
  implicit none

  real(kind=8),dimension(1:3) :: t_rad_spec,h_rad_spec
  real(kind=8) :: J0,z,aexp
  logical :: first=.true.
  save first

  z=1.d0/aexp-1.D0
  if (high_z_realistic_ne) J0min=J0min_ref/(aexp/aexp_ref)**2

  if (teyssier) then
     J0=J0simple(aexp)
     t_rad_spec(HI  ) = taux_rad(HI  ,J0)
     t_rad_spec(HEI ) = taux_rad(HEI ,J0)
     t_rad_spec(HEII) = taux_rad(HEII,J0)
     h_rad_spec(HI  ) = heat_rad(HI  ,J0)
     h_rad_spec(HEI ) = heat_rad(HEI ,J0)
     h_rad_spec(HEII) = heat_rad(HEII,J0)
  elseif (theuns) then
     J0=J0simple(aexp)
     t_rad_spec(HI  ) = taux_rad_theuns(HI  ,J0)
     t_rad_spec(HEI ) = taux_rad_theuns(HEI ,J0)
     t_rad_spec(HEII) = taux_rad_theuns(HEII,J0)
     h_rad_spec(HI  ) = heat_rad_theuns(HI  ,J0)
     h_rad_spec(HEI ) = heat_rad_theuns(HEI ,J0)
     h_rad_spec(HEII) = heat_rad_theuns(HEII,J0)
  elseif (madau) then
     z=1.d0/aexp-1.D0
     t_rad_spec(HI  ) = taux_rad_madau(HI  ,z)
     t_rad_spec(HEI ) = taux_rad_madau(HEI ,z)
     t_rad_spec(HEII) = taux_rad_madau(HEII,z)
     h_rad_spec(HI  ) = heat_rad_madau(HI  ,z)
     h_rad_spec(HEI ) = heat_rad_madau(HEI ,z)
     h_rad_spec(HEII) = heat_rad_madau(HEII,z)     
  elseif (weinberg) then
     if (first) then
        call read_weinberg
        first=.false.
     endif
     t_rad_spec(HI  ) = taux_rad_weinberg(HI  ,z)
     t_rad_spec(HEI ) = taux_rad_weinberg(HEI ,z)
     t_rad_spec(HEII) = taux_rad_weinberg(HEII,z)
     h_rad_spec(HI  ) = heat_rad_weinberg(HI  ,z)
     h_rad_spec(HEI ) = heat_rad_weinberg(HEI ,z)
     h_rad_spec(HEII) = heat_rad_weinberg(HEII,z)     
  elseif (weinbergint) then
     if (first.and.if_read_weinberg) then
        call read_weinberg
        first=.false.
     endif
     t_rad_spec(HI  ) = taux_rad_weinbergint(HI  ,z)
     t_rad_spec(HEI ) = taux_rad_weinbergint(HEI ,z)
     t_rad_spec(HEII) = taux_rad_weinbergint(HEII,z)
     h_rad_spec(HI  ) = heat_rad_weinbergint(HI  ,z)
     h_rad_spec(HEI ) = heat_rad_weinbergint(HEI ,z)
     h_rad_spec(HEII) = heat_rad_weinbergint(HEII,z)     
  endif  
end subroutine set_rates


!=======================================================================
subroutine iterate(i_n,t_rad_spec,h_rad_spec,nbin_T,aexp)
!=======================================================================
  implicit none
  integer :: i_n
  real(kind=8),dimension(1:3)::t_rad_spec,h_rad_spec
  real(kind=8) :: aexp
  
  integer::nbin_T    
  integer::i_T
  real(kind=8) ::T2,nH
  real(kind=8) ::mu
  real(kind=8) ::T2_left,T2_right,err_T2
  real(kind=8) ::cool_tot,heat_tot
  real(kind=8) ::diff
  integer::niter
  real(kind=8) :: n_spec(1:6)
  
  nH=table%nH(i_n)
  
  ! Iteration to find T2_eq
  err_T2=1.
  T2_left=1.d-2
  T2_right=1.d8
  niter=0
  do while (err_T2 > 1.d-10.and.niter <= 100)
     T2=0.5*(T2_left+T2_right)        
     call cmp_cooling(T2,nH,t_rad_spec,h_rad_spec,cool_tot,heat_tot,mu,aexp,n_spec)
     diff = heat_tot-cool_tot
     if(diff>0.)then 
        T2_left =0.5*(T2_left+T2_right)
        T2_right=T2_right
     else
        T2_left =T2_left
        T2_right=0.5*(T2_left+T2_right)
     end if
     err_T2=abs(T2_right-T2_left)/T2_left
     niter=niter+1
  end do
  if (niter > 100) then
     write(*,*) 'ERROR in iterate : too many iterations'
     STOP
  endif
  table%T2eq(i_n)=T2
  
       
  do i_T = 1,nbin_T
     nH=table%nH(i_n)
     T2=table%T2(i_T)*table%T2eq(i_n)
     call cmp_cooling(T2,nH,t_rad_spec,h_rad_spec,cool_tot,heat_tot,mu,aexp,n_spec)
     ! Save to output arrays
     table%cool(i_n,i_T)=cool_tot
     table%heat(i_n,i_T)=heat_tot
     table%mu  (i_n,i_T)=mu
     if (if_species_abundances) table%n_spec(i_n,i_T,1:6)=n_spec(1:6)
  end do
end subroutine iterate
  

!=======================================================================
subroutine cmp_cooling(T2,nH,t_rad_spec,h_rad_spec,cool_tot,heat_tot,mu_out,aexp,n_spec)
!=======================================================================
  implicit none
  
  real(kind=8),dimension(1:3)::t_rad_spec,h_rad_spec
  real(kind=8) ::T2,nH,cool_tot,heat_tot,mu_out,aexp
  
  real(kind=8) ::mu,mu_old,err_mu,mu_left,mu_right
  real(kind=8) ::T
  real(kind=8) ::n_E,n_HI,n_HII,n_HEI,n_HEII,n_HEIII,n_TOT
  
  real(kind=8),dimension(1:6)::n_spec
  real(kind=8) ::cb1,cb2,cb3,ci1,ci2,ci3,cr1,cr2,cr3,cd,ce1,ce2,ce3,ch1,ch2,ch3,coc,coh
  
  integer :: niter
  
  
  ! Iteration to find mu
  err_mu=1.
  mu_left=0.5
  mu_right=1.3
  niter=0
  do while (err_mu > 1.d-4 .and. niter <= 50)
     mu_old=0.5*(mu_left+mu_right)
     T = T2*mu_old
     call cmp_chem_eq(T,nH,t_rad_spec,n_spec,n_TOT,mu)
     err_mu = (mu-mu_old)/mu_old
     if(err_mu>0.)then 
        mu_left =0.5*(mu_left+mu_right)
        mu_right=mu_right
     else
        mu_left =mu_left
        mu_right=0.5*(mu_left+mu_right)
     end if
     err_mu=ABS(err_mu)
     niter=niter+1
  end do
  if (niter > 50) then
     write(*,*) 'ERROR in cmp_cooling : too many iterations.'
     STOP
  endif
    
  ! Get equilibrium abundances
  n_E     = n_spec(1) ! electrons
  n_HI    = n_spec(2) ! H
  n_HII   = n_spec(3) ! H+
  n_HEI   = n_spec(4) ! He
  n_HEII  = n_spec(5) ! He+
  n_HEIII = n_spec(6) ! He++
  ! Bremstrahlung
  cb1 = cool_bre(HI  ,T)*n_E*n_HII  /nH**2
  cb2 = cool_bre(HEI ,T)*n_E*n_HEII /nH**2
  cb3 = cool_bre(HEII,T)*n_E*n_HEIII/nH**2
  ! Ionization cooling
  ci1 = cool_ion(HI  ,T)*n_E*n_HI   /nH**2
  ci2 = cool_ion(HEI ,T)*n_E*n_HEI  /nH**2
  ci3 = cool_ion(HEII,T)*n_E*n_HEII /nH**2
  ! Recombination cooling
  cr1 = cool_rec(HI  ,T)*n_E*n_HII  /nH**2
  cr2 = cool_rec(HEI ,T)*n_E*n_HEII /nH**2
  cr3 = cool_rec(HEII,T)*n_E*n_HEIII/nH**2
  ! Dielectric recombination cooling
  cd  = cool_die(T     )*n_E*n_HEII /nH**2
  ! Line cooling
  ce1 = cool_exc(HI  ,T)*n_E*n_HI   /nH**2
!  ce2 = cool_exc(HEI, T)*n_E*n_HEI  /nH**2 Terme elimine car cubique en nH d'apres Cen
  ce3 = cool_exc(HEII,T)*n_E*n_HEII /nH**2
  ! Compton cooling
  coc  = cool_com(T,aexp)*n_E        /nH**2
  ! Radiative heating
  ch1 = h_rad_spec(HI  )    *n_HI   /nH**2
  ch2 = h_rad_spec(HEI )    *n_HEI  /nH**2
  ch3 = h_rad_spec(HEII)    *n_HEII /nH**2
  ! Compton heating
  coh = heat_com(T,aexp)*n_E        /nH**2
  ! Total cooling and heating rates
  heat_tot = ch1+ch2+ch3+coh
!  cool_tot = cb1+cb2+cb3+ci1+ci2+ci3+cr1+cr2+cr3+cd+ce1+ce2+ce3+coc
  cool_tot = cb1+cb2+cb3+ci1+ci2+ci3+cr1+cr2+cr3+cd+ce1+ce3+coc
  mu_out  =mu
  
  if (if_cooling_functions) then
     cool_out=max(cool_tot,smallnum_cooling)
     heat_out=max(heat_tot,smallnum_cooling)
     cb1s=max(cb1,smallnum_cooling)
     cb2s=max(cb2,smallnum_cooling)
     cb3s=max(cb3,smallnum_cooling)
     ci1s=max(ci1,smallnum_cooling)
     ci2s=max(ci2,smallnum_cooling)
     ci3s=max(ci3,smallnum_cooling)
     cr1s=max(cr1,smallnum_cooling)
     cr2s=max(cr2,smallnum_cooling)
     cr3s=max(cr3,smallnum_cooling)
     cds =max(cd ,smallnum_cooling)
     ce1s=max(ce1,smallnum_cooling)
     ce3s=max(ce3,smallnum_cooling)
     cocs=max(coc,smallnum_cooling)
     cohs=max(coh,smallnum_cooling)
     ch1s=max(ch1,smallnum_cooling)
     ch2s=max(ch2,smallnum_cooling)
     ch3s=max(ch3,smallnum_cooling)
     cohs=max(coh,smallnum_cooling)  
  endif
end subroutine cmp_cooling

!=======================================================================
subroutine cmp_chem_eq(T,n_H,t_rad_spec,n_spec,n_TOT,mu)
!=======================================================================
  implicit none
  real(kind=8)::T,n_H,n_TOT,mu
  real(kind=8),dimension(1:3)::t_rad_spec
  real(kind=8),dimension(1:6)::n_spec
  real(kind=8)::xx,yy
  real(kind=8)::n_HI,n_HII,n_HEI,n_HEII,n_HEIII,n_E
  real(kind=8)::t_rad_HI,t_rad_HEI,t_rad_HEII
  real(kind=8)::t_rec_HI,t_rec_HEI,t_rec_HEII
  real(kind=8)::t_ion_HI,t_ion_HEI,t_ion_HEII
  real(kind=8)::t_ion2_HI,t_ion2_HEI,t_ion2_HEII
  real(kind=8)::x1,err_nE
  
  xx=(1.-Y)
  yy=Y/(1.-Y)/4.
  
  t_rad_HI   = t_rad_spec(HI)
  t_rec_HI   = taux_rec  (HI,T)
  t_ion_HI   = taux_ion  (HI,T)
  
  t_rad_HEI  = t_rad_spec(HEI)
  t_rec_HEI  = taux_rec  (HEI,T)
  t_ion_HEI  = taux_ion  (HEI,T)
  
  t_rad_HEII = t_rad_spec(HEII)
  t_rec_HEII = taux_rec  (HEII,T)
  t_ion_HEII = taux_ion  (HEII,T)
  
  n_E = n_H        
  err_nE = 1.
  
  do while(err_nE > 1.d-8)
     
     t_ion2_HI   = t_ion_HI   + t_rad_HI  /MAX(n_E,1e-15*n_H)
     t_ion2_HEI  = t_ion_HEI  + t_rad_HEI /MAX(n_E,1e-15*n_H)
     t_ion2_HEII = t_ion_HEII + t_rad_HEII/MAX(n_E,1e-15*n_H)
     
     n_HI  = t_rec_HI/(t_ion2_HI+t_rec_HI)*n_H
     n_HII = t_ion2_HI/(t_ion2_HI+t_rec_HI)*n_H
     
     x1 = (t_rec_HEII*t_rec_HEI+t_ion2_HEI*t_rec_HEII+t_ion2_HEII*t_ion2_HEI)
     
     n_HEIII = yy*t_ion2_HEII*t_ion2_HEI/x1*n_H
     n_HEII  = yy*t_ion2_HEI *t_rec_HEII/x1*n_H
     n_HEI   = yy*t_rec_HEII *t_rec_HEI /x1*n_H
     
     err_nE = ABS((n_E - (n_HII + n_HEII + 2.*n_HEIII))/n_H)
     n_E = 0.5*n_E+0.5*(n_HII + n_HEII + 2.*n_HEIII)
     
  end do
    
  n_TOT    =n_E+n_HI+n_HII+n_HEI+n_HEII+n_HEIII
  mu       =n_H/xx/n_TOT
  n_spec(1)=n_E
  n_spec(2)=n_HI
  n_spec(3)=n_HII
  n_spec(4)=n_HEI
  n_spec(5)=n_HEII
  n_spec(6)=n_HEIII
  
end subroutine cmp_chem_eq

!=======================================================================
subroutine read_weinberg
!=======================================================================
  implicit none
  real(kind=8) :: line(7)
  integer :: i

  open(unit=luweinberg,file=table_weinberg_name,form='formatted',status='old',err=2)
  Nweinberg=0
  do
     read(luweinberg,*,end=1) line
     if (line(2) > 0) Nweinberg=Nweinberg+1
  enddo

1 continue
  allocate(table_weinberg(7,Nweinberg))
  rewind(luweinberg)
  read(luweinberg,*) table_weinberg
  close(luweinberg)
  if (weinbergint) call interpolate_weinberg
  return

2 write(*,*) 'ERROR in read_weinberg : I can''t open the input table.'
  STOP  
end subroutine read_weinberg

!=======================================================================
subroutine interpolate_weinberg
!=======================================================================
  implicit none
  integer :: i,j,k,ikind,mp,m
  real(kind=8) :: mat(Norderweinberg,Norderweinberg)
  real(kind=8) :: matin(Norderweinberg,Norderweinberg)
  real(kind=8) :: veca(Norderweinberg),vecb(Norderweinberg)
  real(kind=8) :: redshift_j
  do k=1,Norderweinberg
     do i=1,Norderweinberg
        mat(i,k)=0.d0
        do j=1,Nweinberg
           redshift_j=10.d0**table_weinberg(1,j)-1.
           mat(i,k)=mat(i,k)+redshift_j**(i+k-2)
        enddo
     enddo
  enddo
  do ikind=2,7
     do k=1,Norderweinberg    
        vecb(k)=0.d0
        do j=1,Nweinberg
           redshift_j=10.d0**table_weinberg(1,j)-1.
           vecb(k)=vecb(k)+log(table_weinberg(ikind,j))*redshift_j**(k-1)
        enddo
     enddo
     matin=mat     
     mp=1
     m=1
     call gaussj(matin,Norderweinberg,Norderweinberg,vecb,m,mp)
     coefweinberg(1:Norderweinberg,ikind-1)=vecb(1:Norderweinberg)
  enddo
  if (verbose_cooling) then
     write(*,*) 'Parameters of interpolation in interpolate_weinberg'
     do i=1,Norderweinberg
        write(*,'(6(e10.3))') real(coefweinberg(i,1:6))
     enddo
  endif
end subroutine interpolate_weinberg

!=======================================================================
SUBROUTINE gaussj(a,n,np,b,m,mp)
!=======================================================================
! numrec2, p. 30
!=======================================================================
  implicit real(kind=8) (a-h,o-z)
  INTEGER :: m,mp,n,np,NMAX
  real(kind=8) :: a(np,np),b(np,mp)
  PARAMETER (NMAX=500)
  INTEGER i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
  real(kind=8) big,dum,pivinv
  do j=1,n
     ipiv(j)=0
  enddo
  do i=1,n
     big=0.
     do j=1,n
        if(ipiv(j).ne.1)then
           do k=1,n
              if (ipiv(k).eq.0) then
                 if (abs(a(j,k)).ge.big)then
                    big=abs(a(j,k))
                    irow=j
                    icol=k
                 endif
              else if (ipiv(k).gt.1) then
                 STOP 'singular matrix in gaussj'
              endif
           enddo
        endif
     enddo
     ipiv(icol)=ipiv(icol)+1
     if (irow.ne.icol) then
        do l=1,n
           dum=a(irow,l)
           a(irow,l)=a(icol,l)
           a(icol,l)=dum
        enddo
        do l=1,m
           dum=b(irow,l)
           b(irow,l)=b(icol,l)
           b(icol,l)=dum
        enddo
     endif
     indxr(i)=irow
     indxc(i)=icol
     if (a(icol,icol).eq.0.) STOP 'singular matrix in gaussj'
     pivinv=1./a(icol,icol)
     a(icol,icol)=1.
     do l=1,n
        a(icol,l)=a(icol,l)*pivinv
     enddo
     do l=1,m
        b(icol,l)=b(icol,l)*pivinv
     enddo
     do ll=1,n
        if(ll.ne.icol)then
           dum=a(ll,icol)
           a(ll,icol)=0.
           do l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
           enddo
           do l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
           enddo
        endif
     enddo
  enddo
  do l=n,1,-1
     if(indxr(l).ne.indxc(l))then
        do k=1,n
           dum=a(k,indxr(l))
           a(k,indxr(l))=a(k,indxc(l))
           a(k,indxc(l))=dum
        enddo
     endif
  end do
  return
END SUBROUTINE gaussj

!=======================================================================
function cool_bre(ispec,T)
!=======================================================================
  implicit none
  integer::ispec
  real(kind=8)   ::T,cool_bre
  if(ispec==HI  )cool_bre = 1.42D-27*sqrt(T)*(1.1D0+0.34D0*exp(-(5.5D0-log10(T))**2 /3.D0))
  if(ispec==HEI )cool_bre = 1.42D-27*sqrt(T)*(1.1D0+0.34D0*exp(-(5.5D0-log10(T))**2 /3.D0))
  if(ispec==HEII)cool_bre = 5.68D-27*sqrt(T)*(1.1D0+0.34D0*exp(-(5.5D0-log10(T))**2 /3.D0))
  return
end function cool_bre

!=======================================================================
function cool_exc(ispec,T)
!=======================================================================
  implicit none
  integer::ispec
  real(kind=8)   ::T,cool_exc,T5
  T5=1.d-5*T
  if(ispec==HI  )cool_exc = 7.50D-19/(1.+sqrt(T5))              *exp(-118348.D0/T)
!  Terme elimine car cubique en nH d'apres Cen
!  if(ispec==HEI )cool_exc = 9.10D-27/(1.+sqrt(T5))/(T**0.1687D0)*exp(-13179.D0/T)
  if(ispec==HEII)cool_exc = 5.54D-17/(1.+sqrt(T5))/(T**0.397D0 )*exp(-473638.D0/T)
  return
end function cool_exc

!=======================================================================
function cool_rec(ispec,T)
!=======================================================================
  implicit none
  integer::ispec
  real(kind=8)   ::T,cool_rec
  real(kind=8)   ::T3, T6
  T3 = 1.d-03*T
  T6 = 1.d-06*T
  if(ispec==HI  )cool_rec = 8.70D-27*SQRT(T)/T3**(0.2D0)/(1.D0+T6**0.7D0)
  if(ispec==HEI )cool_rec = 1.55D-26*T**0.3647D0
  if(ispec==HEII)cool_rec = 3.48D-26*SQRT(T)/T3**(0.2D0)/(1.D0+T6**0.7D0)
  return
end function cool_rec

!=======================================================================
function cool_die(T)
!=======================================================================
  implicit none
  real(kind=8) :: T,cool_die
  cool_die=1.24D-13*T**(-1.5D0)*exp(-470000.D0/T)*(1.D0+0.3D0*exp(-94000.D0/T))
  return
end function cool_die

!=======================================================================
function taux_rec(ispec,T)
!=======================================================================
  implicit none
  integer::ispec
  real(kind=8)   ::T,taux_rec
  real(kind=8)   ::T3, T6
  T3 = 1.d-03*T
  T6 = 1.d-06*T
  if(ispec==HI  )taux_rec = dumfac_rec*8.40e-11/SQRT(T)/T3**(0.2)/(1.+T6**0.7)
  if(ispec==HEI )taux_rec = 1.50e-10/T**0.6353+taux_die(T)
  if(ispec==HEII)taux_rec = 3.36e-10/SQRT(T)/T3**(0.2)/(1.+T6**0.7)
  return
end function taux_rec

!=======================================================================
function taux_die(T)
!=======================================================================
  implicit none
  real(kind=8) :: T,taux_die
  taux_die=1.9D-3*T**(-1.5D0)*exp(-470000.D0/T)*(1.D0+0.3D0*exp(-94000.D0/T))
  return
end function taux_die

!=======================================================================
function cool_ion(ispec,T)
!=======================================================================
  implicit none
  integer::ispec
  real(kind=8)   ::T,cool_ion
  real(kind=8)   ::T5
  T5 = 1.d-05*T
  if(ispec==HI  )cool_ion = dumfac_ion*1.27D-21*SQRT(T)/(1.+SQRT(T5))*EXP(-157809.1D0/T)
  if(ispec==HEI )cool_ion = dumfac_ion*9.38D-22*SQRT(T)/(1.+SQRT(T5))*EXP(-285335.4D0/T)
  if(ispec==HEII)cool_ion = dumfac_ion*4.95D-22*SQRT(T)/(1.+SQRT(T5))*EXP(-631515.0D0/T)
  return
end function cool_ion

!=======================================================================
function cool_com(T,aexp)
!=======================================================================
  implicit none
  real(kind=8) ::T,aexp,cool_com
  cool_com=5.406D-36*T/aexp**4 
  return
end function cool_com

!=======================================================================
function heat_com(T,aexp)
!=======================================================================
  implicit none
  real(kind=8) ::T,aexp,heat_com
  heat_com=5.406D-36*2.726D0/aexp**5
  return
end function heat_com

!=======================================================================
function taux_ion(ispec,T)
!=======================================================================
  implicit none
  integer::ispec
  real(kind=8)   :: T,taux_ion
  real(kind=8)   :: T5
  T5 = 1.d-05*T
  if(ispec==HI  )taux_ion = dumfac_ion*5.85D-11*SQRT(T)/(1.+SQRT(T5))*EXP(-157809.1D0/T)
  if(ispec==HEI )taux_ion = dumfac_ion*2.38D-11*SQRT(T)/(1.+SQRT(T5))*EXP(-285335.4D0/T)
  if(ispec==HEII)taux_ion = dumfac_ion*5.68D-12*SQRT(T)/(1.+SQRT(T5))*EXP(-631515.0D0/T)
  return
end function taux_ion

!=======================================================================
function J_nu(e,J0)
!=======================================================================
  implicit none
  real(kind=8) :: e,J_nu,e_L,J0,Jloc
  Jloc = max(J0,J0min) 
  e_L  = 13.598*eV
  J_nu = Jloc*(e_L/e)
  return
end function J_nu

!=======================================================================
function sigma_rad(e,ispec)
!=======================================================================
  implicit none
  integer::ispec
  real(kind=8)   ::sigma_rad,e,e_i,xxx,alph
  if(ispec==HI  )e_i = 13.598D0*eV
  if(ispec==HEI )e_i = 24.587D0*eV
  if(ispec==HEII)e_i = 54.416D0*eV
  xxx = e/e_i
  alph = sqrt(xxx-1.0d0)
  if(ispec==HI  )sigma_rad = 6.30D-18/xxx**4*exp(4.D0-4.D0*atan(alph)/alph) &
       &                             /(1.D0-exp(-twopi/alph))
  if(ispec==HEI )sigma_rad = 7.42D-18*(1.66D0/xxx**2.05D0-0.66D0/xxx**3.05D0)
  if(ispec==HEII)sigma_rad = 1.58D-18/xxx**4*exp(4.D0-4.D0*atan(alph)/alph) &
       &                             /(1.D0-exp(-twopi/alph))
  return
end function sigma_rad

!=======================================================================
function taux_rad(ispec,J0)
!=======================================================================
  implicit none  
  integer::ispec
  real(kind=8) :: J0,taux_rad,e_i,e,de,error,integ
  if(ispec==HI  )e_i = 13.598D0*eV
  if(ispec==HEI )e_i = 24.587D0*eV
  if(ispec==HEII)e_i = 54.416D0*eV
  integ = 0.0d0
  e = e_i
  de = e/100.D0
  error = 1.D0
  do while(error>1.d-6)
     e = e + de
     de = e/100.D0
     error = 2.0d0*twopi*J_nu(e,J0)*sigma_rad(e,ispec)*de/e
     integ = integ + error
     error = error/abs(integ)
  end do
  taux_rad = integ/hplanck
  return
end function taux_rad

!=======================================================================
function taux_rad_madau(ispec,z)
!=======================================================================
  implicit none
  integer :: ispec
  real(kind=8) :: z,taux_rad_madau,tt
  if (z < 15.d0) then
     if (ispec==HI  ) taux_rad_madau=normfacJ0*exp(-31.04D0+2.795D0*z-0.5589D0*z**2)
     if (ispec==HEI ) taux_rad_madau=normfacJ0*exp(-31.08D0+2.822D0*z-0.5664D0*z**2)
     if (ispec==HEII) taux_rad_madau=normfacJ0*exp(-34.30D0+1.826D0*z-0.3899D0*z**2)
  else
     taux_rad_madau=0.d0
  endif
  tt=taux_rad_theuns(ispec,J0min)
  if (taux_rad_madau < tt) taux_rad_madau=tt
  return
end function taux_rad_madau

!=======================================================================
function taux_rad_weinberg(ispec,z)
!=======================================================================
  implicit none
  integer :: ispec
  real(kind=8) :: z,taux_rad_weinberg
  real(kind=8) :: am1log,am1logmin,am1step,weight1,weight2
  integer :: N1,N2,iweinb
  am1log=log10((1.D0+z))
  am1logmin=table_weinberg(1,1)
  am1step=table_weinberg(1,2)-table_weinberg(1,1)
  N1=int((am1log-am1logmin)/am1step)
  N2=N1+1
  if (N1 <= 0) then
     N1=N1+1
     N2=N2+1
  elseif (N1 >= Nweinberg) then
     taux_rad_weinberg=taux_rad_theuns(ispec,J0min)
     return
  endif
  weight2=(am1log-table_weinberg(1,N1))/am1step
  weight1=1.d0-weight2  
  if (ispec==HI  ) iweinb=2
  if (ispec==HEI ) iweinb=3
  if (ispec==HEII) iweinb=4  
  taux_rad_weinberg=normfacJ0*exp(log(table_weinberg(iweinb,N1))*weight1 &
 &                     +log(table_weinberg(iweinb,N2))*weight2)
  return
end function taux_rad_weinberg

!=======================================================================
function taux_rad_weinbergint(ispec,z)
!=======================================================================
  implicit none
  integer :: ispec,i,iweinb
  real(kind=8) :: z,zz,taux_rad_weinbergint,hh,tt
  if (z < 8.5d0) then
     if (ispec==HI  ) iweinb=1
     if (ispec==HEI ) iweinb=2
     if (ispec==HEII) iweinb=3
     hh=0.d0
     zz=max(z,1.0d-15)
     do i=1,Norderweinberg
        hh=hh+coefweinberg(i,iweinb)*zz**(i-1)
     enddo
     taux_rad_weinbergint=normfacJ0*exp(hh)
  else
     taux_rad_weinbergint=0.d0
  endif
  tt=taux_rad_theuns(ispec,J0min)
  if (taux_rad_weinbergint < tt) taux_rad_weinbergint=tt
  return
end function taux_rad_weinbergint

!=======================================================================
function taux_rad_theuns(ispec,J0)
!=======================================================================
  implicit none
  integer :: ispec
  real(kind=8) :: J0,taux_rad_theuns

  if (ispec==HI  ) taux_rad_theuns=1.26D10*J0/(3.D0+alpha)
  if (ispec==HEI ) taux_rad_theuns=1.48D10*J0*0.553D0**alpha &
                     & *(1.66D0/(alpha+2.05D0)-0.66D0/(alpha+3.05D0))
  if (ispec==HEII) taux_rad_theuns=3.34D9*J0*0.249D0**alpha/(3.D0+alpha)
  return
end function taux_rad_theuns

!=======================================================================
function heat_rad(ispec,J0)
!=======================================================================
  implicit none  
  integer::ispec
  real(kind=8) :: J0,heat_rad,e_i,e,de,error,integ
  if(ispec==HI  )e_i = 13.598D0*eV
  if(ispec==HEI )e_i = 24.587D0*eV
  if(ispec==HEII)e_i = 54.416D0*eV
  integ = 0.0d0
  e = e_i
  de = e/100.D0
  error = 1.D0
  do while(error>1.d-6)
     e = e + de
     de = e/100.D0
     error = 2.0d0*twopi*J_nu(e,J0)*sigma_rad(e,ispec)*(e/e_i-1.D0)*de/e
     integ = integ + error
     error=error/abs(integ)
  end do
  heat_rad = integ/hplanck*e_i
  return
end function heat_rad
  
!=======================================================================
function heat_rad_madau(ispec,z)
!=======================================================================
  implicit none
  integer :: ispec
  real(kind=8) :: z,heat_rad_madau,tt
  if (z < 15.d0) then
     if (ispec==HI  ) heat_rad_madau=normfacJ0*exp(-56.62D0+2.788D0*z-0.5594D0*z**2)
     if (ispec==HEI ) heat_rad_madau=normfacJ0*exp(-56.06D0+2.800D0*z-0.5532D0*z**2)
     if (ispec==HEII) heat_rad_madau=normfacJ0*exp(-58.67D0+1.888D0*z-0.3947D0*z**2)
  else
     heat_rad_madau=0.d0
  endif
  tt=heat_rad_theuns(ispec,J0min)
  if (heat_rad_madau < tt) heat_rad_madau=tt
  return
end function heat_rad_madau

!=======================================================================
function heat_rad_weinberg(ispec,z)
!=======================================================================
  implicit none
  integer :: ispec
  real(kind=8) :: z,heat_rad_weinberg
  real(kind=8) :: am1log,am1logmin,am1step,weight1,weight2
  integer :: N1,N2,iweinb
  am1log=log10((1.D0+z))
  am1logmin=table_weinberg(1,1)
  am1step=table_weinberg(1,2)-table_weinberg(1,1)
  N1=int((am1log-am1logmin)/am1step)
  N2=N1+1
  if (N1 <= 0) then
     N1=N1+1
     N2=N2+1
  elseif (N1 >= Nweinberg) then
     heat_rad_weinberg=heat_rad_theuns(ispec,J0min)
     return
  endif
  weight2=(am1log-table_weinberg(1,N1))/am1step
  weight1=1.d0-weight2  
  if (ispec==HI  ) iweinb=5
  if (ispec==HEI ) iweinb=6
  if (ispec==HEII) iweinb=7  
  heat_rad_weinberg=normfacJ0*exp(log(table_weinberg(iweinb,N1))*weight1 &
 &                            +log(table_weinberg(iweinb,N2))*weight2)
  return
end function heat_rad_weinberg

!=======================================================================
function heat_rad_weinbergint(ispec,z)
!=======================================================================
  implicit none
  integer :: ispec,i,iweinb
  real(kind=8) :: z,zz,heat_rad_weinbergint,hh,tt
  if (z < 8.5d0) then
     if (ispec==HI  ) iweinb=4
     if (ispec==HEI ) iweinb=5
     if (ispec==HEII) iweinb=6
     hh=0.d0
     zz=max(z,1.0d-15)
     do i=1,Norderweinberg
        hh=hh+coefweinberg(i,iweinb)*zz**(i-1)
     enddo
     heat_rad_weinbergint=normfacJ0*exp(hh)
  else
     heat_rad_weinbergint=0.d0
  endif
  tt=heat_rad_theuns(ispec,J0min)
  if (heat_rad_weinbergint < tt) heat_rad_weinbergint=tt
  return
end function heat_rad_weinbergint

!=======================================================================
function heat_rad_theuns(ispec,J0)
!=======================================================================
  implicit none
  integer :: ispec
  real(kind=8) :: J0,heat_rad_theuns
  if (ispec==HI  ) heat_rad_theuns=(2.91D-1*J0/(2.D0+alpha))/(3.D0+alpha)
  if (ispec==HEI ) heat_rad_theuns=5.84D-1*J0*0.553D0**alpha* &
                 & (1.66D0/(alpha+1.05D0)-2.32D0/(alpha+2.05D0)+0.66D0/(alpha+3.05D0))
  if (ispec==HEII) heat_rad_theuns=(2.92D-1*J0*0.249D0**alpha/(2.D0+alpha))/(3.D0+alpha)
  return
end function heat_rad_theuns


!=======================================================================
function HsurH0(z,omega0,omegaL,OmegaR)
!=======================================================================
  implicit none
  real(kind=8) :: HsurH0,z,omega0,omegaL,omegaR
  HsurH0=sqrt(Omega0*(1.d0+z)**3+OmegaR*(1.d0+z)**2+OmegaL)
end function HsurH0

end module cooling_module

