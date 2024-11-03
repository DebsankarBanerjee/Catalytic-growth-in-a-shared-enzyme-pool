	program code

	implicit none

	real*8, parameter :: time=1000.d0, t0=0.d0
	real*8, parameter :: er=1.d-6
	real*8, parameter :: sdur=1.d3
	real*8, parameter :: vton=600.d0		! uM to um^3 conversion factor

	real*8,parameter :: c0=1.d0*600.d0		! 0.09d0*600.d0
	real*8,parameter :: c1=0.001d0*600.d0		! 0.9d-3*600.d0

	real*8,parameter :: kp1 = c0/vton, ka1 = c1/vton, km1 = 5.d-3
	real*8,parameter :: kp2 = c0/vton, ka2 = c1/vton, km2 = 5.d-3


	real*8,parameter :: rho = 1.d0/30.d0		! (1.d0/60.d0) concentration in uM (https://doi.org/10.1016/j.cell.2016.08.006)

	real*8, parameter :: dv=10*2.d-5		! (162.d-7) volume of SPD5 molecule (doi 10.1110/ps.04688204)
	real*8, parameter :: v0 = 5.d-3

!---------------------------------
	integer, parameter :: V=5000			! cell volume (http://book.bionumbers.org/how-big-is-a-human-cell/)

	integer, parameter :: N=nint(V*rho*vton)
	integer, parameter :: nstep=12500

	real*8, parameter :: delt=(time-t0)/nstep

	real*8 :: t,r1,r2,pl,ph
	real*8 :: tau
	real*8 :: beta1,beta2,beta3,beta4,beta_sum
	real*8 :: alpha0,alpha_R,alpha_L,alpha_R_sum,alpha_L_sum
	real*8 :: dum, dumran

	integer :: m

	integer :: i,j,k,nt,np,nk
	integer :: l1,l2

	integer :: kk, cnt
	INTEGER :: seed(12), scount

!------------------------------------------------------

	open(unit=100,file='data/m.txt',status='unknown')
	open(unit=200,file='data/l.txt',status='unknown')

!------------------------------------------------------

	write(*,*)'Total subunits =',N

	t=0
	nk=0
  	scount = 1

!----------- initialisation -----------

	dumran = 6772				

  	seed= 9878772+nint(dumran*1091)

        CALL RANDOM_SEED (PUT=seed)

	
	l1=nint((v0+0.01d0)/dv)
	l2=nint(v0/dv)
	m=N-l1-l2

	write(200,*)t,l1*dv,l2*dv,m

!------------------------------------------------------
!------------------- time loop ------------------------

	do while (t.le.time)

!------------------ propensity sum ---------------------

	beta1 = km1*l1
	beta2 = (kp1 + ka1*l1)*m/V
	beta3 = km2*l2
	beta4 = (kp2 + ka2*l2)*m/V
	beta_sum = beta1+beta2+beta3+beta4

	alpha0 = beta_sum

!----------------- next reaction time: tau --------------

	CALL RANDOM_NUMBER(r1)
	if(r1.lt.1.d-8)then		!and hope the next one is not zero too :)
	CALL RANDOM_NUMBER(r1)
	endif

	tau = (1.d0/alpha0)*log(1.d0/r1)

!----------------------- reaction ------------------------

	CALL RANDOM_NUMBER(r2)

!---------------- growth and decay ------------------
	pl=0.d0
	ph=beta1/alpha0
	if (r2.ge.pl.and.r2.lt.ph.and.l1.gt.er) then
	l1 = l1-1
	m = m+1
	endif

	pl=beta1/alpha0
	ph=(beta1+beta2)/alpha0
	if (r2.ge.pl.and.r2.lt.ph) then
	l1 = l1+1
	m = m-1
	endif

	pl=(beta1+beta2)/alpha0
	ph=(beta1+beta2+beta3)/alpha0
	if (r2.ge.pl.and.r2.lt.ph.and.l2.gt.er) then
	l2 = l2-1
	m = m+1
	endif

	pl=(beta1+beta2+beta3)/alpha0
	ph=(beta1+beta2+beta3+beta4)/alpha0
	if (r2.ge.pl.and.r2.lt.ph) then
	l2 = l2+1
	m = m-1
	endif
	
	if((t-t0)/(delt+nk*delt).ge.1.d0) then
	write(200,*)t,l1*dv,l2*dv,m
	nk=nk+1
	endif

	t = t + tau

        if(t.ge.scount*sdur) then
        seed = nint(seed/(1.001d0))
	scount = scount+1
        endif


	enddo

!------------------ time loop ended -----------------------


	stop
	end program code








