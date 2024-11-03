	program code

	implicit none

	real*8, parameter :: time=1000.d0, t0=0.d0
	real*8, parameter :: er=1.d-6
	real*8, parameter :: vton=600.d0		! uM to um^3 conversion factor

	real*8,parameter :: kpe = 10.d0/vton, kpae = 100.d0/vton
	real*8,parameter :: kpa0 = 1.d0/vton, kpa1 = 2000.d0/vton, kma = 5.d-3
	real*8,parameter :: v0 = 5.d-3

!---------------------------------
	real*8, parameter :: dv=10.d0*2.d-5		!162.d-7 volume of SPD5 molecule (doi 10.1110/ps.04688204)
	real*8,parameter :: rhoe = 0.1d0		! enzyme concentration in uM 
	real*8,parameter :: rho = 0.1d0			! concentration in uM (https://doi.org/10.1016/j.cell.2016.08.006)
!---------------------------------
	integer, parameter :: V=5000			! cell volume (http://book.bionumbers.org/how-big-is-a-human-cell/)

	integer, parameter :: Ne=nint(V*rhoe*vton)
	integer, parameter :: N=nint(V*rho*vton)
	integer, parameter :: nstep=200

	real*8, parameter :: delt=(time-t0)/nstep

	real*8 :: t,r1,r2,pl,ph
	real*8 :: tau
	real*8 :: beta1,beta2,beta3,beta4,beta_sum
	real*8 :: peta1,peta2,peta3,peta4
	real*8 :: alpha0
	real*8 :: dum, msize

	integer :: me, m, ms

	integer :: i,j,k,nt,nk
	integer :: le
	integer :: l1,l2

	integer :: kk, cnt

!------------------------------------------------------

	open(unit=100,file='data/ab.txt',status='unknown')
	open(unit=200,file='data/l.txt',status='unknown')

!------------------------------------------------------

	write(*,*)'Total subunits =',N

	t=0
	nk=0
	msize=0.d0

!----------- initialisation -----------

	
	le=0
	l1=nint(v0/dv)
	l2=nint(v0/dv)
	ms=0
	me=Ne-le
	m=N-l1-l2-ms

	if(t0.le.er)then
	write(200,*)t,(l1)*dv,(l2)*dv,m,ms,le,me
	endif

!------------------------------------------------------
!------------------- time loop ------------------------

	do while (t.le.time)

!------------------ propensity sum ---------------------

	beta1 = kpae*m*le/V
	beta2 = kpe*(l1+l2)*me/V

	beta3 = kpa1*ms/V			
	beta4 = kpa1*ms/V 

	peta1 = kma*l1
	peta2 = kpa0*m/V			
	peta3 = kma*l2
	peta4 = kpa0*m/V

	beta_sum = beta1+beta2+beta3+beta4 + peta1+peta2+peta3+peta4

	alpha0 = beta_sum

	if(alpha0.lt.0)then
	write(*,*) t, alpha0, beta1+beta2+beta3+beta4, peta1+peta2+peta3+peta4
	endif

!----------------- next reaction time: tau --------------

	r1 = rand()
	if(r1.lt.1.d-8)then
	r1 = rand()
	endif

	tau = (1.d0/alpha0)*log(1.d0/r1)

!----------------------- reaction ------------------------

	r2 = rand()

!---------------- growth and decay ------------------

!------------------ m+e -> m* ----------------------
	pl=0.d0
	ph=beta1/alpha0
	if (r2.ge.pl.and.r2.lt.ph) then
	ms = ms+1
	le = le-1
	m = m -1
	endif
!------------------ Sn+e -> Sn+e* ----------------------
	pl=beta1/alpha0
	ph=(beta1+beta2)/alpha0
	if (r2.ge.pl.and.r2.lt.ph) then
	le = le+1
	me = me-1
	endif
!------------------ Sn+m* -> Sn+1 ----------------------
	pl=(beta1+beta2)/alpha0
	ph=(beta1+beta2+beta3)/alpha0
	if (r2.ge.pl.and.r2.lt.ph) then
	l1 = l1+1
	ms = ms-1
	me = me+1
	endif

	pl=(beta1+beta2+beta3)/alpha0
	ph=(beta1+beta2+beta3+beta4)/alpha0
	if (r2.ge.pl.and.r2.lt.ph) then
	l2 = l2+1
	ms = ms-1
	me = me+1
	endif

!------------------ non enzymatic PCM ----------------------

!------------------ Sn ->< Sn-1 + m ----------------------
	pl=(beta1+beta2+beta3+beta4)/alpha0
	ph=(beta1+beta2+beta3+beta4+peta1)/alpha0
	if (r2.ge.pl.and.r2.lt.ph.and.l1.gt.er) then
	l1 = l1-1
	m = m+1
	endif

	pl=(beta1+beta2+beta3+beta4+peta1)/alpha0
	ph=(beta1+beta2+beta3+beta4+peta1+peta2)/alpha0
	if (r2.ge.pl.and.r2.lt.ph) then
	l1 = l1+1
	m = m-1
	endif

	pl=(beta1+beta2+beta3+beta4+peta1+peta2)/alpha0
	ph=(beta1+beta2+beta3+beta4+peta1+peta2+peta3)/alpha0
	if (r2.ge.pl.and.r2.lt.ph.and.l2.gt.er) then
	l2 = l2-1
	m = m+1
	endif

	pl=(beta1+beta2+beta3+beta4+peta1+peta2+peta3)/alpha0
	ph=(beta1+beta2+beta3+beta4+peta1+peta2+peta3+peta4)/alpha0
	if (r2.ge.pl.and.r2.lt.ph) then
	l2 = l2+1
	m = m-1
	endif


	
	if((t-t0)/(delt+nk*delt).ge.1.d0) then
	write(200,*)t,l1*dv,l2*dv,m,ms,le,me
	write(100,*)t,l1*dv-l2*dv
	msize = msize + (l1+l2)*dv/2.d0
	nk=nk+1
	endif


	t = t + tau


	enddo

!------------------ time loop ended -----------------------


	write(*,*) V, msize/(nk+1), (l1+l2)*dv/2.d0

	stop
	end program code








