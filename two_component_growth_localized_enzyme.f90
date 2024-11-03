	program code

	implicit none

	real*8, parameter :: time=1200.d0, t0=0.d0
	real*8, parameter :: er=1.d-6
	real*8, parameter :: vton=600.d0							! uM to um^3 conversion factor

	real*8,parameter :: kpa = 10.d0/vton, kma = 5.d-3					!kpa = 1000.d0/vton
	real*8,parameter :: kpae = 5000.d0/vton, kpeb = 1000.d0/vton
	real*8,parameter :: kpb0 = 0.5d0/vton, kpb1 = 10000.d0/vton, kmb = 5.d-3		!kpb = 40.d0/vton
	real*8,parameter :: kmib = 0.01d0

!---------------------------------
	real*8, parameter :: dv=10*2.d-5			! 162.d-7 volume of SPD5 molecule (doi 10.1110/ps.04688204)
	real*8,parameter :: rhoa = 0.25d0			! concentration in uM (https://doi.org/10.1016/j.cell.2016.08.006)
	real*8,parameter :: rhob = 0.5d0			! concentration in uM (https://doi.org/10.1016/j.cell.2016.08.006)
	real*8,parameter :: rhoe = 0.01d0			! concentration in uM
	real*8,parameter :: v0 = 5.d-3
!---------------------------------
	integer, parameter :: V=5000				! cell volume (http://book.bionumbers.org/how-big-is-a-human-cell/)

	integer, parameter :: Na=nint(V*rhoa*vton)
	integer, parameter :: Nb=nint(V*rhob*vton)
	integer, parameter :: Ne=nint(V*rhoe*vton)
	integer, parameter :: nstep=2000

	real*8, parameter :: delt=(time-t0)/nstep

	real*8 :: t,r1,r2,pl,ph
	real*8 :: tau
	real*8 :: beta1,beta2,beta3,beta4
	real*8 :: beta5,beta6,beta7,beta8,beta_sum
	real*8 :: peta00,peta01,peta02,peta03,peta1,peta2,peta3,peta4,peta_sum
	real*8 :: alpha0
	real*8 :: dum, xf

	integer :: ma, mb, me

	integer :: i,j,k,nt,nk
	integer :: la1,la2
	integer :: lb1,lb2, lbi1, lbi2
	integer :: le1,le2
	integer :: leb1, leb2

	integer :: kk, cnt

!------------------------------------------------------

	open(unit=100,file='data/ab.txt',status='unknown')
	open(unit=200,file='data/l.txt',status='unknown')

!------------------------------------------------------

	t=0
	nk=0
	xf = 0.d0	!kpa/kpb

!----------- initialisation -----------

	
	la1=nint((v0+0.01d0)/dv)
	la2=nint(v0/dv)
	lb1=nint(xf*la1)
	lb2=nint(xf*la2)
	le1=0
	le2=0
	leb1=0
	leb2=0
	lbi1=0
	lbi2=0
	ma=Na-la1-la2
	mb=Nb-lb1-lb2
	me=Ne-le1-le2

	write(200,*)t,(la1+lb1)*dv,(la2+lb2)*dv,le1,le2,leb1,leb2
	write(100,*)t,la1*dv,la2*dv,lb1*dv,lb2*dv

!------------------------------------------------------
!------------------- time loop ------------------------

	do while (t.le.time)

!------------------ propensity sum ---------------------

	beta1 = kma*la1
	beta2 = kpa*ma/V

	beta3 = kma*la2
	beta4 = kpa*ma/V

	beta5 = kpae*la1*me/V
	beta6 = kpae*la2*me/V

	beta7 = max(kpeb*le1*lbi1/V,0.d0)
	beta8 = max(kpeb*le2*lbi2/V,0.d0)

	peta00 = kpb0*la1*mb/V
	peta01 = kpb0*la2*mb/V
	peta02 = kmib*lbi1
	peta03 = kmib*lbi2
	peta1 = kmb*lb1
	peta2 = max(kpb1*leb1/V,0.d0)			
	peta3 = kmb*lb2
	peta4 = max(kpb1*leb2/V,0.d0)			

	alpha0 = beta1+beta2+beta3+beta4 + beta5+beta6+beta7+beta8 + peta1+peta2+peta3+peta4 + peta00+peta01+peta02+peta03

	if(alpha0.lt.0)then
	write(*,*) "ERROR!!"
	write(*,*) t, alpha0, beta1+beta2+beta3+beta4, peta1+peta2+peta3+peta4
	endif

!----------------- next reaction time: tau --------------

	r1 = rand()
	if(r1.lt.1.d-8)then
	r1 = rand()
	endif

	tau = (1.d0/alpha0)*log(1.d0/r1)

	if(tau.lt.0)then
	write(*,*) "ERROR!!"
	write(*,*) t, alpha0, tau
	endif

!----------------------- reaction ------------------------

	r2 = rand()

!---------------- growth and decay ------------------

!------------------ scaffold ----------------------
	pl=0.d0
	ph=beta1/alpha0
	if (r2.ge.pl.and.r2.lt.ph.and.la1.gt.er) then
	la1 = la1-1
	ma = ma+1
	endif

	pl=beta1/alpha0
	ph=(beta1+beta2)/alpha0
	if (r2.ge.pl.and.r2.lt.ph) then
	la1 = la1+1
	ma = ma-1
	endif

	pl=(beta1+beta2)/alpha0
	ph=(beta1+beta2+beta3)/alpha0
	if (r2.ge.pl.and.r2.lt.ph.and.la2.gt.er) then
	la2 = la2-1
	ma = ma+1
	endif

	pl=(beta1+beta2+beta3)/alpha0
	ph=(beta1+beta2+beta3+beta4)/alpha0
	if (r2.ge.pl.and.r2.lt.ph) then
	la2 = la2+1
	ma = ma-1
	endif

!----------------- Enzyme --------------------

	pl=(beta1+beta2+beta3+beta4)/alpha0
	ph=(beta1+beta2+beta3+beta4+beta5)/alpha0
	if (r2.ge.pl.and.r2.lt.ph.and.la1.gt.er) then
	le1 = le1+1
	me = me-1
	endif

	pl=(beta1+beta2+beta3+beta4+beta5)/alpha0
	ph=(beta1+beta2+beta3+beta4+beta5+beta6)/alpha0
	if (r2.ge.pl.and.r2.lt.ph.and.la1.gt.er) then
	le2 = le2+1
	me = me-1
	endif

	pl=(beta1+beta2+beta3+beta4+beta5+beta6)/alpha0
	ph=(beta1+beta2+beta3+beta4+beta5+beta6+beta7)/alpha0
	if (r2.ge.pl.and.r2.lt.ph.and.la1.gt.er) then
	leb1 = leb1+1
	lbi1 = lbi1-1
	le1 = le1-1
	endif

	pl=(beta1+beta2+beta3+beta4+beta5+beta6+beta7)/alpha0
	ph=(beta1+beta2+beta3+beta4+beta5+beta6+beta7+beta8)/alpha0
	if (r2.ge.pl.and.r2.lt.ph.and.la1.gt.er) then
	leb2 = leb2+1
	lbi2 = lbi2-1
	le2 = le2-1
	endif

	beta_sum = beta1+beta2+beta3+beta4+beta5+beta6+beta7+beta8

!------------------ PCM ----------------------


	pl=(beta_sum)/alpha0
	ph=(beta_sum+peta1)/alpha0
	if (r2.ge.pl.and.r2.lt.ph.and.la1.gt.er) then
	lb1 = lb1-1
	mb = mb+1
	endif

	pl=(beta_sum+peta1)/alpha0
	ph=(beta_sum+peta1+peta2)/alpha0
	if (r2.ge.pl.and.r2.lt.ph) then
	lb1 = lb1+1
	leb1 = leb1-1
	me = me+1
	endif

	pl=(beta_sum+peta1+peta2)/alpha0
	ph=(beta_sum+peta1+peta2+peta3)/alpha0
	if (r2.ge.pl.and.r2.lt.ph.and.la2.gt.er) then
	lb2 = lb2-1
	mb = mb+1
	endif

	pl=(beta_sum+peta1+peta2+peta3)/alpha0
	ph=(beta_sum+peta1+peta2+peta3+peta4)/alpha0
	if (r2.ge.pl.and.r2.lt.ph) then
	lb2 = lb2+1
	leb2 = leb2-1
	me = me+1
	endif

	pl=(beta_sum+peta1+peta2+peta3+peta4)/alpha0
	ph=(beta_sum+peta1+peta2+peta3+peta4+peta00)/alpha0
	if (r2.ge.pl.and.r2.lt.ph) then
	lbi1 = lbi1+1
	mb = mb-1
	endif

	pl=(beta_sum+peta1+peta2+peta3+peta4+peta00)/alpha0
	ph=(beta_sum+peta1+peta2+peta3+peta4+peta00+peta01)/alpha0
	if (r2.ge.pl.and.r2.lt.ph) then
	lbi2 = lbi2+1
	mb = mb-1
	endif

	peta_sum=beta_sum+peta1+peta2+peta3+peta4+peta00+peta01

	pl=(peta_sum)/alpha0
	ph=(peta_sum+peta02)/alpha0
	if (r2.ge.pl.and.r2.lt.ph) then
	lbi1 = lbi1-1
	mb = mb+1
	endif

	pl=(peta_sum+peta02)/alpha0
	ph=(peta_sum+peta02+peta03)/alpha0
	if (r2.ge.pl.and.r2.lt.ph) then
	lbi2 = lbi2-1
	mb = mb+1
	endif

!------------------------------------------


	
	if((t-t0)/(delt+nk*delt).ge.1.d0) then
	write(200,*)t,(la1+lb1)*dv,(la2+lb2)*dv,le1,le2,leb1,leb2,lbi1,lbi2
	write(100,*)t,la1*dv,la2*dv,lb1*dv,lb2*dv
	nk=nk+1
	endif


	if(la1+la2+ma.ne.Na)then
	write(*,*) "ERROR!! - a not conserved"
	GOTO 100
	endif

	if(lb1+lb2+mb+leb1+leb2+lbi1+lbi2.ne.Nb)then
	write(*,*) "ERROR!! - b not conserved"
	GOTO 100
	endif

	if(le1+le2+me+leb1+leb2.ne.Ne)then
	write(*,*) "ERROR!! - E not conserved"
	GOTO 100
	endif

	t = t + tau

	enddo

!------------------ time loop ended -----------------------

	if (t.ge.time) then
	GOTO 200
	endif

100	write(*,*) "time=",t, "Code stopped due to error"

200	write(*,*) (la1+lb1+lbi1)*dv,(la2+lb2+lbi2)*dv

	stop
	end program code








