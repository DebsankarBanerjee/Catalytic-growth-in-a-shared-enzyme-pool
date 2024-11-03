	program code

	implicit none

	real*8, parameter :: D=1.d0, time=2.d3, t0=0.d0, dt_max=0.05d0
	real*8, parameter :: er=1.d-6
	real*8, parameter :: vton=600.d0		! uM to um^3 conversion factor

	real*8,parameter :: kpe = 10.d0/vton, kpae = 100.d0/vton
	real*8,parameter :: kpa0 = 1.d0/vton, kpa1 = 2000.d0/vton, kma = 5.d-3

!------- system size : [0,L] --------------------------
	real*8,parameter :: L=9.d0

	integer, parameter :: mstar=10, mestar=4, ns=int(L), N=ns*ns*ns*mstar, Ne=ns*ns*ns*mestar
!	integer, parameter :: r=nint((L/ns)*(ns/2))
	integer, parameter :: px1=(ns+1)/2, px2=(ns+1)/2, py1=(ns+1)/2, py2=(ns+1)/2
	integer, parameter :: pz1=(ns+1)/2 - 1, pz2=(ns+1)/2 + 1			! position of centrosomes
	integer, parameter :: nstep=100	

	real*8, parameter :: delt=time/nstep
	real*8, parameter :: dx=(L/ns), kd=D/(dx*dx), dy=dx, dz=dx, vbox=dx*dy*dz

	real*8 :: t,r1,r2,pl,ph
	real*8 :: tau
	real*8 :: beta1,beta21,beta22,beta3,beta4,beta1_sum,beta_sum
	real*8 :: peta1,peta2,peta3,peta4
	real*8 :: alpha0
	real*8 :: alpha_R,alpha_L,alpha_R_sum,alpha_L_sum	! left right
	real*8 :: alpha_F,alpha_B,alpha_F_sum,alpha_B_sum	! front back
	real*8 :: alpha_U,alpha_D,alpha_U_sum,alpha_D_sum	! up down
	real*8 :: dum

	integer :: m(1:ns,1:ns,1:ns), me(1:ns,1:ns,1:ns), ms(1:ns,1:ns,1:ns)
	integer :: le(1:ns,1:ns,1:ns)

	integer :: i,j,k,nt,np,nk
	!integer :: m1,m2
	integer :: l1,l2

	integer :: kk, cnt

	INTEGER :: seed(12)

!------------------------------------------------------

	open(unit=100,file='m.txt',status='unknown')
	open(unit=200,file='l.txt',status='unknown')
	open(unit=600,file='test.txt',status='unknown')

!------------------------------------------------------

	t=0

	l1=1
	l2=5
	nk=0
	pl=0.d0; ph=0.d0
	tau=0.d0
	beta1=0.d0; beta21=0.d0; beta21=0.d0; beta3=0.d0; beta4=0.d0; beta_sum=0.d0
	beta1_sum=0.d0
	alpha0=0.d0; alpha_R=0.d0; alpha_L=0.d0; alpha_R_sum=0.d0; alpha_L_sum=0.d0
	alpha_F=0.d0; alpha_B=0.d0; alpha_F_sum=0.d0; alpha_B_sum=0.d0
	alpha_U=0.d0; alpha_D=0.d0; alpha_U_sum=0.d0; alpha_D_sum=0.d0
	dum=0.d0

  	seed= 29878772
	CALL RANDOM_SEED (PUT=seed)

!----------- uniform monomer profile initialisation -----------

	m=mstar
	le=0
	ms=0
	me=mestar
	

	write(200,*)t,l1,l2,sum(m), sum(ms), sum(le), sum(me)

	write(*,*) N, sum(m), mstar, ns, pz1, pz2

!------------------------------------------------------
!------------------- time loop ------------------------

	do while (t.le.time)

!------------------ propensity sum ---------------------

	!----------- reactions --------------

	beta1_sum = kpae*sum(m*le)/vbox		! everywhere
	beta21 = kpe*(l1)*me(px1,py1,pz1)/vbox	! at two boxes
	beta22 = kpe*(l2)*me(px2,py2,pz2)/vbox	! at two boxes

	beta3 = kpa1*ms(px1,py1,pz1)/vbox			
	beta4 = kpa1*ms(px2,py2,pz2)/vbox 

	peta1 = kma*l1
	peta2 = kpa0*m(px1,py1,pz1)/vbox			
	peta3 = kma*l2
	peta4 = kpa0*m(px2,py2,pz2)/vbox

	beta_sum = beta1_sum+beta21+beta22+beta3+beta4 + peta1+peta2+peta3+peta4



	!----------- diffusion --------------
	alpha_R_sum = kd*(sum(m)+sum(me)+sum(ms)+sum(le))
	alpha_L_sum = kd*(sum(m)+sum(me)+sum(ms)+sum(le))
	alpha_F_sum = kd*(sum(m)+sum(me)+sum(ms)+sum(le))
	alpha_B_sum = kd*(sum(m)+sum(me)+sum(ms)+sum(le))
	alpha_U_sum = kd*(sum(m)+sum(me)+sum(ms)+sum(le))
	alpha_D_sum = kd*(sum(m)+sum(me)+sum(ms)+sum(le))


	alpha0 = beta_sum + alpha_R_sum + alpha_L_sum &
                 + alpha_F_sum + alpha_B_sum + alpha_U_sum + alpha_D_sum

!----------------- next reaction time: tau --------------

	CALL RANDOM_NUMBER(r1)
	if(r1.lt.1.d-8)then
	CALL RANDOM_NUMBER(r1)
	endif

	tau = (1.d0/alpha0)*log(1.d0/r1)

!----------------------- reaction-diffusion ------------------------

	CALL RANDOM_NUMBER(r2)

!---------------- growth and decay ------------------

	!-------- s1 activation -----------
	pl=0.d0
	ph=beta1_sum/alpha0


	if (r2.ge.pl.and.r2.lt.ph) then 		!------S1 Activation if
	
	beta1 = 0.d0

	do j=1,ns
	do k=1,ns
	do i=1,ns
	if (r2.ge.(beta1/alpha0)+pl.and.r2.lt. ((beta1+(kpae*m(i,j,k)*le(i,j,k)/vbox)) /alpha0)+pl) then
	ms(i,j,k) = ms(i,j,k)+1
	le(i,j,k) = le(i,j,k)-1
	m(i,j,k) = m(i,j,k) -1
	endif
	beta1 = beta1 + kpae*m(i,j,k)*le(i,j,k)/vbox
	enddo
	enddo
	enddo

	endif						!------S1 Activation if-end



!------------------ Sn(1)+e -> Sn(1)+e* ----------------------
	pl=beta1_sum/alpha0
	ph=(beta1_sum+beta21)/alpha0
	if (r2.ge.pl.and.r2.lt.ph) then
	le(px1,py1,pz1) = le(px1,py1,pz1)+1
	me(px1,py1,pz1) = me(px1,py1,pz1)-1
	endif

!------------------ Sn(2)+e -> Sn(2)+e* ----------------------
	pl=(beta1_sum+beta21)/alpha0
	ph=(beta1_sum+beta21+beta22)/alpha0
	if (r2.ge.pl.and.r2.lt.ph) then
	le(px2,py2,pz2) = le(px2,py2,pz2)+1
	me(px2,py2,pz2) = me(px2,py2,pz2)-1
	endif

!------------------ Sn(1)+m* -> Sn+1(1) ----------------------
	pl=(beta1_sum+beta21+beta22)/alpha0
	ph=(beta1_sum+beta21+beta22+beta3)/alpha0
	if (r2.ge.pl.and.r2.lt.ph) then
	l1 = l1+1
	ms(px1,py1,pz1) = ms(px1,py1,pz1)-1
	me(px1,py1,pz1) = me(px1,py1,pz1)+1
	endif

!------------------ Sn(2)+m* -> Sn+1(2) ----------------------
	pl=(beta1_sum+beta21+beta22+beta3)/alpha0
	ph=(beta1_sum+beta21+beta22+beta3+beta4)/alpha0
	if (r2.ge.pl.and.r2.lt.ph) then
	l2 = l2+1
	ms(px2,py2,pz2) = ms(px2,py2,pz2)-1
	me(px2,py2,pz2) = me(px2,py2,pz2)+1
	endif


!------------------ non enzymatic PCM ----------------------

!------------------ Sn ->< Sn-1 + m ----------------------
	pl=(beta1_sum+beta21+beta22+beta3+beta4)/alpha0
	ph=(beta1_sum+beta21+beta22+beta3+beta4+peta1)/alpha0
	if (r2.ge.pl.and.r2.lt.ph.and.l1.gt.er) then
	l1 = l1-1
	m(px1,py1,pz1) = m(px1,py1,pz1)+1
	endif

	pl=(beta1_sum+beta21+beta22+beta3+beta4+peta1)/alpha0
	ph=(beta1_sum+beta21+beta22+beta3+beta4+peta1+peta2)/alpha0
	if (r2.ge.pl.and.r2.lt.ph) then
	l1 = l1+1
	m(px1,py1,pz1) = m(px1,py1,pz1)-1
	endif

	pl=(beta1_sum+beta21+beta22+beta3+beta4+peta1+peta2)/alpha0
	ph=(beta1_sum+beta21+beta22+beta3+beta4+peta1+peta2+peta3)/alpha0
	if (r2.ge.pl.and.r2.lt.ph.and.l2.gt.er) then
	l2 = l2-1
	m(px2,py2,pz2) = m(px2,py2,pz2)+1
	endif

	pl=(beta1_sum+beta21+beta22+beta3+beta4+peta1+peta2+peta3)/alpha0
	ph=(beta1_sum+beta21+beta22+beta3+beta4+peta1+peta2+peta3+peta4)/alpha0
	if (r2.ge.pl.and.r2.lt.ph) then
	l2 = l2+1
	m(px2,py2,pz2) = m(px2,py2,pz2)-1
	endif



!---------------- diffusion ------------------------

 	!-------- right hopping---------

	pl=(beta_sum)/alpha0
	ph=(beta_sum+alpha_R_sum)/alpha0

	if (r2.ge.pl.and.r2.lt.ph) then 		!------R-hop if

	alpha_R = 0.d0

	!------------ m ---------------
	do j=1,ns
	do k=1,ns

	do i=1,ns-1
	if (r2.ge.(alpha_R/alpha0)+pl.and.r2.lt.((alpha_R+kd*m(i,j,k))/alpha0)+pl) then
	m(i,j,k)=m(i,j,k)-1
	m(i+1,j,k)=m(i+1,j,k)+1
	endif
	alpha_R = alpha_R + kd*m(i,j,k)
	enddo

	if (r2.ge.(alpha_R/alpha0)+pl.and.r2.lt.((alpha_R+kd*m(ns,j,k))/alpha0)+pl) then
	m(ns,j,k)=m(ns,j,k)-1
	m(1,j,k)=m(1,j,k)+1
	endif
	alpha_R = alpha_R+kd*m(ns,j,k)
	
	enddo
	enddo

	!------------ me ---------------
	do j=1,ns
	do k=1,ns

	do i=1,ns-1
	if (r2.ge.(alpha_R/alpha0)+pl.and.r2.lt.((alpha_R+kd*me(i,j,k))/alpha0)+pl) then
	me(i,j,k)=me(i,j,k)-1
	me(i+1,j,k)=me(i+1,j,k)+1
	endif
	alpha_R = alpha_R + kd*me(i,j,k)
	enddo

	if (r2.ge.(alpha_R/alpha0)+pl.and.r2.lt.((alpha_R+kd*me(ns,j,k))/alpha0)+pl) then
	me(ns,j,k)=me(ns,j,k)-1
	me(1,j,k)=me(1,j,k)+1
	endif
	alpha_R = alpha_R+kd*me(ns,j,k)
	
	enddo
	enddo


	!------------ ms ---------------
	do j=1,ns
	do k=1,ns

	do i=1,ns-1
	if (r2.ge.(alpha_R/alpha0)+pl.and.r2.lt.((alpha_R+kd*ms(i,j,k))/alpha0)+pl) then
	ms(i,j,k)=ms(i,j,k)-1
	ms(i+1,j,k)=ms(i+1,j,k)+1
	endif
	alpha_R = alpha_R + kd*ms(i,j,k)
	enddo

	if (r2.ge.(alpha_R/alpha0)+pl.and.r2.lt.((alpha_R+kd*ms(ns,j,k))/alpha0)+pl) then
	ms(ns,j,k)=ms(ns,j,k)-1
	ms(1,j,k)=ms(1,j,k)+1
	endif
	alpha_R = alpha_R+kd*ms(ns,j,k)
	
	enddo
	enddo


	!------------ le ---------------
	do j=1,ns
	do k=1,ns

	do i=1,ns-1
	if (r2.ge.(alpha_R/alpha0)+pl.and.r2.lt.((alpha_R+kd*le(i,j,k))/alpha0)+pl) then
	le(i,j,k)=le(i,j,k)-1
	le(i+1,j,k)=le(i+1,j,k)+1
	endif
	alpha_R = alpha_R + kd*le(i,j,k)
	enddo

	if (r2.ge.(alpha_R/alpha0)+pl.and.r2.lt.((alpha_R+kd*le(ns,j,k))/alpha0)+pl) then
	le(ns,j,k)=le(ns,j,k)-1
	le(1,j,k)=le(1,j,k)+1
	endif
	alpha_R = alpha_R+kd*le(ns,j,k)
	
	enddo
	enddo


	endif						!------R-hop if-end






	pl=(beta_sum+alpha_R_sum)/alpha0
	ph=(beta_sum+alpha_R_sum+alpha_L_sum)/alpha0	
	if (r2.ge.pl.and.r2.lt.ph) then 		!------L-hop if
	!-------- left hopping---------

	alpha_L = 0.d0

	!------------ m ---------------
	do j=1,ns
	do k=1,ns

	if (r2.ge.(alpha_L/alpha0)+pl.and.r2.lt.((alpha_L+kd*m(1,j,k))/alpha0)+pl) then
	m(1,j,k)=m(1,j,k)-1
	m(ns,j,k)=m(ns,j,k)+1
	endif
	alpha_L = alpha_L+kd*m(1,j,k)

	do i=2,ns
	if (r2.ge.(alpha_L/alpha0)+pl.and.r2.lt.((alpha_L+kd*m(i,j,k))/alpha0)+pl) then
	m(i,j,k)=m(i,j,k)-1
	m(i-1,j,k)=m(i-1,j,k)+1
	endif
	alpha_L = alpha_L + kd*m(i,j,k)
	enddo

	enddo
	enddo


	!------------ me ---------------
	do j=1,ns
	do k=1,ns

	if (r2.ge.(alpha_L/alpha0)+pl.and.r2.lt.((alpha_L+kd*me(1,j,k))/alpha0)+pl) then
	me(1,j,k)=me(1,j,k)-1
	me(ns,j,k)=me(ns,j,k)+1
	endif
	alpha_L = alpha_L+kd*me(1,j,k)

	do i=2,ns
	if (r2.ge.(alpha_L/alpha0)+pl.and.r2.lt.((alpha_L+kd*me(i,j,k))/alpha0)+pl) then
	me(i,j,k)=me(i,j,k)-1
	me(i-1,j,k)=me(i-1,j,k)+1
	endif
	alpha_L = alpha_L + kd*me(i,j,k)
	enddo

	enddo
	enddo



	!------------ ms ---------------
	do j=1,ns
	do k=1,ns

	if (r2.ge.(alpha_L/alpha0)+pl.and.r2.lt.((alpha_L+kd*ms(1,j,k))/alpha0)+pl) then
	ms(1,j,k)=ms(1,j,k)-1
	ms(ns,j,k)=ms(ns,j,k)+1
	endif
	alpha_L = alpha_L+kd*ms(1,j,k)

	do i=2,ns
	if (r2.ge.(alpha_L/alpha0)+pl.and.r2.lt.((alpha_L+kd*ms(i,j,k))/alpha0)+pl) then
	ms(i,j,k)=ms(i,j,k)-1
	ms(i-1,j,k)=ms(i-1,j,k)+1
	endif
	alpha_L = alpha_L + kd*ms(i,j,k)
	enddo

	enddo
	enddo


	!------------ le ---------------
	do j=1,ns
	do k=1,ns

	if (r2.ge.(alpha_L/alpha0)+pl.and.r2.lt.((alpha_L+kd*le(1,j,k))/alpha0)+pl) then
	le(1,j,k)=le(1,j,k)-1
	le(ns,j,k)=le(ns,j,k)+1
	endif
	alpha_L = alpha_L+kd*le(1,j,k)

	do i=2,ns
	if (r2.ge.(alpha_L/alpha0)+pl.and.r2.lt.((alpha_L+kd*le(i,j,k))/alpha0)+pl) then
	le(i,j,k)=le(i,j,k)-1
	le(i-1,j,k)=le(i-1,j,k)+1
	endif
	alpha_L = alpha_L + kd*le(i,j,k)
	enddo

	enddo
	enddo

	endif						!------L-hop if-end






	pl=(beta_sum+alpha_R_sum+alpha_L_sum)/alpha0
	ph=(beta_sum+alpha_R_sum+alpha_L_sum+alpha_B_sum)/alpha0	
	if (r2.ge.pl.and.r2.lt.ph) then 		!------B-hop if
	!-------- back hopping---------

	alpha_B = 0.d0

	!------------ m ---------------
	do i=1,ns
	do k=1,ns

	if (r2.ge.(alpha_B/alpha0)+pl.and.r2.lt.((alpha_B+kd*m(i,1,k))/alpha0)+pl) then
	m(i,1,k)=m(i,1,k)-1
	m(i,ns,k)=m(i,ns,k)+1
	endif
	alpha_B = alpha_B+kd*m(i,1,k)

	do j=2,ns
	if (r2.ge.(alpha_B/alpha0)+pl.and.r2.lt.((alpha_B+kd*m(i,j,k))/alpha0)+pl) then
	m(i,j,k)=m(i,j,k)-1
	m(i,j-1,k)=m(i,j-1,k)+1
	endif
	alpha_B = alpha_B + kd*m(i,j,k)
	enddo

	enddo
	enddo


	!------------ me ---------------
	do i=1,ns
	do k=1,ns

	if (r2.ge.(alpha_B/alpha0)+pl.and.r2.lt.((alpha_B+kd*me(i,1,k))/alpha0)+pl) then
	me(i,1,k)=me(i,1,k)-1
	me(i,ns,k)=me(i,ns,k)+1
	endif
	alpha_B = alpha_B+kd*me(i,1,k)

	do j=2,ns
	if (r2.ge.(alpha_B/alpha0)+pl.and.r2.lt.((alpha_B+kd*me(i,j,k))/alpha0)+pl) then
	me(i,j,k)=me(i,j,k)-1
	me(i,j-1,k)=me(i,j-1,k)+1
	endif
	alpha_B = alpha_B + kd*me(i,j,k)
	enddo

	enddo
	enddo


	!------------ ms ---------------
	do i=1,ns
	do k=1,ns

	if (r2.ge.(alpha_B/alpha0)+pl.and.r2.lt.((alpha_B+kd*ms(i,1,k))/alpha0)+pl) then
	ms(i,1,k)=ms(i,1,k)-1
	ms(i,ns,k)=ms(i,ns,k)+1
	endif
	alpha_B = alpha_B+kd*ms(i,1,k)

	do j=2,ns
	if (r2.ge.(alpha_B/alpha0)+pl.and.r2.lt.((alpha_B+kd*ms(i,j,k))/alpha0)+pl) then
	ms(i,j,k)=ms(i,j,k)-1
	ms(i,j-1,k)=ms(i,j-1,k)+1
	endif
	alpha_B = alpha_B + kd*ms(i,j,k)
	enddo

	enddo
	enddo


	!------------ le ---------------
	do i=1,ns
	do k=1,ns

	if (r2.ge.(alpha_B/alpha0)+pl.and.r2.lt.((alpha_B+kd*le(i,1,k))/alpha0)+pl) then
	le(i,1,k)=le(i,1,k)-1
	le(i,ns,k)=le(i,ns,k)+1
	endif
	alpha_B = alpha_B+kd*le(i,1,k)

	do j=2,ns
	if (r2.ge.(alpha_B/alpha0)+pl.and.r2.lt.((alpha_B+kd*le(i,j,k))/alpha0)+pl) then
	le(i,j,k)=le(i,j,k)-1
	le(i,j-1,k)=le(i,j-1,k)+1
	endif
	alpha_B = alpha_B + kd*le(i,j,k)
	enddo

	enddo
	enddo


	endif						!------B-hop if-end




	pl=(beta_sum+alpha_R_sum+alpha_L_sum+alpha_B_sum)/alpha0
	ph=(beta_sum+alpha_R_sum+alpha_L_sum+alpha_B_sum+alpha_F_sum)/alpha0
	if (r2.ge.pl.and.r2.lt.ph) then 		!------F-hop if
	!-------- front hopping---------

	alpha_F = 0.d0

	!------------ m ---------------
	do i=1,ns
	do k=1,ns

	do j=1,ns-1
	if (r2.ge.(alpha_F/alpha0)+pl.and.r2.lt.((alpha_F+kd*m(i,j,k))/alpha0)+pl) then
	m(i,j,k)=m(i,j,k)-1
	m(i,j+1,k)=m(i,j+1,k)+1
	endif
	alpha_F = alpha_F + kd*m(i,j,k)
	enddo

	if (r2.ge.(alpha_F/alpha0)+pl.and.r2.lt.((alpha_F+kd*m(i,ns,k))/alpha0)+pl) then
	m(i,ns,k)=m(i,ns,k)-1
	m(i,1,k)=m(i,1,k)+1
	endif
	alpha_F = alpha_F+kd*m(i,ns,k)
	
	enddo
	enddo


	!------------ me ---------------
	do i=1,ns
	do k=1,ns

	do j=1,ns-1
	if (r2.ge.(alpha_F/alpha0)+pl.and.r2.lt.((alpha_F+kd*me(i,j,k))/alpha0)+pl) then
	me(i,j,k)=me(i,j,k)-1
	me(i,j+1,k)=me(i,j+1,k)+1
	endif
	alpha_F = alpha_F + kd*me(i,j,k)
	enddo

	if (r2.ge.(alpha_F/alpha0)+pl.and.r2.lt.((alpha_F+kd*me(i,ns,k))/alpha0)+pl) then
	me(i,ns,k)=me(i,ns,k)-1
	me(i,1,k)=me(i,1,k)+1
	endif
	alpha_F = alpha_F+kd*me(i,ns,k)
	
	enddo
	enddo


	!------------ ms ---------------
	do i=1,ns
	do k=1,ns

	do j=1,ns-1
	if (r2.ge.(alpha_F/alpha0)+pl.and.r2.lt.((alpha_F+kd*ms(i,j,k))/alpha0)+pl) then
	ms(i,j,k)=ms(i,j,k)-1
	ms(i,j+1,k)=ms(i,j+1,k)+1
	endif
	alpha_F = alpha_F + kd*ms(i,j,k)
	enddo

	if (r2.ge.(alpha_F/alpha0)+pl.and.r2.lt.((alpha_F+kd*ms(i,ns,k))/alpha0)+pl) then
	ms(i,ns,k)=ms(i,ns,k)-1
	ms(i,1,k)=ms(i,1,k)+1
	endif
	alpha_F = alpha_F+kd*ms(i,ns,k)
	
	enddo
	enddo


	!------------ le ---------------
	do i=1,ns
	do k=1,ns

	do j=1,ns-1
	if (r2.ge.(alpha_F/alpha0)+pl.and.r2.lt.((alpha_F+kd*le(i,j,k))/alpha0)+pl) then
	le(i,j,k)=le(i,j,k)-1
	le(i,j+1,k)=le(i,j+1,k)+1
	endif
	alpha_F = alpha_F + kd*le(i,j,k)
	enddo

	if (r2.ge.(alpha_F/alpha0)+pl.and.r2.lt.((alpha_F+kd*le(i,ns,k))/alpha0)+pl) then
	le(i,ns,k)=le(i,ns,k)-1
	le(i,1,k)=le(i,1,k)+1
	endif
	alpha_F = alpha_F+kd*le(i,ns,k)
	
	enddo
	enddo


	endif						!------F-hop if-end





	pl=(beta_sum+alpha_R_sum+alpha_L_sum+alpha_B_sum+alpha_F_sum)/alpha0
	ph=(beta_sum+alpha_R_sum+alpha_L_sum+alpha_B_sum+alpha_F_sum+alpha_U_sum)/alpha0
	if (r2.ge.pl.and.r2.lt.ph) then 		!------U-hop if
	!-------- Up hopping---------

	alpha_U = 0.d0

	!----------- m --------------
	do i=1,ns
	do j=1,ns

	do k=1,ns-1
	if (r2.ge.(alpha_U/alpha0)+pl.and.r2.lt.((alpha_U+kd*m(i,j,k))/alpha0)+pl) then
	m(i,j,k)=m(i,j,k)-1
	m(i,j,k+1)=m(i,j,k+1)+1
	endif
	alpha_U = alpha_U + kd*m(i,j,k)
	enddo

	if (r2.ge.(alpha_U/alpha0)+pl.and.r2.lt.((alpha_U+kd*m(i,j,ns))/alpha0)+pl) then
	m(i,j,ns)=m(i,j,ns)-1
	m(i,j,1)=m(i,j,1)+1
	endif
	alpha_U = alpha_U+kd*m(i,j,ns)
	
	enddo
	enddo


	!----------- me --------------
	do i=1,ns
	do j=1,ns

	do k=1,ns-1
	if (r2.ge.(alpha_U/alpha0)+pl.and.r2.lt.((alpha_U+kd*me(i,j,k))/alpha0)+pl) then
	me(i,j,k)=me(i,j,k)-1
	me(i,j,k+1)=me(i,j,k+1)+1
	endif
	alpha_U = alpha_U + kd*me(i,j,k)
	enddo

	if (r2.ge.(alpha_U/alpha0)+pl.and.r2.lt.((alpha_U+kd*me(i,j,ns))/alpha0)+pl) then
	me(i,j,ns)=me(i,j,ns)-1
	me(i,j,1)=me(i,j,1)+1
	endif
	alpha_U = alpha_U+kd*me(i,j,ns)
	
	enddo
	enddo


	!----------- ms --------------
	do i=1,ns
	do j=1,ns

	do k=1,ns-1
	if (r2.ge.(alpha_U/alpha0)+pl.and.r2.lt.((alpha_U+kd*ms(i,j,k))/alpha0)+pl) then
	ms(i,j,k)=ms(i,j,k)-1
	ms(i,j,k+1)=ms(i,j,k+1)+1
	endif
	alpha_U = alpha_U + kd*ms(i,j,k)
	enddo

	if (r2.ge.(alpha_U/alpha0)+pl.and.r2.lt.((alpha_U+kd*ms(i,j,ns))/alpha0)+pl) then
	ms(i,j,ns)=ms(i,j,ns)-1
	ms(i,j,1)=ms(i,j,1)+1
	endif
	alpha_U = alpha_U+kd*ms(i,j,ns)
	
	enddo
	enddo


	!----------- le --------------
	do i=1,ns
	do j=1,ns

	do k=1,ns-1
	if (r2.ge.(alpha_U/alpha0)+pl.and.r2.lt.((alpha_U+kd*le(i,j,k))/alpha0)+pl) then
	le(i,j,k)=le(i,j,k)-1
	le(i,j,k+1)=le(i,j,k+1)+1
	endif
	alpha_U = alpha_U + kd*le(i,j,k)
	enddo

	if (r2.ge.(alpha_U/alpha0)+pl.and.r2.lt.((alpha_U+kd*le(i,j,ns))/alpha0)+pl) then
	le(i,j,ns)=le(i,j,ns)-1
	le(i,j,1)=le(i,j,1)+1
	endif
	alpha_U = alpha_U+kd*le(i,j,ns)
	
	enddo
	enddo



	endif						!------U-hop if-end





	pl=(beta_sum+alpha_R_sum+alpha_L_sum+alpha_B_sum+alpha_F_sum+alpha_U_sum)/alpha0
	ph=(beta_sum+alpha_R_sum+alpha_L_sum+alpha_B_sum+alpha_F_sum+alpha_U_sum+alpha_D_sum)/alpha0	
	if (r2.ge.pl.and.r2.lt.ph) then 		!------D-hop if
	!-------- down hopping---------

	alpha_D = 0.d0


	!------------ m --------------
	do i=1,ns
	do j=1,ns

	if (r2.ge.(alpha_D/alpha0)+pl.and.r2.lt.((alpha_D+kd*m(i,j,1))/alpha0)+pl) then
	m(i,j,1)=m(i,j,1)-1
	m(i,j,ns)=m(i,j,ns)+1
	endif
	alpha_D = alpha_D+kd*m(i,j,1)

	do k=2,ns
	if (r2.ge.(alpha_D/alpha0)+pl.and.r2.lt.((alpha_D+kd*m(i,j,k))/alpha0)+pl) then
	m(i,j,k)=m(i,j,k)-1
	m(i,j,k-1)=m(i,j,k-1)+1
	endif
	alpha_D = alpha_D + kd*m(i,j,k)
	enddo

	enddo
	enddo


	!------------ me --------------
	do i=1,ns
	do j=1,ns

	if (r2.ge.(alpha_D/alpha0)+pl.and.r2.lt.((alpha_D+kd*me(i,j,1))/alpha0)+pl) then
	me(i,j,1)=me(i,j,1)-1
	me(i,j,ns)=me(i,j,ns)+1
	endif
	alpha_D = alpha_D+kd*me(i,j,1)

	do k=2,ns
	if (r2.ge.(alpha_D/alpha0)+pl.and.r2.lt.((alpha_D+kd*me(i,j,k))/alpha0)+pl) then
	me(i,j,k)=me(i,j,k)-1
	me(i,j,k-1)=me(i,j,k-1)+1
	endif
	alpha_D = alpha_D + kd*me(i,j,k)
	enddo

	enddo
	enddo

	!------------ ms --------------
	do i=1,ns
	do j=1,ns

	if (r2.ge.(alpha_D/alpha0)+pl.and.r2.lt.((alpha_D+kd*ms(i,j,1))/alpha0)+pl) then
	ms(i,j,1)=ms(i,j,1)-1
	ms(i,j,ns)=ms(i,j,ns)+1
	endif
	alpha_D = alpha_D+kd*ms(i,j,1)

	do k=2,ns
	if (r2.ge.(alpha_D/alpha0)+pl.and.r2.lt.((alpha_D+kd*ms(i,j,k))/alpha0)+pl) then
	ms(i,j,k)=ms(i,j,k)-1
	ms(i,j,k-1)=ms(i,j,k-1)+1
	endif
	alpha_D = alpha_D + kd*ms(i,j,k)
	enddo

	enddo
	enddo


	!------------ le --------------
	do i=1,ns
	do j=1,ns

	if (r2.ge.(alpha_D/alpha0)+pl.and.r2.lt.((alpha_D+kd*le(i,j,1))/alpha0)+pl) then
	le(i,j,1)=le(i,j,1)-1
	le(i,j,ns)=le(i,j,ns)+1
	endif
	alpha_D = alpha_D+kd*le(i,j,1)

	do k=2,ns
	if (r2.ge.(alpha_D/alpha0)+pl.and.r2.lt.((alpha_D+kd*le(i,j,k))/alpha0)+pl) then
	le(i,j,k)=le(i,j,k)-1
	le(i,j,k-1)=le(i,j,k-1)+1
	endif
	alpha_D = alpha_D + kd*le(i,j,k)
	enddo

	enddo
	enddo



	endif						!------D-hop if-end




	
	if(t/(delt+nk*delt).ge.1.d0) then
	write(200,*)t,l1,l2,sum(m), sum(ms), sum(le), sum(me)

	do i=1,ns
	write(100,*) t, i*dx, m(ns/2,ns/2,i), me(ns/2,ns/2,i), ms(ns/2,ns/2,i), le(ns/2,ns/2,i)
	enddo
	write(100,*) ""

	nk=nk+1
	write(*,*) t, nk, (l2-l1)*2.d-4, 0.5d0*(l2+l1)*2.d-4
	endif


	t = t + tau

	enddo

!------------------ time loop ended -----------------------


	stop
	end program code







