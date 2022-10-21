! program for nonparaxial beam propagation in free space

program NBPE_free_space
	implicit none
	
! 	parameters

	complex*16, parameter		::	xi		=	(0.d0,1.d0)
	real*8, parameter			::	pi		=	4.d0*datan(1.d0)	

! 	beam parameters
	real*8						::	kappa

!	space dimenisions
	integer, parameter			::	lr		=	16384               	! no. of points along the transverse direction
	integer, parameter			::	lr2		=	lr/2
	real*8, parameter			::	dr		=	1e-2               	! step size along the transverse direction
	real*8, parameter			::	dz		=	1e-2
	real*8, parameter			::	z_max	=	1.d0
	integer, parameter			::	lz		=	int(z_max/dz)

!	elements of the tridiagonal matrix
	real*8						::	A
	real*8, parameter			::	B	=	1.d0/4.d0
	complex*16					::	zeta, lambda, lambda_tilde, alpha, beta
	complex*16, dimension(lr2)	::	gamma, rho, rhs

!	coordinate variables	
	real*8						::	z, r

!	dummy integers for do loop	
	integer						::	i,j, k

!	matrices for the fields	
	complex*16, dimension(lr2)	::	u_0, u_n, u_n_plus_one, u_z_R_theo

!	field normalization constants
	real*8						::	u_0_max, u_n_plus_one_max, u_z_R_theo_max
	
!	area under the fields
	real*8						::	area_f, area_i, temp
	
!	filenames for saving
	character*200				::	filename

!      Tranverse profile at z = 0
!#####################################
	do i = 1, lr2
		r = (i - 1)*dr
		u_0(i) = dexp(-r**2)
	enddo
	u_0_max  = maxval(abs(u_0))

!      Tranverse profile at z = one Rayleigh length
!#########################################################
	z = z_max
	do i = 1, lr2
		r = (i - 1)*dr
		u_z_R_theo(i) = sqrt(1.d0/(1.d0 + z))*dexp(-r**2/(1.d0 + z))
	enddo
	u_z_R_theo_max  = maxval(abs(u_z_R_theo))
    
!		Plot the transverse profile of the field at z = 0 and the theoretical field at z = z_R
!############################################################################
	do i = lr2, 1, -1
		r = -(i - 1)*dr
		write(1,*) r, abs(u_0(i))/u_0_max, abs(u_z_R_theo(i))/u_0_max, abs(u_z_R_theo(i))/u_z_R_theo_max
	enddo
	do i = 1, lr2
		r = (i - 1)*dr
		write(1,*) r, abs(u_0(i))/u_0_max, abs(u_z_R_theo(i))/u_0_max, abs(u_z_R_theo(i))/u_z_R_theo_max
	enddo

! Repeat the Thomas algorithm for all points along z axis.

!    intialize fields at z = 0
!###############################################

	u_n = u_0

!	Parameters of the tridiagonal matrix element
!##########################################################
	kappa			=	1e-3		! this represents paraxial regime
!	kappa			=	0.1d0		! this represents nonparaxial regime
!	kappa			=	0.5d0		! this represents nonparaxial regime
	A				=	kappa
	zeta			=	B*(A - xi*dz/2.d0)
	lambda 			=	zeta/dr**2
	lambda_tilde	=	1.d0 - 2.d0*lambda

	write(filename,'("z_profile.dat")')		
	open( unit = 5000, file = filename, STATUS='UNKNOWN')

!    Propagation begins here
!###########################################################

	do k = 1, lz
!      Define the rhs column matrix
! #########################################################
		rhs(1)	=	dconjg(lambda_tilde)*u_n(1) - 3.d0*dconjg(lambda)*u_n(2)
		do i = 2, lr2-1
			j = i-1
			alpha	=	lambda*(1.d0 - 0.5d0/j)
			beta	=	lambda*(1.d0 + 0.5d0/j)
			rhs(i)	=	dconjg(alpha)*u_n(i-1) + dconjg(lambda_tilde)*u_n(i) + dconjg(beta)*u_n(i+1)
		enddo
		alpha		=	lambda*(1.d0 - 0.5d0/(lr2-1))
		beta		=	lambda*(1.d0 + 0.5d0/(lr2-1))
		rhs(lr2)	=	dconjg(alpha)*u_n(lr2-1) + dconjg(lambda_tilde)*u_n(lr2)

!       Define the gamma_i
! ########################################################
		gamma(1)	=	2.d0*lambda/lambda_tilde
		do i = 2, lr2-1
			j = i-1
			alpha		=	lambda*(1.d0 - 0.5d0/j)
			beta		=	lambda*(1.d0 + 0.5d0/j)
			gamma(i)	=	beta/(lambda_tilde - alpha*gamma(i-1))
		enddo
		alpha			=	lambda*(1.d0 - 0.5d0/lr2)
		beta			=	lambda*(1.d0 + 0.5d0/lr2)
		gamma(lr2)		=	beta/((lambda_tilde + beta)  - alpha*gamma(i-1))

!       Define the rho column matrix
! ########################################################
		rho(1) = rhs(1)/lambda_tilde
		do i = 2,lr2
			j = i -1
			alpha	=	lambda*(1.d0 - 0.5d0/j)
			rho(i) 	=	(rhs(i) - alpha*rho(i-1))/(lambda_tilde - alpha*gamma(i-1))
		enddo

!       solution at n+1
! ########################################################
		u_n_plus_one(lr2) = rho(lr2)
		do i = lr2-1, 1,-1
			u_n_plus_one(i) = rho(i) - gamma(i)*u_n_plus_one(i+1)
		enddo

		u_n = u_n_plus_one

!   	plot the peak intensity vs longitudinal distance
!##########################################################

		z = (k-1)*dz
		write(5000,*)	z, maxval(abs(u_n_plus_one)/u_0_max)
	enddo

! 				Plot the results
!############################################################
!	Tranverse profile
!############################################################
	u_n_plus_one_max = maxval(abs(u_n_plus_one))
	do i = lr2, 1, -1
		r = -(i - 1)*dr
		write(2,*)	r, abs(u_n_plus_one(i))/u_0_max, abs(u_n_plus_one(i))/u_n_plus_one_max
	enddo
	do i = 1, lr2
		r = (i - 1)*dr
		write(2,*)	r, abs(u_n_plus_one(i))/u_0_max, abs(u_n_plus_one(i))/u_n_plus_one_max
	enddo

! 	Calculate the initial area under the curve
!############################################################
	area_i = 0.d0
	do i = 1, lr2
		r = (i - 1)*dr
		temp 	= r*abs(u_0(i))**2*dr
		area_i 	= area_i + temp
	enddo
	write(*,*) "intial area = ",area_i

! 	Calculate the final area under the curve
!############################################################
	area_f = 0.d0
	do i = 1, lr2
		r = (i - 1)*dr
		temp 	= r*abs(u_n_plus_one(i))**2*dr
		area_f 	= area_f + temp
	enddo
	write(*,*) "final area = ", area_f

	write(*,*) "difference = ", (area_f - area_i)

end program NBPE_free_space