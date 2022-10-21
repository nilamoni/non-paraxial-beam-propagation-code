! program for nonparaxial beam propagation in free space

program NBPE_free_space
	implicit none
	
! 	parameters
	complex*16, parameter		::	xi				=	(0.d0,1.d0)
	real*8, parameter			::	pi				=	4.d0*datan(1.d0)

!	Some constants in CGS
	real*8, parameter			::	speed_of_light	=	3e10

! 	beam parameters
	real*8, parameter			::	ratio			=	10.d0					! this the ratio between wavelength and w_0, i.e., w_0 = 5*wavelength
	real*8, parameter			::	wavelength		=	1.d0/ratio
	real*8, parameter			::	mod_amplitude	=	1e2
	real*8						::	w_z
	real*8						::	kappa

!	medium parameter
	real*8, parameter			::	n_1			=	1.51d0
	real*8, parameter			::	n_2 		= 	6.4e-16
	real*8, parameter			::	n_2_bar		=	n_2*n_1*speed_of_light
	real*8, parameter			::	N_square 	= 	(2.d0*pi*ratio)**2*n_1*n_2_bar*mod_amplitude
	real*8, parameter			::	Z_collapse	=	1.d0/sqrt(N_square - 1.d0)


!	space dimenisions
	integer, parameter			::	lr		=	16384*4               		! no. of points along the transverse direction
	integer, parameter			::	lr2		=	lr/2
	real*8, parameter			::	dr		=	1e-4               			! step size along the transverse direction
	real*8, parameter			::	dz		=	1e-4
	real*8, parameter			::	z_max	=	3.d0*Z_collapse
	integer, parameter			::	lz		=	int(z_max/dz)

!	elements of the tridiagonal matrix
	real*8						::	A
	real*8, parameter			::	B	=	1.d0/4.d0
	complex*16					::	zeta, lambda, lambda_tilde, alpha, beta
	complex*16, dimension(lr2)	::	gamma, rho, rhs

!	coordinate variables	
	real*8						::	z, r

!	dummy integers for do loop	
	integer						::	i, k

!	matrices for the fields	
	complex*16, dimension(lr2)	::	u_0, u_n, u_n_plus_one, dudz

!	field normalization constants
	real*8						::	u_0_max, u_n_plus_one_max
	
!	field power calculation
	real*8						::	input_power, output_power, temp

!	filenames for saving
	character*200				::	filename

!      Tranverse profile at z = 0
!#####################################
	do i = 1, lr2
		r = (i - 1)*dr
		u_0(i) = dexp(-r**2)
	enddo
	u_0_max  = maxval(abs(u_0))

!		Plot the transverse profile of the field at z = 0 and the theoretical field at z = z_R
!############################################################################
	do i = lr2, 1, -1
		r = -(i - 1)*dr
		write(1,*) r, abs(u_0(i))/u_0_max
	enddo
	do i = 1, lr2
		r = (i - 1)*dr
		write(1,*) r, abs(u_0(i))/u_0_max
	enddo

! Repeat the Thomas algorithm for all points along z axis.

!    intialize fields at z = 0
!###############################################

	u_n = u_0

!	Parameters of the tridiagonal matrix element
!##########################################################
	w_z				=	1.d0
	kappa			=	(wavelength/(2.d0*pi*w_z*n_1))**2
	A				=	kappa

	write(filename,'("peak_amplitude_vs_z_adaptive.dat")')		
	open( unit = 5000, file = filename, STATUS='UNKNOWN')

	write(*,*) "ratio = ",		ratio
	write(*,*) "n_1 = ",		n_1
	write(*,*) "n_2 = ",		n_2
	write(*,*) "n_2_bar = ",	n_2_bar
	write(*,*) "kappa = ", 		kappa
	write(*,*) "N_square = ", 	N_square
	write(*,*) "Collapsing length = ", 	Z_collapse
	write(*,*) "Max propagation length = ", z_max

!    Propagation begins here
!###########################################################

	do k = 1, lz
		A 				=	kappa
		zeta			=	B*(A - xi*dz/2.d0)
		lambda 			=	zeta/dr**2
		lambda_tilde	=	1.d0 - 2.d0*lambda
	
! step 1: Propagation through medium
!###########################################################
		call derive(u_n,dudz)
		call rk4(u_n, dudz, u_n_plus_one, dz)
		u_n = u_n_plus_one

! step 2: Propagation through free space
! #########################################################

!      Define the rhs column matrix
! #########################################################
		rhs(1)	=	dconjg(lambda_tilde)*u_n(1) + 2.d0*dconjg(lambda)*u_n(2)
		do i = 2, lr2-1
			alpha	=	lambda*(1.d0 - 0.5d0/(i-1))
			beta	=	lambda*(1.d0 + 0.5d0/(i-1))
			rhs(i)	=	dconjg(alpha)*u_n(i-1) + dconjg(lambda_tilde)*u_n(i) + dconjg(beta)*u_n(i+1)
		enddo
		alpha		=	lambda*(1.d0 - 0.5d0/(lr2-1))
		beta		=	lambda*(1.d0 + 0.5d0/(lr2-1))
		rhs(lr2)	=	dconjg(alpha)*u_n(lr2-1) + dconjg(lambda_tilde)*u_n(lr2)

!       Define the gamma_i
! ########################################################
		gamma(1)	=	2.d0*lambda/lambda_tilde
		do i = 2, lr2
			alpha		=	lambda*(1.d0 - 0.5d0/(i-1))
			beta		=	lambda*(1.d0 + 0.5d0/(i-1))
			gamma(i)	=	beta/(lambda_tilde - alpha*gamma(i-1))
		enddo

!       Define the rho column matrix
! ########################################################
		rho(1) = rhs(1)/lambda_tilde
		do i = 2,lr2
			alpha	=	lambda*(1.d0 - 0.5d0/(i-1))
			rho(i) 	=	(rhs(i) - alpha*rho(i-1))/(lambda_tilde - alpha*gamma(i-1))
		enddo

!       solution at n+1
! ########################################################
		u_n_plus_one(lr2) = rho(lr2)
		do i = lr2-1, 1,-1
			u_n_plus_one(i) = rho(i) - gamma(i)*u_n_plus_one(i+1)
		enddo

		u_n = u_n_plus_one
		
		
!   	For finding beam radius at any given z
!##########################################################
		u_n_plus_one_max = maxval(abs(u_n_plus_one))
		do i = 1,lr2, 1
			if (abs(u_n_plus_one(i)) .lt. u_n_plus_one_max/exp(1.d0)) then
				w_z = abs((i - 1)*dr)
				kappa = (wavelength/(2.d0*pi*w_z*n_1))**2
				goto 10
			endif
		enddo
10		continue

!   	plot the peak intensity vs longitudinal distance
!##########################################################

		z = (k-1)*dz
		write(5000,*)z, maxval(abs(u_n_plus_one)/u_0_max), w_z
	enddo

! 	Calculate the input_power
!############################################################
	input_power = 0.d0
	do i = 1, lr2
		r = (i - 1)*dr
		temp 	= r*abs(u_0(i))**2*dr
		input_power 	= input_power + temp
	enddo
	write(*,*) "input power = ", input_power

! 	Calculate the output_power
!############################################################
	output_power = 0.d0
	do i = 1, lr2
		r = (i - 1)*dr
		temp 			= r*abs(u_n_plus_one(i))**2*dr
		output_power 	= output_power + temp
	enddo
	write(*,*) "Output power = ", output_power

	write(*,*) "difference = ", output_power - input_power

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

contains

	subroutine derive(y, dydz)
		implicit none
		complex*16, dimension(:), intent(in)	::	y
		complex*16, dimension(:), intent(inout)	::	dydz

		dydz = xi*N_square*y*abs(y)**2
		return
		return
	end subroutine derive

!      subroutine for runge kutta
!****************************************************************************************************** !
	subroutine rk4(y, dydz, y_out, dz)
		implicit none
		complex*16, dimension(:), intent(in)		::	y, dydz
		complex*16, dimension(:), intent(inout)		::	y_out
	
		real*8,	intent(in)			::	dz
		
		complex*16, dimension(lr2)	::	yt, dyt, dym
	
		real*8						::	h, hh, h6
	
		h  = dz
		hh = h*0.5d0
		h6 = h/6.d0

		yt  = y + hh*dydz        ! 1st step.
		call derive(yt,dyt)    	! 2nd step.
		yt  = y + hh*dyt

		call derive(yt,dym)    	! 3rd step.
		yt  = y + h*dym

		dym = dyt + dym
		call derive(yt,dyt)     ! 4th step.

		y_out = y + h6*(dydz + dyt + 2.0*dym)

		return
	end subroutine rk4


end program NBPE_free_space