! program for nonparaxial beam propagation in free space

program NBPE_free_space
	implicit none
	
! 	parameters
	complex*16, parameter		::	xi		=	(0.d0,1.d0)
	real*8, parameter			::	pi		=	4.d0*datan(1.d0)

!	laser parameters
	real*8, parameter			::	wavelength	=	1.3e-4				! 1.3 micrometer in cm
	real*8, parameter			::	wavevector	=	2.d0*pi/wavelength


!	transverse space dimenisions
	real*8, parameter			::	X_min	=	-2.d0, 	X_max	=	2.d0
	real*8, parameter			::	X_range	=	abs(X_min - X_max)		! this is the normalized radial distance in units of the normalization constant, k_0*n_core*d_core.
	real*8, parameter			::	dx		=	1e-2/4               	! step size along the transverse direction
	integer, parameter			::	lx		=	X_range/dx            	! no. of points along the transverse direction
	integer, parameter			::	lx2		=	lx/2


!	refractive_index
	real*8, dimension(lx)		::	refractive_index
	real*8, parameter			::	n_cladding	= 1.465, n_core = 1.491		! refractive indices of core and cladding
	real*8, parameter			::	d_core		= 62.5e-4					! diameter of core in cm (62.5 micrometer)
	real*8, parameter			::	r_core 		= d_core/2					! radius of core in cm (62.5 micrometer)

!	extra space normalization constant
	real*8, parameter			::	normalization_constant = wavevector*n_core*d_core


!	longitudinal space dimensions
	real*8, parameter			::	alpha_1	= sqrt((4.d0/normalization_constant**2)*(1.d0 - (n_cladding/n_core)**2))
!	real*8, parameter			::	self_imaging_length = (2.d0*pi/alpha_1)
	real*8, parameter			::	self_imaging_length = (2.d0*pi/alpha_1)*sqrt(1.d0 - alpha_1)
	real*8, parameter			::	dz		=	1e-4
	real*8, parameter			::	Z_max	=	2.d0
	real*8, parameter			::	fiber_length	=	1.d0
	integer, parameter			::	lz		=	int(z_max/dz)
	real*8, parameter			::	dz2		=	dz/2.d0


!	elements of the tridiagonal matrix
	real*8, dimension(lx)			::	A, B
	complex*16, dimension(2, lx)	::	zeta, lambda, lambda_tilde, phi
	complex*16						::	alpha, beta
	complex*16, dimension(lx)		::	rho, rhs
	complex*16, dimension(lx-1)		::	gamma



!	normalized fundamental mode width of the waveguide 
	real*8, parameter				::	w_f = 1.d0/dsqrt(alpha_1)
	
!	normalized Gaussian beam width
	real*8, parameter				::	w_0	= 0.8d0*w_f

!	center of the normalized gaussain beam, also known as missalignment
	real*8, parameter				::	X_0 = 0.2d0*normalization_constant

!	normalized coordinate variables	
	real*8						::	Z, X

!	dummy integers for do loop	
	integer						::	i, j, k, l

!	matrices for the fields	
	complex*16, dimension(lx)	::	u_0, u_n, u_n_plus_one, dudz

!	field normalization constants
	real*8						::	u_0_max, u_n_plus_one_max
	
!	field power calculation
	real*8						::	input_power, output_power, temp

!	filenames for saving
	character*200				::	filename



!      Tranverse profile at z = 0
!#####################################
	do i = 1, lx
		X = (i - lx2)*dx*normalization_constant
		u_0(i) = (1.d0/(w_0*sqrt(2.d0*pi)))*dexp(-((X - X_0)**2/(w_0**2)))
!		u_0(i) = (1.d0/(w_0*sqrt(2.d0*pi)))*dexp(-((X - X_0)**2/(2.d0*w_0**2)))
	enddo
	u_0_max  = maxval(abs(u_0))
	
!      Tranverse refractive_index profile
!#####################################
	do i = 1, lx
		X = ((i - lx2)*dx/(wavevector*n_core))*normalization_constant
		if (X.lt. -r_core .or. X.gt. r_core ) then						! cladding's refractive index
			refractive_index(i) = n_cladding
		elseif(X .ge. -r_core .or. X .le. r_core) then					! core's refractive index
			refractive_index(i) = n_core*dsqrt(1.d0 - 4.d0*(1.d0 - (n_cladding/n_core)**2)*(X/d_core)**2)
!		elseif(X .le. d_core) then					! outside cladding
!			refractive_index(i) = 1.d0
		endif
	enddo


!		Plot the transverse refractive index profile
!############################################################

	write(filename,'("transverse_refractive_index_profile.dat")')		
	open( unit = 5000, file = filename, STATUS='UNKNOWN')

	do i = 1, lx
		X = (i - lx2)*dx
		write(5000,*) X, refractive_index(i)
	enddo
	close(5000)
	
! 		Plot the input transverse field profile
!############################################################
	write(filename,'("input_transverse_profile.dat")')		
	open( unit = 5000, file = filename, STATUS='UNKNOWN')

	do i = 1, lx
		X = (i - lx2)*dx
		write(5000,*) X, abs(u_0(i))/u_0_max
	enddo
	close(5000)



! Repeat the Thomas algorithm for all points along z axis.

!    intialize fields at z = 0
!###############################################

	u_n = u_0

	write(filename,'("peak_amplitude_vs_z_non_paraxial.dat")')		
	open( unit = 5100, file = filename, STATUS='UNKNOWN')
	
	write(filename,'("longitudinal_profille.dat")')		
	open( unit = 5200, file = filename, STATUS='UNKNOWN')
	

	write(*,*) "Max propagation length = ", z_max


!    Propagation begins here
!###########################################################
	write(*,*)"No. of points along z axis: ", lz


!	Parameters of the tridiagonal matrix element
!##########################################################

	do i = 1, lx
		A(i)				=	(n_core/(2.d0*refractive_index(i)))/self_imaging_length
		B(i)				=	(A(i)*self_imaging_length**2)/normalization_constant**2
		zeta(1, i)			=	B(i)*(3.d0*A(i) - xi*dz2/2.d0)
		zeta(2, i)			=	A(i)*B(i)**2*(A(i) - xi*dz2)
		
		phi(2, i)			=	(zeta(1, i) + sqrt(zeta(1,i)**2 - 4.d0*zeta(2,i)))/2.d0
		phi(1, i)			=	zeta(2, i)/phi(2, i)

		lambda(1, i) 		=	phi(1, i)/dx**2
		lambda(2, i) 		=	phi(2, i)/dx**2
		lambda_tilde(1, i)	=	1.d0 - 2.d0*lambda(1, i)
		lambda_tilde(2, i)	=	1.d0 - 2.d0*lambda(2, i)
	enddo


	l = 0
	do k = 1, lz
		Z = (k-1)*dz
!	Parameters of the tridiagonal matrix element
!##########################################################

		if(z .gt. fiber_length) then
			do i = 1, lx
				refractive_index(i) = 1.d0
				A(i)				=	(n_core/(2.d0*refractive_index(i)))/self_imaging_length
				B(i)				=	(A(i)*self_imaging_length**2)/normalization_constant**2
				zeta(1, i)			=	B(i)*(3.d0*A(i) - xi*dz2/2.d0)
				zeta(2, i)			=	A(i)*B(i)**2*(A(i) - xi*dz2)
		
				phi(2, i)			=	(zeta(1, i) + sqrt(zeta(1,i)**2 - 4.d0*zeta(2,i)))/2.d0
				phi(1, i)			=	zeta(2, i)/phi(2, i)

				lambda(1, i) 		=	phi(1, i)/dx**2
				lambda(2, i) 		=	phi(2, i)/dx**2
				lambda_tilde(1, i)	=	1.d0 - 2.d0*lambda(1, i)
				lambda_tilde(2, i)	=	1.d0 - 2.d0*lambda(2, i)
			enddo
		endif

! step 1: Propagation through medium
!###########################################################
		call derive(u_n,dudz, refractive_index)
		call rk4(u_n, dudz, u_n_plus_one, dz, refractive_index)
		u_n = u_n_plus_one

! step 2: Propagation through free space
! #########################################################

		do j = 1, 2

		!      Define the rhs column matrix
		! #########################################################
			rhs(1)	=	dconjg(lambda_tilde(j,1))*u_n(1) + dconjg(lambda(j,1))*u_n(2)
			do i = 2, lx - 1
				alpha	=	lambda(j, i)
				beta	=	lambda(j, i)
				rhs(i)	=	dconjg(alpha)*u_n(i-1) + dconjg(lambda_tilde(j, i))*u_n(i) + dconjg(beta)*u_n(i+1)
			enddo
			alpha		=	lambda(j, lx)
			beta		=	lambda(j, lx)
			rhs(lx)	=	dconjg(alpha)*u_n(lx-1) + dconjg(lambda_tilde(j, lx))*u_n(lx)

		!       Define the gamma_i
		! ########################################################
			gamma(1)	=	lambda(j, 1)/lambda_tilde(j,1)
			do i = 2, lx-1
				alpha		=	lambda(j, i)
				beta		=	lambda(j, i)
				gamma(i)	=	beta/(lambda_tilde(j, i) - alpha*gamma(i-1))
			enddo

!       Define the rho column matrix
! ########################################################
			rho(1) = rhs(1)/lambda_tilde(j,1)
			do i = 2,lx
				alpha	=	lambda(j, i)
				rho(i) 	=	(rhs(i) - alpha*rho(i-1))/(lambda_tilde(j, i) - alpha*gamma(i-1))
			enddo

!       solution at n+1
! ########################################################
			u_n_plus_one(lx) = rho(lx)
			do i = lx-1, 1,-1
				u_n_plus_one(i) = rho(i) - gamma(i)*u_n_plus_one(i+1)
			enddo
			u_n = u_n_plus_one
		enddo
		
!   	plot the peak intensity vs longitudinal distance
!##########################################################
		write(5100,*)z, maxval(abs(u_n_plus_one)/u_0_max)
		
		if (k .eq. 1 + int(lz/1000)*l) then
			write(5200,"(*(g0))") (abs(u_n_plus_one(j))/u_0_max," ",j = 1,lx,1)
			l = l + 1
		endif

	enddo
	close(5100)
	close(5200)


! 				Plot the results
!############################################################
!	Tranverse profile
!############################################################
	write(filename,'("output_transverse_profile_non_paraxial.dat")')		
	open( unit = 5200, file = filename, STATUS='UNKNOWN')

	u_n_plus_one_max = maxval(abs(u_n_plus_one))
	do i = 1, lx
		X = (i - lx2)*dx
		write(5200,*)	X, abs(u_n_plus_one(i))/u_0_max, abs(u_n_plus_one(i))/u_n_plus_one_max
	enddo
	close(5200)

!############################################################
!	Subroutines
!############################################################

contains

	subroutine derive(y, dydz, refractive_index)
		implicit none
		complex*16, dimension(:), intent(in)	::	y
		real*8, dimension(:), intent(in)		::	refractive_index
		complex*16, dimension(:), intent(inout)	::	dydz
		real*8, dimension(lx)					::	susceptibility
		
!		susceptibility	=	(refractive_index**2 - 1.d0)/(4.d0*pi)		! this gives wrong result
		susceptibility	=	(refractive_index - 1.d0)/(2.d0*pi)

		dydz = xi*(2.d0*pi*(refractive_index/n_core)*self_imaging_length)*susceptibility*y

		return
		return
	end subroutine derive

!      subroutine for runge kutta
!****************************************************************************************************** !
	subroutine rk4(y, dydz, y_out, dz, refractive_index)
		implicit none
		complex*16, dimension(:), intent(in)		::	y, dydz
		real*8, dimension(:), intent(in)			::	refractive_index
		complex*16, dimension(:), intent(inout)		::	y_out
	
		real*8,	intent(in)			::	dz
		
		complex*16, dimension(lx)	::	yt, dyt, dym
	
		real*8						::	h, hh, h6
	
		h  = dz
		hh = h*0.5d0
		h6 = h/6.d0

		yt  = y + hh*dydz      						! 1st step.
		call derive(yt,dyt, refractive_index)    	! 2nd step.
		yt  = y + hh*dyt

		call derive(yt,dym, refractive_index)    	! 3rd step.
		yt  = y + h*dym

		dym = dyt + dym
		call derive(yt,dyt, refractive_index)     	! 4th step.

		y_out = y + h6*(dydz + dyt + 2.0*dym)

		return
	end subroutine rk4


end program NBPE_free_space
