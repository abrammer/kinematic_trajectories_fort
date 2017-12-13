
	program main
	real p(4)
	real w 
	p(1) = 1
	p(2) = 3
	p(3) = 4
	p(4) = 7
	w= 0.5
!	print*, cubic_interpolate(p,w)
	print*, cubic_interp(p,w)
	end
