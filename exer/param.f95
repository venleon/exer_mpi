MODULE param

INTEGER(4), PARAMETER :: nx = 201	
INTEGER(4), PARAMETER :: ny = 201	

REAL :: hzero(0:ny+1,0:nx+1), h(0:ny+1,0:nx+1)
REAL :: eta(0:ny+1,0:nx+1),etan(0:ny+1,0:nx+1)
REAL :: u(0:ny+1,0:nx+1), un(0:ny+1,0:nx+1)
REAL :: v(0:ny+1,0:nx+1), vn(0:ny+1,0:nx+1)
REAL :: dt,dx,dy,g
REAL,PARAMETER :: eps = 0.005 ! parameter for Shapiro filter

INTEGER :: i,j,k

INTEGER :: wet(0:ny+1,0:nx+1)
REAL :: hmin

END MODULE param
