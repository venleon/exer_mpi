PROGRAM wave2D

!*****************************************!
!* 2d shallow-Water model                *!
!*                                       *!
!* including:                            *!
!* - horizontal pressure-gradient force  *!
!* - Shapiro filter                      *!
!* - flooding algorithm                  *!
!*                                       *!
!* Author: J. Kaempf, 2008               *!

!* add mpi
!*****************************************!

!USE param
USE sub
USE precision_module
USE timing_module

!* use mpi
use params
use spaceparams
use xmpi_module

! local parameters
REAL :: time,dtmax
REAL :: c,lambda
INTEGER :: n


!* mpi parameters
type(parameters)         :: par
type(spacepars), pointer :: s
type(spacepars), target  :: sglobal
type(spacepars), target  :: slocal
character(len=80)        :: dummystring

!* mpi initialization
#ifdef USEMPI
	call xmpi_initialize
	! t0 = MPI_Wtime()
#endif
!s=>slocal

if (xmaster) then
#ifdef USEMPI
   write(*,*) 'MPI version, running on ',xmpi_size,'processes'
#endif	
	call space_alloc_scalars(sglobal)
	s => sglobal

	!**********
call space_alloc_scalars(sglobal)
	CALL init(s,par)  ! initialisation
	!**********

	! output of initial eta distribution
	OPEN(10,file ='eta0.dat',form='formatted')
	  DO j = 1,s%ny+1
	    WRITE(10,'(201F12.6)')(s%eta(j,k),k=1,s%nx+1)
	  END DO
	CLOSE(10)

	! output of initial layer thickness distribution
	OPEN(10,file ='h0.dat',form='formatted')
	  DO j = 1,s%ny+1
	    WRITE(10,'(201F12.6)')(s%hzero(j,k),k=1,s%nx+1)
	  END DO
	CLOSE(10)

	
	DO j = 1,s%ny
	DO k = 1,s%nx
	  par%hmax = MAX(par%hmax,s%h(j,k))
	END DO
	END DO

	! maximum phase speed
	c = SQRT(2*par%g*par%hmax)

	! determine stability parameter
	lambda = par%dt*SQRT(par%g*par%hmax)/MIN(s%dx,s%dy)

	IF(lambda > 1)THEN
	  WRITE(6,*) "This will not work. Do you know why?"   
	  STOP
	END IF

	! open files for output
	OPEN(10,file ='eta.dat',form='formatted')
!	OPEN(20,file ='h.dat',form='formatted')
!	OPEN(30,file ='u.dat',form='formatted')
!	OPEN(40,file ='v.dat',form='formatted')

	
	call start_timing()
endif


!s => sglobal
#ifdef USEMPI

call xmpi_determine_processor_grid(300,200)
if(xmaster) then
	write(*,*) 'nx',s%nx,'ny',s%ny   ! for debug
  write(*,*) 'processor grid: ',xmpi_m,' X ',xmpi_n
endif
#endif

write (*,*) xmaster
!#ifdef USEMPI
if(xmaster) then
	write(*,*) 'call distribute_par'
		call distribute_par(par)
endif




print 100, time_difference()
100 format (' time = ', f8.3)

END PROGRAM wave2D

subroutine printit(sglobal,slocal,par,it,s)
  use spaceparams
  use params
  use xmpi_module
  IMPLICIT none
  type(spacepars)          :: sglobal,slocal
  type(parameters)         :: par
  integer                  :: it
  character(len=*)         :: s
  integer,save             :: iter=0 
  return
  iter = iter+1
#ifdef USEMPI
  write(*,*) par%t,xmpi_rank,trim(s)
  call space_collect(slocal,sglobal%H,slocal%H)
  call space_collect(slocal,sglobal%zs,slocal%zs)
  call space_collect(slocal,sglobal%zs0,slocal%zs0)
  call space_collect(slocal,sglobal%u,slocal%u)
  call space_collect(slocal,sglobal%uu,slocal%uu)
  call space_collect(slocal,sglobal%ui,slocal%ui)
  call space_collect(slocal,sglobal%hh,slocal%hh)
  call space_collect(slocal,sglobal%vu,slocal%vu)
  call space_collect(slocal,sglobal%v,slocal%v)
  !call space_consistency(slocal,'ALL')
#else
  slocal%nx = slocal%nx  ! to prevent compiler warning about
                         ! unused slocal
#endif
  if(xmaster) call printsum(6,'H',1000*it+iter,sglobal%H)
  if(xmaster) call printsum(6,'zs',1000*it+iter,sglobal%zs)
  if(xmaster) call printsum(6,'zs0',1000*it+iter,sglobal%zs0)
  if(xmaster) call printsum(6,'u',1000*it+iter,sglobal%u)
  if(xmaster) call printsum(6,'uu',1000*it+iter,sglobal%uu)
  if(xmaster) call printsum(6,'ui',1000*it+iter,sglobal%ui)
  if(xmaster) call printsum(6,'hh',1000*it+iter,sglobal%hh)
  if(xmaster) call printsum(6,'vu',1000*it+iter,sglobal%vu)
  if(xmaster) call printsum(6,'v',1000*it+iter,sglobal%v)

  if(xmaster) print *,'par%t:',par%t
  if(xmaster) print *,'par%zs01:',par%zs01
#ifdef USEMPI
  if(xmaster) print *,'s%tideinpt:',slocal%tideinpt
  if(xmaster) print *,'s%tideinpz:',slocal%tideinpz(:,1)
#else
  if(xmaster) print *,'s%tideinpt:',sglobal%tideinpt
  if(xmaster) print *,'s%tideinpz:',sglobal%tideinpz(:,1)
#endif
end subroutine printit
