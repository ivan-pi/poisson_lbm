! 1D POISSON EQUATION SOLVER
! as presented in:
!
! "A novel lattice Boltzmann model for the Poisson equation"
! Zhenhua Chai, Baochang Shi, 2008, Applied Mathematical Modelling
!
! Ivan Pribec, Ljubljana, June 2015
 
! Solves the Poisson-Boltzman equation -> u';'; = k^2 * u
 
!#######################################################################
! PARAMETERS
!#######################################################################
module simParam
        ! precision (real kind)
    integer, parameter :: rKind = 4
 
        ! lattice weights
    real(kind = rKind), parameter ::  w(0:2)  = (/2.0/3.0, 1.0/6.0, 1.0/6.0/)
    real(kind = rKind), parameter :: ws(0:2)  = (/0.0, 0.5, 0.5/)
    real(kind = rKind), parameter ::   alpha  = 1.0/3.0
 
        ! lattice velocities
    integer, parameter :: e(0:2) = (/0,1,-1/)
 
        ! lattice dimension
    integer, parameter :: nx = 100
 
        ! maximum time value and data printing
    integer, parameter :: tMax  = 1000
    integer, parameter :: tPlot = 10
 
        ! lattice spacing and timestep
    real(kind = rKind), parameter :: dx = 0.01
    real(kind = rKind), parameter :: dt = 0.01
    real(kind = rKind), parameter :: c  = dx/dt
 
        ! relaxation time and diffusion coefficient
    real(kind = rKind), parameter :: tau   = 1.0
    real(kind = rKind), parameter :: omega = 1.0/tau
    real(kind = rKind), parameter :: D     = alpha*(c**2)*(0.5-tau)*dt
 
        ! steady state limit
    real(kind = rKind), parameter :: ssc = 0.0001
 
        ! poisson equation coefficient
    real(kind = rKind), parameter :: k = 27.79
 
end module simParam
 
!#######################################################################
! MAIN PROGRAM
!#######################################################################
program poisson1D
    use simParam
    implicit none
 
 
    real(kind = rKind) :: u(0:nx), u0(0:nx), diff(0:nx)
    real(kind = rKind) :: f(0:nx,0:2)
    real(kind = rKind) :: fEq
    integer :: x, i, t
 
        ! initialize u(x)
    u  = 0.0
    u0 = 1.0
    diff = abs(u-u0)
 
        ! boundary values u(0.0) = 1.0, u(1.0) = 1.0
    u(0)  = 1.0e0
    u(nx) = 1.0e0
 
        ! initialize f(x,i) to equilibrium
    do x = 0, nx
        do i = 0, 2
            f(x,i) = fEq(u(x),i)
        end do
    end do
 
        ! initialize time and counter
    t = 0
    call PrintData(f,u,t)
 
    ! ------------------------------------------------------------------
    ! MAIN (TIME) ITERATION LOOP - SOLVES POISSON EQUATION
    !-------------------------------------------------------------------
    do t = 1, tMax
 
            ! STREAMING
        f(1:nx,1)   = f(0:nx-1,1)   ! right e(1)
        f(0:nx-1,2) = f(1:nx,2)     ! left  e(2)
 
            ! BOUNDARY CONDITION
        u(0)  = 1.0e0
        u(nx) = 1.0e0
 
                ! Non-equilibrium extrapolation
            f(0,1)  = fEq(u(0),1)   + (1.0-omega)*(f(1,1)    - fEq(u(1),1))
            f(nx,2) = fEq(u(Nx),2)  + (1.0-omega)*(f(nx-1,2) - fEq(u(nx-1),2))
!~                 ! "Constant concentration"
!~             f(0,1)  = u(0) *(1.0-w(0)) - f(0,2)
!~             f(nx,2) = u(nx)*(1.0-w(0)) - f(nx,1)
 
            ! CALCULATE U (MACRO VARIABLE)
        do x = 0, nx
            u(x) = (f(x,1)+f(x,2))/(1.0-w(0))
        end do
 
        print*, sum(f)
 
            ! PRINT RESULTS
        if (mod(t,tPlot)==0) call PrintData(f,u,t)
 
            ! CALCULATE ABSOLUTE DIFFERENCE BETWEEN ITERATIONS
        diff = abs(u-u0)
        u0   = u
        if (all(diff < ssc)) then
            write(*,*) "Steady state condition was reached."
            exit
        end if
 
            ! COLLISION STEP
        do x = 0, nx
            do i = 0, 2
                f(x,i) = (1.0 - omega)*f(x,i) + omega*fEq(u(x),i) + dt*ws(i)*D*k*k*u(x)
            end do
        end do
 
    end do
 
    call PrintData(f,u,t)
 
    write(*,*) "Program finished."
 
end program poisson1D
 
 
 
!-----------------------------------------------------------------------
! equilibrium distribution function
!-----------------------------------------------------------------------
function fEq(u,i)
    use simParam, only: rKind, w
    implicit none
 
    real(kind = rKind) :: fEq, u
    integer :: i
 
    select case(i)
        case(0)
            fEq = (w(0)-1.0)*u
        case(1)
            fEq = w(1)*u
        case(2)
            fEq = w(2)*u
    end select
 
end function fEq
 
!-----------------------------------------------------------------------
! analytical solution for given Poisson equation
!-----------------------------------------------------------------------
function uSol(x)
    use simParam, only: rKind, k
    implicit none
 
    real(kind = rKind) :: uSol, x
 
    uSol = (exp(k) -1.0)/(exp(k)-exp(-k))*exp(-k*x) + &
         & (1.0-exp(-k))/(exp(k)-exp(-k))*exp(k*x)
 
end function uSol
 
!-----------------------------------------------------------------------
! prints results to file
!-----------------------------------------------------------------------
subroutine PrintData(f,u,t)
    use simParam, only: rKind, nx
    implicit none
 
    real(kind = rKind), intent(in) :: f(0:nx,0:2)
    real(kind = rKind), intent(in) :: u(0:nx)
    integer           , intent(in) :: t
 
    real(kind = rKind) :: uSol, ddx
    integer            :: x
    character(len=100) :: fileNum
 
    ddx = 1.0/real(nx)
 
        ! write iteration number to string
    write(fileNum,*) t
    fileNum = adjustl(fileNum)
 
        ! open file
    open(15,file=';iter_';//trim(fileNum)//';.dat';)
 
        ! write data to file for plotting
    write(15,*) "# x    uSol    uCalc  f(1)    f(2)    f(3)"
    do x = 0, nx
        write(15,*) x, uSol(real(x)*ddx), u(x), f(x,0), f(x,1), f(x,2)
    end do
 
        ! close file
    close(15)
 
    write(*,*) "Iteration = ", t
 
end subroutine PrintData