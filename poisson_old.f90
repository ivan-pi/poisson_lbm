program poisson1D
    implicit none
    integer :: x, t, tMax, old, new, temp
    integer, parameter :: nx = 100
    real :: w0(0:2), w1(0:2), k, tau, omega, alpha, D, dx, c, uA, dt
    real :: u(0:nx), f(0:nx,0:2,2), feq(0:nx,0:2)
    old = 1
    new = 2

    ! weighting factors
    w0(0:2) = (/2./3.,1./6.,1./6./)
    w1(0:2) = (/0.0,0.5,0.5/)

    ! number of nodes and timesteps
    dx = 1.0/real(nx)
    c = 1.0
    dt = dx/c
    tMax = 1000000

    ! coefficient
    k = 27.79

    ! relaxation time
    tau = 1.0
    omega = 1.0/tau

    ! artificial diffusion coefficient
    alpha = 1./3.
    D = alpha*c**2*(0.5-tau)*dt

    !----------------------------------------------------

    ! initial condition
    u(0:nx) = 0.0
    u(0)    = 1.0
    u(nx)   = 1.0

    ! initialize to equilibrium
    f(0:nx,0,old) = (w0(0)-1.0)*u(0:nx)
    f(0:nx,1,old) = w0(1)*u(0:nx)
    f(0:nx,2,old) = w0(2)*u(0:nx)

    do t = 1, tMax

        ! macroscopic variable
        u(0:nx) = (f(0:nx,1,old)+f(0:nx,2,old))/(1.0-w0(0))

        ! equilibrium-
        feq(0:nx,0) = (w0(0)-1.0)*u(0:nx)
        feq(0:nx,1) = w0(1)*u(0:nx)
        feq(0:nx,2) = w0(2)*u(0:nx)

        ! streaming + collision
        f(0:nx,0,new)   = (1-omega)*f(0:nx,0,old) + omega*feq(0:nx,0)
        f(1:nx,1,new)   = (1-omega)*f(0:nx-1,1,old) + omega*feq(0:nx-1,1) + dt*D*w1(1)*k**2*u(0:nx-1)
        f(0:nx-1,2,new) = (1-omega)*f(1:nx,2,old) + omega*feq(1:nx,2) + dt*D*w1(2)*k**2*u(1:nx)

        ! boundaries
        f(0,1,new)  = w0(1)*1.0 + (1-omega)*(f(1,1,new)  - feq(1,1))
        f(nx,2,new) = w0(2)*1.0 + (1-omega)*(f(nx-1,2,new) - feq(nx-1,2))

        temp = old
        old = new
        new = temp

    end do

    open(10,file="results.dat")

    do x = 0, nx
        write(10,*) x*dx, u(x), uA(real(x)*dx,k)
    end do

    close(10)

end program

function uA(x,k)
    implicit none
    real:: uA, x, k

    uA = (exp(k)-1.0)/(exp(k)-exp(-k))*exp(-k*x) + &
       & (1.0-exp(-k))/(exp(k)-exp(-k))*exp(k*x)

end function uA
