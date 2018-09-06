module poisson_example_mod

    use precision_mod, only: wp
    use source_mod, only: source_term
    use poisson_mod, only: Lattice, init_Lattice

    implicit none
    private

    public :: pb_src, init_pb_src, poisson_example, analytical_solution

    type, extends(source_term) :: pb_src ! poisson boltzmann source
        real(wp) :: k
    contains
        procedure :: eval
    end type

contains

    function init_pb_src(k) result(this)
        real(wp), intent(in) :: k
        type(pb_src) :: this
        this%k = k
    end function

    real(wp) function eval(this,u,x)
        class(pb_src), intent(in) :: this
        real(wp), intent(in) :: u
        real(wp), intent(in) :: x
        eval = this%k**2*u
    end function

!>  Returns analytical solution of $\Delta u = k^2 u$ on
!   the domain $\[0,1\]$ with constant Dirichlet boundary conditions
!   $u(0) = 1$, $u(1) = 1$.
!
!   See Eq. (3.4) in paper by Chai and Shi (2008)
!
    elemental function analytical_solution(x,k) result(u)
        real(wp), intent(in) :: x
        real(wp), intent(in) :: k
        real(wp) :: u

        u = (exp(k) - 1.0_wp)/(exp(k) - exp(-k))*exp(-k*x) + &
           & (1.0_wp - exp(-k))/(exp(k) - exp(-k))*exp(k*x)
    end function


    subroutine poisson_example(k,nx,steps,check,tol)
        real(wp), intent(in) :: k
        integer, intent(in) :: nx, steps, check
        real(wp) :: tol

        real(wp), allocatable :: uold(:), u(:)
        real(wp) :: gre, omega
        type(Lattice) :: latt
        type(pb_src), target :: mysrc
        integer :: i, step

        ! Initalize density arrays
        allocate(uold(nx),u(nx))
        uold = 0.0_wp; uold([1,nx]) = 1._wp

        ! Create source term
        mysrc = init_pb_src(k)
        
        ! Create lattice
        omega = 1.0_wp
        latt = init_Lattice(nx=nx,uinit=uold,omega=omega,src=mysrc)

        time_loop: do step = 1, steps

            call latt%collide_and_stream()

            call latt%set_left_neem(1.0_wp)
            call latt%set_right_neem(1.0_wp)
            
            ! Check convergence
            if (mod(step,check) == 0) then
                call latt%get_macroscopic_values(u)
                gre = sum(abs(u-uold))/sum(abs(u)) ! global relative error
                uold = u
                write(*,*) "step = ", step, "error = ", gre
                if (gre < tol) exit time_loop
            end if

            call latt%prep_next_step()
        end do time_loop

        call latt%output("poisson_results.txt")

    end subroutine

end module


module reaction_example_mod

    use precision_mod, only: wp
    use source_mod, only: source_term
    use poisson_mod, only: Lattice, init_Lattice

    implicit none
    private

    public :: reaction_example

    type, extends(source_term) :: rc_src ! poisson boltzmann source
        real(wp) :: Th
    contains
        procedure :: eval
    end type

contains

    function init_rc_src(Th) result(this)
        real(wp), intent(in) :: Th
        type(rc_src) :: this
        this%Th = Th
    end function

    real(wp) function eval(this,u,x)
        class(rc_src), intent(in) :: this
        real(wp), intent(in) :: u
        real(wp), intent(in) :: x
        eval = this%Th**2*u
    end function

    subroutine reaction_example(Th,nx,steps,check,tol)
        real(wp), intent(in) :: Th
        integer, intent(in) :: nx
        integer, intent(in) :: steps
        integer, intent(in) :: check
        real(wp), intent(in) :: tol

        real(wp), allocatable :: uold(:), u(:)
        real(wp) :: interval(2), gre, omega, rval
        type(Lattice) :: latt

        type(rc_src), target :: mysrc

        integer :: i, step

        ! Initalize density arrays
        allocate(uold(nx),u(nx))
        uold = 0.0_wp;    
        u = uold

        ! Create source term
        mysrc = init_rc_src(Th)
        
        ! Create lattice
        omega = 1.0_wp
        latt = init_Lattice(nx=nx,uinit=uold,omega=omega,src=mysrc)

        time_loop: do step = 1, steps

            call latt%collide_and_stream(uout=u)

            ! Boundaries
            call latt%set_left_neem(1.0_wp)
            rval = (4.0_wp*u(nx-1) - u(nx-2))/3.0_wp ! one-sided finite difference
            call latt%set_right_neem(rval)
            
            ! Convergence check
            if (mod(step,check) == 0) then
                gre = sum(abs(u-uold))/sum(abs(u)) ! global relative error
                uold = u
                write(*,*) "step = ", step, "error = ", gre
                if (gre < tol) exit time_loop
            end if

            call latt%prep_next_step()
        end do time_loop

        call latt%output("reaction_results.txt")

    end subroutine

end module


program main

    use precision_mod, only: wp
    use poisson_example_mod, only: poisson_example
    use reaction_example_mod, only: reaction_example

    real(wp) :: tol

    tol = epsilon(tol) ! Set tolerance

    write(*,*) "Poisson example"
    call poisson_example(k=27.79_wp,nx=101,steps=1000000,check=5000,tol=tol)
    
    write(*,*) ""
    write(*,*) "Reaction example"
    call reaction_example(Th=1.0_wp,nx=101,steps=1000000,check=5000,tol=tol)

end program