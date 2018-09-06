module poisson_mod

    use precision_mod, only: wp
    use source_mod, only: source_term

    implicit none
    private

    public :: wp, source_term 
    public :: Lattice, init_Lattice

    ! D1Q3 Model Parameters
    integer, parameter :: e(3) = [0,1,-1]
    real(wp), parameter :: w(3) = [2._wp/3._wp,1._wp/6._wp,1._wp/6._wp]
    real(wp), parameter :: wbar(3) = [0._wp,0.5_wp,0.5_wp]
    real(wp), parameter :: alpha = 1._wp/3._wp ! Speed of sound squared $c_\mathrm{s}^2$

    real(wp), parameter :: coef = 1._wp/(1._wp - w(1)) ! coefficient in macroscopic value sum, see Eq. (2.5)


!>  Derived type to encapsulate the lattice Boltzmann algorithm
!   
    type :: Lattice
        real(wp) :: interval(2)
        integer :: nx
        real(wp) :: dx, dt, diff, omega
        real(wp), allocatable :: f(:,:,:)
        integer :: old, new
        class(source_term), pointer :: src
    contains
        procedure :: collide_and_stream
        procedure :: get_macroscopic_values
        procedure :: prep_next_step
        procedure :: output
        
        ! Non-equilibrium extrapolation
        procedure :: set_left_neem
        procedure :: set_right_neem

        ! Non-equilibrium bounceback
        ! procedure :: set_left_neeb
        ! procedure :: set_right_neeb

        ! Inamuro's boundary condition
        procedure :: set_left_inamuro
        procedure :: set_right_inamuro
    end type


contains

    function init_Lattice(ab,nx,uinit,omega,src) result(this)
        real(wp), intent(in), optional :: ab(2)
        integer, intent(in) :: nx
        real(wp), intent(in) :: uinit(nx)
        real(wp), intent(in) :: omega
        class(source_term), intent(in), pointer :: src

        type(Lattice) :: this

        real(wp) :: ab_(2), dx, tau
        integer :: i

        ! Set interval
        ab_ = [0.0_wp,1.0_wp] ! default interval
        if (present(ab)) ab_ = ab
        this%interval = ab_

        ! Set number of lattice sites
        this%nx = nx

        ! Set spatial and time step
        this%dx = (ab_(2) - ab_(1))/real(nx-1,wp)
        this%dt = this%dx

        ! Set relaxation time and diffusivity
        this%omega = omega
        tau = 1.0_wp/omega
        this%diff = alpha*(0.5_wp - tau)*this%dt ! Not the same as in standard LBM!

        ! Set flags
        this%old = 1
        this%new = 2

        ! Allocate space
        allocate(this%f(0:nx+1,3,2))
        
        ! Initialize to equilibrium values
        do i = 1, nx
            this%f(i,1,this%old) = (w(1) - 1._wp)*uinit(i)
            this%f(i,2,this%old) = w(2)*uinit(i)
            this%f(i,3,this%old) = w(3)*uinit(i)
        end do

        ! Associate source term functor
        this%src => src

    end function


    subroutine collide_and_stream(this,uout)
        class(Lattice), intent(inout) :: this
        real(wp), intent(out), optional ::uout(this%nx)
        real(wp) :: one_minus_omega, feq(3), u, source, x
        integer :: i

        associate(f => this%f, nx => this%nx, new => this%new, old => this%old, & 
            dt => this%dt, omega => this%omega, D => this%diff)

        one_minus_omega = 1._wp - omega

        spatial_loop: do i = 1, this%nx
            ! Macroscopic value, Eq. (2.5)
            u = coef*(f(i,2,old) + f(i,3,old))
            if (present(uout)) uout(i) = u

            ! Equilibirum distribution
            feq(1) = (w(1) - 1._wp)*u
            feq(2) = w(2)*u
            feq(3) = w(3)*u

            ! Evaluate source term
            x = this%interval(1) + (i-1)*this%dx
            source = this%src%eval(u,x)

            ! Collide and stream
            f(i,1,new)   = one_minus_omega*f(i,1,old) + omega*feq(1)    ! wbar(1) = 0
            f(i+1,2,new) = one_minus_omega*f(i,2,old) + omega*feq(2) + dt*wbar(2)*D*source
            f(i-1,3,new) = one_minus_omega*f(i,3,old) + omega*feq(3) + dt*wbar(3)*D*source
        end do spatial_loop

        ! Perform bounceback by default
        f(1,2,new) = f(0,3,new)
        f(nx,3,new) = f(nx+1,2,new)

        end associate

    end subroutine

!>  Sets a left Dirichlet boundary condition using the
!   non-equilibrium extrapolation method.
!
    subroutine set_left_neem(this,value)
        class(Lattice), intent(inout) :: this
        real(wp), intent(in) :: value
        real(wp) :: u
        u = coef*sum(this%f(2,2:3,this%old))
        this%f(1,1,this%new)  = (w(1) - 1._wp)*value + (this%f(2,1,this%new) - (w(1) - 1._wp)*u)
        this%f(1,2,this%new)  = w(2)*value + (this%f(2,2,this%new) - w(2)*u)
        this%f(1,3,this%new)  = w(3)*value + (this%f(2,3,this%new) - w(3)*u)
    end subroutine

!>  Sets a right Dirichlet boundary condition using the
!   non-equilibrium extrapolation method.
!
    subroutine set_right_neem(this,value)
        class(Lattice), intent(inout) :: this
        real(wp), intent(in) :: value
        real(wp) :: u
        u = coef*sum(this%f(this%nx-1,2:3,this%old))
        this%f(this%nx,1,this%new) = (w(1) - 1._wp)*value + (this%f(this%nx-1,3,this%new) - (w(1) - 1._wp)*u)
        this%f(this%nx,2,this%new) = w(2)*value + (this%f(this%nx-1,2,this%new) - w(2)*u)
        this%f(this%nx,3,this%new) = w(3)*value + (this%f(this%nx-1,3,this%new) - w(3)*u)
    end subroutine


    subroutine set_left_inamuro(this,value)
        class(Lattice), intent(inout) :: this
        real(wp), intent(in) :: value
        this%f(1,2,this%new)  = value*(1.0_wp-w(1)) - this%f(1,3,this%new)
    end subroutine

    subroutine set_right_inamuro(this,value)
        class(Lattice), intent(inout) :: this
        real(wp), intent(in) :: value
        this%f(this%nx,3,this%new) = value*(1.0_wp-w(1)) - this%f(this%nx,2,this%new)
    end subroutine

!>  Swaps flags between old and new arrays.
!
!   This routine must be called at the end of the time loop!
!
    subroutine prep_next_step(this)
        class(Lattice), intent(inout) :: this
        integer :: temp
        temp = this%old
        this%old = this%new
        this%new = temp
    end subroutine


!>  Calculates the physical quantity $u$ (macroscopic value) 
!   according to Eq. (2.5) in the paper by Chai and Shi (2008).
!
    subroutine get_macroscopic_values(this,u)
        class(Lattice), intent(in) :: this
        real(wp), intent(out) :: u(this%nx)
        u(:) = coef*sum(this%f(1:this%nx,2:3,this%old),dim=2)
    end subroutine


    subroutine output(this,fname)
        class(Lattice), intent(in) :: this
        character(len=*), intent(in) :: fname

        integer :: iunit, i
        real(wp) :: x(this%nx), u(this%nx)

        call this%get_macroscopic_values(u)

        x = [(this%interval(1) + (i-1)*this%dx,i=1,this%nx)]

        open(newunit=iunit,file=fname)

        do i = 1, this%nx
            write(iunit,*) x(i), u(i)
        end do

        close(iunit)
    end subroutine

end module

