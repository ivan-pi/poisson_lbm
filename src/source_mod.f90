module source_mod

    use precision_mod, only: wp

    implicit none
    private

    public :: source_term

    type, abstract :: source_term
    contains
        procedure(eval_interface), deferred :: eval
    end type

    interface
        real(wp) function eval_interface(this,u,x)
            import wp, source_term
            class(source_term), intent(in) :: this
            real(wp), intent(in) :: u
            real(wp), intent(in) :: x
        end function
    end interface

end module