module precision_mod

    implicit none
    private

    public :: wp

    integer, parameter :: sp = kind(1.0e0)  ! single precision
    integer, parameter :: dp = kind(1.0d0)  ! double precision

    integer, parameter :: wp = dp   ! working precison

end module