module timing_module

  implicit none
  integer, dimension (8), private :: dt
  real, private :: h, m, s, ms, tt
  real, private :: last_tt

contains

  subroutine start_timing()
    implicit none

    call date_and_time(values=dt)
    print 100, dt(1:3), dt(5:8)
    h = real(dt(5))
    m = real(dt(6))
    s = real(dt(7))
    ms = real(dt(8))
    last_tt = 60*(60*h+m) + s + ms/1000.0
100 format (1x, i4, '/', i2, '/', i2, 1x, i2, ':', i2, ':', i2, 1x, i3)
  end subroutine start_timing

  subroutine end_timing()
    implicit none

    call date_and_time(values=dt)
    print 100, dt(1:3), dt(5:8)
100 format (1x, i4, '/', i2, '/', i2, 1x, i2, ':', i2, ':', i2, 1x, i3)
  end subroutine end_timing

  real function time_difference()
    implicit none

    tt = 0.0
    call date_and_time(values=dt)
    h = real(dt(5))
    m = real(dt(6))
    s = real(dt(7))
    ms = real(dt(8))
    tt = 60*(60*h+m) + s + ms/1000.0
    time_difference = tt - last_tt
    last_tt = tt
  end function time_difference

end module timing_module
