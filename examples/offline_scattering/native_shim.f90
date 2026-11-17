! Copyright (c) 2026 Stardust Labs. MIT license; see repository LICENSE.
! Test-only C bridge to the repository's unmodified vertical adding routine.
subroutine evaluate_adding(n, r, rd, t, td, b, fd, fu) bind(c)
  use iso_c_binding
  use rrtmg_sw_vrtqdr, only: vrtqdr_sw
  implicit none
  integer(c_int), value :: n
  real(c_double), intent(in) :: r(n+1), rd(n+1), t(n+1), td(n+1), b(n+1)
  real(c_double), intent(out) :: fd(n+1), fu(n+1)
  real(c_double) :: down(n+1), up(n+1), upd(n+1), beam(n+1)
  real(c_double) :: d(n+1,1), u(n+1,1)
  integer :: j
  down = 0; up = 0; upd = 0; beam = 1
  up(n+1) = r(n+1); upd(n+1) = rd(n+1)
  do j = 1, n
    beam(j+1) = beam(j)*b(j)
  end do
  call vrtqdr_sw(n, 1, r, rd, t, td, b, down, up, upd, beam, d, u)
  fd = d(:,1); fu = u(:,1)
end subroutine evaluate_adding
