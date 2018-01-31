      module mod_qcdloop_c
        implicit none

      interface
        function cln(x,isig) bind(C,name="cln")
          use iso_c_binding
          complex(c_double_complex), intent(in) :: x
          real(c_double), intent(in) :: isig
          complex(c_double_complex) :: cln
        end function

        function qlzero(x) bind(C,name="qlzero")
          use iso_c_binding
          real(c_double), intent(in) :: x
          logical(c_bool) :: qlzero
        end function

        function qlnonzero(x) bind(C,name="qlnonzero")
          use iso_c_binding
          real(c_double), intent(in) :: x
          logical(c_bool) :: qlnonzero
        end function

        function qli1(m1,mu2,ep) bind(C,name="qli1")
          use iso_c_binding
          real(c_double), intent(in) :: m1,mu2
          integer(c_int), intent(in) :: ep
          complex(c_double_complex) :: qli1
        end function

        function qli2(p1,m1,m2,mu2,ep) bind(C,name="qli2")
          use iso_c_binding
          real(c_double), intent(in) :: p1,m1,m2,mu2
          integer(c_int), intent(in) :: ep
          complex(c_double_complex) :: qli2
        end function

        function qli3(p1,p2,p3,m1,m2,m3,mu2,ep) bind(C,name="qli3")
          use iso_c_binding
          real(c_double), intent(in) :: p1,p2,p3,m1,m2,m3,mu2
          integer(c_int), intent(in) :: ep
          complex(c_double_complex) :: qli3
        end function

        function qli4(p1,p2,p3,p4,s12,s23,m1,m2,m3,m4,mu2,ep)
     &      bind(C,name="qli4")
          use iso_c_binding
          real(c_double), intent(in) :: p1,p2,p3,p4,s12,s23
          real(c_double), intent(in) :: m1,m2,m3,m4,mu2
          integer(c_int), intent(in) :: ep
          complex(c_double_complex) :: qli4
        end function

        ! I consider the quad precision interface to be broken:
        ! it probably only works with gfortran since there are no
        ! standard c bindings for quad precision types
        ! (on gfortran there is c_float128 and c_float128_complex)

        subroutine qli1q(m1,mu2,ep,ret) bind(C,name="qli1q")
          use iso_c_binding
          real*16, intent(in) :: m1,mu2
          integer(c_int), intent(in) :: ep
          complex*32, intent(out) :: ret
        end subroutine

        subroutine qli2q(p1,m1,m2,mu2,ep,ret) bind(C,name="qli2q")
          use iso_c_binding
          real*16, intent(in) :: p1,m1,m2,mu2
          integer(c_int), intent(in) :: ep
          complex*32, intent(out) :: ret
        end subroutine

        subroutine qli3q(p1,p2,p3,m1,m2,m3,mu2,ep,ret)
     &      bind(C,name="qli3q")
          use iso_c_binding
          real*16, intent(in) :: p1,p2,p3,m1,m2,m3,mu2
          integer(c_int), intent(in) :: ep
          complex*32, intent(out) :: ret
        end subroutine

        subroutine qli4q(p1,p2,p3,p4,s12,s23,m1,m2,m3,m4,mu2,ep,ret)
     &      bind(C,name="qli4q")
          use iso_c_binding
          real*16, intent(in) :: p1,p2,p3,p4,s12,s23
          real*16, intent(in) :: m1,m2,m3,m4,mu2
          integer(c_int), intent(in) :: ep
          complex*32, intent(out) :: ret
        end subroutine

      end interface

      end module
