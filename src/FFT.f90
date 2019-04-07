module FFT_mod
  use, intrinsic :: iso_c_binding 
  implicit none
  include 'fftw3.f03'
contains
  subroutine FFT(in,out,n)
    double complex::in(:),out(:)
    integer::n,plan,i
    call dfftw_plan_dft_1d(plan,n,in,out,FFTW_FORWARD,FFTW_ESTIMATE)
    call dfftw_execute_dft(plan, in, out)
    call dfftw_destroy_plan(plan)
  end subroutine FFT
end module FFT_mod
