module experiment
  use global
  use mesh_mod
  use IO
contains
  ! ----------------------------------------------------------------------------
  subroutine calc_sph_harm(mesh,l,m,SH)
    implicit none
    integer::l,m,i,j
    type(t_mesh)::mesh
    double precision::theta,phi,x,y,z,r
    !    character (len=1024) :: filename
    double complex:: SH(:)
    do i=1,mesh%Ntot
       x=mesh%node(i)%q(1)-mesh%box%center(1)
       y=mesh%node(i)%q(2)-mesh%box%center(2)
       z=mesh%node(i)%q(3)-mesh%box%center(3)
       r=sqrt(x*x+y*y+z*z)
       phi=calc_phi(x,y)
       theta=calc_theta(x,y,z)
       call sph_harm(l,m,theta,phi,SH(i))
       if((z.eq.0.0).and.(l.eq.1).and.(m.eq.1)) then
          print *,'calc_sph_harm(',l,',',m,')=',SH(i),&
               -sqrt(3.0*0.5*0.25/pi)*sin(theta)*cos(phi)-dreal(SH(i)),&
               -sqrt(3.0*0.5*0.25/pi)*sin(theta)*sin(phi)-dimag(SH(i))!,&               theta*180/pi,phi*180/pi
       end if
       if((z.eq.0.0).and.(l.eq.1).and.(m.eq.-1)) then
          print *,'calc_sph_harm(',l,',',m,')=',SH(i),&
               sqrt(3.0*0.5*0.25/pi)*sin(theta)*cos(phi)-dreal(SH(i)),&
               -sqrt(3.0*0.5*0.25/pi)*sin(theta)*sin(phi)-dimag(SH(i))!,&               theta*180/pi,phi*180/pi
       end if
!       SH(i)=SH(i)/r
!       print *,x,y,z,r,theta,phi,SH(i)
!       print *,x,y,z,180*theta/pi,180*phi/pi
!       print *,x,y,z,180*theta/pi,cos(theta)
    end do
!    print *,SH,abs(SH)


  end subroutine calc_sph_harm
end module experiment
