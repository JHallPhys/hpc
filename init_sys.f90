module initialise 
use constants

implicit none 


contains 

! initalise potentials on referance plane and inner cylinder

subroutine init(z_max,r_max,r_in,z_in,phi_0,grid)

! args 

integer,intent(in) :: z_max,r_max,z_in,r_in

real(kind=dp),intent(in) :: phi_0
real(kind=dp),dimension(0:z_max,0:r_max),intent(out) :: grid 

! local vars 

real(kind=dp) :: phi_ref
integer :: i,j 

!---------------------------------------------------------------
! set up grid with refernce plane potentials and inner potential 
!---------------------------------------------------------------


grid(:,:) = 0.0_dp

! set up potential on surface on the inner cylinder

do i = 0,r_in
  do j = 0,z_in

  grid(j,i) = phi_0

  end do 
end do  

! set up potential on the refernce plane 

do i = r_in+1,r_max

  phi_ref = phi_0*( (log(real(r_max,kind=dp)) - log(real(i,kind=dp)))/(log(real(r_max,kind=dp)) - log(real(r_in,kind=dp))) )
  
  grid(0,i) = phi_ref

end do 

end subroutine 

end module 
