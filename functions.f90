module SOR_functions

use constants

implicit none 

contains 


! r>0 case

function f_1(i_plus,i_minus,k_plus,k_minus,pos_current,h_step)

real(kind=dp) :: f_1 
real(kind=dp),intent(in) :: i_plus,i_minus,k_plus,k_minus,pos_current,h_step

f_1 = 0.25_dp*(i_plus+i_minus+k_plus+k_minus) + (h_step/(8.0_dp*pos_current))*(k_plus-k_minus)

end function f_1


! r=0 case 

function f_0(i_plus,i_minus,i)

real(kind=dp) :: f_0
real(kind=dp),intent(in) :: i,i_plus,i_minus

f_0 = (2.0_dp/3.0_dp)*(i) + (i_plus + i_minus)/6.0_dp

end function f_0

! New potential 

function new_pot(phi_old,U,w) 

real(kind=dp) :: new_pot 
real(kind=dp) :: phi_old,U,w

new_pot = phi_old + w*(U-phi_old) 

end function new_pot

end module 
