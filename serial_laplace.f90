program serial_laplace

use constants
use initialise
use SOR_functions 
use pgm_maker

implicit none 

real(kind=dp) :: phi_0,h,k_algor,phi_grid,phi_new,phi_1,phi_2,w
real(kind=dp),allocatable,dimension(:,:) :: grid

integer :: i,j,k,r_in,r_out,l_in,l_out,ierr,counter,l,m	

open(file='convergence.dat',unit=20,status='replace')


!-------------------------------------------------
!               System Geometry 
!-------------------------------------------------

h = 0.1_dp

!do m = 1,2



!print*, 'h=',h 

l_out = int(20.0_dp/h)
r_out = int(10.0_dp/h)

l_in = int(5.0_dp/h)
r_in = int(1.0_dp/h)
print*, l_out,r_out,l_in,r_in
phi_0 = 1000
allocate(grid(1:L_out,1:r_out),stat=ierr)
if (ierr.ne.0) stop 'error allocating grid'
!print*, l_out,r_out,l_in,r_in

! step-size



!-------------------------------------------------
!                 Initalise Grid
!-------------------------------------------------
w = 1.9_dp





! check effect of w on convergence rate


counter = 0



call init(l_out,r_out,r_in,l_in,phi_0,grid)
print*, 'init sum:',sum(grid)
call write_image(l_out/2,r_out,grid)
stop
!-------------------------------------------------
!            SOR Algorithm 
!-------------------------------------------------
 

! number of grid itterations

do k = 1,1
  
  phi_1 = grid(int(6.0_dp/h),int(5.0_dp/h))

  do j = 1,l_out

    do i = 1,r_out

      

      ! are we inside the conductor 

      if ((i.lt.r_in+1).and.(j.lt.l_in+1)) then 

        grid(j,i) = grid(j,i)      

      ! not inside

      else  

        if (i.eq.1) then ! r=0 case

          ! potential at current grid point

          phi_grid = f_0(grid(j+1,1),grid(j-1,1),grid(j,2))
        

          ! new potential 

          phi_new = new_pot(grid(j,i),phi_grid,w)


          ! update grid point with new potential

          grid(j,i) = phi_new

        else ! r > 0 
     
          ! potential at current grid point
        
          phi_grid = f_1(grid(j+1,i),grid(j-1,i),grid(j,i+1),grid(j,i-1),real(i,kind=dp),h)
      
          ! new potential 

          phi_new = new_pot(grid(j,i),phi_grid,w)
    
          ! track convergence 

          ! update grid point with new potential

          grid(j,i) = phi_new

        end if 

      end if 

    end do 

  end do

  !phi_2 = grid(int(6.0_dp/h),int(5.0_dp/h))
  !if (k.gt.5000) then
   ! if ( (phi_1.gt.0).and.(phi_2.gt.0) ) then
      
    !  if (abs(phi_1-phi_2).lt.0.01) then
     !  print*, 'converged'! at k =',k
       ! print*, 'using a value of w =',w
       ! print*, 'with a value of phi =',grid(int(6/h),int(5/h))
      !  write(20,*) h,w,k
       ! exit
      !end if 
    !end if
  !end if
!print*, abs(phi_1-phi_2)
end do 
print*, sum(grid)
!h = h/10.0_dp

!end do 
!-------------------------------------------------
!            Data Output 
!-------------------------------------------------

! comparison with parallel 

!print*, sum(grid(99:199,0:99))


! output image as pgm 

call write_image(l_out/2,r_out,grid)


close(20)

deallocate(grid,stat=ierr)
if (ierr.ne.0) stop 'error deallocating grid'

end program 
