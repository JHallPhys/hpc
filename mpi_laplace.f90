program serial_laplace

use constants
use SOR_functions 
use pgm_maker
use mpi

implicit none 

real(kind=dp) :: phi_0,h,k_algor,phi_grid,phi_old,phi_new,w
real(kind=dp) :: time_start,time_stop,eps,a,b
real(kind=dp),allocatable,dimension(:,:) :: grid
integer :: i,j,k,r_in,r_out,l_in,l_out,ierror,max_itteration,split

! MPI Variables 


integer :: my_rank,num_procs,send_req_l,send_req_r, &
           & recv_req_l,recv_req_r,imax,imin


integer :: new_type

integer,dimension(1:MPI_STATUS_SIZE) :: send_rstat,recv_rstat, &
                                       & send_lstat,recv_lstat



!##############################################################
!-------------------------------------------------------------
!                    MPI - Stuff
!-------------------------------------------------------------
!##############################################################

! create the world

call MPI_init(ierror)
if (ierror.ne.0) stop 'Error in MPI_init'

! get rank on each processor in world: rank=0,nprocs-1 

call MPI_comm_rank(MPI_comm_world,my_rank,ierror)
if (ierror.ne.0) stop 'Error in MPI_comm_rank'


! Get total number of processors

call MPI_comm_size(MPI_comm_world,num_procs,ierror)
if (ierror.ne.0) stop 'Error in MPI_comm_size'


! MPI info for user 

!print*,
!print*, '--------- MPI INFO --------------'
!print*,
!print*, 'Number of processors:', num_procs
!print*, 'Rank of current:', my_rank
!print*, 
!print*, '---------------------------------'

!-------------------------------------------------
!            INITISIALISE SYSTEM 
!-------------------------------------------------


! System Geometry 


h = 0.1_dp ! step-size

l_out = int(20.0_dp/h)
r_out = int(10.0_dp/h)

l_in = int(5.0_dp/h)
r_in = int(1.0_dp/h)

phi_0 = 1000.0_dp
!print*, l_out,r_out,l_in,r_in
! SOR algor variables 

w = 1.9_dp 
phi_old = 0.0_dp
eps = 0.001_dp


! progress tracking params 

time_start = 0.0_dp
time_stop = 0.0_dp


! start timer

time_start = MPI_wtime()

!-------------------------------------------------------------
!       Initilaise  parallel qunatities
!-------------------------------------------------------------

! split the data across the processors 

! we require the number of processors to be an integer devisor of l_out 

! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! come back to me 
!if (
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 

split = (l_out)/num_procs

! assign index range in terms of rank 

imin = 1 + (my_rank*split)

imax = split + imin  

! check its working 

!print*, 
!print*, '#############################################'
!print*, '----------- Data Spliiting info -------------'
!print*, '============================================='
!write(*,'(1x,a8,1x,i2,a1,1x,a6,1x,i3,2x,a1,1x,a6,1x,i3)') 'for rank',my_rank,':','imin =',imin,';','imax =',imax
!print*, '============================================='
!print*, 

! allocate mpi grid and split data - including ghost cells  

allocate(grid(0:split+1,0:r_out+1),stat=ierror)
if (ierror.ne.0) stop 'error allocating mpi grid'



!--------------------------------------------------------------
! set up charge distribution within system 
!--------------------------------------------------------------

! outer surface grounded 

grid(:,:) = 0.0_dp

! set up reference plane charge using potential of infinite coax 

if (my_rank.eq.0) then 

  do i = r_in+1,r_out  

    grid(1,i) = phi_0*( (log(real(r_out,kind=dp)) - log(real(i,kind=dp)))/(log(real(r_out,kind=dp)) - log(real(r_in,kind=dp))) )

  end do 

end if 

! set up inner cylinder - may cross more than one grid slice 

if (imax.lt.l_in) then 

   grid(1:split,1:r_in) = phi_0

else if ( (imin.lt.l_in).and.(imax.gt.l_in) ) then 

   grid(1:(l_in-imin+1),1:r_in) = phi_0

end if 


! define new data type to allow passing of array rows

call MPI_Type_vector(100,1,split+2,MPI_DOUBLE_PRECISION,new_type,ierror)
if (ierror.ne.0) stop 'problem in mpi_type_vector'

call MPI_TYPE_COMMIT(new_type,ierror)
if (ierror.ne.0) stop 'problem with type commit'

print*, sum(grid)

!##############################################################
!--------------------------------------------------------------
!               SOR Algorithm - main loop 
!--------------------------------------------------------------
!##############################################################


! number of grid itterations

do k = 1,1
  
  !-------------------------------------------------
  !        split data across the world 
  !-------------------------------------------------
 
  if (my_rank < num_procs-1) then
    
    
    call MPI_issend(grid(split,1), 1, new_type, &
         & (my_rank+1), 1, MPI_comm_world, send_req_r, ierror)
    if (ierror.ne.0) stop 'Error in send right'

    
    call MPI_irecv(grid(split+1,1),1, new_type, &
         & (my_rank+1), 2, MPI_comm_world, recv_req_r, ierror)
    if (ierror.ne.0) stop 'Error in recev right'

  end if

  if (my_rank > 0) then

    call MPI_irecv(grid(0,1),1,new_type, &
          & (my_rank-1), 1, MPI_comm_world,recv_req_l, ierror)
    if (ierror.ne.0) stop 'Error in recev left'

    

    call MPI_issend(grid(1,1),1, new_type, &
          & (my_rank-1), 2, MPI_comm_world, send_req_l, ierror)
    if (ierror.ne.0) stop 'Error in send left'
    
  end if

  
  !##############################################################
  !          update non-halo region 
  !-------------------------------------------------
  ! Cases: 1) Conducter in slice.
  !        2) Slice free of conductor 
  !-------------------------------------------------
  !##############################################################

  ! free of conducter slice is << Yoda comment   
  
  if (imin.gt.l_in) then
    
    do j = 2,split-1
      do i = 1,r_out-1
       
        if (i.gt.1) then ! r =/= case

        ! potential at current grid point
        phi_grid = f_1(grid(j+1,i),grid(j-1,i),grid(j,i+1),grid(j,i-1),real(i,kind=dp),h)
 
        ! new potential 
        grid(j,i) = new_pot(grid(j,i),phi_grid,w)
  

        else ! r=0 case

          phi_grid = f_0(grid(j+1,1),grid(j-1,1),grid(j,2)) 

          grid(j,1) = new_pot(grid(j,1),phi_grid,w)

        end if 

      end do 
    end do   


  ! conductor in slice
  
  else if ((imin.le.l_in).and.(imax.gt.l_in)) then 

    do j = 2,split-1
      do i = 1,r_out-1
        
        if ( (j.le.(l_in-imin)+1).and.(i.le.r_in) ) then  

          grid(j,i)=grid(j,i)             

        else if (j.gt.(l_in-imin)+1) then 

          if(i.gt.1) then
            phi_grid = f_1(grid(j+1,i),grid(j-1,i),grid(j,i+1),grid(j,i-1),real(i,kind=dp),h) 
            grid(j,i) = new_pot(grid(j,i),phi_grid,w)
            
          else
            phi_grid = f_0(grid(j+1,1),grid(j-1,1),grid(j,2)) 
            grid(j,i) = new_pot(grid(j,i),phi_grid,w)
          
          end if 
        
        else if ( (j.le.(l_in-imin)+1).and.(i.gt.r_in) ) then  

          if(i.gt.1) then
            phi_grid = f_1(grid(j+1,i),grid(j-1,i),grid(j,i+1),grid(j,i-1),real(i,kind=dp),h) 
            grid(j,i) = new_pot(grid(j,i),phi_grid,w)            
          else
            phi_grid = f_0(grid(j+1,1),grid(j-1,1),grid(j,2)) 
            grid(j,i) = new_pot(grid(j,i),phi_grid,w)           
          end if 

        end if 
  
      end do 
    end do 

  end if 


  !##############################################################
  !-------------------------------------------------
  !            Swap the halo data
  !-------------------------------------------------
  !##############################################################


  if (my_rank < num_procs-1) then

    call MPI_wait(send_req_r, send_rstat, ierror)
    if (ierror.ne.0) stop 'Error in wait: send right'

    call MPI_wait(recv_req_r, recv_rstat, ierror)
    if (ierror.ne.0) stop 'Error in wait: recieve right'

  end if

     
  if (my_rank > 0) then

    call MPI_wait(recv_req_l, recv_lstat, ierror)
    if (ierror.ne.0) stop 'Error in wait: recieve left'


    call MPI_wait(send_req_l, send_lstat, ierror)
    if (ierror.ne.0) stop 'Error in wait: send left'
    
    
  end if

 
  ! conducter not in data slice
  if (imin.gt.l_in) then
    j=1
    phi_grid = f_0(grid(j+1,1),grid(j-1,1),grid(j,2)) 
    grid(j,1) = new_pot(grid(j,1),phi_grid,w)

    j=split
    phi_grid = f_0(grid(split+1,1),grid(split-1,1),grid(split,2)) 
    grid(j,1) = new_pot(grid(j,1),phi_grid,w)
    
    
    do i = 2,r_out-1


      j=1
      phi_grid = f_1(grid(j+1,i),grid(j-1,i),grid(j,i+1),grid(j,i-1),real(i,kind=dp),h) 
      grid(j,i) = new_pot(grid(j,i),phi_grid,w)


      j=split
      phi_grid = f_1(grid(j+1,i),grid(j-1,i),grid(j,i+1),grid(j,i-1),real(i,kind=dp),h) 
      grid(j,i)  = new_pot(grid(j,i),phi_grid,w)

    end do 


  ! conducter in data slice
  else if ((imin.le.l_in).and.(imax.gt.l_in)) then

    do j = 1,split,split-1

      do i = 1,r_out-1

        if ( (j.le.l_in-imin+1).and.(i.gt.r_in) ) then

            if ((my_rank.eq.0).and.(j.eq.1)) then ! maintain reference B.C
              grid(j,i) = grid(j,i)
            else            
              phi_grid = f_1(grid(j+1,i),grid(j-1,i),grid(j,i+1),grid(j,i-1),real(i,kind=dp),h) 
              grid(j,i) = new_pot(grid(j,i),phi_grid,w)
            end if   

        else if (j.gt.l_in-imin+1) then 
    
          phi_grid = f_1(grid(j+1,i),grid(j-1,i),grid(j,i+1),grid(j,i-1),real(i,kind=dp),h) 
          grid(j,i) = new_pot(grid(j,i),phi_grid,w)
           
        end if 

      end do 

      
    end do 
   
  else if (imax.lt.l_in) then

    do i = r_in+1,r_out-1

      j=1
      phi_grid = f_1(grid(j+1,i),grid(j-1,i),grid(j,i+1),grid(j,i-1),real(i,kind=dp),h) 
      grid(j,i) = new_pot(grid(j,i),phi_grid,w)


      j=split
      phi_grid = f_1(grid(j+1,i),grid(j-1,i),grid(j,i+1),grid(j,i-1),real(i,kind=dp),h) 
      grid(j,i)  = new_pot(grid(j,i),phi_grid,w)
    
    end do  

  end if 

end do 

! stop timer 

time_stop = MPI_wtime()

print*, 'total time on proc', my_rank, 'is', time_stop - time_start

!==================================================
!--------------------------------------------------
!            Data Output 
!--------------------------------------------------
!==================================================



if(my_rank.eq.0) then 

  !call write_image(100,100,grid)
  print*, 'sum of rank 0 array:',sum(grid)

end if 

if (my_rank.eq.1) then 

  print*, 'sum of rank 1 array:',sum(grid)
  !call write_image(100,100,grid)
  
end if

!++++++++++++++++++++++++++++++++++++++++++++++++++
!==================================================


deallocate(grid,stat=ierror)
if (ierror.ne.0) stop 'error deallocating grid'


! end mpi
call MPI_finalize(ierror)
if (ierror.ne.0) stop 'Error in MPI_finalize'




end program 
