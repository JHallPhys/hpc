module pgm_maker

use constants

implicit none
  
contains

subroutine write_image(Nx,Ny,grid)
   
implicit none

real (kind=dp), dimension(:,:), intent(in) :: grid
integer,intent(in) :: Nx,Ny

integer, dimension(:,:), allocatable :: pixels
real (kind=dp) :: Tmax, Tmin

integer :: ierr, max_greys, i, j, k, out_unit 
  
allocate (pixels(1:Nx,1:Ny),stat=ierr)

if (ierr/=0) stop 'Error in allocating pixels'

max_greys=255
Tmax=maxval(grid)
Tmin=minval(grid)

do j=1,Ny
   do i=1,Nx
     pixels(i,j)=int((grid(i,j)-Tmin)*max_greys/(Tmax-Tmin)) !Tmin<T<Tmax
   end do
end do

out_unit=10

open (file='laplace_pgm',unit=out_unit,status='replace')

  write (out_unit,11) 'P2'                 !pgm magic number
  write (out_unit,12) Nx,Ny                !width, height
  write (out_unit,13) max_greys            !max gray value

do j=1,Ny

  do i=1,Nx-15,15
    write (out_unit,14) (pixels(i+k-1,j),k=1,15)  !each line < 70 chars
  end do

  write (out_unit,14) (pixels(k,j),k=i,Nx)

end do

close (unit=out_unit)

11  format(a2)
12  format(i10,1x,i10)
13  format(i10)
14  format(15(1x,i3))


end subroutine write_image
  
end module 
