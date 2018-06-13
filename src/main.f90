program main

use read
use fingerprint

implicit none

real, dimension(3,3)                            :: cell                         
real, dimension(:,:), allocatable               :: coordinates, forces
real, dimension(:,:,:), allocatable             :: fprint
integer, dimension(:), allocatable              :: ntype
integer                                         :: io
character(len=30)                               :: strucname
integer                                         :: i, k

open(unit=123,file='test.vasp',status='old',action='read')

call read_vasp(123,cell,ntype,coordinates,forces,strucname,io)


allocate(fprint(sum(ntype),0:2*50,3))
call fp(coordinates,cell,ntype,2,50,10.0,fprint)

do i = 1, size(ntype)
    print*, "atom", i
    do k = 0, 50
        print*, fprint(i,k,:)
    enddo
enddo

end program
