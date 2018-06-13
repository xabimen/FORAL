module read

public :: read_vasp






contains

subroutine read_vasp(fileUnit,cell,ntype,coordinates,forces,strucname,io)
implicit none
integer, intent(in)                             :: fileUnit
real, dimension(:,:), intent(out)               :: cell
integer, dimension(:), allocatable, intent(out) :: ntype
real, dimension(:,:), allocatable, intent(out):: coordinates, forces
integer, intent(out)                            :: io
real, dimension(3,3)                          :: cell_inv
character(len=15), intent(out)                  :: strucname
integer                                         :: i
integer, dimension(2)                           :: atomtype, numIons
character(len=2), dimension(size(atomType))     :: atomName



read(unit=fileUnit,iostat=io,fmt='(a15)') strucname
if (io == -1) then
    RETURN
endif
read(unit=fileUnit,fmt=*) 

do i = 1, 3
    read(unit=fileUnit,fmt=*) cell(i,:)
enddo

read(unit=fileUnit,fmt=*) atomName(:)

if (atomName(1)=="Si") then
    atomType(1) = 14
elseif (atomName(1) == "Ca") then
    atomType(1) = 20
else
    print*, "Error: reading line 6 of vasp file", atomName(1)
endif
if (atomName(2)=="O ") then
    atomType(2) = 8
else
    print*, "Error: reading line 6 of vasp file", atomName(1),"  and   ", atomName(2)
endif

read(unit=fileUnit,fmt=*) numIons(:)
read(unit=fileUnit,fmt=*)
allocate(coordinates(sum(numIons),3),forces(sum(numIons),3))

call inverse(cell,cell_inv,3)

allocate(ntype(sum(numIons)))
ntype = 0
do i = 1, numions(1)
    ntype(i) = atomType(1)
    read(unit=fileUnit,fmt=*) coordinates(i,:), forces(i,:)
    coordinates(i,:) = matmul(cell_inv,coordinates(i,:))
enddo
do i = numIons(1)+1,numIons(1)+numIons(2)
    ntype(i) = atomType(2)
    read(unit=fileUnit,fmt=*) coordinates(i,:), forces(i,:)
    coordinates(i,:) = matmul(cell_inv,coordinates(i,:))
enddo

end subroutine

  subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
implicit none 
integer n
real a(n,n), c(n,n)
real L(n,n), U(n,n), b(n), d(n), x(n)
real coeff
integer i, j, k

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0

! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
end subroutine inverse



end module
