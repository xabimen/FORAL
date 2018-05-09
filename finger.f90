!****************************************************
!*                                                  *    
!*                  FINGERPRINT                     *                
!*                   F.O.R.A.L                      *
!*                                                  *
!*   Author:    Xabier M. Aretxabaleta              *
!*                 2018, EHU/UPV                    *
!*                                                  *
!****************************************************

! Programa honek input moduan VASP fitxategi bat hartzen
! du kartestar kordenatuetan eta egitura horren
! fingerprint-a kalkulatzen du.

program fingerp


implicit none

real*8, parameter                               :: pi = acos(-1.0d0)
real*8                                          :: Rmax, V, r,start, w
real*8, dimension(3)                            :: r_vec
real*8, dimension(3,3)                          :: cell
real*8, dimension(:), allocatable               :: a
real*8, dimension(:,:), allocatable             :: coordinates, dist_matrix, forces
real*8, dimension(:,:,:), allocatable           :: fingerprint
integer                                         :: n_inp, i, natoms, d, k, j, j2, cont, io
integer, dimension(:), allocatable              :: atomType, numIons, N, typ_i, typ_j
integer, dimension(:,:), allocatable            :: dist_atoms, atom_list
character(len=15)                               :: strucname
character(len=35)                               :: input_file, output_file, arg1, arg2, arg3

call get_command_argument(1,arg1)
read(arg1,*) n_inp
call get_command_argument(2,arg2)
read(arg2,*) d
call get_command_argument(3,arg3)
read(arg3,*) Rmax
call get_command_argument(4,input_file)
call get_command_argument(5,output_file)
allocate(numions(n_inp),atomtype(n_inp))

open(unit=123, file=input_file, action='read', status='old')

open(unit=124,file=output_file,status='replace',action='write')

do

call read_vasp(123,cell,atomType,numIons,coordinates,forces,strucname,io)

if (io == -1) then
    exit
endif

natoms = sum(numIons)
call makeMatrices(cell,coordinates,numIons,atomType,Rmax,N,V,dist_matrix,dist_atoms)


!"do i = 1, size(atomType)
!    print*, numIons(i)
!enddo
!do i = 1, size(dist_atoms,1)
!    print*, dist_atoms(i,:), dist_matrix(i,:), N(i)
!enddo


allocate(fingerprint(sum(numIons),0:d,3),a(d),atom_list(sum(numIons),2))

cont = 1
do i = 1, 2
    do j = 1, numIons(i)
        atom_list(cont,1) = cont
        atom_list(cont,2) = atomType(i)
        cont = cont+1
    enddo
enddo

start = minval(sqrt(sum(dist_matrix(2:,:)**2,2)))/2
w = (Rmax - start)/real(d-1,8)

a = 0.0d0
do k = 1, d
       a(k) = (Rmax - start)*real(k-1)/real(d-1,8) + start
enddo

fingerprint = 0.0d0
do i = 1, size(dist_atoms,1)
    r_vec = dist_matrix(i,:)
    r = sqrt(sum(r_vec**2))
    j = dist_atoms(i,2)
 !   print*, "r", r_vec, r
    do k = 1, d
        fingerprint(j,k,:) = fingerprint(j,k,:) + &
        (r_vec/r)*exp(-(r-a(k))**2/(2.0d0*w**2))/(sqrt(2.0d0*pi)*w)*(cos(pi*r/Rmax)+1.0d0)/2.0d0*N(i)
    enddo
enddo


write(*,fmt='(a,a)') strucName,"-(r)en fingerprint-a kalkulatzen"
do i = 1, sum(numIons)
    write(unit=124, fmt='(a,a,i3)') strucname, "atom ", i
    fingerprint(i,0,:) = atom_list(i,2)
    do k = 0, d
        write(unit=124, fmt='(4f20.10)') a(k), fingerprint(i,k,:)
    enddo
    write(unit=124, fmt='(a,3f20.10)') "FORCE (eV/Angstom): ", forces(i,:)
enddo

deallocate(fingerprint,a,atom_list,coordinates,forces)
enddo

print*,
print*, "Fingerprint-a ",output_file," fitxategian gorde da."

contains


subroutine read_vasp(fileUnit,cell,atomType,numIons,coordinates,forces,strucname,io)
integer, intent(in)                             :: fileUnit
real*8, dimension(:,:), intent(out)             :: cell
integer, dimension(:), intent(out)              :: atomType
integer, dimension(:), intent(out)              :: numIons
real*8, dimension(:,:), allocatable, intent(out):: coordinates, forces
integer, intent(out)                            :: io
real*8, dimension(3,3)                          :: cell_inv
character(len=15), intent(out)                  :: strucname
integer                                         :: i
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

do i = 1, sum(numIons)
    read(unit=fileUnit,fmt=*) coordinates(i,:), forces(i,:)
    coordinates(i,:) = matmul(cell_inv,coordinates(i,:))
enddo

end subroutine

subroutine makeMatrices(cell,coordinates,numIons,atomType,Rmax,N,V,dist_matrix,dist_atoms)

implicit none

real*8, dimension(:,:), intent(in)              :: cell, coordinates
integer, dimension(:), intent(in)               :: atomtype, numIons
real*8, intent(in)                              :: Rmax
integer, dimension(:), allocatable, intent(out) :: N
real*8, intent(out)                             :: V
real*8, dimension(:,:), allocatable, intent(out):: dist_matrix
integer                                         :: species, i, j, c, lengthX,&
                                                   lengthY, lengthZ,&
                                                   quit_marker_x, quit_marker_y, &
                                                   quit_marker_z, k, quadrants,&
                                                   current_cell, natoms,&
                                                   marker, basic_cell, Nfull, &
                                                   c2, c3
integer, dimension(sum(numIons))                :: atomType1
real*8, dimension(sum(numIons))                 :: distances
real*8, dimension(13,3)                         :: vect
real*8, dimension(13)                           :: abs_vect
integer, dimension(8,3)                         :: condition, signum
real*8                                          :: x, y, z, Rij
integer, dimension(100000)                      :: N2
integer, dimension(100000,2)                    :: dist_atoms2
integer, dimension(:,:), allocatable, intent(out)   :: dist_atoms
real*8, dimension(100000,3)                     :: dist_matrix2

V = abs(det(cell))
species = size(numIons)
natoms = sum(numIons)
N2 = 0
distances = 0
dist_matrix2 = 0
c = 0
c2 = 1

do i = 1, species
    do j = 1, numIons(i)
        c = c + 1
        atomType1(c) = atomType(i)
    enddo
enddo

c = 0

condition(1,:) = (/0, 0, 0/)
condition(2,:) = (/1, 0, 0/)
condition(3,:) = (/0, 1, 0/)
condition(4,:) = (/0, 0, 1/)
condition(5,:) = (/1, 1, 0/)
condition(6,:) = (/0, 1, 1/)
condition(7,:) = (/1, 0, 1/)
condition(8,:) = (/1, 1, 1/)
signum(1,:) = (/1, 1, 1/)
signum(2,:) = (/-1, 1, 1/)
signum(3,:) = (/1, -1, 1/)
signum(4,:) = (/1, 1, -1/)
signum(5,:) = (/-1, -1, 1/)
signum(6,:) = (/1, -1, -1/)
signum(7,:) = (/-1, 1, -1/)
signum(8,:) = (/-1, -1, -1/)
vect(1,:) = cell(1,:);
vect(2,:) = cell(2,:);
vect(3,:) = cell(3,:);
vect(4,:) = vect(1,:)+vect(2,:);
vect(5,:) = vect(1,:)-vect(2,:);
vect(6,:) = vect(1,:)+vect(3,:);
vect(7,:) = vect(1,:)-vect(3,:);
vect(8,:) = vect(3,:)+vect(2,:);
vect(9,:) = vect(3,:)-vect(2,:);
vect(10,:) = vect(1,:)+vect(2,:)+vect(3,:);
vect(11,:) = vect(1,:)+vect(2,:)-vect(3,:);
vect(12,:) = vect(1,:)-vect(2,:)+vect(3,:);
vect(13,:) = -vect(1,:)+vect(2,:)+vect(3,:);

do i = 1, 13
    abs_vect(i) = sqrt(dot_product(vect(i,:),vect(i,:)))
enddo

lengthX = ceiling((Rmax+maxval(abs_vect))/minval(abs_vect)) + 1
lengthY = lengthX
lengthZ = lengthX
do i = 0, lengthX
    quit_marker_x = 1
    do j = 0, lengthY
        quit_marker_y = 1
        do k = 0, lengthZ
            quit_marker_z = 1
            do quadrants = 1, 8
                if (condition(quadrants,1)*iszero(i) + condition(quadrants,2)&
                    *iszero(j) + condition(quadrants,3)*iszero(k) == 0) then
                    do current_cell = 1, natoms
                        distances = 0.0
                        marker = 0
                        do basic_cell = 1, natoms
                            x = coordinates(current_cell,1) + &
                            signum(quadrants,1)*i - &
                            coordinates(basic_cell,1)
                            y = coordinates(current_cell,2) + &
                            signum(quadrants,2)*j - &
                            coordinates(basic_cell,2)
                            z = coordinates(current_cell,3) + &
                            signum(quadrants,3)*k - &
                            coordinates(basic_cell,3)

                            Rij = (x*cell(1,1)+y*cell(2,1)+z*cell(3,1))**2
                            Rij = Rij + (x*cell(1,2)+y*cell(2,2)+z*cell(3,2))**2
                            Rij = Rij + (x*cell(1,3)+y*cell(2,3)+z*cell(3,3))**2
                            !print*, k, quadrants, current_cell, basic_cell, Rij
                            !print*, x, y, z
                            if (Rij < Rmax**2 .and. Rij > 0.0001) then
                                quit_marker_z = 0
                                quit_marker_y = 0
                                quit_marker_x = 0


                                !print*, i, j, k, current_cell, basic_cell, matmul((/x,y,z/),cell), atomType1(current_cell)
                                !if (marker == 0 .and. Nfull >= natoms) then
                                c = c + 1 
                                    !print*, "barruan"
                                N2(c) = atomType1(current_cell) 
                                dist_atoms2(c,1) = current_cell
                                dist_atoms2(c,2) = basic_cell
                                dist_matrix2(c,:) = matmul((/x,y,z/),cell)

                                !endif
                                !Nfull = Nfull + 1 - marker
                                marker = 1
                                distances(basic_cell) = sqrt(Rij)
                            endif
                        enddo
                        !if (i+j+k+current_cell == 1) then
                        !    dist_matrix2(1,:) = distances(:)
                        !    typ_j2(1) = typ_j_coef
                        !    c3 = 1
                        !elseif (marker == 1) then
                        !    c2 = c2 + 1
                        !    dist_matrix2(c2,:) = distances(:)
                        !    c3 = c3 + 1
                        !    typ_j2(c3) = typ_j_coef
                        !endif
                    enddo
                endif
            enddo
            if (quit_marker_z == 1) then
                exit
            endif
        enddo
        if (quit_marker_y == 1) then
            exit
        endif
    enddo
    if (quit_marker_x == 1) then
        exit
    endif
enddo

allocate(N(c),dist_matrix(c,3),dist_atoms(c,2))

do i = 1, c
    N(i) = N2(i)
    dist_matrix(i,:) = dist_matrix2(i,:)
    dist_atoms(i,:) = dist_atoms2(i,:)
enddo


end subroutine

SUBROUTINE COSINE_DISTANCE(FING,FING2,NUMIONS,COS_DIST)
IMPLICIT NONE
REAL*8, DIMENSION(:,:), INTENT(IN)          :: FING, FING2
INTEGER, DIMENSION(:), INTENT(IN)           :: NUMIONS
REAL*8, INTENT(OUT)                         :: COS_DIST
INTEGER                                     :: I, J, K, BINS
REAL*8                                      :: A, B, C

A = 0.0D0
B = 0.0D0
c = 0.0D0
DO I = 1, SIZE(NUMIONS)
    DO J = 1, SIZE(NUMIONS)
        DO BINS = 1, SIZE(FING,2)
            K = (I-1)*2 + J
            A = A + FING(K,BINS)*FING2(K,BINS)*NUMIONS(I)*NUMIONS(J)/REAL(SUM(NUMIONS))      
            B = B + FING(K,BINS)**2*NUMIONS(I)*NUMIONS(J)/REAL(SUM(NUMIONS))
            C = C + FING2(K,BINS)**2*NUMIONS(I)*NUMIONS(J)/REAL(SUM(NUMIONS))
        ENDDO
    ENDDO
ENDDO

COS_DIST = (1 - A/(SQRT(B)*SQRT(C)))/2.0D0

END SUBROUTINE

SUBROUTINE FINGERPRINT_CALC(N,V,DIST_MATRIX,TYP_I,TYP_J,NUMIONS,RMAX,SIGMA,DELTA,&
ORDER,FING,ATOM_FING)
IMPLICIT NONE
INTEGER, DIMENSION(:), INTENT(IN)               :: N, TYP_I, TYP_J, NUMIONS
REAL*8, INTENT(IN)                              :: RMAX, SIGMA, DELTA, V
REAL*8, DIMENSION(:,:), INTENT(IN)              :: DIST_MATRIX
REAL*8, DIMENSION(:), INTENT(OUT)               :: ORDER
REAL*8, DIMENSION(:,:), INTENT(OUT)             :: FING
REAL*8, DIMENSION(:,:,:), INTENT(OUT)             :: ATOM_FING
REAL*8, DIMENSION(803)                          :: ERF_TABLE
INTEGER                                         :: NUMBINS, BINS, I, J
REAL*8, DIMENSION(2)                              :: INTERVAL, MY_ERF
REAL*8                                          :: SIGM, R0, DELT,WEIGHT
REAL*8, PARAMETER                               :: PI = ACOS(-1.0D0)

DO I = 1, 803
    ERF_TABLE(I) = ERF((I - 402.0D0)/100.0D0);
ENDDO
NUMBINS = NINT(RMAX/DELTA)
FING = 0
FING(:,1) = -1.0d0
ORDER = 0
ATOM_FING = 0
SIGM = SIGMA/SQRT(2.0*LOG(2.0))
INTERVAL = 0

DO BINS = 2, NUMBINS
    DO I = 1, SUM(NUMIONS)
        DO J = 1, SIZE(DIST_MATRIX,1)
            IF (DIST_MATRIX(J,I) > 0 .AND.&
                ABS(DIST_MATRIX(J,I)-DELTA*(BINS-0.5D0)) < 4*SIGM) THEN

                R0 = DIST_MATRIX(J,I)
                INTERVAL(2) = SIGN(5.0D0,DELTA*BINS-R0)
                IF (ABS((DELTA*BINS-R0)/SQRT(2.0D0)) <= 5.0*SIGM) THEN
                    INTERVAL(2) = (DELTA*BINS-R0)/(SQRT(2.0D0)*SIGM)
                ENDIF
                INTERVAL(1) = SIGN(5.0D0,DELTA*(BINS-1)-R0)
                IF (ABS((DELTA*(BINS-1)-R0)/SQRT(2.0D0)) <= 5.0*SIGM) THEN
                    INTERVAL(1) = (DELTA*(BINS-1)-R0)/(SQRT(2.0D0)*SIGM)
                ENDIF
                MY_ERF = APPROX_ERF(INTERVAL,ERF_TABLE)
                DELT = 0.5*(MY_ERF(2)-MY_ERF(1))
                ATOM_FING(I,TYP_J(J),BINS) = ATOM_FING(I,TYP_J(J),BINS) +&
                                            DELT/(N(J)*R0**2)
                FING((TYP_I(I)-1)*SIZE(NUMIONS)+TYP_J(J),BINS) = &    
                    FING((TYP_I(I)-1)*SIZE(NUMIONS)+TYP_J(J),BINS) + DELT/R0**2
             ENDIF
        ENDDO
        ATOM_FING(I,:,BINS) = V*ATOM_FING(I,:,BINS)/(4.0D0*PI*DELTA) - 1.0D0
        DO J = 1, SIZE(NUMIONS)
            WEIGHT = NUMIONS(J)/REAL(SUM(NUMIONS))
            ORDER(I) = ORDER(I) + &
                     WEIGHT*DELTA*ATOM_FING(I,J,BINS)**2/(V/REAL(SUM(NUMIONS)))**(1.0/3.0)
        ENDDO
    ENDDO
    DO I = 1, SIZE(NUMIONS)
        DO J = 1, SIZE(NUMIONS)
            FING((I-1)*SIZE(NUMIONS)+J,BINS) = V*FING((I-1)*SIZE(NUMIONS)+J,BINS)&
                            /(4.0D0*PI*NUMIONS(I)*NUMIONS(J)*DELTA) - 1
        ENDDO
    ENDDO
ENDDO
ORDER = SQRT(ORDER)
END SUBROUTINE

FUNCTION DET (A) RESULT (ans)
IMPLICIT NONE
REAL*8, DIMENSION(:,:), INTENT(IN)  :: A
REAL*8  :: ans

ans =   A(1,1)*A(2,2)*A(3,3)  &
    - A(1,1)*A(2,3)*A(3,2)  &
    - A(1,2)*A(2,1)*A(3,3)  &
    + A(1,2)*A(2,3)*A(3,1)  &
    + A(1,3)*A(2,1)*A(3,2)  &
    - A(1,3)*A(2,2)*A(3,1)

END FUNCTION

FUNCTION ISZERO (N) RESULT (ANS)
IMPLICIT NONE
INTEGER, INTENT(IN)     :: N
INTEGER                 :: ANS

IF (N == 0) THEN
    ANS = 1
ELSE
    ANS = 0
ENDIF

END FUNCTION

FUNCTION APPROX_ERF(DATAPOINTS,ERF_TABLE) RESULT(MY_ERF)
IMPLICIT NONE
REAL*8, DIMENSION(:), INTENT(IN)        :: DATAPOINTS, ERF_TABLE
REAL*8, DIMENSION(SIZE(DATAPOINTS))     :: MY_ERF
INTEGER                                 :: I, X1
REAL*8                                  :: Y1, Y2, Y3, Y0, A0, A1, A2, A3, Q
DO I = 1, SIZE(DATAPOINTS)
    IF (DATAPOINTS(I) >= 4) THEN
        MY_ERF(I) = 1
    ELSEIF (DATAPOINTS(I) <= -4) THEN
        MY_ERF(I) = -1
    ELSE
        X1 = FLOOR(DATAPOINTS(I)/0.01D0+402.0D0)
        Y3 = ERF_TABLE(X1+2)
        Y2 = ERF_TABLE(X1+1)
        Y1 = ERF_TABLE(X1)
        Y0 = ERF_TABLE(X1-1)
        A0 = -Y0/6.0D0 + Y1/2.0D0 - Y2/2.0D0 + Y3/6.0D0
        A1 = Y0/2.0D0 - Y1 + Y2/2.0D0
        A2 = -Y0/3.0D0 - Y1/2.0D0 + Y2 - Y3/6.0D0
        A3 = Y1
        Q  = DATAPOINTS(I)/0.01D0 + 402.0D0 - X1
        MY_ERF(I) = ((A0*Q + A1)*Q + A2)*Q+A3
    ENDIF
ENDDO

END FUNCTION

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
double precision a(n,n), c(n,n)
double precision L(n,n), U(n,n), b(n), d(n), x(n)
double precision coeff
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

end program

