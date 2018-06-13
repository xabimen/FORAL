module fingerprint

public :: fp

contains

subroutine fp(coor,cell,ntype,ndif,d,rmax,fprint)
implicit none
real, dimension(:,:), intent(in)        :: coor, cell
integer, dimension(:), intent(in)       :: ntype
integer, intent(in)                     :: d, ndif
real, intent(in)                        :: rmax
real, dimension(:,0:,:), intent(out)    :: fprint
integer, dimension(:), allocatable      :: n
real, dimension(:,:), allocatable       :: dist_matrix
integer, dimension(:,:), allocatable    :: dist_atoms
real                                    :: start, w, r, v
real, dimension(3)                      :: r_vec
integer                                 :: i, j, k
integer, dimension(:), allocatable      :: natoms, atomtype
real, dimension(:), allocatable         :: a
real, parameter                         :: pi = acos(-1.0)

allocate(natoms(ndif), atomtype(ndif),a(d))

atomtype(1) = 14
atomtype(2) = 8
natoms = 0
do i = 1, size(ntype)
    if (ntype(i) == 14) then
        natoms(1) = natoms(1) + 1
    elseif (ntype(i) == 8) then
        natoms(2) = natoms(2) + 1
    else
        print*, "Atomoren bat ez da silizioa edo oxigenoa"
    endif
enddo

call makeMatrices(cell,coor,natoms,atomtype,rmax,n,v,dist_matrix,dist_atoms)

start = minval(sqrt(sum(dist_matrix(2:,:)**2,2)))/2.0
w = (Rmax - start)/real(d-1,8)
do k = 1, d
       a(k) = (Rmax - start)*real(k-1)/real(d-1,4) + start
enddo

fprint = 0.0d0

do j = 1, size(ntype)
    fprint(j,0,:) = ntype(j)
enddo

do i = 1, size(dist_atoms,1)
    r_vec = dist_matrix(i,:)
    r = sqrt(sum(r_vec**2))
    j = dist_atoms(i,2)
    if (N(i) == 14) then
        do k = 1, d
            fprint(j,k,:) = fprint(j,k,:) + &
            (r_vec/r)*exp(-(r-a(k))**2/(2.0d0*w**2))/(sqrt(2.0d0*pi)*w)*(cos(pi*r/Rmax)+1.0d0)/2.0d0
        enddo
    elseif (N(i) == 8) then
        do k = 1, d
            fprint(j,d+k,:) = fprint(j,d+k,:) + &
            (r_vec/r)*exp(-(r-a(k))**2/(2.0d0*w**2))/(sqrt(2.0d0*pi)*w)*(cos(pi*r/Rmax)+1.0d0)/2.0d0
        enddo
    endif

enddo

end subroutine

subroutine makeMatrices(cell,coordinates,numIons,atomType,Rmax,N,V,dist_matrix,dist_atoms)

implicit none

real, dimension(:,:), intent(in)              :: cell, coordinates
integer, dimension(:), intent(in)               :: atomtype, numIons
real, intent(in)                              :: Rmax
integer, dimension(:), allocatable, intent(out) :: N
real, intent(out)                             :: V
real, dimension(:,:), allocatable, intent(out):: dist_matrix
integer                                         :: species, i, j, c, lengthX,&
                                                   lengthY, lengthZ,&
                                                   quit_marker_x, quit_marker_y, &
                                                   quit_marker_z, k, quadrants,&
                                                   current_cell, natoms,&
                                                   marker, basic_cell, Nfull, &
                                                   c2, c3
integer, dimension(sum(numIons))                :: atomType1
real, dimension(sum(numIons))                 :: distances
real, dimension(13,3)                         :: vect
real, dimension(13)                           :: abs_vect
integer, dimension(8,3)                         :: condition, signum
real                                          :: x, y, z, Rij
integer, dimension(100000)                      :: N2
integer, dimension(100000,2)                    :: dist_atoms2
integer, dimension(:,:), allocatable, intent(out)   :: dist_atoms
real, dimension(100000,3)                     :: dist_matrix2

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

FUNCTION APPROX_ERF(DATAPOINTS,ERF_TABLE) RESULT(MY_ERF)
IMPLICIT NONE
REAL, DIMENSION(:), INTENT(IN)        :: DATAPOINTS, ERF_TABLE
REAL, DIMENSION(SIZE(DATAPOINTS))     :: MY_ERF
INTEGER                                 :: I, X1
REAL                                  :: Y1, Y2, Y3, Y0, A0, A1, A2, A3, Q
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

FUNCTION DET (A) RESULT (ans)
IMPLICIT NONE
REAL, DIMENSION(:,:), INTENT(IN)  :: A
REAL  :: ans

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
end module
