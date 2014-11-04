program defect_new

use atoms
use tetrahedra
use octahedra
use cell
 
implicit none

type species
    type (atom), dimension(:), allocatable :: ion
end type species  

type (species), dimension(:), allocatable :: spec
class (tetrahedron), dimension(:), allocatable :: tetra
class (octahedron), dimension(:), allocatable :: octa

character(len=20) :: posfile, inptfile, tetfile, octfile, poly_out, atoms_out, npolyocc
integer :: i, j, k, l, m, polyswitch
integer :: nconfigs, nstep, thispoly, nspec, lattice_spec
integer :: natomsout ! number of mobile atoms 
integer, allocatable, dimension(:) :: nsp, polylist, sitelist
character(len=20) :: fmtout, fmtout2
integer :: fin, fout1, fout2, fout3

interface

    function diagonal( square_matrix )
        double precision, dimension(:,:), intent(in) :: square_matrix
        double precision, dimension( size(square_matrix, 1) ) :: diagonal
    end function diagonal

end interface

inptfile = 'defect_new.inpt'
poly_out = 'sites_atoms.dat'
atoms_out = 'atoms_sites.dat'
npolyocc = 'npolyocc.dat'

open(file=inptfile, status='old', newunit=fin)
read(fin,*) posfile
read(fin,*) tetfile
read(fin,*) octfile
read(fin,*) nconfigs
read(fin,*) nspec
allocate (nsp(nspec))
allocate (spec(nspec))
do i=1, nspec
    read(fin,*) nsp(i)
    allocate (spec(i)%ion(nsp(i)))
    forall (j=1:nsp(i)) spec(i)%ion(j)%id = j
enddo
read(fin,*) lattice_spec
read(fin,*) ntet
read(fin,*) noct
read(fin,*) cboxlen(:)
read(fin,*) h(:,:)
close(fin)

!generate cell lengths for orthorhombic cell from boxlen(:)
boxlen = cboxlen*diagonal(h)

halfcboxlen = cboxlen/2.0
halfboxlen = boxlen/2.0

call setup_tet(tetra, ntet)
call setup_oct(octa, noct)

natomsout = sum(nsp) - nsp(lattice_spec) ! mobile atoms are assumed to be those not defining lattice sites

allocate(polylist(natomsout))
allocate(sitelist(ntet+noct))
 
! read in ids for tetrahedra vertices
open(file=tetfile, status='old', newunit=fin)
do i=1, ntet 
    read(fin,*)  tetra(i)%vertex_ids
enddo
close(fin)

! read in ids for octahedra vertices
open(file=octfile, status='old', newunit=fin)
do i=1, noct
    read(fin,*)  octa(i)%vertex_ids
enddo
close(fin)

open(file=posfile, status='old', newunit=fin)
open(file=npolyocc, newunit=fout1)
open(file=atoms_out, newunit=fout2)
open(file=poly_out, newunit=fout3)

do nstep=1, nconfigs

    tetra%occnum = 0
    tetra%occupied = .false.
    octa%occnum = 0
    octa%occupied = .false.

    forall (i=1:nspec) 
        spec(i)%ion%polyid = 0
        spec(i)%ion%inoct = .false.
        spec(i)%ion%intet = .false.
    end forall

    do i=1, nspec
        do j=1, nsp(i)
            read(fin,*) spec(i)%ion(j)%r(1:3)
        enddo
    enddo

    do thispoly=1, ntet
        associate( tet => tetra(thispoly) )
            call tet%set_vertices_from_ids(spec(lattice_spec)%ion)
            call tet%enforce_pbc
            call tet%assign_faces
            do i=1, nspec
                if (i /= lattice_spec) then
                    do j=1, nsp(i) 
                        call tet%occupied_by(spec(i)%ion(j))
                    enddo
                endif
            enddo
        end associate 
    enddo

    do thispoly=1, noct
        associate( oct => octa(thispoly) )
            call oct%set_vertices_from_ids(spec(lattice_spec)%ion)
            call oct%enforce_pbc
            call oct%assign_faces
            do i=1, nspec
                if(i /= lattice_spec) then
                    do j=1, nsp(i)
                        call oct%occupied_by(spec(i)%ion(j))
                    enddo
                endif
            enddo
        end associate
    enddo

    write(6,*) "step ", nstep, count(tetra%occupied), count(octa%occupied)
    write(fout1,*) nstep, count(tetra%occupied), count(octa%occupied)

    natomsout = 0
    polylist = 0
    do i=1, nspec
        if (i /= lattice_spec) then
            do j=1, nsp(i)
                natomsout = natomsout + 1
                if( spec(i)%ion(j)%inoct )then
                    polyswitch = -1
                else
                    polyswitch = +1
                endif
                polylist(natomsout) = spec(i)%ion(j)%polyid * polyswitch
            enddo
        endif
    enddo

    sitelist(1:ntet) = tetra%occnum
    sitelist(ntet+1:ntet+noct) = octa%occnum

    write(fmtout, '(A4,I4,A7)') "(I5,",natomsout,"(I5,X))" ! internal write to define output formatting
    write(fmtout2, '(A4,I4,A7)') "(I5,",ntet+noct,"(I5,X))" ! internal write to define output formatting
    write(fout2, fmtout) nstep, polylist(:)
    write(fout3, fmtout2) nstep, sitelist(:)

enddo !ends loop over nconfig steps

close(fin)
close(fout1)
close(fout2)
close(fout3)

stop

end program defect_new

function diagonal( square_matrix )
    implicit none
    double precision, dimension(:,:), intent(in) :: square_matrix
    double precision, dimension(size(square_matrix)) :: temp_array
    double precision, dimension(size(square_matrix, 1)) :: diagonal
    temp_array = pack(square_matrix, .true.)
    diagonal = temp_array(1::size(diagonal)+1)
end function diagonal

