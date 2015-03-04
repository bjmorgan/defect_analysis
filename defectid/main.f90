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

character(len=30) :: posfile, cellfile, inptfile, tetfile, octfile, poly_out, atoms_out, npolyocc
integer :: i, j, polyswitch
integer :: nconfigs, nstep, thispoly, nspec, lattice_spec, mobile_spec
integer :: natomsout ! number of mobile atoms 
integer, allocatable, dimension(:) :: nsp, polylist, sitelist
character(len=30) :: fmtout, fmtout2
integer :: fin, fout1, fout2, fout3, fcell

interface

    function diagonal( square_matrix )
        double precision, dimension(:,:), intent(in) :: square_matrix
        double precision, dimension( size(square_matrix, 1) ) :: diagonal
    end function diagonal

end interface

inptfile  = 'defect_new.inpt'
poly_out  = 'sites_atoms.dat'
atoms_out = 'atoms_sites.dat'
npolyocc  = 'npolyocc.dat'

open( file=inptfile, status='old', newunit=fin)
read( fin, * ) posfile
read( fin, * ) cellfile
read( fin, * ) tetfile
read( fin, * ) octfile
read( fin, * ) nconfigs
read( fin, * ) nspec
allocate (nsp(nspec))
allocate (spec(nspec))
do i=1, nspec
    read(fin,*) nsp(i)
    allocate (spec(i)%ion(nsp(i)))
    forall (j=1:nsp(i)) spec(i)%ion(j)%id = j
enddo
read( fin, * ) lattice_spec
read( fin, * ) mobile_spec
read( fin, * ) ntet
read( fin, * ) noct
read( fin, * ) boxlen(:)
read( fin, * ) h(:,:)
close(fin)

call setup_tet(tetra, ntet)
call setup_oct(octa, noct)

natomsout = nsp( mobile_spec ) 

allocate( polylist( natomsout ) )
allocate( sitelist( ntet+noct ) )
 
! read in ids for tetrahedra vertices
open(file=tetfile, status='old', newunit=fin)
do i=1, ntet 
    read(fin,*) tetra(i)%vertex_ids
enddo
close(fin)

! read in ids for octahedra vertices
open(file=octfile, status='old', newunit=fin)
do i=1, noct
    read(fin,*) octa(i)%vertex_ids
enddo
close(fin)

open( file = posfile,   status='old', newunit=fin )
open( file = cellfile,  status='old', newunit=fcell )
open( file = npolyocc,  newunit=fout1 )
open( file = atoms_out, newunit=fout2 )
open( file = poly_out,  newunit=fout3 )

do concurrent (i=1:nspec)
    spec(i)%ion%polyid = 0
    spec(i)%ion%inoct = .false.
    spec(i)%ion%intet = .false.
    spec(i)%ion%prev_intet = .false.
    spec(i)%ion%prev_inoct = .false.
end do

do nstep=1, nconfigs

    tetra%occnum = 0
    tetra%occupied = .false.
    octa%occnum = 0
    octa%occupied = .false.

    do j=1, nsp( mobile_spec )
        spec(mobile_spec)%ion(j)%previous_polyid = spec(mobile_spec)%ion(j)%polyid
        spec(mobile_spec)%ion(j)%prev_inoct = spec(mobile_spec)%ion(j)%inoct
        spec(mobile_spec)%ion(j)%prev_intet = spec(mobile_spec)%ion(j)%intet
        spec(mobile_spec)%ion(j)%polyid = 0
        spec(mobile_spec)%ion(j)%inoct  = .false.
        spec(mobile_spec)%ion(j)%intet  = .false.
    end do

    read( fcell, * ) h(:,:) ! cell unit vectors are rows in h
    read( fcell, * ) boxlen(:) 
    halfboxlen = boxlen / 2.0
    cboxlen = diagonal( h ) * boxlen
    halfcboxlen = cboxlen / 2.0
    
    do i=1, nspec
        do j=1, nsp(i)
            associate( this_ion => spec(i)%ion(j) )
                read(fin,*) this_ion%r(1:3) 
!                 if( spec(i) == mobile_spec .and. j==329 ) write(6,*) this_ion%r        
                this_ion%r = move_inside_cell( this_ion%r ) ! apply periodic boundary conditions
!                 if( spec(i) == mobile_spec .and. j==329 ) write(6,*) this_ion%r 
            end associate
        enddo
    enddo
  
    do thispoly = 1, ntet ! construct tetrahedra
        call tetra(thispoly)%set_vertices_from_ids( spec( lattice_spec )%ion )
        call tetra(thispoly)%enforce_pbc
        call tetra(thispoly)%assign_faces
    end do

    do thispoly = 1, noct ! construct octahedra
        call octa(thispoly)%set_vertices_from_ids( spec( lattice_spec )%ion )
        call octa(thispoly)%enforce_pbc
        call octa(thispoly)%assign_faces
    end do

    ionloop: do j=1, nsp( mobile_spec )
        ! test whether ions are in the same polyhedra as the previous step
        if ( spec( mobile_spec )%ion(j)%prev_intet ) then
            call tetra( spec( mobile_spec )%ion(j)%previous_polyid )%occupied_by( spec( mobile_spec )%ion(j) )
        else if ( spec( mobile_spec )%ion(j)%prev_inoct ) then
            call  octa( spec( mobile_spec )%ion(j)%previous_polyid )%occupied_by( spec( mobile_spec )%ion(j) )
        end if
        if ( spec(mobile_spec)%ion(j)%intet .or. spec(mobile_spec)%ion(j)%inoct ) cycle ionloop
        ! ion has moved. search over remaining tetrahedra
        tetloop: do thispoly = 1, ntet
            if ( spec( mobile_spec )%ion(j)%prev_intet .and. thispoly .eq. spec( mobile_spec)%ion(j)%previous_polyid ) cycle tetloop
            call tetra( thispoly )%occupied_by( spec(mobile_spec)%ion(j) )
            if ( spec(mobile_spec)%ion(j)%intet ) cycle ionloop
        end do tetloop
        ! search over remaining octahedra
        octloop: do thispoly = 1, noct
            if ( spec( mobile_spec )%ion(j)%prev_inoct .and. thispoly .eq. spec( mobile_spec)%ion(j)%previous_polyid ) cycle octloop
            call octa( thispoly )%occupied_by( spec(mobile_spec)%ion(j) )
            if ( spec(mobile_spec)%ion(j)%inoct ) cycle ionloop
        end do octloop
        ! if we reach here without cycling ionloop, this ion has not been located in any polyhedron
        write(6,*) 'Ion ', j,' not in any polyhedra, at ', spec(mobile_spec)%ion(j)%r
        stop
    end do ionloop

    write(6,*) "step ", nstep, count(tetra%occupied), count(octa%occupied)
    write(fout1,*) nstep, count(tetra%occupied), count(octa%occupied)

    natomsout = 0
    polylist = 0
    do j=1, nsp( mobile_spec )
        associate( this_ion => spec(mobile_spec)%ion(j) )
            natomsout = natomsout + 1
            if ( this_ion%inoct ) then
                polyswitch = -1
            else
                polyswitch = +1
            endif
            polylist(natomsout) = this_ion%polyid * polyswitch
        end associate
    enddo

    sitelist( 1:ntet ) = tetra%occnum
    sitelist( ntet+1:ntet+noct ) = octa%occnum

    write( fmtout,  '(A4,I5,A7)') "(I5,", natomsout, "(I6,X))" ! internal write to define output formatting
    write( fmtout2, '(A4,I5,A7)') "(I5,", ntet+noct, "(I6,X))" ! internal write to define output formatting
    write( fout2, fmtout )  nstep, polylist(:)
    write( fout3, fmtout2 ) nstep, sitelist(:)

enddo !ends loop over nconfig steps

close( fin )
close( fout1 )
close( fout2 )
close( fout3 )

stop

end program defect_new

function diagonal( square_matrix )
    implicit none
    double precision, dimension(:,:), intent(in) :: square_matrix
    double precision, dimension(size(square_matrix)) :: temp_array
    double precision, dimension(size(square_matrix, 1)) :: diagonal
    temp_array = pack( square_matrix, .true. )
    diagonal = temp_array( 1::size(diagonal) + 1 )
end function diagonal
