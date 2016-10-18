program defect_sort

implicit none

type ions
    double precision, dimension(3) :: r
    double precision, dimension(3) :: rref
    integer :: site
    integer :: siteold
    integer :: siteref
    integer :: sitedum
end type ions

type sites
    logical :: occupied
    integer :: pair
end type sites

type iofile
    character(len=20) :: filename
    integer :: iounit    
end type iofile

type(ions), allocatable, dimension(:) :: ion
type(sites), allocatable, dimension(:) :: site
character(len=20) :: fmtin
integer :: nbins, ncoords, nions, nfions, ntet, noct
double precision :: maxz
integer :: i, j, dum
integer, allocatable, dimension(:,:) :: zbin
integer :: thisbin, nocc1, nocc2, ndocc
character(len=1) :: zdir
integer :: zindex
integer :: nskip
logical :: dynflag ! only consider defects if present for sequential frames

type(iofile) :: inptfile, defectfile, posfile, pairfile, outfile1, outfile2, outfile3

inptfile%filename = "defect_sort.inpt"
outfile1%filename = "defectnum.dat"
outfile2%filename = "defectdist.dat"
outfile3%filename = "defectpos.dat"

open( file=outfile1%filename, newunit=outfile1%iounit )
open( file=outfile2%filename, newunit=outfile2%iounit )
open( file=outfile3%filename, newunit=outfile3%iounit )
open( file=inptfile%filename, status='old', newunit=inptfile%iounit )

read( inptfile%iounit, * ) posfile%filename
read( inptfile%iounit, * ) defectfile%filename
read( inptfile%iounit, * ) pairfile%filename
read( inptfile%iounit, * ) ncoords
read( inptfile%iounit, * ) nions
read( inptfile%iounit, * ) nfions
read( inptfile%iounit, * ) ntet
read( inptfile%iounit, * ) noct
read( inptfile%iounit, * ) maxz
read( inptfile%iounit, * ) nbins
read( inptfile%iounit, * ) zdir
read( inptfile%iounit, * ) dynflag
if (dynflag) then
    read( inptfile%iounit, *) nskip
    if ( nskip .eq. 0 ) dynflag = .false.
else
    nskip = 0
end if

close( inptfile%iounit )

allocate(ion(nions))
allocate(zbin(4,nbins))
allocate(site(-noct:ntet))

zbin = 0

select case (zdir)
    case ('x')
        zindex = 1
    case ('y')
        zindex = 2
    case ('z')
        zindex = 3
    case default
        write(6,*) "Not a valid string for close-packed stacking direction"
        stop
end select

open(10, file=pairfile%filename, status="old", form="formatted")
do i=1, ntet ! tetrahedra can be paired if we have hcp stacking
    read(10,*) dum ,j
    site(i)%pair = j
end do

open(11, file=defectfile%filename, status="old", form="formatted")
open(12, file=posfile%filename, status="old")

write(fmtin, '(A4,I4,A7)') "(I5,",nions,"(I5,X))" ! internal write to define output formatting

do i=1, nskip-1
    do j=1, nions
        read(12,*) ion(j)%r(1:3)
    end do
    do j=1, nfions
        read(12,*) dum, dum, dum
    end do
    read(11,*) dum, ion(:)%site
end do

if (dynflag) then
    write(6,*) "initial step"
    do j=1, nions
        read(12,*) ion(j)%rref(1:3)
    end do
    do j=1, nfions
        read(12,*) dum, dum, dum
    end do
    read(11,*) dum, ion(:)%site
    ion(:)%siteref = ion(:)%site
    ion(:)%siteold = ion(:)%site
end if

write(outfile2%iounit,*) "# step  noct  ntet1  ntet2  ndtet"

steploop: do i=nskip+1, ncoords
    do j=1, nions
        read(12,*) ion(j)%r(1:3) 
    end do
    do j=1, nfions
        read(12,*) dum, dum, dum
    end do
    
    if (dynflag) then
        if (mod(i,nskip) /= 0) then
            read(11,*) dum, ion(:)%sitedum
            cycle steploop
        end if
    end if

    read(11,*) dum, ion(:)%site

    write(6,*) "step ", i

    site(:)%occupied = .false.

    if (dynflag) then
        do j=1, nions
            if ( ion(j)%site == ion(j)%siteold ) then
                ion(j)%siteref = ion(j)%site
                ion(j)%rref = ion(j)%r
            end if
        end do
        ion(:)%siteold = ion(:)%site
    else
        ion(:)%siteref = ion(:)%site
        do j=1, nions
            ion(j)%rref = ion(j)%r
        end do
    end if

    do j=1, nions
        site(ion(j)%siteref)%occupied = .true.
    end do
    site(0)%occupied = .false.

    nocc1 = 0
    nocc2 = 0
    noct = 0
    ndocc = 0

    do j=1, nions
        thisbin = int((ion(j)%rref(zindex) / maxz) * nbins) + 1 
        if (ion(j)%siteref < 1) then
            zbin( 1, thisbin ) = zbin( 1, thisbin ) + 1
            noct = noct + 1
            write(outfile3%iounit,*) "oct", ion(j)%rref(:)
        else if (ion(j)%siteref <= ntet/2) then
            zbin(2,thisbin) = zbin(2,thisbin) + 1
            nocc1 = nocc1 + 1
            write(outfile3%iounit,*) "tet1", ion(j)%rref(:)
        else
            if (site(site(ion(j)%siteref)%pair)%occupied) then
                zbin(4,thisbin) = zbin(4,thisbin) + 1
                ndocc = ndocc + 1
                write(6,*) ion(j)%siteref
                write(6,*) site(ion(j)%siteref)%pair
                write(6,*) site(site(ion(j)%siteref)%pair)%occupied
                write(outfile3%iounit,*) "tetdd", ion(j)%rref(:)
            else
                zbin(3,thisbin) = zbin(3,thisbin) + 1
                nocc2 = nocc2 + 1
                write(outfile3%iounit,*) "tet2", ion(j)%rref(:)
            end if
        end if
    end do
    write(outfile1%iounit,*) i, noct, nocc1, nocc2, ndocc
end do steploop

do i=1, nbins
    write(outfile2%iounit,*) i, zbin(:,i)
end do
end program
