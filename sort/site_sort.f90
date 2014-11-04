program vacancy_sort

implicit none

type thision
    double precision, dimension(3) :: r
    integer :: tet
end type thision

type, abstract :: polyhedron
    double precision, dimension(3) :: centre
    integer :: occupied
    integer :: occold
    integer :: occref
    integer :: occdum
end type polyhedron

type, extends(polyhedron) :: tetrahedron
    type(thision), dimension(4) :: ion 
    integer, dimension(4) :: num
end type tetrahedron

type iofile
    character(len=20) :: filename
    integer :: iounit
end type iofile

type(thision), allocatable, dimension(:) :: ion
type(tetrahedron), allocatable, dimension(:) :: tetra
type(iofile) :: inptfile, tetlistf, tetoccf, posfile, outfile1, outfile2, outfile3
character(len=20) :: fmtin
double precision, dimension(3) :: boxlen, hboxlen
double precision, dimension(3) :: r1, r2, r3, r4
double precision :: tetspread, vaccent
double precision :: maxz
character(len=1) :: zdir
integer :: nbins, ntet, nions, nfions, ncoords
integer :: i, j, k, l, dum
integer, allocatable, dimension(:,:) :: zbin
integer :: thisbin, kstart
integer :: zindex
integer :: nskip
integer :: nvac1, nvac2, nocc1, nocc2
logical :: vaclog
logical :: dynflag ! only consider defects if present for sequential frames

inptfile%filename = "vacancy_sort.inpt"
outfile1%filename = "vacancynum.dat"
outfile2%filename = "vacancydist.dat"
outfile3%filename = "vaccentpos.dat"

open( file=outfile1%filename, newunit=outfile1%iounit )
open( file=outfile2%filename, newunit=outfile2%iounit )
open( file=outfile3%filename, newunit=outfile3%iounit )
open( file=inptfile%filename, status="old", newunit=inptfile%iounit )

read( inptfile%iounit, * ) posfile%filename
read( inptfile%iounit, * ) ncoords
read( inptfile%iounit, * ) nions
read( inptfile%iounit, * ) nfions
read( inptfile%iounit, * ) tetlistf%filename
read( inptfile%iounit, * ) ntet
read( inptfile%iounit, * ) maxz
read( inptfile%iounit, * ) nbins
read( inptfile%iounit, * ) tetoccf%filename
read( inptfile%iounit, * ) boxlen(1:3)
read( inptfile%iounit, * ) zdir
read( inptfile%iounit, * ) dynflag
if (dynflag) then
    read(inptfile%iounit,*) nskip
else
    nskip = 0
endif
read(inptfile%iounit,*) vaclog
if (vaclog) then
    read(inptfile%iounit,*) vaccent
endif
close(inptfile%iounit)

select case (zdir)
    case ("x")
        zindex = 1
    case ("y")
        zindex = 2
    case ("z")
        zindex = 3
    case default
        write(6,*) "Not a valid string for close-packed stacking direction:", zdir
        stop
end select

allocate(ion(nfions))
allocate(zbin(5,nbins))
allocate(tetra(ntet))

hboxlen = boxlen/2.0

zbin = 0

open(file=tetlistf%filename, status="old", newunit=tetlistf%iounit)
open(file=posfile%filename, status="old", newunit=posfile%iounit)
open(file=tetoccf%filename, status="old", newunit=tetoccf%iounit)

write(fmtin, '(A4,I4,A7)') "(I5,",nions,"(I5,X))" ! internal write to define output formatting

do i=1, ntet
    read(tetlistf%iounit,*) tetra(i)%num(1:4)
enddo

do i=1, nskip-1
    call read_positions( posfile, nions, ion )
    read(tetoccf%iounit,*)dum, tetra(:)%occupied
enddo

if (dynflag) then
    write(6,*) "initial step"
    call read_positions( posfile, nions, ion )
    read(tetoccf%iounit,*) dum, tetra%occupied
    tetra%occref = tetra%occupied
    tetra%occold = tetra%occupied
endif

write(outfile2%iounit,*) "#    step    nvac1    nocc1    nvac2    nocc2"

steploop: do i=nskip+1, ncoords
    call read_positions( posfile, nions, ion )
    
    if (dynflag) then
        if (mod(i,nskip) /= 0) then
            read(tetoccf%iounit,*) dum, tetra%occdum
            cycle steploop
        endif
    endif
 
    nvac1 = 0
    nocc1 = 0
    nvac2 = 0
    nocc2 = 0
  
    read(tetoccf%iounit,*) dum, tetra%occupied

    write(6,*) "step ", i

    if (dynflag) then
        do j=1, ntet
            if (tetra(j)%occupied == tetra(j)%occold) then
                tetra(j)%occref = tetra(j)%occupied
             endif
        enddo
        tetra%occold = tetra%occupied
    else
        tetra%occref = tetra%occupied
    endif

! allocate ion positions to lattice sites

    if (vaclog) then ! only apply periodic boundaries along xz
        kstart = 2
    else
        kstart = 1
    endif
   
    do j=1, ntet
        r1 = ion(tetra(j)%num(1))%r
        r2 = ion(tetra(j)%num(2))%r
        r3 = ion(tetra(j)%num(3))%r
        r4 = ion(tetra(j)%num(4))%r

! enforce periodic boundary conditions
! only for orthorhombic cells at the moment
        if (vaclog) then
            if (r1(1) > vaccent)then
                r1(1) = r1(1) - boxlen(1)
            endif

            if (r2(1) > vaccent)then
                r2(1) = r2(1) - boxlen(1)
            endif

            if (r3(1) > vaccent)then
                r3(1) = r3(1) - boxlen(1)
            endif

            if (r4(1) > vaccent)then
                r4(1) = r4(1) - boxlen(1)
            endif
        endif

        do k=kstart, 3
            tetspread = max(r1(k), r2(k), r3(k), r4(k)) - min(r1(k), r2(k), r3(k), r4(k))
            if (tetspread > hboxlen(k)) then
                if (r1(k) < hboxlen(k)) then 
                    r1(k) = r1(k) + boxlen(k)
                endif
                if (r2(k) < hboxlen(k)) then 
                    r2(k) = r2(k) + boxlen(k)
                endif
                if (r3(k) < hboxlen(k)) then 
                    r3(k) = r3(k) + boxlen(k)
                endif
                if (r4(k) < hboxlen(k)) then 
                    r4(k) = r4(k) + boxlen(k)
                endif
            endif
        enddo
        tetra(j)%ion(1)%r = r1
        tetra(j)%ion(2)%r = r2
        tetra(j)%ion(3)%r = r3
        tetra(j)%ion(4)%r = r4
        
        tetra(j)%centre = (r1+r2+r3+r4)/4
    
        thisbin = int((tetra(j)%centre(zindex) / maxz) * nbins)+1 
        if (thisbin > nbins) then
            thisbin = thisbin - nbins
        endif
        zbin(1,thisbin) = zbin(1,thisbin) + 1
        if (j <= ntet/2) then
            if (tetra(j)%occref == 0) then
                zbin(2,thisbin) = zbin(2,thisbin) + 1
                nvac1 = nvac1 + 1
                write( outfile3%iounit, * ) "vac1", tetra(j)%centre
             else
                zbin(3,thisbin) = zbin(3,thisbin) + 1
                nocc1 = nocc1 + 1
                write( outfile3%iounit, * ) "occ1", tetra(j)%centre
            endif
        else
            if (tetra(j)%occref == 0) then
                zbin(4,thisbin) = zbin(4,thisbin) + 1
                nvac2 = nvac2 + 1
                write( outfile3%iounit, * ) "vac2", tetra(j)%centre
            else
                zbin(5,thisbin) = zbin(5,thisbin) + 1
                nocc2 = nocc2 + 1
                write( outfile3%iounit, * ) "occ2", tetra(j)%centre
            endif
        endif
    enddo
    write(outfile1%iounit,*) i, nvac1, nocc1, nvac2, nocc2
enddo steploop

do i=1, nbins
    write(outfile2%iounit,*) i,zbin(:,i)
enddo

stop

contains

subroutine read_positions( posfile, nskip, ions )
    type(iofile), intent(in) :: posfile
    type(thision), dimension(:), intent(inout) :: ions
    integer, intent(in) :: nskip
    integer :: i
    do i=1, nskip
        read(posfile%iounit, *)
    end do
    do i=1, size(ions)
        read(posfile%iounit, *) ions(i)%r
    end do
end subroutine read_positions

end program
