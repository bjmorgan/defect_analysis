module tetrahedra
    
!    use atoms
!    use cell
    use polyhedra

    implicit none
    
    integer :: ntet = 0
    integer :: id_offset
    integer, parameter :: ntet_vertices = 4
    integer, parameter :: ntet_faces = 4
    character(len=string_length), parameter :: tet_string = 'tetrahedron'

    type, extends(polyhedron) :: tetrahedron
        integer :: orientation
    contains
        procedure :: init => init_tet
        procedure :: assign_faces => tet_assign_faces
        procedure, nopass :: contains_atom => tet_contains_atom
        procedure :: unique_id => tet_unique_id
    end type tetrahedron

    contains

    subroutine tet_contains_atom(this_atom)
        class(atom) :: this_atom
        this_atom%intet = .true.
        this_atom%inoct = .false.
    end subroutine tet_contains_atom

    subroutine tet_assign_faces( this )
        class(tetrahedron) :: this
        integer :: i, j
        do i=1, 4 
            forall(j=1:3) this%face(i)%vertex(j) = this%vertex(mod(i+j,4)+1)
            call this%face(i)%define_normal( this%centre )
        end do
    end subroutine tet_assign_faces

    subroutine setup_tet( tetra, ntet )
        class(tetrahedron), dimension(:), allocatable :: tetra
        integer, intent(in) :: ntet
        integer :: i
        id_offset = n_sites
        allocate(tetra(ntet))
        do i=1, ntet
            call tetra(i)%init
            tetra(i)%id = i
        enddo 
        n_sites = n_sites + ntet
        end subroutine setup_tet

     subroutine init_tet( this )
        class(tetrahedron) :: this
        this%num_vert = ntet_vertices
        this%num_faces = ntet_faces
        this%string = tet_string
        call this%alloc_vertices
        call this%alloc_faces
    end subroutine init_tet

    integer function tet_unique_id( this )
        class(tetrahedron), intent(in) :: this
        tet_unique_id = id_offset + this%id
    end function tet_unique_id

end module tetrahedra
