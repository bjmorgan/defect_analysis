module octahedra
    
    use polyhedra

    implicit none
    
    integer :: noct = 0
    integer :: id_offset
    integer, parameter :: noct_vertices = 6
    integer, parameter :: noct_faces = 8
    integer :: occoct
    character(len=string_length), parameter :: oct_string = 'octahedron'

    type, extends(polyhedron) :: octahedron
    contains
        procedure :: init => init_oct
        procedure :: assign_faces => oct_assign_faces
        procedure, nopass :: contains_atom => oct_contains_atom 
        procedure :: unique_id => oct_unique_id
    end type octahedron

    contains

    subroutine oct_contains_atom(this_atom)
        class(atom) :: this_atom
        this_atom%intet = .false.
        this_atom%inoct = .true.
        this_atom%insph = .false.
    end subroutine oct_contains_atom

    subroutine oct_assign_faces( this ) ! vertices are arranged in opposite pairs: (1,2)(3,4)(5,6)
        class(octahedron) :: this
        integer :: i, j, k, n
        do i=0, 1
            do j=0, 1
                do k=0, 1
                    n = 1 + i + j*2 + k*4
                    this%face(n)%vertex(1) = this%vertex(i+1)
                    this%face(n)%vertex(2) = this%vertex(j+3)
                    this%face(n)%vertex(3) = this%vertex(k+5)
                    call this%face(n)%define_normal( this%centre )
                enddo
            enddo
        enddo  
    end subroutine oct_assign_faces

    subroutine setup_oct( octa, noct )
        class(octahedron), dimension(:), allocatable :: octa
        integer, intent(in) :: noct
        integer :: i
        id_offset = n_sites
        allocate(octa(noct))
        do i=1, noct
            call octa(i)%init
            octa(i)%id = i
        enddo
        n_sites = n_sites + noct
    end subroutine setup_oct

    subroutine init_oct( this )
        class(octahedron) :: this
        this%num_vert = noct_vertices
        this%num_faces = noct_faces
        this%string = oct_string
        call this%alloc_vertices
        call this%alloc_faces
    end subroutine init_oct

    logical function oct_exists( vertex_ids, octs )
        integer, dimension(6), intent(in) :: vertex_ids
        class(octahedron), dimension(:), intent(in) :: octs
        integer i
        do i=1, size(octs)
            if (count(vertex_ids .eq. octs(i)%vertex_ids) == 6) then
                oct_exists = .true.
                return
            end if
        end do
        oct_exists = .false.
    end function oct_exists

    integer function oct_unique_id( this )
        class(octahedron), intent(in) :: this
        oct_unique_id = id_offset + this%id
    end function oct_unique_id

end module octahedra
