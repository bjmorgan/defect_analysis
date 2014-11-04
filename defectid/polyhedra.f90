module polyhedra
    
    use atoms
    use cell
    use faces

    implicit none
   
    integer, parameter :: string_length = 11
 
    type, abstract :: polyhedron
        double precision, dimension(3) :: centre
        type(atom), dimension(:), allocatable :: vertex
        integer, dimension(:), allocatable :: vertex_ids
        integer :: num_vert
        integer :: num_faces
        type(face), dimension(:), allocatable :: face
        logical :: occupied 
        integer :: occnum
        type(atom) :: occ_atom
        integer :: id
        character(len=string_length) :: string
    contains
        procedure :: set_vertex
        procedure :: set_vertices
        procedure :: set_vertices_from_ids
        procedure :: alloc_vertices
        procedure :: alloc_faces
        procedure :: enforce_pbc
        procedure :: occupied_by
        procedure, private :: set_centre
        procedure(contains_atom_sub), deferred :: contains_atom
        procedure(assign_faces_sub), deferred :: assign_faces
    end type polyhedron

    abstract interface
        subroutine contains_atom_sub(this, this_atom)
            import :: polyhedron, atom
            class(polyhedron) :: this
            class(atom) :: this_atom
        end subroutine contains_atom_sub
    
        subroutine assign_faces_sub( this )
            import :: polyhedron
            class( polyhedron ) :: this
        end subroutine assign_faces_sub
    end interface

    contains

    subroutine occupied_by( this, this_atom )
        implicit none
        class(polyhedron) :: this
        type(atom) :: this_atom
        integer :: dotsum, k, l
        double precision, dimension(3) :: testp, shiftvec
        double precision :: dotprod
        do l=0, 7 ! loop over periodic images
            shiftvec(:) = mod(l,2)        * cboxlen(1) * h(1,:) &
                        + mod(int(l/2),2) * cboxlen(2) * h(2,:) &
                        + mod(int(l/4),2) * cboxlen(3) * h(3,:)
            dotsum = 0
            do k=1, this%num_faces ! loop over each face of the polyhedron
                testp = (this_atom%r + shiftvec ) - this%face(k)%vertex(1)%r
                dotprod = sum(testp * this%face(k)%normal)
                dotsum = dotsum - sign(1.0,dotprod)
            enddo
            if (dotsum == this%num_faces) then
                if (this%occupied) then
                    write(6,*) "warning! double occupation of ", this%string, this%id
                    write(6,*) "by ions ", this%occnum, this_atom%id
                endif
                this%occupied = .true.
                this%occnum = this_atom%id
                this_atom%polyid = this%id ! this is going to cause problems with identically numbered tetrahedra and octahedra
                call this%contains_atom(this_atom)
                return
            endif
        enddo
    end subroutine occupied_by

    subroutine alloc_vertices( this )
        class(polyhedron) :: this
        allocate( this%vertex(this%num_vert) )
        allocate( this%vertex_ids(this%num_vert) )
    end subroutine alloc_vertices

    subroutine alloc_faces( this )
        class(polyhedron) :: this
        allocate( this%face(this%num_faces) )
    end subroutine alloc_faces

    subroutine set_vertex( this, v, this_atom )
        class(polyhedron) :: this
        integer :: v
        type(atom) :: this_atom

        this%vertex(v) = this_atom
        this%vertex_ids(v) = this_atom%id
    end subroutine set_vertex

    subroutine set_vertices( this, atoms )
        class(polyhedron) :: this
        integer :: i
        type(atom), dimension(:) :: atoms
    
        do i=1, this%num_vert
            this%vertex(i) = atoms(i)
            this%vertex_ids(i) = atoms(i)%id
        end do
        call this%set_centre
    end subroutine set_vertices

    subroutine set_centre( this )
        class(polyhedron) :: this
        integer :: l
        forall (l=1:3)
            this%centre(l) = sum( this%vertex(:)%r(l) ) / this%num_vert
        end forall
    end subroutine set_centre

    subroutine set_vertices_from_ids( this, atoms )
        class(polyhedron) :: this
        type(atom), dimension(:), intent(in) :: atoms
        type(atom), dimension(this%num_vert) :: temp_atoms
        integer :: i
        forall (i=1:this%num_vert) temp_atoms(i) = atoms(this%vertex_ids(i))
        call this%set_vertices( temp_atoms )
    end subroutine set_vertices_from_ids

    subroutine enforce_pbc( this )
        class( polyhedron ) :: this
        double precision, dimension(3) :: minimum_r
        integer :: i

        do i=1, 3
            minimum_r(i) = minval( this%vertex(:)%r(i) )
        end do
        do i=1, this%num_vert
            this%vertex(i)%r = minimum_r - r_as_minimum_image( dr( minimum_r, this%vertex(i)%r ) )
        end do
        call this%set_centre ! update centre point
    end subroutine enforce_pbc

end module polyhedra
