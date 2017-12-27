module polyhedra
    
    use atoms
    use cell
    use faces

    implicit none
   
    integer, parameter :: string_length = 11
    integer :: n_polyhedra = 0
 
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
        procedure(contains_atom_sub), deferred, nopass :: contains_atom
        procedure(assign_faces_sub), deferred :: assign_faces
        procedure(unique_id_func), deferred :: unique_id
    end type polyhedron

    abstract interface
        subroutine contains_atom_sub(this_atom)
            import :: atom
            class(atom) :: this_atom
        end subroutine contains_atom_sub
    
        subroutine assign_faces_sub( this )
            import :: polyhedron
            class( polyhedron ) :: this
        end subroutine assign_faces_sub

        integer function unique_id_func(this)
            import :: polyhedron
            class(polyhedron), intent(in) :: this
        end function unique_id_func
    end interface

    contains

    subroutine occupied_by( this, this_atom )
        implicit none
        class(polyhedron) :: this
        type(atom) :: this_atom
        integer :: dotsum, k, l
        double precision, dimension(3) :: testp, shiftvec
        
        do l=0, 7 ! loop over periodic images
            shiftvec(:) = mod(l,2)        * boxlen(1) * h(1,:) &
                        + mod(int(l/2),2) * boxlen(2) * h(2,:) &
                        + mod(int(l/4),2) * boxlen(3) * h(3,:)
            dotsum = 0
            do k=1, this%num_faces ! loop over each face of the polyhedron
                testp = (this_atom%r + shiftvec ) - this%face(k)%vertex(1)%r
                dotsum = dotsum - int( sign( dble(1.0), dot_product( testp, this%face(k)%normal ) ) )
            enddo
            if (dotsum == this%num_faces) then
                if (this%occupied) then
                    write(6,*) "warning! double occupation of ", this%string, this%unique_id()
                    write(6,*) "by ions ", this%occnum, this_atom%id
                    write(6,*) "at =>", this_atom%r
                endif
                this%occupied = .true.
                this%occnum = this_atom%id
                this_atom%polyid = this%id 
                call this%contains_atom( this_atom )
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

    pure subroutine set_vertex( this, v, this_atom )
        class(polyhedron), intent(inout) :: this
        integer, intent(in) :: v
        type(atom), intent(in) :: this_atom

        this%vertex(v) = this_atom
        this%vertex_ids(v) = this_atom%id
    end subroutine set_vertex

    subroutine set_vertices( this, atoms )
        class(polyhedron), intent(inout) :: this
        integer :: i
        type(atom), dimension(:), intent(in) :: atoms
    
        do i=1, this%num_vert
            this%vertex(i) = atoms(i)
            this%vertex_ids(i) = atoms(i)%id
        end do
        call this%set_centre
    end subroutine set_vertices

    pure subroutine set_centre( this )
        class(polyhedron), intent(inout) :: this
        integer :: l
        forall (l=1:3)
            this%centre(l) = sum( this%vertex(:)%r(l) ) / this%num_vert
        end forall
    end subroutine set_centre

    subroutine set_vertices_from_ids( this, atoms )
        class(polyhedron), intent(inout) :: this
        type(atom), dimension(:), intent(in) :: atoms
        type(atom), dimension(this%num_vert) :: temp_atoms
        integer :: i
        do i=1, this%num_vert
            temp_atoms(i) = atoms( this%vertex_ids( i ) )
        end do
        call this%set_vertices( temp_atoms )
    end subroutine set_vertices_from_ids

    pure subroutine enforce_pbc( this )
        class(polyhedron), intent(inout) :: this
        type point
            double precision, dimension(3) :: r
        end type point
        integer :: i, k
        type(point), dimension(this%num_vert) :: p
        double precision :: tetspread
        double precision, dimension(this%num_vert) :: dr

        forall (k=1:this%num_vert) p(k)%r = this%vertex(k)%r
        do i=1, 3
            forall (k=1:this%num_vert) dr(k) = relr( p(k)%r, i )
            tetspread = maxval(dr) - minval(dr)
            if (tetspread > halfboxlen(i)) then
                do k=1, this%num_vert
                    if ( relr(p(k)%r, i) < halfboxlen(i) ) then
                        p(k)%r = p(k)%r + h(i,:) * boxlen(i)
                    end if
                end do
            end if
        end do
        forall (k=1:this%num_vert) this%vertex(k)%r = p(k)%r
        call this%set_centre ! update centre point
    end subroutine enforce_pbc

end module polyhedra
