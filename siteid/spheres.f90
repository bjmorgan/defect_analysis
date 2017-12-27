module spheres

    use atoms
    use cell
    use sites

    implicit none

    integer :: nspheres = 0
    integer :: id_offset
    character(len=string_length), parameter :: sph_string = 'sphere'

    type, extends(site) :: spherical_site
        double precision :: cutoff
    contains
        procedure :: occupied_by => sphere_occupied_by
        procedure :: init => init_sphere
        procedure, nopass :: contains_atom => sph_contains_atom
        procedure :: unique_id => sph_unique_id
    end type spherical_site

    contains

    subroutine sph_contains_atom(this_atom)
        class(atom) :: this_atom
        this_atom%intet = .false.
        this_atom%inoct = .false.
        this_atom%insph = .true.
    end subroutine sph_contains_atom

    subroutine init_sphere( this )
        class(spherical_site) :: this
        this%string = sph_string
    end subroutine init_sphere

    integer function sph_unique_id( this )
        class(spherical_site), intent(in) :: this
        sph_unique_id = id_offset + this%id
    end function sph_unique_id

    subroutine sphere_occupied_by( this, this_atom )
        implicit none
        class(spherical_site) :: this
        type(atom) :: this_atom
        double precision, dimension(3) :: dr
        dr = r_as_minimum_image( this_atom%r - this%centre )
        if ( dot_product( dr, dr ) < this%cutoff**2 ) then ! TODO should really store cutoffsq
            this%occupied = .true.
            this%occnum = this_atom%id
            this_atom%polyid = this%id
            call this%contains_atom( this_atom )
            return
        endif
    end subroutine sphere_occupied_by

    subroutine setup_sph( spheres, nspheres )
        class(spherical_site), dimension(:), allocatable :: spheres
        integer, intent(in) :: nspheres
        integer :: i
        id_offset = n_sites
        allocate(spheres(nspheres))
        do i=1, nspheres
            call spheres(i)%init
            spheres(i)%id = i
        enddo
        n_sites = n_sites + nspheres
    end subroutine setup_sph

end module spheres
   

