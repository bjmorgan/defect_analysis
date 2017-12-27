module sites
    
    use atoms
    use cell

    implicit none
   
    integer, parameter :: string_length = 11
    integer :: n_sites = 0
 
    type, abstract :: site
        double precision, dimension(3) :: centre
        logical :: occupied 
        integer :: occnum
        type(atom) :: occ_atom
        integer :: id
        character(len=string_length) :: string
    contains
        procedure(occupied_by_sub), deferred :: occupied_by
        procedure(contains_atom_sub), deferred, nopass :: contains_atom
        procedure(unique_id_func), deferred :: unique_id
    end type site

    abstract interface
        subroutine contains_atom_sub(this_atom)
            import :: atom
            class(atom) :: this_atom
        end subroutine contains_atom_sub
    
        integer function unique_id_func(this)
            import :: site
            class(site), intent(in) :: this
        end function unique_id_func

        subroutine occupied_by_sub(this, this_atom)
            import :: atom, site
            class(site) :: this
            type(atom) :: this_atom
        end subroutine occupied_by_sub
    
    end interface

end module sites
