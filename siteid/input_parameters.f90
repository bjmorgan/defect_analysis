module input_parameters

    implicit none

    type input_parameter_set
        character(len=30) :: posfile, cellfile, inptfile, tetfile, octfile, sphfile
        integer :: nconfigs, nspec, lattice_spec, mobile_spec
        integer, allocatable, dimension(:) :: nsp
        integer :: ntet, noct, nspheres
        logical :: variable_cell
        double precision, dimension(3) :: boxlen
        double precision, dimension(3,3) :: h
    contains
        procedure :: read_from_file
    end type input_parameter_set

    contains

    subroutine read_from_file( params, inptfile )
        class(input_parameter_set) :: params
        character(len=30), intent(in) :: inptfile
        integer :: fin, i

        open( file=inptfile, status='old', newunit=fin )
        read( fin, * ) params%posfile
        read( fin, * ) params%cellfile
        read( fin, * ) params%tetfile
        read( fin, * ) params%octfile
        read( fin, * ) params%sphfile
        read( fin, * ) params%nconfigs
        read( fin, * ) params%nspec
        allocate (params%nsp(params%nspec))
        do i=1, params%nspec
            read(fin,*) params%nsp(i)
        enddo
        read( fin, * ) params%lattice_spec
        read( fin, * ) params%mobile_spec
        read( fin, * ) params%ntet
        read( fin, * ) params%noct
        read( fin, * ) params%nspheres
        read( fin, * ) params%variable_cell
        read( fin, * ) params%boxlen(:)
        read( fin, * ) params%h(:,:)
        close(fin)
    end subroutine read_from_file

end module input_parameters

        

