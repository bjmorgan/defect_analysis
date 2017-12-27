module cell

    double precision, dimension(3,3) :: h
    double precision, dimension(3) :: boxlen, halfboxlen, cboxlen, halfcboxlen, cpplane

contains

    pure double precision function relr(r,i)
        ! returns the relative distance along direction i of point r from the (lower) cell boundary with index i
        integer, intent(in) :: i
        double precision, dimension(3), intent(in) :: r
    
        select case (i)
            case (1)
                relr = r(1) - r(2)*h(2,1)/h(2,2) - r(3)*h(3,1)/h(3,3)
            case (2)
                relr = r(2) - r(1)*h(1,2)/h(1,1) - r(3)*h(3,2)/h(3,3)
            case (3)
                relr = r(3) - r(1)*h(1,3)/h(1,1) - r(2)*h(2,3)/h(2,2)
        end select
    end function relr

    pure function r_as_minimum_image( r )
        double precision, dimension(3) :: temp_r
        double precision, dimension(3), intent(in) :: r
        double precision, dimension(3) :: r_as_minimum_image
        integer :: k

        temp_r = r

        do k=1, 3
            if ( temp_r(k) < -halfcboxlen(k) ) then
                temp_r(:) = temp_r(:) + boxlen(k) * h(k,:)
            else if ( temp_r(k) > halfcboxlen(k) ) then
                temp_r(:) = temp_r(:) - boxlen(k) * h(k,:)
            end if
        end do

        r_as_minimum_image = temp_r
    end function r_as_minimum_image

    pure function dr( r1, r2 )
        double precision, dimension(3), intent(in) :: r1, r2
        double precision, dimension(3) :: dr

        dr = r1 - r2

    end function dr

    pure function move_inside_cell( r ) result( new_r )
        double precision, dimension(3), intent(in) :: r
        double precision, dimension(3) :: new_r, temp_r, mid_cell_r
        mid_cell_r = matmul( h, halfboxlen )
        temp_r = r - mid_cell_r
        new_r = r_as_minimum_image( temp_r ) + mid_cell_r
    end function move_inside_cell

end module cell
