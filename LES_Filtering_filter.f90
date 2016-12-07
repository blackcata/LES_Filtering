!-----------------------------------------------------------------------------------!
!
!   PROGRAM : LES_Filtering_filter.f90
!
!   PURPOSE : Filter datas from the DNS data for turbulent channel flow
!             using Gaussian filter.
!
!                                                             2016.12.07 K.Noh
!
!-----------------------------------------------------------------------------------!

        SUBROUTINE FILTER
          USE LES_FILTERING_module,                                           &
              ONLY : Nx, Ny, Nz, dx, dz, Del

          USE LES_FILTERING_module,                                           &
              ONLY : X, Y, Z, U, V, W, dy, U_Fil, V_Fil, W_Fil

          IMPLICIT NONE
          INTEGER :: i,j,k
          REAL(KIND=8) :: r


        END SUBROUTINE FILTER
