!-----------------------------------------------------------------------------------!
!
!   PROGRAM : LES_Filtering_main.f90
!
!   PURPOSE : To LES-filter the DNS data for turbulent channel flow.
!             (a) : obtain a filtered velocity field by applying Gaussian filter
!             (b) : obtain the residual-stress tensor
!             (c) : obtain the residual viscosity
!             (d) : obtain averaged Smagorinsky coefficient
!
!                                                             2016.12.07 K.Noh
!
!-----------------------------------------------------------------------------------!

        PROGRAM LES_FILTERING

          IMPLICIT NONE
          REAL(KIND=8) :: time_sta, time_end

          CALL SETUP
          CALL READ_DNS
          CALL FILTER
          CALL OUTPUT

        END PROGRAM LES_FILTERING
