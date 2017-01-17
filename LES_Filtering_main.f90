!------------------------------------------------------------------------------!
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
!------------------------------------------------------------------------------!

        PROGRAM LES_FILTERING
          USE LES_FILTERING_module,                                             &
              ONLY : FILTER_OX, VS_ONLY

          IMPLICIT NONE
          REAL(KIND=8) :: time_sta, time_end

          CALL SETUP
          CALL READ_DNS

          IF ( FILTER_OX > 0) CALL FILTER
          IF ( FILTER_OX > 1) CALL SECOND_FILTER

          IF ( mod(VS_ONLY,2) == 0 ) CALL VORTICAL_STRUCTURE

          CALL OUTPUT

        END PROGRAM LES_FILTERING
