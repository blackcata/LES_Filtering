!-----------------------------------------------------------------------------------!
!
!   PROGRAM : LES_Filtering_module.f90
!
!   PURPOSE : Module for LES-filter the DNS data for turbulent channel flow.
!
!                                                             2016.12.07 K.Noh
!
!-----------------------------------------------------------------------------------!

        MODULE LES_FILTERING_module

          IMPLICIT NONE
          INTEGER :: Nx, Ny, Nz
          REAL(KIND=8) :: Del,dx,dz,FW,pi
          CHARACTER(LEN=65) :: file_name, dir_name, path_name

          REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: X,Y,Z,dy
          REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: U,V,W,U_Fil,V_Fil,W_Fil

          CONTAINS
            !------------------------------------------------------!
            !                Gaussian Filter Function
            !------------------------------------------------------!
            FUNCTION G(Del,r)
              REAL(KIND=8) :: G
              REAL(KIND=8),INTENT(IN) :: Del, r

              G = sqrt(6/(pi*Del**2))*exp(-6*(r**2/Del**2))

            END FUNCTION G

        END MODULE
