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
          REAL(KIND=8) :: Del,dx,dz,FW
          CHARACTER(LEN=65) :: file_name, dir_name, path_name

          REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: X,Y,Z,dy
          REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: U,V,W,U_Fil,V_Fil,W_Fil

        END MODULE
