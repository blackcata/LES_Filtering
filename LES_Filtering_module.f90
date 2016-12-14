!------------------------------------------------------------------------------!
!
!   PROGRAM : LES_Filtering_module.f90
!
!   PURPOSE : Module for LES-filter the DNS data for turbulent channel flow.
!
!                                                             2016.12.07 K.Noh
!
!------------------------------------------------------------------------------!

        MODULE LES_FILTERING_module

          IMPLICIT NONE
          INTEGER :: Nx, Ny, Nz, Nx_fil, Nz_fil
          REAL(KIND=8) :: Del,dx,dz,FW,pi,tol
          CHARACTER(LEN=65) :: file_name, dir_name, path_name

          REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: X,Y,Z,dy
          REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: U,V,W,U_Fil,V_Fil,W_Fil
          REAL(KIND=8),DIMENSION(:,:,:,:,:),ALLOCATABLE :: Resi_T

          CONTAINS
            !------------------------------------------------------------------!
            !                     Gaussian Filter Function                     !
            !------------------------------------------------------------------!
            FUNCTION G(Del,r)
              REAL(KIND=8) :: G
              REAL(KIND=8),INTENT(IN) :: Del, r

              G = sqrt(6/(pi*Del**2))*exp(-6*(r**2/Del**2))

            END FUNCTION G

            !------------------------------------------------------------------!
            !                     Finding Y point Function                     !
            !------------------------------------------------------------------!
            FUNCTION J_det(y_val)
              INTEGER :: j, J_det
              REAL(KIND=8),INTENT(IN) :: y_val

              DO j = 1,Ny
                IF (Y(j)>y_val) THEN
                  IF ( (Y(j) - y_val) < (y_val - Y(j-1))) THEN
                    J_det = j
                  ELSE
                    J_det = j-1
                  END IF
                  EXIT
                END IF

              END DO

            END FUNCTION J_det

            !------------------------------------------------------------------!
            !                             U Selection                          !
            !------------------------------------------------------------------!
            FUNCTION FIND_U(i_tmp,j,k_tmp,v_i)
              INTEGER,INTENT(IN) :: i_tmp,j,k_tmp,v_i
              REAL(KIND=8) :: FIND_U

              SELECT CASE(v_i)
                CASE(1)
                  FIND_U = U(i_tmp,j,k_tmp)
                CASE(2)
                  FIND_U = V(i_tmp,j,k_tmp)
                CASE(3)
                  FIND_U = W(i_tmp,j,k_tmp)
              END SELECT

            END FUNCTION FIND_U

            !------------------------------------------------------------------!
            !                        Filtered U Selection                      !
            !------------------------------------------------------------------!
            FUNCTION FIND_U_Fil(i,j,k,v_i)
              INTEGER,INTENT(IN) :: i,j,k,v_i
              REAL(KIND=8) :: FIND_U_Fil

              SELECT CASE(v_i)
                CASE(1)
                  FIND_U_Fil = U_Fil(i,j,k)
                CASE(2)
                  FIND_U_Fil = V_Fil(i,j,k)
                CASE(3)
                  FIND_U_Fil = W_Fil(i,j,k)
              END SELECT

            END FUNCTION FIND_U_Fil

        END MODULE
