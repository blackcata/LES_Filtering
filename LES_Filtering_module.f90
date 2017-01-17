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
          INTEGER :: N, Nx, Ny, Nz, Nx_fil, Nz_fil, VS_CASE, FILTER_OX, VS_ONLY,&
                     Y_ORDER
          REAL(KIND=8) :: Del,dx,dz,FW,pi,tol
          CHARACTER(LEN=65) :: file_name, dir_name, path_name

          REAL(KIND=8),DIMENSION(:),ALLOCATABLE :: X,Y,Z,dy,YP
          REAL(KIND=8),DIMENSION(:,:,:),ALLOCATABLE :: U,V,W,U_Fil,V_Fil,W_Fil, &
                                                       U_Fil_2,V_Fil_2,W_Fil_2, &
                                                       VS,Cs
          REAL(KIND=8),DIMENSION(:,:,:,:,:),ALLOCATABLE :: Resi_T,S_T,S_T_Fil,  &
                                                           Nu_R,O_T,O_T_Fil,    &
                                                           S_T_Fil_2,L_T,M_T

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
            FUNCTION FIND_U(i,j,k,v_i)
              INTEGER,INTENT(IN) :: i,j,k,v_i
              REAL(KIND=8) :: FIND_U

              SELECT CASE(v_i)
                CASE(1); FIND_U = U(i,j,k)
                CASE(2); FIND_U = V(i,j,k)
                CASE(3); FIND_U = W(i,j,k)
              END SELECT

            END FUNCTION FIND_U

            !------------------------------------------------------------------!
            !                        Filtered U Selection                      !
            !------------------------------------------------------------------!
            FUNCTION FIND_U_Fil(i,j,k,v_i)
              INTEGER,INTENT(IN) :: i,j,k,v_i
              REAL(KIND=8) :: FIND_U_Fil

              SELECT CASE(v_i)
                CASE(1); FIND_U_Fil = U_Fil(i,j,k)
                CASE(2); FIND_U_Fil = V_Fil(i,j,k)
                CASE(3); FIND_U_Fil = W_Fil(i,j,k)
              END SELECT

            END FUNCTION FIND_U_Fil


            !------------------------------------------------------------------!
            !                  Second Filtered U Selection                     !
            !------------------------------------------------------------------!
            FUNCTION FIND_U_Fil_2(i,j,k,v_i)
              INTEGER,INTENT(IN) :: i,j,k,v_i
              REAL(KIND=8) :: FIND_U_Fil_2

              SELECT CASE(v_i)
                CASE(1); FIND_U_Fil_2 = U_Fil_2(i,j,k)
                CASE(2); FIND_U_Fil_2 = V_Fil_2(i,j,k)
                CASE(3); FIND_U_Fil_2 = W_Fil_2(i,j,k)
              END SELECT

            END FUNCTION FIND_U_Fil_2

            !------------------------------------------------------------------!
            !                            dU Selection                          !
            !------------------------------------------------------------------!
            FUNCTION FIND_dU(i,j,k,x_i,x_j)
              INTEGER,INTENT(IN) :: i,j,k,x_i,x_j
              REAL(KIND=8) :: FIND_dU, U_ij(0:2,0:2,0:2)

              SELECT CASE(x_i)
                CASE(1); U_ij(0:2,0:2,0:2) = U(i-1:i+1,j-1:j+1,k-1:k+1)
                CASE(2); U_ij(0:2,0:2,0:2) = V(i-1:i+1,j-1:j+1,k-1:k+1)
                CASE(3); U_ij(0:2,0:2,0:2) = W(i-1:i+1,j-1:j+1,k-1:k+1)
              END SELECT

              SELECT CASE(x_j)
                CASE(1); FIND_dU = (U_ij(2,1,1) - U_ij(0,1,1))
                CASE(2); FIND_dU = (U_ij(1,2,1) - U_ij(1,0,1))
                CASE(3); FIND_dU = (U_ij(1,1,2) - U_ij(1,1,0))
              END SELECT

            END FUNCTION FIND_dU

            !------------------------------------------------------------------!
            !                          dU_Fil Selection                        !
            !------------------------------------------------------------------!
            FUNCTION FIND_dU_Fil(i,j,k,x_i,x_j)
              INTEGER,INTENT(IN) :: i,j,k,x_i,x_j
              REAL(KIND=8) :: FIND_dU_Fil, U_ij(0:2,0:2,0:2)

              SELECT CASE(x_i)
                CASE(1); U_ij(0:2,0:2,0:2) = U_Fil(i-1:i+1,j-1:j+1,k-1:k+1)
                CASE(2); U_ij(0:2,0:2,0:2) = V_Fil(i-1:i+1,j-1:j+1,k-1:k+1)
                CASE(3); U_ij(0:2,0:2,0:2) = W_Fil(i-1:i+1,j-1:j+1,k-1:k+1)
              END SELECT

              SELECT CASE(x_j)
                CASE(1); FIND_dU_Fil = (U_ij(2,1,1) - U_ij(0,1,1))
                CASE(2); FIND_dU_Fil = (U_ij(1,2,1) - U_ij(1,0,1))
                CASE(3); FIND_dU_Fil = (U_ij(1,1,2) - U_ij(1,1,0))
              END SELECT

            END FUNCTION FIND_dU_Fil

            !------------------------------------------------------------------!
            !                        dU_Fil_2 Selection                        !
            !------------------------------------------------------------------!
            FUNCTION FIND_dU_Fil_2(i,j,k,x_i,x_j)
              INTEGER,INTENT(IN) :: i,j,k,x_i,x_j
              REAL(KIND=8) :: FIND_dU_Fil_2, U_ij(0:2,0:2,0:2)

              SELECT CASE(x_i)
                CASE(1); U_ij(0:2,0:2,0:2) = U_Fil_2(i-1:i+1,j-1:j+1,k-1:k+1)
                CASE(2); U_ij(0:2,0:2,0:2) = V_Fil_2(i-1:i+1,j-1:j+1,k-1:k+1)
                CASE(3); U_ij(0:2,0:2,0:2) = W_Fil_2(i-1:i+1,j-1:j+1,k-1:k+1)
              END SELECT

              SELECT CASE(x_j)
                CASE(1); FIND_dU_Fil_2 = (U_ij(2,1,1) - U_ij(0,1,1))
                CASE(2); FIND_dU_Fil_2 = (U_ij(1,2,1) - U_ij(1,0,1))
                CASE(3); FIND_dU_Fil_2 = (U_ij(1,1,2) - U_ij(1,1,0))
              END SELECT

            END FUNCTION FIND_dU_Fil_2

            !------------------------------------------------------------------!
            !                            dx Selection                          !
            !------------------------------------------------------------------!
            FUNCTION FIND_dx(i,j,k,x_j)
              INTEGER,INTENT(IN) :: i,j,k,x_j
              REAL(KIND=8) :: FIND_dx

              SELECT CASE(x_j)
                CASE(1); FIND_dx = 2*dx
                CASE(2); FIND_dx = (dy(j)+dy(j+1))
                CASE(3); FIND_dx = 2*dz
              END SELECT

            END FUNCTION FIND_dx

        END MODULE
