!------------------------------------------------------------------------------!
!
!   PROGRAM : LES_Filtering_filter.f90
!
!   PURPOSE : Filter datas from the DNS data for turbulent channel flow
!             using Gaussian filter.
!
!                                                             2016.12.07 K.Noh
!
!------------------------------------------------------------------------------!

        SUBROUTINE FILTER
          USE LES_FILTERING_module,                                             &
              ONLY : G, FIND_U, FIND_U_Fil, FIND_dU_Fil, FIND_dx

          USE LES_FILTERING_module,                                             &
              ONLY : Nx, Ny, Nz, dx, dz, Del, Nx_fil, Nz_fil

          USE LES_FILTERING_module,                                             &
              ONLY : X, Y, Z, U, V, W, dy, U_Fil, V_Fil, W_Fil,                 &
                     Resi_T, S_T_Fil, O_T_Fil, NU_R

          IMPLICIT NONE
          INTEGER :: i,j,k, i_loc, k_loc, i_tmp, k_tmp, x_i, x_j
          REAL(KIND=8) :: r, G_tot, time_sta, time_end, U_i, U_j, S,            &
                          dU_i, dU_j, dx_i, dx_j

          WRITE(*,*) '----------------------------------------------------'
          WRITE(*,*) '             FILTERING PROCESS STARTED              '
          WRITE(*,"(7X,A9,I4,5X,A9,I4)")"Nx : ",Nx,"Nz : ",Nz
          WRITE(*,"(7X,A9,I4,5X,A9,I4)")"Nx_fil : ",Nx_fil,"Nz_fil : ",Nz_fil
          CALL CPU_TIME(time_sta)

          !--------------------------------------------------------------------!
          !                  Main loop of filtering velocities                 !
          !--------------------------------------------------------------------!
          !$OMP PARALLEL DO private(k,j,i,k_loc,i_loc,i_tmp,k_tmp,x_i,x_j,U_i,U_j,G_tot,r,S)
          DO j = 1,Ny
            DO k = 1,Nz
              DO i = 1,Nx
                G_tot = 0.0

                !--------------------------------------------------------------!
                !                Sub loop of horizontal filtering              !
                !--------------------------------------------------------------!
                DO k_loc = -Nz_fil,Nz_fil
                  DO i_loc = -Nx_fil,Nx_fil
                    i_tmp = i + i_loc
                    k_tmp = k + k_loc

                    !----------------------------------------------------------!
                    !               Periodic boundary conditions               !
                    !----------------------------------------------------------!
                    IF (i_tmp > Nx) THEN
                      i_tmp = i_tmp - Nx
                    ELSEIF( i_tmp <= 0) THEN
                      i_tmp = i_tmp + Nx
                    END IF

                    IF (k_tmp > Nz) THEN
                      k_tmp = k_tmp - Nz
                    ELSEIF( k_tmp <= 0) THEN
                      k_tmp = k_tmp + Nz
                    END IF

                    !----------------------------------------------------------!
                    !                Gaussian Filter Calculation               !
                    !----------------------------------------------------------!
                    r = sqrt( (X(i)-X(i_tmp))**2 + (Z(k)-Z(k_tmp))**2 )
                    G_tot = G_tot + G(Del,r)

                    !----------------------------------------------------------!
                    !                        Filtering U_i                     !
                    !----------------------------------------------------------!
                    U_Fil(i,j,k) = U_Fil(i,j,k) + G(Del,r)* U(i_tmp,j,k_tmp)
                    V_Fil(i,j,k) = V_Fil(i,j,k) + G(Del,r)* V(i_tmp,j,k_tmp)
                    W_Fil(i,j,k) = W_Fil(i,j,k) + G(Del,r)* W(i_tmp,j,k_tmp)

                    !----------------------------------------------------------!
                    !                   Filtering (U_i) X (U_j)                !
                    !----------------------------------------------------------!
                    DO x_i = 1,3
                      DO x_j = 1,3
                        U_i = FIND_U(i_tmp,j,k_tmp,x_i)
                        U_j = FIND_U(i_tmp,j,k_tmp,x_j)
                        Resi_T(i,j,k,x_i,x_j) = Resi_T(i,j,k,x_i,x_j)           &
                                              + G(Del,r)*U_i*U_j
                      END DO
                    END DO

                  END DO
                END DO

                !--------------------------------------------------------------!
                !                        Weight adjustment                     !
                !--------------------------------------------------------------!
                U_Fil(i,j,k) = U_Fil(i,j,k)/G_tot
                V_Fil(i,j,k) = V_Fil(i,j,k)/G_tot
                W_Fil(i,j,k) = W_Fil(i,j,k)/G_tot

                DO x_i = 1,3
                  DO x_j = 1,3
                    U_i = FIND_U_Fil(i,j,k,x_i)
                    U_j = FIND_U_Fil(i,j,k,x_j)
                    Resi_T(i,j,k,x_i,x_j) = Resi_T(i,j,k,x_i,x_j)/G_tot         &
                                          - U_i*U_j
                  END DO
                END DO

                ! WRITE(*,"(3(I10),3(F15.9))")                                  &
                !                i,j,k,U_Fil(i,j,k),V_Fil(i,j,k),W_Fil(i,j,k)
              END DO
            END DO
          END DO
          !OMP END PARALLEL

          !--------------------------------------------------------------------!
          !                         Filtering Strain rate                      !
          !--------------------------------------------------------------------!
          !$OMP PARALLEL DO private(k,j,i,x_i,x_j,U_i,U_j,S)
          DO k = 2,Nz-1
            DO j = 2,Ny-1
              DO i = 2,Nx-1

                S = 0.0
                DO x_i = 1,3
                  DO x_j = 1,3
                    dU_i = FIND_dU_Fil(i,j,k,x_i,x_j)
                    dU_j = FIND_dU_Fil(i,j,k,x_j,x_i)
                    dx_i = FIND_dx(i,j,k,x_i)
                    dx_j = FIND_dx(i,j,k,x_j)

                    S_T_Fil(i,j,k,x_i,x_j) = ( dU_i / dx_j + dU_j / dx_i ) / 2.0
                    O_T_Fil(i,j,k,x_i,x_j) = ( dU_i / dx_j - dU_j / dx_i ) / 2.0
                    S  = S + ( S_T_Fil(i,j,k,x_i,x_j) )**2
                  END DO
                END DO
                S = sqrt(2*S)
                NU_R(i,j,k,1:3,1:3) = - Resi_T(i,j,k,1:3,1:3)                    &
                                      * S_T_Fil(i,j,k,1:3,1:3)/S**2

              END DO
            END DO
          END DO
          !OMP END PARALLEL

          CALL CPU_TIME(time_end)
          WRITE(*,*) '              FILTERING PROCESS ENDED               '
          WRITE(*,*) ' Total Filtering time : ',time_end - time_sta,' s'
          WRITE(*,*) '----------------------------------------------------'
          WRITE(*,*) ''

        END SUBROUTINE FILTER
