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
              ONLY : G

          USE LES_FILTERING_module,                                             &
              ONLY : Nx, Ny, Nz, dx, dz, Del, Nx_fil, Nz_fil

          USE LES_FILTERING_module,                                             &
              ONLY : X, Y, Z, U, V, W, dy, U_Fil, V_Fil, W_Fil

          IMPLICIT NONE
          INTEGER :: i,j,k, i_loc, k_loc, i_tmp, k_tmp, myid, OMP_GET_THREAD_NUM
          REAL(KIND=8) :: r, G_tot, time_sta, time_end

          WRITE(*,*) '----------------------------------------------------'
          WRITE(*,*) '             FILTERING PROCESS STARTED              '
          WRITE(*,"(2X,A9,I2,1X,A9,I2)")"Nx : ",Nx,"Nz : ",Nz
          WRITE(*,"(2X,A9,I2,1X,A9,I2)")"Nx_fil : ",Nx_fil,"Nz_fil : ",Nz_fil
          CALL CPU_TIME(time_sta)

          !$OMP PARALLEL DO private(k,i,k_loc,i_loc,i_tmp,k_tmp,G_tot,myid)
          DO j = 1,Ny
            DO k = 1,Nz
              DO i = 1,Nx
                G_tot = 0.0
                DO k_loc = -Nz_fil,Nz_fil!-Nz/2,Nz/2
                  DO i_loc = -Nx_fil,Nx_fil!-Nx/2,Nx/2
                    i_tmp = i + i_loc
                    k_tmp = k + k_loc

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

                    r = sqrt( (X(i)-X(i_tmp))**2 + (Z(k)-Z(k_tmp))**2 )
                    G_tot = G_tot + G(Del,r)
                    U_Fil(i,j,k) = U_Fil(i,j,k) + G(Del,r)* U(i_tmp,j,k_tmp)
                    V_Fil(i,j,k) = V_Fil(i,j,k) + G(Del,r)* V(i_tmp,j,k_tmp)
                    W_Fil(i,j,k) = W_Fil(i,j,k) + G(Del,r)* W(i_tmp,j,k_tmp)
                  END DO
                END DO
                U_Fil(i,j,k) = U_Fil(i,j,k)/G_tot
                V_Fil(i,j,k) = V_Fil(i,j,k)/G_tot
                W_Fil(i,j,k) = W_Fil(i,j,k)/G_tot

                ! WRITE(*,"(3(I10),3(F15.9))")                                  &
                !                i,j,k,U_Fil(i,j,k),V_Fil(i,j,k),W_Fil(i,j,k)
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
