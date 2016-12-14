!------------------------------------------------------------------------------!
!
!   PROGRAM : LES_Filtering_output.f90
!
!   PURPOSE : Write each variables in the RESULT folder.
!
!                                                             2016.12.07 K.Noh
!
!------------------------------------------------------------------------------!
          SUBROUTINE OUTPUT
              USE LES_FILTERING_module,                                         &
                  ONLY : Nx, Ny, Nz, file_name, dir_name, path_name

              USE LES_FILTERING_module,                                         &
                  ONLY : X, Y, Z, U, V, W, U_Fil, V_Fil, W_Fil, dy

              IMPLICIT NONE
              INTEGER :: i,j,k
              REAL(KIND=8) :: time_sta, time_end, U_ave, V_ave, W_ave

              WRITE(*,*) '----------------------------------------------------'
              WRITE(*,*) '              WRITING PROCESS STARTED               '
              CALL CPU_TIME(time_sta)

              dir_name = 'RESULT'

              !----------------------------------------------------------------!
              !                Outputs for U_Fil(Filtered velocity)
              !----------------------------------------------------------------!
              file_name = '/U_filtered.plt'
              path_name = TRIM(dir_name)//TRIM(file_name)
              OPEN(100,FILE=path_name,FORM='FORMATTED',POSITION='APPEND')
              WRITE(100,*) 'VARIABLES = X,Y,Z,U_fil,V_Fil,W_Fil'
              WRITE(100,"(3(A,I3,2X))")' ZONE  I = ',Nx,' J = ',Ny, ' K = ', Nz
              DO k = 1,Nz
                DO j = 1,Ny
                  DO i = 1,Nx
                    WRITE(100,"(6F15.9)") X(i),Y(j),Z(k),                       &
                                          U_fil(i,j,k),V_Fil(i,j,k),W_Fil(i,j,k)
                  END DO
                END DO
              END DO
              CLOSE(100)

              !----------------------------------------------------------------!
              !            Outputs for Averaged U on x,z directions            !
              !----------------------------------------------------------------!
              file_name = '/U_averaged_profile.plt'
              path_name = TRIM(dir_name)//TRIM(file_name)
              OPEN(100,FILE=path_name,FORM='FORMATTED',POSITION='APPEND')
              WRITE(100,*) 'VARIABLES = Y,U_ave,V_ave,W_ave'
              WRITE(100,"(A,I3)")' J = ',Ny,

              DO j = 1,Ny
                U_ave = 0.0
                V_ave = 0.0
                W_ave = 0.0

                DO k = 1,Nz
                  DO i = 1,Nx
                    U_ave = U_ave + U(i,j,k)
                    V_ave = V_ave + V(i,j,k)
                    W_ave = W_ave + W(i,j,k)
                  END DO
                END DO
                U_ave = U_ave/(Nx*Nz)
                V_ave = V_ave/(Nx*Nz)
                W_ave = W_ave/(Nx*Nz)

                WRITE(100,"(4F15.9)")Y(j),U_ave,V_ave,W_ave
              END DO
              CLOSE(100)

              CALL CPU_TIME(time_end)

              WRITE(*,*) '           WRITING PROCESS IS COMPLETED            '
              WRITE(*,*) '  Total Writing time : ',time_end - time_sta,' s'
              WRITE(*,*) '----------------------------------------------------'
              WRITE(*,*) ''

              DEALLOCATE(X,Y,Z,U,V,W,U_Fil,V_Fil,W_Fil,dy)

          END SUBROUTINE OUTPUT
