!-----------------------------------------------------------------------------------!
!
!   PROGRAM : LES_Filtering_output.f90
!
!   PURPOSE : Write each variables in the RESULT folder.
!
!                                                             2016.12.07 K.Noh
!
!-----------------------------------------------------------------------------------!
          SUBROUTINE OUTPUT
              USE LES_FILTERING_module,                                         &
                  ONLY : Nx, Ny, Nz, file_name, dir_name, path_name

              USE LES_FILTERING_module,                                         &
                  ONLY : X, Y, Z, U, V, W, U_Fil, V_Fil, W_Fil, dy

              IMPLICIT NONE
              INTEGER :: i,j,k
              REAL(KIND=8) :: time_sta, time_end

              WRITE(*,*) '----------------------------------------------------'
              WRITE(*,*) '              WRITING PROCESS STARTED               '
              CALL CPU_TIME(time_sta)

              dir_name = 'RESULT'
              !----------------------------------------------------------!
              !            Outputs for U_Fil(Filtered velocity)
              !----------------------------------------------------------!
              file_name = '/U_filtered.plt'
              path_name = TRIM(dir_name)//TRIM(file_name)
              OPEN(100,FILE=file_name,FORM='FORMATTED',POSITION='APPEND')
              DO k = 1,Nz
                DO j = 1,Ny
                  DO i = 1,Nx
                    WRITE(100,"(6F15.9)") X(i),Y(j),Z(k),                       &
                                          U_fil(i,j,k),V_Fil(i,j,k),W_Fil(i,j,k)
                  END DO
                END DO
              END DO
              CLOSE(100)

              WRITE(*,*) '            WRITING PROCESS IS COMPLETED            '
              WRITE(*,*) '  Total Writing time : ',time_end - time_sta,' s'
              WRITE(*,*) '----------------------------------------------------'
              WRITE(*,*) ''

              DEALLOCATE(X,Y,Z,U,V,W,U_Fil,V_Fil,W_Fil,dy)

          END SUBROUTINE OUTPUT
