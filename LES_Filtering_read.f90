!-----------------------------------------------------------------------------------!
!
!   PROGRAM : LES_Filtering_read.f90
!
!   PURPOSE : Reading datas from the DNS data for turbulent channel flow.
!
!                                                             2016.12.07 K.Noh
!
!-----------------------------------------------------------------------------------!

        SUBROUTINE READ_DNS
          USE LES_FILTERING_module,                                           &
              ONLY : Nx, Ny, Nz, file_name, dir_name, path_name

          USE LES_FILTERING_module,                                           &
              ONLY : X, Y, Z, U, V, W, dx, dy, dz

          IMPLICIT NONE

          INTEGER :: i,j,k
          REAL(KIND=8) :: tmp_x, tmp_y, tmp_z
          CHARACTER(20) :: header

          dir_name = 'Data'
          path_name = TRIM(dir_name)//'/'//TRIM(file_name)

          OPEN(100,FILE=path_name,FORM='FORMATTED',STATUS='OLD')
          READ(100,*) header
          READ(100,*) header

          DO k = 1,Nz
            DO j = 1,Ny
              DO i = 1,Nx
                READ(100,*) tmp_x, tmp_y, tmp_z, U(i,j,k), V(i,j,k), W(i,j,k)

                IF (j==1 .AND. k==1) X(i)      = tmp_x
                IF (i==1 .AND. k==1) Y(Ny-j+1) = tmp_y
                IF (i==1 .AND. j==1) Z(k)      = tmp_z
              END DO
            END DO
          END DO

          CLOSE(100)

          DO i = 1,Nx-1
            dx(i) = X(i+1) - X(i)
          END DO

          DO j = 1,Ny-1
            dy(j) = Y(j+1) - Y(j)
          END DO

          DO k = 1,Nz-1
            dz(k) = Z(k+1) - Z(k)
          END DO

        END SUBROUTINE READ_DNS
