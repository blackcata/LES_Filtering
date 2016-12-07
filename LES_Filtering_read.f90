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
              ONLY : X, Y, Z, U, V, W

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

        END SUBROUTINE READ_DNS
