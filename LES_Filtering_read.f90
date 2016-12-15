!------------------------------------------------------------------------------!
!
!   PROGRAM : LES_Filtering_read.f90
!
!   PURPOSE : Reading datas from the DNS data for turbulent channel flow.
!
!                                                             2016.12.07 K.Noh
!
!------------------------------------------------------------------------------!

        SUBROUTINE READ_DNS

          USE LES_FILTERING_module,                                             &
              ONLY : G, FIND_dU, FIND_dx

          USE LES_FILTERING_module,                                             &
              ONLY : Nx, Ny, Nz, dx, dz, Del, FW, tol, Nx_fil, Nz_fil,          &
                     file_name, dir_name, path_name

          USE LES_FILTERING_module,                                             &
              ONLY : X, Y, Z, U, V, W, dy, S_T

          IMPLICIT NONE

          INTEGER :: i,j,k,it,x_i,x_j
          REAL(KIND=8) :: tmp_x, tmp_y, tmp_z, time_sta, time_end, dg, r,       &
                          dU_i, dU_j, dx_i, dx_j
          CHARACTER(20) :: header

          WRITE(*,*) '----------------------------------------------------'
          WRITE(*,*) '              READING PROCESS STARTED               '
          CALL CPU_TIME(time_sta)

          dir_name = 'Data'
          path_name = TRIM(dir_name)//'/'//TRIM(file_name)

          OPEN(100,FILE=path_name,FORM='FORMATTED',STATUS='OLD')
          READ(100,*) header
          READ(100,*) header

          !--------------------------------------------------------------------!
          !                   Main loop of reading DNS datas                   !
          !--------------------------------------------------------------------!
          DO k = 1,Nz
            DO j = Ny,1,-1
              DO i = 1,Nx
                READ(100,*) tmp_x, tmp_y, tmp_z, U(i,j,k), V(i,j,k), W(i,j,k)

                IF (j==1 .AND. k==1) X(i)      = tmp_x
                IF (i==1 .AND. k==1) Y(j)      = tmp_y
                IF (i==1 .AND. j==1) Z(k)      = tmp_z
              END DO
            END DO
          END DO

          CLOSE(100)

          !--------------------------------------------------------------------!
          !                         Calculating dx,dy,dz                       !
          !--------------------------------------------------------------------!
          DO j = 1,Ny-1
            dy(j) = Y(j+1) - Y(j)
          END DO

          dx  = X(2) - X(1)
          dz  = Z(2) - Z(1)
          Del = FW*sqrt(dx * dz)

          !--------------------------------------------------------------------!
          !                      Calculating Stress Tensor                     !
          !--------------------------------------------------------------------!
          !$OMP PARALLEL DO private(k,j,i,x_i,x_j,dU_i,dU_j,dx_i,dx_j)
          DO k = 2,Nz-1
            DO j = 2,Ny-1
              DO i = 2,Nx-1

                DO x_i = 1,3
                  DO x_j = 1,3
                    dU_i = FIND_dU(i,j,k,x_i)
                    dU_j = FIND_dU(i,j,k,x_j)
                    dx_i = FIND_dx(i,j,k,x_i)
                    dx_j = FIND_dx(i,j,k,x_j)

                    S_T(i,j,k,x_i,x_j) = ( dU_i / dx_j + dU_j / dx_i ) / 2.0
                  END DO
                END DO

              END DO
            END DO
          END DO
          !OMP END PARALLEL

          S_T(1,1:Ny,1:Nz,1:3,1:3)  = S_T(2,1:Ny,1:Nz,1:3,1:3)
          S_T(Nx,1:Ny,1:Nz,1:3,1:3) = S_T(Nx-1,1:Ny,1:Nz,1:3,1:3)

          S_T(1:Nx,1,1:Nz,1:3,1:3)  = S_T(1:Nx,2,1:Nz,1:3,1:3)
          S_T(1:Nx,Ny,1:Nz,1:3,1:3) = S_T(1:Nx,Ny-1,1:Nz,1:3,1:3)

          S_T(1:Nx,1:Ny,1,1:3,1:3)  = S_T(1:Nx,1:Ny,2,1:3,1:3)
          S_T(1:Nx,1:Ny,Nz,1:3,1:3) = S_T(1:Nx,1:Ny,Nz-1,1:3,1:3)

          !--------------------------------------------------------------------!
          !       Determining the number of filtering node in x,z dirctions    !
          !--------------------------------------------------------------------!
          DO it = 1,1000
            r = REAL(it,KIND=8)
            dg = G(Del,r) - G(Del,r+1)
            IF (dg < tol) THEN
              Nx_fil = CEILING(REAL(it,KIND=8)/(2*dx))
              Nz_fil = CEILING(REAL(it,KIND=8)/(2*dz))
              EXIT
            END IF
          END DO

          CALL CPU_TIME(time_end)
          WRITE(*,*) '            READING PROCESS IS COMPLETED            '
          WRITE(*,*) '  Total Reading time : ',time_end - time_sta,' s'
          WRITE(*,*) '----------------------------------------------------'
          WRITE(*,*) ''

        END SUBROUTINE READ_DNS
