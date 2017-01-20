!------------------------------------------------------------------------------!
!                                                                              !
!   PROGRAM : Vortical_Structure.f90                                           !
!                                                                              !
!   PURPOSE : Make the vortircal structures by using 3 methods                 !
!               (a) Q criteria                                                 !
!               (b) Lambda_2 criteria                                          !
!               (c) Lambda_ci criteria                                         !
!                                                                              !
!                                                             2016.12.17 K.Noh !
!                                                                              !
!------------------------------------------------------------------------------!

        SUBROUTINE VORTICAL_STRUCTURE

          USE LES_FILTERING_module,                                             &
              ONLY : FIND_dU, FIND_dx

          USE LES_FILTERING_module,                                             &
              ONLY : Nx, Ny, Nz, VS_CASE, tol

          USE LES_FILTERING_module,                                             &
              ONLY : S_T, O_T, VS

          IMPLICIT NONE
          EXTERNAL DSYEV

          INTEGER :: i,j,k,x_i,x_j,INFO,LWORK,LDA,LWMAX
          REAL(KIND=8) :: time_sta, time_end, dU_i, dU_j, dx_i, dx_j, P, Q, R   &
                                            , Q_hat, R_hat, Delta, det, tr
          REAL(KIND=8) :: D_T_tmp(3,3), Q_loc(3,3), DD(3,3), WORK(1000), EIG_R(3)&
                         ,S_T_tmp(3,3), O_T_tmp(3,3)
          COMPLEX(KIND=8) :: eig(3)

          WRITE(*,*) '----------------------------------------------------'
          WRITE(*,*) '     VORTICAL STRUCTURE MAKING PROCESS STARTED      '
          CALL CPU_TIME(time_sta)

          SELECT CASE(VS_CASE)
            CASE(1)
          !--------------------------------------------------------------------!
          !             Main loop of making vortical structures - Q            !
          !--------------------------------------------------------------------!
          !$OMP PARALLEL DO private(k,j,i,Q_loc,x_i,x_j,dU_i,dU_j,dx_i,dx_j)
            DO k = 2,Nz-1
              DO j = 2,Ny-1
                DO i = 2,Nx-1

                  DO x_j = 1,3
                    DO x_i = 1,3
                      dU_i = FIND_dU(i,j,k,x_i,x_j)
                      dU_j = FIND_dU(i,j,k,x_j,x_i)
                      dx_i = FIND_dx(i,j,k,x_i)
                      dx_j = FIND_dx(i,j,k,x_j)

                      Q_loc(x_i,x_j) = -0.5*(dU_i/dx_j)*(dU_j/dx_i)
                    END DO
                  END DO

                  VS(i,j,k) = SUM(Q_loc(1:3,1:3))
                END DO
              END DO
            END DO
            !OMP END PARALLEL

          CASE(2)
          !--------------------------------------------------------------------!
          !         Main loop of making vortical structures - Lambda_2         !
          !--------------------------------------------------------------------!
            LWMAX = 1000

            !$OMP PARALLEL DO private(k,j,i,D_T_tmp,eig,x_i,x_j,LDA,WORK,LWORK,INFO,S_T_tmp,O_T_tmp)
            DO k = 2,Nz-1
              DO j = 2,Ny-1
                DO i = 2,Nx-1

                  CALL MMDOT(S_T(i,j,k,1:3,1:3),S_T(i,j,k,1:3,1:3),S_T_tmp(1:3,1:3),3)
                  CALL MMDOT(O_T(i,j,k,1:3,1:3),O_T(i,j,k,1:3,1:3),O_T_tmp(1:3,1:3),3)
                  DO x_j = 1,3
                    DO x_i = 1,3
                      D_T_tmp(x_i,x_j) = S_T_tmp(x_i,x_j) + O_T_tmp(x_i,x_j)
                    END DO
                  END DO

                  ! LWORK = -1
                  ! CALL DSYEV('N','U',3,D_T_tmp,3,EIG_R,WORK,LWORK,INFO)
                  ! LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
                  ! CALL DSYEV('N','U',3,D_T_tmp,3,EIG_R,WORK,LWORK,INFO)
                  ! VS(i,j,k) = EIG_R(2)

                  ! WRITE(*,"(3I,6F15.9)") i,j,k,REAL(eig(1:3)),AIMAG(eig(1:3))
                  CALL EIG33(D_T_tmp,eig)
                  VS(i,j,k) = REAL(eig(2))

                END DO
              END DO
            END DO
            !OMP END PARALLEL

          CASE(3)
          !--------------------------------------------------------------------!
          !         Main loop of making vortical structures - Lambda_ci        !
          !--------------------------------------------------------------------!
            !$OMP PARALLEL DO private(k,j,i,x_i,x_j,D_T_tmp,eig)
            DO k = 2,Nz-1
              DO j = 2,Ny-1
                DO i = 2,Nx-1

                  DO x_j = 1,3
                    DO x_i = 1,3
                      D_T_tmp(x_i,x_j) = S_T(i,j,k,x_i,x_j) + O_T(i,j,k,x_i,x_j)
                    END DO
                  END DO

                  CALL EIG33(D_T_tmp,eig)
                  IF (abs( IMAG(eig(2) + eig(3)) ) < tol) VS(i,j,k) = AIMAG(eig(2))**2

                END DO
              END DO
            END DO
            !OMP END PARALLEL
            END SELECT

          CALL CPU_TIME(time_end)
          WRITE(*,*) '   VORTICAL STRUCTURE MAKING PROCESS IS COMPLETED   '
          WRITE(*,*) '  Total Making time : ',time_end - time_sta,' s'
          WRITE(*,*) '----------------------------------------------------'
          WRITE(*,*) ''

        END SUBROUTINE VORTICAL_STRUCTURE
