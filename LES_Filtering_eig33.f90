!------------------------------------------------------------------------------!
!
!   PROGRAM : LES_Filtering_eig33.f90
!
!   PURPOSE : Getting eigenvalue of 3 X 3 matrix
!
!                                                             2016.12.17 K.Noh
!
!------------------------------------------------------------------------------!

        SUBROUTINE EIG33(A,eig)

          USE LES_FILTERING_module,                                             &
              ONLY : VS_CASE

          IMPLICIT NONE

          REAL(KIND=8):: A(3,3)
          REAL(KIND=8) :: B(3,3),I(3,3),a0,a1,a2,a3,det,tr
          COMPLEX(KIND=8) :: eig(3)

          CALL MMDOT(A,A,B,3)
          !---------------------------------------------------------------------!
          !                   Trace function of 3 X 3 matrix                     !
          !----------------------------------------------------------------------!
          a0 = det(A)
          a1 = 0.5*(tr(B) - tr(A)**2)
          a2 = tr(A)
          a3 = -1

          IF (VS_CASE == 2) CALL CUBIC_2(a3,a2,a1,a0,eig)
          IF (VS_CASE == 3) CALL CUBIC(a3,a2,a1,a0,eig)

        END SUBROUTINE EIG33

        !----------------------------------------------------------------------!
        !                 Determinant Function of 3 X 3 matrix                 !
        !----------------------------------------------------------------------!
        FUNCTION det(A)
            REAL(KIND=8) :: det
            REAL(KIND=8) :: A(3,3)

            det = A(1,1)*A(2,2)*A(3,3) + A(1,2)*A(2,3)*A(3,1) &
                + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1) &
                - A(1,2)*A(2,1)*A(3,3) - A(1,1)*A(2,3)*A(3,2)

        END FUNCTION det

        !----------------------------------------------------------------------!
        !                   Trace function of 3 X 3 matrix                     !
        !----------------------------------------------------------------------!
        FUNCTION tr(A)
            REAL(KIND=8) :: tr
            REAL(KIND=8) :: A(3,3)

            tr   = A(1,1) + A(2,2) + A(3,3)

        END FUNCTION tr

        !----------------------------------------------------------------------!
        !                 Matrix,matrix multiplication function                !
        !----------------------------------------------------------------------!
        SUBROUTINE MMDOT(A,B,C,n)
            INTEGER,INTENT(IN) :: n
            REAL(KIND=8),INTENT(INOUT) :: A(n,n), B(n,n), C(n,n)
            INTEGER :: i,j,k

            C(1:n,1:n) = 0.0
            DO i = 1,n
              DO j = 1,n
                DO k = 1,n
                  C(i,j) = C(i,j) + A(i,k) * B(k,j)
                END DO
              END DO
            END DO
        END SUBROUTINE MMDOT

        !----------------------------------------------------------------------!
        !                   General cubic equation solver                      !
        !----------------------------------------------------------------------!
        SUBROUTINE CUBIC(a,b,c,d,x)

            IMPLICIT NONE
            REAL(KIND=8),INTENT(IN) :: a,b,c,d
            REAL(KIND=8) :: S1_tmp,S2_tmp,Q,R
            COMPLEX(KIND=8) ::i,x(1:3),S1,S2
            i = CMPLX(0,1)

            Q  = 2.*b**3. - 9.*a*b*c + 27.*a**2.*d
            R  = b**2 - 3.*a*c

            x(1) = - b/(3.*a)                                                   &
                 - 1./(3.*a)*((Q + sqrt(Q**2. - 4.*R**3.))/2.)**(1.0d0/3.0d0)   &
                 - 1./(3.*a)*((Q - sqrt(Q**2. - 4.*R**3.))/2.)**(1.0d0/3.0d0)

            x(2) = - b/(3.*a)                                                   &
                 + (1.-i*sqrt(3.))/(6.*a)*((Q + sqrt(Q**2. - 4.*R**3.))/2.)**(1.0d0/3.0d0)         &
                 + (1.+i*sqrt(3.))/(6.*a)*((Q - sqrt(Q**2. - 4.*R**3.))/2.)**(1.0d0/3.0d0)

            x(3) = - b/(3.*a)                                                   &
                 + (1.+i*sqrt(3.))/(6.*a)*((Q + sqrt(Q**2. - 4.*R**3.))/2.)**(1.0d0/3.0d0)         &
                 + (1.-i*sqrt(3.))/(6.*a)*((Q - sqrt(Q**2. - 4.*R**3.))/2.)**(1.0d0/3.0d0)

        END SUBROUTINE CUBIC

        SUBROUTINE CUBIC_2(a,b,c,d,x)

          IMPLICIT NONE
          REAL(KIND=8),INTENT(IN) :: a,b,c,d
          COMPLEX(KIND=8) ::i,x(1:3),Q,R,S1,S2
          i = CMPLX(0,1)

          Q = 2*b**3 - 9*a*b*c + 27*a**2*d
          R = sqrt(Q**2 -4*(b**2-3*a*c)**3)
          S1 = ( 0.5*(Q+R) )**(1./3.)
          S2 = ( 0.5*(Q-R) )**(1./3.)

          x(1) = -b/(3*a) - 1./(3.*a)*S1- 1./(3.*a)*S2
          x(2) = -b/(3*a) +(1+i*3**0.5)/(6*a)*S1 +(1-i*3**0.5)/(6*a)*S2
          x(3) = -b/(3*a) +(1-i*3**0.5)/(6*a)*S1 +(1+i*3**0.5)/(6*a)*S2

      END SUBROUTINE CUBIC_2
