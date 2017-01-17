!------------------------------------------------------------------------------!
!
!   PROGRAM : LES_Filtering_setup.f90
!
!   PURPOSE : Setup for LES-filter the DNS data for turbulent channel flow.
!
!                                                             2016.12.07 K.Noh
!
!------------------------------------------------------------------------------!

        SUBROUTINE SETUP

            USE LES_FILTERING_module,                                           &
                ONLY : Nx, Ny, Nz, dx, dz, FW, pi, tol,                         &
                       file_name, dir_name, path_name, VS_CASE, FILTER_OX,      &
                       VS_ONLY, Y_ORDER

            USE LES_FILTERING_module,                                           &
                ONLY : X, Y, Z, dy, U, V, W, U_Fil, V_Fil, W_Fil,               &
                       Resi_T, S_T, S_T_Fil, NU_R, O_T, O_T_Fil, VS,            &
                       L_T, M_T, U_Fil_2, V_Fil_2, W_Fil_2, S_T_Fil_2, Cs, YP

            IMPLICIT NONE
            INTEGER :: i,j,k,N

            pi = atan(1.0)*4

            !------------------------------------------------------------------!
            !                  Make & Initialize Result folder                 !
            !------------------------------------------------------------------!
            file_name = 'instantaneous_velocity_field_re644.plt'
            ! file_name = 'INSU_XYZ.plt'
            dir_name  = 'RESULT'

            CALL SYSTEM('mkdir '//TRIM(dir_name))
            CALL SYSTEM('mkdir '//TRIM(dir_name)//'/U')
            CALL SYSTEM('mkdir '//TRIM(dir_name)//'/RESIDUAL_STRESS')
            CALL SYSTEM('mkdir '//TRIM(dir_name)//'/STRAIN_RATE')
            CALL SYSTEM('mkdir '//TRIM(dir_name)//'/ROTATION_RATE')
            CALL SYSTEM('mkdir '//TRIM(dir_name)//'/EDDY_VISCOSITY')
            CALL SYSTEM('mkdir '//TRIM(dir_name)//'/VORTICAL_STRUCTURE')
            CALL SYSTEM('mkdir '//TRIM(dir_name)//'/SMARGORINSKY_COEFFICIENT')

            CALL SYSTEM('rm -rf ./'//TRIM(dir_name)//'/*.plt')
            CALL SYSTEM('rm -rf ./'//TRIM(dir_name)//'/U'//'/*.plt')
            CALL SYSTEM('rm -rf ./'//TRIM(dir_name)//'/RESI'//'/*.plt')
            CALL SYSTEM('rm -rf ./'//TRIM(dir_name)//'/STRAIN_RATE'//'/*.plt')
            CALL SYSTEM('rm -rf ./'//TRIM(dir_name)//'/ROTATION_RATE'//'/*.plt')
            CALL SYSTEM('rm -rf ./'//TRIM(dir_name)//'/EDDY_VISCOSITY'//'/*.plt')
            CALL SYSTEM('rm -rf ./'//TRIM(dir_name)//'/VORTICAL_STRUCTURE'//'/*.plt')
            CALL SYSTEM('rm -rf ./'//TRIM(dir_name)//'/SMARGORINSKY_COEFFICIENT'//'/*.plt')

            !------------------------------------------------------------------!
            !                          Statistic type                          !
            !                                                                  !
            !   (a) All Statistic             : 0                              !
            !   (b) Excpet Vortical Structure : 1                              !
            !   (c) Only Vortical Structure   : 2                              !
            !                                                                  !
            !------------------------------------------------------------------!
            VS_ONLY = 1

            !------------------------------------------------------------------!
            !                    Vortical Structure methods                    !
            !                                                                  !
            !   (a) Q-criteria         : 1                                     !
            !   (b) Lambda_2 criteria  : 2                                     !
            !   (c) Lambda_ci criteria : 3                                     !
            !                                                                  !
            !------------------------------------------------------------------!
            VS_CASE = 2

            !------------------------------------------------------------------!
            !                       The number of Filter                       !
            !                                                                  !
            !   (a) No-Filter     : 0                                          !
            !   (b) First Filter  : 1                                          !
            !   (c) Second Filter : 2                                          !
            !                                                                  !
            !------------------------------------------------------------------!
            FILTER_OX = 0

            !------------------------------------------------------------------!
            !                         Y sorting order                          !
            !                                                                  !
            !   (a) Ascending Order  : 0                                       !
            !   (b) Descending Order : 1                                       !
            !                                                                  !
            !------------------------------------------------------------------!
            Y_ORDER = 0

            !------------------------------------------------------------------!
            !                       Statistic Slice Point                      !
            !                                                                  !
            !    N  : Total Number of Slice Point                              !
            !    YP : Slice Point Position                                     !
            !                                                                  !
            !------------------------------------------------------------------!
            N = 3
            ALLOCATE( YP(1:N) )
            YP = [5,30,200]

            !------------------------------------------------------------------!
            !                    Constants for LES filtering                   !
            !------------------------------------------------------------------!
            Nx = 288
            Ny = 257
            Nz = 288

            FW = 4     ! Filter width constant
            tol = 1e-8 ! Tolerance for the number of nodes in x,z directions

            !------------------------------------------------------------------!
            !                         Allocate variables                       !
            !------------------------------------------------------------------!
            ALLOCATE( X(1:Nx),Y(1:Ny),Z(1:Nz),dy(1:Ny-1) )
            ALLOCATE( U(1:Nx,1:Ny,1:Nz), V(1:Nx,1:Ny,1:Nz), W(1:Nx,1:Ny,1:Nz) )
            ALLOCATE( U_Fil(1:Nx,1:Ny,1:Nz), V_Fil(1:Nx,1:Ny,1:Nz),             &
                      W_Fil(1:Nx,1:Ny,1:Nz), VS(1:Nx,1:Ny,1:Nz),                &
                      U_Fil_2(1:Nx,1:Ny,1:Nz), V_Fil_2(1:Nx,1:Ny,1:Nz),         &
                      W_Fil_2(1:Nx,1:Ny,1:Nz), Cs(1:Nx,1:Ny,1:Nz)             )
            ALLOCATE(Resi_T(1:Nx,1:Ny,1:Nz,1:3,1:3),NU_R(1:Nx,1:Ny,1:Nz,1:3,1:3))
            ALLOCATE(S_T(1:Nx,1:Ny,1:Nz,1:3,1:3),S_T_Fil(1:Nx,1:Ny,1:Nz,1:3,1:3))
            ALLOCATE(O_T(1:Nx,1:Ny,1:Nz,1:3,1:3),O_T_Fil(1:Nx,1:Ny,1:Nz,1:3,1:3))
            ALLOCATE(S_T_Fil_2(1:Nx,1:Ny,1:Nz,1:3,1:3))
            ALLOCATE(M_T(1:Nx,1:Ny,1:Nz,1:3,1:3),L_T(1:Nx,1:Ny,1:Nz,1:3,1:3))

            !------------------------------------------------------------------!
            !                         Initial Conditions                       !
            !------------------------------------------------------------------!
            X(1:Nx) = 0.0
            Y(1:Ny) = 0.0
            Z(1:Nz) = 0.0

            dy(1:Ny-1) = 0.0

            U(1:Nx,1:Ny,1:Nz) = 0.0
            V(1:Nx,1:Ny,1:Nz) = 0.0
            W(1:Nx,1:Ny,1:Nz) = 0.0

            VS(1:Nx,1:Ny,1:Nz)    = 0.0
            U_Fil(1:Nx,1:Ny,1:Nz) = 0.0
            V_Fil(1:Nx,1:Ny,1:Nz) = 0.0
            W_Fil(1:Nx,1:Ny,1:Nz) = 0.0

            Cs(1:Nx,1:Ny,1:Nz)      = 0.0
            U_Fil_2(1:Nx,1:Ny,1:Nz) = 0.0
            V_Fil_2(1:Nx,1:Ny,1:Nz) = 0.0
            W_Fil_2(1:Nx,1:Ny,1:Nz) = 0.0

            Resi_T(1:Nx,1:Ny,1:Nz,1:3,1:3)  = 0.0
            S_T(1:Nx,1:Ny,1:Nz,1:3,1:3)     = 0.0
            S_T_Fil(1:Nx,1:Ny,1:Nz,1:3,1:3) = 0.0
            O_T(1:Nx,1:Ny,1:Nz,1:3,1:3)     = 0.0
            O_T_Fil(1:Nx,1:Ny,1:Nz,1:3,1:3) = 0.0
            NU_R(1:Nx,1:Ny,1:Nz,1:3,1:3)    = 0.0

            S_T_Fil_2(1:Nx,1:Ny,1:Nz,1:3,1:3)     = 0.0
            L_T(1:Nx,1:Ny,1:Nz,1:3,1:3)           = 0.0
            M_T(1:Nx,1:Ny,1:Nz,1:3,1:3)           = 0.0

        END SUBROUTINE SETUP
