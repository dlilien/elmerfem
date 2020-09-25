!/*****************************************************************************/
! *
! *  Elmer/Ice, a glaciological add-on to Elmer
! *  http://elmerice.elmerfem.org
! *
! * 
! *  This program is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) any later version.
! * 
! *  This program is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *  GNU General Public License for more details.
! *
! *  You should have received a copy of the GNU General Public License
! *  along with this program (in file fem/GPL-2); if not, write to the 
! *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
! *  Boston, MA 02110-1301, USA.
! *
! *****************************************************************************/
! ******************************************************************************
! *
! *  Authors:  Nicholas Rathmann and David Lilien
! *  Email:   nicholas.rathmann@nbi.ku.dk and dlilien90@gmail.com
! *  Web:     http://elmerice.elmerfem.org
! *  Address: Niels Bohr Institute, University of Copenhagen
! *           Tagensvej 16
! *           Copenhagen 2200, Denmark
! *
! *       Date of modification: 09/20
! *
! *****************************************************************************/
      PROGRAM TestSpectralModel
        IMPLICIT NONE
        INTEGER, PARAMETER :: N_TSTEPS=200,SpectralDim=6,SpectralOrder=3
        INTEGER :: i, j, k, INFO
        INTEGER :: IPIV(SpectralDim)
        REAL(KIND=8), PARAMETER :: TSTEP=0.1
        REAL(KIND=8) :: C(SpectralDim), DCDt(SpectralDim,SpectralDim)
        REAL(KIND=8) :: Stiffness(SpectralDim, SpectralDim), Force(SpectralDim)
        REAL(KIND=8) :: StrainRate(3,3), Spin(3,3)
        REAL(KIND=8) :: OverlapMatrix
        CHARACTER(len=64) :: Message
        EXTERNAL         DGESV
        !!!!
        ! Test that the code works for a single point. Imposed strainrate for pure shear,
        ! no spin, implicit timestepping.
        ! Relies on Lapack for the matrix solve, with no assumptions about
        ! Symmetry, Positive definiteness, etc.
        ! Initially, this is just a very basic version where we call the subroutines
        !!!!!

        ! Vertical compression
        StrainRate = 0.0
        StrainRate(1,1) = 0.0005
        StrainRate(2,2) = 0.0005
        StrainRate(3,3) = -0.001
        Spin = 0.0

        ! Initialize the fabric
        C = 0.0
        C(1) = 1.0

        WRITE(*,'(A,F8.4,A,F7.4,F7.4,F7.4,F7.4,A)') 'At time: ', 0.0, ' C is ', C(1), &
               C(2), C(3), C(4), '...'
        DO i = 1,N_TSTEPS
          Stiffness = 0.0

          ! Implicit scheme, so diagonal starts as 1
          DO j = 1,SpectralDim
            Stiffness(j, j) = 1.0
            Force(j) = C(j)
          END DO
          CALL SpectralModel(SpectralDim, SpectralOrder, C, StrainRate, Spin, &
                             OverlapMatrix, DCDt)
          ! Be careful of sign errors in this...
          DO j = 1,SpectralDim
            DO k = 1,SpectralDim
              Stiffness(j, k) = Stiffness(j, k) - TSTEP * DCDt(j, k)
            END DO
          END DO

          CALL DGESV( SpectralDim, 1, Stiffness, SpectralDim, IPIV, FORCE, SpectralDim, INFO )

          DO j = 1,SpectralDim
            C(j) = Force(j)
          END DO
          IF (MODULO(i, 10).EQ.0) THEN
            WRITE(*,'(A,F8.4,A,F7.4,F7.4,F7.4,F7.4,A)') 'At time: ', i * TSTEP, ' C is ', C(1), &
                  C(2), C(3), C(4), '...'
          END IF
          CALL PostProcessFabric(C, SpectralDim)
          IF (MODULO(i, 10).EQ.0) THEN
            WRITE(*,'(A,F7.4,F7.4,F7.4,F7.4,A)') 'After PostProcess C is ', C(1), &
                 C(2), C(3), C(4), '...'
          END IF
        END DO
      END PROGRAM TestSpectralModel


      SUBROUTINE SpectralModel(ProbDim, SpectralOrder, C, StrainRate, Spin, &
                               OverlapMatrix, DCDt)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: SpectralOrder, ProbDim
        ! SpectralOrder is the maximum L, ProbDim is the resultant length of C
        REAL(KIND=8), INTENT(IN) :: C(ProbDim), StrainRate(3, 3), Spin(3,3)
        ! I have not assumed anything about the ordering of C. Define as you wish,
        ! then maybe document?
        ! I am just passing the StrainRate and Spin because the calculation of these
        ! quantities from the velocity is a bit annoying with basis functions, etc.
        REAL(KIND=8), INTENT(IN) :: OverlapMatrix
        ! Need to relate the size of OverlapMatrix
        ! to the SpectralOrder. You can size this however you like,
        ! and I will modify the main code to match
        REAL(KIND=8), INTENT(OUT) :: DCDt(ProbDim, ProbDim)
        ! I have the advection piece handled, so all that needs to be done is to deal
        ! with the Lagrangian fabric evolution. Elmer does the timestepping, so you just
        ! need to return a matrix such that the matrix DCDt * C = DC/Dt in a Lagrangian sense.
            
        INTEGER :: i
        REAL :: MeanC
        ! Set to dummy values so that I can check if the overall code works
        ! Here, just diffusion amongst the coefficients since it is easy
        DCDt = 0.0
        DO i = 1,ProbDim - 1
          DCDt(i, i) = -2
          DCDt(i, i + 1) = 1
          DCDt(i + 1, i) = 1
        END DO
        DCDt(1,1) = -1.0
        DCDt(ProbDim, ProbDim) = -1.0
      END SUBROUTINE SpectralModel


      SUBROUTINE PostProcessFabric(C, ProbDim)
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: ProbDim
        REAL(kind=8), INTENT(INOUT):: C(ProbDim)
        ! Use this function for anything that has to happen after the matrix solve.
        ! Modify C in place. This is just a dummy loop for now so that I have something
        ! To test with
        C(1) = C(1) - 0.01
        c(3) = C(3) + 0.01
      END SUBROUTINE
