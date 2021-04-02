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
! *  Authors: Juha Ruokolainen, Fabien Gillet-Chaulet, Olivier Gagliardini,
! *  David Lilien
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 08 Jun 1997
! *  Date of modification: 13/10/05 from version 1.5
! *  Date of modification: 26/03/21 from version 1.5
! *
! *****************************************************************************
!> Module containing a solver for computing deformational heat
!> Differs from the canonical version since stress is calculated using the
!> anisotropic effective viscosity. This solver will calculate the values of the
!> deformational heat on the integration points, which is the advantage compared
!> to simply taking sigma:epsilon_dot e.g. via MATC
!------------------------------------------------------------------------------
RECURSIVE SUBROUTINE DeformationalHeatSolverAI( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
  
  USE DefUtils
  
  IMPLICIT NONE
  
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Solve stress equations for one timestep
!
!  ARGUMENTS:
!
!  TYPE(Model_t) :: Model,  
!     INPUT: All model information (mesh,materials,BCs,etc...)
!
!  TYPE(Solver_t) :: Solver
!     INPUT: Linear equation solver options
!
!  REAL(KIND=dp) :: dt,
!     INPUT: Timestep size for time dependent simulations (NOTE: Not used
!            currently)
!
!******************************************************************************
  
  TYPE(Model_t)  :: Model
  TYPE(Solver_t), TARGET :: Solver
  
  LOGICAL ::  TransientSimulation
  REAL(KIND=dp) :: dt
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
  TYPE(Solver_t), POINTER :: PSolver
  
  TYPE(Matrix_t),POINTER :: StiffMatrix
  
  INTEGER :: i, j, k, l, n, t, iter, NDeg
  INTEGER :: dim, STDOFs, LocalNodes, istat, spoofdim
  
  TYPE(ValueList_t),POINTER :: Material, BC, BodyForce
  TYPE(Nodes_t) :: ElementNodes
  TYPE(Element_t),POINTER :: CurrentElement
  
  REAL(KIND=dp) :: RelativeChange, UNorm=0.0, PrevUNorm=0.0
  
  REAL(KIND=dp), ALLOCATABLE :: Basis(:),ddBasisddx(:,:,:)
  REAL(KIND=dp), ALLOCATABLE :: dBasisdx(:,:)
  REAL(KIND=dp) :: u, v, w, detJ
  
  LOGICAL :: stat, CSymmetry 
  
  INTEGER :: NewtonIter, NonlinearIter
  
  TYPE(Variable_t), POINTER :: FlowVariable 
  
  REAL(KIND=dp), POINTER ::  Solution(:), &
       ForceVector(:), FlowValues(:) 
  
  INTEGER, POINTER :: NodeIndexes(:), &
       FlowPerm(:)
  
  INTEGER :: body_id
  INTEGER :: old_body = -1
  
  LOGICAL :: Isotropic, AllocationsDone = .FALSE.,  &
       Requal0
  LOGICAL :: GotIt, UnFoundFatal=.TRUE.,OutOfPlaneFlow
  
  REAL(KIND=dp), ALLOCATABLE:: LocalMassMatrix(:,:), &
       LocalStiffMatrix(:,:), LocalForce(:), &
       LocalP(:),  &
       LocalVelo(:,:)
  
  INTEGER :: NumberOfBoundaryNodes
  INTEGER, POINTER :: BoundaryReorder(:)
  
  REAL(KIND=dp), POINTER :: BoundaryNormals(:,:), &
       BoundaryTangent1(:,:), BoundaryTangent2(:,:)
  CHARACTER(LEN=MAX_NAME_LEN) :: FlowSolverName

  ! Variables needed for material defs
  TYPE(ValueList_t),POINTER :: SolverParams
  CHARACTER(LEN=MAX_NAME_LEN) :: viscosityFile, TempVar
  REAL(KIND=dp) :: FabricGrid(4878), Wn(7), MinSRInvariant
  TYPE(Variable_t), POINTER :: FabricVariable, TempSol
  REAL(KIND=dp), POINTER :: FabricValues(:), Temperature(:)
  INTEGER, POINTER :: FabricPerm(:), TempPerm(:)
  REAL(KIND=dp), ALLOCATABLE:: LocalFluidity(:), &
       LocalTemperature(:), K1(:), K2(:), E1(:), &
       E2(:), E3(:)

  ! And variables to spoof the dimension
  CHARACTER(LEN=MAX_NAME_LEN) :: OOPlaneRotVar13, OOPLaneRotVar23
  REAL(KIND=dp), ALLOCATABLE :: LocalOOP13(:), LocalOOP23(:)
  LOGICAL :: OOPlaneRot13, OOPlaneRot23
  TYPE(Variable_t), POINTER :: OOPlaneRotSol13, OOPlaneRotSol23
  REAL(KIND=dp), POINTER :: OOPlaneRotValues13(:), OOPlaneRotValues23(:)


  INTEGER, POINTER :: OOPlaneRotPerm13(:), OOPlaneRotPerm23(:)
  
#ifdef USE_ISO_C_BINDINGS
  REAL(KIND=dp) :: at, at0
#else
  REAL(KIND=dp) :: at, at0, CPUTime, RealTime
#endif
!------------------------------------------------------------------------------
  SAVE NumberOfBoundaryNodes, BoundaryReorder, BoundaryNormals, &
       BoundaryTangent1, BoundaryTangent2, ViscosityFile, Wn, FabricGrid, &
       MinSRInvariant, LocalTemperature, K1, K2, E1, E2, E3, N, UnFoundFatal, &
       LocalFluidity

  SAVE Spoofdim, LocalOOP13, LocalOOP23
  
  SAVE Basis, dBasisdx, ddBasisddx
  SAVE LocalMassMatrix, LocalStiffMatrix, LocalForce, &
       ElementNodes,  &
       AllocationsDone,  &
       old_body, Unorm, PrevUnorm
  
  SAVE LocalVelo, LocalP, dim


      Wn(7) = GetConstReal( Model % Constants, 'Gas Constant', GotIt )
      IF (.NOT.GotIt) THEN
        WRITE(Message,'(A)') 'VariableGas Constant  not found. &
                     &Setting to 8.314'
        CALL INFO('DeformationalHeatAI', Message, level=20)
        Wn(7) = 8.314
      ELSE
        WRITE(Message,'(A,F10.4)') 'Gas Constant = ',   Wn(7)
        CALL INFO('DeformationalHeatAI', Message , level = 20)
      END IF

!------------------------------------------------------------------------------
!  Read the name of the Flow Solver (NS)
!------------------------------------------------------------------------------
  FlowSolverName = GetString( Solver % Values, 'Flow Solver Name', GotIt )    
  IF (.NOT.Gotit) FlowSolverName = 'AIFlow'
  FlowVariable => VariableGet( Solver % Mesh % Variables, FlowSolverName,&
       UnFoundFatal=UnFoundFatal )
  FlowPerm    => FlowVariable % Perm
  FlowValues  => FlowVariable % Values
  OutOfPlaneFlow = .FALSE.

  SolverParams => GetSolverParams()
  TempVar = ListGetString( SolverParams,'Temperature Solution Name',GotIt,UnFoundFatal )
  IF (.NOT.GotIt) THEN
      TempVar = 'Temperature'
  END IF
  TempSol => VariableGet( Solver % Mesh % Variables, TempVar )
  IF ( ASSOCIATED( TempSol) ) THEN
    TempPerm    => TempSol % Perm
    Temperature => TempSol % Values
  END IF
  WRITE(Message,'(A,A)') 'Temperature variable = ', TempVar
  CALL INFO('DeformationalHeatAI', Message , level = 20)

  OOPlaneRotVar23 = ListGetString( SolverParams,'OOPlane23 Strain Name',GotIt,UnFoundFatal=.FALSE.)
  IF (.NOT.GotIt) THEN
      OOPlaneRot23 = .FALSE.
      OOPlaneRotVar23 = 'Dummy'
  ELSE
      OOPlaneRot23 = .TRUE.
      OOPlaneRotSol23 => VariableGet( Solver % Mesh % Variables, OOPLaneRotVar23 )
      IF ( ASSOCIATED( OOPlaneRotSol23) ) THEN
        OOPlaneRotPerm23 => OOPlaneRotSol23 % Perm
        OOPlaneRotValues23 => OOPlaneRotSol23 % Values
      END IF
  END IF
  WRITE(Message,'(A,A)') 'OOPlane23 variable = ', OOPlaneRotVar23
  CALL INFO('DeformationalHeatAI', Message , level = 20)

  OOPlaneRotVar13 = ListGetString( SolverParams,'OOPlane13 Strain Name',GotIt,UnFoundFatal=.FALSE.)
  IF (.NOT.GotIt) THEN
      OOPlaneRot13 = .FALSE.
      OOPlaneRotVar13 = 'Dummy'
  ELSE
      OOPlaneRot13 = .TRUE.
      OOPlaneRotSol13 => VariableGet( Solver % Mesh % Variables, OOPLaneRotVar13 )
      IF ( ASSOCIATED( OOPlaneRotSol13) ) THEN
        OOPlaneRotPerm13 => OOPlaneRotSol13 % Perm
        OOPlaneRotValues13 => OOPlaneRotSol13 % Values
      END IF
  END IF
!------------------------------------------------------------------------------
!  Read constants from constants section of SIF file
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
  IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN
  
  Solution => Solver % Variable % Values
  STDOFs   =  Solver % Variable % DOFs
  
  IF ( STDOFs /=1 ) THEN
     CALL Fatal( 'DeformationalHeatAI', 'DOF must be equal to 1' )
  END IF
  
  StiffMatrix => Solver % Matrix
  ForceVector => StiffMatrix % RHS
!------------------------------------------------------------------------------
!     Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
  IF ( .NOT. AllocationsDone .OR. Solver % MeshChanged) THEN
     N = Model % MaxElementNodes
     dim = CoordinateSystemDimension()
      IF ( OOPlaneRot13 .OR. OOPlaneRot23 ) THEN
        spoofdim = dim + 1
      ELSE
        spoofdim = dim
      END IF
     IF ( AllocationsDone ) THEN
        DEALLOCATE( ElementNodes % x,     &
             ElementNodes % y,     &
             ElementNodes % z,     &
             LocalVelo, LocalP,    &                      
             Basis, ddBasisddx,    &
             dBasisdx,             &
             LocalMassMatrix,      &
             LocalStiffMatrix,     &
             LocalForce,           &
             K1, K2, E1, E2, E3,   &
             LocalTemperature,     &
             LocalOOP23,LocalOOP13,&
             LocalFluidity )
     END IF

     ALLOCATE( ElementNodes % x( N ), &
          ElementNodes % y( N ), &
          ElementNodes % z( N ), &
          LocalVelo( 3,N ), LocalP( N ), &          
          Basis( 2*N ),ddBasisddx(1,1,1), dBasisdx( 2*N,3 ), &
          LocalMassMatrix( 2*STDOFs*N,2*STDOFs*N ),  &
          LocalStiffMatrix( 2*STDOFs*N,2*STDOFs*N ),  &
          LocalForce( 2*STDOFs*N ),  &
          K1( N ), K2( N ), E1( N ), E2( N ), E3( N ), &
          LocalTemperature(N), &
          LocalOOP13(N), LocalOOP23(N), &
          LocalFluidity(N), STAT=istat )

     IF ( istat /= 0 ) THEN
        CALL Fatal( 'DeformationalHeatAI', 'Memory allocation error.' )
     END IF
!------------------------------------------------------------------------------

     AllocationsDone = .TRUE.
  END IF

  FabricVariable => VariableGet(Solver % Mesh % Variables, 'Fabric')
  IF ( ASSOCIATED( FabricVariable ) ) THEN
    FabricPerm    => FabricVariable % Perm
    FabricValues => FabricVariable % Values
  END IF

!------------------------------------------------------------------------------
  NonlinearIter = 1
  DO iter=1,NonlinearIter

     at  = CPUTime()
     at0 = RealTime()

     CALL Info( 'DeformationalHeatAI', ' ', Level=4 )
     CALL Info( 'DeformationalHeatAI', ' ', Level=4 )
     CALL Info( 'DeformationalHeatAI', ' ', Level=4 )
     CALL Info( 'DeformationalHeatAI', ' ', Level=4 )
     CALL Info( 'DeformationalHeatAI', 'Starting assembly...',Level=4 )

!------------------------------------------------------------------------------
        CALL DefaultInitialize()
!------------------------------------------------------------------------------

        DO t=1,Solver % NumberOFActiveElements
           IF ( RealTime() - at0 > 1.0 ) THEN
              WRITE(Message,'(a,i3,a)' ) '   Assembly: ',  &
                   INT(100.0 - 100.0 * (Solver % NumberOfActiveElements-t) / &
                   (1.0*Solver % NumberOfActiveElements)), ' % done'
              CALL Info( 'DeformationalHeatAI', Message, Level=5 )
              at0 = RealTime()
           END IF

           CurrentElement => GetActiveElement(t)
           n = GetElementNOFNodes()
           NodeIndexes => CurrentElement % NodeIndexes

           ElementNodes % x(1:n) = Model % Nodes % x(NodeIndexes(1:n))
           ElementNodes % y(1:n) = Model % Nodes % y(NodeIndexes(1:n))
           ElementNodes % z(1:n) = Model % Nodes % z(NodeIndexes(1:n))

           Material => GetMaterial()
           body_id = CurrentElement % BodyId
           IF (body_id /= old_body) Then 
             old_body = body_id
             Call  GetMaterialDefs()
           END IF
           LocalFluidity(1:n) = ListGetReal( Material, &
                            'Fluidity Parameter', n, NodeIndexes, GotIt,&
                            UnFoundFatal=UnFoundFatal)

            LocalTemperature = 0.0D0
            IF ( ASSOCIATED(TempSol) ) THEN
              DO i=1,n
               k = TempPerm(NodeIndexes(i))
                 LocalTemperature(i) = Temperature(k)
              END DO
            ELSE
               LocalTemperature(1:n) = 0.0d0
            END IF
!------------------------------------------------------------------------------
!    Read in material constants from Material section
!------------------------------------------------------------------------------


!!!! Restricted to the Power Law case

           IF(.NOT.Isotropic) Then
             K1(1:n) = FabricValues( 5 * (FabricPerm(NodeIndexes(1:n))-1) + 1 ) 
             K2(1:n) = FabricValues( 5 * (FabricPerm(NodeIndexes(1:n))-1) + 2 )
             E1(1:n) = FabricValues( 5 * (FabricPerm(NodeIndexes(1:n))-1) + 3 )
             E2(1:n) = FabricValues( 5 * (FabricPerm(NodeIndexes(1:n))-1) + 4 )
             E3(1:n) = FabricValues( 5 * (FabricPerm(NodeIndexes(1:n))-1) + 5 )
           END IF

           LocalVelo = 0.0_dp
           DO i=1, dim
             LocalVelo(i,1:n) = FlowValues((dim+1)*(FlowPerm(NodeIndexes(1:n))-1) + i)
           END DO
           BodyForce => GetBodyForce()
           
           LocalP(1:n) = FlowValues((dim+1)*FlowPerm(NodeIndexes(1:n)))

          If ( OOPlaneRot13 ) THEN
             LocalOOP13 = 0.0_dp
             LocalOOP13(1:n) = OOPlaneRotValues13(OOPlaneRotPerm13(NodeIndexes))
          END IF
          If ( OOPlaneRot23 ) THEN
             LocalOOP23 = 0.0_dp
             LocalOOP23(1:n) = OOPlaneRotValues23(OOPlaneRotPerm23(NodeIndexes))
          END IF

           CALL LocalAIMatrix(LocalMassMatrix, LocalStiffMatrix, LocalForce, &
                K1, K2, E1, E2, E3, Wn, MinSRInvariant, LocalVelo, LocalP, &
                LocalFluidity, LocalTemperature, CurrentElement, n, &
                ElementNodes, Isotropic, OOPlaneRot13, OOPlaneRot23, LocalOOP13, LocalOOP13)

!------------------------------------------------------------------------------
!        Update global matrices from local matrices 
!------------------------------------------------------------------------------
           CALL DefaultUpdateEquations( LocalStiffMatrix, LocalForce )

        END DO

        CALL Info( 'DeformationalHeatAI', 'Assembly done', Level=4 )


        CALL DefaultFinishAssembly()

!------------------------------------------------------------------------------
!     Dirichlet boundary conditions
!------------------------------------------------------------------------------
        CALL DefaultDirichletBCs()

!------------------------------------------------------------------------------

        CALL Info( 'DeformationalHeatAI', 'Set boundaries done', Level=4 )

!------------------------------------------------------------------------------
!     Solve the system and check for convergence
!------------------------------------------------------------------------------
        PrevUNorm = UNorm

        UNorm = DefaultSolve()

     Solver % Variable % Norm = Unorm  

     IF ( PrevUNorm + UNorm /= 0.0d0 ) THEN
        RelativeChange = 2.0d0 * ABS( PrevUNorm - UNorm) / ( PrevUnorm + UNorm)
     ELSE
        RelativeChange = 0.0d0
     END IF

     WRITE( Message, * ) 'Result Norm   : ',UNorm, PrevUNorm
     CALL Info( 'DeformationalHeatAI', Message, Level=4 )
     WRITE( Message, * ) 'Relative Change : ',RelativeChange
     CALL Info( 'DeformationalHeatAI', Message, Level=4 )


!------------------------------------------------------------------------------
  END DO ! of nonlinear iter
!------------------------------------------------------------------------------


CONTAINS


    ! Copied this one from AIflowSolve_nlS2
SUBROUTINE GetMaterialDefs()
      ! check if we are isotropic or not
      Isotropic = ListGetLogical( Material , 'Isotropic',Gotit )
      IF (.NOT.Gotit) Then
          Isotropic = .False.
           WRITE(Message,'(A)') 'Isotropic set to False'
           CALL INFO('AIFlowSolve', Message, Level = 20)
      ELSE
           IF ( (ASSOCIATED( FabricVariable )).AND.Isotropic ) Then
              WRITE(Message,'(A)') 'Be careful Isotropic is true &
                           & and Fabric is defined!'
              CALL INFO('AIFlowSolve', Message, Level = 1)
           END IF
      END IF

      IF (.NOT.Isotropic) Then
        ! Get the viscosity file and store the viscosities into FabricGrid
         viscosityFile = ListGetString( Material ,'Viscosity File',GotIt,UnFoundFatal )
         OPEN( 1, File = viscosityFile)
         DO i=1,813
             READ( 1, '(6(e14.8))' ) FabricGrid( 6*(i-1)+1:6*(i-1)+6 )
         END DO
         CLOSE(1)
      ENDIF

      Wn(2) = ListGetConstReal( Material , 'Powerlaw Exponent', GotIt,UnFoundFatal=UnFoundFatal)
         !Previous default value: Wn(2) = 1.0
      WRITE(Message,'(A,F10.4)') 'Powerlaw Exponent = ',   Wn(2)
      CALL INFO('AIFlowSolve', Message, Level = 20)

      Wn(3) = ListGetConstReal( Material, 'Activation Energy 1', GotIt,UnFoundFatal=UnFoundFatal)
         !Previous default value: Wn(3) = 1.0
      WRITE(Message,'(A,F10.4)') 'Activation Energy 1 = ',   Wn(3)
      CALL INFO('AIFlowSolve', Message, Level = 20)

      Wn(4) = ListGetConstReal( Material, 'Activation Energy 2', GotIt,UnFoundFatal=UnFoundFatal)
         !Previous default value: Wn(4) = 1.0
      WRITE(Message,'(A,F10.4)') 'Activation Energy 2 = ',   Wn(4)
      CALL INFO('AIFlowSolve', Message, Level = 20)

      Wn(5) = ListGetConstReal(Material, 'Reference Temperature', GotIt,UnFoundFatal=UnFoundFatal)
         !Previous default value: Wn(5) = -10.0 (Celsius)
      WRITE(Message,'(A,F10.4)') 'Reference Temperature = ',   Wn(5)
      CALL INFO('AIFlowSolve', Message, Level = 20)

      Wn(6) = ListGetConstReal( Material, 'Limit Temperature', GotIt,UnFoundFatal=UnFoundFatal)
         !Previous default value: Wn(6) = -10.0 (Celsius)
      WRITE(Message,'(A,F10.4)') 'Limit Temperature = ',   Wn(6)
      CALL INFO('AIFlowSolve', Message, Level = 20)

! Get the Minimum value of the Effective Strain rate 
      MinSRInvariant = 100.0*AEPS

      IF ( Wn(2) > 1.0  ) THEN
        MinSRInvariant =  &
             ListGetConstReal( Material, 'Min Second Invariant', GotIt )
        IF (.NOT.GotIt) THEN
          WRITE(Message,'(A)') 'Variable Min Second Invariant not &
                    &found. Setting to 100.0*AEPS )'
          CALL INFO('AIFlowSolve', Message, Level = 20)
        ELSE
          WRITE(Message,'(A,E14.8)') 'Min Second Invariant = ', MinSRInvariant
          CALL INFO('AIFlowSolve', Message, Level = 20)
        END IF
      END IF

!------------------------------------------------------------------------------
      END SUBROUTINE GetMaterialDefs

!------------------------------------------------------------------------------
  SUBROUTINE LocalAIMatrix(MassMatrix, StiffMatrix, ForceVector, &
       NodalK1, NodalK2, NodalE1, NodalE2, NodalE3, Wn, MinSRInvariant, &
       NodalVelo, NodalP, NodalFluidity, &
       NodalTemp, Element, n, Nodes, isotropic, OOPlaneRot13, OOPlaneRot23, &
       NodalOOP13, NodalOOP23 )
!------------------------------------------------------------------------------
    
    USE MaterialModels
    
    REAL(KIND=dp) :: StiffMatrix(:,:), MassMatrix(:,:)
    REAL(KIND=dp) ::  NodalVelo(:,:)
    REAL(KIND=dp), DIMENSION(:) :: ForceVector, NodalP
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element
    INTEGER :: n

    REAL(KIND=dp) :: MinSRInvariant, Wn(7)
    REAL(KIND=dp) :: NodalK1(:), NodalK2(:)
    REAL(KIND=dp) :: NodalE1(:), NodalE2(:), NodalE3(:)
    REAL(KIND=dp) :: NodalFluidity(:), NodalTemp(:)

    REAL(KIND=dp) :: NodalOOP13(:), NodalOOP23(:)
    LOGICAL :: OOPlaneRot13, OOPlaneRot23
!------------------------------------------------------------------------------
!
    REAL(KIND=dp) :: Basis(2*n),ddBasisddx(1,1,1)
    REAL(KIND=dp) :: dBasisdx(2*n,3),detJ, pBasis(n)
    
    REAL(KIND=dp) :: Stress, epsi
    
    REAL(KIND=dp) :: Pressure
    REAL(KIND=dp) :: LGrad(3,3), StrainRate(3,3), StressTensor(3,3)
    
    INTEGER :: i, j, k, p, q, t, dim, NBasis,  LinearBasis
    
    REAL(KIND=dp) :: s, u, v, w, eta
    
    REAL(KIND=dp) :: dDispldx(3,3)
    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
    
    INTEGER, POINTER :: EdgeMap(:,:)
    INTEGER :: N_Integ, nd
    INTEGER, DIMENSION(6), PARAMETER :: indx = (/1, 2, 3, 1, 2, 3/), &
                                         indy = (/1, 2, 3, 2, 3, 1/)
    
    REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ
    
    LOGICAL :: stat, CSymmetry

    REAL(KIND=dp) :: C(6,6)
    REAL(KIND=dp) :: Temp, ai(3), Angle(3), a2(6), D(6)
    REAL(KIND=dp) :: Bg, BGlenT, ss, nn
    LOGICAL :: Isotropic

    INTEGER :: spoofdim

     INTERFACE
      Subroutine R2Ro(a2,dim,spoofdim,ai,angle)
         USE Types
         REAL(KIND=dp),intent(in) :: a2(6)
         Integer :: dim,spoofdim
         REAL(KIND=dp),intent(out) :: ai(3), Angle(3)
      End Subroutine R2Ro
                 
      Subroutine OPILGGE_ai_nl(ai,Angle,etaI,eta36)
          USE Types
          REAL(kind=dp), INTENT(in),  DIMENSION(3)   :: ai
          REAL(kind=dp), INTENT(in),  DIMENSION(3)   :: Angle
          REAL(kind=dp), INTENT(in),  DIMENSION(:)   :: etaI
          REAL(kind=dp), INTENT(out), DIMENSION(6,6) :: eta36
        END SUBROUTINE OPILGGE_ai_nl
      END INTERFACE
    
!------------------------------------------------------------------------------
    dim = CoordinateSystemDimension()
    
    IF ((OOPlaneRot13).OR.( OOPlaneRot23 )) THEN
          SpoofDim = dim + 1
    ELSE
          SpoofDim = dim
    END IF

    ForceVector = 0.0_dp
    StiffMatrix = 0.0_dp
    MassMatrix  = 0.0_dp
    
    IntegStuff = GaussPoints( Element )
     
    U_Integ => IntegStuff % u
    V_Integ => IntegStuff % v
    W_Integ => IntegStuff % w
    S_Integ => IntegStuff % s
    N_Integ =  IntegStuff % n
!
!   Now we start integrating
!
    DO t=1,N_Integ

       u = U_Integ(t)
       v = V_Integ(t)
       w = W_Integ(t)

!------------------------------------------------------------------------------
!     Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
       stat = ElementInfo(Element,Nodes,u,v,w,detJ, &
            Basis,dBasisdx,ddBasisddx,.FALSE.,.FALSE.)
       
       s = detJ * S_Integ(t)
       
       LGrad = MATMUL( NodalVelo(:,1:n), dBasisdx(1:n,:) )
       StrainRate = 0.5 * ( LGrad + TRANSPOSE(LGrad) )
      ! If we are adding in manual spins, do it here before the stress
      ! calculations
      IF ( OOPlaneRot23 ) THEN
          StrainRate(2, 3) = SUM( NodalOOP23(1:n) * Basis(1:n) )
          StrainRate(3, 2) = StrainRate(2, 3)
      END IF
      IF ( OOPlaneRot13 ) THEN
          StrainRate(3, 1) = SUM( NodalOOP13(1:n) * Basis(1:n) )
          StrainRate(1, 3) = StrainRate(3, 1)
      END IF
       StressTensor = 0.0_dp

       Temp = SUM( NodalTemp(1:n)*Basis(1:n) )
       Wn(1) = SUM( NodalFluidity(1:n)*Basis(1:n) )

!    ----------------------------
      IF (.Not.Isotropic) then
	    C = 0.0_dp
!    Material parameters at that point
        a2(1) = SUM( NodalK1(1:n) * Basis(1:n) ) 
        a2(2) = SUM( NodalK2(1:n) * Basis(1:n) ) 
        a2(3) = 1.d0 - a2(1) - a2(2)
        a2(4) = SUM( NodalE1(1:n) * Basis(1:n) )
        a2(5) = SUM( NodalE2(1:n) * Basis(1:n) )
        a2(6) = SUM( NodalE3(1:n) * Basis(1:n) )
      
        CALL R2Ro(a2,dim,spoofdim,ai,Angle)
        CALL OPILGGE_ai_nl(ai,Angle,FabricGrid,C)
         
!
!    Compute deviatoric stresses: 
!    ----------------------------
      D(1) = StrainRate(1,1)
      D(2) = StrainRate(2,2)
      D(3) = StrainRate(3,3)
      D(4) = 2. * StrainRate(1,2)
      D(5) = 2. * StrainRate(2,3)
      D(6) = 2. * StrainRate(3,1)
      
      DO k = 1, 2*spoofdim
        DO j = 1, 2*spoofdim
          StressTensor( INDx(k),INDy(k) ) = &
          StressTensor( INDx(k),INDy(k) ) + C(k,j) * D(j)
        END DO
        IF (k > 3)  StressTensor( INDy(k),INDx(k) ) = StressTensor( INDx(k),INDy(k) )
      END DO
	ELSE  ! ISOTROPIC CASE
	  StressTensor=2._dp * StrainRate
	END IF
	   

! non relative viscosities
    ! Glen fluidity       
	 Bg=BGlenT(Temp,Wn)
	 ss=1.0_dp
     ! Case Non linear
	 IF (Wn(2) > 1.0) THEN 
	   Bg=Bg**(1.0/Wn(2))
       ss = 0.0_dp
       DO i = 1, 3
         DO j = 1, 3
           ss = ss + StressTensor(i,j)**2
         END DO
       END DO
       nn = (1.0 - Wn(2))/(2.0*Wn(2))
       ss = (ss / 2.0)**nn
       IF (ss < MinSRInvariant ) ss = MinSRInvariant
    END IF
       
	 StressTensor=StressTensor*ss/Bg

     ! Pressure = SUM( NodalP(1:n)*Basis(1:n) )
     ! Stress = Stress - Pressure
       
       DO p=1,n         
          DO q=1,n        
             StiffMatrix(p,q) =  &
                  StiffMatrix(p,q) + s*Basis(q)*Basis(p)
          END DO
          DO i=1,3
            DO j=1,3
              ForceVector(p) =  &
               ForceVector(p) + s*StressTensor(i, j)*StrainRate(i,j)*Basis(p) 
           END DO
         END DO
       END DO
     END DO
       
!------------------------------------------------------------------------------
  END SUBROUTINE LocalAIMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
END SUBROUTINE DeformationalHeatSolverAI
!------------------------------------------------------------------------------
