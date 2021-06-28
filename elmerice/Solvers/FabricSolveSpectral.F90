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
! *  Authors:  David Lilien modified from standard fabric by
! *                Juha Ruokolainen, Olivier Gagliardini, Fabien Gillet-Chaulet
! *  Email:   david.lilien@nbi.ku.dk and Juha.Ruokolainen@csc.fi
! *  Web:     http://elmerice.elmerfem.org
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *       Date of modification: 06/21
! *
! *****************************************************************************/
!>  Solver for fabric parameter equations 
!------------------------------------------------------------------------------
      RECURSIVE SUBROUTINE FabricSolverSpectral( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------
      USE SpecFab
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

     INTEGER :: dim,n1,n2,i,j,k,l,n,t,iter,NDeg,STDOFs,LocalNodes,istat,spoofdim

     TYPE(ValueList_t),POINTER :: Material, BC, SolverParams
     TYPE(Nodes_t) :: ElementNodes
     TYPE(Element_t),POINTER :: CurrentElement, Element, &
              ParentElement, LeftParent, RightParent, Edge

     REAL(KIND=dp) :: RelativeChange,UNorm,PrevUNorm,Gravity(3), &
         Tdiff,Normal(3),NewtonTol,NonlinearTol,s,Wn(12)


     INTEGER :: NewtonIter,NonlinearIter

     TYPE(Variable_t), POINTER :: FabricSol, TempSol, FabricVariable, FlowVariable, &
                                  MeshVeloVariable,TensorFabricVariable,&
                                  OOPlaneRotSol13,OOPlaneRotSol23,GradSol

     REAL(KIND=dp), POINTER :: Temperature(:),Fabric(:), &
           FabricValues(:), FlowValues(:),TensorFabricValues(:),&
           MeshVeloValues(:), Solution(:), Ref(:), &
           OOPlaneRotValues13(:), OOPlaneRotValues23(:)


     INTEGER, POINTER :: TempPerm(:),FabricPerm(:),NodeIndexes(:), &
                        FlowPerm(:),MeshVeloPerm(:),TensorFabricPerm(:),&
                        OOPlaneRotPerm13(:),OOPlaneRotPerm23(:),GradPerm(:)

     REAL(KIND=dp) :: rho,lambda0,gamma0   !Interaction parameter,diffusion parameter
     REAL(KIND=dp) :: A1plusA2
     Real(KIND=dp), parameter :: Rad2deg=180._dp/Pi
     REAL(KIND=dp) :: a2(6), a2short(5)
     REAL(KIND=dp) :: ai(3), Angle(3)

     LOGICAL :: GotForceBC,GotIt,NewtonLinearization = .FALSE.,UnFoundFatal=.TRUE.
     LOGICAL :: OOPlaneRot13
     LOGICAL :: OOPlaneRot23

     INTEGER :: body_id,bf_id,eq_id, comp, Indexes(128)
!
     INTEGER :: old_body = -1

     REAL(KIND=dp) :: FabricGrid(4879)
                        
     LOGICAL :: AllocationsDone = .FALSE., FirstTime = .TRUE., FreeSurface

     TYPE(Variable_t), POINTER :: TimeVar

     REAL(KIND=dp), ALLOCATABLE:: MASS(:,:), STIFF(:,:),  &
       LocalFluidity(:), LOAD(:,:),Force(:), LocalTemperature(:), &
       Alpha(:,:),Beta(:), &
       Velocity(:,:), MeshVelocity(:,:), LocalOOP13(:), LocalOOP23(:), &
       LocalA4(:,:),ElGradVals(:,:,:),LocalGrad(:, :),LocalFabric(:,:),&
       LocalLHS(:, :), ElLHSVals(:, :, :)
     COMPLEX(KIND=dp), ALLOCATABLE:: nlm(:)

     SAVE MASS, STIFF, LOAD, Force,ElementNodes,Alpha,Beta, & 
          LocalTemperature, LocalFluidity,  AllocationsDone, &
          Wn,  FabricGrid, rho, lambda0, Velocity, &
          MeshVelocity, old_body, dim, comp, LocalOOP13, LocalOOP23, &
          spoofdim, FabVarName, FirstTime, ElGradVals, LocalGrad,gamma0,&
          LocalFabric, LocalLHS, ElLHSVals, nlm
!------------------------------------------------------------------------------
     CHARACTER(LEN=MAX_NAME_LEN) :: viscosityFile, TempVar, &
     OOPlaneRotVar13, OOPLaneRotVar23, FabVarName

     REAL(KIND=dp) :: Bu,Bv,Bw,RM(3,3), SaveTime = -1
     REAL(KIND=dp), POINTER :: PrevFabric(:),CurrFabric(:),TempFabVal(:)

     INTEGER, PARAMETER :: LCap=4, fab_len=15

     SAVE  ViscosityFile, PrevFabric, CurrFabric,TempFabVal
#ifdef USE_ISO_C_BINDINGS
     REAL(KIND=dp) :: at, at0
#else
     REAL(KIND=dp) :: at, at0, CPUTime, RealTime
#endif
!------------------------------------------------------------------------------
      INTERFACE
        Subroutine R2Ro(a2,dim,spoofdim,ai,angle)
        USE Types
        REAL(KIND=dp),intent(in) :: a2(6)
        Integer :: dim, spoofdim
        REAL(KIND=dp),intent(out) :: ai(3), Angle(3)
       End Subroutine R2Ro
      End Interface                                                       
!------------------------------------------------------------------------------
!  Read constants from constants section of SIF file
!------------------------------------------------------------------------------
      IF (FIRSTTIME) THEN
        CALL initspecfab(Lcap)
        FIRSTTIME = .FALSE.
      END IF

      Wn(7) = ListGetConstReal( Model % Constants, 'Gas Constant', GotIt,UnFoundFatal=UnFoundFatal )
      WRITE(Message,'(A,F10.4)')'Gas Constant =',Wn(7)
      CALL INFO('FabricSolveRecryst',Message,Level=4)
!------------------------------------------------------------------------------
!    Get variables needed for solution
!------------------------------------------------------------------------------
      IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN

      Solution => Solver % Variable % Values
      STDOFs   =  Solver % Variable % DOFs

      SolverParams => GetSolverParams()
      FabVarName = ListGetString( SolverParams,'Fabric Name',GotIt,UnFoundFatal=.FALSE. )
      IF (.NOT.GotIt) THEN
          FabVarName = 'Fabric'
          WRITE(Message,'(A,A)') 'Fabric name unfound, using ', FabVarName
          CALL INFO('FabricSolveRecryst', Message , level = 3)
      END IF
      FabricSol    => VariableGet( Solver % Mesh % Variables, FabVarName )
      IF ( ASSOCIATED( FabricSol) ) THEN
        FabricPerm   => FabricSol % Perm
        FabricValues => FabricSol % Values
      ELSE
        WRITE(Message,'(A,A,A)') 'Fabric variable ', FabVarName, &
            ' not found, quitting'
        CALL FATAL('FabricSolveRecryst', Message)
      END IF

      IF (FabricSol % DOFs .NE. fab_len) THEN
        WRITE(Message,'(A,A,A)') 'Fabric variable ', FabVarName,&
                                 ' must have fab_len DoFs'
        CALL FATAL('FabricSolveRecryst', Message)
      END IF

      TempVar = ListGetString( SolverParams,'Temperature Solution Name',&
            GotIt,UnFoundFatal=.FALSE. )
      IF (.NOT.GotIt) THEN
          TempVar = 'Temperature'
      END IF
      TempSol => VariableGet( Solver % Mesh % Variables, TempVar )
      IF ( ASSOCIATED( TempSol) ) THEN
        TempPerm    => TempSol % Perm
        Temperature => TempSol % Values
      END IF
      WRITE(Message,'(A,A)') 'Temperature variable = ', TempVar
      CALL INFO('FabricSolveRecryst', Message , level = 20)
      OOPlaneRotVar23 = ListGetString( SolverParams, &
        'OOPlane23 Strain Name',GotIt,UnFoundFatal=.FALSE.)
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
      CALL INFO('FabricSolveRecryst', Message , level = 20)

      OOPlaneRotVar13 = ListGetString( SolverParams,&
        'OOPlane13 Strain Name',GotIt,UnFoundFatal=.FALSE.)
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
      WRITE(Message,'(A,A)') 'OOPlane13 variable = ', OOPlaneRotVar13
      CALL INFO('FabricSolveRecryst', Message , level = 20)
      FlowVariable => VariableGet( Solver % Mesh % Variables, 'AIFlow' )
      IF ( ASSOCIATED( FlowVariable ) ) THEN
       FlowPerm    => FlowVariable % Perm
       FlowValues  => FlowVariable % Values
      END IF
      
!!!!! Mesh Velo
     MeshVeloVariable => VariableGet( Solver % Mesh % Variables, &
            'Mesh Velocity' )

     IF ( ASSOCIATED( MeshVeloVariable ) ) THEN
       MeshVeloPerm    => MeshVeloVariable % Perm
       MeshVeloValues  => MeshVeloVariable % Values
     END IF
       
                                       
      StiffMatrix => Solver % Matrix
      ! UNorm = Solver % Variable % Norm
      Unorm = SQRT( SUM( FabricValues**2 ) / SIZE(FabricValues) )
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!     Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
      IF ( .NOT. AllocationsDone .OR. Solver % Mesh % Changed) THEN
        N = Model % MaxElementNodes

       dim = CoordinateSystemDimension()
       IF ( OOPlaneRot13 .OR. OOPlaneRot23 ) THEN
           spoofdim = dim + 1
       ELSE
           spoofdim = dim
       END IF
       
       IF ( AllocationsDone ) THEN
         DEALLOCATE( LocalTemperature, &
                     LocalFabric, &
                     ElGradVals, &
                     ElLHSVals, &
                     LocalGrad,&
                     LocalLHS,&
                     Force, LocalFluidity, &
                     Velocity,MeshVelocity, &
                     MASS,STIFF,      &
                     LOAD, Alpha, Beta, &
                     CurrFabric, TempFabVal, &
                     nlm, &
                     LocalOOP13, LocalOOP23 )
       END IF

       ALLOCATE( LocalTemperature( N ), LocalFluidity( N ), &
                 LocalFabric(fab_len, N),&
                 nlm(fab_len),&
                 ElGradVals(Solver % NumberOFActiveElements, fab_len, N), &
                 ElLHSVals(Solver % NumberOFActiveElements, fab_len, N), &
                 LocalGrad(fab_len,N), &
                 LocalLHS(fab_len,N), &
                 Force( 2*STDOFs*N ), &
                 Velocity(4, N ),MeshVelocity(3,N), &
                 MASS( 2*STDOFs*N,2*STDOFs*N ),  &
                 STIFF( 2*STDOFs*N,2*STDOFs*N ),  &
                 LOAD( 4,N ), Alpha( 3,N ), Beta( N ), &
                 CurrFabric( fab_len*SIZE(Solver % Variable % Values)), &
                 TempFabVal( SIZE(FabricValues)), &
                 LocalOOP13(N), LocalOOP23(N), &
                 STAT=istat )


       IF ( istat /= 0 ) THEN
          CALL Fatal( 'FabricSolveRecryst', 'Memory allocation error.' )
       END IF

       CurrFabric = 0.
       TempFabVal = 0.
       IF ( TransientSimulation ) THEN
          IF (AllocationsDone ) DEALLOCATE (PrevFabric)
          ALLOCATE( PrevFabric( fab_len*SIZE(Solver % Variable % Values)) )
         PrevFabric = 0.
       END IF

       DO i=1,Solver % NumberOFActiveElements
          CurrentElement => GetActiveElement(i)   
          n = GetElementDOFs( Indexes )
          n = GetElementNOFNodes()
          NodeIndexes => CurrentElement % NodeIndexes
          Indexes(1:n) = Solver % Variable % Perm( Indexes(1:n) )
          DO COMP=1,fab_len
            IF ( TransientSimulation ) THEN
               PrevFabric(fab_len*(Indexes(1:n)-1)+COMP) = &
                   FabricValues(fab_len*(FabricPerm(NodeIndexes(1:n))-1)+COMP)
            END IF
               CurrFabric(fab_len*(Indexes(1:n)-1)+COMP) = &
                   FabricValues(fab_len*(FabricPerm(NodeIndexes(1:n))-1)+COMP)
          END DO
       END DO

       AllocationsDone = .TRUE.
      END IF

      IF( TransientSimulation ) THEN
        TimeVar => VariableGet( Solver % Mesh % Variables, 'Time' )
        IF ( SaveTime /= TimeVar % Values(1) ) THEN
           SaveTime = TimeVar % Values(1)
           PrevFabric = CurrFabric
        END IF
      END IF

!------------------------------------------------------------------------------
!    Do some additional initialization, and go for it
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
      NonlinearTol = ListGetConstReal( Solver % Values, &
        'Nonlinear System Convergence Tolerance' )

      NewtonTol = ListGetConstReal( Solver % Values, &
        'Nonlinear System Newton After Tolerance' )

      NewtonIter = ListGetInteger( Solver % Values, &
        'Nonlinear System Newton After Iterations' )

      NonlinearIter = ListGetInteger( Solver % Values, &
         'Nonlinear System Max Iterations',GotIt )

      IF ( .NOT.GotIt ) NonlinearIter = 1
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
      DO iter=1,NonlinearIter
!------------------------------------------------------------------------------
       at  = CPUTime()
       at0 = RealTime()
       CALL Info( 'FabricSolve', ' ', Level=4 )
       CALL Info( 'FabricSolve', ' ', Level=4 )
       CALL Info( 'FabricSolve', &
                    '-------------------------------------',Level=4 )
       WRITE( Message, * ) 'Fabric solver  iteration', iter
       CALL Info( 'FabricSolve', Message,Level=4 )
       CALL Info( 'FabricSolve', &
                     '-------------------------------------',Level=4 )
       CALL Info( 'FabricSolve', ' ', Level=4 )
       CALL Info( 'FabricSolve', 'Starting assembly...',Level=4 )
       PrevUNorm = UNorm

       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !! Calculate and save gradients once per iteration
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       DO t=1,Solver % NumberOFActiveElements
         CurrentElement => GetActiveElement(t)
         CALL GetElementNodes( ElementNodes )
         n = GetElementDOFs( Indexes )
         n = GetElementNOFNodes()
         NodeIndexes => CurrentElement % NodeIndexes
         Material => GetMaterial()
         body_id = CurrentElement % BodyId
         IF (body_id /= old_body) Then 
           old_body = body_id
           CALL GetMaterialDefs()
         END IF
         LocalFluidity(1:n) = ListGetReal( Material, &
                         'Fluidity Parameter', n, NodeIndexes, GotIt,&
                         UnFoundFatal=UnFoundFatal)
!------------------------------------------------------------------------------
!        Get element local stiffness & mass matrices
!------------------------------------------------------------------------------
         LocalTemperature = 0.0D0
         IF ( ASSOCIATED(TempSol) ) THEN
            DO i=1,n
               k = TempPerm(NodeIndexes(i))
               LocalTemperature(i) = Temperature(k)
            END DO
         ELSE
            LocalTemperature(1:n) = 0.0d0
         END IF

         DO i = 1,fab_len
            LocalFabric(i, 1:n) = CurrFabric( fab_len*(Solver % Variable %Perm(Indexes(1:n))-1)+i)
         END DO

         ! Two variables needed for velocity
         k = FlowVariable % DOFs
         Velocity = 0.0d0
         DO i=1,k-1
            Velocity(i,1:n) = FlowValues(k*(FlowPerm(NodeIndexes)-1)+i)
         END DO
         MeshVelocity=0._dp
         IF (ASSOCIATED(MeshVeloVariable)) Then
           k = MeshVeloVariable % DOFs
           DO i=1,k
              MeshVelocity(i,1:n) = MeshVeloValues(k*(MeshVeloPerm(NodeIndexes)-1)+i)
           END DO
         EndIF

         LocalOOP13 = 0.0_dp
         LocalOOP23 = 0.0_dp
         IF (spoofdim.gt.dim) THEN
            LocalOOP13(1:n) = OOPlaneRotValues13(OOPlaneRotPerm13(NodeIndexes))
            LocalOOP23(1:n) = OOPlaneRotValues23(OOPlaneRotPerm23(NodeIndexes))
        END IF
         CALL FabGrad( LocalGrad, LocalLHS, fab_len, LocalFabric, &
                       LocalTemperature, LocalFluidity,  Velocity, &
                       MeshVelocity, CurrentElement, n, ElementNodes, &
                       Wn, rho, lambda0, gamma0, LocalOOP23, LocalOOP13)
         ElGradVals(t, :, :) = LocalGrad(:, :)
         ElLHSVals(t, :, :) = LocalGrad(:, :)
       END DO  ! active elements
 
       
       ! Note that the first component of the fabric is unchanged
       ! (normalization)
       DO COMP=1,fab_len
         IF ((SPOOFDIM.LE.2).AND.ANY(COMP==(/1, 3, 5, 8, 10, 12, 14/))) CYCLE
       Solver % Variable % Values = CurrFabric( COMP::fab_len )
       IF ( TransientSimulation ) THEN
          Solver % Variable % PrevValues(:,1) = PrevFabric( COMP::fab_len )
       END IF

       CALL DefaultInitialize()
!------------------------------------------------------------------------------
       DO t=1,Solver % NumberOFActiveElements
!------------------------------------------------------------------------------

         IF ( RealTime() - at0 > 1.0 ) THEN
           WRITE(Message,'(a,i3,a)' ) '   Assembly: ', INT(100.0 - 100.0 * &
            (Solver % NumberOfActiveElements-t) / &
               (1.0*Solver % NumberOfActiveElements)), ' % done'
                       
           CALL Info( 'FabricSolve', Message, Level=5 )
           at0 = RealTime()
         END IF

         CurrentElement => GetActiveElement(t)
         CALL GetElementNodes( ElementNodes )
         n = GetElementDOFs( Indexes )
         n = GetElementNOFNodes()
         NodeIndexes => CurrentElement % NodeIndexes

         
         Material => GetMaterial()
         body_id = CurrentElement % BodyId
!------------------------------------------------------------------------------
!        Read in material constants from Material section
!------------------------------------------------------------------------------
         IF (body_id /= old_body) Then 
           old_body = body_id
           CALL GetMaterialDefs()
         END IF
      
         LocalFluidity(1:n) = ListGetReal( Material, &
                         'Fluidity Parameter', n, NodeIndexes, GotIt,&
                         UnFoundFatal=UnFoundFatal)
!------------------------------------------------------------------------------
!        Get element local stiffness & mass matrices
!------------------------------------------------------------------------------
         LocalTemperature = 0.0D0
         IF ( ASSOCIATED(TempSol) ) THEN
            DO i=1,n
               k = TempPerm(NodeIndexes(i))
               LocalTemperature(i) = Temperature(k)
            END DO
         ELSE
            LocalTemperature(1:n) = 0.0d0
         END IF

         k = FlowVariable % DOFs
         Velocity = 0.0d0
         DO i=1,k-1
            Velocity(i,1:n) = FlowValues(k*(FlowPerm(NodeIndexes)-1)+i)
         END DO

!------------meshvelocity
         MeshVelocity=0._dp
         IF (ASSOCIATED(MeshVeloVariable)) Then
           k = MeshVeloVariable % DOFs
           DO i=1,k
              MeshVelocity(i,1:n) = MeshVeloValues(k*(MeshVeloPerm(NodeIndexes)-1)+i)
           END DO
         EndIF
!----------------------------------

        LocalGrad(:, :) = ElGradVals(t, :, :)
        LocalLHS(:, :) = ElLHSVals(t, :, :)
        CALL LocalMatrix( comp, MASS, STIFF, FORCE, LOAD, LocalGrad, &
                          LocalLHS, Velocity, MeshVelocity, &
                          CurrentElement, n, ElementNodes, spoofdim)

!------------------------------------------------------------------------------
!        Update global matrices from local matrices 
!------------------------------------------------------------------------------
         IF ( TransientSimulation )  CALL Default1stOrderTime(MASS,STIFF,FORCE)
         CALL DefaultUpdateEquations( STIFF, FORCE )
!------------------------------------------------------------------------------
      END DO ! Active Elements
      CALL Info( 'FabricSolve', 'Assembly done', Level=4 )
!------------------------------------------------------------------------------
   
!------------------------------------------------------------------------------
!     Assembly of the edge terms
!------------------------------------------------------------------------------
!
!3D => Edges => Faces
      If (dim.eq.3) then 
      DO t=1,Solver % Mesh % NumberOfFaces
         Edge => Solver % Mesh % Faces(t)
         IF ( .NOT. ActiveBoundaryElement(Edge) ) CYCLE
       
         LeftParent  => Edge % BoundaryInfo % Left
         RightParent => Edge % BoundaryInfo % Right
         IF ( ASSOCIATED(RightParent) .AND. ASSOCIATED(LeftParent) ) THEN
            n  = GetElementNOFNodes( Edge )
            n1 = GetElementNOFNodes( LeftParent )
            n2 = GetElementNOFNodes( RightParent )

            k = FlowVariable % DOFs
            Velocity = 0.0d0
            DO i=1,k
               Velocity(i,1:n) = FlowValues(k*(FlowPerm(Edge % NodeIndexes)-1)+i)
            END DO

     !-------------------mesh velo
         MeshVelocity=0._dp
      IF ( ASSOCIATED( MeshVeloVariable ) ) THEN
            k = MeshVeloVariable % DOFs
            DO i=1,k
               MeshVelocity(i,1:n) = MeshVeloValues(k*(MeshVeloPerm(Edge % NodeIndexes)-1)+i)
            END DO
      END IF
     !--------------------------

            FORCE = 0.0d0
            MASS  = 0.0d0
            CALL LocalJumps( STIFF,Edge,n,LeftParent,n1,RightParent,n2,Velocity,MeshVelocity )
            IF ( TransientSimulation )  CALL Default1stOrderTime(MASS, STIFF, FORCE)
            CALL DefaultUpdateEquations( STIFF, FORCE, Edge )
         END IF
      END DO
!
!  2D
      Else

      DO t=1,Solver % Mesh % NumberOfEdges
         Edge => Solver % Mesh % Edges(t)
         IF ( .NOT. ActiveBoundaryElement(Edge) ) CYCLE
       
         LeftParent  => Edge % BoundaryInfo % Left
         RightParent => Edge % BoundaryInfo % Right
         IF ( ASSOCIATED(RightParent) .AND. ASSOCIATED(LeftParent) ) THEN
            n  = GetElementNOFNodes( Edge )
            n1 = GetElementNOFNodes( LeftParent )
            n2 = GetElementNOFNodes( RightParent )

            k = FlowVariable % DOFs
            Velocity = 0.0d0
            DO i=1,k
               Velocity(i,1:n) = FlowValues(k*(FlowPerm(Edge % NodeIndexes(1:n))-1)+i)
            END DO

     !-------------------mesh velo
         MeshVelocity=0._dp
      IF ( ASSOCIATED( MeshVeloVariable ) ) THEN
            k = MeshVeloVariable % DOFs
            DO i=1,k
               MeshVelocity(i,1:n) = MeshVeloValues(k*(MeshVeloPerm(Edge % NodeIndexes)-1)+i)
            END DO
      END IF
     !--------------------------

            FORCE = 0.0d0
            MASS  = 0.0d0
            CALL LocalJumps( STIFF,Edge,n,LeftParent,n1,RightParent,n2,Velocity,MeshVelocity )
            IF ( TransientSimulation )  CALL Default1stOrderTime(MASS, STIFF, FORCE)
            CALL DefaultUpdateEquations( STIFF, FORCE, Edge )
         END IF
       END DO

      END IF

      CALL DefaultFinishBulkAssembly()
!------------------------------------------------------------------------------
!     Loop over the boundary elements
!------------------------------------------------------------------------------
      DO t = 1, Solver % Mesh % NumberOfBoundaryElements
!------------------------------------------------------------------------------

         Element => GetBoundaryElement(t)
         IF( .NOT. ActiveBoundaryElement() )  CYCLE
         IF( GetElementFamily(Element) == 1 ) CYCLE

         ParentElement => Element % BoundaryInfo % Left
         IF ( .NOT. ASSOCIATED( ParentElement ) ) &
            ParentElement => Element % BoundaryInfo % Right
          
         n  = GetElementNOFNodes( Element )
         n1 = GetElementNOFnodes( ParentElement )
       
         k = FlowVariable % DOFs
         Velocity = 0.0d0
         DO i=1,k
            Velocity(i,1:n) = FlowValues(k*(FlowPerm(Element % NodeIndexes(1:n))-1)+i)
         END DO

!-------------------mesh velo
         MeshVelocity=0._dp
      IF ( ASSOCIATED( MeshVeloVariable ) ) THEN
        k = MeshVeloVariable % DOFs
        DO i=1,k
         MeshVelocity(i,1:n) = MeshVeloValues(k*(MeshVeloPerm(Element % NodeIndexes)-1)+i)
        End do
      END IF
!--------------------------


         BC => GetBC()
         LOAD = 0.0d0
         GotIt = .FALSE.
         IF ( ASSOCIATED(BC) ) THEN
            LOAD(1,1:n) = GetReal( BC, ComponentName('Fabric', Comp) , GotIt )
         END IF

         MASS = 0.0d0
         CALL LocalMatrixBoundary(  STIFF, FORCE, LOAD(1,1:n), &
                              Element, n, ParentElement, n1, Velocity,MeshVelocity, GotIt )

         IF ( TransientSimulation )  CALL Default1stOrderTime(MASS, STIFF, FORCE)
         CALL DefaultUpdateEquations( STIFF, FORCE )
      END DO

      CALL DefaultFinishAssembly()
!------------------------------------------------------------------------------
      CALL Info( 'FabricSolve', 'Set boundaries done', Level=4 )

!------------------------------------------------------------------------------
!     Solve the system and check for convergence
!------------------------------------------------------------------------------
      Unorm = DefaultSolve()
      WRITE(Message,*) 'solve done', minval( solver % variable % values), maxval( Solver % variable % values)
      CALL Info( 'FabricSolve', Message, Level=4 )
      
      n1 = Solver % Mesh % NumberOfNodes
      ALLOCATE( Ref(n1) )
      Ref = 0
      TempFabVal(COMP::fab_len ) = 0. !fab
      
      DO t=1,Solver % NumberOfActiveElements
         Element => GetActiveElement(t) 
         n = GetElementDOFs( Indexes )
         n = GetElementNOFNodes()
         
         DO i=1,n
            k = Element % NodeIndexes(i)
            TempFabVal( fab_len*(FabricPerm(k)-1) + COMP ) =    & 
            TempFabVal( fab_len*(FabricPerm(k)-1) + COMP ) + &
            Solver % Variable % Values( Solver % Variable % Perm(Indexes(i)) )
            FabricValues( fab_len*(FabricPerm(k)-1) + COMP ) = &
                          TempFabVal(fab_len*(FabricPerm(k)-1) + COMP ) 
            Ref(k) = Ref(k) + 1
         END DO
      END DO

      DO i=1,n1
         j=FabricPerm(i)
         IF (j < 1) CYCLE
         IF ( Ref(i) > 0 ) THEN
            FabricValues( fab_len*(j-1)+COMP ) = &
                   FabricValues( fab_len*(j-1)+COMP ) / Ref(i)
         END IF
      END DO

      DEALLOCATE( Ref )
      END DO ! End DO Comp

       DO i=1,Solver % NumberOFActiveElements
          CurrentElement => GetActiveElement(i)   
          n = GetElementDOFs( Indexes )
          n = GetElementNOFNodes()
          NodeIndexes => CurrentElement % NodeIndexes
          Indexes(1:n) = Solver % Variable % Perm( Indexes(1:n) )
          DO COMP=1,fab_len
            CurrFabric(fab_len*(Indexes(1:n)-1)+COMP) = &
                        FabricValues(fab_len*(FabricPerm(NodeIndexes(1:n))-1)+COMP)
          END DO
       END DO


      Unorm = SQRT( SUM( FabricValues**2 ) / SIZE(FabricValues) )
      Solver % Variable % Norm = Unorm  

      IF ( PrevUNorm + UNorm /= 0.0d0 ) THEN
         RelativeChange = 2.0d0 * ABS( PrevUNorm - UNorm) / ( PrevUnorm + UNorm)
      ELSE
         RelativeChange = 0.0d0
     END IF

      WRITE( Message, * ) 'Result Norm   : ',UNorm
      CALL Info( 'FabricSolve', Message, Level=4 )
      WRITE( Message, * ) 'Relative Change : ',RelativeChange
      CALL Info( 'FabricSolve', Message, Level=4 )

      
!------------------------------------------------------------------------------
      IF ( RelativeChange < NewtonTol .OR. &
            iter > NewtonIter ) NewtonLinearization = .TRUE.

      IF ( RelativeChange < NonLinearTol ) EXIT
      TensorFabricVariable => &
       VariableGet( Solver % Mesh % Variables, 'TensorFabric' )
     IF ( ASSOCIATED( TensorFabricVariable ) ) THEN
         TensorFabricPerm  => TensorFabricVariable % Perm
         TensorFabricValues => TensorFabricVariable % Values
      
         DO t=1,Solver % NumberOFActiveElements

           CurrentElement => GetActiveElement(t)
           n = GetElementNOFNodes()
           NodeIndexes => CurrentElement % NodeIndexes

           DO i = 1,fab_len
             LocalFabric(i, 1:n) = CurrFabric( fab_len*(Solver % Variable %Perm(Indexes(1:n))-1)+i)
           END DO

           Do i=1,n
             nlm = CMPLX(PACK(LocalFabric(:, i),.TRUE.), KIND=dp)
             a2short(:) = a2_to_ae2(a2_ij(nlm))
             TensorFabricValues(5*(TensorFabricPerm(NodeIndexes(i))-1)+1:&
                                5*(TensorFabricPerm(NodeIndexes(i))-1)+5)=&
               a2short(:)
           END DO
         END DO
      END IF
!------------------------------------------------------------------------------
    END DO ! of nonlinear iter
!------------------------------------------------------------------------------
CONTAINS

      SUBROUTINE GetMaterialDefs()

      viscosityFile = ListGetString( Material ,'Viscosity File',GotIt, UnFoundFatal)
      OPEN( 1, File = viscosityFile)
      DO i=1,813
         READ( 1, '(6(e14.8))' ) FabricGrid( 6*(i-1)+1:6*(i-1)+6 )
      END DO
      READ(1 , '(e14.8)' ) FabricGrid(4879)
      CLOSE(1)

       rho = ListGetConstReal(Material, 'Interaction Parameter', GotIt )
       IF (.NOT.GotIt) THEN
           WRITE(Message,'(A)') 'Interaction  Parameter notfound. &
                         &Setting to the value in ViscosityFile'
           CALL INFO('AIFlowSolve', Message, Level = 20)
           rho = FabricGrid(4879)
       ELSE
           WRITE(Message,'(A,F10.4)') 'Interaction Parameter = ', rho
           CALL INFO('AIFlowSolve', Message, Level = 20)
       END IF

       lambda0 = ListGetConstReal( Material, 'Diffusion Intercept', GotIt,UnFoundFatal=UnFoundFatal)
           !Previous default value: lambda = 0.0_dp
      WRITE(Message,'(A,F10.4)') 'Diffusion Intercept = ', lambda0
      CALL INFO('AIFlowSolve', Message, Level = 20)

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
           !Previous default value: Wn(5) = -10.0
      WRITE(Message,'(A,F10.4)') 'Reference Temperature = ',   Wn(5)
      CALL INFO('AIFlowSolve', Message, Level = 20)

      Wn(6) = ListGetConstReal( Material, 'Limit Temperature', GotIt,UnFoundFatal=UnFoundFatal)
           !Previous default value: Wn(6) = -10.0
      WRITE(Message,'(A,F10.4)') 'Limit Temperature = ',   Wn(6)
      CALL INFO('AIFlowSolve', Message, Level = 20)


      Wn(8) = ListGetConstReal( Material, 'Migration A', GotIt,UnFoundFatal=.FALSE.)
           !Previous default value: Wn(6) = -10.0
      IF (.NOT.GotIt) THEN
        Wn(8) = 0.0_dp
        WRITE(Message,'(A,F10.4)') &
            'Migration A unfound, taken to be ',   Wn(8)
        CALL INFO('AIFlowSolve', Message, Level = 3)
      ELSE
        WRITE(Message,'(A,F10.4)') 'Migration Gamma = ',   Wn(8)
        CALL INFO('AIFlowSolve', Message, Level = 20)
      END IF

      Wn(9) = ListGetConstReal( Material, 'Lattice Rotation', GotIt,UnFoundFatal=.FALSE.)
      IF (.NOT.GotIt) THEN
        Wn(9) = 1.0_dp
        WRITE(Message,'(A,F10.4)') &
            'Lattice Rotation unfound, assumed active'
        CALL INFO('AIFlowSolve', Message, Level = 3)
      ELSE
        WRITE(Message,'(A,F10.4)') 'Lattice Rotation = ',   Wn(9)
        CALL INFO('AIFlowSolve', Message, Level = 20)
      END IF

      Wn(10) = ListGetConstReal( Material, 'Diffusion Temp Dependence', GotIt,UnFoundFatal=.FALSE.)
      IF (.NOT.GotIt) THEN
        Wn(10) = 0.0_dp
        WRITE(Message,'(A,F10.4)') &
            'Diffusion temp dependence unfound, assumed zero'
        CALL INFO('AIFlowSolve', Message, Level = 3)
      END IF

      Wn(11) = ListGetConstReal( Material, 'Max Diffusion', GotIt,UnFoundFatal=.FALSE.)
      IF (.NOT.GotIt) THEN
        Wn(11) = 1.0e8_dp
      END IF

      Wn(12) = ListGetConstReal( Material, 'Max Migration', GotIt,UnFoundFatal=.FALSE.)
      IF (.NOT.GotIt) THEN
        Wn(12) = 1.0e8_dp
      END IF

      gamma0 = ListGetConstReal( Material, 'Migration Prefactor',GotIt,UnFoundFatal=.TRUE.)
      WRITE(Message,'(A,F10.4)') 'Migration prefactor = ', gamma0
      CALL INFO('AIFlowSolve', Message, Level = 20)

!------------------------------------------------------------------------------
      END SUBROUTINE GetMaterialDefs
!------------------------------------------------------------------------------

      SUBROUTINE FabGrad( Gradient, LHS, fab_len, NodalFabric, &
                          NodalTemperature, NodalFluidity, NodalVelo, &
                          NodMeshVel, Element, n, Nodes, Wn, rho, &
                          lambda0,gamma0,LocalOOP23,LocalOOP13)
!------------------------------------------------------------------------------
     INTEGER :: fab_len
     REAL(KIND=dp) :: NodalVelo(:,:),NodMeshVel(:,:),NodalFabric(:,:)
     REAL(KIND=dp), DIMENSION(:) :: NodalTemperature, NodalFluidity
     REAL(KIND=dp), DIMENSION(:), OPTIONAL :: LocalOOP23, LocalOOP13
     REAL(KIND=dp), Intent(OUT) :: Gradient(:,:), LHS(:,:)

     TYPE(Nodes_t) :: Nodes
     TYPE(Element_t) :: Element
     INTEGER :: n, Comp
!------------------------------------------------------------------------------
!
     REAL(KIND=dp) :: Basis(2*n),ddBasisddx(1,1,1),NodalGradient(fab_len)
     REAL(KIND=dp) :: dBasisdx(2*n,3),SqrtElementMetric

     REAL(KIND=dp) :: Theta
     REAL(KIND=dp) :: OOP23

     REAL(KIND=dp) :: A,M, tau,pe1,pe2,unorm,C0, SU(n), SW(n)
     REAL(KIND=dp) :: LoadAtIp, Temperature
     REAL(KIND=dp) :: rho,Deq,ai(6),hmax
     REAL(KIND=dp) :: lambda0, gamma0, lambda, gammav

     INTEGER :: i,j,k,p,q,t,dim,NBasis,ind(3),spoofdim,DOFs = 1

     REAL(KIND=dp) :: s,u,v,w, Radius, B(6,3), G(3,6), a2full(3, 3)
     REAL(KIND=dp) :: Wn(:),Velo(3),DStress(6),StrainR(6),Spin(3),SD(6)

     REAL(KIND=dp) :: LGrad(3,3),StrainRate(3,3),D(6),angle(3),epsi
     REAL(KIND=dp) :: ap(3),C(6,6),Spin1(3,3),Stress(3,3),eps(3,3)
     REAL(KIND=dp) :: SStar(3,3), SStarMean, TrS, Fabric(fab_len)
     COMPLEX(KIND=dp) :: dndt(fab_len, fab_len), dndt_ROT(fab_len, fab_len),&
                      dndt_DDRX(fab_len, fab_len), dndt_CDRX(fab_len, fab_len)
     Integer :: INDi(6),INDj(6)
     INTEGER :: N_Integ
     REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ

     LOGICAL :: stat
     TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

     INTERFACE
        FUNCTION BGlenT( Tc, W)
           USE Types
           REAL(KIND=dp) :: BGlenT,Tc,W(8)
        END FUNCTION
         Subroutine R2Ro(ai,dim,spoofdim,a2,angle)
         USE Types
         REAL(KIND=dp),intent(in) :: ai(6)
         Integer :: dim, spoofdim
         REAL(KIND=dp),intent(out) :: a2(3), Angle(3)
        End Subroutine R2Ro

        Subroutine OPILGGE_ai_nl(a2,Angle,etaI,eta36)
          USE Types
          REAL(kind=dp), INTENT(in),  DIMENSION(3)   :: a2
          REAL(kind=dp), INTENT(in),  DIMENSION(3)   :: Angle
          REAL(kind=dp), INTENT(in),  DIMENSION(:)   :: etaI
          REAL(kind=dp), INTENT(out), DIMENSION(6,6) :: eta36
        END SUBROUTINE OPILGGE_ai_nl
        
      END INTERFACE

      dim = CoordinateSystemDimension()
      IF (( PRESENT(LocalOOP13) ).OR.( PRESENT(LocalOOP23))) THEN
          SpoofDim = dim + 1
      ELSE
          SpoofDim = dim
      END IF

!    Integration stuff:
!    ------------------
      NBasis = n
      IntegStuff = GaussPoints( Element  )

      U_Integ => IntegStuff % u
      V_Integ => IntegStuff % v
      W_Integ => IntegStuff % w
      S_Integ => IntegStuff % s
      N_Integ =  IntegStuff % n


      LHS = 0.0_dp
      Gradient = 0.0_dp
!
!   Now we start integrating:
!   -------------------------
      DO t=1,N_Integ

      u = U_Integ(t)
      v = V_Integ(t)
      w = W_Integ(t)

!------------------------------------------------------------------------------
!     Basis function values & derivatives at the integration point
!------------------------------------------------------------------------------
      stat = ElementInfo( Element,Nodes,u,v,w,SqrtElementMetric, &
               Basis, dBasisdx, ddBasisddx, .FALSE. )

      s = SqrtElementMetric * S_Integ(t)
!------------------------------------------------------------------------------
!
!     Orientation parameters at the integration point:
!     ------------------------------------------------
      Temperature = SUM( NodalTemperature(1:n) * Basis(1:n) ) 

      DO i=1,fab_len
       Fabric(i) = SUM(NodalFabric(i,1:n)*Basis(1:n))
      END DO
      Theta = 1._dp / ( FabricGrid(5) + FabricGrid(6) )
      
      Stress = 0.0
      StrainRate = 0.0
      Spin1 = 0.0

      a2full = a2_ij(CMPLX(Fabric, KIND=dp))
!
!    Material parameters at that point
!    ---------------------------------
      ai(1) = a2full(1, 1)
      ai(2) = a2full(2, 2)
      ai(3) = a2full(3, 3)
      ai(4) = a2full(1, 2)
      ai(5) = a2full(2, 3)
      ai(6) = a2full(1, 3)

      call R2Ro(ai,dim,spoofdim,ap,angle)
      CALL OPILGGE_ai_nl(ap, Angle, FabricGrid, C)
      LGrad = MATMUL( NodalVelo(1:3,1:n), dBasisdx(1:n,1:3) )

      StrainRate = 0.5 * ( LGrad + TRANSPOSE(LGrad) )

      Spin1 = 0.5 * ( LGrad - TRANSPOSE(LGrad) )
      epsi = StrainRate(1,1)+StrainRate(2,2)+StrainRate(3,3)
        DO i=1,dim 
          StrainRate(i,i) = StrainRate(i,i) - epsi/dim
        END DO

      ! If we are adding in manual spins, do it here before the stress
      ! calculations
      IF ( PRESENT(LocalOOP23) ) THEN
          StrainRate(2, 3) = SUM( LocalOOP23(1:n) * Basis(1:n) )
          StrainRate(3, 2) = StrainRate(2, 3)
      END IF
      IF ( PRESENT(LocalOOP13) ) THEN
          StrainRate(3, 1) = SUM( LocalOOP13(1:n) * Basis(1:n) )
          StrainRate(1, 3) = StrainRate(3, 1)
      END IF

!
!    Compute deviatoric stresses: 
!    ----------------------------
      D(1) = StrainRate(1,1)
      D(2) = StrainRate(2,2)
      D(3) = StrainRate(3,3)
      D(4) = 2. * StrainRate(1,2)
      D(5) = 2. * StrainRate(2,3)
      D(6) = 2. * StrainRate(3,1)
      
      INDi(1:6) = (/ 1, 2, 3, 1, 2, 3 /)
      INDj(1:6) = (/ 1, 2, 3, 2, 3, 1 /)
      DO k = 1, 2*spoofdim
       DO j = 1, 2*spoofdim
        Stress( INDi(k),INDj(k) ) = &
        Stress( INDi(k),INDj(k) ) + C(k,j) * D(j)
       END DO
       IF (k > 3)  Stress( INDj(k),INDi(k) ) = Stress( INDi(k),INDj(k) )
      END DO

      SStar = 0.0
      TrS = 0.0
      DO k = 1, 2 * spoofdim
        TrS = TrS + Stress(INDi(k), INDj(k)) * Stress(INDj(k), INDi(k))
        IF (k > 3) TrS = TrS + Stress(INDi(k), INDj(k)) * Stress(INDj(k), INDi(k))
      END DO
      IF (TrS.GE.1.0e-16) THEN
        DO k = 1, 2*spoofdim
          DO j = 1, 2*spoofdim
            SStar( INDi(k),INDj(k) ) = 5.0_dp * (Stress( INDi(k),INDj(k) )) / TrS
          END DO
          IF (k > 3)  SStar( INDj(k),INDi(k) ) = SStar( INDi(k),INDj(k) )
        END DO
      END IF

!     SD=(1-r)D + r psi/2 S :
!     -----------------------
      SD=0._dp
      DO i=1,2*spoofdim
        SD(i)= (1._dp - rho)*StrainRate(INDi(i),INDj(i)) + rho *&
                                   Theta *  Stress(INDi(i),INDj(i))
      END DO
      eps(1, 1) = SD(1)
      eps(1, 2) = SD(4)
      eps(1, 3) = SD(6)
      eps(2, 1) = SD(4)
      eps(2, 2) = SD(2)
      eps(2, 3) = SD(5)
      eps(3, 1) = SD(6)
      eps(3, 2) = SD(5)
      eps(3, 3) = SD(3)
    
      Deq=sqrt((SD(1)*SD(1)+SD(2)*SD(2)+SD(3)*SD(3)+2._dp* &
                             (SD(4)*SD(4)+SD(5)*SD(5)+SD(6)*SD(6)))/3._dp)

      ! simplest to do this in celcius
      lambda = MIN((lambda0 + Temperature * Wn(10)) * Deq, Wn(11))
      ! Arrhenius relations are in Kelvin
      gammav = Min((gamma0 * EXP(-Wn(8) / (Temperature + 273.15))) * Deq, Wn(12))

      dndt_ROT = dndt_ij_LATROT(EPS, Spin1, 0.0_dp * Strainrate,&
                                0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp)
      dndt_DDRX = dndt_ij_DDRX(CMPLX(Fabric, kind=dp), EPS)
      dndt_CDRX = f_nu_eps(lambda, StrainRate) * dndt_ij_REG()

      dndt = gammav * dndt_DDRX + dndt_CDRX + Wn(9) * dndt_ROT
      NodalGradient = REAL(MATMUL(dndt, Fabric))
      DO i=1,fab_len
        Gradient(i, t) = NodalGradient(i) - REAL(dndt_ROT(i, i) * Fabric(i))
        LHS(i, t) = REAL(dndt_ROT(i, i))
      END DO
      END DO ! N_Integ
       END SUBROUTINE FabGrad



!------------------------------------------------------------------------------
      SUBROUTINE LocalMatrix( Comp, MASS, STIFF, FORCE, LOAD, &
          Grad, LHS, NodalVelo, NodMeshVel, &
          Element, n, Nodes, dim)
!------------------------------------------------------------------------------
      ! Note that dim is increased by 1 for out-of-page stress/strain

      REAL(KIND=dp) :: STIFF(:,:),MASS(:,:)
      REAL(KIND=dp) :: LOAD(:,:),NodalVelo(:,:),NodMeshVel(:,:),&
                       Grad(:,:),LHS(:,:)
      REAL(KIND=dp), DIMENSION(:) :: FORCE
      TYPE(Nodes_t) :: Nodes
      TYPE(Element_t) :: Element
      INTEGER :: n, Comp

      REAL(KIND=dp) :: Basis(2*n),ddBasisddx(1,1,1)
      REAL(KIND=dp) :: dBasisdx(2*n,3),SqrtElementMetric
      REAL(KIND=dp) :: A,M,unorm,C0
      REAL(KIND=dp) :: LoadAtIp
      INTEGER :: i,j,k,p,q,t,dim,NBasis,ind(3),N_Integ
      REAL(KIND=dp) :: s,u,v,w, Radius, Velo(3)
      REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ,W_Integ,S_Integ
      LOGICAL :: stat
      TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff

      FORCE = 0.0D0
      MASS  = 0.0D0
      STIFF = 0.0D0
!    
!    Integration stuff:
!    ------------------
      NBasis = n
      IntegStuff = GaussPoints( Element  )

      U_Integ => IntegStuff % u
      V_Integ => IntegStuff % v
      W_Integ => IntegStuff % w
      S_Integ => IntegStuff % s
      N_Integ =  IntegStuff % n

      ! Start integrating:
      DO t=1,N_Integ
        u = U_Integ(t)
        v = V_Integ(t)
        w = W_Integ(t)

        !------------------------------------------------------------------------------
        !     Get basis function values & derivatives at the integration point
        !------------------------------------------------------------------------------
        stat = ElementInfo( Element,Nodes,u,v,w,SqrtElementMetric, &
               Basis, dBasisdx, ddBasisddx, .FALSE. )

        s = SqrtElementMetric * S_Integ(t)
        Velo = 0.0d0
        DO i=1,dim
         Velo(i) = SUM( Basis(1:n) * (NodalVelo(i,1:n) - NodMeshVel(i,1:n)) )
        END DO
        Unorm = SQRT( SUM( Velo**2._dp ) )

        C0 = LHS(COMP, t)
        LoadAtIp = Grad(COMP, t)
        DO p=1,NBasis
          DO q=1,NBasis
            A = 0.0d0
            M = Basis(p) * Basis(q)
            ! Reaction terms:
            A = A - C0 * Basis(q) * Basis(p)

            ! Advection terms:
            DO j=1,dim
               A = A - Velo(j) * Basis(q) * dBasisdx(p,j)
            END DO

            ! Add nodal matrix to element matrix:
            MASS( p,q )  = MASS( p,q )  + s * M
            STIFF( p,q ) = STIFF( p,q ) + s * A
          END DO
          FORCE(p) = FORCE(p) + s*LoadAtIp*Basis(p)
        END DO
      END DO 

 1000 FORMAT((a),x,i2,x,6(e13.5,x)) 
 1001 FORMAT(6(e13.5,x))
!------------------------------------------------------------------------------
       END SUBROUTINE LocalMatrix
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
      SUBROUTINE FindParentUVW( Edge, nEdge, Parent, nParent, U, V, W, Basis )
!------------------------------------------------------------------------------
        IMPLICIT NONE
        TYPE(Element_t), POINTER :: Edge, Parent
        INTEGER :: nEdge, nParent
        REAL( KIND=dp ) :: U, V, W, Basis(:)
        !------------------------------------------------------------------------------
        INTEGER :: i, j,l
        REAL(KIND=dp) :: NodalParentU(nEdge),NodalParentV(nEdge),NodalParentW(nEdge)
        !------------------------------------------------------------------------------
        DO i = 1,nEdge
        DO j = 1,nParent
          IF ( Edge % NodeIndexes(i) == Parent % NodeIndexes(j) ) THEN
            NodalParentU(i) = Parent % Type % NodeU(j)
            NodalParentV(i) = Parent % Type % NodeV(j)
            NodalParentW(i) = Parent % Type % NodeW(j)
            EXIT
          END IF
        END DO
        END DO
        U = SUM( Basis(1:nEdge) * NodalParentU(1:nEdge) )
        V = SUM( Basis(1:nEdge) * NodalParentV(1:nEdge) )
        W = SUM( Basis(1:nEdge) * NodalParentW(1:nEdge) )
!------------------------------------------------------------------------------      
      END SUBROUTINE FindParentUVW
!------------------------------------------------------------------------------      


!------------------------------------------------------------------------------
      SUBROUTINE LocalJumps( STIFF,Edge,n,LeftParent,n1,RightParent,n2,Velo,MeshVelo )
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: STIFF(:,:), Velo(:,:),MeshVelo(:,:)
      INTEGER :: n,n1,n2
      TYPE(Element_t), POINTER :: Edge, LeftParent, RightParent
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: EdgeBasis(n), EdgedBasisdx(n,3), EdgeddBasisddx(n,3,3)
      REAL(KIND=dp) :: LeftBasis(n1), LeftdBasisdx(n1,3), LeftddBasisddx(n1,3,3)
      REAL(KIND=dp) :: RightBasis(n2), RightdBasisdx(n2,3), RightddBasisddx(n2,3,3)
      REAL(KIND=dp) :: Jump(n1+n2), Average(n1+n2)
      REAL(KIND=dp) :: detJ, U, V, W, S, Udotn, xx, yy
      LOGICAL :: Stat
      INTEGER :: i, j, p, q, dim, t, nEdge, nParent
      TYPE(GaussIntegrationPoints_t) :: IntegStuff
      REAL(KIND=dp) :: hE, Normal(3), cu(3), LeftOut(3)

      TYPE(Nodes_t) :: EdgeNodes, LeftParentNodes, RightParentNodes

      Save EdgeNodes, LeftParentNodes, RightParentNodes
!------------------------------------------------------------------------------
      dim = CoordinateSystemDimension()
      STIFF = 0.0d0

      CALL GetElementNodes( EdgeNodes, Edge )
      CALL GetElementNodes( LeftParentNodes,  LeftParent )
      CALL GetElementNodes( RightParentNodes, RightParent )
!------------------------------------------------------------------------------
!     Numerical integration over the edge
!------------------------------------------------------------------------------
      IntegStuff = GaussPoints( Edge )

      LeftOut(1) = SUM( LeftParentNodes % x(1:n1) ) / n1
      LeftOut(2) = SUM( LeftParentNodes % y(1:n1) ) / n1
      LeftOut(3) = SUM( LeftParentNodes % z(1:n1) ) / n1
      LeftOut(1) = SUM( EdgeNodes % x(1:n) ) / n - LeftOut(1)
      LeftOut(2) = SUM( EdgeNodes % y(1:n) ) / n - LeftOut(2)
      LeftOut(3) = SUM( EdgeNodes % z(1:n) ) / n - LeftOut(3)

      DO t=1,IntegStuff % n
        U = IntegStuff % u(t)
        V = IntegStuff % v(t)
        W = IntegStuff % w(t)
        S = IntegStuff % s(t)

        ! Basis function values & derivatives at the integration point:
        !--------------------------------------------------------------
        stat = ElementInfo( Edge, EdgeNodes, U, V, W, detJ, &
             EdgeBasis, EdgedBasisdx, EdgeddBasisddx, .FALSE. )

        S = S * detJ

        Normal = NormalVector( Edge, EdgeNodes, U, V, .FALSE. )
        IF ( SUM( LeftOut*Normal ) < 0 ) Normal = -Normal

        ! Find basis functions for the parent elements:
        ! ---------------------------------------------
        CALL FindParentUVW( Edge,n,LeftParent,n1,U,V,W,EdgeBasis )
        stat = ElementInfo( LeftParent, LeftParentNodes, U, V, W, detJ, &
                LeftBasis, LeftdBasisdx, LeftddBasisddx, .FALSE. )

        CALL FindParentUVW( Edge,n,RightParent,n2,U,V,W,EdgeBasis )
        stat = ElementInfo( RightParent, RightParentNodes, U, V, W, detJ, &
              RightBasis, RightdBasisdx, RightddBasisddx, .FALSE. )

        ! Integrate jump terms:
        ! ---------------------
        Jump(1:n1) = LeftBasis(1:n1)
        Jump(n1+1:n1+n2) = -RightBasis(1:n2)

        Average(1:n1) = LeftBasis(1:n1) / 2
        Average(n1+1:n1+n2) = RightBasis(1:n2) / 2

        cu = 0.0d0
        DO i=1,dim
          cu(i) = SUM( (Velo(i,1:n)-MeshVelo(i,1:n)) * EdgeBasis(1:n) )
        END DO
        Udotn = SUM( Normal * cu )

        DO p=1,n1+n2
          DO q=1,n1+n2
            STIFF(p,q) = STIFF(p,q) + s * Udotn * Average(q) * Jump(p)
            STIFF(p,q) = STIFF(p,q) + s * ABS(Udotn)/2 * Jump(q) * Jump(p)
          END DO
        END DO
     END DO
!------------------------------------------------------------------------------
      END SUBROUTINE LocalJumps
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
   SUBROUTINE LocalMatrixBoundary( STIFF, FORCE, LOAD, &
        Element, n, ParentElement, np, Velo,MeshVelo, InFlowBC )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: STIFF(:,:),  FORCE(:), LOAD(:), Velo(:,:),MeshVelo(:,:)
     INTEGER :: n, np
     LOGICAL :: InFlowBC
     TYPE(Element_t), POINTER :: Element, ParentElement
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3)
     REAL(KIND=dp) :: ParentBasis(np), ParentdBasisdx(np,3), ParentddBasisddx(np,3,3)
     INTEGER :: i,j,p,q,t,dim

     REAL(KIND=dp) :: Normal(3), g, L, Udotn, UdotnA, cu(3), detJ,U,V,W,S
     LOGICAL :: Stat
     TYPE(GaussIntegrationPoints_t) :: IntegStuff

     TYPE(Nodes_t) :: Nodes, ParentNodes
     SAVE Nodes, ParentNodes
!------------------------------------------------------------------------------
     dim = CoordinateSystemDimension()
     FORCE = 0.0d0
     STIFF = 0.0d0

     CALL GetElementNodes( Nodes, Element )
     CALL GetElementNodes( ParentNodes, ParentElement )

     ! Numerical integration:
     !-----------------------
     IntegStuff = GaussPoints( Element )
!
! Compute the average velocity.dot.Normal        
!        
     UdotnA = 0.0   
     DO t=1,IntegStuff % n
       U = IntegStuff % u(t)
       V = IntegStuff % v(t)
       W = IntegStuff % w(t)
       S = IntegStuff % s(t)

       Normal = NormalVector( Element, Nodes, U, V, .TRUE. ) 

       ! Basis function values & derivatives at the integration point:
       ! -------------------------------------------------------------
       stat = ElementInfo( Element, Nodes, U, V, W, detJ, &
               Basis, dBasisdx, ddBasisddx, .FALSE. )
       S = S * detJ
       cu = 0.0d0
       DO i=1,dim
          cu(i) = SUM( (Velo(i,1:n)-MeshVelo(i,1:n)) * Basis(1:n) )
       END DO
       UdotnA = UdotnA + s*SUM( Normal * cu )

     END DO

     DO t=1,IntegStuff % n
       U = IntegStuff % u(t)
       V = IntegStuff % v(t)
       W = IntegStuff % w(t)
       S = IntegStuff % s(t)

       Normal = NormalVector( Element, Nodes, U, V, .TRUE. ) 

       ! Basis function values & derivatives at the integration point:
       ! -------------------------------------------------------------
       stat = ElementInfo( Element, Nodes, U, V, W, detJ, &
               Basis, dBasisdx, ddBasisddx, .FALSE. )
       S = S * detJ

       CALL FindParentUVW( Element, n, ParentElement, np, U, V, W, Basis )
       stat = ElementInfo( ParentElement, ParentNodes, U, V, W, &
           detJ,  ParentBasis, ParentdBasisdx, ParentddBasisddx, .FALSE. )

       L = SUM( LOAD(1:n) * Basis(1:n) )
       cu = 0.0d0
       DO i=1,dim
          cu(i) = SUM( (Velo(i,1:n)-MeshVelo(i,1:n)) * Basis(1:n) )
       END DO
       Udotn = SUM( Normal * cu )

       DO p = 1,np
        IF (InFlowBC .And. (UdotnA < 0.) ) THEN
            FORCE(p) = FORCE(p) - s * Udotn*L*ParentBasis(p)
         ELSE
           DO q=1,np
             STIFF(p,q) = STIFF(p,q) + s*Udotn*ParentBasis(q)*ParentBasis(p)
           END DO
         END IF
       END DO
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE LocalMatrixBoundary
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
 END SUBROUTINE FabricSolverSpectral
!------------------------------------------------------------------------------
