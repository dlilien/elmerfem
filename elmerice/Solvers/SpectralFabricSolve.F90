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
! *  Authors:  David A. Lilien
! *  Email:   dlilien90@gmail.com
! *  Based on FabricSolve.F90, by Juha Ruokolainen, Olivier Gagliardini,
! *      and Fabien Gillet-Chaulet
! *  Web:     http://elmerice.elmerfem.org
! *       Date of modification: 09/20
! *
! *****************************************************************************/
!>  Solver for fabric parameter equations 
!------------------------------------------------------------------------------
RECURSIVE SUBROUTINE SpectralFabricSolver( Model,Solver,dt,TransientSimulation )
!------------------------------------------------------------------------------

      USE DefUtils
      USE specfab

      IMPLICIT NONE
!------------------------------------------------------------------------------
!******************************************************************************
!
!  Solve Fabric equations at one timestep
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
      TYPE(Nodes_t) :: ElementNodes
      TYPE(Solver_t), POINTER :: PSolver
      TYPE(Matrix_t), POINTER :: StiffMatrix
      TYPE(ValueList_t), POINTER :: Material, BC, SolverParams
      TYPE(Element_t), POINTER :: CurrentElement, Element, &
           ParentElement, LeftParent, RightParent, Edge

      TYPE(Variable_t), POINTER :: &
           FlowVariable, EigenFabricVariable, MeshVeloVariable, &
           TimeVar, FabricVariable, C2Variable

      REAL(KIND=dp), POINTER :: &
           FlowValues(:), EigenFabricValues(:), &
           MeshVeloValues(:), FabricValues(:), &
           PrevFabric(:),CurrFabric(:),TempFabVal(:),&
           Solution(:), Ref(:), C2Values(:)

      INTEGER, POINTER :: NodeIndexes(:), &
           FlowPerm(:), MeshVeloPerm(:), C2Perm(:), EigenFabricPerm(:), FabricPerm(:)

      INTEGER :: body_id,bf_id,eq_id, Indexes(128),SpectralOrder,&
                 old_body = -1, NewtonIter,NonlinearIter, Node,&
                 dim,n1,n2,i,j,k,l,n,nd,t,iter,NDeg,STDOFs,LocalNodes,istat,&
                 comp, SpectralDim,INDi(6),INDj(6)
      REAL(KIND=dp) :: rho, lambda, A1plusA2, Bu, Bv, Bw, &
           a2(6), ai(3), Angle(3), RM(3,3), &
           SaveTime = -1, RelativeChange,UNorm,PrevUNorm,Gravity(3), &
           Tdiff,Normal(3),NewtonTol,NonlinearTol,s,OverlapMatrix,&
           C2(3,3)

      REAL(KIND=dp), parameter :: Rad2deg=180._dp/Pi

      CHARACTER(LEN=MAX_NAME_LEN) :: SolverName='SpectralFabric',&
                                     FlowName, ComponentName, &
                                     FabricName, OverlapMatrixFile,&
                                     C2Name, EigName

      LOGICAL :: GotForceBC,GotIt,NewtonLinearization=.FALSE.,&
                 UnFoundFatal=.TRUE.,AllocationsDone = .FALSE.,&
                 FreeSurface,ExportFabric=.False.,FirstTime=.TRUE.,&
                 ExportEigV=.False.

      REAL(KIND=dp), ALLOCATABLE:: MASS(:,:), STIFF(:,:), LOAD(:,:),Force(:), &
          Alpha(:,:),Beta(:),Velocity(:,:), MeshVelocity(:,:), LocalFabric(:)

      COMPLEX(KIND=dp), ALLOCATABLE:: FabOut(:)

      SAVE MASS, STIFF, LOAD, Force,ElementNodes,Alpha,Beta, & 
           AllocationsDone, rho, lambda, Velocity, &
           MeshVelocity, old_body, dim, comp, SolverName, &
           CurrFabric, TempFabVal, PrevFabric, &
           SpectralOrder, OverlapMatrix, LocalFabric, &
           ExportEigV, ExportFabric, FabOut

#ifdef USE_ISO_C_BINDINGS
      REAL(KIND=dp) :: at, at0
#else
      REAL(KIND=dp) :: at, at0, CPUTime, RealTime
#endif

      INTERFACE
        SUBROUTINE PostProcessFabric(C, ProbDim)
            INTEGER, INTENT(IN) :: ProbDim
            REAL(kind=8), INTENT(INOUT):: C(ProbDim)
        END SUBROUTINE PostProcessFabric
      END INTERFACE

      INTERFACE
        Subroutine R2Ro_mat(A,ai,angle)
        USE Types
        REAL(KIND=dp),intent(in) :: A(3,3)
        REAL(KIND=dp),intent(out) :: ai(3), Angle(3)
       End Subroutine R2Ro_mat
      End Interface  
!------------------------------------------------------------------------------
!  Read constants from constants section of SIF file
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!    Get variables needed for solution
!    All of these can be changed manually except for Mesh Velocity
!------------------------------------------------------------------------------
      IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN

      Solution => Solver % Variable % Values
      STDOFs   =  Solver % Variable % DOFs

      SolverParams => GetSolverParams()

      SpectralOrder = ListGetInteger( SolverParams,'Fabric Order', &
          GotIt,UnFoundFatal=.TRUE. )
      ! SpectralDim = SpectralOrder * (SpectralOrder + 1) / 2
      SpectralDim = sum([(1+i*2, i=0, SpectralOrder,2)])


      FlowName = ListGetString( SolverParams, 'Flow Solution Name', &
            GotIt,UnFoundFatal=.FALSE. )
      IF (.NOT.GotIt) THEN
          FlowName = 'AIFlow'
      END IF
      FlowVariable => VariableGet( Solver % Mesh % Variables, FlowName)
      IF ( ASSOCIATED( FlowVariable ) ) THEN
        FlowPerm    => FlowVariable % Perm
        FlowValues  => FlowVariable % Values
        WRITE(Message,'(A,A)') 'Flow variable = ', FlowName
        CALL INFO(SolverName, Message , level = 20)
      ELSE
        WRITE(Message,'(A,A)') 'Flow variable not found:', FlowName
        CALL FATAL(SolverName, message)
      END IF

      FabricName = ListGetString( SolverParams, 'Fabric Name', &
            GotIt,UnFoundFatal=.FALSE. )
      IF (.NOT.GotIt) THEN
          FabricName = 'SpectralFabric'
      END IF
      FabricVariable => VariableGet( Solver % Mesh % Variables, FabricName)
      IF ( ASSOCIATED( FabricVariable ) ) THEN
        FabricPerm    => FabricVariable % Perm
        FabricValues  => FabricVariable % Values
        WRITE(Message,'(A,A)') 'Fabric variable = ', FabricName
        CALL INFO(SolverName, Message , level = 20)
      ELSE
        WRITE(Message,'(A,A)') 'Fabric variable not found:', FabricName
        CALL FATAL(SolverName, message)
      END IF
      IF (SpectralDim .NE. FabricVariable % DOFs) THEN
          WRITE(Message,*) 'Fabric DOF mismatch', SpectralDim, 'is not', FabricVariable % DOFs
        CALL FATAL(SolverName, message)
      END IF

      MeshVeloVariable => VariableGet( Solver % Mesh % Variables, &
            'Mesh Velocity' )
      IF ( ASSOCIATED( MeshVeloVariable ) ) THEN
        MeshVeloPerm    => MeshVeloVariable % Perm
        MeshVeloValues  => MeshVeloVariable % Values
      END IF

      C2Name = ListGetString( SolverParams, 'C2 Name', &
            GotIt,UnFoundFatal=.FALSE. )
      IF (GotIt) THEN
        ExportFabric = .TRUE.
        C2Variable => VariableGet( Solver % Mesh % Variables, C2Name)
        IF ( ASSOCIATED( C2Variable ) ) THEN
          C2Perm    => C2Variable % Perm
          C2Values  => C2Variable % Values
          WRITE(Message,'(A,A)') 'C2 variable = ', C2Name
          CALL INFO(SolverName, Message , level = 20)
        ELSE
          WRITE(Message,'(A,A)') 'C2 variable called for but not found:', C2Name
          CALL FATAL(SolverName, message)
        END IF
      END IF

      EigName = ListGetString( SolverParams, 'EigenV Name', &
            GotIt,UnFoundFatal=.FALSE. )
      IF (GotIt) THEN
        ExportEigV = .TRUE.
        EigenFabricVariable => VariableGet( Solver % Mesh % Variables, EigName)
        IF ( ASSOCIATED( C2Variable ) ) THEN
          EigenFabricPerm    => EigenFabricVariable % Perm
          EigenFabricValues  => EigenFabricVariable % Values
          WRITE(Message,'(A,A)') 'EigenV variable = ', EigName
          CALL INFO(SolverName, Message , level = 20)
        ELSE
          WRITE(Message,'(A,A)') 'EigenV variable called for but not found:', EigName
          CALL FATAL(SolverName, message)
        END IF
      END IF

      !!!!!!!!!!!!!
      !! Must be handled later
      !!!!!!!!!!!!!
      IF (FirstTime) THEN
        FirstTime = .FALSE.
        OverlapMatrix = 0.0
      END IF
                                       
!------------------------------------------------------------------------------
      StiffMatrix => Solver % Matrix
      Unorm = SQRT( SUM( FabricValues**2 ) / SIZE(FabricValues) )
!------------------------------------------------------------------------------
!     Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
      IF ( .NOT. AllocationsDone .OR. Solver % Mesh % Changed) THEN
        N = Model % MaxElementNodes

       dim = CoordinateSystemDimension()
       
       IF ( AllocationsDone ) THEN
         DEALLOCATE( Force, &
                     Velocity,MeshVelocity, &
                     MASS,STIFF,&
                     LOAD, Alpha, Beta,&
                     LocalFabric,&
                     CurrFabric, TempFabVal,&
                     FabOut )
       END IF

       ALLOCATE( Force( 2*STDOFs*N ), &
                 Velocity(4, N ),MeshVelocity(3,N), &
                 MASS( 2*STDOFs*N,2*STDOFs*N ),  &
                 STIFF( 2*STDOFs*N,2*STDOFs*N ),  &
                 LOAD( 4,N ), Alpha( 3,N ), Beta( N ), &
                 CurrFabric( SpectralDim*SIZE(Solver % Variable % Values)), &
                 TempFabVal( SIZE(FabricValues)), &
                 LocalFabric(N * SpectralDim), &	
                 FabOut(SpectralDim), &	
                 STAT=istat )


       IF ( istat /= 0 ) THEN
         CALL FATAL(SolverName, 'Memory allocation error.' )
       END IF

       CurrFabric = 0.
       TempFabVal = 0.
       IF ( TransientSimulation ) THEN
          IF (AllocationsDone ) DEALLOCATE (PrevFabric)
          ALLOCATE( PrevFabric( SpectralDim*SIZE(Solver % Variable % Values)) )
         PrevFabric = 0.
       END IF

       DO i=1,Solver % NumberOFActiveElements
          CurrentElement => GetActiveElement(i)   
          n = GetElementDOFs( Indexes )
          n = GetElementNOFNodes()
          NodeIndexes => CurrentElement % NodeIndexes
          Indexes(1:n) = Solver % Variable % Perm( Indexes(1:n) )
          DO COMP=1,SpectralDim
            IF ( TransientSimulation ) THEN
               PrevFabric(SpectralDim*(Indexes(1:n)-1)+COMP) = &
                   FabricValues(SpectralDim*(FabricPerm(NodeIndexes(1:n))-1)+COMP)
            END IF
               CurrFabric(SpectralDim*(Indexes(1:n)-1)+COMP) = &
                   FabricValues(SpectralDim*(FabricPerm(NodeIndexes(1:n))-1)+COMP)
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
      NonlinearTol = ListGetConstReal( Solver % Values, &
        'Nonlinear System Convergence Tolerance' )

      NewtonTol = ListGetConstReal( Solver % Values, &
        'Nonlinear System Newton After Tolerance' )

      NewtonIter = ListGetInteger( Solver % Values, &
        'Nonlinear System Newton After Iterations' )

      NonlinearIter = ListGetInteger( Solver % Values, &
         'Nonlinear System Max Iterations',GotIt )

      IF ( .NOT.GotIt ) NonlinearIter = 1
      DO iter=1,NonlinearIter
        at  = CPUTime()
        at0 = RealTime()

        CALL Info( SolverName, ' ', Level=4 )
        CALL Info( SolverName, ' ', Level=4 )
        CALL Info( SolverName, &
                    '-------------------------------------',Level=4 )
        WRITE( Message, * ) 'Spectral Fabric solver  iteration', iter
        CALL Info( SolverName, Message,Level=4 )
        CALL Info( SolverName, &
                     '-------------------------------------',Level=4 )
        CALL Info( SolverName, ' ', Level=4 )
        CALL Info( SolverName, 'Starting assembly...',Level=4 )

        PrevUNorm = UNorm
       
        DO COMP=1,SpectralDim
            Solver % Variable % Values = CurrFabric( COMP::SpectralDim )
            IF ( TransientSimulation ) THEN
              Solver % Variable % PrevValues(:,1) = PrevFabric( COMP::SpectralDim )
            END IF

            CALL DefaultInitialize()
    !------------------------------------------------------------------------------
            DO t=1,Solver % NumberOFActiveElements
    !------------------------------------------------------------------------------

                IF ( RealTime() - at0 > 1.0 ) THEN
                  WRITE(Message,'(a,i3,a)' ) '   Assembly: ', INT(100.0 - 100.0 * &
                   (Solver % NumberOfActiveElements-t) / &
                      (1.0*Solver % NumberOfActiveElements)), ' % done'
                               
                  CALL Info(SolverName, Message, Level=6 )
                  at0 = RealTime()
                END IF

             CurrentElement => GetActiveElement(t)
             CALL GetElementNodes( ElementNodes )
             n = GetElementDOFs( Indexes )
             n = GetElementNOFNodes()
             NodeIndexes => CurrentElement % NodeIndexes
             Material => GetMaterial()
             body_id = CurrentElement % BodyId

             DO i=1,SpectralDim
               LocalFabric(i::SpectralDim) = CurrFabric(SpectralDim*(&
                 Solver % Variable % Perm(Indexes(1:n)) - 1) + i )
             END DO

    !------------------------------------------------------------------------------
    !        Get element local stiffness & mass matrices
    !------------------------------------------------------------------------------
             k = FlowVariable % DOFs
             Velocity = 0.0d0
             DO i=1,k
                Velocity(i,1:n) = FlowValues(k*(FlowPerm(NodeIndexes)-1)+i)
             END DO

             MeshVelocity=0._dp
             IF (ASSOCIATED(MeshVeloVariable)) Then
               k = MeshVeloVariable % DOFs
               DO i=1,k
                  MeshVelocity(i,1:n) = MeshVeloValues(k*(MeshVeloPerm(NodeIndexes)-1)+i)
               END DO
             EndIF

             CALL LocalMatrix(MASS, STIFF, FORCE, LOAD, Comp, LocalFabric, Velocity, &
               MeshVelocity, CurrentElement, n, ElementNodes, rho, &
               lambda, SpectralOrder, SpectralDim, OverlapMatrix)

    !------------------------------------------------------------------------------
    !        Update global matrices from local matrices 
    !------------------------------------------------------------------------------
             IF ( TransientSimulation )  CALL Default1stOrderTime(MASS,STIFF,FORCE)
             CALL DefaultUpdateEquations( STIFF, FORCE )
    !------------------------------------------------------------------------------
          END DO
          CALL Info( SolverName, 'Assembly done', Level=4 )
    !------------------------------------------------------------------------------
    !------------------------------------------------------------------------------
    !     Assembly of the edge terms
    !------------------------------------------------------------------------------
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
            STIFF = 0.0d0
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
            STIFF = 0.0d0
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
            LOAD(1,1:n) = GetReal( BC, ComponentName(FabricName, Comp) , GotIt )
         END IF

         MASS = 0.0d0
         STIFF = 0.0d0
         FORCE = 0.0d0
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
      WRITE(Message,*) 'Solve Done', minval( solver % variable % values), maxval( Solver % variable % values)
      CALL Info( SolverName, Message, Level=4 )
      
      n1 = Solver % Mesh % NumberOfNodes
      ALLOCATE( Ref(n1) )
      Ref = 0
      TempFabVal(COMP::SpectralDim ) = 0. !fab
      
      DO t=1,Solver % NumberOfActiveElements
         Element => GetActiveElement(t) 
         n = GetElementDOFs( Indexes )
         n = GetElementNOFNodes()
         
         DO i=1,n
            k = Element % NodeIndexes(i)
            TempFabVal( SpectralDim*(FabricPerm(k)-1) + COMP ) =    & 
            TempFabVal( SpectralDim*(FabricPerm(k)-1) + COMP ) + &
            Solver % Variable % Values(Solver % Variable % Perm(Indexes(i)) )
            FabricValues( SpectralDim*(FabricPerm(k)-1) + COMP ) = &
                          TempFabVal(SpectralDim*(FabricPerm(k)-1) + COMP ) 
            Ref(k) = Ref(k) + 1
         END DO
      END DO

        DO i=1,n1
          j=FabricPerm(i)
          IF (j < 1) CYCLE
          IF ( Ref(i) > 0 ) THEN
            FabricValues( SpectralDim*(j-1)+COMP ) = &
                   FabricValues( SpectralDim*(j-1)+COMP ) / Ref(i)
          END IF
        END DO
        DEALLOCATE( Ref )
        
      END DO ! End DO Comp

      INDi(1:5) = (/ 1, 3, 1, 3, 2 /)
      INDj(1:5) = (/ 1, 3, 3, 2, 1 /)

      ! Need to post-process every node to make sure we are normalized
      DO i=1,Solver % Mesh % NumberOfNodes
        IF (ExportFabric.OR.ExportEigV) THEN
          FabOut = CMPLX(FabricValues(SpectralDim*(i - 1) + 1:SpectralDim*i))
          C2 = A2_ij(FabOut)
          ! C2 = reshape((/ 0.33, 0.0, 0.0, 0.0, 0.33, 0.0, 0.0, 0.0, 0.33 /), shape(C2))

        END IF
        IF (ExportFabric) THEN
          DO j=1,5
            C2Values(5 * (C2Perm(i) - 1) + j) = C2(INDi(j), INDj(j))
          END DO
        END IF
        IF (ExportFabric) THEN
          call R2Ro_mat(C2,ai,angle)
          angle(:)=angle(:)*rad2deg
          If (angle(1).gt.90._dp) angle(1)=angle(1)-180._dp
          If (angle(1).lt.-90._dp) angle(1)=angle(1)+180._dp
          DO j=1,3
            EigenFabricValues(6 * (EigenFabricPerm(i) - 1) + j) = ai(j)
            EigenFabricValues(6 * (EigenFabricPerm(i) - 1) + 3 + j) = angle(j)
          END DO
        END IF
        CALL PostProcessFabric(FabricValues(SpectralDim*(i - 1) + 1:SpectralDim*i), SpectralDim)
      END DO

      ! And finally need to reset the current fabric
       DO i=1,Solver % NumberOFActiveElements
          CurrentElement => GetActiveElement(i)   
          n = GetElementDOFs( Indexes )
          n = GetElementNOFNodes()
          NodeIndexes => CurrentElement % NodeIndexes
          Indexes(1:n) = Solver % Variable %Perm( Indexes(1:n) )
          DO COMP=1,SpectralDim
            CurrFabric(SpectralDim*(Indexes(1:n)-1)+COMP) = &
                        FabricValues(SpectralDim*(FabricPerm(NodeIndexes(1:n))-1)+COMP)
          END DO
       END DO

      ! Normalization will have altered the norm, so recalculate and override
      Unorm = SQRT( SUM( Solution**2 ) / SIZE(Solution) )
      Solver % Variable % Norm = Unorm  

      !------------------------------------------------------------------------------
      ! Check convergence
      IF ( PrevUNorm + UNorm /= 0.0d0 ) THEN
         RelativeChange = 2.0d0 * ABS( PrevUNorm - UNorm) / ( PrevUnorm + UNorm)
      ELSE
         RelativeChange = 0.0d0
      END IF

      WRITE( Message, * ) 'Result Norm   : ',UNorm
      CALL Info( SolverName, Message, Level=4 )
      WRITE( Message, * ) 'Relative Change : ',RelativeChange
      CALL Info( SolverName, Message, Level=4 )
      IF ( RelativeChange < NewtonTol .OR. &
            iter > NewtonIter ) NewtonLinearization = .TRUE.
      IF ( RelativeChange < NonLinearTol ) EXIT
      !------------------------------------------------------------------------------

!------------------------------------------------------------------------------
      END DO ! of nonlinear iter
!------------------------------------------------------------------------------
CONTAINS

!------------------------------------------------------------------------------
      SUBROUTINE LocalMatrix(MASS, STIFF, FORCE, LOAD, Comp, LocalFabric,& 
          NodalVelo, NodMeshVel, Element, n, Nodes, rho, lambda, SpectralOrder,&
          SpectralDim, OverlapMatrix)
!------------------------------------------------------------------------------
! Inputs and Outputs
!------------------------------------------------------------------------------
      REAL(KIND=dp), Target, INTENT(INOUT) :: STIFF(:,:),MASS(:,:)
      REAL(KIND=dp) :: LOAD(:,:), NodalVelo(:,:),NodMeshVel(:,:)
      REAL(KIND=dp), DIMENSION(:) :: FORCE, LocalFabric

      TYPE(Nodes_t) :: Nodes
      TYPE(Element_t) :: Element
      INTEGER :: n, Comp
!------------------------------------------------------------------------------
!    Local variables
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: Basis(2*n),ddBasisddx(1,1,1)
      REAL(KIND=dp) :: dBasisdx(2*n,3),SqrtElementMetric

      REAL(KIND=dp) :: Theta, LoadAtIp

      REAL(KIND=dp) :: A,M, hK,tau,pe1,pe2,unorm,C0, SU(n), SW(n)
      REAL(KIND=dp) :: rho,lambda,Deq, ai(6),a4(9),hmax

      INTEGER :: i,j,k,p,q,t,dim,ind(3),&
          SpectralOrder,SPectralDim, C1, C2

      REAL(KIND=dp) :: s,u,v,w, Radius, B(6,3), G(3,6)
      REAL(KIND=dp) :: Velo(3),StrainR(6),Spin(3),SD(6)

      REAL(KIND=dp) :: LGrad(2,2),SR(2,2),StrainRate(3,3),D(6),angle(3),epsi
      REAL(KIND=dp) :: ap(3),Spin1(3,3),&
        ThisNodeFabric(SpectralDim)
      LOGICAL :: CSymmetry
              
      INTEGER :: N_Integ
      REAL(KIND=dp), DIMENSION(:), POINTER :: U_Integ,V_Integ, &
                                              W_Integ,S_Integ
      REAL(KIND=dp), POINTER :: LocalMass(:,:),LocalStiff(:,:)

      LOGICAL :: stat
      TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
      REAL(kind=dp) :: OverlapMatrix, DCDt(SpectralDim, SpectralDim)

      INTERFACE
        SUBROUTINE SpectralModel(ProbDim, SpectralOrder, C, StrainRate, Spin, &
                                 OverlapMatrix, DCDt)
            INTEGER, intent(in) :: SpectralOrder, ProbDim
            REAL(kind=8), intent(in) :: C(ProbDim), StrainRate(3, 3), Spin(3,3)
            REAL(kind=8), intent(in) :: OverlapMatrix
            REAL(kind=8), intent(out) :: DCDt(ProbDim, ProbDim)
        END SUBROUTINE SpectralModel
      END INTERFACE
!------------------------------------------------------------------------------
      
      dim = CoordinateSystemDimension()

      FORCE = 0.0D0
      MASS  = 0.0D0
      STIFF = 0.0D0
!    
!    Integration stuff:
!    ------------------
      IntegStuff = GaussPoints( Element  )

      U_Integ => IntegStuff % u
      V_Integ => IntegStuff % v
      W_Integ => IntegStuff % w
      S_Integ => IntegStuff % s
      N_Integ =  IntegStuff % n

      hk = ElementDiameter( Element, Nodes )
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

        !      Strain-Rate, Stresses and Spin
        StrainRate = 0.0
        Spin1 = 0.0

        !    Compute strainRate and Spin :
        !    -----------------------------
        LGrad = MATMUL( NodalVelo(1:2,1:n), dBasisdx(1:n,1:2) )
        SR = 0.5 * ( LGrad + TRANSPOSE(LGrad) )
        StrainRate(1,1) = SR(1, 1)
        StrainRate(3,1) = SR(2, 1)
        StrainRate(1,3) = SR(1, 2)
        StrainRate(3,3) = SR(2, 2)
        SR = 0.5 * ( LGrad - TRANSPOSE(LGrad) )
        Spin1(1,1) = SR(1, 1)
        Spin1(3,1) = SR(2, 1)
        Spin1(1,3) = SR(1, 2)
        Spin1(3,3) = SR(2, 2)

        IF ( CSymmetry ) THEN
          StrainRate(1,3) = 0.0
          StrainRate(2,3) = 0.0
          StrainRate(3,1) = 0.0
          StrainRate(3,2) = 0.0
          StrainRate(3,3) = 0.0

          Radius = SUM( Nodes % x(1:n) * Basis(1:n) )

          IF ( Radius > 10*AEPS ) THEN
            StrainRate(3,3) = SUM( Nodalvelo(1,1:n) * Basis(1:n) ) / Radius
          END IF

          epsi = StrainRate(1,1)+StrainRate(2,2)+StrainRate(3,3)
          DO i=1,3   
            StrainRate(i,i) = StrainRate(i,i) - epsi/3.0
          END DO

        ELSE
          epsi = StrainRate(1,1)+StrainRate(2,2)+StrainRate(3,3)
          DO i=1,dim 
            StrainRate(i,i) = StrainRate(i,i) - epsi/dim
          END DO

        END IF

        ! Velocity must be corrected by the mesh velocity
        Velo = 0.0d0
        DO i=1,dim
           Velo(i) = SUM( Basis(1:n) * (NodalVelo(i,1:n) - NodMeshVel(i,1:n)) )
        END DO

        ! Need the fabric vector
        ThisNodeFabric = 0.0
        DO i = 1,SpectralDim
          ThisNodeFabric(i) = SUM(Basis(1:n) * LocalFabric(i::STDOFs))
        END DO


        !    Plug in Nicholas's model
        CALL SpectralModel(SpectralDim, SpectralOrder,ThisNodeFabric, StrainRate, Spin1, &
                          OverlapMatrix, DCDt)

                      
!     Loop over basis functions (of both unknowns and weights):
!     ---------------------------------------------------------
      DO p=1,n
         DO q=1,n
            A = 0.0d0
            M = Basis(p) * Basis(q)
!
!           Reaction terms:
!           ---------------
            A = A - DCDt(Comp,Comp) * Basis(q) * Basis(p)

            !
            ! Advection terms:
            ! ----------------
            DO j=1,dim
               A = A - Velo(j) * Basis(q) * dBasisdx(p,j)
            END DO

!           Add nodal matrix to element matrix:
!           -----------------------------------
            MASS( p,q )  = MASS( p,q )  + s * M
            STIFF( p,q ) = STIFF( p,q ) + s * A
         END DO

        ! The righthand side
        ! ----------------------
        ! Need to subtract off the component that is caused
        ! by the diagonal element, since we solve for that above
        LoadAtIp = SUM(DCDt(comp, :) * ThisNodeFabric) - DCDt(comp, comp) * ThisNodeFabric(comp)

        LoadAtIp= LoadAtIp * Basis(p)
        FORCE(p) = FORCE(p) + s*LoadAtIp
        END DO
!-----------------------------------------------------------------------------
      END DO  ! Of integration points
!------------------------------------------------------------------------------

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
 END SUBROUTINE SpectralFabricSolver
!------------------------------------------------------------------------------
