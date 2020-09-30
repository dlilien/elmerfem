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
           TimeVar

      REAL(KIND=dp), POINTER :: &
           FlowValues(:), EigenFabricValues(:), &
           MeshVeloValues(:), Solution(:)

      INTEGER, POINTER :: Perm(:),NodeIndexes(:), &
           FlowPerm(:), MeshVeloPerm(:), EigenFabricPerm(:)

      INTEGER :: body_id,bf_id,eq_id, Indexes(128),SpectralOrder,&
                 old_body = -1, NewtonIter,NonlinearIter, Node,&
                 dim,n1,n2,i,j,k,l,n,nd,t,iter,NDeg,STDOFs,LocalNodes,istat,&
                 comp

      REAL(KIND=dp) :: rho, lambda, A1plusA2, Bu, Bv, Bw, &
           a2(6), ai(3), Angle(3), RM(3,3), &
           SaveTime = -1, RelativeChange,UNorm,PrevUNorm,Gravity(3), &
           Tdiff,Normal(3),NewtonTol,NonlinearTol,s,OverlapMatrix

      REAL(KIND=dp), parameter :: Rad2deg=180._dp/Pi

      CHARACTER(LEN=MAX_NAME_LEN) :: SolverName='SpectralFabric',&
                                     FlowName, ComponentName, &
                                     FabricName, OverlapMatrixFile

      LOGICAL :: GotForceBC,GotIt,NewtonLinearization=.FALSE.,&
                 UnFoundFatal=.TRUE.,AllocationsDone = .FALSE.,&
                 FreeSurface,FirstTime=.TRUE.

      REAL(KIND=dp), ALLOCATABLE:: MASS(:,:), STIFF(:,:), LOAD(:,:),Force(:), &
          Alpha(:,:),Beta(:),Velocity(:,:), MeshVelocity(:,:), LocalFabric(:)

      SAVE MASS, STIFF, LOAD, Force,ElementNodes,Alpha,Beta, & 
           AllocationsDone, rho, lambda, Velocity, &
           MeshVelocity, old_body, dim, comp, SolverName, &
           SpectralOrder, LocalFabric, &
           OverlapMatrix

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

!------------------------------------------------------------------------------
!  Read constants from constants section of SIF file
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
!    Get variables needed for solution
!    All of these can be changed manually except for Mesh Velocity
!------------------------------------------------------------------------------
      IF ( .NOT. ASSOCIATED( Solver % Matrix ) ) RETURN

      Solution => Solver % Variable % Values
      Perm => Solver % Variable % Perm
      STDOFs   =  Solver % Variable % DOFs
      FabricName   =  Solver % Variable % Name

      SolverParams => GetSolverParams()

      SpectralOrder = ListGetInteger( SolverParams,'Fabric Order', &
          GotIt,UnFoundFatal=.TRUE. )
      IF (STDOFs .NE. (SpectralOrder * (SpectralOrder + 1)) / 2) THEN
          WRITE(Message,'(A,I4,A,I4)') 'Fabric variable order ', &
              STDOFs, ' does not match order ', SpectralOrder
          CALL FATAL(SolverName, message)
      END IF
      WRITE(Message,'(A,I4,A,I4)') 'Fabric dimension is', &
            STDOFs, ', matching order ', SpectralOrder
      CALL INFO(SolverName, Message , level = 20)

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

      MeshVeloVariable => VariableGet( Solver % Mesh % Variables, &
            'Mesh Velocity' )
      IF ( ASSOCIATED( MeshVeloVariable ) ) THEN
        MeshVeloPerm    => MeshVeloVariable % Perm
        MeshVeloValues  => MeshVeloVariable % Values
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
      Unorm = SQRT( SUM( Solution**2 ) / SIZE(Solution) )
!------------------------------------------------------------------------------
!     Allocate some permanent storage, this is done first time only
!------------------------------------------------------------------------------
      IF ( .NOT. AllocationsDone .OR. Solver % Mesh % Changed) THEN
        N = Model % MaxElementNodes

       dim = CoordinateSystemDimension()
       
       IF ( AllocationsDone ) THEN
         DEALLOCATE( Force, &
                     Velocity,MeshVelocity, &
                     MASS,STIFF,      &
                     LOAD, Alpha, Beta, &
                     LocalFabric)
       END IF

       ALLOCATE( Force( 2*STDOFs*N ), &
                 Velocity(4, N ),MeshVelocity(3,N), &
                 MASS( 2*STDOFs*N,2*STDOFs*N ),  &
                 STIFF( 2*STDOFs*N,2*STDOFs*N ),  &
                 LOAD( 4,N ), Alpha( 3,N ), Beta( N ), &
                 LocalFabric(N * STDOFs), &
                 STAT=istat )


       IF ( istat /= 0 ) THEN
         CALL FATAL(SolverName, 'Memory allocation error.' )
       END IF

       AllocationsDone = .TRUE.
      END IF

      IF( TransientSimulation ) THEN
        TimeVar => VariableGet( Solver % Mesh % Variables, 'Time' )
        IF ( SaveTime /= TimeVar % Values(1) ) THEN
           SaveTime = TimeVar % Values(1)
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
       
        CALL DefaultInitialize()

!------------------------------------------------------------------------------
        DO t=1,Solver % NumberOFActiveElements
!------------------------------------------------------------------------------

            IF ( RealTime() - at0 > 1.0 ) THEN
              WRITE(Message,'(a,i3,a)' ) '   Assembly: ', INT(100.0 - 100.0 * &
               (Solver % NumberOfActiveElements-t) / &
                  (1.0*Solver % NumberOfActiveElements)), ' % done'
                           
              CALL Info(SolverName, Message, Level=5 )
              at0 = RealTime()
            END IF

         CurrentElement => GetActiveElement(t)
         CALL GetElementNodes( ElementNodes )
         n = GetElementNOFNodes()
         nd = GetElementDOFs( Indexes )
         NodeIndexes => CurrentElement % NodeIndexes
         Material => GetMaterial()
         body_id = CurrentElement % BodyId


!------------------------------------------------------------------------------
!        Get element local stiffness & mass matrices
!------------------------------------------------------------------------------
         DO i=1,STDOFs
            LocalFabric(((i - 1) * n + 1):i * n) = Solution(STDOFs*(&
                Solver % Variable % Perm(Indexes(1:n))-1) + i )
         END DO

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

         CALL LocalMatrix(MASS, STIFF, FORCE, LOAD, LocalFabric, Velocity, &
           MeshVelocity, CurrentElement, n, nd, ElementNodes, rho, &
           lambda, SpectralOrder, STDOFs, OverlapMatrix)

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
!
!3D => Edges => Faces
      IF (dim.eq.3) THEN 
        CALL FATAL(SolverName, "This solver cannot actually handle 3d...")
      ELSE
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
            CALL LocalJumps( STIFF,Edge,n,STDOFs,LeftParent,n1,RightParent,n2,Velocity,MeshVelocity )
            IF ( TransientSimulation )  CALL Default1stOrderTime(MASS, STIFF, FORCE)
            CALL DefaultUpdateEquations( STIFF, FORCE, Edge )
         END IF
       END DO

      END IF

      CALL DefaultFinishBulkAssembly()

!------------------------------------------------------------------------------
!     Loop over the boundary elements
!         Cannot use DefaultDirichletBCs because of advection effects
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
        MeshVelocity=0._dp
        IF ( ASSOCIATED( MeshVeloVariable ) ) THEN
          k = MeshVeloVariable % DOFs
          DO i=1,k
            MeshVelocity(i,1:n) = MeshVeloValues(k*(MeshVeloPerm(Element % NodeIndexes)-1)+i)
          End do
        END IF

         BC => GetBC()
         MASS = 0.0d0
         LOAD = 0.0d0
         GotIt = .FALSE.

         ! Does not necessarily need to be found...
         IF ( ASSOCIATED(BC) ) THEN
           ! Different if we only have one component
           IF (STDOFs.EQ.1) THEN
             LOAD(1,1:n) = GetReal( BC, FabricName, GotIt )
           ELSE
             DO i=1,STDOFs
               LOAD(1,i:i + n * STDOFs:STDOFs ) = GetReal( BC, ComponentName(FabricName, i),&
                                                           GotIt )
             END DO
           END IF
         END IF

         MASS = 0.0d0
         CALL LocalMatrixBoundary(  STIFF, FORCE, LOAD(1,1:n * STDOFs), &
                              Element, n, STDOFs, ParentElement, n1, Velocity,MeshVelocity, GotIt )


      IF ( TransientSimulation )  CALL Default1stOrderTime(MASS, STIFF, FORCE)
        CALL DefaultUpdateEquations( STIFF, FORCE )
      END DO
      CALL DefaultFinishAssembly()

!------------------------------------------------------------------------------
      CALL Info( SolverName, 'Set boundaries done', Level=4 )

!------------------------------------------------------------------------------
!     Solve the system and check for convergence
!------------------------------------------------------------------------------
      Unorm = DefaultSolve()

      WRITE(Message,*) 'solve done', minval( solver % variable % values), maxval( Solver % variable % values)
      CALL Info( SolverName, Message, Level=4 )
      
      ! Need to post-process every node to make sure we are normalized
      DO i=1,Solver % Mesh % NumberOfNodes
        Node = Perm(i)
        CALL PostProcessFabric(Solution(STDOFS*(Node-1) + 1:STDOFS*Node), STDOFs)
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
      SUBROUTINE LocalMatrix(MASS, STIFF, FORCE, LOAD, NodalFabric, & 
          NodalVelo, NodMeshVel, Element, n, nd, Nodes, rho, lambda, SpectralOrder,&
          SpectralDim, OverlapMatrix)
!------------------------------------------------------------------------------
! Inputs and Outputs
!------------------------------------------------------------------------------
      REAL(KIND=dp), Target, INTENT(INOUT) :: STIFF(:,:),MASS(:,:)
      REAL(KIND=dp) :: LOAD(:,:), NodalVelo(:,:),NodMeshVel(:,:)
      REAL(KIND=dp), DIMENSION(:) :: FORCE, NodalFabric

      TYPE(Nodes_t) :: Nodes
      TYPE(Element_t) :: Element
      INTEGER :: n, nd, Comp
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

      REAL(KIND=dp) :: s,u,v,w, Radius, B(6,3), G(3,6), C44,C55,C66
      REAL(KIND=dp) :: Velo(3),DStress(6),StrainR(6),Spin(3),SD(6)

      REAL(KIND=dp) :: LGrad(3,3),StrainRate(3,3),D(6),angle(3),epsi
      REAL(KIND=dp) :: ap(3),C(6,6),Spin1(3,3),&
        ThisNodeFabric(SpectralDim)
      LOGICAL :: CSymmetry, Fond
              
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
      Fond=.False.
      
      hmax = maxval (Nodes % y(1:n))
        
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
        LGrad = MATMUL( NodalVelo(1:3,1:n), dBasisdx(1:n,1:3) )
        StrainRate = 0.5 * ( LGrad + TRANSPOSE(LGrad) )
        Spin1 = 0.5 * ( LGrad - TRANSPOSE(LGrad) )

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


        !    Plug in Nicholas's model
        CALL SpectralModel(STDOFs, SpectralOrder,ThisNodeFabric, StrainRate, Spin1, &
                          OverlapMatrix, DCDt)

        ! I think we need two extra loops here, since we need to handle the
        ! extra dimensions. And we still need to do all the basis functions
        ! at each point too.
        DO p=1,nd
          DO q=1,nd
            i = SpectralDim * (p -1) + 1
            j = SpectralDim * (q -1) + 1
            LocalMass => MASS(i:i + SpectralDim, j:j + SpectralDim)
            LocalStiff => Stiff(i:i + SpectralDim, j:j + SpectralDim)

            ! Loop over basis functions (of both unknowns and weights):
            ! ---------------------------------------------------------
            DO i=1,SpectralDim
              ! Advection terms:
              ! ----------------
              A = 0.0_dp
              DO k=1,dim
                A = A - Velo(k) * Basis(q) * dBasisdx(p,k)
              END DO
              LocalStiff( i,i ) = LocalStiff( i,i ) + s * A
              M = Basis(p) * Basis(q)
              LocalMass( i,i )  = LocalMass( i,i )  + s * M

              DO j=1,SpectralDim
                A = 0.0_dp

                ! Reaction term:
                ! ----------------
                A = A - DCDt(i, j) * Basis(q) * Basis(p)

                ! Add nodal matrix to element matrix:
                !  -----------------------------------
                LocalStiff( i,j ) = LocalStiff( i,j ) + s * A
              END DO
            END DO
          END DO


        ! The righthand side...is zero:
        ! ----------------------
        LoadAtIp = 0.0_dp
        LoadAtIp= LoadAtIp * Basis(p)
        FORCE(p) = FORCE(p) + s*LoadAtIp
        END DO
!-----------------------------------------------------------------------------
      END DO  ! Of integration points
!------------------------------------------------------------------------------

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
      SUBROUTINE LocalJumps( STIFF,Edge,n,SpectralDim,LeftParent,n1,RightParent,n2,Velo,MeshVelo )
!------------------------------------------------------------------------------
      REAL(KIND=dp), TARGET :: STIFF(:,:)
      REAL(KIND=dp) :: Velo(:,:),MeshVelo(:,:)
      INTEGER :: n,n1,n2,SpectralDim
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
      REAL(KIND=dp), POINTER :: LocalStiff(:, :)

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

        DO i = 1,SpectralDim
          LocalStiff => STIFF(i:i + (n1 + n2) * SpectralDim:SpectralDim,&
                              i:i + (n1 + n2) * SpectralDim:SpectralDim)
          DO p=1,n1+n2
            DO q=1,n1+n2
              STIFF(p,q) = STIFF(p,q) + s * Udotn * Average(q) * Jump(p)
              STIFF(p,q) = STIFF(p,q) + s * ABS(Udotn)/2 * Jump(q) * Jump(p)
            END DO
          END DO
        END DO
      END DO
!------------------------------------------------------------------------------
      END SUBROUTINE LocalJumps
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
      SUBROUTINE LocalMatrixBoundary( STIFF, FORCE, LOAD, &
          Element, n, SpectralDim, ParentElement, np, Velo,MeshVelo, InFlowBC )
!------------------------------------------------------------------------------
      REAL(KIND=dp), TARGET :: STIFF(:,:),  FORCE(:)
      REAL(KIND=dp) :: LOAD(:), Velo(:,:),MeshVelo(:,:)
      INTEGER :: n, np, spectraldim
      LOGICAL :: InFlowBC
      TYPE(Element_t), POINTER :: Element, ParentElement
!------------------------------------------------------------------------------
      REAL(KIND=dp) :: Basis(n), dBasisdx(n,3), ddBasisddx(n,3,3)
      REAL(KIND=dp) :: ParentBasis(np), ParentdBasisdx(np,3), ParentddBasisddx(np,3,3)
      REAL(KIND=dp) :: Normal(3), g, L, Udotn, UdotnA, cu(3), detJ,S

      INTEGER :: i,j,p,q,t,dim

      LOGICAL :: Stat

      TYPE(GaussIntegrationPoints_t) :: IP
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
      IP = GaussPoints( Element )

      ! Compute the average velocity.dot.Normal        
      UdotnA = 0.0   
      DO t=1,IP % n
        S = IP % s(t)

        Normal = NormalVector( Element, Nodes, IP % U(t), IP % V(t), .TRUE. ) 

        ! Basis function values & derivatives at the integration point:
        ! -------------------------------------------------------------
        stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), IP % W(t), detJ, &
               Basis, dBasisdx, ddBasisddx, .FALSE. )
        S = S * detJ
        cu = 0.0d0
        DO i=1,dim
          cu(i) = SUM( (Velo(i,1:n)-MeshVelo(i,1:n)) * Basis(1:n) )
        END DO
        UdotnA = UdotnA + s*SUM( Normal * cu )
      END DO

      ! Now switch to the actual integration
      DO t=1,IP % n
        S = IP % s(t)

        Normal = NormalVector( Element, Nodes, IP % U(t), IP % V(t), .TRUE. ) 

        ! Basis function values & derivatives at the integration point:
        ! -------------------------------------------------------------
        stat = ElementInfo( Element, Nodes, IP % U(t), IP % V(t), IP % W(t), detJ, &
                Basis, dBasisdx, ddBasisddx, .FALSE. )
        S = S * detJ

        CALL FindParentUVW( Element, n, ParentElement, np, IP % U(t), IP % V(t), IP % W(t), Basis )
        stat = ElementInfo( ParentElement, ParentNodes, &
            IP % U(t), IP % V(t), IP % W(t), &
            detJ,  ParentBasis, ParentdBasisdx, ParentddBasisddx, .FALSE. )
        
        ! Calculate the normal component of the velocity
        cu = 0.0d0
        DO i=1,dim
          cu(i) = SUM( (Velo(i,1:n)-MeshVelo(i,1:n)) * Basis(1:n) )
        END DO
        Udotn = SUM( Normal * cu )

        ! Now go through each component
        DO p = 1,np
          DO j=1,SpectralDim
            L = SUM( LOAD(j:j + n * SpectralDim:SpectralDim ) * Basis(1:n) )
            IF (InFlowBC .And. (UdotnA < 0.) ) THEN
              FORCE(SpectralDim * (p - 1) + j) = FORCE(SpectralDim * (p - 1) + j) - s * Udotn*L*ParentBasis(p)
            ELSE
              DO q=1,np
                STIFF(SpectralDim * (p - 1) + j, SpectralDim * (q - 1) + j) = &
                    STIFF(SpectralDim * (p - 1) + j, SpectralDim * (q - 1) + j) &
                    + s*Udotn*ParentBasis(q)*ParentBasis(p)
              END DO
            END IF
          END DO
        END DO
      END DO
!------------------------------------------------------------------------------
      END SUBROUTINE LocalMatrixBoundary
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
 END SUBROUTINE SpectralFabricSolver
!------------------------------------------------------------------------------
