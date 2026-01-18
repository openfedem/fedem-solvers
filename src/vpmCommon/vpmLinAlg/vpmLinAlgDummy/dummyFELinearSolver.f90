!! SPDX-FileCopyrightText: 2023 SAP SE
!!
!! SPDX-License-Identifier: Apache-2.0
!!
!! This file is part of FEDEM - https://openfedem.org
!!==============================================================================

module FELinearSolver

  !!============================================================================
  !! This module is a dummy implementation of the FE linear solver from DNVS.
  !! It is only used on platforms where this equation solver is not supported.
  !!============================================================================

  implicit none

  type Error_Flag
     integer :: dummy
  end type Error_Flag

  type MessageLevel
     integer :: dummy
  end type MessageLevel

  type SAM
     integer :: dummy
  end type SAM

  type CAM
     integer :: dummy
  end type CAM

  type GSFColSup
     integer :: dummy
  end type GSFColSup

  interface FEAnalyze
     module procedure FEAnalyzeA, FEAnalyzeB, FEAnalyzeC, FEAnalyzeD
  end interface

  interface FEEndAssembly
     module procedure FEEndAssemblyA, FEEndAssemblyB
  end interface

  interface FEExtractDiagA
     module procedure FEExtractDiagAA, FEExtractDiagAB
  end interface

  interface FEExtractK12
     module procedure FEExtractK12A, FEExtractK12B
  end interface

  interface FEPermuteVectors
     module procedure FEPermuteVectorsA, FEPermuteVectorsB, FEPermuteVectorsC
  end interface

  interface FEDestroy
     module procedure FEDestroyA, FEDestroyB
  end interface

  private :: FEAnalyzeA, FEAnalyzeB, FEAnalyzeC, FEAnalyzeD
  private :: FEEndAssemblyA, FEEndAssemblyB
  private :: FEExtractDiagAA, FEExtractDiagAB
  private :: FEExtractK12A, FEExtractK12B
  private :: FEPermuteVectorsA, FEPermuteVectorsB
  private :: FEDestroyA, FEDestroyB
  private :: errorMsg


contains

  subroutine errorMsg (prgnam,INFO,Msg)
    use reportErrorModule, only : error_p, reportError
    character(len=*)            , intent(in) :: prgnam
    type(Error_Flag)            , pointer    :: INFO
    type(MessageLevel), optional, intent(in) :: Msg
    call reportError (error_p, &
         'The GSF equation solver is not available in this version.', &
         'Toggle off the "Equation solver out-of-core" option in the GUI,', &
         'or remove -gsfSolver as additional option.', addString=prgnam)
    allocate(INFO)
    if (present(Msg)) then
       INFO%dummy = -max(1,abs(Msg%dummy))
    else
       INFO%dummy = -999
    end if
  end subroutine errorMsg

  subroutine SetOutputUnit (Msg,iNum)
    type(MessageLevel), intent(inout) :: Msg
    integer           , intent(in)    :: iNum
    print *,' ** SetOutputUnit dummy: ',Msg%dummy,iNum
  end subroutine SetOutputUnit

  subroutine SetMemorySize (Msg,iNum)
    type(MessageLevel), intent(inout) :: Msg
    integer           , intent(in)    :: iNum
    print *,' ** SetMemorySize dummy: ',Msg%dummy,iNum
  end subroutine SetMemorySize

  subroutine SetLinearSolverType (Msg,iNum)
    type(MessageLevel), intent(inout) :: Msg
    integer           , intent(in)    :: iNum
    print *,' ** SetLinearSolverType dummy: ',Msg%dummy,iNum
  end subroutine SetLinearSolverType

  subroutine SetScratchFileLocation (Msg,tmpDir)
    type(MessageLevel), intent(inout) :: Msg
    character(len=*)  , intent(in)    :: tmpDir
    print *,' ** SetScratchFileLocation dummy: ',Msg%dummy,tmpDir
  end subroutine SetScratchFileLocation

  subroutine WriteMessageData (Msg)
    type(MessageLevel), intent(in) :: Msg
    print *,' ** WriteMessageData dummy: ',Msg%dummy
  end subroutine WriteMessageData

  subroutine FEInitialize (Msg,INFO)
    type(MessageLevel), intent(in) :: Msg
    type(Error_Flag)  , pointer    :: INFO
    call errorMsg ('FEInitialize',INFO,Msg)
  end subroutine FEInitialize

  subroutine FEAnalyzeA (Msg,P,Q,T,INFO)
    use FEData        , only          :  FEDataInput
    type(MessageLevel), intent(inout) :: Msg
    type(FEDataInput) , intent(in)    :: P
    type(SAM)         , pointer       :: Q
    type(GSFColSup)   , pointer       :: T
    type(Error_Flag)  , pointer       :: INFO
    print *,' ** FEAnalyzeA dummy: ',associated(P%sam), &
         &                           associated(Q),associated(T)
    call errorMsg ('FEAnalyze',INFO,Msg)
  end subroutine FEAnalyzeA

  subroutine FEAnalyzeB (Msg,Q,T,INFO)
    type(MessageLevel), intent(inout) :: Msg
    type(SAM)         , intent(inout) :: Q
    type(GSFColSup)   , pointer       :: T
    type(Error_Flag)  , pointer       :: INFO
    print *,' ** FEAnalyzeB dummy: ',Q%dummy,associated(T)
    call errorMsg ('FEAnalyze',INFO,Msg)
  end subroutine FEAnalyzeB

  subroutine FEAnalyzeC (Msg,P,Q,Qsub,T,INFO)
    use FEData        , only          :  FEDataInput
    type(MessageLevel), intent(inout) :: Msg
    type(FEDataInput) , intent(in)    :: P
    type(SAM)         , intent(in)    :: Q
    type(SAM)         , pointer       :: Qsub
    type(GSFColSup)   , pointer       :: T
    type(Error_Flag)  , pointer       :: INFO
    print *,' ** FEAnalyzeC dummy: ',associated(P%sam),Q%dummy, &
         &                           associated(Qsub),associated(T)
    call errorMsg ('FEAnalyze',INFO,Msg)
  end subroutine FEAnalyzeC

  subroutine FEAnalyzeD (Msg,P,Q,Qsub,INFO)
    use FEData        , only          :  FEDataInput
    type(MessageLevel), intent(inout) :: Msg
    type(FEDataInput) , intent(in)    :: P
    type(SAM)         , intent(in)    :: Q
    type(SAM)         , pointer       :: Qsub
    type(Error_Flag)  , pointer       :: INFO
    print *,' ** FEAnalyzeD dummy: ',associated(P%sam),Q%dummy,associated(Qsub)
    call errorMsg ('FEAnalyze',INFO,Msg)
  end subroutine FEAnalyzeD

  function FEMatrixOrder (Q)
    type(SAM), intent(in) :: Q
    integer               :: FEMatrixOrder
    print *,' ** FEMatrixOrder dummy: ',Q%dummy
    FEMatrixOrder = 0
  end function FEMatrixOrder

  subroutine FEGetSuperelementOrder (Q,neq1,neq2)
    type(SAM), intent(in)  :: Q
    integer  , intent(out) :: neq1, neq2
    print *,' ** FEGetSuperelementOrder dummy: ',Q%dummy
    neq1 = 0
    neq2 = 0
  end subroutine FEGetSuperelementOrder

  subroutine FEGetActiveOrder (Q,MEQN,INFO)
    type(SAM)       , intent(in)  :: Q
    integer         , intent(out) :: MEQN(:)
    type(Error_Flag), pointer     :: INFO
    print *,' ** FEGetActiveOrder dummy: ',Q%dummy,associated(INFO)
    MEQN = 0
  end subroutine FEGetActiveOrder

  subroutine FEBeginAssembly (Q,S,INFO)
    type(SAM)       , intent(in) :: Q
    type(CAM)       , pointer    :: S
    type(Error_Flag), pointer    :: INFO
    print *,' ** FEBeginAssembly dummy: ',Q%dummy,associated(S)
    call errorMsg ('FEBeginAssembly',INFO)
  end subroutine FEBeginAssembly

  subroutine FEAssembleElement (OPT,iE,k,P,Q,S,NRHS,B,INFO)
    use kindModule   , only          :  dp
    use FEData       , only          :  FEDataInput
    integer          , intent(in)    :: OPT, iE
    real(dp)         , intent(in)    :: k(:,:)
    type(FEDataInput), intent(in)    :: P
    type(SAM)        , intent(in)    :: Q
    type(CAM)        , intent(inout) :: S
    integer          , intent(in)    :: NRHS
    real(dp)         , intent(inout) :: B(:)
    type(Error_Flag) , pointer       :: INFO
    print *,' ** FEAssembleElement dummy: ',OPT,iE,size(k),associated(P%sam), &
         &                                  Q%dummy,S%dummy,NRHS,size(B)
    call errorMsg ('FEAssembleElement',INFO)
  end subroutine FEAssembleElement

  subroutine FEEndAssemblyA (Msg,Q,S,T,INFO)
    Type(MessageLevel), intent(in)    :: Msg
    Type(SAM)         , intent(inout) :: Q
    Type(CAM)         , intent(inout) :: S
    Type(GSFColSup)   , intent(inout) :: T
    type(Error_Flag)  , pointer       :: INFO
    print *,' ** FEEndAssemblyA dummy: ',Q%dummy,S%dummy,T%dummy
    call errorMsg ('FEEndAssembly',INFO,Msg)
  end subroutine FEEndAssemblyA

  subroutine FEEndAssemblyB (Q,S)
    Type(SAM)         , intent(inout) :: Q
    Type(CAM)         , intent(inout) :: S
    print *,' ** FEEndAssemblyB dummy: ',Q%dummy,S%dummy
  end subroutine FEEndAssemblyB

  subroutine FEFactorize (Msg,S,T,INFO)
    type(MessageLevel), intent(inout) :: Msg
    type(CAM)         , intent(inout) :: S
    type(GSFColSup)   , intent(inout) :: T
    type(Error_Flag)  , pointer       :: INFO
    print *,' ** FEFactorize dummy: ',S%dummy,T%dummy
    call errorMsg ('FEFactorize',INFO,Msg)
  end subroutine FEFactorize

  subroutine FESolve (Msg,Action,N,NRHS,Q,T,B,X,INFO,KEEPB)
    use kindModule    , only          :  dp
    type(MessageLevel), intent(in)    :: Msg
    character(len=4)  , intent(in)    :: Action
    integer           , intent(in)    :: N, NRHS
    type(SAM)         , intent(in)    :: Q
    type(GSFColSup)   , intent(inout) :: T
    real(dp)          , intent(in)    :: B(:,:)
    real(dp)          , intent(out)   :: X(:,:)
    type(Error_Flag)  , pointer       :: INFO
    logical, optional , intent(in)    :: KEEPB
    X = 0.0_dp
    print *,' ** FESolve dummy: ',Action,N,NRHS,Q%dummy,T%dummy, &
         &                        size(B),size(X),present(KEEPB)
    call errorMsg ('FESolve',INFO,Msg)
  end subroutine FESolve

  subroutine FEExtractS (Q,T,NRHS,B,R,INFO)
    use kindModule  , only        :  dp
    type(SAM)       , intent(in)  :: Q
    type(GSFColSup) , intent(in)  :: T
    integer         , intent(in)  :: NRHS
    real(dp)        , intent(in)  :: B(:,:)
    real(dp)        , intent(out) :: R(:,:)
    type(Error_Flag), pointer     :: INFO
    R = 0.0_dp
    print *,' ** FEExtractS dummy: ',Q%dummy,T%dummy,NRHS,size(B),size(R)
    call errorMsg ('FEExtractS',INFO)
  end subroutine FEExtractS

  subroutine FEFindIsolatedZeroPivots (S,indices,INFO)
    type(CAM)       , intent(inout) :: S
    integer         , pointer       :: indices(:)
    type(Error_Flag), pointer       :: INFO
    nullify(indices)
    print *,' ** FEFindIsolatedZeroPivots dummy: ',S%dummy
    call errorMsg ('FEFindIsolatedZeroPivots',INFO)
  end subroutine FEFindIsolatedZeroPivots

  subroutine FEInsertPivots (S,indices,values,INFO)
    use kindModule  , only          :  dp
    type(CAM)       , intent(inout) :: S
    integer         , intent(in)    :: indices(:)
    real(dp)        , intent(in)    :: values(:)
    type(Error_Flag), pointer       :: INFO
    print *,' ** FEInsertPivots dummy: ',S%dummy,size(indices),size(values)
    call errorMsg ('FEInsertPivots',INFO)
  end subroutine FEInsertPivots

  subroutine FEExtractDiagAA (M,Q,S,A,INFO)
    use kindModule  , only          :  dp
    integer         , intent(in)    :: M
    type(SAM)       , intent(in)    :: Q
    type(CAM)       , intent(inout) :: S
    real(dp)        , intent(out)   :: A(:)
    type(Error_Flag), pointer       :: INFO
    A = 0.0_dp
    print *,' ** FEExtractDiagAA dummy: ',M,Q%dummy,S%dummy
    call errorMsg ('FEExtractDiagA',INFO)
  end subroutine FEExtractDiagAA

  subroutine FEExtractDiagAB (S,A,INFO)
    use kindModule  , only          :  dp
    type(CAM)       , intent(inout) :: S
    real(dp)        , intent(out)   :: A(:)
    type(Error_Flag), pointer       :: INFO
    A = 0.0_dp
    print *,' ** FEExtractDiagAB dummy: ',S%dummy
    call errorMsg ('FEExtractDiagA',INFO)
  end subroutine FEExtractDiagAB

  subroutine FEExtractK12A (S,T,A,IntExt,INFO)
    use kindModule  , only          :  dp
    type(SAM)       , intent(in)    :: S
    type(CAM)       , intent(inout) :: T
    real(dp)        , intent(out)   :: A(:,:)
    integer         , intent(out)   :: IntExt(:)
    type(Error_Flag), pointer       :: INFO
    A = 0.0_dp
    IntExt = 0
    print *,' ** FEExtractK12A dummy: ',S%dummy,T%dummy
    call errorMsg ('FEExtractK12',INFO)
  end subroutine FEExtractK12A

  subroutine FEExtractK12B (S,T,IROW,JCOL,A,IntExt,INFO)
    use kindModule  , only          :  dp
    type(SAM)       , intent(in)    :: S
    type(CAM)       , intent(inout) :: T
    integer         , pointer       :: IROW(:), JCOL(:)
    real(dp)        , pointer       :: A(:)
    integer         , intent(out)   :: IntExt(:)
    type(Error_Flag), pointer       :: INFO
    nullify(IROW,JCOL,A)
    IntExt = 0
    print *,' ** FEExtractK12B dummy: ',S%dummy,T%dummy
    call errorMsg ('FEExtractK12',INFO)
  end subroutine FEExtractK12B

  subroutine FEExtractK22 (LowAll,S,T,A,INFO)
    use kindModule  , only          :  dp
    character(len=1), intent(in)    :: LowAll
    type(SAM)       , intent(in)    :: S
    type(CAM)       , intent(inout) :: T
    real(dp)        , intent(out)   :: A(:,:)
    type(Error_Flag), pointer       :: INFO
    A = 0.0_dp
    print *,' ** FEExtractK22 dummy: ',LowAll,S%dummy,T%dummy
    call errorMsg ('FEExtractK22',INFO)
  end subroutine FEExtractK22

  subroutine FEExtractB12 (Msg,Q,T,B,FstCol,LstCol,IntExt,INFO)
    use kindModule    , only          :  dp
    type(MessageLevel), intent(in)    :: Msg
    type(SAM)         , intent(in)    :: Q
    type(GSFColSup)   , intent(inout) :: T
    real(dp)          , intent(out)   :: B(:,:)
    integer           , intent(in)    :: FstCol, LstCol
    integer           , intent(out)   :: IntExt(:)
    type(Error_Flag)  , pointer       :: INFO
    B = 0.0_dp
    IntExt = 0
    print *,' ** FEExtractB12 dummy: ',Q%dummy,T%dummy,FstCol,LstCol
    call errorMsg ('FEExtractB12',INFO,Msg)
  end Subroutine FEExtractB12

  subroutine FEExtractL22 (LowAll,Q,T,A,INFO)
    use kindModule  , only          :  dp
    character(len=1), intent(in)    :: LowAll
    type(SAM)       , intent(in)    :: Q
    type(GSFColSup) , intent(inout) :: T
    real(dp)        , intent(out)   :: A(:,:)
    type(Error_Flag), pointer       :: INFO
    A = 0.0_dp
    print *,' ** FEExtractL22 dummy: ',LowAll,Q%dummy,T%dummy
    call errorMsg ('FEExtractL22',INFO)
  end subroutine FEExtractL22

  subroutine FEPRM (AllInt,perm,M,N,Q,alpha,A,B,beta,C,INFO)
    use kindModule  , only          :  dp
    character(len=1), intent(in)    :: AllInt, perm
    integer         , intent(in)    :: M, N
    type(SAM)       , intent(in)    :: Q
    real(dp)        , intent(in)    :: alpha
    type(CAM)       , intent(inout) :: A
    real(dp)        , intent(in)    :: B(:,:)
    real(dp)        , intent(in)    :: beta
    real(dp)        , intent(inout) :: C(:,:)
    type(Error_Flag), pointer       :: INFO
    print *,' ** FEPRM dummy: ',AllInt,perm,M,N, &
         &                      Q%dummy,alpha,A%dummy,size(B),beta,size(C)
    call errorMsg ('FEPRM',INFO)
  end subroutine FEPRM

  subroutine FETRA (AllInt,perm,M,N,Q,alpha,A,B,beta,C,INFO)
    use kindModule  , only          :  dp
    character(len=1), intent(in)    :: AllInt, perm
    integer         , intent(in)    :: M, N
    type(SAM)       , intent(in)    :: Q
    real(dp)        , intent(in)    :: alpha
    type(CAM)       , intent(inout) :: A
    real(dp)        , intent(in)    :: B(:,:)
    real(dp)        , intent(in)    :: beta
    real(dp)        , intent(inout) :: C(:,:)
    type(Error_Flag), pointer       :: INFO
    print *,' ** FETRA dummy: ',AllInt,perm,M,N, &
         &                      Q%dummy,alpha,A%dummy,size(B),beta,size(C)
    call errorMsg ('FETRA',INFO)
  end subroutine FETRA

  subroutine FETRMM (AllInt,transa,diag,perm,M,N,Q,alpha,A,B,C,INFO)
    use kindModule  , only          :  dp
    character(len=1), intent(in)    :: AllInt, transa, diag, perm
    integer         , intent(in)    :: M, N
    type(SAM)       , intent(in)    :: Q
    type(GSFColSup) , intent(inout) :: A
    real(dp)        , intent(in)    :: alpha, B(:,:)
    real(dp)        , intent(out)   :: C(:,:)
    type(Error_Flag), pointer       :: INFO
    C = 0.0_dp
    print *,' ** FETRMM dummy: ',AllInt,transa,diag,perm,M,N, &
         &                       Q%dummy,A%dummy,alpha,size(B),size(C)
    call errorMsg ('FETRMM',INFO)
  end subroutine FETRMM

  subroutine FEPermuteVectorsA (AllInt,Ext2Int,M,Q,B,INFO)
    use kindModule  , only          :  dp
    character(len=1), intent(in)    :: AllInt, Ext2Int
    integer         , intent(in)    :: M
    type(SAM)       , intent(in)    :: Q
    real(dp)        , intent(inout) :: B(:)
    type(Error_Flag), pointer       :: INFO
    print *,' ** FEPermuteVectorsA dummy: ',AllInt,Ext2Int,M,Q%dummy,size(B)
    call errorMsg ('FEPermuteVectors',INFO)
  end subroutine FEPermuteVectorsA

  subroutine FEPermuteVectorsB (AllInt,Ext2Int,M,N,Q,B,INFO)
    use kindModule  , only          :  dp
    character(len=1), intent(in)    :: AllInt, Ext2Int
    integer         , intent(in)    :: M, N
    type(SAM)       , intent(in)    :: Q
    real(dp)        , intent(inout) :: B(:,:)
    type(Error_Flag), pointer       :: INFO
    print *,' ** FEPermuteVectorsB dummy: ',AllInt,Ext2Int,M,N,Q%dummy,size(B)
    call errorMsg ('FEPermuteVectors',INFO)
  end subroutine FEPermuteVectorsB

  subroutine FEPermuteVectorsC (AllInt,Ext2Int,M,N,Q,B,C,INFO)
    use kindModule  , only          :  dp
    character(len=1), intent(in)    :: AllInt, Ext2Int
    integer         , intent(in)    :: M, N
    type(SAM)       , intent(in)    :: Q
    real(dp)        , intent(in)    :: B(:,:)
    real(dp)        , intent(out)   :: C(:,:)
    type(Error_Flag), pointer       :: INFO
    C = 0.0_dp
    print *,' ** FEPermuteVectorsC dummy: ',AllInt,Ext2Int,M,N,Q%dummy,size(B)
    call errorMsg ('FEPermuteVectors',INFO)
  end subroutine FEPermuteVectorsC

  subroutine FEPermuteIndices (Ext2Int,Q,indices,INFO)
    character(len=1), intent(in)    :: Ext2Int
    type(SAM)       , intent(in)    :: Q
    integer         , intent(inout) :: indices(:)
    type(Error_Flag), pointer       :: INFO
    print *,' ** FEPermuteIndices dummy: ',Ext2Int,Q%dummy,size(indices)
    call errorMsg ('FEPermuteIndices',INFO)
  end subroutine FEPermuteIndices

  subroutine FEDestroyA (Msg,Q,S,T,INFO)
    type(MessageLevel)       , intent(in) :: Msg
    type(SAM)      , optional, pointer    :: Q
    type(CAM)      , optional, pointer    :: S
    type(GSFColSup), optional, pointer    :: T
    type(Error_Flag)         , pointer    :: INFO
    print *,' ** FEDestroyA dummy: ',present(Q),present(S),present(T)
    call errorMsg ('FEDestroy',INFO,Msg)
  end subroutine FEDestroyA

  subroutine FEDestroyB (Q,S,T,INFO)
    type(SAM)      , optional, pointer :: Q
    type(CAM)      , optional, pointer :: S
    type(GSFColSup), optional, pointer :: T
    type(Error_Flag)         , pointer :: INFO
    print *,' ** FEDestroyB dummy: ',present(Q),present(S),present(T)
    call errorMsg ('FEDestroy',INFO)
  end subroutine FEDestroyB

  subroutine FEErrorHandler (INFO,P,Q,output)
    use FEData, only : FEDataInput
    type(Error_Flag) , pointer              :: INFO
    type(FEDataInput), optional, intent(in) :: P
    type(SAM)        , optional, intent(in) :: Q
    integer          , optional, intent(in) :: output
    if (present(output) .and. associated(INFO)) then
       write(output,*) ' ** FEErrorHandler: ',present(P),present(Q)
    end if
  end subroutine FEErrorHandler

  function FEPivotError (INFO,checkNegative)
    Logical                       :: FEPivotError
    Type(Error_Flag) , pointer    :: INFO
    Logical, optional, intent(in) :: checkNegative
    if (associated(INFO)) then
       FEPivotError = .true.
    else if (present(checkNegative)) then
       FEPivotError = checkNegative
    else
       FEPivotError = .false.
    end if
  end function FEPivotError

  subroutine FESetError (INFO,ProcName,ErrorNumber)
    Type(Error_Flag) , pointer    :: INFO
    Character(len=*) , intent(in) :: ProcName
    Integer, optional, intent(in) :: ErrorNumber
    call errorMsg (ProcName,INFO)
    if (present(ErrorNumber)) print *,'    Error =',ErrorNumber
  end subroutine FESetError

  subroutine GSFDump (Q,S,T,lpu,INFO)
    type(SAM)       , intent(in) :: Q
    type(CAM)       , intent(in) :: S
    type(GSFColSup) , intent(in) :: T
    integer         , intent(in) :: lpu
    type(Error_Flag), pointer    :: INFO
    write(lpu,*) ' ** GSFDump: ',Q%dummy,S%dummy,T%dummy,associated(INFO)
  end subroutine GSFDump

end module FELinearSolver
