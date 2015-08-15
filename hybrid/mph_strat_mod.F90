#include "mph_defs_f.h"
#include "mph_macros_f.h"  
! Handle MaPHyS resolution partitioning strategy
Module mph_strat_mod
  Use MPH_error_mod
  Implicit None

  Public :: mph_strat_init
  Public :: mph_strat_set

  ! Private constants
  Character(len=MAPHYS_STRL), Parameter, Private :: &
       FLNAME= "mph_strat_mod.F90"

  

Contains

  Subroutine mph_strat_init(ikeep,rkeep)
    Use MPH_strat_utils_mod
    Implicit None
    Integer      , Intent(inout) :: ikeep(MAPHYS_IKEEP_SIZE)
    Real(Kind=8) , Intent(inout) :: rkeep(MAPHYS_RKEEP_SIZE)
    Integer :: i

    Do i= 1, MAPHYS_IKEEP_SIZE
       Call unseti(ikeep(i))
    End Do

    Do i= 1, MAPHYS_RKEEP_SIZE
       Call unsetr(rkeep(i))
    End Do

  End Subroutine mph_strat_init

  !> set the maphys strategy 
  !! 
  !! @param [in,out] ikeep  the strategy + state of maphys (integers).
  !!             Assumes on input that several component must be set :
  !!             - IKEEP_CURJOB
  !!             - IKEEP_NBDOMAINS           
  !!             - IKEEP_SYMMETRY if IKEEP_CURJOB > CURJOB_Is_Analysis 
  !!             - IKEEP_MPIRANK, IKEEP_HOSTRANK
  !! 
  !! @param [in,out] rkeep  the strategy + state of maphys (reals).
  !! @param [in    ] icntl  the user controls (integers)
  !! @param [in    ] rcntl  the user controls (reals)
  !! @param [in,out] comm  the user controls (reals)
  Subroutine mph_strat_set(ikeep,rkeep,icntl,rcntl,comm,info)
    Use MPH_maphys_enum
    Use MPH_strat_utils_mod
    Implicit None
    Integer      , Intent(inout) :: ikeep(MAPHYS_IKEEP_SIZE)
    Real(Kind=8) , Intent(inout) :: rkeep(MAPHYS_RKEEP_SIZE)
    Integer      , Intent(in)    :: icntl(MAPHYS_ICNTL_SIZE)
    Real(Kind=8) , Intent(in)    :: rcntl(MAPHYS_RCNTL_SIZE)
    Integer      , Intent(inout) :: comm
    Integer      , Intent(out)   :: info

    info = MPH_SUCCESS
    Call set_print
    Call set_ana
    Call set_fac
    Call set_pcd
    Call set_slv

  Contains
    ! 
    ! Warning : this routine do not follow standard coding style
    ! and use several construction that may handicap your debugger.
    !
    ! This special coding/construction are here to reduce 
    ! maintenance cost/headaches of the strategy selection. 
    ! This is done by loosely coupling each strategy while 
    ! providing enough parameter checking.
    ! 
    ! Each strategy have its own subroutine. General rules in 
    ! a subroutine is as follows:
    !
    !  1. What should I set ? 
    !   You should set yourself and your children, but :
    !   - if a previous error occured -> Return (do not log)
    !   - if I was already set. Set my children 
    !     (because they may not be setted)
    !
    !  2. How to set myself ?
    !   - If decision is selfmade, simply check the ?cntl value
    !     and perform copy from ?cntl to ?keep.
    !     The "generic setters" may help you to perform this operation.
    !
    !   - If decision depends on other strategy(ies).
    !     simply call the setters of the other strategy(ies),
    !     then take the decision. Try to check the value if possible.
    !
    !     - If you cannot decide, leave stuff as it was,
    !       so that the caller may use "is_unset()".
    !
    !     - If a dependency to other strategy is optional, 
    !       call its setter and use "is_unset()". 
    !
    !     - If a dependency is compulsory, 
    !       call its setter. If its setter failed, return.
    !       If it passed, try to check the return value.
    !
    !  3. How to set my children ?
    !     Simply call their setters.
    !
    !                        Written by Yohan Lee-tin-yien (2011-11-26)
    !

    !
    !*** generic setters ***!
    !
    Subroutine set_enum(k,c,min,max)
      Integer, Intent(in) :: c,k,min,max
      !* copy an enumeration after checking its range *
      ! avoid on failure or if already set
      If(info < 0) Return
      If(is_seti(ikeep(k)))Return
      ! check icntl
      If(.Not.checkRange(icntl,c,min,max))Then
         info=-__LINE__
         Return
      EndIf;
      ! copy
      ikeep (k) = icntl(c)
    End Subroutine set_enum

    Subroutine set_positiveI(k,c)
      Integer, Intent(in) :: c,k
      ! avoid on failure or if already set
      If(info < 0) Return
      If(is_seti(ikeep(k)))Return
      ! check icntl
      If(.Not.checkPositiveI(icntl,c))Then
         info=-__LINE__
         Return
      EndIf;
      ! copy
      ikeep (k) = icntl(c)
    End Subroutine set_positiveI

    Subroutine set_positiveR(k,c)
      Integer, Intent(in) :: c,k
      ! avoid on failure or if already set
      If(info < 0) Return
      If(is_setr(rkeep(k)))Return
      ! check icntl
      If(.Not.checkPositiveR(rcntl,c))Then
         info=-__LINE__
         Return
      EndIf;
      ! copy
      rkeep (k) = rcntl(c)
    End Subroutine set_positiveR
    !
    !*** Actual implementations by order of appearance ***!
    !
    Subroutine set_print
      ! set me
      Call set_enum(IKEEP_PRINT_Cntls,ICNTL_PRINT_Cntls,0,2)
      Call set_enum(IKEEP_PRINT_Infos,ICNTL_PRINT_Infos,0,2)
    End Subroutine set_print
    !****
    Subroutine set_ana
      Call set_insystem
    End Subroutine set_ana
    !****
    Subroutine set_insystem
      ! set me
      Call set_enum(IKEEP_INSYSTEM,ICNTL_INSYSTEM,0,2)
      MPH_ONFAILURE_RETURN(info)
      ! set my children
      If (ikeep(IKEEP_INSYSTEM ) == INSYSTEM_IsCentralized ) Call set_part
    End Subroutine set_insystem
    !****
    Subroutine set_part
      ! set me 
      Call set_enum(IKEEP_PART_Strategy,ICNTL_PART_Strategy,0,4)
      MPH_ONFAILURE_RETURN(info)
      ! specific checks
      If(.Not. is_aPowerOf2(ikeep(IKEEP_NBDOMAINS)))Then
         Call MPH_LogWithInfo(MSG_ERROR,ikeep(IKEEP_NBDOMAINS),&
              "MaPHyS partitioner only support *a power of 2* subdomains. nb subdomains=")
         CHCKASSRT(.False.,info)
      End If
    End Subroutine set_part
    !****
    Subroutine set_fac
      Call set_schur
    End Subroutine set_fac
    !****
    Subroutine set_schur
      ! set me 
      Call set_enum(IKEEP_SCHUR_Strategy,ICNTL_SCHUR_Strategy,0,2)
      ! set my children
      Select Case(ikeep(IKEEP_SCHUR_Strategy))
      Case(SCHUR_STRATEGY_ApprxWithIluT)
         Call set_facsds
         Call set_pilut
      Case(SCHUR_STRATEGY_EXACT)
         Call set_facsds
      Case(SCHUR_STRATEGY_ExactFromFacto)
         Call MPH_Log(MSG_ERROR,"strategy not implemented yet.")
         CHCKASSRT(.False.,info)
      Case Default
         Call MPH_LogWithInfo(MSG_ERROR,__LINE__,&
              "Unsupported case in file '"//Trim(flname)//"' at line: ") 
         ASSRT(.False.)
      End Select
    End Subroutine set_schur
    !****
    Subroutine set_facsds
      Call set_enum(IKEEP_SDS_Facto,ICNTL_SDS_Default,1,3)
      If (ikeep(IKEEP_SDS_Facto) == SDS_DEFAULT_UsrDefined)&
           Call set_enum(IKEEP_SDS_Facto,ICNTL_SDS_Facto,1,2)
    End Subroutine set_facsds
    !**** 
    Subroutine set_pilut
      ! set me
      Call set_positiveI(IKEEP_ILU_LUFill         ,ICNTL_ILU_LUFill)
      Call set_positiveI(IKEEP_ILU_SCHURFill      ,ICNTL_ILU_SCHURFill)
      Call set_positiveR(RKEEP_ILU_LUThreshold    ,RCNTL_ILU_LUThreshold)
      Call set_positiveR(RKEEP_ILU_SchurThreshold ,RCNTL_ILU_SchurThreshold)
    End Subroutine set_pilut
    !****
    Subroutine set_pcd
      Integer, Parameter :: kp = IKEEP_PCD_Strategy
      Logical :: tmp
      ! set 
      Call set_enum(kp,ICNTL_PCD_Strategy,-1,5)
      ! autodetect
      If(ikeep(kp) == PCD_STRATEGY_isAutodetected )Then
         Call set_schur ! dependency
         Select Case(ikeep(IKEEP_SCHUR_Strategy))
         Case( SCHUR_STRATEGY_Exact)
            ikeep(kp) = PCD_STRATEGY_isLocalExact
         Case( SCHUR_STRATEGY_ExactFromFacto );
            ikeep(kp) = PCD_STRATEGY_isLocalExact
         Case( SCHUR_STRATEGY_ApprxWithIluT )
            ikeep(kp) = PCD_STRATEGY_isForcedByILUT
         Case Default   
            ikeep(kp) = PCD_STRATEGY_isNone
         End Select
      End If
      ! check
      Call set_schur ! dependency
      Select Case(ikeep(kp))
      Case (PCD_STRATEGY_isLocalApprox,PCD_STRATEGY_isLocalExact)
         tmp = checkValue(icntl,ICNTL_SCHUR_Strategy,SCHUR_STRATEGY_Exact) 
         If(.Not. tmp) Call printInconsistencyI(&
              icntl, ICNTL_PCD_Strategy, ICNTL_SCHUR_Strategy )
         ASSRT(tmp)
      Case (PCD_STRATEGY_isForcedByILUT)
         tmp = checkValue(icntl,ICNTL_SCHUR_Strategy,SCHUR_STRATEGY_ApprxWithIluT) 
         If(.Not. tmp) Call printInconsistencyI(&
              icntl, ICNTL_PCD_Strategy, ICNTL_SCHUR_Strategy )
         ASSRT(tmp)
      Case (PCD_STRATEGY_isNone)
         Continue
      Case Default ;
         Call MPH_LogWithInfo(MSG_ERROR,__LINE__,&
              "Unsupported case in file '"//Trim(flname)//"' at line: ") 
         ASSRT(.False.)
      End Select
      ! children
      Select Case(ikeep(kp))
      Case (PCD_STRATEGY_isLocalApprox); call set_pcdlocalapprox
      Case (PCD_STRATEGY_isForcedByILUT); call set_pcdilut
      Case (PCD_STRATEGY_isLocalExact);  Call set_pcdlocalexact
      Case (PCD_STRATEGY_isNone); Continue
      Case Default ;
         Call MPH_LogWithInfo(MSG_ERROR,__LINE__,&
              "Unsupported case in file '"//Trim(flname)//"' at line: ") 
         ASSRT(.False.)
      End Select
    End Subroutine set_pcd
    !****
    Subroutine set_pcdlocalapprox
      Call set_pcdsds
      Call set_positiveR(&
           RKEEP_PCD_SparsifyThreshold,RCNTL_PCD_SparsifyThreshold)
    End Subroutine set_pcdlocalapprox
    !****
    Subroutine set_pcdilut
      Logical :: tmp
      tmp = checkValue(icntl,&
           ICNTL_SCHUR_Strategy,SCHUR_STRATEGY_ApprxWithIluT)
      If(.Not. tmp) Call MPH_LogWithInfo(MSG_ERROR,ICNTL_PCD_Strategy, &
           "This error occured while processing ICNTL:")
      CHCKASSRT( tmp , info)
      Call set_pcdsds
    End Subroutine set_pcdilut
    !****
    Subroutine set_pcdsds
      Call set_enum(IKEEP_SDS_Precond,ICNTL_SDS_Default,1,3)
      If (ikeep(IKEEP_SDS_Precond) == SDS_DEFAULT_UsrDefined)&
           Call set_enum(IKEEP_SDS_Precond,ICNTL_SDS_Precond,1,2)
    End Subroutine set_pcdsds
    !****
    Subroutine set_pcdlocalexact
      Call set_pcddds
      Call set_schurUselessAfter
    End Subroutine set_pcdlocalexact
    !****
    Subroutine set_pcddds
      ! TODO add controls on 
      ! statistic activation (rcond,norm1) here
    End Subroutine set_pcddds
    !****
    Subroutine set_schurUselessAfter

      ikeep(IKEEP_SCHUR_UselessAfter) = CURJOB_isSolve

      Call set_itsmatvect
      If (is_unseti(ikeep(IKEEP_ITS_MatVect))) Return
      If(ikeep(IKEEP_ITS_MatVect) == MAT_VECT_isImplicit)&
           ikeep(IKEEP_SCHUR_UselessAfter) = CURJOB_isPrecond

    End Subroutine set_schurUselessAfter
    !**** 
    Subroutine set_slv
      Call set_itsmatvect
      Call set_itssolver
    End Subroutine set_slv
    !****
    Subroutine set_itsmatvect
      Call set_enum(IKEEP_ITS_MatVect,ICNTL_ITS_MatVect,0,2)
      If (ikeep(IKEEP_ITS_MatVect) == MAT_VECT_isAutodetected  )&
           ikeep(IKEEP_ITS_MatVect) = MAT_VECT_isExplicit
    End Subroutine set_itsmatvect
    !****
    Subroutine set_itssolver
      Use MPH_sparse_matrix_enum
      Integer, Parameter :: kp = IKEEP_ITS_Solver
      ! set me
      Call set_enum(kp,ICNTL_ITS_Solver,1,4)
      ! depend on how the schur was computed
      Call set_pcd ! dependency
      If (is_unsetI(ikeep(IKEEP_PCD_Strategy))) Goto 9
      If ( ikeep(IKEEP_PCD_Strategy) == &
           PCD_STRATEGY_isForcedByILUT)Then
         If (ikeep(kp) == ITS_Solver_isPackCG) & 
              Call mph_log(MSG_WARNING, "Using GMRES instead of CG: &
              & CG does not support preconditioner from PILUT")  
         ikeep(kp) = ITS_Solver_isPackGMRES
      End If
      ! autodetect
      If (is_unsetI(ikeep(IKEEP_SYMMETRY))) Goto 9
      If (ikeep(kp) == ITS_Solver_isAuto  ) Then
         Select Case (ikeep(IKEEP_SYMMETRY)) 
         Case (SM_SYM_isGeneral, SM_SYM_isSymmetric)
            ikeep(kp) = ITS_Solver_isPackGMRES
         Case(SM_SYM_isSPD)
            ikeep(kp) = ITS_Solver_isPackCG
         End Select
      End If
      ! children
9     Continue
      Call set_initGuess
      Call set_itsTolerance
      Call set_itsmaxiter
      Select Case(ikeep(kp))
      Case( ITS_Solver_isPackCG    ); Call  set_packcg
      Case( ITS_Solver_isPackGMRES ); Call  set_packgmres
      Case( ITS_Solver_isPackFGMRES ); Call set_packfgmres      
      Case( ITS_Solver_isAuto      ); Continue
      Case Default; 
         Call MPH_LogWithInfo(MSG_ERROR,__LINE__,&
              "Unsupported case in file '"//Trim(flname)//"' at line: ") 
         ASSRT(.False.) 
      End Select
    End Subroutine set_itssolver
    !****
    Subroutine set_initGuess
      ! only set it during the solve
      If(ikeep(IKEEP_CURJOB) /= CURJOB_IsSolve) Return
      Call set_enum(IKEEP_ITS_InitGuess,ICNTL_ITS_InitGuess,0,1)
    End Subroutine set_initGuess
    !****
    Subroutine set_itsTolerance
      ! only set it during the solve
      If(ikeep(IKEEP_CURJOB) /= CURJOB_IsSolve) Return
      Call set_positiveR(RKEEP_ITS_Tolerance,RCNTL_ITS_Tolerance)
    End Subroutine set_itsTolerance
    !****
    Subroutine set_itsmaxiter
      Integer, Parameter :: kp = IKEEP_ITS_MaxIters
      ! only set it during the solve
      If(ikeep(IKEEP_CURJOB) /= CURJOB_IsSolve) Return
      Call set_positivei(kp,ICNTL_ITS_MaxIters)
      MPH_ONFAILURE_RETURN(info)
      ! autoset
      If (ikeep(kp) == 0)Then
         ikeep(kp) = 100
         ! log
         If (ikeep(IKEEP_MPIrank) /= ikeep(IKEEP_HOSTRANK)) Return
         Call MPH_logWithInfo(MSG_WARNING,ICNTL_ITS_MaxIters,&
              "Iterative solver -> Maximal Nb of Iterations unset. ICNTL:")
         Call MPH_LogWithInfo(MSG_WARNING, ikeep(kp),&
              "Iterative solver -> Maximal Nb of Iterations unset. Take :")
      End If
    End Subroutine set_itsmaxiter
    !****
    Subroutine set_packcg
      ! nothing specific here
    End Subroutine set_packcg
    !****
    Subroutine set_packgmres
      Call set_enum(IKEEP_ITS_OrtStrategy,ICNTL_ITS_OrtStrategy,0,3)
      Call set_enum(IKEEP_ITS_ResStrategy,ICNTL_ITS_ResStrategy,0,1)
      Call set_itsRestart
    End Subroutine set_packgmres
    !****
    Subroutine set_packfgmres
      Call set_enum(IKEEP_ITS_OrtStrategy,ICNTL_ITS_OrtStrategy,0,3)
      Call set_enum(IKEEP_ITS_ResStrategy,ICNTL_ITS_ResStrategy,0,1)
      Call set_itsRestart
    End Subroutine set_packfgmres
    
    Subroutine set_itsrestart
      Integer, Parameter :: kp = IKEEP_ITS_Restart
      ! only set it during the solve
      If(ikeep(IKEEP_CURJOB) /= CURJOB_IsSolve) Return
      Call set_positivei(kp,ICNTL_ITS_Restart)
      MPH_ONFAILURE_RETURN(info)
      ! autoset
      If (ikeep(kp) == 0)Then
         ikeep(kp) = 100
         ! log
         If (ikeep(IKEEP_MPIrank) /= ikeep(IKEEP_HOSTRANK)) Return
         Call MPH_logWithInfo(MSG_WARNING,ICNTL_ITS_Restart,&
              "Iterative solver -> Restart every X iterations unset. ICNTL:")
         Call MPH_LogWithInfo(MSG_WARNING, ikeep(kp),&
              "Iterative solver -> Restart every X iterations unset. Take :")
      End If
    End Subroutine set_itsrestart

End Subroutine mph_strat_set



End Module mph_strat_mod

