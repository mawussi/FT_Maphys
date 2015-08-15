! Warning: XMPH_GENFILE_COMMENT
#include "mph_defs_f.h"
#include "mph_macros_f.h"

! [+] module : XMPH_state_print_mod --------------------------------------------------
!
!> module to handle statistics in maphys
!!
Module XMPH_state_print_mod
  !* Modules *!
  Use MPH_log_mod
  Use MPH_error_mod
  Use MPH_maphys_enum
  Use XMPH_maphys_type
  
  !* No Implicit Typing *!
  Implicit None

  !* Access specifiers *!

  Public :: XMPH_state_print_stepstart
  Public :: XMPH_state_print_stepend
  Private :: XMPH_state_print_step
  Public :: XMPH_state_print_listCntl
  Public :: XMPH_state_print_listInfo
  Public :: XMPH_state_print

  Private :: PrintSH  ! section header
  Private :: PrintTH1 ! table header (1 value)
  Private :: PrintTH4 ! table header (4 values)
  Private :: PrintTF1 ! table footer (1 value)
  Private :: PrintTF4 ! table footer (4 values)
  Private :: PrintTL  ! table line 
  Private :: PrintINTERP  ! table line   
  Public :: XMPH_state_print_legend
  Public :: XMPH_state_print_about

  !* Private constants *!

  Integer, Private , Parameter :: F_ICNTL   =0 
  Integer, Private , Parameter :: F_RCNTL   =1
  Integer, Private , Parameter :: F_IINFOG  =2
  Integer, Private , Parameter :: F_RINFOG  =3
  Integer, Private , Parameter :: F_IINFO1  =4
  Integer, Private , Parameter :: F_RINFO1  =5
  Integer, Private , Parameter :: F_IINFO   =6
  Integer, Private , Parameter :: F_RINFO   =7
  Integer, Private , Parameter :: F_IsSYM   =8 ! print symmetry
  Integer, Private , Parameter :: F_IsARITH =9 ! print arithmetic 
  Integer, Private , Parameter :: F_IsVERSION =10 ! print the version

  Character(len=MAPHYS_STRL), Private, Parameter ::&
       FLNAME = "XMPH_state_print_mod.F90"


  !* routines *!

Contains

  Subroutine XMPH_state_print_stepstart(mphs)
    Use XMPH_maphys_type
    Implicit None
    Type(XMPH_maphys_t), Intent(in) :: mphs
    Call XMPH_state_print_step(mphs,1)
  End Subroutine XMPH_state_print_stepstart

  Subroutine XMPH_state_print_stepend  (mphs)
    Use XMPH_maphys_type
    Implicit None
    Type(XMPH_maphys_t), Intent(in) :: mphs
    Call XMPH_state_print_step(mphs,2)
  End Subroutine XMPH_state_print_stepend

  Subroutine XMPH_state_print_step(mphs,startOrEnd)
    Implicit None

    !* Arguments *!
    Type(XMPH_maphys_t), Intent(in) :: mphs
    Integer            , Intent(in) :: startOrEnd

    !* Local variables *!
    Integer :: verb
    Integer :: unit
    Integer :: i
    Character(len=MAPHYS_STRL) :: tmp, tmp2
    
    !- End of header -----------------------------------------------------------

    ! Init
    If (mphs%ikeep(IKEEP_MPIRANK) /= mphs%ikeep(IKEEP_HOSTRANK) ) Return 
    !
    Call mph_log_GetVerbosity(verb)
    If ( verb < MSG_VERBOSE ) Return

    Call mph_log_GetUnit(MSG_VERBOSE,unit)
    If ( unit < 0) Return

    ! select comment
    Select Case(mphs%ikeep(IKEEP_CURJOB))
    Case(CURJOB_IsInit    ); tmp = "(Initialization)"
    Case(CURJOB_IsAnalysis); tmp = "(Analysis)"
    Case(CURJOB_IsFacto)   ; tmp = "(Factorisation)" 
    Case(CURJOB_IsPrecond) ; tmp = "(Preconditioning)" 
    Case(CURJOB_IsSolve)   ; tmp = "(Solving)" 
    Case(CURJOB_IsFinalize); tmp = "(Finalizing)" 
    Case Default           ; tmp = "(Unknown)" 
    End Select



    ! print the message
    Select Case (startOrEnd)
    Case(1)

       ! Avoid step that do not print something
       Select Case(mphs%ikeep(IKEEP_CURJOB))
       Case(CURJOB_IsInit    ); Return
       Case(CURJOB_IsAnalysis); Continue
       Case(CURJOB_IsFacto)   ; Return
       Case(CURJOB_IsPrecond) ; Return
       Case(CURJOB_IsSolve)   ; Continue
       Case(CURJOB_IsFinalize); Return
       Case Default           ; Continue
       End Select

       Write(tmp2,'(A,I1,1X,A)') "* Starting Job",mphs%ikeep(IKEEP_CURJOB),Trim(tmp)
    Case(2)
       Write(tmp2,'(A,I1,1X,A)') "* Ending Job",mphs%ikeep(IKEEP_CURJOB),Trim(tmp)
    End Select

    Write(unit,*)
    Write(unit,*)
    Write(unit,*) Trim(tmp2)
    Write(unit,*) ('=',i=1,Len_Trim(tmp2))
    Write(unit,*)

  End Subroutine XMPH_state_print_step


  !
  !-----------------------------------------------------------------------------
  !
  !> List the maphys controls
  !!
  !! @param[in ] mphs  The maphys instance containing the stats.
  !!          - comm
  !!          - {icntl,rcntl}
  !! @param[out] info  the routine status
  !!
  !! @author Yohan Lee-tin-yien
  !!
  Subroutine XMPH_state_print_listCntl( mphs, info )
    Implicit None

    !* Arguments *!

    Type(XMPH_maphys_t), Intent(in ) :: mphs
    Integer       , Intent(out) :: info

    !* Local variables *!

    ! scalars
    Integer :: unit
    Integer :: i

    ! string
    Character(len=MAPHYS_STRL) :: iform ! format for the integers
    Character(len=MAPHYS_STRL) :: rform ! format for the reals
    Character(len=MAPHYS_STRL) :: lform ! format for the legends

    !- End of header -----------------------------------------------------------

    ! Exit early
    info = 0

    If (mphs%ikeep(IKEEP_MPIRANK) /= mphs%ikeep(IKEEP_HOSTRANK) ) Return 
    !
    Select Case(mphs%ikeep(IKEEP_PRINT_Cntls))
    Case( PRINT_Not )
       Return
    Case( PRINT_Once )
       If (mphs%ikeep(IKEEP_Curjob) /= CURJOB_IsAnalysis ) Return 
    Case( PRINT_EveryStep )
       Continue
    Case Default
       CHCKASSRT(.False., info )
       Return
    End Select
    !
    Call mph_log_GetUnit(MSG_STD,unit)
    If ( unit < 0) info = 1
    If ( unit < 0) Return

    !

    Call PrintSH(unit,"List of controls") 

    lform='(A,A5,1A11)'
    iform='(A,I5,1I11)'
    rform='(A,I5,1PE11.3)'

    Write(unit,FMT=lform) "FIELD ", "INDEX", &
         "VALUE"
    Write(unit,*) "--"
    Do i=1,MAPHYS_ICNTL_SIZE
       Write(unit,FMT=iform) "ICNTL ", i, mphs%icntl(i)
    End Do
    Write(unit,*) "--"
    Do i=1,MAPHYS_RCNTL_SIZE
       Write(unit,FMT=rform) "RCNTL ",i, mphs%rcntl(i)
    End Do

        Call PrintSH(unit,"List of controls") 

    lform='(A,A5,1A11)'
    iform='(A,I5,1I11)'
    rform='(A,I5,1PE11.3)'

  End Subroutine XMPH_state_print_listCntl

  !
  !-----------------------------------------------------------------------------
  !
  !> List the maphys infos
  !!
  !! @param[in ] mphs  The maphys instance containing the stats.
  !!          - {i,r}info*
  !! @param[out] info  the routine status
  !!
  !! @author Yohan Lee-tin-yien
  !!
  Subroutine XMPH_state_print_listInfo( mphs, info )
    Implicit None

    !* Arguments *!

    Type(XMPH_maphys_t), Intent(in ) :: mphs
    Integer       , Intent(out) :: info

    !* Local variables *!

    ! scalars
    Integer :: unit
    Integer :: i

    ! string
    Character(len=MAPHYS_STRL) :: iform ! format for the integers
    Character(len=MAPHYS_STRL) :: rform ! format for the reals
    Character(len=MAPHYS_STRL) :: lform ! format for the legends

    !- End of header -----------------------------------------------------------

    ! Exit early
    info = 0
    If (mphs%ikeep(IKEEP_MPIRANK) /= mphs%ikeep(IKEEP_HOSTRANK) ) Return 

    Select Case(mphs%ikeep(IKEEP_PRINT_Infos))
    Case( PRINT_Not )
       Return
    Case( PRINT_Once )
       If (mphs%ikeep(IKEEP_Curjob) /= CURJOB_IsAnalysis ) Return 
    Case( PRINT_EveryStep )
       Continue
    Case Default
       CHCKASSRT(.False., info )
       Return
    End Select

    !
    Call mph_log_GetUnit(MSG_STD,unit)
    If ( unit < 0) info = 1
    If ( unit < 0) Return

    Call PrintSH(unit,"List of infos") 

    !-- iinfog / rinfog

    lform='(A,A5,1A11)'
    iform='(A,I5,1I11)'
    rform='(A,I5,1PE11.3)'
    
    Write(unit,*)
    Write(unit,FMT=lform) "FIELD ", "INDEX", &
         "VALUE"
    Write(unit,*) "--"
    Do i=1,MAPHYS_IINFOG_SIZE
       Write(unit,FMT=iform) "IINFOG", i, mphs%iinfog(i)
    End Do
    Write(unit,*) "--"
    Do i=1,MAPHYS_RINFOG_SIZE
       Write(unit,FMT=rform) "RINFOG",i, mphs%rinfog(i)
    End Do

#if MAPHYS_DEBUG
    Write(unit,*)
    Write(unit,FMT=lform) "FIELD ", "INDEX", &
         "VALUE"
    Write(unit,*) "--"
    Do i=1,MAPHYS_IKEEP_SIZE
       Write(unit,FMT=iform) "IKEEP", i, mphs%ikeep(i)
    End Do
    Write(unit,*) "--"
    Do i=1,MAPHYS_RKEEP_SIZE
       Write(unit,FMT=rform) "RKEEP",i, mphs%rkeep(i)
    End Do
#endif

    !-- iinfo / rinfo

    lform='(A,A5,4A11)'
    iform='(A,I5,2I11,2(1PE11.3))'
    rform='(A,I5,4(1PE11.3))'

    Write(unit,*)
    Write(unit,FMT=lform) "FIELD ", "INDEX", &
         "MIN","MAX", &
         "AVERAGE", "STD.DEV."
    Write(unit,*) "--"
    Do i=1,MAPHYS_IINFO_SIZE
       Write(unit,FMT=iform) "IINFO ", i, &
            mphs%iinfomin(i), mphs%iinfomax(i), &
            mphs%iinfoavg(i), mphs%iinfosig(i)
    End Do
    Write(unit,*) "--"
    Do i=1,MAPHYS_RINFO_SIZE
       Write(unit,FMT=rform) "RINFO ",i, &
            mphs%rinfomin(i), mphs%rinfomax(i), &
            mphs%rinfoavg(i), mphs%rinfosig(i)
    End Do

  End Subroutine XMPH_state_print_listInfo


  ! [+] routine : XMPH_state_print ---------------------------------------------
  !
  !> Print the maphys statistics
  !!
  !! @param[in ] unit  The unit where to print the statistics
  !! @param[in ] mphs  The maphys instance containing the stats.
  !!          - comm
  !!          - {iinfo,rinfo}{,min,max,avg}
  !! @param[out] info  the routine status
  !!
  !! @author Yohan Lee-tin-yien
  !!
  Subroutine XMPH_state_print( mphs, info )
    Use MPH_log_mod
    Use XMPH_maphys_type
    Use MPH_maphys_enum
    Implicit None

    !* Arguments *!

    Type(XMPH_maphys_t) , Intent(in ) :: mphs
    Integer             , Intent(out) :: info

    !* Local variables *!

    ! scalars
    Integer :: verb
    Integer :: unit
    Integer :: lend            ! lenght of the description
    Logical :: isAfterIni
    Logical :: isAfterAna
    Logical :: isAfterFac
    Logical :: isAfterPcd
    Logical :: isAfterSlv
    Logical :: isVerb
    Logical :: isVerb2

    !- End of header -----------------------------------------------------------

    !---------------------------------------------------------------------------
    ! [0] Init
    !---------------------------------------------------------------------------

    info = 0
    If (mphs%ikeep(IKEEP_MPIRANK) /= mphs%ikeep(IKEEP_HOSTRANK) ) Return 

    Call mph_log_GetVerbosity(verb)
    If (verb < MSG_STD ) Return
    isVerb  = (verb >= MSG_VERBOSE  )   
    isVerb2 = (verb >= MSG_VERBOSE2 )   

    Call mph_log_GetUnit(MSG_STD,unit)
    If ( unit < 0) info = 1
    If ( unit < 0) Return

    ! Only print in at during solve step in standard verbosity
    If ( (verb == MSG_STD) .And. &
         ( mphs%ikeep(IKEEP_CURJOB) < CURJOB_IsSolve )) Return

    ! Do not stage the output on verbose 2
    isAfterIni = (mphs%ikeep(IKEEP_CURJOB) >= CURJOB_IsInit    ).Or.(isVerb2)
    isAfterAna = (mphs%ikeep(IKEEP_CURJOB) >= CURJOB_IsAnalysis).Or.(isVerb2)
    isAfterFac = (mphs%ikeep(IKEEP_CURJOB) >= CURJOB_IsFacto   ).Or.(isVerb2)
    isAfterPcd = (mphs%ikeep(IKEEP_CURJOB) >= CURJOB_IsPrecond ).Or.(isVerb2)
    isAfterSlv = (mphs%ikeep(IKEEP_CURJOB) >= CURJOB_IsSolve   ).Or.(isVerb2)

    
#define PARITH(x,y) 
#define PICNTL(x,y) Call PrintTL(F_ICNTL,x,lend,mphs,unit,y)
#define PRCNTL(x,y) Call PrintTL(F_RCNTL,x,lend,mphs,unit,y)
#define PIINFOG(x,y) Call PrintTL(F_IINFOG,x,lend,mphs,unit,y)
#define PRINFOG(x,y) Call PrintTL(F_RINFOG,x,lend,mphs,unit,y)
#define PIINFO1(x,y) Call PrintTL(F_IINFO1,x,lend,mphs,unit,y)
#define PIINFO(x,y) Call PrintTL(F_IINFO,x,lend,mphs,unit,y)
#define PRINFO(x,y) Call PrintTL(F_RINFO,x,lend,mphs,unit,y)
#define PRINTERP(x,y) Call PrintInterp(x,lend,mphs,unit,y)


    !---------------------------------------------------------------------------
    ! [1] Print the statistics about the solution
    !---------------------------------------------------------------------------

    !--
    If( isAfterPcd )Then
       Call PrintSH(unit,'Solution Quality')
    End If

    If( isAfterSlv )Then
       lend=1+Len_Trim("Schur. system backward error")
       Call PrintTH1(unit,lend)
       PRINFOG(RINFOG_BckwrdErrorDist,"|A.x-b|/|b|")
       PRCNTL(RCNTL_ITS_Tolerance,"Wanted precision")             
       PIINFOG(IINFOG_ITS_NbIters,"Nb. of iters.")
       PICNTL(ICNTL_ITS_MaxIters,"Max Nb. of iters.")
       PRINFOG(RINFOG_ITS_BckErr,"Schur. system backward error")
       Call PrintTF1(unit,lend)
    End If

    If( isAfterPcd )Then
       lend=1+Len_Trim("Schur -> rcond1 (estimation)")
       Call PrintTH4(unit,lend)
       PRINFO(RINFO_DSCHUR_NORM1,"Schur -> norm1")
       PRINFO(RINFO_DSCHUR_RCOND1,"Schur -> rcond1 (estimation)")
       Call PrintTF4(unit,lend)
    End If

    !---------------------------------------------------------------------------
    ! [2] Print the strategy used
    !---------------------------------------------------------------------------

    If (isVerb)Then
       Call PrintSH(unit,'Resolution Strategy')
       
       lend = 1+ Len_Trim("Precond Sparse Direct Solver")
       Call PrintTH1(unit,lend)
       Call PrintTL(F_IsARITH,0,lend,mphs,unit," ")
       Call PrintTL(F_IsSYM  ,0,lend,mphs,unit," ")
       If( isAfterAna ) PIINFO1(IINFO_STRAT_PART  , "Partitioning                 ")
       If( isAfterFac ) PIINFO1(IINFO_STRAT_SCHUR , "Schur Computation            ")
       If( isAfterFac ) PIINFO1(IINFO_STRAT_FACTO , "Facto Sparse Direct Solver   ")
       If( isAfterPcd ) PIINFO1(IINFO_STRAT_PCD   , "Precond Strategy             ")
       If( isAfterPcd ) PIINFO1(IINFO_STRAT_PCDSDS, "Precond Sparse Direct Solver ")
       If( isAfterSlv ) PIINFO1(IINFO_STRAT_ITS   , "Iterative Solver             ")
       If( isAfterSlv ) PIINFO1(IINFO_STRAT_ITSSMV, "Iterative Solver - MatVect   ")
       Call PrintTF4(unit,lend)
    End If

    !---------------------------------------------------------------------------
    ! [3] Print the differents sizes in maphys
    !---------------------------------------------------------------------------

    Call PrintSH(unit,'Sizes')

    ! [3.1] Global datas

    lend=Len_Trim("input matrix -> nb. entries") + 1
    Call PrintTH1(unit,lend)
    PIINFOG(IINFOG_MAT_ORDER    ,"input matrix -> order")
    PIINFOG(IINFOG_MAT_NBENTRIES,"input matrix -> nb. entries")
    Call PrintTF1(unit,lend)

    ! [3.1] Local datas
    lend=Len_Trim("prcnd.     -> % of kept entries") + 1
    Call PrintTH4(unit,lend)
    If (isAfterAna)Then
       PIINFO(IINFO_LCMAT_ORDER     ,"lc. matrix -> order")
       PIINFO(IINFO_LCMAT_NBENTRIES ,"lc. matrix -> nb. entries")
       PIINFO(IINFO_LCMAT_SIZEOF    ,"lc. matrix -> mem. cost")
       PIINFO(IINFO_SCHUR_ORDER     ,"schur.     -> order")
    End If
    If (isAfterFac)Then
       PIINFO(IINFO_LCFACTO_SIZEOF  ,"lc. matrix facto. -> mem. cost")
       PIINFO(IINFO_SCHUR_SIZEOF    ,"schur.     -> mem. cost")
       PIINFO(IINFO_SCHUR_NBENTRIES ,"schur.     -> nb. entries")
    End If
    If (isAfterPcd)Then
       PIINFO(IINFO_PCD_ORDER       ,"prcnd.     -> order")
       PIINFO(IINFO_PCD_PERKEPT     ,"prcnd.     -> % of kept entries")
       PIINFO(IINFO_PCD_PERKEPT     ,"prcnd.     -> % of kept entries")
    End If
    Call PrintTF4(unit,lend)

    !---------------------------------------------------------------------------
    ! [4] Print the timings
    !---------------------------------------------------------------------------

    lend=1 + Len_Trim("Timing -> Analyze -> Pretreat Input Mtx.") 

    Call PrintSH(unit,'Timings')

    Call PrintTH4(unit,lend)
    If(isAfterSlv)Then
       PRINFO(RINFO_TIMING_Total,"Timing -> Total")
       PRINFO(RINFO_TIMING_AllSteps,"Timing -> Spent in all steps")
    End If

    If(isAfterAna) PRINFO(RINFO_TIMING_Analysis,"Timing -> Analyze")
    If(isAfterFac) PRINFO(RINFO_TIMING_Facto,"Timing -> Factorise")
    If(isAfterPcd) PRINFO(RINFO_TIMING_Precond,"Timing -> Precond")
    If(isAfterSlv) PRINFO(RINFO_TIMING_Solve,"Timing -> Solve")

    If((isAfterAna).And.(isVerb)) Then
       PRINFO(RINFO_TIMING_PrcInMat,"Timing -> Analyze -> Pretreat Input Mtx.")
       PRINFO(RINFO_TIMING_IN2LCSYS,"Timing -> Analyze -> Form Lc. Systems")
       PRINFO(RINFO_TIMING_EXTBLOCS,"Timing -> Analyze -> Form Subblocs")
    End If

    If((isAfterFac).And.(isVerb)) Then
       PRINFO(RINFO_TIMING_FacLcMat,"Timing -> Factorise -> lc. matrix")
       PRINFO(RINFO_TIMING_GetSchur,"Timing -> Factorise -> compute schur.")
    End If

    If((isAfterPcd).And.(isVerb)) Then
       PRINFO(RINFO_TIMING_AsmsSchur,"Timing -> Precond -> Assemble matrix")
       PRINFO(RINFO_TIMING_PcdFacto,"Timing -> Precond -> Facto matrix")
    End If

    If((isAfterSlv).And.(isVerb)) Then
       PRINFO(RINFO_TIMING_SolveDistRHS,"Timing -> Solve -> distribute b")
       PRINFO(RINFO_TIMING_SolveGenRHS,"Timing -> Solve -> compute Schur RHS")
       PRINFO(RINFO_TIMING_SolveITS,"Timing -> Solve -> the interface (SCHUR)")
       PRINFO(RINFO_TIMING_SolveSDS,"Timing -> Solve -> the interiors")
       PRINFO(RINFO_TIMING_SolveGathSol,"Timing -> Solve -> Gather x")
       Call PrintSH(unit,'Interpolation')
       
       PRINTERP(RINFO_TIMING_Interp,"Timing -> Solve -> Intepolation")
       PRINTERP(RINFO_TIMING_Fault,"Timing ->  Choice of faulty proc")
       PRINTERP(RINFO_RATIO_Lost,"Ratio of lost data")
       PRINTERP(RINFO_RATIO_Interp,"Ratio of interpolated data")
       PRINTERP(RINFO_PROC_INVOLVED,"Nb of procs involved in recovery")
    End If
    Call PrintTF4(unit,lend)

#undef PICNTL
#undef PRCNTL
#undef PIINFOG
#undef PRINFOG
#undef PIINFO1
#undef PIINFO
#undef PRINFO
#undef PRINTERP


    info = 0

  End Subroutine XMPH_state_print

  !> Print a section header
  Subroutine PrintSH(unit,msg)
    Integer, Intent(in) :: unit
    Character(len=*) :: msg
    Integer :: i

    Write(unit,*)
    Write(unit,*) '     * ',Trim(msg)
    Write(unit,*)
    ! Write(unit,*) '* ',Trim(msg)
    ! Write(unit,*) ('=',i=1,Len_Trim(msg)+2)

  End Subroutine PrintSH

  !> Print a table's header 
  !! where the table's line consisting have only 1 value
  Subroutine PrintTH1(unit,lend)
    Integer, Intent(in) :: unit
    Integer, Intent(in) :: lend
    
    Character(len=MAPHYS_STRL) :: stmp
    Integer :: i
    !-

    stmp = "Description"
    Write(unit,*)
    Write(unit,FMT='(A,A9,A5,A,A11)') &
         stmp(1:lend),"  FIELD (","INDEX",")","VALUE"
    Write(unit,*)  ('-', i=1,lend+15+11)


  End Subroutine PrintTH1

  !> Print a table's header 
  !! where the table's line consisting have
  !! 4 values (min,max,avg,std)
  Subroutine PrintTH4(unit,lend)
    Integer, Intent(in) :: unit
    Integer, Intent(in) :: lend

    Character(len=MAPHYS_STRL) :: stmp
    Integer :: i
    !-

    Write(unit,*)
    stmp = "Description"
    Write(unit,FMT='(A,A9,A5,A1,4A11)') &
         stmp(1:lend),"  FIELD (","INDEX",")",&
         'MIN','MAX','AVERAGE','STD.DEV.'
    Write(unit,*) ('-', i=1,lend+15+4*11)


  End Subroutine PrintTH4



  !> Print a table's footer
  !! where the table's line consisting have 1 value
  !!  values (min,max,avg,std)
  Subroutine PrintTF1(unit,lend)
    Integer, Intent(in) :: unit
    Integer, Intent(in) :: lend

    Write(unit,*)
    ! Integer :: i
    ! Write(unit,*) ('-', i=1,lend+15+11)


  End Subroutine PrintTF1

  !> Print a table's footer
  !! where the table's line consisting have
  !! 4 values (min,max,avg,std)
  Subroutine PrintTF4(unit,lend)
    Integer, Intent(in) :: unit
    Integer, Intent(in) :: lend

    Write(unit,*)
    
    ! Integer :: i
    ! Write(unit,*) ('-', i=1,lend+15+4*11)

  End Subroutine PrintTF4



  ! Print a table's line 
  Subroutine PrintTL(field,ind,lend,mphs,unit,msg)
    Use XMPH_maphys_type
    Implicit None

    Integer, Intent(in) :: field
    Integer, Intent(in) :: ind
    Integer, Intent(in) :: lend
    Type(XMPH_maphys_t), Intent(in) :: mphs
    Integer, Intent(in) :: unit
    Character(len=*), Intent(in) :: msg

    Character(len=MAPHYS_STRL) :: stmp

    !-
    stmp = Trim(msg) 
    Select Case(field)
       !- val
    Case(F_ICNTL); Write(unit,'(A,A9,I5,A,I11)')&
         stmp(1:lend),"  ICNTL (",ind,")",mphs%icntl(ind)

    Case(F_RCNTL); Write(unit,'(A,A9,I5,A,1PE11.3)')&
         stmp(1:lend),"  RCNTL (",ind,")",mphs%rcntl(ind)

    Case(F_IINFOG); Write(unit,'(A,A9,I5,A,I11)')&
         stmp(1:lend),"  IINFOG(",ind,")",mphs%iinfog(ind)

    Case(F_RINFOG); Write(unit,'(A,A9,I5,A,1PE11.3)') &
         stmp(1:lend),"  RINFOG(",ind,")",mphs%rinfog(ind)

    Case(F_IINFO1); Write(unit,'(A,A9,I5,A,I11)') &
         stmp(1:lend),"  IINFO (",ind,")",mphs%iinfo(ind)

       !- min, max, avg, sig 
    Case(F_IINFO ); Write(unit,'(A,A9,I5,A,2I11,2(1PE11.3))')&
         stmp(1:lend),"  IINFO (",ind,")",&
         mphs%iinfomin(ind), mphs%iinfomax(ind),&
         mphs%iinfoavg(ind), mphs%iinfosig(ind)

    Case(F_RINFO); Write(unit,'(A,A9,I5,A,4(1PE11.3))') &
         stmp(1:lend),"  RINFO (",ind,")",&
         mphs%rinfomin(ind), mphs%rinfomax(ind),&
         mphs%rinfoavg(ind), mphs%rinfosig(ind)

       ! special stuff
    Case(F_IsSYM   ); stmp = "Symmetry"; Write(unit,'(A,A9,6X,I11)') &
         stmp(1:lend),"  SYM    ", mphs%SYM
       
    Case(F_IsARITH); stmp = "Arithmetic"; Write(unit,'(A,A9,5X,A20)')&
         stmp(1:lend),"  .....  ",' XMPH_FLOAT '

    Case(F_IsVERSION); stmp = "Version"; Write(unit,'(A,A9,5X,A20)') &
         stmp(1:lend),"  VERSION",Trim(mphs%version)

    End Select
    
  End Subroutine PrintTL



  Subroutine PrintInterp(ind,lend,mphs,unit,msg)
    Use XMPH_maphys_type
    Implicit None

    Integer, Intent(in) :: ind
    Integer, Intent(in) :: lend
    Type(XMPH_maphys_t), Intent(in) :: mphs
    Integer, Intent(in) :: unit
    Character(len=*), Intent(in) :: msg

    Character(len=MAPHYS_STRL) :: stmp

    !-
    stmp = Trim(msg) 
    Write(unit,'(A,A9,I5,A,1PE11.3)') &
         stmp(1:lend),"  RINFO (",ind,")",&
         mphs%rinfomax(ind)

    
  End Subroutine PrintInterp



  !
  !-----------------------------------------------------------------------------
  !
  !> Print the legend
  !!
  !! the legend is printed 
  !! - on the solve step
  !! - by the host,
  !! - in verbose 2 and above.
  Subroutine XMPH_state_print_legend(mphs)
    Implicit None
    Type(XMPH_maphys_t), Intent(in) :: mphs

    Integer :: verb
    Integer :: unit
    !- End of header ----------------------------------------------------------- 

    If (mphs%ikeep(IKEEP_CURJOB)  /= CURJOB_IsSolve ) Return
    If (mphs%ikeep(IKEEP_MPIRANK) /= mphs%ikeep(IKEEP_HOSTRANK) ) Return 

    Call mph_log_GetVerbosity(verb)
    If (verb < MSG_VERBOSE2 ) Return
    Call mph_log_GetUnit(MSG_VERBOSE2,unit)

    Write(unit,*)
    Write(unit,*) "== ABREVIATIONS == "
    Write(unit,*)
    Write(unit,*) " A  denotes the input matrix"
    Write(unit,*) " b  denotes the input right-hand-side"
    Write(unit,*) " x  denotes the computed solution"
    Write(unit,*) "| | denotes the norm2"
    Write(unit,*) "STD.DEV. stands for standard deviation"
    Write(unit,*) "iters. stands for iterations"
    Write(unit,*) "nb.    stands for number"
    Write(unit,*) "lc.    stands for local"
    Write(unit,*) "mtx.   stands for matrix"
    Write(unit,*) "mem.   stands for memory"
    Write(unit,*) "facto. stands for factorisation"
    Write(unit,*) "schur. stands for the schur complement matrix"
    Write(unit,*) "prcnd. stands for the preconditioner"
    Write(unit,*)
    Write(unit,*) "== NOTES == "
    Write(unit,*)
    Write(unit,*) "All memory costs are in megabytes"
    Write(unit,*) "All times        are in seconds"
    Write(unit,*) " ---- denotes an unaccessible field to the user"
    Write(unit,*) " -1   denotes an unavailable  data"
    Write(unit,*) "MIN, MAX, AVERAGE and STD.DEV refer to"
    Write(unit,*) "the minimum, the maximum, the average of a value"
    Write(unit,*) "and its standard deviation"
    Write(unit,*) "on the MPI processes of mphs%COMM "


  End Subroutine XMPH_state_print_legend


  !
  !-----------------------------------------------------------------------------
  !
  !> Print the legend
  !!
  !! the legend is printed 
  !! - on the analysis step
  !! - by the host,
  !! - in verbose 2 and above.
  Subroutine XMPH_state_print_about(mphs)
    Implicit None
    Type(XMPH_maphys_t), Intent(in) :: mphs

    Integer :: verb
    Integer :: unit
    Integer :: i, len
    Character(len=MAPHYS_STRL) :: msg

    !- End of header ----------------------------------------------------------- 

    If (mphs%ikeep(IKEEP_CURJOB)  /= CURJOB_IsAnalysis ) Return
    If (mphs%ikeep(IKEEP_MPIRANK) /= mphs%ikeep(IKEEP_HOSTRANK) ) Return 

    Call mph_log_GetVerbosity(verb)
    If (verb < MSG_STD ) Return
    Call mph_log_GetUnit(MSG_STD,unit)
    If (unit < 0) Return


    Write(msg,'(20X,5A)') "* MaPHyS ", Trim(mphs%Version),' [', Trim(' XMPH_FLOAT '), ' ] *'
    len = Len_Trim(msg) - 20

    Write(unit,*)
    Write(unit,*)
    Write(unit,*) (' ',i=1,20 ),('*',i=1,len )
    Write(unit,*) Trim(msg)
    Write(unit,*) (' ',i=1,20 ),('*',i=1,len )
    Write(unit,*)

  End Subroutine XMPH_state_print_about



End Module XMPH_state_print_mod

