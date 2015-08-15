!> utilitary module for unit testing
!!
!! This module defines useful routines and variables for unit testing.
!!
!! Provides the routine "Test" which :
!!   - launch a procedure (without arguments)
!!   - log the result of the procedures
!!
!! The procedure talks to "Test" by setting 2 variables namely :
!!   - test_result, the result of the procedure
!!                  (integer, valid values are in enumeration TEST_*) 
!!   - test_msg, the message of the procedure (string of size 1024)
!!
!!----
!! 
!! @par History
!!
!! - Date     : version : comments
!! - 18/01/11 : v0.1a   : add test_result
!! - 1*/01/11 : v0.1a   : creation
!!
Module test_utils_mod

  Implicit None

  !* Shared Public data *!
  Integer, Public, Save  :: test_result
  Character(len=1024)    :: test_msg

  !* Enumerations *!
  Integer, Parameter     :: TEST_PASS  =  0
  Integer, Parameter     :: TEST_FAIL  = -1
  Integer, Parameter     :: TEST_KFAIL =  1
  Integer, Parameter     :: TEST_KPASS =  2
  Integer, Parameter     :: TEST_XFAIL =  3
  Integer, Parameter     :: TEST_XPASS = -2
  Integer, Parameter     :: TEST_SKIP  =  5
  Integer, Parameter     :: TEST_UNSET =  6

  !* Local data *!
  Integer, Private, Save :: notest = 0
  Integer, Private, Save :: countpass = 0
  Integer, Private, Save :: countfail = 0
  Integer, Private, Save :: countothr = 0
  
  !* Access specifiers *!
  Public  :: test 

Contains

  !> launch a test
  !! 
  !! Launch the test "sub" (which should set "test_result" to TEST_* )
  !! Print a report in the form 
  !! Testing : <no>  "text of the test " [<test_result>]
  !! 
  !!--- 
  !!
  !! @param [in,out] sub the subroutine to call (should set test_status)
  !! @param [in    ] txt the text of the test.
  !!
  Subroutine test( sub, text )

    ! Arguments
    External sub
    Character(len=*)  :: text

    ! explicit interface 
    Interface 
       Subroutine sub
       End Subroutine sub
    End Interface
    
    ! Local variables
    Character(len=15), Parameter  :: fmt = '(a,i4.3,a60,a)'
    Integer :: status
    Logical :: haveMPI, performPrinting
    Integer :: mpirank, mpisize
    Include "mpif.h"

    !---
    Call MPI_Initialized(haveMPI,status)
    If (haveMPI)Then
       Call MPI_Comm_rank(MPI_COMM_WORLD,mpirank,status)
       Call MPI_Comm_size(MPI_COMM_WORLD,mpisize,status)
       performPrinting = .False.
       If (mpirank == 0) performPrinting = .True.
    Else
       Write(*,*) "****KO****"
       performPrinting = .True.
    End If

    notest = notest+1

    test_msg = ""
    test_result = TEST_UNSET

    ! Print the test that will be launched
    If (performPrinting) &
         Write(6,FMT='(a,i4.3,a60)',advance="no") "Test : ",notest, Trim(text)

    ! Exercise the test
    Call sub

    ! Print the result of the test
    If (.Not. performPrinting) Return

    Select Case(test_result)
    Case( TEST_PASS  ) ; Write(6,FMT='(a)',advance="no") " [PASS ]"
    Case( TEST_FAIL  ) ; Write(6,FMT='(a)',advance="no") " [FAIL ]"
    Case( TEST_KFAIL ) ; Write(6,FMT='(a)',advance="no") " [KFAIL]"
    Case( TEST_KPASS ) ; Write(6,FMT='(a)',advance="no") " [KPASS]"
    Case( TEST_XFAIL ) ; Write(6,FMT='(a)',advance="no") " [XFAIL]"
    Case( TEST_XPASS ) ; Write(6,FMT='(a)',advance="no") " [XPASS]"
    Case( TEST_SKIP  ) ; Write(6,FMT='(a)',advance="no") " [SKIPPED]"
    Case Default       ; Write(6,FMT='(a)',advance="no") " [UNK  ]"
    End Select

    ! Increment global counters
    Select Case(test_result)
    Case( TEST_PASS  ) ; countpass = countpass +1
    Case( TEST_FAIL  ) ; countfail = countfail +1
    Case( TEST_KFAIL ) ; countothr = countothr +1
    Case( TEST_KPASS ) ; countothr = countothr +1
    Case( TEST_XFAIL ) ; countothr = countothr +1
    Case( TEST_XPASS ) ; countothr = countothr +1
    Case( TEST_SKIP  ) ; countothr = countothr +1
    Case Default       ; countothr = countothr +1
    End Select

    ! Print additional information about its execution
    Write(6,FMT='(a)') Trim(test_msg)
    
  End Subroutine test

  Subroutine test_summarize
    Integer :: countall
    countall = countpass + countfail + countothr
    Write(6,FMT='(a)') ""
    Write(6,FMT='(a)') "Summary of Launched tests"
    Write(6,FMT='(a,I10,a,I10)') "PASSED = ", countpass, "/", countall
    Write(6,FMT='(a,I10,a,I10)') "FAILED = ", countfail, "/", countall
    Write(6,FMT='(a,I10,a,I10)') "OTHERS = ", countothr, "/", countall
    Write(6,FMT='(a)')"OTHERS = SKIPPED tests or those returning {K,X}{PASS,FAIL}"
  End Subroutine test_summarize


End Module test_utils_mod
