! Warning: XMPH_GENFILE_COMMENT
#include "mph_defs_f.h"
#include "mph_macros_f.h"  
! utilitary routines for strategy manipulation
Module mph_strat_utils_mod
  Use MPH_Log_mod
  Implicit None

  Public :: unsetI
  Public :: unsetR
  Public :: is_seti
  Public :: is_setr
  Public :: is_unseti
  Public :: is_unsetr
  Public :: is_apowerOf2
  Public :: checkRange
  Public :: checkValue
  Public :: checkPositiveI
  Public :: checkPositiveR
  Public :: printInconsistencyI

  ! Private constants
  Character(len=MAPHYS_STRL), Parameter, Private :: &
       FLNAME= "xmph_strat_utils_mod.F90"
  Integer, Parameter, Private :: IUNSET = -99999 
  Real(Kind=8), Parameter, Private :: RUNSET = -99999 

Contains 
  !*****
  Subroutine unseti(param)
    Integer, Intent(out) :: param
    param = IUNSET
  End Subroutine unseti
  !*****
  Subroutine unsetr(param)
    Real(Kind=8), Intent(out) :: param
    param = RUNSET
  End Subroutine unsetr
  !****
  Logical Function is_seti( param )
    Integer, Intent(in) :: param
    is_seti = ( param /= IUNSET )
  End Function is_seti
  !****
  Logical Function is_setr( param )
    Real(Kind=8), Intent(in) :: param
    is_setr = ( param /= RUNSET )
  End Function is_setr
  !****
  Logical Function is_unseti( param )
    Integer, Intent(in) :: param
    is_unseti = ( param == IUNSET )
  End Function is_unseti
  !****
  Logical Function is_unsetr( param )
    Real(Kind=8), Intent(in) :: param
    is_unsetr = ( param == RUNSET )
  End Function is_unsetr
  !****
  Logical Function is_aPowerOf2(int)
    Implicit None
    Integer, Intent(in) :: int
    ! check if "int" is a power of 2
    Integer :: i
    is_APowerOf2 = .False. 
    i=1
    Do While ( i < int )
       i = i*2
       If ( i == int ) Then
          is_APowerOf2 = .True. 
          Return
       End If
    End Do
  End Function is_APowerOf2
  !****
  Logical Function checkRange(icntl,ind, min, max )
    Implicit None
    Integer, Intent(in) :: icntl(MAPHYS_ICNTL_SIZE)
    Integer, Intent(in) :: ind, min, max
    !  check if the incl(ind) is in the range [min,max]
    Character(len=MAPHYS_STRL) :: msg
    checkRange = ( min <= icntl(ind) ).And.(icntl(ind) <= max )
    If (.Not. checkRange )Then
       Write(msg,'(A,I2,A,3(A,I8))') "ICNTL(",ind,") is out of range:", &
            " val =", icntl(ind), " min=", min, " max=", max
       Call MPH_Log(MSG_ERROR,msg)
    End If

  End Function checkRange
  !**** 
  Logical Function checkValue(icntl,ind, val )
    Implicit None
    Integer, Intent(in) :: icntl(MAPHYS_ICNTL_SIZE)
    Integer, Intent(in) :: ind, val
    !  check if the intcl(ind) has the value val
    Character(len=MAPHYS_STRL) :: msg
    checkValue = ( icntl(ind) == val )
    If (.Not. checkValue )Then
       Write(msg,'(A,I2,A,2(A,I8))') "ICNTL(",ind,") wrong value:", &
            " val =", icntl(ind), " must be =", val
       Call MPH_Log(MSG_ERROR,msg)
    End If

  End Function checkValue
  !***
  Logical Function checkPositivei(icntl,ind)
    Implicit None
    Integer, Intent(in) :: icntl(MAPHYS_ICNTL_SIZE)
    Integer, Intent(in) :: ind
    !  check if the incl(ind) is in the range [min,max]
    Character(len=MAPHYS_STRL) :: msg
    checkPositivei = ( 0 <= icntl(ind) )
    If (.Not. checkPositivei )Then
       Write(msg,'(A,I2,A,1(A,I8))') "ICNTL(",ind,") must be positive:",&
            " val =", icntl(ind)
       Call MPH_Log(MSG_ERROR,msg)
    End If

  End Function checkPositivei
  !*** 
  Logical Function checkPositiver(rcntl,ind)
    Implicit None
    Real(Kind=8), Intent(in) :: rcntl(MAPHYS_RCNTL_SIZE)
    Integer, Intent(in) :: ind
    !  check if the incl(ind) is in the range [min,max]
    Character(len=MAPHYS_STRL) :: msg
    checkPositiveR = ( 0.d0 <= rcntl(ind) )
    If (.Not. checkPositiver )Then
       Write(msg,'(A,I2,A,A,G15.3))') "RCNTL(",ind,") must be positive:",&
            " val =", rcntl(ind)
       Call MPH_Log(MSG_ERROR,msg)
    End If

  End Function checkPositiver
  !***
  Subroutine printInconsistencyI(icntl,ind1,ind2)
    Implicit None
    Integer, Intent(in) :: icntl(MAPHYS_ICNTL_SIZE)
    Integer, Intent(in) :: ind1,ind2
    Character(len=MAPHYS_STRL) :: msg

    Write(msg,'(A,I2,A,I2,A)') "Inconsitency detected between &
         &ICNTL(",ind1,") and ICNTL(",ind2,")."
    Call MPH_Log(MSG_ERROR,msg)
    msg = ""
    Write(msg,'(A,I5,A,I5)') "Please correct their value. Respective values are: ",&
         icntl(ind1)," and ", icntl(ind2)
    Call MPH_Log(MSG_ERROR,msg)

  End Subroutine printInconsistencyI
End Module mph_strat_utils_mod
