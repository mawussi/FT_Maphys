#include "mph_defs_f.h"
#include "mph_macros_f.h"
!> module containing the routines to partition/distribute the system. 
!!
!! @author Azzam Haidar
!! @author Yohan Lee-tin-yien
!!
Module MPH_domain_mod

  !* Imported types *!
  Use MPH_error_mod
  Use MPH_part_type, Only : &
       maphys_domain_t

  !* Defined routines *!
  Public :: MPH_domain_nullify
  Public :: MPH_domain_free
  Public :: MPH_domain_write
  Public :: MPH_domain_read
  Public :: MPH_domain_check
  Public :: MPH_domain_getBoundary

  !* Private constants *!
  Character(len=MAPHYS_STRL), Private, Parameter :: &
       FLNAME = "mph_domain_mod.F90"

  Contains

    !****
    
    Subroutine MPH_domain_nullify(this)
      Type(maphys_domain_t), Intent(inout) :: this
      Nullify(this%myindexVi , this%myindexintrf ,&
           this%myinterface  , this%mylogicintrf , &
           this%myptrindexVi , this%weight )
    End Subroutine MPH_domain_nullify

    Subroutine MPH_domain_free(this)
      Type(maphys_domain_t), Intent(inout) :: this
      If(Associated(this%myindexVi    )) Deallocate(this%myindexVi )
      If(Associated(this%myindexintrf )) Deallocate(this%myindexintrf)
      If(Associated(this%myinterface  )) Deallocate(this%myinterface  ) 
      If(Associated(this%mylogicintrf )) Deallocate(this%mylogicintrf)
      If(Associated(this%myptrindexVi )) Deallocate(this%myptrindexVi)
      If(Associated(this%weight       )) Deallocate(this%weight)
    End Subroutine MPH_domain_free


    !****

    Subroutine MPH_domain_write(domain, wunit, info)
      Implicit None

      !* Arguments *!

      Type(maphys_domain_t), Intent(in) :: domain
      Integer, Intent(in) :: wunit
      Integer, Intent(out) :: info
      
      !- End of header ---------------------------------------------------------
      
      Write(wunit,*) "ndof=", domain%myndof
      Write(wunit,*) "ndofinterior=", domain%myndofinterior
      Write(wunit,*) "ndofintrf=", domain%myndofintrf
      Write(wunit,*) "gballintrf=", domain%gballintrf
      Write(wunit,*) 
      Write(wunit,*) "int_Lg=",domain%myint_Lg    
      Write(wunit,*) 
      Write(wunit,*) "IndexIntfMaxSize=",domain%lenindintrf 
      Write(wunit,*) "sizeIntrf=",domain%mysizeIntrf 
      Write(wunit,*) "nbVi=",domain%mynbvi      
      Write(wunit,*) "indexVi(:)=",domain%myindexVi(1:domain%mynbvi)
      Write(wunit,*) "ptrIndexVi(:)=",domain%myptrindexVi(1:domain%mynbvi+1)
      Write(wunit,*) "IndexIntrf(:)=",&
           domain%myindexintrf(  1  : domain%myptrindexVi(domain%mynbvi+1) -1  )
      ! Nullify(domain%myinterface) ! list_Intrfc ! (should be removed)

      Write(wunit,*) 
      Write(wunit,*) "ndofLogicIntf=", domain%myndoflogicintrf
      Write(wunit,*) "LogicIntf(:)=", domain%mylogicintrf(1:domain%myndoflogicintrf)
      Write(wunit,*) "weight(:)=", domain%weight(1:domain%mysizeIntrf)

      info = 0

    End Subroutine MPH_domain_write

    Subroutine MPH_domain_read(this, runit, info)
      Implicit None

      !* Arguments *!

      Type(maphys_domain_t), Intent(out) :: this
      Integer, Intent(in) :: runit
      Integer, Intent(out) :: info
      
      !* Local variables *!
      Character(len=MAPHYS_STRL) :: dummy
      
      
      !- End of header ---------------------------------------------------------
      
      Read(runit,*) dummy, this%myndof
      Read(runit,*) dummy, this%myndofinterior
      Read(runit,*) dummy, this%myndofintrf
      Read(runit,*) dummy, this%gballintrf
      Read(runit,*) 
      Read(runit,*) dummy,this%myint_Lg    
      Read(runit,*) ! line jump
      Read(runit,*) dummy,this%lenindintrf 
      Read(runit,*) dummy,this%mysizeIntrf 
      Read(runit,*) dummy,this%mynbvi     

      Call allocandread(this%myindexVi,this%mynbvi,info) 
      MPH_ONFAILURE_RETURN(info)

      Call allocandread(this%myptrindexVi,this%mynbvi+1,info) 
      MPH_ONFAILURE_RETURN(info)

      Call allocandread(&
           this%myindexintrf,this%myptrindexVi(this%mynbvi+1)-1,info) 
      MPH_ONFAILURE_RETURN(info)
      ! Nullify(this%myinterface) ! list_Intrfc ! (should be removed)

      Read(runit,*) ! line jump
      Read(runit,*) dummy, this%myndoflogicintrf
      Call allocandread(this%mylogicintrf,this%myndoflogicintrf,info)
      MPH_ONFAILURE_RETURN(info)

      Call allocandreadr8(this%weight,this%mysizeIntrf,info)
      MPH_ONFAILURE_RETURN(info)

      ! small internal routines
      Contains

        Subroutine allocandread(ptr,ptrsize,info)
          MPH_INT, Pointer :: ptr(:)
          MPH_INT :: ptrsize
          Integer :: info
          !
          Allocate(ptr(ptrsize),STAT=info)
          CHCKALLOC(info)
          Read(runit,*) dummy, ptr(1:ptrsize)

        End Subroutine allocandread

        Subroutine allocandreadr8(ptr,ptrsize,info)
          Real(Kind=8), Pointer :: ptr(:)
          MPH_INT :: ptrsize
          Integer :: info
          !
          Allocate(ptr(ptrsize),STAT=info)
          CHCKALLOC(info)
          Read(runit,*) dummy, ptr(1:ptrsize)
        End Subroutine allocandreadr8

    End Subroutine MPH_domain_read


    
    
    Subroutine MPH_domain_check(domain, info)
      Implicit None
      
      !* Arguments *!
      
      Type(maphys_domain_t), Intent(in) :: domain
      Integer, Intent(out) :: info

      !* Local variables *! 

      Integer, Parameter :: nchecks = 6
      Integer :: checks(nchecks)

      !- End of header ---------------------------------------------------------

      checks = 0

      CHCKASSRT( domain%mysizeIntrf > 0            , checks(1) )
      CHCKASSRT( domain%mynbvi > 0                 , checks(2) )
      CHCKASSRT( domain%mynbvi > 0                 , checks(3) )
      CHCKASSRT( Associated( domain%myindexVi )    , checks(4) )
      CHCKASSRT( Associated( domain%myptrindexVi ) , checks(5) )
      CHCKASSRT( Associated( domain%myindexintrf ) , checks(6) )

      info = MinVal(checks)

    End Subroutine MPH_domain_check


    !> Get the boundary with a neighbor
    !!
    !! Associate the pointer "Bound" to the boundary with neighbor "neigh".
    !! On error, return the bound => Null() , and boundsize = -1
    !! 
    !! @param [in ] domain    the domain description
    !! @param [in ] neigh     the neighbor
    !! @param [out] bound     the boundary
    !! @param [out] boundsize the boundary size
    !!
    Subroutine MPH_domain_getBoundary(domain, neigh,bound, boundsize )

      Implicit None

      !* Arguments *!

      Type(maphys_domain_t), Intent(in) :: domain
      Integer              , Intent(in) :: neigh
      MPH_INT, Pointer  , Intent(out) :: bound(:)
      MPH_INT           , Intent(out) :: boundsize
      
      !* Local variables *!

      MPH_INT, Pointer :: ptr2bounds(:)
      MPH_INT, Pointer :: bounds(:)

      !- End of header ---------------------------------------------------------
      
      ! Exit early
      boundsize = -1
      Nullify(bound)
      If ((neigh < 1).OR.(neigh>domain%mynbvi)) Return

      !
      Nullify( ptr2Bounds, Bounds )

      ptr2Bounds => domain%myPtrIndexVi
      Bounds     => domain%myIndexIntrf

      boundsize = ptr2Bounds(neigh+1) - ptr2Bounds(neigh)
      bound => Bounds(ptr2Bounds(neigh):ptr2Bounds(neigh+1)-1)

    End Subroutine MPH_domain_getBoundary
      

    !> Get the size of the boundary with a neighbor
    !!
    !! On error, return boundsize = -1
    !! 
    !! @param [in ] domain    the domain description
    !! @param [in ] neigh     the neighbor
    !! @return  boundsize the boundary size
    !!
    MPH_INT Function MPH_domain_getBoundarySize(domain, neigh ) &
         Result (boundsize)
      
      Implicit None
      
      !* Arguments *!
      
      Type(maphys_domain_t), Intent(in) :: domain
      Integer              , Intent(in) :: neigh
      
      boundsize = -1
      
      If ((neigh < 1).OR.(neigh>domain%mynbvi)) Return
      
      boundsize = domain%myptrindexVi(neigh+1) - domain%myptrindexVi(neigh)

    End Function MPH_domain_getBoundarySize
    


End Module MPH_domain_mod
