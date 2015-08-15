#include "mph_defs_f.h"
#include "mph_macros_f.h"
!> module to manipulate the mph_rhs_partition structure
Module MPH_rhs_partition_mod
  
  Use mph_error_mod
  Use mph_part_type, Only : &
       maphys_rhs_partition_t

  Implicit None
  
  !* Public routines *!
  Public :: mph_rhspart_null
  Public :: mph_rhspart_free

  !* Private routines *!
  Private :: ReadI
  Private :: ReadIptr

  !* Private constants *!
  Character(len=MAPHYS_STRL), Private, Parameter :: &
       FLNAME = "mph_rhs_partition_mod.F90"


  Contains

    Subroutine MPH_rhspart_null(this)
      Implicit None
      Type(maphys_rhs_partition_t), Intent(inout) :: this

      Nullify(this%procintdisp )      
      Nullify(this%procintrfdisp  )   
      
      Nullify(this%procintdisp )      
      Nullify(this%procintrfdisp  )   
      Nullify(this%domLg )            
      Nullify(this%domintdof )        
      Nullify(this%domintrfdof )     
      Nullify(this%metperm )          
      
      Nullify(this%scatindices )
      Nullify(this%scatlogicindices )
      Nullify(this%domlogicintrfdof)

      Nullify(this%gbtoloc)

    End Subroutine MPH_rhspart_null
    !
    !-------------------------------------------------------------------------
    !    
    Subroutine MPH_rhspart_free(this)
      Implicit None
      Type(maphys_rhs_partition_t), Intent(inout) :: this

      If(Associated(this%procintdisp      ))Deallocate(this%procintdisp)      
      If(Associated(this%procintrfdisp    ))Deallocate(this%procintrfdisp )   
      If(Associated(this%domLg            ))Deallocate(this%domLg)            
      If(Associated(this%domintdof        ))Deallocate(this%domintdof)        
      If(Associated(this%domintrfdof      ))Deallocate(this%domintrfdof)     
      If(Associated(this%metperm          ))Deallocate(this%metperm)          
      If(Associated(this%scatindices      ))Deallocate(this%scatindices)      
      If(Associated(this%scatlogicindices ))Deallocate(this%scatlogicindices)
      If(Associated(this%domlogicintrfdof ))Deallocate(this%domlogicintrfdof)
      If(Associated(this%gbtoloc          ))Deallocate(this%gbtoloc)
      
    End Subroutine MPH_rhspart_free
    !
    !-------------------------------------------------------------------------
    !    
    Subroutine MPH_rhspart_read(this,funit,info)
      Implicit None
      
      Type(maphys_rhs_partition_t), Intent(out) :: this
      Integer, Intent(in) :: funit
      Integer, Intent(out) :: info

      !- End of header ---------------------------------------------------------
      Call MPH_rhspart_null(this)

      Call ReadI(this%nbdom,funit,info) ! line 1
      MPH_ONFAILURE_GOTO9999(info)
      Call ReadI(this%gballndof,funit,info) ! line 2
      MPH_ONFAILURE_RETURN(info)
      Call ReadI(this%gballintrf,funit,info) ! line 3
      MPH_ONFAILURE_RETURN(info)
      Call ReadI(this%combivtxpcsz,funit,info) ! line 4
      MPH_ONFAILURE_GOTO9999(info)
      Call ReadIptr(this%procintdisp  ,this%nbdom,funit,info) ! line 5
      MPH_ONFAILURE_GOTO9999(info)
      Call ReadIptr(this%procintrfdisp,this%nbdom,funit,info) ! line 6
      MPH_ONFAILURE_GOTO9999(info)
      ! Call ReadIptr(this%domLg        ,this%nbdom,funit,info_ignore)
      ! 
      Call ReadIptr(this%domintdof   ,this%nbdom,funit,info) ! line 7
      MPH_ONFAILURE_GOTO9999(info)
      Call ReadIptr(this%domintrfdof ,this%nbdom,funit,info) ! line 8
      MPH_ONFAILURE_GOTO9999(info)
      Call ReadIptr(this%metperm      ,this%gballndof,funit,info) !  line 9
      MPH_ONFAILURE_GOTO9999(info)
      Call ReadI(this%rhsway,funit,info) ! line 10
      Select Case (this%rhsway)
      Case(1) 

         Call ReadIptr(this%scatindices,this%combivtxpcsz,funit,info) ! line 11
         MPH_ONFAILURE_GOTO9999(info)

      Case(2)

           Call ReadIptr(this%scatlogicindices,this%gballintrf,funit,info) ! line 11
           MPH_ONFAILURE_GOTO9999(info)
           Call ReadIptr(this%domlogicintrfdof,this%nbdom,funit,info) ! line 12
           MPH_ONFAILURE_GOTO9999(info)

      Case Default

         info = -1
         MPH_ONFAILURE_GOTO9999(info)

      End Select

      Call ReadIptr(this%gbtoloc,this%gballintrf,funit,info) ! line 12/13
      MPH_ONFAILURE_RETURN(info)

9999 Continue

      If (info < 0) Call MPH_rhspart_free(this)
      
    End Subroutine MPH_rhspart_read
    !
    !-------------------------------------------------------------------------
    !
    Subroutine ReadI(i,funit,info)
      Implicit None
      Integer,  Intent(out) :: i
      Integer,  Intent(in) :: funit
      Integer,  Intent(out) :: info

      Character(len=MAPHYS_STRL) :: dummy

      Read(funit,*,IOSTAT=info) dummy, i
      CHCKIO(info)

    End Subroutine ReadI
    !
    !-------------------------------------------------------------------------
    !        
    Subroutine ReadIptr(ptr,ptrsize,funit,info)
      Implicit None
      Integer, Pointer, Intent(out) :: ptr(:)
      Integer, Intent(in) :: ptrsize
      Integer,  Intent(in) :: funit
      Integer,  Intent(out) :: info

      Character(len=MAPHYS_STRL) :: dummy

      Allocate(ptr(ptrsize),STAT=info)
      CHCKALLOC(info)
      If (info<0) Return
      
      Read(funit,*,IOSTAT=info) dummy, ptr(1:ptrsize)
      CHCKIO(info)
      If (info<0) Deallocate(ptr)

    End Subroutine ReadIptr
    !
    !-------------------------------------------------------------------------
    !    
    Subroutine MPH_rhspart_write(this,funit,info)
      Implicit None
      
      Type(maphys_rhs_partition_t), Intent(in) :: this
      Integer, Intent(in) :: funit
      Integer, Intent(out) :: info

      !- End of header ---------------------------------------------------------
      
      Call WriteI(this%nbdom,funit,info,"nbdom")
      MPH_ONFAILURE_RETURN(info)

      Call WriteI(this%gballndof,funit,info,"gballndof")
      MPH_ONFAILURE_RETURN(info)

      Call WriteI(this%gballintrf,funit,info,"gballintrf")
      MPH_ONFAILURE_RETURN(info)

      Call WriteI(this%combivtxpcsz,funit,info,"combivtxpcsz")
      MPH_ONFAILURE_RETURN(info)

      Call WriteIptr(this%procintdisp  ,this%nbdom,funit,info,"procintdisp")
      MPH_ONFAILURE_RETURN(info)

      Call WriteIptr(this%procintrfdisp,this%nbdom,funit,info,"procintrfdisp")
      MPH_ONFAILURE_RETURN(info)

      ! Call WriteIptr(this%domLg        ,this%nbdom,funit,info_ignore,"domLg")
      ! MPH_ONFAILURE_DO_NOTHING(info_ignore)


      Call WriteIptr(this%domintdof   ,this%nbdom,funit,info,"domintdof")
      MPH_ONFAILURE_RETURN(info)

      Call WriteIptr(this%domintrfdof ,this%nbdom,funit,info,"domintrfdof")
      MPH_ONFAILURE_RETURN(info)

      Call WriteIptr(this%metperm      ,this%gballndof,funit,info,"metperm")
      MPH_ONFAILURE_RETURN(info)
      !**
      Call WriteI(this%rhsway,funit,info,"rhsway")
      Select Case (this%rhsway)
      Case(1) 
         
         Call WriteIptr(this%scatindices,this%combivtxpcsz,funit,info,&
              "scatindices")
         MPH_ONFAILURE_RETURN(info)

      Case(2)
         Call WriteIptr(this%scatlogicindices,this%gballintrf,funit,info,&
              "scatlogicindices")         
         MPH_ONFAILURE_RETURN(info)
         Call WriteIptr(this%domlogicintrfdof,this%nbdom,funit,info,&
              "domlogicintrfdof")
         MPH_ONFAILURE_RETURN(info)
         
      Case Default

         info = -1
         MPH_ONFAILURE_RETURN(info)

      End Select

      Call WriteIptr(this%gbtoloc,this%gballintrf,funit,info,"gbtoloc")
      MPH_ONFAILURE_RETURN(info)


    End Subroutine MPH_rhspart_write
    !
    !-------------------------------------------------------------------------
    !
    Subroutine WriteIptr(ptr,ptrsize,funit,info,desc)
      Implicit None
      Integer, Pointer, Intent(in) :: ptr(:)
      Integer, Intent(in) :: ptrsize
      Integer,  Intent(in) :: funit
      Integer,  Intent(out) :: info
      Character(len=*), Intent(in) :: desc
      
      If (size(ptr) >= ptrsize)Then
         Write(funit,*,IOSTAT=info) Trim(desc)//":", ptr(1:ptrsize)
         CHCKIO(info)
      Else
         Write(funit,*,IOSTAT=info) Trim(desc)//":"
         CHCKIO(info)
         info = Max(info, - __LINE__ )
      End If

    End Subroutine WriteIptr
    !
    !-------------------------------------------------------------------------
    !
    Subroutine WriteI(i,funit,info,desc)
      Implicit None
      Integer,  Intent(in) :: i
      Integer,  Intent(in) :: funit
      Integer,  Intent(out) :: info
      Character(len=*), Intent(in) :: desc

      Character(len=MAPHYS_STRL) :: dummy

      Write(funit,*,IOSTAT=info) Trim(desc)//":", i
      CHCKIO(info)

    End Subroutine WriteI




End Module MPH_rhs_partition_mod


