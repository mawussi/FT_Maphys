#include "mph_defs_f.h"
#include "mph_macros_f.h"

! [+] Module : MPH_part_scotch_mod ----------------------------------------------------
!
!> Module for partitioning with scotch
!!
Module MPH_part_scotch_mod

  !* Modules *!
  Use MPH_part_type
  
  !* No implicit typing *!
  Implicit None

#ifdef HAVE_LIBSCOTCH
  !* Import scotch constants *!
  Include 'scotchf.h'
#endif

  !* Enumerations *!
  Integer, Parameter :: MAPHYS_NUMKIND = MPH_INTKIND
!  Integer, Parameter :: SCOTCH_NUMKIND = 8  ! SCOTCH with 64 bits integer
  Integer, Parameter :: SCOTCH_NUMKIND = 4   ! SCOTCH with 32 bits integers
  Integer, Parameter :: SCOTCH_IDXKIND = 8
  Integer, Parameter :: SCOTCH_SYSKIND = 4

  !* structures *!

#ifndef HAVE_LIBSCOTCH
  Type maphys_scotch_graph_t
     sequence
     Integer :: dummy
  End type maphys_scotch_graph_t
#else
  Type maphys_scotch_graph_t
     sequence
     ! scalars
     Integer(kind=SCOTCH_NUMKIND) :: baseval
     Integer(kind=SCOTCH_NUMKIND) :: vertnbr
     Integer(kind=SCOTCH_NUMKIND) :: edgenbr
     
     ! array
     Integer(kind=SCOTCH_NUMKIND), Pointer :: verttab(:)
     Integer(kind=SCOTCH_NUMKIND), Pointer :: vendtab(:)
     Integer(kind=SCOTCH_NUMKIND), Pointer :: velotab(:)
     Integer(kind=SCOTCH_NUMKIND), Pointer :: edgetab(:)
     
     ! opaque scotch structure 
     Real(kind=8) :: data(SCOTCH_GRAPHDIM)

  End type maphys_scotch_graph_t
#endif

  !* Access specifiers *!
  
  Public  :: Part_scotch_partgraph
  Public  :: Maphys2scotch_graph
  Public  :: scotch2maphys_part_tree
  Private :: get_sons


  !* Procedures *!

  Contains

    ! [+] routine : PART_scotch_partgraph ---------------------------------
    !
    !> partition the matrix graph with scotch
    !!
    !!----
    !!
    !! @param metalgo     choosen algorithm
    !! @param nbdom       wanted number of domain
    !! @param azz         the binary tree
    !! @param metperm     the permutation
    !! @param metiperm    the inverse permutation
    !! @param info        status of the subroutine (0=SUCCESS)
    !!
    !!----
    !!
    !! @author Yohan Lee-tin-yien
    !!
    !!
    subroutine PART_scotch_partgraph              &  ! intents
         ( metalgo, nbdom, graph_global_A,      &  ! in
         & azz, metperm,metiperm, info )           ! out

      !* Modules *!

      Use MPH_part_type
      Use mph_log_mod
      Implicit None

      !* Arguments *!
      integer   , intent(in) :: metalgo
      integer   , intent(in) :: nbdom
      type(maphys_matrix_graph_t), intent(in) :: graph_global_A
      type(maphys_binary_node_t) , pointer, intent(  out) :: azz (:)
      integer, intent(out), pointer :: metperm  (:)
      integer, intent(out), pointer :: metiperm (:)
      integer, intent(out)          :: info

#ifndef HAVE_LIBSCOTCH
      Call MPH_Log(MSG_ERROR,"MaPHyS was not compiled with SCOTCH (-DHAVE_LIBSCOTCH)")
      info = -1
      Return
#else
      !* Local variables *!
      
      ! constants
      Integer, Parameter :: STRAT_LEN = 10000
!       Character(len=STRAT_LEN), Parameter :: strattemplate = "n{&
!            &sep=(/((levl)<(MAPHYSMAXLEVEL))?((m{asc=b{bnd=f{move=200,pass=1000,bal=0.5},&
!            &org=(|h{pass=10})f{move=200,pass=1000,bal=0.5},width=3},low=h{pass=10},&
!            &type=h,vert=100,rat=0.7}|m{asc=b{bnd=f{move=200,pass=1000,bal=0.5},&
!            &org=(|h{pass=10})f{move=200,pass=1000,bal=0.5},width=3},low=h{pass=10},&
!            &type=h,vert=100,rat=0.7}));)&
!            &,ole=g{pass=3}&
!            &,ose=g{pass=3}}"

      Character(len=STRAT_LEN), Parameter :: strattemplate = "n{&
           &sep=(/((levl)<(MAPHYSMAXLEVEL))?((m{asc=b{bnd=f{move=200,pass=1000,bal=0.01},&
           &org=(|h{pass=10})f{move=200,pass=1000,bal=0.01},width=3},low=h{pass=10},&
           &type=h,vert=100,rat=0.7}|m{asc=b{bnd=f{move=200,pass=1000,bal=0.01},&
           &org=(|h{pass=10})f{move=200,pass=1000,bal=0.01},width=3},low=h{pass=10},&
           &type=h,vert=100,rat=0.7}));)&
           &,ole=g{pass=3}&
           &,ose=g{pass=3}}"


      ! scalars
      Integer :: istep
      Integer :: msg_class
      Integer :: maphysMaxLevel
      Integer(kind=SCOTCH_SYSKIND) :: iinfo 
      Integer(kind=SCOTCH_NUMKIND) :: vertnbr
      Integer(kind=SCOTCH_NUMKIND) :: cblknbr
      Integer(kind=SCOTCH_NUMKIND) :: i,j

      ! arrays     
      Integer(kind=SCOTCH_NUMKIND), Pointer :: permtab (:)
      Integer(kind=SCOTCH_NUMKIND), Pointer :: peritab (:)
      Integer(kind=SCOTCH_NUMKIND), Pointer :: rangtab (:)
      Integer(kind=SCOTCH_NUMKIND), Pointer :: treetab (:)

      ! types
      Type(maphys_scotch_graph_t) :: scotch_graph
      Real(kind=8) :: scotch_strat(SCOTCH_STRATDIM)

      ! Strings
      Character(len=MAPHYS_STRL) :: rname = "PART_scotch_partgraph"
      Character(len=STRAT_LEN) :: string, strat
      !- End of header ---------------------------------------------------------

      !-------------------------------------------------------------------------
      ! [1.0] Init
      !-------------------------------------------------------------------------

      ! check algorithm
      istep = 1
      iinfo = 0
      If (metalgo /= 4) iinfo = -1
      If (iinfo /= 0 ) Goto 9999
      
      ! scalars
      vertnbr = graph_global_A%ndof
      maphysMaxLevel = INT(Log(REAL(nbdom))/Log(2.0)) ! = log2(nbdom)

      ! arrays
      Nullify(permtab)
      Nullify(peritab)
      Nullify(rangtab)
      Nullify(treetab)
      Nullify(metperm)
      Nullify(metiperm)

      istep = 11
      Allocate(permtab(vertnbr  ), stat=iinfo)
      Allocate(peritab(vertnbr  ), stat=iinfo)
      Allocate(rangtab(vertnbr+1), stat=iinfo)
      Allocate(treetab(vertnbr  ), stat=iinfo)
      If (iinfo /= 0 ) iinfo = -iinfo
      If (iinfo /= 0 ) Goto 9999

      ! strings

      ! form "strat" from strattemplate, by
      ! substituing "MAPHYSMAXLEVEL" by its value
      ! i= starting index of "MAPHYSMAXLEVEL"
      ! j= ending index of "MAPHYSMAXLEVEL"
      ! string = contains the value of MAPHYSMAXLEVEL (in string format)
      string=""
      Write(string,"(I5)") maphysMaxLevel
      string=ADJUSTL(string)
      
      ! Replace all occurences of MAPHYSMAXLEVEL with its value.
      i=0
      strat= strattemplate
      Do While ( i /= -1 )
         ! find "MAPHYSMAXLEVEL"
         i=INDEX(strat,"MAPHYSMAXLEVEL") - 1
         ! Stop if no occurence were found
         If ( i == -1 ) Cycle 
         ! Compute the ending index
         j=i+1+Len_Trim("MAPHYSMAXLEVEL")
         ! Replace
         strat=strat(1:i)//Trim(string)//strat(j:)
      End Do
      
      !-------------------------------------------------------------------------
      ! [2.0] Init SCOTCH structures
      !-------------------------------------------------------------------------
      
      ! convert graph
      istep = 21
      Call MAPHYS2SCOTCH_Graph(graph_global_A, scotch_graph, iinfo)
      If (iinfo < 0 ) Goto 9999

      ! init the strategy
      istep = 22
      Call SCOTCHFStratInit(scotch_strat(1),iinfo)
      If (iinfo /= 0 ) Goto 9999

      istep = 23
      Call SCOTCHFStratGraphOrder(scotch_strat(1), Trim(strat), iinfo)
      If (iinfo /= 0 ) iinfo = -iinfo
      If (iinfo /= 0 ) Goto 9999

      !-------------------------------------------------------------------------
      ! [3.0] Partition the SCOTCH GRAPH
      !-------------------------------------------------------------------------

      ! call scotch
      istep = 31
      Call SCOTCHFGraphOrder(                     &
           scotch_graph%data(1), scotch_strat(1), &
           permtab(1), peritab(1),                &
           cblknbr, rangtab(1), treetab(1),       &
           iinfo )
      If (iinfo /= 0 ) iinfo = -iinfo
      If (iinfo /= 0 ) Goto 9999
      
      ! verify
      istep = 32
      If ( cblknbr /= (2*nbdom - 1)) iinfo = -1
      If (iinfo < 0 ) Goto 9999

      !-------------------------------------------------------------------------
      ! [4.0] Convert SCOTCH data into MAPHYS ones
      !-------------------------------------------------------------------------

      ! Convert SCOTCH treetab/rangtab into MaPHyS ones
      istep = 41
      Call SCOTCH2MAPHYS_part_Tree(nbdom,cblknbr, rangtab, treetab, azz, iinfo )
      If (iinfo < 0 ) Goto 9999

      ! permutations 
      ! - METIS and SCOTCH permuations are inversed 
      ! - we perform copy for 32/64 bits convertions 
      istep = 42
      Allocate(metperm(vertnbr), stat= iinfo)
      Allocate(metiperm(vertnbr), stat= iinfo)
      If (iinfo /= 0 ) iinfo = -iinfo
      If (iinfo /= 0) Goto 9999

      Do i=1, vertnbr
         metperm (i) = INT(peritab(i), kind=MAPHYS_NUMKIND)
         metiperm(i) = INT(permtab(i), kind=MAPHYS_NUMKIND)
      End Do

      !-------------------------------------------------------------------------
      ! [5.0] End
      !-------------------------------------------------------------------------

9999  Continue
      ! print error/warning messages
      If ( iinfo /=  0 ) Then
         
         If ( iinfo > 0) msg_class = MSG_WARNING
         If ( iinfo < 0) msg_class = MSG_ERROR
         
         Select Case(istep) 
         Case(11); Call mph_logWithInfo (msg_class,vertnbr+1,Trim(rname)//&
              "Failed to allocate a few arrays, max size =")
         Case(21); Call mph_logWithInfo (msg_class,iinfo,Trim(rname)//&
              "while converting maphys graph into scotch graph")
         Case(22); Call mph_logWithInfo (msg_class,iinfo,Trim(rname)//&
              "while initializing scotch strategy")
         Case(31); Call mph_logWithInfo (msg_class,iinfo,Trim(rname)//&
              "while calling scotch")
         Case(34); Call mph_logWithInfo (msg_class,cblknbr,Trim(rname)//&
              "while asserting cblknbr == 2*nbdom -1, cblknbr =")
         Case(41); Call mph_logWithInfo (msg_class,iinfo,Trim(rname)//&
              "while converting scotch partitioning tree into maphys one")
         Case(42); Call mph_logWithInfo (msg_class,vertnbr,Trim(rname)//&
              "while allocating permutation vectors, size=")

         End Select
         
      End If
      
      ! report success, failure or warning
      If ( iinfo == 0 ) info =  0
      If ( iinfo <  0 ) info = - istep
      If ( iinfo >  0 ) info = + istep

      ! Free Memory
      If (info == 0 ) Call SCOTCHfGraphExit(scotch_graph%data, iinfo)
      If (info == 0 ) Call SCOTCHFStratExit(scotch_strat,iinfo)

      If( Associated(permtab)) DeAllocate(permtab)
      If( Associated(peritab)) DeAllocate(peritab)
      If( Associated(rangtab)) DeAllocate(rangtab)
      If( Associated(treetab)) DeAllocate(treetab)

      If( Associated(scotch_graph%vendtab)) Nullify   (scotch_graph%vendtab)
      If( Associated(scotch_graph%verttab)) DeAllocate(scotch_graph%verttab)
      If( Associated(scotch_graph%edgetab)) DeAllocate(scotch_graph%edgetab)
      If( Associated(scotch_graph%velotab)) DeAllocate(scotch_graph%velotab)
#endif

    End Subroutine PART_scotch_partgraph


    ! [+] routine : Maphys2scotch_graph ----------------------------------------
    !
    !> Convert a maphys matrix graph into a scotch graph
    !!
    !!----
    !!
    !! @param [in ] maphys_graph   the maphys matrix graph
    !! @param [out] scotch_graph   the scotch matrix graph
    !! @param [out] info           the status
    !!
    !!----
    !!
    !! @author Yohan Lee-tin-yien
    Subroutine Maphys2scotch_graph(maphys_graph, scotch_graph,info)

      !* Modules *!
      Use MPH_part_type, Only : &
           maphys_matrix_graph_t ! struct
      Use mph_log_mod
      Implicit None

      !* Arguments *!
      Type(maphys_matrix_graph_t) , Intent(in)  :: maphys_graph
      Type(maphys_scotch_graph_t) , Intent(out) :: scotch_graph
      Integer                     , Intent(out) :: info

#ifndef HAVE_LIBSCOTCH
      Call MPH_Log(MSG_ERROR,"MaPHyS was not compiled with SCOTCH (-DHAVE_LIBSCOTCH)")
      info = -1
      Return
#else

      !* Local variables *!

      ! scalars
      Integer :: istep, msg_class
      Logical :: GRAPH_VERTEX_is_weighted
      Integer(kind=SCOTCH_SYSKIND) :: iinfo
      Integer(kind=SCOTCH_NUMKIND) :: baseval, vertnbr, edgenbr
      Integer(kind=SCOTCH_NUMKIND) :: i

      ! string
      Character(len=MAPHYS_STRL) :: rname = "Maphys2scotch_graph"

      !- End of header ---------------------------------------------------------

      !-------------------------------------------------------------------------
      ! [1.0] Init
      !-------------------------------------------------------------------------

      ! warning : 
      ! - we need to copy data for 32/64 bits integer conversions
      ! - vendtab is not defined (see SCOTCH documentation about optional arg.)
      ! - edlotab is not defined
      ! - vlbltab is not defined
      ! - velotab is only defined when the vertices are weighted
      

      ! scalars
      GRAPH_VERTEX_is_weighted = ( maphys_graph%algo == 3 )
      scotch_graph%baseval =  1
      scotch_graph%vertnbr =  maphys_graph%ndof
      scotch_graph%edgenbr =  maphys_graph%xadj(maphys_graph%ndof+1)-1
      scotch_graph%data(:) =  0.d0 ! unsignificant init

      baseval = scotch_graph%baseval
      vertnbr = scotch_graph%vertnbr
      edgenbr = scotch_graph%edgenbr

      ! arrays
      istep = 11

      Nullify(scotch_graph%verttab)
      Nullify(scotch_graph%vendtab)
      Nullify(scotch_graph%edgetab)
      Nullify(scotch_graph%velotab)

      Allocate(scotch_graph%verttab(vertnbr+1), stat=iinfo)
      Allocate(scotch_graph%edgetab(edgenbr), stat=iinfo)

      If (GRAPH_VERTEX_is_weighted) &
           Allocate(scotch_graph%velotab(vertnbr), stat=iinfo)
      If(iinfo /= 0) iinfo = -iinfo
      If(iinfo /= 0) Goto 9999

      Do i=1, vertnbr+1
         scotch_graph%verttab(i) = INT(maphys_graph%xadj(i),SCOTCH_NUMKIND)
      End Do

      Do i=1, edgenbr
         scotch_graph%edgetab(i) = INT(maphys_graph%adjncy(i),SCOTCH_NUMKIND)
      End Do

      If (GRAPH_VERTEX_is_weighted) Then
         Do i=1,vertnbr
            scotch_graph%velotab(i) = INT(maphys_graph%vwgt(i),SCOTCH_NUMKIND)
         End Do
      End If

      scotch_graph%vendtab => scotch_graph%verttab(2:vertnbr+1)

      ! structures
      istep = 12
      Call SCOTCHfGraphInit(scotch_graph%data(1), iinfo)
      If(iinfo /= 0) iinfo = -iinfo
      If(iinfo /= 0) Goto 9999
      
      !-------------------------------------------------------------------------
      ! [2.0] Build
      !-------------------------------------------------------------------------

      If (GRAPH_VERTEX_is_weighted) Then
         istep = 21
         Call SCOTCHfGraphBuild(scotch_graph%data(1), &
              scotch_graph%baseval    , scotch_graph%vertnbr, &
              scotch_graph%verttab(1) , scotch_graph%verttab(2), &
              scotch_graph%velotab(1), scotch_graph%verttab(1) , &
              scotch_graph%edgenbr, scotch_graph%edgetab(1) , &
              scotch_graph%edgetab(1), iinfo)
      Else
         istep = 22
         Call SCOTCHfGraphBuild(scotch_graph%data(1), &
              scotch_graph%baseval    , scotch_graph%vertnbr, &
              scotch_graph%verttab(1) , scotch_graph%verttab(2), &
              scotch_graph%verttab(1), scotch_graph%verttab(1) , &
              scotch_graph%edgenbr, scotch_graph%edgetab(1) , &
              scotch_graph%edgetab(1), iinfo)
      End If
      If (iinfo /= 0) iinfo = -iinfo
      If (iinfo /= 0) Goto 9999

      !-------------------------------------------------------------------------
      ! [3.0] Check
      !-------------------------------------------------------------------------

      istep = 30
      Call mph_logWithInfo(MSG_DEBUG,vertnbr,"vertnbr =")
      Call mph_logWithInfo(MSG_DEBUG,edgenbr,"edgenbr =")
      Call SCOTCHfGraphCheck(scotch_graph%data(1), iinfo )
      If (iinfo /= 0) iinfo = -iinfo
      If (iinfo /= 0) Goto 9999

      !-------------------------------------------------------------------------
      ! [4.0] End
      !-------------------------------------------------------------------------

9999  Continue

      ! print error/warning messages
      If ( iinfo /=  0 ) Then
         
         If ( iinfo > 0) msg_class = MSG_WARNING
         If ( iinfo < 0) msg_class = MSG_ERROR
         
         Select Case(istep) 
         Case(11); Call mph_logWithInfo (msg_class,&
              Max(vertnbr+1,edgenbr),Trim(rname)//&
              "Failed to allocate a few arrays, max size =")
         Case(12); Call mph_logWithInfo (msg_class,iinfo,Trim(rname)//&
              "while initializing into scotch graph")
         Case(21); Call mph_logWithInfo (msg_class,iinfo,Trim(rname)//&
              "while building scotch graph (weighted vertices)")
         Case(22); Call mph_logWithInfo (msg_class,iinfo,Trim(rname)//&
              "while building scotch graph (no weighted vertices)")
         Case(30); Call mph_logWithInfo (msg_class,iinfo,Trim(rname)//&
              "while checking scotch graph")
         End Select
         
      End If
      
      ! report success, failure or warning
      If ( iinfo == 0 ) info =  0
      If ( iinfo <  0 ) info = - istep
      If ( iinfo >  0 ) info = + istep

      ! free memory on error
      If ( iinfo < 0 )Then
         ! msg
         write(6,*) trim(rname), " : Internal Error at step = ", istep
         info = -istep

         Call SCOTCHfGraphExit(scotch_graph%data, iinfo)
         If( Associated(scotch_graph%vendtab)) Nullify   (scotch_graph%vendtab)
         If( Associated(scotch_graph%verttab)) DeAllocate(scotch_graph%verttab)
         If( Associated(scotch_graph%edgetab)) DeAllocate(scotch_graph%edgetab)
         If( Associated(scotch_graph%velotab)) DeAllocate(scotch_graph%velotab)
      Endif
#endif

    End Subroutine Maphys2scotch_graph

    ! [+] routine : Scotch2maphys_part_tree -----------------------------------
    !
    !> Convert a scotch partitioning tree into a maphys one
    !!
    !!----
    !!
    !! @param [in ] scotch_tree    the scotch partitioning tree
    !! @param [out] maphys_tree    the maphys partitioning tree
    !! @param [out] info           the status
    !!
    !!----
    !!
    !! @author Yohan Lee-tin-yien
    Subroutine scotch2maphys_part_tree &
         (nbdom, cblknbr, rangtab, treetab, azz, info )

      !* Modules *!
      Use MPH_part_type, Only : &
           maphys_binary_tree_t ! struct
      Use mph_log_mod
      Implicit None

      !* Arguments *!
      Integer                     , Intent(in) :: nbdom
      Integer(kind=SCOTCH_NUMKIND), Intent(in) :: cblknbr
      Integer(kind=SCOTCH_NUMKIND), Intent(in), Pointer :: rangtab(:)
      Integer(kind=SCOTCH_NUMKIND), Intent(in), Pointer :: treetab(:)
      Type(maphys_binary_node_t) , Intent(out), Pointer :: azz(:)
      Integer                    , Intent(out) :: info

#ifndef HAVE_LIBSCOTCH
      Call MPH_Log(MSG_ERROR,"MaPHyS was not compiled with SCOTCH (-DHAVE_LIBSCOTCH)")
      info = -1
      Return
#else

      !* Local variables *!
      ! constants
      Integer, Parameter :: NOT_COMPUTED_HERE = -1
      Integer, Parameter :: is_LEAF = -1
      
      ! scalars
      Integer :: msg_class
      Integer(kind=SCOTCH_SYSKIND) :: istep
      Integer(kind=SCOTCH_SYSKIND) :: iinfo
      Integer(kind=SCOTCH_NUMKIND) :: i
      Integer(kind=SCOTCH_NUMKIND) :: rblk, lblk
      Integer(kind=SCOTCH_NUMKIND) :: rblk2, lblk2
      Integer(kind=SCOTCH_NUMKIND) :: sep

      !> convert a SCOTCH block index
      !! into a MAPHYS separator index.
      !!
      !! @warning the SCOTCH block index of root (= -1) is omitted.
      Integer(kind=SCOTCH_NUMKIND) :: blk2sepidx(cblknbr) 

      ! Strings
      Character(len=MAPHYS_STRL) :: rname = "SCOTCH2MAPHYS_PART_TREE"

      !- End of header ---------------------------------------------------------

      !-------------------------------------------------------------------------
      ! [1.0] Initialize local data
      !-------------------------------------------------------------------------

      !---------------------------------------------------------------
      ! [1.1] allocate memory
      !---------------------------------------------------------------
      istep =11
      iinfo = 0
      Nullify(azz)
      Allocate(azz(nbdom-1), stat= iinfo )
      If (iinfo /= 0) iinfo = -iinfo
      If (iinfo /= 0 ) Goto 9999

      !---------------------------------------------------------------
      ! Construct the index converters  (blk2sepidx, blk2leafidx)
      !---------------------------------------------------------------
      !
      ! - If the block "iblk" is not a separator, blk2sepidx(iblk) = -1
      ! - and start from 1.
      ! - We do not store block = -1 (root) 
      istep=12

      ! init
      sep  = 0
      blk2sepidx (:) = -1

      ! fill
      Do i = cblknbr, 1, -1
         ! get the left and right sons blocks 
         istep = 121
         Call get_sons(i, cblknbr, rangtab, treetab, lblk, rblk, iinfo )
         If (iinfo < 0 ) Goto 9999

         ! Check sons
         istep = 122
         If ((lblk  == is_LEAF) .and. (rblk  /= is_LEAF)) iinfo = -1 
         If ((lblk  /= is_LEAF) .and. (rblk  == is_LEAF)) iinfo = -1 
         If (iinfo < 0 ) Goto 9999

         ! handle separators
         If ((lblk  /= is_LEAF) .and. (rblk  /= is_LEAF)) Then
            sep = sep+1
            blk2sepidx(i) = sep
         End If

      End Do

      ! verify
      istep = 123
      If ( sep  /= nbdom - 1) iinfo = -1
      If (iinfo < 0 ) Goto 9999

      !-------------------------------------------------------------------------
      ! [2.0] Init the later computed fields
      !-------------------------------------------------------------------------

      azz(:)%level  = NOT_COMPUTED_HERE
      azz(:)%rdomid = NOT_COMPUTED_HERE
      azz(:)%ldomid = NOT_COMPUTED_HERE

      !-------------------------------------------------------------------------
      ! [3.0] Fill the other fields
      !-------------------------------------------------------------------------
      istep = 30

      ! init counters
      sep = 0
      ! fill the structures
      Do i = cblknbr, 1, -1

         ! get the left and right sons 
         istep = 31
         Call get_sons(i, cblknbr, rangtab, treetab, lblk, rblk, iinfo )
         If (iinfo < 0 ) Goto 9999

         ! Check sons
         istep = 32
         If ((lblk  == is_LEAF) .and. (rblk  /= is_LEAF)) iinfo = -1 
         If ((lblk  /= is_LEAF) .and. (rblk  == is_LEAF)) iinfo = -1 
         If (iinfo < 0 ) Goto 9999

         ! jump the block "i" is a leaf
         If ((lblk  == is_LEAF) .and. (rblk  == is_LEAF)) Cycle 

         ! increment the separator counter
         sep = sep +1

         ! save the father index
         If ( treetab(i) == -1 ) azz(sep)%father = 0 ! root
         If ( treetab(i) /= -1 ) azz(sep)%father = blk2sepidx( treetab(i) )

         ! save the sons indexes
         azz(sep)%rson  = min(blk2sepidx(rblk),blk2sepidx(lblk))         
         azz(sep)%lson  = max(blk2sepidx(rblk),blk2sepidx(lblk))
         
         ! get the range of each domains
         azz(sep)%sepst = rangtab(i)
         azz(sep)%seped = rangtab(i+1) -1

         ! fill lgst, lged if the sons of lson are leafs
         Call get_sons(lblk, cblknbr, rangtab, treetab, lblk2, rblk2, iinfo )
         If (iinfo /= 0 ) iinfo = 0 ! ignore errors
         If ((lblk2  == is_LEAF) .and. (rblk2  == is_LEAF)) Then
            azz(sep)%lgst = rangtab(lblk)
            azz(sep)%lged = rangtab(lblk+1) -1
         Else
            azz(sep)%lgst = 0
            azz(sep)%lged = 0
         End If

         ! fill rgst, rged if the sons of rson are leafs
         Call get_sons(rblk, cblknbr, rangtab, treetab, lblk2, rblk2, iinfo )
         If (iinfo /= 0 ) iinfo = 0 ! ignore errors
         If ((lblk2  == is_LEAF) .and. (rblk2  == is_LEAF)) Then
            azz(sep)%rgst = rangtab(rblk)
            azz(sep)%rged = rangtab(rblk+1) -1
         Else
            azz(sep)%rgst = 0
            azz(sep)%rged = 0
         End If

      End Do

      ! verify
      istep = 33
      If (sep /= nbdom - 1) iinfo = -1
      If (iinfo < 0 ) Goto 9999
      
      !-------------------------------------------------------------------------
      ! [4.0] Return
      !-------------------------------------------------------------------------

9999 Continue
      
      ! print error/warning messages
      If ( iinfo /=  0 ) Then
         
         If ( iinfo > 0) msg_class = MSG_WARNING
         If ( iinfo < 0) msg_class = MSG_ERROR
         
         Select Case(istep) 
         Case(11); Call mph_logWithInfo (msg_class, nbdom-1,Trim(rname)//&
              "Failed to allocate a structure, nb elements =")
         Case(121); Call mph_logWithInfo (msg_class,iinfo,Trim(rname)//&
              "while getting the right and left son blocs")
         Case(122); Call mph_logWithInfo (msg_class,iinfo,Trim(rname)//&
              "while checking the sons blocs")
         Case(123); Call mph_logWithInfo (msg_class,sep,Trim(rname)//&
              "while asserting nbdom-1 == sep, sep= ")
         Case(31); Call mph_logWithInfo (msg_class,iinfo,Trim(rname)//&
              "while getting the right and left son blocs")
         Case(32); Call mph_logWithInfo (msg_class,iinfo,Trim(rname)//&
              "while checking the sons blocs")
         Case(33); Call mph_logWithInfo (msg_class,sep,Trim(rname)//&
              "while asserting nbdom-1 == sep, sep= ")
         End Select
         
      End If
      
      ! report success, failure or warning
      If ( iinfo == 0 ) info =  0
      If ( iinfo <  0 ) info = - istep
      If ( iinfo >  0 ) info = + istep

#endif

    End Subroutine scotch2maphys_part_tree


    ! [+] routine : get_sons ---------------------------------------------------
    !
    !> Get the left and right sons ('LBLK','RBLK')
    !! of the block "father".
    !!
    !! LBLK is the "first" domain which have "DOM" has father
    !! RBLK is the "last" domain which have "DOM" has father
    !! If they are leafs, return LBLK = RBLK = -1
    !! This routine do not print any outputs.
    !!
    !!----
    !!
    !! @param [in] father   specifies which block 
    !! @param [in] cblcknbr specifies the number of blocks
    !! @param [in] rangtab  specifies the range of each blocks
    !! @param [in] treetab  specifies the father of each block
    !! @param [out] rblk    specifies the right son block index
    !! @param [out] lblk    specifies the left  son block index
    !! @param [out] info    specifies the routine status
    !!
    !!-----
    !!
    !! @author Yohan Lee-tin-yien
    Subroutine get_sons( father, &
         cblknbr, rangtab, treetab, &
         rblk, lblk, info )
      Implicit None
      
      !* Arguments *!
      Integer(kind=SCOTCH_NUMKIND), Intent(in) :: father
      Integer(kind=SCOTCH_NUMKIND), Intent(in) :: cblknbr
      Integer(kind=SCOTCH_NUMKIND), Intent(in), Pointer :: rangtab(:)
      Integer(kind=SCOTCH_NUMKIND), Intent(in), Pointer :: treetab(:)
      Integer(kind=SCOTCH_NUMKIND), Intent(out) :: rblk
      Integer(kind=SCOTCH_NUMKIND), Intent(out) :: lblk
      Integer                     , Intent(out) :: info

#ifndef HAVE_LIBSCOTCH
      info = -1
      Return
#else

      !* Local Variables *!
      Integer :: istep
      Integer(kind=SCOTCH_NUMKIND) :: i

      !- End of header ---------------------------------------------------------
      !-------------------------------------------------------------------------
      ! [1.0] Init 
      !-------------------------------------------------------------------------

      info = 0

      !yohan>
      ! Trick to avoid the compilation warning "rangtab unused".
      ! As "treetab","rangtab" are closely linked,
      ! I don't want to remove it from the list of arguments.
      ! Hence I trick the compiler with the following line :
      !yohan<
      i = rangtab(1)

      !-------------------------------------------------------------------------
      ! [2.0] check arguments
      !-------------------------------------------------------------------------
      istep = 2
      If ( father < 1       ) info = -1
      If ( father > cblknbr ) info = -1
      If (info /= 0) Goto 9999

      !-------------------------------------------------------------------------
      ! [3.0] Compute lblk
      !-------------------------------------------------------------------------

      lblk = -1
      Do i=1, cblknbr
         If ( treetab(i) == father ) Then
            lblk = i
            Exit
         End If
      End Do

      !-------------------------------------------------------------------------
      ! [4.0] Compute rblk
      !-------------------------------------------------------------------------

      rblk = -1
      Do i=cblknbr, 1, -1
         If ( treetab(i) == father ) Then
            rblk = i
            Exit
         End If
      End Do

      !-------------------------------------------------------------------------
      ! [5.0] if found, check that they are different
      !-------------------------------------------------------------------------
      istep = 5
      If ( (rblk == lblk) .and. (lblk /= -1)) info = -1 
      If (info /= 0) Goto 9999

      !-------------------------------------------------------------------------
      ! [6.0] Return
      !-------------------------------------------------------------------------

      ! Handle errors
      ! On errors print nothing as
      ! caller may ignore errors.
9999 Continue
      If (info /= 0 ) info = - istep
#endif
    End Subroutine get_sons


  End Module MPH_part_scotch_mod
