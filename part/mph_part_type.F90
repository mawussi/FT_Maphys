#include "mph_defs_f.h"

!> derived types for MAPHYS partitioner.
!!
Module MPH_part_type

  !* Modules *!

  ! C interopability
#if defined(HAVE_ISO_C_BINDING)
  Use, Intrinsic :: ISO_C_BINDING
#endif

  !* No Implicit typing *!
  Implicit None


  !lty> begin unsure comments

  ! [+] Type : maphys_rhs_partition_t ---------------------------------------
  !
  !> data used to split the global system into local system 
    Type maphys_rhs_partition_t
       sequence

       !* global data *! (only on master = 0)

       !> number of domains
       Integer :: nbdom

       !> order of the global matrix 
       Integer :: gballndof          
       
       !> order of the schur complement
       !! @Note  it is a copy of maphys_domain_t%gballintrf.
       Integer :: gballintrf  

       !> ? vertex ? size
       INTEGER :: combivtxpcsz                 

       !> processus displacement for domain's interior 
       INTEGER, pointer :: procintdisp (:)      
       !> processus displacement for domain's interface
       INTEGER, pointer :: procintrfdisp  (:)   
       !> domain langrage  
       INTEGER, pointer :: domLg (:)            
       !> degree of freedom of each domain's interior 
       INTEGER, pointer :: domintdof (:)        

       !> degree of freedom of each domain's interface
       INTEGER, pointer :: domintrfdof (:)     
       
       !> permutation done on the global matrix
       INTEGER, pointer :: metperm (:)          

       !* Strategy *!
       INTEGER :: rhsway ! selected strategy (1..2) , forced to 2
       ! 1 rst strategy 
       INTEGER, pointer :: scatindices (:)      ! allow to scatter data ?
       ! 2 nd strategy   
       INTEGER, pointer :: scatlogicindices (:) ! allow to scatter data ?
       INTEGER, pointer :: domlogicintrfdof(:)  ! degree of freedom 

       !> global to local numbering in the interface (size = gballintrf)
       !! Note: (Yohan Lee-tin-yien 2011-10-03)
       !! gbtoloc is only useful with a global preconditioner.
       !! It is unused in the current version,
       !! and not exported in the read/write routine,
       !! but are left there for future use.
       Integer, pointer :: gbtoloc      (:)     

       
    end Type maphys_rhs_partition_t

    ! [+] Type : maphys_domain_t -----------------------------------------------
    !
    !> defines the data on each domain
    Type maphys_domain_t
       sequence
       !------------------------------------------
       !> @par [ local matrix description ]
       !------------------------------------------

       !> number of rows (columns) in the local matrix
       Integer :: myndof           

       !> number of rows (columns) in the local matrix' interior  
       Integer :: myndofinterior   

       !> number of rows (columns) in the local matrix' interface
       Integer :: myndofintrf      

       !------------------------------------------
       !> @par [ local matrix treatment ]
       !------------------------------------------

       !> lagrange related
       Integer :: myint_Lg       

       ! -------------------------------------------
       !> @par [ Interface Description ]
       ! -------------------------------------------

       !> size of array myindexintrf
       Integer :: lenindintrf                   

       !> size of the interface
       Integer :: mysizeIntrf  

       !> sum on all processes of "mysizeIntrf"
       !! @Note  it is here to avoid a synchrone communication 
       !!         during the intialization of the iterative solver.
       !! @Note  it is a copy of maphys_part_rhs_t%gballintrf.
       Integer :: gballintrf  

       ! -------------------------------------------
       !> @par [ Neigbor Description ]
       ! -------------------------------------------

       !> number of neighbors 
       !!
       !! number of subdomains neighbor of the subdomain.
       !! no (me+1) as me is the MPI process number starting
       !! at 0.
       Integer :: mynbvi            

       !> index of the neighbors (size mynbvi)
       !!
       !! no of the neighboring subdomain 
       Integer, pointer :: myindexVi    (:)     

       !> pointer on the interface points (size mynbvi ?? +1 ?? )
       !!
       !! pointer on the index of the first point on the interface
       !! with the sudomain. The indices of the interface points are
       !! stored in the array IndexIntrf.
       Integer, pointer :: myptrindexVi (:)

       !> list of the interface points (size lenindintrf)
       !!
       !! list of the interface points, pointed by the
       !! array ptr_Index_Intrfc.
       !! The indices are given in the local ordering.
       Integer, pointer :: myindexintrf (:)

       !> array of size mysizeIntrf
       !!
       !! Convert local indexing to the global indexing
       !!
       Integer, pointer :: myinterface  (:)
       
       ! -------------------------------------------
       !> @par [ Interface rows/columns responsability ]
       ! -------------------------------------------
       !
       ! Each processor is responsible of
       ! several rows on the interface.
       ! This is namely used to know how to gather the solution. 
       !

       !> number of rows/columns in the interface
       !! that this domain is responsible
       Integer :: myndoflogicintrf
 
       !> list of rows/columns in the interface
       !! that this domain is responsible. (in local indexing)
       Integer, pointer :: mylogicintrf (:)

       !> The weight/contribution of each node on the interface.
       !!
       !! Array of size mysizeIntrf.
       !! value must be between 0.d0 and 1.d0.
       !!
       !! Currently, for each node on the interface with global index "i". 
       !! - weight(gbtoloc(i)) = 1.d0,
       !!   if the processor is the one responsible of the node.
       !! - weight(gbtoloc(i)) = 0.d0, 
       !!   if not.
       !!
       Real(kind=8), Pointer :: weight (:)

    end Type maphys_domain_t

    !< end unsure comments

    ! [+] Type : maphys_matrix_graph_t -----------------------------------------
    !
    !> the corresponding graph to a sparse square matrix
    !! matrix graph (optionally weighted) is represented by an adjacency list.
    !!
    !! @see metis manual section 5.1
    Type maphys_matrix_graph_t
       sequence
       !> order of the sparse matrix
       Integer          :: ndof       
       !> number of entries of sparse matrix (number of non-zeros)
       MPH_INT       :: nnz 

       !> see metis manual section 5.1 (size ndof+1)
       Integer, Pointer :: xadj   (:) 
       !> see metis manual section 5.1 (size nnz   )
       Integer, Pointer :: adjncy (:) 
       !> specifies the weight on vertices (size ndof)
       Integer, Pointer :: vwgt   (:) 

       !> algorithm used to split the graph
       Integer          :: algo       

       !* statistics *!

       !> maximum number of adjacent vertices in the graph
       Integer :: maxadj    
    end Type maphys_matrix_graph_t

    ! [+] Type : maphys_binary_node_t -------------------------------------------
    !
    !> a node in the binary tree
    !!
#ifndef HAVE_ISO_C_BINDING
      Type maphys_binary_node_t
         sequence

         !-------------------------------
         !> @par [ binary tree description ]
         !-------------------------------

         !> the parent node      (1..nbsep)
         !> if root, value = 0
         Integer  :: father 

         !> the right child node (1..nbsep)
         !> if no child, value = -1
         Integer  :: rson

         !> the left  child node (1..nbsep)
         !> if no child, value = -1
         Integer  :: lson   

         !-------------------------------
         !> @par [ ordering description ]
         !-------------------------------

         !> row (column) index at which the separator starts (1..ndof)
         Integer  :: sepst  
         !> row (column) index at which the separator ends   (1..ndof) 
         Integer  :: seped  
         !> row (column) index at which the right interior domain starts 
         Integer  :: rgst   
         !> row (column) index at which the right interior domain ends
         Integer  :: rged   
         !> row (column) index at which the left  interior domain starts
         Integer  :: lgst   
         !> row (column) index at which the left  interior domain ends
         !>  if no interior domain, value = 0
         Integer  :: lged   

         !-------------------------------
         !> @par [ appended attributes ]
         !-------------------------------

         !> depth of the node    
         Integer  :: level   
         !> number of the right interior domain (1..nbdomain)
         Integer  :: rdomid  
         !> number of the  left interior domain (1..nbdomain)
         !>  if no interior domain, value = -1
         Integer  :: ldomid  
         
      end Type maphys_binary_node_t

#else
      Type, bind(C) :: maphys_binary_node_t
         Integer(C_INT)  :: father
         Integer(C_INT)  :: rson
         Integer(C_INT)  :: lson
         
         Integer(C_INT)  :: sepst
         Integer(C_INT)  :: seped

         Integer(C_INT)  :: rgst
         Integer(C_INT)  :: rged
         Integer(C_INT)  :: lgst
         Integer(C_INT)  :: lged

         Integer(C_INT)  :: level  

         Integer(C_INT)  :: rdomid 
         Integer(C_INT)  :: ldomid 

      end Type maphys_binary_node_t

#endif

      ! [+] Type : maphys_binary_tree_t ----------------------------------------
      !
      !> describe maphys partition binary tree
      Type maphys_binary_tree_t
#ifndef HAVE_ISO_C_BINDING
         sequence
#endif
         !----------------------------
         !> @par [ the tree elements ]
         !----------------------------

         !> size of the binary tree (number of separator)
         Integer :: nbsep    
         !> the elements of the binary tree 
         Type(maphys_binary_node_t), Pointer :: node (:) 

         !----------------------------
         !> @par [ related data ]
         !----------------------------

         !> level of the root of the the tree
         Integer :: toplevel            
         !> number of vertexes in the separators
         Integer :: totinterface         
         !> pointer to domain's separator (level=1) (size nbsep + 1) 
         Integer, Pointer :: domptr   (:) 

         Integer, Pointer :: domstptr (:) ! 
         Integer, Pointer :: domintdof(:) ! 

         !----------------------------
         !> @par [ other ]
         !----------------------------
         
         !> permutation applied on the matrix
         INTEGER, Pointer :: metperm (:) 
         !> inverse of the permutation applied on the matrix
         INTEGER, Pointer :: metiperm (:) 
         
       end Type maphys_binary_tree_t

       ! [+] Type : maphys_domains_t ----------------------------------------
       !
       !> data used in the analysis phase to describes
       !! all the domains.
       !!
       Type maphys_domains_t
         sequence

         !> number of domains
         Integer :: nbdom 

         !-------------------------------
         !> @par [ interior description ] 
         !-------------------------------

         !> nnz in interiors (size nbdom)
         MPH_INT, Pointer :: domintnnz  (:) 
         !> dof of each interior (size nbdom)
         Integer , Pointer :: domintdof (:)
         !> Pointer to partitioning tree's leafs
         Integer, Pointer :: domstptr (:)  

         !-------------------------------
         !> @par [ interface description ] 
         !-------------------------------

         !> nnz in interfaces (size nbdom)
         MPH_INT, Pointer :: domintrfnnz(:)
         !> nnz on all interfaces
         MPH_INT :: nnzallinterface     
         !> degree of freedom of each domain's interface
         Integer, Pointer  :: domintrfdof (:)
         !> weight (lty : ? which one ?) on vertex(size totinterface)
         Integer, Pointer :: vtxweight (:) 
         !> lty :? vertex position ? (size totinterface+1)
         Integer, Pointer :: vtxpos    (:) 

         ! To construct each domain, we need to identify 
         ! the interfaces by analyzing the binary tree.
         ! To each node "i" on an interface, we associate 
         ! a unique identifier and its domain id,
         ! which are respectively stored into 
         ! 'intrindices(i)' and 'intrfproc(i)'
         ! (range for i is 1..combivtxpcsz) 
         
         !> lty: ? sum on interfaces of "myintrfndof" ?
         Integer :: combivtxpcsz               
         !> interface's vertex index  (size >= combivtxpcsz) 
         Integer, Pointer  :: intrfindices (:) 
         !> interface's vertex domain (size >= combivtxpcsz) 
         Integer, Pointer  :: intrfproc    (:) 

         !-------------------------------
         !> @par [ statistics ] 
         !-------------------------------

         !> maximum of interfaces' degrees of freedom 
         Integer :: maxprocintrf               
         !> minimum of interfaces' degrees of freedom 
         Integer :: minprocintrf  
         !> see maphys_binary_tree%totinterface             
         Integer :: totinterface  

         !-------------------------------
         !> @par [ other ] 
         !-------------------------------

         !> Lagrange related
         Integer, Pointer :: domLg    (:) 
         !> permutation on rows/columns
         Integer, Pointer :: metperm  (:) 

      End Type maphys_domains_t

    End Module MPH_part_Type
