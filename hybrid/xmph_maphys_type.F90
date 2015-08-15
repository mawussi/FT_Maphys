! Warning: XMPH_GENFILE_COMMENT
#include "mph_defs_f.h"

!> Module defining the derived type "XMPH_maphys_t"
  Module XMPH_maphys_type
    
    !* Modules *!
    Use XMPH_sparse_matrix_type
    Use XMPH_sls_mod
    Use XMPH_dls_mod
    Use MPH_part_type, Only : &
         maphys_rhs_partition_t, &
         maphys_domain_t
    Use MPH_mem_mod, Only : &
         MPH_mem_t

    !* No implicit typing *!
    Implicit none
    Include "mph_env_type_f.inc"

    !* Type defintions *!
    Type XMPH_maphys_t; Sequence
     ! ===============
     ! USER PARAMETERS
     ! ===============
     !   
     !> MPI communicator
     !!   
     !! The MPI Communicator used by maphys.\n
     !! It must be set with a valid MPI Communicator before any call
     !! to maphys.
     !!    
     Integer :: comm 

     ! Input matrix (in coordinate format)
     ! ---------------------------------------
     
     !> matrix symmetry
     !! - 0 = unsymmetric
     !! - 1 = symmetric
     !! - 2 = SPD
     !! .
     Integer                :: sym 
     !> order of the matrix
     Integer                :: n
     !> number of entries in the sparse matrix.
     !! If the matrix is symmetric or SPD,
     !! user only gives half of the entries.
     Integer                :: nnz
     Integer      , pointer :: rows (:)
     Integer      , pointer :: cols (:)
     XMPH_FLOAT , pointer :: values    (:)

     !> ask to write the input matrix to a file
     Character(len = MAPHYS_STRL) :: write_matrix

     ! Right-hand-side (in dense format ordered by columns)
     ! ----------------------------------------------------
     Integer  :: nrhs  
     XMPH_FLOAT , pointer :: rhs (:) 

     ! Solution (in dense format ordered by columns)
     ! ----------------------------------------------------
     XMPH_FLOAT , pointer :: sol (:) 

     ! Controls
     ! --------
     Integer      :: job 
     Integer(kind=4) :: icntl  ( MAPHYS_ICNTL_SIZE  ) 
     Real   (kind=8) :: rcntl  ( MAPHYS_RCNTL_SIZE  ) 
     
     ! Statistics
     ! -----------
     ! version
     Character(len=MAPHYS_STRL) :: version

     ! on this process (MPI)

     Integer(kind=4) :: iinfo ( MAPHYS_IINFO_SIZE ) 
     Real   (kind=8) :: rinfo ( MAPHYS_RINFO_SIZE ) 
     
     ! on all process (MPI)

     Integer(kind=4) :: iinfomin ( MAPHYS_IINFO_SIZE ) 
     Integer(kind=4) :: iinfomax ( MAPHYS_IINFO_SIZE ) 
     Real   (kind=8) :: iinfoavg ( MAPHYS_IINFO_SIZE ) 
     Real   (kind=8) :: iinfosig ( MAPHYS_IINFO_SIZE ) 

     Real   (kind=8) :: rinfomin ( MAPHYS_RINFO_SIZE ) 
     Real   (kind=8) :: rinfomax ( MAPHYS_RINFO_SIZE ) 
     Real   (kind=8) :: rinfoavg ( MAPHYS_RINFO_SIZE )
     Real   (kind=8) :: rinfosig ( MAPHYS_RINFO_SIZE ) 

     Integer(kind=4) :: iinfog   ( MAPHYS_IINFOG_SIZE ) 
     Real   (kind=8) :: rinfog   ( MAPHYS_RINFOG_SIZE ) 

     ! =====================
     ! Internal working data
     ! =====================
     
     ! internal controls
     ! -----------------
     Integer(kind=4) :: ikeep ( MAPHYS_IKEEP_SIZE )    
     Real   (kind=8) :: rkeep ( MAPHYS_RKEEP_SIZE )    

     ! memory
     Type(MPH_mem_t) :: mem

     ! environement
     Type(mph_env_t) :: env

     ! Description of the Non-overlapping domain Decomposition 
     ! -------------------------------------------------------
     ! Local domain description (interface + interior)
     type(maphys_domain_t)     :: lc_domain 

     ! Necessary data to part the right hand side
     Type(maphys_rhs_partition_t) :: part_rhs
  
     ! Blocs on the matrix of the domain 
     ! [  Aii     Aib ] 
     ! [  Abi     Abb ] 
     ! - i: interior  related (nodes inside the domain)
     ! - b: boundary  related (nodes on its interface )
     Type(XMPH_sparse_matrix_t) :: sm_Aii , sm_Aib   
     Type(XMPH_sparse_matrix_t) :: sm_Abi , sm_Abb   

     ! saved scalings on the local matrix
     XMPH_FLOAT, pointer :: row_scaling (:)  
     XMPH_FLOAT, pointer :: col_scaling (:)  
      
     ! linear system on a domain
     ! Its interior is solved with a sparse direct solver
     ! --------------------------------------------------
     Type(XMPH_sls_t) :: sls_domain  

     ! On the interface : the schur complement linear system
     ! (solved with an iterative method)
     ! -----------------------------------------------------
     ! local schur complements matrices (exact or approximate)
     Type(XMPH_dense_matrix_t)  :: dm_schur ! exact   
     Type(XMPH_sparse_matrix_t) :: sm_schur ! approximate   

     ! preconditioners for the iterative method
     Type(XMPH_dls_t)  :: dls_precond_schur
     Type(XMPH_sls_t)  :: sls_precond_schur

  End Type XMPH_maphys_t


End Module XMPH_maphys_type
