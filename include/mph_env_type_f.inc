  Type mph_env_t
     sequence
     Integer :: comm            ! Global MPI communicator
     Integer :: nodecomm        ! MPI communicator on the node.
     Integer :: nodeid       ! color of the node.
  End type mph_env_t
