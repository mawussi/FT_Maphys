! fortran definitions for maphys, common to all arithmetics
#ifndef MPH_DEFS_F_H__
#define MPH_DEFS_F_H__

#ifndef MAPHYS_VERSION
#define MAPHYS_VERSION '0.1.5'
#endif

! Status
#define MPH_SUCCESS         0

! Logical
#define MPH_LOGICAL   INTEGER
#define MPH_TRUE            0
#define MPH_FALSE           1
#define MPH_LOGICAL_UNSET   2

! Integers
#define MPH_INT        Integer
#define MPH_INTKIND         4
#define MPH_INTBYTESIZE     4
#define MPH_INTMPI     MPI_INTEGER


! global size
#define MAPHYS_STRL        1024
#define MAPHYS_ICNTL_SIZE  45
#define MAPHYS_RCNTL_SIZE  40 
#define MAPHYS_IINFO_SIZE  40
#define MAPHYS_RINFO_SIZE  40
#define MAPHYS_IINFOG_SIZE 5
#define MAPHYS_RINFOG_SIZE 20 
#define MAPHYS_IKEEP_SIZE  41
#define MAPHYS_RKEEP_SIZE  40   

#define MTHREAD_ICNTL_SIZE        5
#define ANA_TIMING_SIZE          10

#define SDS_MAX_INDEX      3

#endif 
