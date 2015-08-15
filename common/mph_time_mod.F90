  ! handle timing 
  Module MPH_time_mod
    Implicit None
    
    Public :: MPH_time_start
    Public :: MPH_time_stop
    Public :: MPH_time_add
    Public :: MPH_time_reset

    Contains

      Subroutine MPH_time_start(time)
        Implicit None
        Real(Kind=8), Intent(out) :: time
        Real(Kind=8), External    :: MPI_Wtime
        time = MPI_Wtime()
      End Subroutine MPH_time_start

      Subroutine MPH_time_stop(time)
        Implicit None
        Real(Kind=8), Intent(inout) :: time
        Real(Kind=8), External    :: MPI_Wtime
        time = MPI_Wtime() - time
      End Subroutine MPH_time_stop

      Subroutine MPH_time_add(time,timetoadd)
        Implicit None
        Real(Kind=8), Intent(inout) :: time
        Real(Kind=8), Intent(in   ) :: timetoadd
        time = time + timetoadd
      End Subroutine MPH_time_add

      Subroutine MPH_time_reset(time)
        Implicit None
        Real(Kind=8), Intent(out) :: time
        time = 0.d0
      End Subroutine MPH_time_reset

  End Module MPH_time_mod
