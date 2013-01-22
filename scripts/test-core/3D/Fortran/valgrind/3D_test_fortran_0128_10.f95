! @file   3D_test_fortran_0128_10.f95
! @author Markovich Dmitry <dmmrkovich at gmail (.) com>
! @copyright 2013 Markovich Dmitry
PROGRAM MAIN
  ! Variables declaration
  integer depth
  integer repeats, size, steps
  integer sum_flag, print_flag, time_flag, size_check_flag
  integer Hx, Chxe, Chxh, Hy, Chyh, Chye, Hz, Chzh, Chze, Ex, Cexe, Cexh, Ey, Ceye, Ceyh, Ez, Ceze, Cezh, SrcEz
  integer i, j, k, c, d, t, task
  integer printflag
  real(4) :: start_time, finish_time
  double precision time, x, sum
  double precision initial, coefficient
  parameter (depth = 2, initial = 0.0, coefficient = 0.1, size = 128, steps = 10)
  double precision, dimension(0:(size+1),0:(size+1),0:(size+1),19,0:1) :: data
  ! End of Variables declaration
  
  ! Variables initialization
  ! Enumerating data components
  ! Fields and coefficients
  Hx = 1
  Chxh = 2
  Chxe = 3
  Hy = 4
  Chyh = 5
  Chye = 6
  Hz = 7
  Chzh = 8
  Chze = 9
  Ex = 10
  Cexe = 11
  Cexh = 12
  Ey = 13
  Ceye = 14
  Ceyh = 15
  Ez = 16
  Ceze = 17
  Cezh = 18
  SrcEz = 19
  ! End Fields and coefficients
  ! Setting Fields
  data(0:(size + 1), 0:(size + 1), 0:(size + 1), Hx, 0:1) = initial
  data(0:(size + 1), 0:(size + 1), 0:(size + 1), Hy, 0:1) = initial
  data(0:(size + 1), 0:(size + 1), 0:(size + 1), Hz, 0:1) = initial
  data(0:(size + 1), 0:(size + 1), 0:(size + 1), Ex, 0:1) = initial
  data(0:(size + 1), 0:(size + 1), 0:(size + 1), Ey, 0:1) = initial
  data(0:(size + 1), 0:(size + 1), 0:(size + 1), Ez, 0:1) = initial
  ! End Setting Fields
  ! Setting Coefficients
  data(0:(size + 1), 0:(size + 1), 0:(size + 1), Chxh, 0:1) = coefficient
  data(0:(size + 1), 0:(size + 1), 0:(size + 1), Chxe, 0:1) = coefficient
  data(0:(size + 1), 0:(size + 1), 0:(size + 1), Chyh, 0:1) = coefficient
  data(0:(size + 1), 0:(size + 1), 0:(size + 1), Chye, 0:1) = coefficient
  data(0:(size + 1), 0:(size + 1), 0:(size + 1), Chzh, 0:1) = coefficient
  data(0:(size + 1), 0:(size + 1), 0:(size + 1), Chze, 0:1) = coefficient
  data(0:(size + 1), 0:(size + 1), 0:(size + 1), Cexe, 0:1) = coefficient
  data(0:(size + 1), 0:(size + 1), 0:(size + 1), Cexh, 0:1) = coefficient
  data(0:(size + 1), 0:(size + 1), 0:(size + 1), Ceye, 0:1) = coefficient
  data(0:(size + 1), 0:(size + 1), 0:(size + 1), Ceyh, 0:1) = coefficient
  data(0:(size + 1), 0:(size + 1), 0:(size + 1), Ceze, 0:1) = 1.0
  data(0:(size + 1), 0:(size + 1), 0:(size + 1), Cezh, 0:1) = coefficient
  ! End Setting Coefficients
  ! End of Variables initialization
  
  call cpu_time(start_time)
  DO t = 0, (steps - 1) ! Calculation cycle 
     ! Setting the source
     x = dble(t) * dble(t) / (dble(steps) * dble(steps))
     data(:, :, :, SrcEz, 0) = exp(-x)
     ! Hx update
     data(1:size,1:size,1:size,Hx,1) = data(1:size,1:size,1:size,Chxh,0) * data(1:size,1:size,1:size,Hx,0) + &
          data(1:size,1:size,1:size,Chxe,0) &
          * ((data(1:size,1:size,2:(size+1),Ey,0) - data(1:size,1:size,1:size,Ey,0)) &
          - (data(1:size,2:(size+1),1:size,Ez,0) - data(1:size,1:size,1:size,Ez,0)))
     ! End Hx update
     ! Hy update
     data(1:size,1:size,1:size,Hy,1) = data(1:size,1:size,1:size,Chyh,0) * data(1:size,1:size,1:size,Hy,0) + &
          data(1:size,1:size,1:size,Chye,0) &
          * ((data(2:(size+1),1:size,1:size,Ez,0) - data(1:size,1:size,1:size,Ez,0)) &
          - (data(1:size,1:size,2:(size+1),Ex,0) - data(1:size,1:size,1:size,Ex,0)))
     ! End Hy update
     ! Hz update
     data(1:size,1:size,1:size,Hz,1) = data(1:size,1:size,1:size,Chzh,0) * data(1:size,1:size,1:size,Hz,0) + &
          data(1:size,1:size,1:size,Chze,0) &
          * ((data(1:size,2:(size+1),1:size,Ex,0) - data(1:size,1:size,1:size,Ex,0)) &
          - (data(2:(size+1),1:size,1:size,Ey,0) - data(1:size,1:size,1:size,Ey,0)))
     ! End Hz update
     ! Ex update
     data(1:size,1:size,1:size,Ex,1) = data(1:size,1:size,1:size,Cexe,0) * data(1:size,1:size,1:size,Ex,0) + &
          data(1:size,1:size,1:size,Cexh,0) &
          * ((data(1:size,1:size,1:size,Hz,0) - data(1:size,0:(size-1),1:size,Hz,0)) &
          - (data(1:size,1:size,1:size,Hy,0) - data(1:size,1:size,0:(size-1),Hy,0)))
     ! End Ex update
     ! Ey update
     data(1:size,1:size,1:size,Ey,1) = data(1:size,1:size,1:size,Ceye,0) * data(1:size,1:size,1:size,Ey,0) + &
          data(1:size,1:size,1:size,Ceyh,0) &
          * ((data(1:size,1:size,1:size,Hx,0) - data(1:size,1:size,0:(size-1),Hx,0)) &
          - (data(1:size,1:size,1:size,Hz,0) - data(0:(size-1),1:size,1:size,Hz,0)))
     ! End Ey update
     ! Ez update
     data(1:size,1:size,1:size,Ez,1) = data(1:size,1:size,1:size,Ceze,0) * data(1:size,1:size,1:size,Ez,0) + &
          data(1:size,1:size,1:size,Cezh,0) &
          * ((data(1:size,1:size,1:size,Hy,0) - data(0:(size-1),1:size,1:size,Hy,0)) &
          - (data(1:size,1:size,1:size,Hx,0) - data(1:size,0:(size-1),1:size,Hx,0))) &
          + data(1:size,1:size,1:size,SrcEz,0)
     ! End Ez update
     ! Cycle arrays
     data(1:size,1:size,1:size,Hx,0) = data(1:size,1:size,1:size,Hx,1)
     data(1:size,1:size,1:size,Hy,0) = data(1:size,1:size,1:size,Hy,1)
     data(1:size,1:size,1:size,Hz,0) = data(1:size,1:size,1:size,Hz,1)
     data(1:size,1:size,1:size,Ex,0) = data(1:size,1:size,1:size,Ex,1)
     data(1:size,1:size,1:size,Ey,0) = data(1:size,1:size,1:size,Ey,1)
     data(1:size,1:size,1:size,Ez,0) = data(1:size,1:size,1:size,Ez,1)
     ! End Cycle arrays
  END DO  ! End calculation cycle
  call cpu_time(finish_time)
  time = finish_time-start_time
  PRINT *,'>> Task',size,'x',size,'x',size,'x',steps,'(steps)'
  PRINT *, '>> Fortran 3D took',time,'s'
  ! Calculating Ez control sum
  sum = 0.0
  DO k = 0, (size + 1)
     DO j = 0, (size + 1)
        DO i = 0, (size + 1)
           sum = sum + dabs(data(i,j,k,Ez,0))
        END DO
     END DO
  END DO
  !PRINT *,'>> Ez control sum is',sum
  ! End Calculating Ez control sum
  STOP
END PROGRAM MAIN

   
