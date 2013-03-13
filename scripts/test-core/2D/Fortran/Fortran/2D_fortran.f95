! @file   2D_fortran.f95
! @author Markovich Dmitry <dmmrkovich at gmail (.) com>
! @copyright 2013 Markovich Dmitry
PROGRAM MAIN
  ! Variables declaration
  integer depth
  integer repeats, size, steps
  integer sum_flag, print_flag, time_flag, size_check_flag
  integer Ez, Ceze, Cezh, Hy, Chyh, Chye, Hx, Chxe, Chxh, SrcEz
  integer i, j, c, d, t, task
  integer printflag
  real(4) :: start_time, finish_time
  double precision time, x, sum
  double precision initial, coefficient
  parameter (depth = 2, initial = 0.0, coefficient = 0.1d+0, size = sz, steps = st)
  ! Fields and coefficients
  parameter (Hx = 1, Chxh = 2, Chxe = 3, Ez = 4, Ceze = 5, Cezh = 6, Hy = 7, Chyh = 8, Chye = 9, SrcEz = 10)
  ! End Fields and coefficients
  ! End of Variables declaration
  double precision, dimension(0:(size+1),0:(size+1),SrcEz,0:1) :: data
  ! Variables initialization
  ! Setting Fields
  data(:, :, Hx, :) = initial
  data(:, :, Hy, :) = initial
  data(:, :, Ez, :) = initial
  ! End Setting Fields
     ! Setting Coefficients
  data(:,:, Chxh, :) = coefficient
  data(:,:, Chxe, :) = coefficient
  data(:,:, Ceze, :) = 1.0
  data(:,:, Cezh, :) = coefficient
  data(:,:, Chyh, :) = coefficient
  data(:,:, Chye, :) = coefficient
  ! End Setting Coefficients
  ! End Enumerating data components
  ! End of Variables initialization
  call cpu_time(start_time)
     DO t = 0, (-1 + steps/2) ! Calculation cycle  
        ! Setting the source
	x = dble(t*2) * dble(t*2) / (dble(steps) * dble(steps))
     	data(:, :, SrcEz, 0) = exp(-x)
        ! End Setting the source
        ! Hx update
        data(1:size,1:size,Hx,1) = data(1:size,1:size,Chxh,0) * data(1:size,1:size,Hx,0) &
             - data(1:size,1:size,Chxe,0) * (data(1:size,2:(size+1),Ez,0) - data(1:size,1:size,Ez,0))
        ! End Hx update
        ! Hy update
        data(1:size,1:size,Hy,1) = data(1:size,1:size,Chyh,0) * data(1:size,1:size,Hy,0) &
             + data(1:size,1:size,Chye,0) * (data(2:(size+1),1:size,Ez,0) - data(1:size,1:size,Ez,0))
        ! End Hy update
        ! Ez update
        data(1:size,1:size,Ez,1) = data(1:size,1:size,Ceze,0) * data(1:size,1:size,Ez,0) &
             + data(1:size,1:size,Cezh,0) * ((data(1:size,1:size,Hy,1) - data(0:(size-1),1:size,Hy,1)) &
             - (data(1:size,1:size,Hx,1) - data(1:size,0:(size-1),Hx,1))) + data(1:size,1:size,SrcEz,0)
        ! End Ez update

       ! Setting the source
       x = dble(t*2+1) * dble(t*2+1) / (dble(steps) * dble(steps))
       data(:, :, SrcEz, 1) = exp(-x)
       ! Hx update
       data(1:size,1:size,Hx,0) = data(1:size,1:size,Chxh,1) * data(1:size,1:size,Hx,1) &
             - data(1:size,1:size,Chxe,1) * (data(1:size,2:(size+1),Ez,1) - data(1:size,1:size,Ez,1))	
     ! End Hx update
     ! Hy update
        data(1:size,1:size,Hy,0) = data(1:size,1:size,Chyh,1) * data(1:size,1:size,Hy,1) &
             + data(1:size,1:size,Chye,1) * (data(2:(size+1),1:size,Ez,1) - data(1:size,1:size,Ez,1))
     ! End Hy update
        ! Ez update
        data(1:size,1:size,Ez,0) = data(1:size,1:size,Ceze,1) * data(1:size,1:size,Ez,1) &
             + data(1:size,1:size,Cezh,1) * ((data(1:size,1:size,Hy,0) - data(0:(size-1),1:size,Hy,0)) &
             - (data(1:size,1:size,Hx,0) - data(1:size,0:(size-1),Hx,0))) + data(1:size,1:size,SrcEz,1)
        ! End Ez update
     END DO ! End calculation cycle
     call cpu_time(finish_time)
     time = finish_time-start_time
     PRINT *,
     PRINT *,'>> 2D Fortran',size
	PRINT *,'>> 2D Fortran took',time,'s'
     IF (cs == 1) THEN
      ! Calculating Ez control sum
      sum = 0.0
      DO j = 0, (size + 1)
         DO i = 0, (size + 1)
            sum = sum + data(i,j,Ez,0)
         END DO
      END DO
         PRINT *,'>> 2D Fortran control sum for Ez is',sum
      ! End Calculating Ez control sum 
     END IF
     STOP
   END PROGRAM MAIN
   
