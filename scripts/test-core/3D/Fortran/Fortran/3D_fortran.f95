! @file   3D_fortran.f95
! @author Markovich Dmitry <dmmrkovich at gmail (.) com>
! @copyright 2013 Markovich Dmitry
PROGRAM MAIN
  ! Variables declaration
  integer depth
  integer repeats, size, steps
  integer sum_flag, print_flag, time_flag, size_check_flag
  integer Hx, Chxe, Chxh, Hy, Chyh, Chye, Hz, Chzh, Chze, Ex, Cexe, Cexh, Ey, Ceye, Ceyh, Ez, Ceze, Cezh, SrcEz, toPrint
  integer i, j, k, c, d, t, task
  integer printflag
  real(4) :: start_time, finish_time
  double precision time, x, sum
  double precision initial, coefficient, boundary
  parameter (depth = 2, initial = 0.0d+0, boundary = 0.0, coefficient = 0.1d+0, size = sz, steps = st, full_output = 0)
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
  ! Setting Boundary
  data(0, 0:(size + 1), 0:(size + 1), Hx, 0:1) = boundary
  data(size + 1, 0:(size + 1), 0:(size + 1), Hx, 0:1) = boundary
  data(0:(size + 1), 0, 0:(size + 1), Hx, 0:1) = boundary
  data(0:(size + 1), size + 1, 0:(size + 1), Hx, 0:1) = boundary
  data(0:(size + 1), 0:(size + 1), 0, Hx, 0:1) = boundary
  data(0:(size + 1), 0:(size + 1), size + 1, Hx, 0:1) = boundary

  data(0, 0:(size + 1), 0:(size + 1), Hy, 0:1) = boundary
  data(size + 1, 0:(size + 1), 0:(size + 1), Hy, 0:1) = boundary
  data(0:(size + 1), 0, 0:(size + 1), Hy, 0:1) = boundary
  data(0:(size + 1), size + 1, 0:(size + 1), Hy, 0:1) = boundary
  data(0:(size + 1), 0:(size + 1), 0, Hy, 0:1) = boundary
  data(0:(size + 1), 0:(size + 1), size + 1, Hy, 0:1) = boundary

  data(0, 0:(size + 1), 0:(size + 1), Hz, 0:1) = boundary
  data(size + 1, 0:(size + 1), 0:(size + 1), Hz, 0:1) = boundary
  data(0:(size + 1), 0, 0:(size + 1), Hz, 0:1) = boundary
  data(0:(size + 1), size + 1, 0:(size + 1), Hz, 0:1) = boundary
  data(0:(size + 1), 0:(size + 1), 0, Hz, 0:1) = boundary
  data(0:(size + 1), 0:(size + 1), size + 1, Hz, 0:1) = boundary

  data(0, 0:(size + 1), 0:(size + 1), Ex, 0:1) = boundary
  data(size + 1, 0:(size + 1), 0:(size + 1), Ex, 0:1) = boundary
  data(0:(size + 1), 0, 0:(size + 1), Ex, 0:1) = boundary
  data(0:(size + 1), size + 1, 0:(size + 1), Ex, 0:1) = boundary
  data(0:(size + 1), 0:(size + 1), 0, Ex, 0:1) = boundary
  data(0:(size + 1), 0:(size + 1), size + 1, Ex, 0:1) = boundary

  data(0, 0:(size + 1), 0:(size + 1), Ey, 0:1) = boundary
  data(size + 1, 0:(size + 1), 0:(size + 1), Ey, 0:1) = boundary
  data(0:(size + 1), 0, 0:(size + 1), Ey, 0:1) = boundary
  data(0:(size + 1), size + 1, 0:(size + 1), Ey, 0:1) = boundary
  data(0:(size + 1), 0:(size + 1), 0, Ey, 0:1) = boundary
  data(0:(size + 1), 0:(size + 1), size + 1, Ey, 0:1) = boundary

  data(0, 0:(size + 1), 0:(size + 1), Ez, 0:1) = boundary
  data(size + 1, 0:(size + 1), 0:(size + 1), Ez, 0:1) = boundary
  data(0:(size + 1), 0, 0:(size + 1), Ez, 0:1) = boundary
  data(0:(size + 1), size + 1, 0:(size + 1), Ez, 0:1) = boundary
  data(0:(size + 1), 0:(size + 1), 0, Ez, 0:1) = boundary
  data(0:(size + 1), 0:(size + 1), size + 1, Ez, 0:1) = boundary
  ! End Setting boundary
  ! Setting Fields
  data(1:size, 1:size, 1:size, Hx, 0:1) = initial
  data(1:size, 1:size, 1:size, Hy, 0:1) = initial
  data(1:size, 1:size, 1:size, Hz, 0:1) = initial
  data(1:size, 1:size, 1:size, Ex, 0:1) = initial
  data(1:size, 1:size, 1:size, Ey, 0:1) = initial
  data(1:size, 1:size, 1:size, Ez, 0:1) = initial
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
  DO t = 0, (steps - 1), 2 ! Calculation cycle 
     !Setting the source
     x = dble(t) * dble(t) / (dble(steps) * dble(steps))
     ! x = dble(t*2) * dble(t*2) / (dble(steps) * dble(steps))
     data(:, :, :, SrcEz, 0) = exp(-x)
     ! Hx update
 data(1:size,1:size,1:size,Hx,1) = data(1:size,1:size,1:size,Chxh,0) &
      * data(1:size,1:size,1:size,Hx,0) + data(1:size,1:size,1:size,Chxe,0)&
      * ((data(1:size,1:size,2:(size+1),Ey,0) - data(1:size,1:size,1:size,Ey,0))&
       - (data(1:size,2:(size+1),1:size,Ez,0) - data(1:size,1:size,1:size,Ez,0)))
 ! End Hx update
 ! Hy update
 data(1:size,1:size,1:size,Hy,1) = data(1:size,1:size,1:size,Chyh ,0) &
      * data(1:size,1:size,1:size,Hy,0) +  data(1:size,1:size,1:size,Chye,0)&
      * ((data(2:(size+1),1:size,1:size,Ez,0) - data(1:size ,1:size,1:size,Ez,0))&
       - (data(1:size,1:size ,2:(size+1),Ex,0) - data(1:size,1:size,1:size,Ex,0))) 
 ! End Hy update
 ! Hz update
 data(1:size,1:size,1:size,Hz,1) = data(1:size,1:size,1:size,Chzh,0)&
      * data(1:size,1:size,1:size,Hz,0) + data(1:size,1:size,1:size,Chze,0)&
      * ((data(1:size,2:(size+1),1:size,Ex,0) - data(1:size,1:size,1:size,Ex,0))&
       - (data(2:(size+1),1:size,1:size,Ey,0) - data(1:size,1:size,1:size,Ey,0)))
 ! End Hz update
 ! Ex update
 data(1:size,1:size,1:size,Ex,1) = data(1:size,1:size,1:size,Cexe,0)&
      * data(1:size,1:size,1:size,Ex,0) + data(1:size,1:size,1:size,Cexh,0)&
      * ((data(1:size,1:size,1:size,Hz,1) - data(1:size,0:(size-1),1:size,Hz,1))&
       - (data(1:size,1:size,1:size,Hy,1) - data(1:size,1:size,0:(size-1),Hy,1)))
 ! End Ex update
 ! Ey update
 data(1:size,1:size,1:size,Ey,1) = data(1:size,1:size,1:size,Ceye,0)&
      * data(1:size,1:size,1:size,Ey,0) + data(1:size,1:size,1:size,Ceyh,0)&
      * ((data(1:size,1:size,1:size,Hx,1) - data(1:size,1:size,0:(size-1),Hx,1))&
       - (data(1:size,1:size,1:size,Hz,1) - data(0:(size-1),1:size,1:size,Hz,1)))
 ! End Ey update
 ! Ez update
 data(1:size,1:size,1:size,Ez,1) = data(1:size,1:size,1:size,Ceze,0)&
      * data(1:size,1:size,1:size,Ez,0) + data(1:size,1:size,1:size,Cezh,0)&
      * ((data(1:size,1:size,1:size,Hy,1) - data(0:(size-1),1:size,1:size,Hy,1))&
       - (data(1:size,1:size,1:size,Hx,1) - data(1:size,0:(size-1),1:size,Hx,1)))&
       + data(1:size,1:size,1:size,SrcEz,0)
 !End DO update
     ! Setting the source
     x = dble(t+1) * dble(t+1) / (dble(steps) * dble(steps))
     data(:, :, :, SrcEz, 1) = exp(-x)
     ! Hx update
     data(1:size,1:size,1:size,Hx,0) = data(1:size,1:size,1:size,Chxh,1) * data(1:size,1:size,1:size,Hx,1) + &
          data(1:size,1:size,1:size,Chxe,1) &
          * ((data(1:size,1:size,2:(size+1),Ey,1) - data(1:size,1:size,1:size,Ey,1)) &
          - (data(1:size,2:(size+1),1:size,Ez,1) - data(1:size,1:size,1:size,Ez,1)))
     ! End Hx update
     ! Hy update
     data(1:size,1:size,1:size,Hy,0) = data(1:size,1:size,1:size,Chyh,1) * data(1:size,1:size,1:size,Hy,1) + &
          data(1:size,1:size,1:size,Chye,1) &
          * ((data(2:(size+1),1:size,1:size,Ez,1) - data(1:size,1:size,1:size,Ez,1)) &
          - (data(1:size,1:size,2:(size+1),Ex,1) - data(1:size,1:size,1:size,Ex,1)))
     ! End Hy update
     ! Hz update
     data(1:size,1:size,1:size,Hz,0) = data(1:size,1:size,1:size,Chzh,1) * data(1:size,1:size,1:size,Hz,1) + &
          data(1:size,1:size,1:size,Chze,1) &
          * ((data(1:size,2:(size+1),1:size,Ex,1) - data(1:size,1:size,1:size,Ex,1)) &
          - (data(2:(size+1),1:size,1:size,Ey,1) - data(1:size,1:size,1:size,Ey,1)))
     ! End Hz update
     ! Ex update
     data(1:size,1:size,1:size,Ex,0) = data(1:size,1:size,1:size,Cexe,1) * data(1:size,1:size,1:size,Ex,1) + &
          data(1:size,1:size,1:size,Cexh,1) &
          * ((data(1:size,1:size,1:size,Hz,0) - data(1:size,0:(size-1),1:size,Hz,0)) &
          - (data(1:size,1:size,1:size,Hy,0) - data(1:size,1:size,0:(size-1),Hy,0)))
     ! End Ex update
     ! Ey update
     data(1:size,1:size,1:size,Ey,0) = data(1:size,1:size,1:size,Ceye,1) * data(1:size,1:size,1:size,Ey,1) + &
          data(1:size,1:size,1:size,Ceyh,1) &
          * ((data(1:size,1:size,1:size,Hx,0) - data(1:size,1:size,0:(size-1),Hx,0)) &
          - (data(1:size,1:size,1:size,Hz,0) - data(0:(size-1),1:size,1:size,Hz,0)))
     ! End Ey update
     ! Ez update
     data(1:size,1:size,1:size,Ez,0) = data(1:size,1:size,1:size,Ceze,1) * data(1:size,1:size,1:size,Ez,1) + &
          data(1:size,1:size,1:size,Cezh,1) &
          * ((data(1:size,1:size,1:size,Hy,0) - data(0:(size-1),1:size,1:size,Hy,0)) &
          - (data(1:size,1:size,1:size,Hx,0) - data(1:size,0:(size-1),1:size,Hx,0))) &
          + data(1:size,1:size,1:size,SrcEz,1)
     ! End Ez update
 END DO  ! End calculation cycle
  call cpu_time(finish_time)
  time = finish_time-start_time
  PRINT *,
  PRINT *,'>> 3D Fortran ',size,' steps ',steps
  PRINT *, '>> 3D Fortran took',time,'s'
   ! Calculating Ez control sum
   sum = 0.0
   DO k = 0, (size+1)
      DO j = 0, (size+1)
         DO i = 0, (size+1)
            sum = sum + data(k,j,i,Ez,0)
         END DO
      END DO
   END DO
   PRINT *,'>> 3D Fortran control sum for Ez is',sum
   IF (cs == 1) THEN
      !toPrint = Cezh
      toPrint = Ez
      !toPrint = Ex
      !toPrint = Ey
      !toPrint = Hz
      !toPrint = Hx
      !toPrint = Hy
      DO k = 0, (size+1)
         DO j = 0, (size+1)
           ! print *, data(k,j,:,toPrint,0)
         END DO
        ! print *,
      END DO
      ! End Calculating Ez control sum
   END IF
   STOP
 END PROGRAM MAIN

   
