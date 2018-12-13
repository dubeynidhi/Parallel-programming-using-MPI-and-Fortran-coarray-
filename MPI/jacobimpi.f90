!  PROGRAMMER
!  Nidhi Dubey, ndubey2@asu.edu
!  jacobian iteration implementation using mpi
use, intrinsic::iso_fortran_env  ! Fortran 2003 and later
use mpi
implicit none
integer, parameter:: ROOT=0
integer, parameter:: INTERNAL_ERROR=1  ! shell error code; must be nonzero
integer, parameter:: DOUBLE=kind(1.0d0)  ! double precision

!  Local variables

integer:: ierr,npes,M,maxiter,boundary
integer:: me  ! my MPI rank (0 = root process)
real(DOUBLE)::f(0:255,0:255),t(0:255,0:255)
real x,i,j,y
M=255
maxiter=5000

! change boundary here
boundary=0

call mpi_init(ierr)  ! should be first executable statement
call mpi_comm_rank(MPI_COMM_WORLD, me, ierr)  ! my rank
call mpi_comm_size(MPI_COMM_WORLD, npes, ierr) ! number of processes

if(me.eq.ROOT) then

!calculate true solution
    t = boundary
    do i=1,m-1
        x=i/(m)
        do j=1,m-1
            y=j/(m)
            t(i,j)=((exp(x+y))*((x ** 2) - x)*((y ** 2) - y))
        enddo
    enddo
    sync all

!calculate f matrix
    f = boundary
    do i=1,m-1
        x=i/(m)
        do j=1,m-1
            y=j/(m)
            f(i,j)=((exp(x+y))*((((x ** 2) + 3*x)*((y ** 2) - y)) + (((y ** 2) + 3*y)*((x ** 2) - x))))
        enddo
    enddo
end if

! broadcast true solution and f array to all 16 processes
call mpi_bcast(t, size(t), MPI_DOUBLE, ROOT, MPI_COMM_WORLD, ierr)
call mpi_bcast(f, size(f), MPI_DOUBLE, ROOT, MPI_COMM_WORLD, ierr)

! call function to compute u
call jacobian(me,boundary,maxiter,f,t)
call mpi_finalize(ierr)
! should be last executable statement before STOP
stop
!---------   End of normal computation

911 continue
call mpi_abort(MPI_COMM_WORLD, INTERNAL_ERROR, ierr)
end


subroutine jacobian(me,boundary,maxiter,f,t)
use mpi
integer, parameter:: DOUBLE=kind(1.0d0)

integer, intent(in)::me,maxiter,boundary
real(DOUBLE), intent(in)::f(0:255,0:255)
real(DOUBLE), intent(in)::t(0:255,0:255)
integer:: status(MPI_STATUS_SIZE), ierr
integer::i,j,stop_iter,send_up,receive_from_down,send_down,receive_from_up,send_right,receive_from_left,send_left,receive_from_right
integer:: k,beg_x, end_x, beg_y, end_y, f_x_beg,f_y_beg,f_x_end,f_y_end

! buffers to send data up, down and left and right
real(DOUBLE):: receive_buf_up_down(66)
real(DOUBLE):: send_buf_up_down(66)

real(DOUBLE):: receive_buf_up_down1(66)
real(DOUBLE):: send_buf_up_down1(66)

real(DOUBLE):: receive_buf_left_right(66)
real(DOUBLE):: send_buf_left_right(66)

real(DOUBLE):: receive_buf_left_right1(66)
real(DOUBLE):: send_buf_left_right1(66)

! to store final results and reduced final results
real(DOUBLE):: results(2),final_results(2)

! define u_old, u_new and temp(for delta calculation)
real(DOUBLE), dimension(:,:), allocatable:: u_new,u_old,temp

! each subgrid of u has 66*66 elements which includes ghost points for all four directions
1 the boundaries are kept constant throughout the code and inner elements are updates
allocate(u_new(0:65,0:65),u_old(0:65,0:65),temp(0:65,0:65))

! to indicate buffer send receive direction.
send_up = 0
receive_from_down = 0

send_down = 1
receive_from_up = 1

send_right = 2
receive_from_left = 2

send_left = 3
receive_from_right = 3

u_old=0
temp=0

! we need beginning and end index for all the 4*4 subgrids also we need to initialize boundary. We also need the corresponding elements in
! f and true solution matrix

! update first row of top processors with boundary and find start and end index for !these processors
if((me==1).or.(me==2).or.(me==3).or.(me==4)) then
    beg_x=2
    end_x=64
    f_x_beg=1
    f_x_end=63

    u_old(1,:)=boundary
end if

! update last row of bottom processors with boundary and find start and end index for !these processors
if((me==12).or.(me==13).or.(me==14).or.(me==15)) then
    beg_x=1
    end_x=63
    f_x_beg=192
    f_x_end=254

    u_old(64,:)=boundary
end if

! find start and end index for these processors
if((me==4).or.(me==5).or.(me==6).or.(me==7)) then
    beg_x=1
    end_x=64
    f_x_beg=64
    f_x_end=127
end if

! find start and end index for these processors
if((me==8).or.(me==9).or.(me==10).or.(me==11)) then
    beg_x=1
    end_x=64
    f_x_beg=128
    f_x_end=191
end if

! update first column of left processors and find start and end index for these processors
if((me==0).or.(me==4).or.(me==8).or.(me==12)) then
    u_old(:,1)=boundary
    beg_y=2
    end_y=64
    f_y_beg=1
    f_y_end=63
end if

! update last column of right processors and find start and end index for these processors
if((me==3).or.(me==7).or.(me==11).or.(me==15)) then
    u_old(:,64)=boundary
    beg_y=1
    end_y=63
    f_y_beg=192
    f_y_end=254
end if

! find start and end index for these processors
if((me==1).or.(me==5).or.(me==9).or.(me==13)) then
    beg_y=1
    end_y=64
    f_y_beg=64
    f_y_end=127
end if

! find start and end index for these processors
if((me==2).or.(me==6).or.(me==10).or.(me==14)) then
    beg_y=1
    end_y=64
    f_y_beg=128
    f_y_end=191
end if


! initialise u_new and temp too
u_new = u_old
temp = u_old

! now run till maximum iterations and calculate the value in u_new and also send receive data
do k=1,maxiter

    ! implement jacobi scheme
    do i=beg_x,end_x
        do j=beg_y,end_y
            u_new(i,j) = ((u_old(i+1,j) + u_old(i-1,j) + u_old(i,j+1) + u_old(i,j-1) - ((f(f_x_beg+i-beg_x,f_y_beg+j-beg_y))/(255 ** 2)))/4)
        enddo
    enddo

    ! check for convergence after every 10 iterations
    if(mod(k,10).eq.0) then
        ! call this function to calculate delta and c
        call calc_delta_c(temp,u_new,t,results, beg_x, beg_y, end_x, end_y,f_y_beg,f_x_beg)
        call mpi_barrier(MPI_COMM_WORLD, ierr)  ! wait for everyone to finish to sync
        ! compute final results by all the processes
        call mpi_allreduce(results,final_results,2,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD)

        ! compare with set delta to continue or exit
        if(final_results(1).lt.0.0001) then
            ! note iteration number and then  break
            stop_iter=k
            goto 911
        else
            ! else store u_new values in temp
            temp=u_new
        endif
    endif

    ! send data up for all except processors 0,1,2,3
    if ((me .ne. 0).and.(me .ne. 1).and.(me .ne. 2).and.(me .ne. 3)) then
        send_buf_up_down=u_old(beg_y,:)
        call MPI_Sendrecv( send_buf_up_down, 66, MPI_DOUBLE,(me-4),send_up, receive_buf_up_down,66, MPI_DOUBLE,(me-4),receive_from_up, MPI_COMM_WORLD, status, ierr)
        u_old(65,:)=receive_buf_up_down
    endif

    ! send data down for all except processors 12,13,14,15
    if ((me .ne. 12).and.(me .ne. 13).and.(me .ne. 14).and.(me .ne. 15)) then
        send_buf_up_down1=u_old(end_y,:)
        call MPI_Sendrecv( send_buf_up_down1, 66,MPI_DOUBLE,(me+4),send_down, receive_buf_up_down1,66,MPI_DOUBLE,(me+4), receive_from_down, MPI_COMM_WORLD, status, ierr)
        u_old(0,:)=receive_buf_up_down1
    endif

    ! send data left for all except processors 0,4,8,12
    if((me.ne.0).and.(me.ne.4).and.(me.ne.8).and.(me.ne.12)) then
        send_buf_left_right=u_old(:,beg_y)
        call MPI_Sendrecv( send_buf_left_right,66,MPI_DOUBLE,(me-1),send_left, receive_buf_left_right,66,MPI_DOUBLE,(me-1), receive_from_left, MPI_COMM_WORLD, status, ierr)
        u_new(:,65)=receive_buf_left_right
    endif

     ! send data right for all except processors 3,7,11,15
     if((me.ne.3).and.(me.ne.7).and.(me.ne.11).and.(me.ne.15)) then
        receive_buf_left_right1=u_old(:,end_y)
        call MPI_Sendrecv( send_buf_left_right1,66,MPI_DOUBLE,(me+1),send_right, receive_buf_left_right1,66,MPI_DOUBLE,(me+1), receive_from_right, MPI_COMM_WORLD, status, ierr)
        u_old(:,0)=receive_buf_left_right1
    endif

    ! old will be updated now
    u_old=u_new

enddo



! print results
911 continue
if(me.eq.0) then
    !write(*,*) "Stopped after ", stop_iter, "iterations"
    write(*,*),'delta is :' , final_results(1)
    write(*,*) "c is :", final_results(2)
end if
ierr=0
return

end subroutine jacobian

! to calculate delta and c
subroutine calc_delta_c(temp,u_new,t,results,beg_x, beg_y, end_x, end_y,f_y_beg,f_x_beg)
use, intrinsic::iso_fortran_env
integer, parameter:: DOUBLE=kind(1.0d0)  ! double precision
real(DOUBLE),intent(in):: u_new(0:65,0:65)
real(DOUBLE),intent(in):: temp(0:65,0:65)
real(DOUBLE),intent(in):: t(0:255,0:255)
real(DOUBLE),intent(out)::results(2)
integer, intent(in)::beg_x, beg_y, end_x, end_y,f_y_beg,f_x_beg
real(DOUBLE)::dist
real(DOUBLE)::c
real(DOUBLE)::delta
integer::j
integer::k
dist = 0
c = 0
delta = 0

! calculation for delta. Subtract from temp and calculate maximum
do j=beg_x,end_x
    do k=beg_y,end_y
        dist = u_new(j,k) - temp(j,k)
        dist = (dist ** 2)
        dist = sqrt(dist)
        if(dist.gt.delta) then
            delta = dist
        endif
    enddo
enddo

! calculation for c. Subtract from true solution and calculate maximum
do j=beg_x,end_x
    do k=beg_y,end_y
        dist = u_new(j,k) - t(f_x_beg+j-beg_x,f_y_beg+k-beg_y)
        dist = (dist ** 2)
        dist = sqrt(dist)
        if(dist.gt.c) then
            c = dist
        endif
    enddo
enddo

!store and return results
results(1) = delta
results(2) = c

return
end subroutine calc_delta_c
