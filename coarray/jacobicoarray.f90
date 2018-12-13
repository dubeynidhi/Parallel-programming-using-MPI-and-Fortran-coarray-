!  PROGRAMMER
!  Nidhi Dubey, ndubey2@asu.edu
! jacobian implementation using coarray.
program jacobicoarray
    use, intrinsic::iso_fortran_env  ! Fortran 2003 and later
    implicit none

    integer m,k,c1,c2,q,beg_x, end_x, beg_y, end_y, f_x_beg,f_y_beg,f_x_end,f_y_end,stop_iter
    ! u is 3-D, as I am storing new and old with same variable, and temp and f are 2D 256*256 matrix
    ! each subgrid of u has 66*66 elements which includes ghost points for all four directions
    ! the boundaries are kept constant throughout the code and inner elements are updates
    real u(0:65,0:65,0:1)[*], f(0:255,0:255),t(0:255,0:255),temp(0:65,0:65), d[*],del[*]
    integer total_imgs, curr_img, old, new, maxiter, boundary
    real x,i,j,c,y,delta

    c1=1
    c2=64
    q=10
    m=255
    new = 1
    old = 1-new
    maxiter=5000
    d=0
    delta=0
    c=0

    ! change boundary here
    boundary=0

    !get total number of images and current image
    total_imgs = NUM_IMAGES()
    curr_img = THIS_IMAGE()

    ! calculate true solution
    t = boundary ! set to boundary and compute for inner indices
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
    sync all
    u=0
    temp=0

    ! we need beginning and end index for all the 4*4 subgrids also we need to initialize boundary. We also need the corresponding elements in
    ! f and true solution matrix


   ! update first row of top processors with boundary and find start and end index for !these processors
    if((curr_img==1).or.(curr_img==2).or.(curr_img==3).or.(curr_img==4)) then
        beg_x=2
        end_x=64
        f_x_beg=1
        f_x_end=63

        u(1,:,:)=boundary
    end if

  ! update last row of bottom processors with boundary and find start and end index for !these processors
    if((curr_img==13).or.(curr_img==14).or.(curr_img==15).or.(curr_img==16)) then
        beg_x=1
        end_x=63
        f_x_beg=192
        f_x_end=254
        u(64,:,:)=boundary
    end if

    ! find start and end index for these processors
    if((curr_img==5).or.(curr_img==6).or.(curr_img==7).or.(curr_img==8)) then
        beg_x=1
        end_x=64
         f_x_beg=64
         f_x_end=127
    end if

    ! find start and end index for these processors
    if((curr_img==9).or.(curr_img==10).or.(curr_img==11).or.(curr_img==12)) then
        beg_x=1
        end_x=64
        f_x_beg=128
        f_x_end=191
    end if

    ! update first column of left processors and find start and end index for these processors
    if((curr_img==1).or.(curr_img==5).or.(curr_img==9).or.(curr_img==13)) then
        u(:,1,:)=boundary
        beg_y=2
        end_y=64
        f_y_beg=1
        f_y_end=63

    end if

    ! update last column of right processors and find start and end index for these processors
    if((curr_img==4).or.(curr_img==8).or.(curr_img==12).or.(curr_img==16)) then
        u(:,64,:)=boundary
        beg_y=1
        end_y=63
        f_y_beg=192
        f_y_end=254
    end if

    ! find start and end index for these processors
    if((curr_img==2).or.(curr_img==6).or.(curr_img==10).or.(curr_img==14)) then
        beg_y=1
        end_y=64
         f_y_beg=64
          f_y_end=127
    end if

    ! find start and end index for these processors
    if((curr_img==3).or.(curr_img==7).or.(curr_img==11).or.(curr_img==15)) then
        beg_y=1
        end_y=64
         f_y_beg=128
          f_y_end=191
    end if


! now run till maximum iterations and calculate the value in u_new and also send receive data
do k=1, maxiter

      ! send data right for all except processors 4,8,12,16
    if ((curr_img .ne. 4).and.(curr_img .ne. 8).and.(curr_img .ne. 12).and.(curr_img .ne. 16)) then
        u(:,c1-1,old)[curr_img+1] = u(:,c2,old)
    end if

    ! send data left for all except processors 1,5,9,13
    if ((curr_img .ne. 1).and.(curr_img .ne. 5).and.(curr_img .ne. 9).and.(curr_img .ne. 13)) then
        u(:,c2+1,old)[curr_img-1] = u(:,c1,old)
    end if

      ! send data down for all except processors 13,14,15,16
    if ((curr_img .ne. 13).and.(curr_img .ne. 14).and.(curr_img .ne. 15).and.(curr_img .ne. 16)) then
        u(c1-1,:,old)[curr_img+4] = u(c2,:,old)
    end if

     ! send data up for all except processors 1,2,3,4
    if ((curr_img .ne. 1).and.(curr_img .ne. 2).and.(curr_img .ne. 3).and.(curr_img .ne. 4)) then
        u(c2+1,:,old)[curr_img-4] = u(c1,:,old)
    end if

sync all

! update u_new using jacobian equation
do i=beg_x,end_x
    do j=beg_y,end_y
    u(i,j,new)=((u(i+1,j,old)+u(i-1,j,old)+u(i,j+1,old)+u(i,j-1,old)-((f(f_x_beg+i-beg_x,f_y_beg+j-beg_y))/(m ** 2)))/4)
    enddo
enddo

! check for convergence after every 10 iterations
if(mod(k,q).eq.0) then

    !calculate delta for all 16 images
    del[curr_img]=maxval(SQRT(((u(beg_x:end_x,beg_y:end_y,new)-temp(beg_x:end_x,beg_y:end_y))**2)))
    sync all

    ! in image 1 find maximum delta from all images
    if (this_image() == 1) then
        do i=1,total_imgs
            if (delta .lt. del[i]) then
                delta=del[i]
            end if
        end do
    end if

    ! check for convergence, if delta less than 0.0001 break. Else save u into temp
    if(delta.lt. 0.0001) then
        stop_iter=k
        goto 911
    else
        temp=u(:,:,new)
    end if

end if

    ! swap old and new
     new = old
     old = 1-new
end do

911 continue

! after the iterations are done calculate c for all images ( compare with true solution
d[curr_img]=maxval(SQRT(((u(beg_x:end_x,beg_y:end_y,new)-t(f_x_beg:f_x_end,f_y_beg:f_y_end))**2)))
sync all

! in image 1 find maximum delta from all images
if (this_image() == 1) then
    do i=1,total_imgs
        if (c .lt. d[i]) then
            c=d[i]
        end if
    end do
end if

! print results in image 1
if(curr_img .eq. 1) then
!write(*,*) "Printing results in image 1"
!write(*,*) "Stopped after ",k, "iterations"
write(*,*),'delta is :',delta
write(*,*) "c is:",c
end if


end program jacobicoarray
