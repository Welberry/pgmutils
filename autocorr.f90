program autocorr

  ! A program to calculate the autocorrelation coefficient of an image
  ! for a number of vectors

  use iso_varying_string
  use cmdline_arguments, only: assignment(=), have_args, get_options, next_arg, option_exists, get_value, has_value
  use file_functions, only: stderr
  use pnm_class
  use precision, only: rd_kind

  implicit none

  character(len=1000) :: fname

  type (pnm_object) :: pnm

  integer :: i, j, k, l, m, error, nx, ny, nxny, weight

  real :: datamean, sumsqr
  real(kind=rd_kind) :: moran, geary

  integer, allocatable, dimension(:,:) :: data
  real, allocatable, dimension(:,:) :: shifted, datanorm
  logical, allocatable, dimension(:,:) :: mask, shiftedmask

  type (varying_string) :: myoptions(3)

  logical :: debug = .FALSE.

  ! Populate our 
  myoptions(1) = "help"
  myoptions(2) = "debug"
  myoptions(3) = "num"

  call get_options(myoptions, error)
  if (error > 0) stop 'Whoops error processing options!'

  if ( option_exists('help') ) then
     call usage
     stop
  end if

  if ( .not. have_args() ) then
     write(stderr,*) "ERROR! Must provide a pgm file as command line argument"
     call usage
     stop
  end if

  debug = option_exists("debug")

  do while (have_args())

     ! The files are on the command line
     fname = next_arg()

     call read(pnm,trim(fname))

     nx = size(pnm,1)
     ny = size(pnm,2)

     nxny = nx*ny

     ! Allocate the space required for the data
     allocate(data(nx,ny),shifted(nx,ny),datanorm(nx,ny),mask(nx,ny),shiftedmask(nx,ny))
     
     data = pnm

     datamean = mean(data, data > 0)

     datanorm = real(data) - datamean

     mask = data > 0

     sumsqr = sum( datanorm**2, mask = mask ) 

     if (debug) print *,'mean = ', datamean
     if (debug) print *,count(mask)
     if (debug) print *,'sumsqr = ',sumsqr
     if (debug) print *,'nxny/sumsqr = ',real(nxny)/sumsqr

     if (debug) print *,'sum(data**2)', sum( data**2, mask = mask ) 

     do i = 0, nx - 1

        shifted = eoshift(datanorm,i,dim=1)
        shiftedmask = eoshift(mask,i,dim=1)
        if (debug) print *,count(shiftedmask)

        moran = 0
        geary = 0
        weight = 1

!!$        do j = -1, 1
!!$           do k = -1, 1
!!$              moran = moran + sum( eoshift(eoshift(datanorm,j,dim=1),k,dim=2)*shifted, mask = shiftedmask )
!!$              ! if (debug) print *,moran
!!$              weight = weight + 1
!!$           end do
!!$        end do
!!$
!!$        moran = (real(nxny) / sumsqr * real(count(shiftedmask)) * real(weight) ) * moran

        weight = 0

        do j = 1, nx
           do k = 1, ny
              if (.not. mask(j,k)) cycle
              ! do l = -1, 1
              do l = 0, 0
                 if (j+l > nx .or. j+l < 1) cycle
                 ! do m = -1, 1
                 do m = 0, 0
                    if (k+m > ny .or. k+m < 1) cycle
                    if (.NOT. shiftedmask(j+l,k+m) ) cycle
                    moran = moran + datanorm(j,k) * shifted(k+m,j+l)
                    geary = geary + (data(j,k) - shifted(k+m,j+l))**2
                    weight = weight + 1
                 end do
              end do
           end do
        end do

        if (debug) print *,i,weight,moran

        moran = (real(nxny) / (sumsqr * real(weight) ) ) * moran

        if (debug) print *,i,moran

        if (debug) print *,i,weight,geary

        geary = 0.5 * (real(nxny - 1) / (sumsqr * real(weight) ) ) * geary

        if (debug) print *,i,geary

     end do

     deallocate(data)

  end do

contains

  real function mean (array, mask)

    integer, intent(in) :: array(:,:)
    logical, intent(in) :: mask(:,:)

    integer(kind=8) :: total, npixels

    total = 0
    npixels = 0 
    do i = 1, size(array,2)
       total = total + sum(array(:,i),mask=mask(:,i))
       npixels = npixels + count(mask(:,i))
       ! print *,i,total
    end do
    mean = real(total,8) / real(npixels,8)

  end function mean

  subroutine usage

    write(stderr,*)
    write(stderr,*) 'autocorr calculates the autocorrelation coefficient of an image file'
    write(stderr,*)
    write(stderr,*) 'Usage: autocorr [OPTIONS] mol2file [group-spec]'
    write(stderr,*)
    write(stderr,*) ' OPTIONS:'
    write(stderr,*)
    write(stderr,*) '  --help    - print this message and exit'
    write(stderr,*) '  --debug   - suppress diagnostic output'
    write(stderr,*) '  --num     - number of vectors'
    write(stderr,*)

  end subroutine usage

end program autocorr
