program braggextract

  use cluster_functions
  use iso_varying_string
  use cmdline_arguments, only: get_options, bad_options, have_args, option_exists, &
      has_value, get_value, assignment(=), next_arg, num_args
  use string_functions, only: join, real
  use file_functions
  use sort_functions, only: sort, unique
  use pnm_class
  use precision
  use fundamental_constants, only: pi
  use variable_array
  use image_transforms, only: rotate_image

  implicit none

  integer :: i, j, error, nx, ny, unit, x, y, largest_peak, nclusters

  character(len=2000) :: fname, outname

  logical :: initialised = .FALSE., fixed_name = .FALSE., have_missing = .FALSE.

  type (pnm_object) :: pnm, newpnm
  type (varying_string) :: myoptions(7)

  integer(i8_kind), dimension(:,:), allocatable :: total

  integer, dimension(:,:), allocatable :: data, number, peakdata, clusters
  real, dimension(:,:), allocatable :: smoothbackground, smoothbacksquared, smoothstddev
  logical, dimension(:,:), allocatable :: data_mask
  integer, dimension(:), allocatable :: indices, onedclusters
  real, dimension(:), allocatable :: xcentres, ycentres
  integer, dimension(:), pointer :: clustertmp

  real :: sigma, braggthreshold, datamean, datastddev, maskvalue
  real :: xcentre, ycentre, radius
  integer :: nxny, npixels, x0, x1, y0, y1, n_xc, n_yc

  logical :: fixed_sigma, fixed_braggthreshold, verbose, fixed_radius, centre_max
  
  ! These are our accepted command line options (see subroutine usage for
  ! an explanation)

  myoptions(1) = 'help'
  myoptions(2) = 'sigma'
  myoptions(3) = 'braggthreshold'
  myoptions(4) = 'verbose'
  myoptions(5) = 'radius'
  myoptions(6) = 'maskimage'
  myoptions(7) = 'centremax'

  ! This call parses the command line arguments for command line options
  call get_options(myoptions, error)

  ! Check we weren't passed duff options -- spit the dummy if we were
  if (error > 0) then
     write(stderr,*) 'ERROR! Unknown options: ',join(bad_options()," ")
     call usage
     STOP
  end if

  ! Check if we just want to print the usage
  if (option_exists('help')) then
     call usage
     STOP
  end if

  verbose = option_exists('verbose')
  centre_max = option_exists('centremax')

  if (num_args() < 1) then
     write(stderr,*) 'ERROR! Must supply pgm file(s) as command-line arguments'
     call usage
     STOP
  end if

  ! See if we have specified a sigma. If I/sigma(I) > sigma we 
  ! will treat the pixel as being part of a peak.
  if (option_exists('sigma')) then
     ! Make sure we have a value
     if (.NOT. has_value('sigma')) then
        write(stderr,*) 'Option sigma must have a value!'
        call usage
        stop
     end if
     ! Get the maximum value we will scale to
     sigma = get_value('sigma')
     fixed_sigma = .TRUE.
  else
     sigma = 5.
  end if

  if (option_exists('braggthreshold')) then
     ! Make sure we have a value
     if (.NOT. has_value('braggthreshold')) then
        write(stderr,*) 'Option braggthreshold must have a value!'
        call usage
        stop
     end if
     ! Get the maximum value we will scale to
     braggthreshold = get_value('braggthreshold')
     fixed_braggthreshold = .TRUE.
  else
     braggthreshold = 0.5
  end if

  if (option_exists('radius')) then
     ! Make sure we have a value
     if (.NOT. has_value('radius')) then
        write(stderr,*) 'Option radius must have a value!'
        call usage
        stop
     end if
     ! Get the maximum value we will scale to
     radius = get_value('radius')
     radius = radius - 1.
     fixed_radius = .TRUE.
  end if

  do while (have_args())

     ! The files are on the command line
     fname = next_arg()

     call read(pnm,trim(fname))

     nx = size(pnm,1)
     ny = size(pnm,2)

     ! Allocate the space required for the data
     allocate(data(nx,ny),peakdata(nx,ny),smoothbackground(nx,ny),smoothbacksquared(nx,ny),smoothstddev(nx,ny))
     allocate(clusters(nx,ny))
     allocate(data_mask(nx,ny))
     
     ! Get the data out of the pnm object
     data = pnm

     ! We make the assumption that pixels with no counts are
     ! outside the detector range and can be regarded as having
     ! no data (this is mostly the area round the outside of the
     ! circular image plate)
     data_mask = (data /= 0)

     ! Calculate the mean
     datamean = mean(data,data_mask)

     ! Calculate the standard deviation (use N-1 population). Don't 
     ! use the 'mean' routine as the squared data gets 'clipped' by
     ! the limited precision and gives a bogus answer. Instead we
     ! promote to high precision and use inplace sum routine.
     datastddev = sqrt(real(sum(int(data,8)**2,data_mask),8)/real(count(data_mask)-1,8) - datamean**2)

     if (verbose) print *,'Mean = ', datamean
     if (verbose) print *,'stddev = ', datastddev

     smoothbackground = 0.
     smoothbackground = rolling_ball(int(data,8), data_mask)
     smoothbacksquared = rolling_ball(int(data,8)**2, data_mask)
     smoothstddev = sqrt(smoothbacksquared/((nx+1)*(ny+1)/(50.**2)) - smoothbackground/((nx+1)*(ny+1)/(50.**2)))
     
     if (verbose) then 
        newpnm = nint(sqrt(smoothbacksquared/((nx+1)*(ny+1)/(50.**2)) - smoothbackground/((nx+1)*(ny+1)/(50.**2))))
        call write(newpnm,'localvariance.pgm')
     end if

     if (verbose) print *,'Initial number of true values: ',count((data >= smoothbackground + smoothstddev*sigma) .and. data_mask)

     ! This function from the cluster_functions module will return an
     ! integer array where all the pixels which are 'true' in the input
     ! are labelled by cluster number. In this case all the true pixels
     ! will be those which exceed our mask value
     clusters = cluster_image( (data >= smoothbackground + smoothstddev*sigma) .and. data_mask )

     ! Make a one dimensional version of the cluster image
     allocate(onedclusters(nx*ny))
     allocate(indices(nx*ny))

     if (verbose) then
        print *,'Setting up onedclusters ...'
        pnm = clusters
        call write(pnm,'clusters.pgm')
     end if

     onedclusters = pack(clusters,mask=.TRUE.)

     if (verbose) print *,'Setting up indices ...'

     ! The much tidier inplace array notation was horrible slow for
     ! values of nxny >= 1000000. Use loop instead.
     ! indices = (/ (i, i=1,nxny) /)
     do i = 1, nx*ny
        indices(i) = i
     end do

     ! The clusters are numbered from 1 .. n, so the maximum value will
     ! be the number of different clusters
     nclusters = maxval(clusters)

     allocate(xcentres(nclusters), ycentres(nclusters))

     peakdata = 0

     do i = 1, nclusters

        ! We make a temporary array which is just the indices of the
        ! cluster we are interested in (the push function returns the
        ! size of the array after the push, which is the number of 
        ! pixels in our cluster)
        npixels = push(clustertmp, pack(indices, mask=(onedclusters==i)))

        ! print *,clustertmp
        ! print *,pack(onedclusters, mask=(onedclusters==i))
        ! stop

        ! Find the centre of the cluster
        if (centre_max) then
           largest_peak = -1
           do j = 1, npixels
              x = mod(clustertmp(j)-1,nx)+1
              y = ((clustertmp(j)-1)/ny)+1
              if (data(x,y) > largest_peak) then
                 ! print *,i,j,clustertmp(j),x,y,largest_peak,data(x,y)
                 largest_peak = data(x,y)
                 xcentre = x
                 ycentre = y
              end if
           end do
        else
           xcentre = real(sum(mod(clustertmp-1,nx)+1))/real(npixels)
           ycentre = real(sum(((clustertmp-1)/ny)+1))/real(npixels)
        end if

        if (verbose) print '("Cluster ",I4," is centred at x=",F0.2," y=",F0.2," and contains ",I5," pixels")'&
             ,i,xcentre,ycentre,npixels

        ! Determine the radius of the peak
        if (.not. fixed_radius) radius = sqrt(real(npixels) / pi)

        n_xc = nint(xcentre)
        n_yc = nint(ycentre)
        x0 = n_xc-nint(radius)
        x1 = n_xc+nint(radius)
        y0 = n_yc-nint(radius)
        y1 = n_yc+nint(radius)
        peakdata(x0:x1,y0:y1) = circular_mask(peakdata(x0:x1,y0:y1),data(n_xc, n_yc))

        ! xcentres(i) = n_xc
        ! ycentres(i) = n_yc
        xcentres(i) = xcentre
        ycentres(i) = ycentre

        ! Remove all the elements from the clustertmp array
        npixels = splice(clustertmp, 0)
        
     end do

     call sort(xcentres)
     call sort(ycentres)

     print *,unique(xcentres)
     print *,unique(ycentres)
     
     ! Now we mask off the data ...
     where (clusters > 0) data = 0

     ! Use the last input filename as a template for the output filename.
     ! Search for a ".pgm" suffix .. if we find one than 
     ! replace it with "_bragg.pgm", otherwise just slap "_bragg.pgm" on the end
     i = index(fname,".pgm")
     if (i == 0) then
        ! Couldn't find a matching suffix, so look for the end of the
        ! filename string in the character variable
        i = verify(fname," ",back=.TRUE.) + 1
     end if
     outname = fname
     outname(i:) = "_bragg.pgm"

     pnm = peakdata

     if (verbose) print *,'Writing ',trim(outname)
     
     call write(pnm, outname)

     if (option_exists('maskimage')) then
        ! Make sure we have a value
        if (.NOT. has_value('maskimage')) then
           write(stderr,*) 'Option maskimage must have a value!'
           call usage
           stop
        end if
        ! Get the maximum value we will scale to
        fname = get_value('maskimage')

        call read(pnm, trim(fname))

        call rotate_image(as_array_2d(pnm), 0., data, 0)

        where (peakdata /= 0) data = peakdata

        i = index(fname,".pgm")
        if (i == 0) then
           ! Couldn't find a matching suffix, so look for the end of the
           ! filename string in the character variable
           i = verify(fname," ",back=.TRUE.) + 1
        end if
        outname = fname
        outname(i:) = "_bragg.pgm"

        pnm = data
        call write(pnm, trim(outname))

     end if

     deallocate(data)
     deallocate(peakdata)
     deallocate(smoothbackground)
     deallocate(smoothbacksquared)
     deallocate(smoothstddev)
     deallocate(clusters)
     deallocate(data_mask)
     deallocate(onedclusters)
     deallocate(indices)
     deallocate(xcentres)
     deallocate(ycentres)
     
  end do

contains

  subroutine usage

    write(stderr,*)
    write(stderr,*) 'Convert pgm files to kuplot NIPL format'
    write(stderr,*)
    write(stderr,*) 'Usage: maskpeaks [--help] pgmfile(s)'
    write(stderr,*)
    write(stderr,*) '  --help    - print this message'
    write(stderr,*)

  end subroutine usage

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
    if (npixels > 0) then
       ! print *,total,npixels
       mean = real(total,8) / real(npixels,8)
    else
       mean = -9999
    end if

  end function mean

  real function var (array, mask)

    integer, intent(in) :: array(:,:)
    logical, intent(in) :: mask(:,:)

    integer(kind=8) :: total, npixels
    integer(kind=8), dimension(size(array,1),size(array,2)) :: array8

    total = 0
    npixels = 0 

    array8 = array**2

    do i = 1, size(array,2)
       total = total + sum(array8(:,i),mask=mask(:,i))
       npixels = npixels + count(mask(:,i))
       ! print *,i,total
    end do
    if (npixels > 0) then
       ! print *,total,npixels
       var = real(total,8) / real(npixels,8) - (mean(array,mask))**2
    else
       var = -9999
    end if

  end function var

!!$  function block_background (array,mask) result(arrout)
!!$
!!$    real, dimension(size(array,1),size(array,2)) :: arrout
!!$    integer, intent(in) :: array(:,:)
!!$
!!$    integer :: x0, x1, y0, y1, xinc, yinc
!!$
!!$    yinc = nint(size(array,2)/10.)
!!$
!!$    do i = 1, size(array,2) - yinc
!!$    end do
!!$
!!$  end function block_background

  function rolling_ball (array,mask) result(arrout)

    integer(kind=8), intent(in) :: array(:,:)
    logical, intent(in) :: mask(:,:)

    real, dimension(size(array,1),size(array,2)) :: arrout

    integer :: x0, x1, y0, y1, dx, dy, nx, ny

    nx = size(array,1)
    ny = size(array,1)

    dx = nint(nx/50.)+1
    dy = nint(ny/50.)+1

    x0 = 1
    x1 = dx
    y0 = 1
    y1 = dy

    arrout = real(array)

    yloop: do 
       xloop: do 
          ! print *,x0,x1,y0,y1
          if ((count(mask(x0:x1,y0:y1)) /= 0) .and. mask(x0+(dx/2),y0+(dy/2)) ) then
             arrout(x0+(dx/2),y0+(dy/2)) = sum(array(x0:x1,y0:y1),mask=mask(x0:x1,y0:y1))/real(count(mask(x0:x1,y0:y1)))
             ! print *,x0+dx/2,y0+dy/2,mean(array(x0:x1,y0:y1),mask(x0:x1,y0:y1))
          end if
          x0 = x0 + 1
          x1 = x1 + 1
          if (x1 > nx) exit xloop
       end do xloop
       x0 = 1
       x1 = dx
       y0 = y0 + 1
       y1 = y1 + 1
       if (y1 > ny) exit yloop
    end do yloop

  end function rolling_ball

  function circular_mask (square, mask_value) result(circle)
      
    ! A small square section of an image is the input
    integer, dimension(:,:), intent(in) :: square
    integer, intent(in)                 :: mask_value
    
    ! And a circular masked version of this is the output
    integer, dimension(size(square,1),size(square,2)) :: circle
    
    integer :: i, nx, ny, radius
       
    nx = size(square,1)
    ny = size(square,2)
    
    radius = (nx-1)/2
    
    circle = (nint(nx/2.) - spread((/ (i, i=1,nx) /),2,ny))**2 + (nint(ny/2.) - spread((/ (i, i=1,ny) /),1,ny) )**2

    ! totalmasked = totalmasked + count(circle > radius**2)
       
    where (circle <= radius**2) 
       circle = mask_value
    elsewhere
       circle = square
    end where
    
    ! do i = ny, 1, -1
    !    print '(100I5:)',distance(:,i)
    ! end do
    
  end function circular_mask

end program braggextract
