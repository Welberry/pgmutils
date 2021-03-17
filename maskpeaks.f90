program maskpeaks

  use cluster_functions
  use iso_varying_string
  use cmdline_arguments, only: get_options, bad_options, have_args, option_exists, &
      has_value, get_value, assignment(=), next_arg, num_args
  use string_functions, only: join, real
  use file_functions
  use pnm_class
  use precision
  use fundamental_constants, only: pi
  use variable_array

  implicit none

  integer :: i, error, nx, ny, unit

  character(len=2000) :: fname, outname

  logical :: initialised = .FALSE., fixed_name = .FALSE., have_missing = .FALSE.

  type (pnm_object) :: pnm
  type (varying_string) :: myoptions(3)

  integer(i8_kind), dimension(:,:), allocatable :: total

  integer, dimension(:,:), allocatable :: data, number, clusters
  logical, dimension(:,:), allocatable :: data_mask
  integer, dimension(:), allocatable :: indices, onedclusters
  integer, dimension(:), pointer :: clustertmp

  real :: peaksigma, braggthreshold, datamean, datastddev, maskvalue
  real :: xcentre, ycentre, radius
  integer :: nxny, npixels

  logical :: fixed_peaksigma, fixed_braggthreshold, verbose
  
  ! These are our accepted command line options (see subroutine usage for
  ! an explanation)

  myoptions(1) = 'help'
  myoptions(2) = 'peaksigma'
  myoptions(3) = 'braggthreshold'

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

  if (num_args() < 1) then
     write(stderr,*) 'ERROR! Must supply pgm file(s) as command-line arguments'
     call usage
     STOP
  end if

  ! See if we have specified a peaksigma. If I/sigma(I) > peaksigma we 
  ! will treat the pixel as being part of a peak.
  if (option_exists('peaksigma')) then
     ! Make sure we have a value
     if (.NOT. has_value('peaksigma')) then
        write(stderr,*) 'Option peaksigma must have a value!'
        call usage
        stop
     end if
     ! Get the maximum value we will scale to
     peaksigma = get_value('peaksigma')
     fixed_peaksigma = .TRUE.
  else
     peaksigma = 3.
  end if

  ! See if we have specified a peaksigma. If I/sigma(I) > peaksigma we 
  ! will treat the pixel as being part of a peak.
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

  do while (have_args())

     ! The files are on the command line
     fname = next_arg()

     call read(pnm,trim(fname))

     nx = size(pnm,1)
     ny = size(pnm,2)

     ! Allocate the space required for the data
     allocate(data(nx,ny))
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
     datastddev = sqrt(mean(data**2,data_mask) - datamean**2)

     print *,'Mean = ', datamean
     print *,'stddev = ', datastddev

     stop

     if (verbose) print *,'Initial number of true values: ',count(data >= maskvalue)

     ! This function from the cluster_functions module will return an
     ! integer array where all the pixels which are 'true' in the input
     ! are labelled by cluster number. In this case all the true pixels
     ! will be those which exceed our mask value
     clusters = cluster_image( data >= maskvalue )

     ! Make a one dimensional version of the cluster image
     allocate(onedclusters(nx*ny))
     allocate(indices(nx*ny))

     if (verbose) print *,'Setting up onedclusters ...'

     onedclusters = pack(clusters,mask=.TRUE.)

     if (verbose) print *,'Setting up indices ...'

     ! The much tidier inplace array notation was horrible slow for
     ! values of nxny >= 1000000. Use loop instead.
     ! indices = (/ (i, i=1,nxny) /)
     do i = 1, nxny
        indices(i) = i
     end do

     if (verbose) print *,'Masking data ...'
     
     ! The clusters are numbered from 1 .. n, so the maximum value will
     ! be the number of different clusters
     do i = 1, maxval(clusters)

        ! Remove all the elements from the clustertmp array
        npixels = splice(clustertmp, 0)
        
        ! We make a temporary array which is just the indices of the
        ! cluster we are interested in (the push function returns the
        ! size of the array after the push, which is the number of 
        ! pixels in our cluster)
        npixels = push(clustertmp, pack(indices, mask=(onedclusters==i)))

        ! Find the centre of the cluster
        xcentre = real(sum(mod(clustertmp,nx)))/real(npixels)
        ycentre = real(sum((clustertmp/ny)+1))/real(npixels)

        if (verbose) print '("Cluster ",I4," is centred at x=",F0.2," y=",F0.2," and contains ",I5," pixels")'&
             ,i,xcentre,ycentre,npixels

        ! Determine the radius of overexposed peak
        radius = sqrt(real(npixels) / pi)

     end do
     
     ! Now we mask off the data ...
     where (clusters > 0) data = 0

     ! Use the last input filename as a template for the output filename.
     ! Search for a ".pgm" suffix .. if we find one than 
     ! replace it with "_ave.pgm", otherwise just slap "_ave.pgm" on the end
     i = index(fname,".pgm")
     if (i == 0) then
        ! Couldn't find a matching suffix, so look for the end of the
        ! filename string in the character variable
        i = verify(fname," ",back=.TRUE.) + 1
     end if
     outname = fname
     outname(i:) = "_masked.pgm"
     
     ! call write(pnm, outname)
     unit = open(trim(outname))

     deallocate(data)

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
    mean = real(total,8) / real(npixels,8)

  end function mean

end program maskpeaks
