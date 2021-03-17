program ni2pgm

  use iso_varying_string
  use cmdline_arguments, only: get_options, bad_options, have_args, option_exists, &
      has_value, get_value, assignment(=), next_arg, num_args
  use string_functions, only: join, real
  use file_functions
  use pnm_class
  use precision
  use image_transforms
  use fundamental_constants, only: radian

  implicit none

  integer :: i, error, nx, ny, unit, nfold

  character(len=2000) :: fname, outname

  logical :: verbose = .FALSE., autoscale, twofold
  logical :: rotaverage, hmirror, vmirror, fixed_max

  type (pnm_object) :: pnm
  type (varying_string) :: myoptions(8)

  real :: maximum, xnorm

  real, dimension(:,:), allocatable :: data

  real, parameter :: no_data_value = -9999
  
  ! These are our accepted command line options (see subroutine usage for
  ! an explanation)

  myoptions(1) = 'help'
  myoptions(2) = 'noscale'
  myoptions(3) = 'verbose'
  myoptions(4) = 'rotave'
  myoptions(5) = 'twofold'
  myoptions(6) = 'norm'
  myoptions(7) = 'hmirror'
  myoptions(8) = 'vmirror'

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
     write(stderr,*) 'ERROR! Must supply nipl file(s) as command-line arguments'
     call usage
     STOP
  end if

  autoscale = .NOT. option_exists('noscale')
  verbose = option_exists('verbose')
  twofold = option_exists('twofold')
  hmirror = option_exists('hmirror')
  vmirror = option_exists('vmirror')

  if (option_exists('rotave')) then
     if (twofold) then 
        write(stderr,*) 'Requested rotave and twofold simultaneously. Ignoring twofold option.'
        twofold = .FALSE.
     end if
     if (.NOT. has_value('rotave')) then
        write(stderr,*) 'Option rotave must have a value!'
        call usage
        stop
     end if
     ! Get the n-fold rotation ...
     nfold = get_value('rotave')
     rotaverage = .TRUE.
  end if

  ! See if we have specified a maximum we wish to scale the data to
  if (option_exists('norm')) then
     autoscale = .FALSE.
     ! Make sure we have a value
     if (.NOT. has_value('norm')) then
        write(stderr,*) 'Option norm must have a value!'
        call usage
        stop
     end if
     ! Get the maximum value we will scale to
     xnorm = get_value('norm')
     fixed_max = .TRUE.
  end if

  do while (have_args())

     ! The files are on the command line
     fname = next_arg()

     unit = open(trim(fname), status='old')
     read(unit,*) nx, ny

     ! Allocate the space required for the data
     allocate(data(nx,ny))
     
     read(unit,*)

     do i = 1, ny
        read(unit,*) data(:,i)
     end do
     close(unit)

     if (verbose) then
        print *,'Read in ',trim(fname)
        print *,'Max value =  ',maxval(data)
        print *,'Min value =  ',minval(data)
     end if

     ! 2-fold average the data ? (This is hokey two-fold rotation averaging
     ! where we run the x and y directions in reverse order, i.e. a flip
     ! and a flop)
     if (twofold) call average(data,data(nx:1:-1,ny:1:-1))
     
     ! n-fold average the data ?
     if (rotaverage) call rotation_average(data, nfold)
     
     ! Apply horizontal mirror?
     if (hmirror) call average(data,data(nx:1:-1,:))
     
     ! Apply vertical mirror?
     if (vmirror) call average(data,data(:,ny:1:-1))
        
     if (autoscale .or. fixed_max) then
        if (autoscale) maximum = maxval(data)
        if (fixed_max) maximum = xnorm
        data = data * (real(2**16 - 1) / maximum) 
        if (verbose) then
           print *,'Normalising to ',maximum
        end if
     end if

     pnm = nint(data)

     ! Use the last input filename as a template for the output filename.
     ! Search for a ".ni" suffix .. if we find one than 
     ! replace it with ".pgm", otherwise just slap ".pgm" on the end
     i = index(fname,".ni")
     if (i == 0) then
        ! Couldn't find a matching suffix, so look for the end of the
        ! filename string in the character variable
        i = verify(fname," ",back=.TRUE.) + 1
     end if
     outname = fname
     outname(i:) = ".pgm"
     
     if (verbose) print *,'Wrote ',trim(outname)

     ! Write the pgm file
     call write(pnm, trim(outname))

     deallocate(data)

  end do

contains

  subroutine usage

    write(stderr,*)
    write(stderr,*) 'Convert kuplot NIPL formatted files to pgm format'
    write(stderr,*)
    write(stderr,*) 'Usage: ni2pgm [--help] [--verbose] [--noscale] nifile(s)'
    write(stderr,*)
    write(stderr,*) '  --help     - print this message'
    write(stderr,*) '  --verbose  - verbose output'
    write(stderr,*) '  --noscale  - do not rescale data to 65535'
    write(stderr,*) '  --norm     - scale the data to this value'
    write(stderr,*) '  --twofold  - two-fold average the data.'
    write(stderr,*) '  --rotave=n - n-fold average the data (n is any integer).'
    write(stderr,*) '  --hmirror  - apply horizonal mirror averaging'
    write(stderr,*) '  --vmirror  - apply vertical mirror averaging'
    write(stderr,*)

  end subroutine usage

  subroutine rotation_average(input, nfold)

    ! Interface variables
    real, dimension(:,:), intent(inout) :: input
    integer, intent(in)                 :: nfold

    ! Local variables
    real, dimension(size(input,1),size(input,2))    :: dsitot, dsibuffer
    integer, dimension(size(input,1),size(input,2)) :: bufcount, totcount
    real(kind=rd_kind) :: step
    integer :: i, j, k

    dsitot = input
    totcount = 0
    where(dsitot /= no_data_value) totcount = 1
    
    step = (360.d0/real(nfold,rd_kind))/radian
    
    do i = 1, nfold-1
       ! print *,'number ',i,real(i,rd_kind)*step*radian
       dsibuffer = input
       call rotate_image(dsibuffer,real(i,rd_kind)*step,no_data_value)
       do j = 1, size(input,1)
          do k = 1, size(input,2)
             if (dsibuffer(j,k) /= no_data_value) then
                if (dsitot(j,k) /= no_data_value) then
                   dsitot(j,k) = dsitot(j,k) + dsibuffer(j,k)
                else
                   dsitot(j,k) = dsibuffer(j,k)
                end if
                totcount(j,k) = totcount(j,k) + 1
             end if
          end do
       end do
    end do

    where(totcount > 0) dsitot = dsitot/real(totcount)
    input = dsitot

  end subroutine rotation_average

  subroutine average (original, transform)

    ! Interface variables
    real, dimension(:,:), intent(inout) :: original
    real, dimension(:,:), intent(in)    :: transform

    ! Local variables
    real, dimension(size(original,1),size(original,2)) :: flipflopcopy

    flipflopcopy = transform

    where (flipflopcopy == no_data_value)
       ! we have no data in our copy, get data from the original
       flipflopcopy = original
    end where
           
    where (original == no_data_value)
       ! we have no data in our original, get data from the copy
       original = flipflopcopy
    elsewhere
       ! we have data in both now -- take average
       original = (original + flipflopcopy) / 2.
    end where

  end subroutine average
    
end program ni2pgm
