program raw2pgm

  use cmdline_arguments
  use string_functions, only: join, real
  use file_functions, only: stderr
  use pnm_class
  use precision, only: i8_kind, rd_kind 
  use binary_io, only: read, open, close, binary_filehandle
  use iso_varying_string

  implicit none

  integer :: i, j, error, nx, ny, unit, depth, tmp, tmp2

  character(len=2000) :: fname, outname

  logical :: verbose = .FALSE., autoscale, fixed_max, fixed_nx, fixed_ny, fixed_depth
  logical :: little_endian = .TRUE.

  type (pnm_object) :: pnm
  type (varying_string) :: myoptions(9)

  real :: maximum, xnorm

  integer(i8_kind), dimension(:,:), allocatable :: data
  
  type (binary_filehandle) :: fh

  ! These are our accepted command line options (see subroutine usage for
  ! an explanation)

  myoptions(1) = 'help'
  myoptions(2) = 'noscale'
  myoptions(3) = 'verbose'
  myoptions(4) = 'depth'
  myoptions(5) = 'nx'
  myoptions(6) = 'ny'
  myoptions(7) = 'norm'
  myoptions(8) = 'bigendian'

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
     write(stderr,*) 'ERROR! Must supply raw file(s) as command-line arguments'
     call usage
     STOP
  end if

  verbose = option_exists('verbose')
  little_endian = .NOT. option_exists('bigendian')

  ! See if we have specified a maximum we wish to scale the data to
  if (option_exists('norm')) then
     ! Make sure we have a value
     if (.NOT. has_value('norm')) then
        write(stderr,*) 'Option norm must have a value!'
        call usage
        stop
     end if
     ! Get the maximum value we will scale to
     xnorm = get_value('norm')
     fixed_max = .TRUE.
  else
     autoscale = .TRUE.
     fixed_max = .FALSE.
  end if

  ! See if we have specified a width
  if (option_exists('nx')) then
     ! Make sure we have a value
     if (.NOT. has_value('nx')) then
        write(stderr,*) 'Option nx must have a value!'
        call usage
        stop
     end if
     ! Get the maximum value we will scale to
     nx = get_value('nx')
     fixed_nx = .TRUE.
  else
     fixed_nx = .FALSE.
  end if

  ! See if we have specified a height
  if (option_exists('ny')) then
     ! Make sure we have a value
     if (.NOT. has_value('ny')) then
        write(stderr,*) 'Option ny must have a value!'
        call usage
        stop
     end if
     ! Get the maximum value we will scale to
     ny = get_value('ny')
     fixed_ny = .TRUE.
  else
     fixed_ny = .FALSE.
  end if

  ! See if we have specified a depth (bits per pixel)
  if (option_exists('depth')) then
     ! Make sure we have a value
     if (.NOT. has_value('depth')) then
        write(stderr,*) 'Option depth must have a value!'
        call usage
        stop
     end if
     ! Get the maximum value we will scale to
     depth = get_value('depth')
     fixed_depth = .TRUE.
  else
     depth = 2
     fixed_depth = .FALSE.
  end if

  if (fixed_nx) then
     if (.not. fixed_ny) then
     else
     end if
  else
     if (.not. fixed_ny) then
        write(stderr,*) 'Must specify at least nx or ny as command line options'
        call usage
        stop
     end if
  end if

  do while (have_args())

     ! The files are on the command line
     fname = next_arg()

     call open(fh, trim(fname))
     ! unit = open(trim(fname), status='old')
     ! read(unit,*) nx, ny

     ! Allocate the space required for the data
     allocate(data(nx,ny))
     
     do i = 1, ny
        do j = 1, nx
           call read(fh, tmp,  depth, littleendian=little_endian)
           ! call read(fh, tmp2, 2, littleendian=little_endian)
           ! data(j,i) = ishft(tmp2,16) + tmp
           data(j,i) = ishft(tmp2,16) + tmp
        end do
        ! print *,i
     end do

     call close(fh)

     ! autoscale = .false.

     if (autoscale) then
        maximum = maxval(data)
        print *,'scaling to ',maximum
        data = nint(data * (real(2**16 - 1) / maxval(data)))
     end if

     pnm = int(data)

     ! Use the last input filename as a template for the output filename.
     ! Search for a ".***" suffix .. if we find one than 
     ! replace it with ".pgm", otherwise just slap ".pgm" on the end
     i = index(fname,".")
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
    write(stderr,*) '  --help    - print this message'
    write(stderr,*) '  --verbose - verbose output'
    write(stderr,*) '  --noscale - do not rescale data to 65535'
    write(stderr,*)

  end subroutine usage

end program raw2pgm
