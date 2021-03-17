program pgmave

  use iso_varying_string
  use cmdline_arguments, only: get_options, bad_options, have_args, option_exists, &
      has_value, get_value, assignment(=), next_arg, num_args
  use string_functions, only: join, real
  use file_functions, only: exists, freeunit, stderr, stdout
  use pnm_class
  use precision

  implicit none

  integer :: i, error, nx, ny, missing_value

  character(len=2000) :: fname, outname

  logical :: initialised = .FALSE., fixed_name = .FALSE., have_missing = .FALSE.

  type (pnm_object) :: pnm
  type (varying_string) :: myoptions(3)

  integer(i8_kind), dimension(:,:), allocatable :: total

  integer, dimension(:,:), allocatable :: data, number
  
  ! These are our accepted command line options (see subroutine usage for
  ! an explanation)

  myoptions(1) = 'help'
  myoptions(2) = 'outfile'
  myoptions(3) = 'missing'

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

  if (option_exists('outfile')) then
     ! We have specified an output filename as command line option
     ! Make sure we have a value
     if (.NOT. has_value('outfile')) then
        write(stderr,*) 'Option outfile must specify a filename!'
        call usage
        stop
     end if
     outname = ""
     outname = get_value('outfile')
     fixed_name = .TRUE.
  end if

  ! See if we have specified a maximum we wish to scale the data to
  if (option_exists('missing')) then
     ! Make sure we have a value
     if (.NOT. has_value('missing')) then
        write(stderr,*) 'Option missing must have a value!'
        call usage
        stop
     end if
     ! Get the maximum value we will scale to
     missing_value = get_value('missing')
     have_missing = .TRUE.
  end if

  if (num_args() < 1) then
     write(stderr,*) 'ERROR! Must supply pgm file(s) as command-line arguments'
     call usage
     STOP
  end if

  do while (have_args())

     ! The files are on the command line
     fname = next_arg()

     call read(pnm,trim(fname))

     if (.NOT. initialised) then

        nx = size(pnm,1)
        ny = size(pnm,2)

        ! Allocate the space required for the data
        allocate(data(nx,ny), total(nx,ny), number(nx,ny))

        initialised = .TRUE.

        number = 0

     end if

     data = pnm

     total = total + data

     if (have_missing) then
        where(data /= missing_value) number = number + 1
     else
        number = number + 1
     end if

  end do

  total = total / real(number)

  if (maxval(total) > (2**16 - 1)) then
     total = (real(total)/real(maxval(total))) * (2**16 - 1)
  end if

  pnm = int(total)

  if (.NOT. fixed_name) then
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
     outname(i:) = "_ave.pgm"
  end if
     
  ! Write the pgm file.
  call write(pnm, outname)

  deallocate(data,total)

contains

  subroutine usage

    write(stderr,*)
    write(stderr,*) 'Usage: pgmave --help --outfile=<filename>'
    write(stderr,*)
    write(stderr,*) '  --help    - print this message.'
    write(stderr,*) '  --outfile - output file (default is replace ".pgm" with "_ave.pgm").'
    write(stderr,*) '  --missing - missing value (no data) which will be ignored when averaging.'
    write(stderr,*)

  end subroutine usage

end program pgmave
