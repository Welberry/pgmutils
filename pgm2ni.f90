program pgm2ni

  use iso_varying_string
  use cmdline_arguments, only: get_options, bad_options, have_args, option_exists, &
      has_value, get_value, assignment(=), next_arg, num_args
  use string_functions, only: join, real
  use file_functions, only: exists, freeunit, stderr, stdout, open
  use pnm_class
  use precision

  implicit none

  integer :: i, error, nx, ny, unit

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

  if (num_args() < 1) then
     write(stderr,*) 'ERROR! Must supply pgm file(s) as command-line arguments'
     call usage
     STOP
  end if

  do while (have_args())

     ! The files are on the command line
     fname = next_arg()

     call read(pnm,trim(fname))

     nx = size(pnm,1)
     ny = size(pnm,2)

     ! Allocate the space required for the data
     allocate(data(nx,ny))
     
     data = pnm

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
     outname(i:) = ".ni"
     
     ! Write the nipl file.
     ! call write(pnm, outname)
     unit = open(trim(outname))
     write(unit,*) nx, ny
     write(unit,*) 0, nx-1, 0, ny-1
     do i = 1, ny
        write(unit,'(20000I8)') data(:,i)
     end do
     close(unit)

     deallocate(data)

  end do

contains

  subroutine usage

    write(stderr,*)
    write(stderr,*) 'Convert pgm files to kuplot NIPL format'
    write(stderr,*)
    write(stderr,*) 'Usage: pgm2ni [--help] pgmfile(s)'
    write(stderr,*)
    write(stderr,*) '  --help    - print this message'
    write(stderr,*)

  end subroutine usage

end program pgm2ni
