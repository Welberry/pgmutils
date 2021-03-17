program pgmcombine
  
  use iso_varying_string
  use cmdline_arguments
  use string_functions, only: join, real
  use pnm_class
  use file_functions, only: stderr

  implicit none

  type (pnm_object), dimension(:), allocatable :: pnm_array

  type (pnm_object) :: pnm

  integer, dimension(:,:), allocatable :: data, datatmp, count

  type (varying_string) :: myoptions(3)

  integer :: error, maxx, maxy, i

  character(len=400) :: fname, outfile

  logical :: verbose

  ! These are our accepted command line options (see subroutine usage for
  ! an explanation)

  myoptions(1) = 'help'
  myoptions(2) = 'outfile'
  myoptions(3) = 'verbose'

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

  ! Check if we have specified an output file name
  if (option_exists('outfile')) then
     if (.not. has_value('outfile')) then
       write (stderr,*) 'Option outfile missing a filename'
       call usage
       STOP
     end if
     ! Get the horizontal shift
     outfile = get_value('outfile')
  else
     outfile = 'combined.pgm'
  end if

  if (num_args() < 2) then
     write(stderr,*) 'ERROR! Must supply two or more pgm files to combine'
     call usage
     STOP
  end if

  verbose = option_exists('verbose')

  allocate(pnm_array(num_args()))

  i = 0

  do while (have_args())

     i = i + 1

     ! The files are on the command line
     fname = next_arg()

     if (verbose) print *,'Reading ',trim(fname)

     ! Read in a pgm
     call read(pnm_array(i),trim(fname))

     maxx = max(maxx,size(pnm_array(i),1))
     maxy = max(maxy,size(pnm_array(i),1))

  end do

  if (verbose) write(*,'(A,I0,A,I0)') 'Combined image will be ',maxx,' x ',maxy

  allocate(data(maxx,maxy),datatmp(maxx, maxy),count(maxx,maxy))

  do i = 1, size(pnm_array)

     datatmp = 0

     ! Grab the data from the pnm object
     datatmp = pnm_array(i)


     ! Move it about a bit if it is smaller than our maximum area.
     ! Note that this is integer arithmetic, so should produce no
     ! move for an image off by only one pixel

     datatmp = eoshift(eoshift(datatmp, (maxx-size(pnm_array(i),1))/2, dim=1), (maxy-size(pnm_array(i),2))/2, dim=2)

     ! Add it to our combined output
     where (datatmp /= 0)
        data = data + datatmp
        count = count + 1
     end where

     ! Write it out

  end do

  where (count /= 0) data = nint(real(data)/real(count))

  pnm = data

  call write(pnm,trim(outfile))

  deallocate(pnm_array,data)

contains

  subroutine usage

    write(stderr,*)
    write(stderr,*) 'Usage: pgmcombine [OPTIONS] pgmfile1 pgmfile2 ...'
    write(stderr,*)
    write(stderr,*) '  --help             - print this message.'
    write(stderr,*) '  --outfile=filename - specify an output file (default: combine.pgm ).'
    write(stderr,*)
    
  end subroutine usage

end program pgmcombine
