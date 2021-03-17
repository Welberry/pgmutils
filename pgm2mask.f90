program pgm2mask

  use iso_varying_string
  use cmdline_arguments, only: get_options, bad_options, have_args, option_exists, &
      has_value, get_value, assignment(=), next_arg, num_args
  use string_functions, only: join, real
  use file_functions, only: exists, freeunit, stderr, stdout, open
  use pnm_class
  use precision

  implicit none

  integer :: i, j, num, error, nx, ny, unit, maskvalue

  character(len=2000) :: fname, outname

  type (pnm_object) :: pnm
  type (varying_string) :: myoptions(3)

  integer(i8_kind), dimension(:,:), allocatable :: total

  integer, dimension(:,:), allocatable :: data, indices
  
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

  maskvalue = 0

  do while (have_args())

     ! The files are on the command line
     fname = next_arg()

     call read(pnm,trim(fname))

     nx = size(pnm,1)
     ny = size(pnm,2)

     ! Allocate the space required for the data
     allocate(data(nx,ny),indices(nx,ny))
     
     data = pnm

     ! Reverse the order of the y-axis
     data = data(:,ny:1:-1)

     num = 0
     do j = 1, ny
        do i = 1, nx
           num = num+1
           indices(i,j) = num
        end do
     end do

     ! Use the last input filename as a template for the output filename.
     ! Search for a ".pgm" suffix .. if we find one than 
     ! replace it with "_ave.pgm", otherwise just slap "_ave.pgm" on the end
     i = index(fname,".pgm",back=.TRUE.)
     if (i == 0) then
        ! Couldn't find a matching suffix, so look for the end of the
        ! filename string in the character variable
        i = verify(fname," ",back=.TRUE.) + 1
     end if
     outname = fname
     outname(i:) = ".mask"
     print *,i,trim(outname)

     i = rle_encode(trim(outname), pack(indices, data == maskvalue))
     
     deallocate(data,indices)

  end do

contains

  subroutine usage

    write(stderr,*)
    write(stderr,*) 'Convert pgm files to mask file format'
    write(stderr,*)
    write(stderr,*) 'Usage: pgm2mask [--help] pgmfile(s)'
    write(stderr,*)
    write(stderr,*) '  --help    - print this message'
    write(stderr,*)

  end subroutine usage

  integer function rle_encode (filename, indices) result(nrecords)

    implicit none
    
    character(len=*), intent(in)  :: filename
    integer, intent(in)           :: indices(:)

    ! Local variables
    integer :: i, num, index, unitnum, masklength

    logical :: masking

    unitnum = freeunit()
    num = size(indices)

    open (unit=unitnum, file=filename, form='unformatted', status='replace')
    
    masking = .false.
    masklength = 1
    index = indices(1)
    
    write(unitnum) 'RLE'
    write(unitnum) num
    
    ! Return value
    nrecords = 0

    if (num /= 0) then
       do i = 2, num
          if ((indices(i) - indices(i-1)) == 1) then
             masklength=masklength+1
          else
             write(unitnum)index,masklength
             nrecords = nrecords + 1
             index = indices(i)
             masklength = 1
          end if
       end do

       nrecords = nrecords + 1
       ! Must write out the last record
       write(unitnum)index,masklength

    end if
    
    close(unitnum)
    
  end function rle_encode

end program pgm2mask
