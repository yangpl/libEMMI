program main
  implicit none

  integer :: ir1, ir2, is1, is2, isrc, irec
  integer :: nsrc1, nsrc2, nrec1, nrec2
  real :: xr1, xr2, xr3, xs1, xs2, xs3
  real :: drec1, drec2, dsrc1, dsrc2
  real :: xs1min, xs2min, xs1max, xs2max
  real :: xr1min, xr1max, xr2min, xr2max, tmp
  real :: xmax, xmin
  real, parameter :: PI=3.14159265

  
  !survey src range
  xs1min = -6000
  xs1max = 6000
  xs2min = -6000
  xs2max = 6000
  !number of survey in xr1 and xr2
  nsrc1 = 5
  nsrc2 = 5
  ! survey src separation dsrc1 and dsrc2
  dsrc1 = (xs1max-xs1min)/(nsrc1-1)
  dsrc2 = (xs2max-xs2min)/(nsrc2-1)

  
  !receiver range
  xr1min = -8000
  xr1max = 8000
  xr2min = -8000
  xr2max = 8000
  !number of receivers in xr1 and xr2
  nrec1 = 161
  nrec2 = 161
  ! receiver separation drec1 and drec2
  drec1 = (xr1max-xr1min)/(nrec1 - 1)
  drec2 = (xr2max-xr2min)/(nrec2 - 1)

  !------------------------------------------------
  open(10, file='src_rec_table.txt', status='replace')
  write(10,'(A12,A12)') 'isrc', 'irec'
  do is2=1,nsrc2
     do is1=1,nsrc1
        isrc = is1 + nsrc1*(is2-1)
        
        irec = 0
        do ir2=1,nsrc2
           do ir1=1,nrec1
              irec = irec + 1

              write(10,*) isrc, irec
           enddo
        enddo
        do ir2=1,nrec2
           do ir1=1,nsrc1
              irec = irec + 1

              write(10,*) isrc, irec
           enddo
        enddo
     enddo
  enddo
  close(10)

  !---------------------------------------------
  open(10, file='receivers.txt', status='replace')
  write(10,*)'            x                y                z         azimuth        dip        iRx'
  xr3 = 700
  irec = 0
  do ir2=1,nsrc2
     xr2 = xs2min + (ir2-1)*dsrc2
     do ir1=1,nrec1
        xr1 = xr1min + (ir1-1)*drec1
        irec = irec + 1

        write(10,*) xr1,xr2,xr3,0,0,irec !receiver src
     enddo
  enddo
  do ir2=1,nrec2
     xr2 = xr2min + (ir2-1)*drec2
     do ir1=1,nsrc1
        xr1 = xs1min + (ir1-1)*dsrc1
        irec = irec + 1

        write(10,*) xr1,xr2,xr3,0,0,irec !receiver src
     enddo
  enddo
  close(10)


  
  xmax = 9000
  xmin = -9000
  !-------------------------------------------
  open(10, file='sources.txt', status='replace')
  write(10,*)'            x                y                z         azimuth        dip        iTx'
  do is2=1,nsrc2
     xs2 = xs2min + (is2-1)*dsrc2
     do is1=1,nsrc1
        xs1 = xs1min + (is1-1)*dsrc1

        xs3 = 800 - 100* cos(2.*PI*xs1/(xmax-xmin)) !at seabed
        
        isrc = is1 + nsrc1*(is2-1)
        write(10,*) xs1, xs2, xs3, 0, 0, isrc !source src
     enddo
  enddo
  close(10)


end program main
