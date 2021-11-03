c     low wave number forcing      
      integer NU                  ! number of grid points per direction
      parameter (NU = 58)         ! important: NU has to be even
      integer plan(8)             ! FFTW
      real TL                     ! forcing time scale
      parameter (TL = 1.7)
      real sigma                  ! "std deviation" of rand. walk
      parameter (sigma = 0.282)
      real udx                    ! uniform spacing of comp. grid
      real k0                     ! zero wave number
      real kF                     ! forcing sphere radius
      real abv(NU,NU,NU)          ! squared abs value of wavenumber vec.
      real flwx(lx1,ly1,lz1,lelv) ! forcing in x direction
      real flwy(lx1,ly1,lz1,lelv) ! forcing in y direction
      real flwz(lx1,ly1,lz1,lelv) ! forcing in z direction
c     at previous time steps (used in TLES)
      real flwx_old(lx1,ly1,lz1,lelv)
      real flwy_old(lx1,ly1,lz1,lelv)
      real flwz_old(lx1,ly1,lz1,lelv)
c     wavenumber vector components
      real wnvx(NU,NU,NU),wnvy(NU,NU,NU),wnvz(NU,NU,NU)             
c     random walk (UO) components
      double complex bx(NU,NU,NU),by(NU,NU,NU),bz(NU,NU,NU)
c     components of forcing acceleration aF
      real aFxr(NU,NU,NU),aFyr(NU,NU,NU),aFzr(NU,NU,NU)

c     point source to initiate turbulent transition
      real fpsx(lx1,ly1,lz1,lelv) ! forcing in x direction
      real fpsy(lx1,ly1,lz1,lelv) ! forcing in y direction
      real fpsz(lx1,ly1,lz1,lelv) ! forcing in z direction 
     
c     common blocks
      common /forcing_real/ udx,k0,kF,abv,flwx,flwy,flwz,
     + wnvx,wnvy,wnvz,aFxr,aFyr,aFzr
      common /forcing_cmplx/ bx,by,bz
      common /ps_real/ fpsx,fpsy,fpsz
