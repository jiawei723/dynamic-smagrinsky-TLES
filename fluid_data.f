      real rlsim ! target Taylor-microscale Reynolds number
      parameter (rlsim = 53.5)
      integer nksim ! wavenumbers
      parameter(nksim = 8) 
      real an(3,nksim), bn(3,nksim) ! Fourier vector 
      real kint(3,nksim) ! discrete wavenumber vector 
      real kdbl(nksim) ! continuous wavenumber
      real fvf ! constant volume flow forcing
      real utildeold  ! utilde of prev. time step
      real vxold(lx1,ly1,lz1,lelv),vyold(lx1,ly1,lz1,lelv),
     + vzold(lx1,ly1,lz1,lelv) ! vel. field at t^(n-3)

c     common blocks
      common /fluid_real/ an,bn,kint,kdbl,fvf,utildeold,
     + vxold,vyold,vzold
