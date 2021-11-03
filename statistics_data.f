c     velocity components averaged in the z-direction
      real ux_avgz(lx1,ly1,lz1,lelv)
      real uy_avgz(lx1,ly1,lz1,lelv)
      real uz_avgz(lx1,ly1,lz1,lelv)
      real uxdconv_avgz(lx1,ly1,lz1,lelv)
      real uydconv_avgz(lx1,ly1,lz1,lelv)
      real uzdconv_avgz(lx1,ly1,lz1,lelv)

c     interpolation integer
      integer igs_z
      integer nprobes
      integer maxsamples
      parameter (maxsamples = 1200)

c     sample points and profiles
      real probe_cords(3,maxsamples)
      real u_profile(3,maxsamples)
      real udconv_profile(3,maxsamples)

c     common block
      common /statistics_int/ igs_z,nprobes
      common /statistics_real/ ux_avgz,uy_avgz,uz_avgz,uxdconv_avgz,
     + uydconv_avgz,uzdconv_avgz,probe_cords,u_profile,udconv_profile
