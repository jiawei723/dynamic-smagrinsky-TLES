c-----------------------------------------------------------------------
c>    initialize shared data
      subroutine shared_init(xm1,ym1,zm1)
      implicit none
c     interface
      include 'SIZE'
      include 'lib/shared_data.f'
      include 'lib/slf_data.f'
      real xm1(lx1,ly1,lz1,lelt),ym1(lx1,ly1,lz1,lelt),
     + zm1(lx1,ly1,lz1,lelt)
      real glmax
c     implementation

c     get boundaries of global and local domain
      ldmx = glmax(xm1,lx1*ly1*lz1*lelt) 
      ldmy = glmax(ym1,lx1*ly1*lz1*lelt) 
      ldmz = glmax(zm1,lx1*ly1*lz1*lelt) 

c     initialize the forcing parameter
      Af = Afr

      return
      end

