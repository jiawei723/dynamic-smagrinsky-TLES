      real tij(lx1,ly1,lz1,lelv,3,3) ! stress components <uiuj>-<ui><uj>
      real tij_old(lx1,ly1,lz1,lelv,3,3) ! of previous time step
      real tij_lap(lx1,ly1,lz1,lelv,3,3) ! laplacian of tij
      real tij_xdiv(lx1,ly1,lz1,lelv) ! d(tau_1j)/(dxj)
      real tij_ydiv(lx1,ly1,lz1,lelv) ! d(tau_2j)/(dxj)
      real tij_zdiv(lx1,ly1,lz1,lelv) ! d(tau_3j)/(dxj)
      real uxdconv(lx1,ly1,lz1,lelv),uydconv(lx1,ly1,lz1,lelv),
     + uzdconv(lx1,ly1,lz1,lelv) ! deconvoluted velocity
      real dvdt(lx1,ly1,lz1,lelv,3) ! filt. velocity time derivative
      real dvdt_old(lx1,ly1,lz1,lelv,3) ! time der. of prev. time step
      real vxtemp(lx1,ly1,lz1,lelv),vytemp(lx1,ly1,lz1,lelv),
     + vztemp(lx1,ly1,lz1,lelv) ! velocity of the first time step
      real FTS ! filter time scale T
      real qbar(3,lx1,ly1,lz1,lelv) ! filtered deconv. vel. 2ndreg
      real FTS2 ! filter time scale of the secondary regularization
      real freg(3,lx1,ly1,lz1,lelv) ! regularization term

c     velocity fields of 3rd order approx. deconvolution
      real u2(3,lx1,ly1,lz1,lelv)
      real u3(3,lx1,ly1,lz1,lelv)
      real u4(3,lx1,ly1,lz1,lelv)
      real u5(3,lx1,ly1,lz1,lelv)
      real wxuf(lx1,ly1,lz1,lelv)
      real wyuf(lx1,ly1,lz1,lelv)
      real wzuf(lx1,ly1,lz1,lelv)
      real wxbar(lx1,ly1,lz1,lelv)
      real wybar(lx1,ly1,lz1,lelv)
      real wzbar(lx1,ly1,lz1,lelv)

c     coefficients of 3rd order approx. deconvolution
c     for q = 3
      real DC0
      parameter (DC0 = 35.0/16.0)
      real DC1
      parameter (DC1 = -29.0/16.0)
      real DC2
      parameter (DC2 = 3.0/4.0)
      real DC3
      parameter (DC3 = -1.0/8.0)
c     for q = 4
      real DC0_q4
      parameter (DC0_q4 = 315.0/128.0)
      real DC1_q4
      parameter (DC1_q4 = -325.0/128.0)
      real DC2_q4
      parameter (DC2_q4 = 190.0/128.0)
      real DC3_q4
      parameter (DC3_q4 = -60.0/128.0)
      real DC4_q4
      parameter (DC4_q4 = 8.0/128.0)


c     common blocks
      common /filter_real/ tij,tij_old,tij_lap,tij_xdiv,tij_ydiv,
     + tij_zdiv,uxdconv,uydconv,uzdconv,dvdt,dvdt_old,
     + vxtemp,vytemp,vztemp,FTS,qbar,FTS2,freg,u2,u3,u4,u5,
     + wxuf,wyuf,wzuf,wxbar,wybar,wzbar

