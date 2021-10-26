      real vxv, vyv, vzv ! volume average of velocity
      real Af ! actual forcing parameter
      real kin, eps ! turbulent kinetic energy and dissipation
      real velrms ! rms velocity
      real lambda ! transverse Taylor-scale
      real relambda ! Taylor-scale Reynolds number
      real rel ! turbulence Reynolds number
      real teddy ! eddy time turnover time
      real leddy ! lengthscale
      real teta, leta ! Kolmogorov time- and lengthscale

      real tauw ! average wall shear stress
      real utau ! friction velocity
      real Retau ! Reynolds tau
      real ycoords(ly1*lely) ! y-coordinates
      real uprof(ly1*lely) ! velocity profile u(y)

c     time averages
      real uavg(lx1,ly1,lz1,lelv),urms(lx1,ly1,lz1,lelv),
     + vrms(lx1,ly1,lz1,lelv),wrms(lx1,ly1,lz1,lelv),
     + uvms(lx1,ly1,lz1,lelv)
      real uavg_dconv(lx1,ly1,lz1,lelv),urms_dconv(lx1,ly1,lz1,lelv),
     + vrms_dconv(lx1,ly1,lz1,lelv),wrms_dconv(lx1,ly1,lz1,lelv),
     + uvms_dconv(lx1,ly1,lz1,lelv)
      real tijavg_xx(lx1,ly1,lz1,lelv),tijavg_yy(lx1,ly1,lz1,lelv),
     + tijavg_zz(lx1,ly1,lz1,lelv),tijavg_xy(lx1,ly1,lz1,lelv)

c     planar averages of the time averages
      real uavg_y(ly1*lely),urms_y(ly1*lely),
     + vrms_y(ly1*lely),wrms_y(ly1*lely),
     + uvms_y(ly1*lely)
      real uavg_dconv_y(ly1*lely),urms_dconv_y(ly1*lely),
     + vrms_dconv_y(ly1*lely),wrms_dconv_y(ly1*lely),
     + uvms_dconv_y(ly1*lely)
      real tijavg_xx_y(ly1*lely),tijavg_yy_y(ly1*lely),
     + tijavg_zz_y(ly1*lely),tijavg_xy_y(ly1*lely)

c     Reynolds stresses
      real Ruu(ly1*lely),Rvv(ly1*lely),Rww(ly1*lely),Ruv(ly1*lely)
      real Ruu_dconv(ly1*lely),Rvv_dconv(ly1*lely),Rww_dconv(ly1*lely),
     + Ruv_dconv(ly1*lely)

c     common block
      common /slf_real/ vxv, vyv, vzv, Af, kin, eps, velrms,
     + lambda, relambda, rel, teddy, leddy, teta, leta, tauw, utau,
     + Retau, ycoords, uprof, uavg, urms, vrms, wrms, uvms, uavg_dconv,
     + urms_dconv, vrms_dconv, wrms_dconv, uvms_dconv, tijavg_xx,
     + tijavg_yy, tijavg_zz, tijavg_xy, uavg_y, urms_y, vrms_y, wrms_y,
     + uvms_y, uavg_dconv_y, urms_dconv_y, vrms_dconv_y, wrms_dconv_y,
     + uvms_dconv_y, tijavg_xx_y, tijavg_yy_y, tijavg_zz_y, tijavg_xy_y,
     + Ruu, Rvv, Rww, Ruv, Ruu_dconv, Rvv_dconv, Rww_dconv, Ruv_dconv