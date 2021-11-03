c-----------------------------------------------------------------------
c>    set variable properties
      subroutine usr_uservp(ifield,istep,utrans,udiff,temp)
      implicit none
c     interface
      include 'SIZE'
      include 'lib/shared_data.f'
      include 'lib/dysmg_data.f'
      integer ifield,istep,ix,iy,iz,eg
      real utrans,udiff,temp
c     implementation

      udiff  = ediff(ix,iy,iz,eg)
      utrans = 1.0

      return
      end

c-----------------------------------------------------------------------
c>    set acceleration term <BR>
c>    note: this is an acceleration term, not a force!
c>    thus, ffx will subsequently be multiplied by rho(x,t). <BR>
c>    the linear forcing scheme is from
c>    "Linear forcing in numerical simulations of isotropic turbulence:
c>    Physical space implementations and convergence properties",
c>    Rosales et al., 2005
c>    the low wave number forcing is from
c>    "An examination of forcing in direct numerical simulation of
c>    turbulence", Eswaran and Pope,1988
      subroutine usr_userf(istep,ix,iy,iz,e,ux,uy,uz,ffx,ffy,ffz)
      implicit none
c     interface
      include 'SIZE'
      include 'lib/shared_data.f'
      include 'lib/filtering_data.f'
      include 'lib/fluid_data.f'
      include 'lib/forcing_data.f'
      include 'lib/dysmg_data.f'
      integer istep,ix,iy,iz,e
      real ux,uy,uz,ffx,ffy,ffz
c     implementation

c     linear forcing scheme force
      if (iflfs) then
         ffx = Af*(ux - vxv)
         ffy = Af*(uy - vyv)
         ffz = Af*(uz - vzv)
c     low wave number forcing
      elseif (iflwf) then
         ffx = flwx(ix,iy,iz,e)
         ffy = flwy(ix,iy,iz,e)
         ffz = flwz(ix,iy,iz,e)
c     forcing field consisting point sources
      elseif (ifps .and. (istep .lt. ips)) then
         ffx = fpsx(ix,iy,iz,e)*(-1.0)**istep
         ffy = fpsy(ix,iy,iz,e)*(-1.0)**istep
         ffz = fpsz(ix,iy,iz,e)*(-1.0)**istep
      else
         ffx = 0.0
         ffy = 0.0
         ffz = 0.0
      endif

c     contant volume flow in x-direction
      if (ifcvf) then
         ffx = ffx + fvf
      endif

c     TLES: apply the res. stress
      if (iftles) then
         ffx = ffx - tij_xdiv(ix,iy,iz,e)
         ffy = ffy - tij_ydiv(ix,iy,iz,e)
         ffz = ffz - tij_zdiv(ix,iy,iz,e)
      endif

c     TLES: secondary regularization
      if (if2ndreg .and. iftles) then
         ffx = ffx - freg(1,ix,iy,iz,e)
         ffy = ffy - freg(2,ix,iy,iz,e)
         ffz = ffz - freg(3,ix,iy,iz,e)
      endif

      return
      end

c-----------------------------------------------------------------------
c>    set source term
      subroutine usr_userq(ix,iy,iz,e,ifield,qvol)
      implicit none
c     interface
      include 'SIZE'
      include 'lib/shared_data.f'
      integer ix,iy,iz,e,ifield
      real qvol
c     implementation

      qvol = 0.0

      return
      end

c-----------------------------------------------------------------------
c>    set up initial conditions
      subroutine usr_useric(ifield,pi,x,y,z,ux,uy,uz,temp)
      implicit none
c     interface
      include 'SIZE'
      include 'lib/shared_data.f'
      integer ifield
      real pi,x,y,z,ux,uy,uz,temp
c     implementation

      ux = 0.0
      uy = 0.0
      uz = 0.0
      temp = 0.0

      return
      end

c-----------------------------------------------------------------------
c>    initialize and compute the flow
      subroutine usr_chk1(istep,iostep,dt,torder,volvm1,pi,bm1,binvm1,
     + xm1,ym1,zm1,vx,vxlag,vy,vylag,vz,vzlag,pr,rho,mu)
      implicit none
c     interface
      include 'SIZE'
      include 'lib/shared_data.f'
      include 'lib/statistics_data.f'
      integer istep,iostep
      real dt,torder,volvm1,pi,rho,mu
      real bm1(lx1,ly1,lz1,lelt),binvm1(lx1,ly1,lz1,lelt),
     + xm1(lx1,ly1,lz1,lelt),ym1(lx1,ly1,lz1,lelt),zm1(lx1,ly1,lz1,lelt)
      real vx(lx1,ly1,lz1,lelv),vy(lx1,ly1,lz1,lelv),
     + vz(lx1,ly1,lz1,lelv)
      real pr(lx2,ly2,lz2,lelv)
      real vxlag(lx1,ly1,lz1,lelv,int(torder)-1),
     + vylag(lx1,ly1,lz1,lelv,int(torder)-1),
     + vzlag(lx1,ly1,lz1,lelv,int(torder)-1)
c     implementation
      integer atorder
      character*10 fieldName

      atorder=int(torder)-1

c     initialize data
      if (istep .eq. 0) then
         call shared_init(xm1,ym1,zm1)
         if (iflwf) call forcing_init(pi)
         if (ifps) call forcing_ps(xm1,ym1,zm1)
         if (ifksm) then
            call fluid_ksim(pi)
            call fluid_init(xm1,ym1,zm1,vx,vy,vz,pi,pr)
         else
            call slf_init()
         endif

c     Compute eddy viscosity using dynamic smagorinsky model
      if(ifuservp) then
        if(nid.eq.0) write(6,*) 'Calculating eddy visosity'
        do e=1,nelv
           call eddy_visc(ediff,e)
        enddo
        call copy(t,ediff,n)
      endif

c        TLES
         if (iftles) then
            if (iffullrestart .and. ifreTLES) then
               call filter_restart(vx,vy,vz)
            else
               call filter_init(vx,vy,vz,dt)
            endif
         endif

c        initialize the statistics
         if (ifstat) then
            call statistics_init()
         endif
      endif

c     low wave number forcing
      if (iflwf) then
         call forcing_setzero()
         if (nid .eq. 0) call forcing_computeterm(dt)
         call forcing_interpolate(xm1,ym1,zm1,dt)
      endif

c     constant volume flow
      if (ifcvf) then
         call set_forcing(vx,dt,istep,bm1,volvm1)
      endif

c     TLES
      if (iftles) then
         if (ifvalid) then
            call filter_valid(xm1,ym1,zm1,vx,vy,vz,istep,dt,pi)
         endif
         call filter_diffterm(binvm1)
         call filter_stress(istep,vx,vxlag,vy,vylag,vz,vzlag,dt,atorder)
         call filter_div(binvm1)
         call filter_deconv(vx,vy,vz)
         if (if2ndreg) then ! compute freg - used in usr_userf
            if (ifadreg) then ! approx. deconv. regularization (Pruett)
               call filter_ADregularization(istep,vx,vy,vz,dt)
            else ! own regularization
               call filter_regularization(istep,vx,vy,vz,dt)
            endif
         endif
         if (ifdiag .and. nid .eq. 0) then
            call filter_diag(xm1,ym1,zm1,vx,vy,vz,dt)
         endif
      endif

c     compute statistics
      call slf_vk_compute(istep,iostep,bm1,vx,vy,vz,volvm1)
      if (ifchannel) then
         call slf_channel_compute(istep,bm1,vx,vy,vz,volvm1,rho,mu/rho)
         if (istep*dt .ge. 10.0) then ! start after 10 seconds
            call slf_channel_compute_stresses(istep*dt,vx,vy,vz)
         endif
      elseif (istep .gt. 0) then
         call slf_trb_compute(bm1,vx,vy,vz,volvm1,mu/rho)
      endif

c     periodic hill statistics
      if (ifstat .and. (mod(istep,profilestep) .eq. 0)) then
c        compute the average slice (contraction in z)
         call statistics_computeaverage_z()

c        interpolates the 'u_avgz' field to probe points
         call statistics_interpolate_samples(nprobes)
      endif

      return
      end

c-----------------------------------------------------------------------
c>    write out data
      subroutine usr_chk2(istep,iostep,nsteps,ifreguo,ifxyo,time,pr,t,
     + xm1,ym1,zm1,vx,vy,vz)
      implicit none
c     interface
      include 'SIZE'
      include 'lib/shared_data.f'
      include 'lib/filtering_data.f'
      include 'lib/forcing_data.f'
      include 'lib/statistics_data.f'
      integer istep,iostep,nsteps
      logical ifreguo,ifxyo
      real time,pr(lx2,ly2,lz2,lelv),t(lx1,ly1,lz1,lelt,ldimt)
      real xm1(lx1,ly1,lz1,lelt),ym1(lx1,ly1,lz1,lelt),
     + zm1(lx1,ly1,lz1,lelt),vx(lx1,ly1,lz1,lelt),vy(lx1,ly1,lz1,lelt),
     + vz(lx1,ly1,lz1,lelt)
c     implementation

c     periodic hill stats
      if (ifstat .and. (mod(istep,profilestep) .eq. 0)) then

c        write out the average slice u(x,y)
         call statistics_write_slice(xm1,ym1,zm1,time)

c        write out the average vel. profiles
         if (nid .eq. 0) then
            call statistics_write_profiles(time)
         endif
      endif

c     write out data to regular and irregular grid
      if (mod(istep,iostep) .eq. 0) then
c        coordinates only in first file
         ifxyo = .true.
         if (istep .gt. 0) ifxyo = .false.
         ifreguo = .false.
c        call prepost(.true.,'his') ! b3d_reg

         if (iftles) then
c            call outpost(uxdconv,uydconv,uzdconv,pr,t,'udc')
         endif
         if (iftles) then
c            call outpost(tij(:,:,:,:,1,1),tij_xdiv,
c     + dvdt(:,:,:,:,1),pr,t,'txx')
c         call outpost(tij(:,:,:,:,1,2),tij_xdiv,tij_lap(:,:,:,:,1,1),
c     + pr,t,'txx')
         !call outpost(flwx,flwy,flwz,pr,t,'flw')
         endif
         ifreguo = .false.
      endif

      if (mod(istep,statisticstep) .eq. 0) then
         if (nid .eq. 0) then
            if (ifchannel) then
               call slf_channel_write(time)
            else
               call slf_write(time)
            endif
         endif
      endif

      if ((istep .eq. nsteps) .and. (iftles)) then
         call filter_write_restart()
      endif
c     write out u(t) each time step after istep >= 0.5*nsteps
      if ((nid .eq. 0) .and. (iftimeseries) .and.
     + (istep .le. 400)) then
         if ((ifadreg) .and. (istep .le. 400)) then
            call slf_ADRvalidation(time,xm1,ym1,zm1,vx,vy,vz)
         elseif (.not. ifadreg) then
            call slf_timeseries(time,xm1,ym1,zm1,vx,vy,vz)
         endif
      endif

      return
      end

c-----------------------------------------------------------------------
c>     Compute Lij for dynamic Smagorinsky model:
c                    _   _      _______
c          L_ij  :=  u_i u_j  - u_i u_j
c
      subroutine comp_lij(lij,u,v,w,fu,fv,fw,fh,fht,e)

      include 'SIZE'
c
      integer e
c
      real lij(lx1*ly1*lz1,3*ldim-3)
      real u  (lx1*ly1*lz1,lelv)
      real v  (lx1*ly1*lz1,lelv)
      real w  (lx1*ly1*lz1,lelv)
      real fu (1) , fv (1) , fw (1)
     $   , fh (1) , fht(1)

      call tens3d1(fu,u(1,e),fh,fht,nx1,nx1)  ! fh x fh x fh x u
      call tens3d1(fv,v(1,e),fh,fht,nx1,nx1)
      call tens3d1(fw,w(1,e),fh,fht,nx1,nx1)

      n = nx1*ny1*nz1
      do i=1,n
         lij(i,1) = fu(i)*fu(i)
         lij(i,2) = fv(i)*fv(i)
         lij(i,3) = fw(i)*fw(i)
         lij(i,4) = fu(i)*fv(i)
         lij(i,5) = fv(i)*fw(i)
         lij(i,6) = fw(i)*fu(i)
      enddo

      call col3   (fu,u(1,e),u(1,e),n)    !  _______
      call tens3d1(fv,fu,fh,fht,nx1,nx1)  !  u_1 u_1
      call sub2   (lij(1,1),fv,n)

      call col3   (fu,v(1,e),v(1,e),n)    !  _______
      call tens3d1(fv,fu,fh,fht,nx1,nx1)  !  u_2 u_2
      call sub2   (lij(1,2),fv,n)

      call col3   (fu,w(1,e),w(1,e),n)    !  _______
      call tens3d1(fv,fu,fh,fht,nx1,nx1)  !  u_3 u_3
      call sub2   (lij(1,3),fv,n)

      call col3   (fu,u(1,e),v(1,e),n)    !  _______
      call tens3d1(fv,fu,fh,fht,nx1,nx1)  !  u_1 u_2
      call sub2   (lij(1,4),fv,n)

      call col3   (fu,v(1,e),w(1,e),n)    !  _______
      call tens3d1(fv,fu,fh,fht,nx1,nx1)   !  u_2 u_3
      call sub2   (lij(1,5),fv,n)

      call col3   (fu,w(1,e),u(1,e),n)    !  _______
      call tens3d1(fv,fu,fh,fht,nx1,nx1)  !  u_3 u_1
      call sub2   (lij(1,6),fv,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine comp_mij(mij,sij,dg2,fs,fi,fh,fht,nt,e)
c
c     Compute Mij for dynamic Smagorinsky model:
c
c                     2 _  ____     _______
c          M_ij  :=  a  S  S_ij  -  S  S_ij
c
      include 'SIZE'
c
      integer e
c
      real mij(lx1*ly1*lz1,3*ldim-3)
      real dg2(lx1*ly1*lz1,lelv)
      real fs (1) , fi (1) , fh (1) , fht(1)

      real magS(lx1*ly1*lz1)
      real sij (lx1*ly1*lz1*ldim*ldim)

      integer imap(6)
      data imap / 0,4,8,1,5,2 /

      n = nx1*ny1*nz1

      call mag_tensor_e(magS,sij)
      call cmult(magS,2.0,n)

c     Filter S
      call tens3d1(fs,magS,fh,fht,nx1,nx1)  ! fh x fh x fh x |S|

c     a2 is the test- to grid-filter ratio, squared

      a2 = nx1-1       ! nx1-1 is number of spaces in grid
      a2 = a2 /(nt-1)  ! nt-1 is number of spaces in filtered grid

      do k=1,6
         jj = n*imap(k) + 1
         call col3   (fi,magS,sij(jj),n)
         call tens3d1(mij(1,k),fi,fh,fht,nx1,nx1)  ! fh x fh x fh x (|S| S_ij)
         call tens3d1(fi,sij(jj),fh,fht,nx1,nx1)  ! fh x fh x fh x S_ij
         do i=1,n
            mij(i,k) = (a2**2 * fs(i)*fi(i) - mij(i,k))*dg2(i,e)
         enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine eddy_visc(ediff,e)
c
c     Compute eddy viscosity using dynamic smagorinsky model
c
      include 'SIZE'
      include 'TOTAL'
      include 'ZPER'

      real ediff(nx1*ny1*nz1,nelv)
      integer e

      common /dynsmg/ sij (lx1*ly1*lz1,ldim,ldim)
     $              , mij (lx1*ly1*lz1,3*ldim-3)
     $              , lij (lx1*ly1*lz1,3*ldim-3)
     $              , dg2 (lx1*ly1*lz1,lelv)
     $              , num (lx1*ly1*lz1,lelv)
     $              , den (lx1*ly1*lz1,lelv)
     $              , snrm(lx1*ly1*lz1,lelv)
     $              , numy(ly1*lely),deny(ly1*lely),yy(ly1*lely)
      real sij,mij,lij,dg2,num,den,snrm,numy,deny,yy

      parameter(lxyz=lx1*ly1*lz1)
      common /xzmp0/ ur (lxyz) , us (lxyz) , ut (lxyz)
      real           vr (lxyz) , vs (lxyz) , vt (lxyz)
     $     ,         wr (lxyz) , ws (lxyz) , wt (lxyz)
      common /xzmp1/ w1(lx1*lelv),w2(lx1*lelv)

      !! NOTE CAREFUL USE OF EQUIVALENCE HERE !!
      equivalence (vr,lij(1,1)),(vs,lij(1,2)),(vt,lij(1,3))
     $          , (wr,lij(1,4)),(ws,lij(1,5)),(wt,lij(1,6))

      common /sgsflt/ fh(lx1*lx1),fht(lx1*lx1),diag(lx1)

      integer nt
      save    nt
      data    nt / -9 /

      ntot = nx1*ny1*nz1

      if (nt.lt.0) call
     $   set_ds_filt(fh,fht,nt,diag,nx1)! dyn. Smagorinsky filter

      call comp_gije(sij,vx(1,1,1,e),vy(1,1,1,e),vz(1,1,1,e),e)
      call comp_sije(sij)

      call mag_tensor_e(snrm(1,e),sij)
      call cmult(snrm(1,e),2.0,ntot)

      call set_grid_spacing(dg2)
      call comp_mij   (mij,sij,dg2,ur,us,fh,fht,nt,e)

      call comp_lij   (lij,vx,vy,vz,ur,us,ut,fh,fht,e)

c     Compute numerator (ur) & denominator (us) for Lilly contraction

      n = nx1*ny1*nz1
      do i=1,n
         ur(i) = mij(i,1)*lij(i,1)+mij(i,2)*lij(i,2)+mij(i,3)*lij(i,3)
     $      + 2*(mij(i,4)*lij(i,4)+mij(i,5)*lij(i,5)+mij(i,6)*lij(i,6))
         us(i) = mij(i,1)*mij(i,1)+mij(i,2)*mij(i,2)+mij(i,3)*mij(i,3)
     $      + 2*(mij(i,4)*mij(i,4)+mij(i,5)*mij(i,5)+mij(i,6)*mij(i,6))
      enddo

c     smoothing numerator and denominator in time
      call copy (vr,ur,nx1*nx1*nx1)
      call copy (vs,us,nx1*nx1*nx1)

      beta1 = 0.0                   ! Temporal averaging coefficients
      if (istep.gt.1) beta1 = 0.9   ! Retain 90 percent of past
      beta2 = 1. - beta1

      do i=1,n
         num (i,e) = beta1*num(i,e) + beta2*vr(i)
         den (i,e) = beta1*den(i,e) + beta2*vs(i)
      enddo


      if (e.eq.nelv) then  ! planar avg and define nu_tau

         call dsavg(num)   ! average across element boundaries
         call dsavg(den)

         call planar_average_s      (numy,num,w1,w2)
c        call wall_normal_average_s (numy,ny1,nely,w1,w2)
         call planar_fill_s         (num,numy)

         call planar_average_s      (deny,den,w1,w2)
c        call wall_normal_average_s (deny,ny1,nely,w1,w2)
         call planar_fill_s         (den,deny)

         call planar_average_s(yy,ym1,w1,w2)
      endif

      return
      end