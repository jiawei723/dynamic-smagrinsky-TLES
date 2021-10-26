c-----------------------------------------------------------------------
c>    set variable properties
      subroutine usr_uservp(ifield,istep,utrans,udiff,temp)
      implicit none
c     interface
      include 'SIZE'
      include 'lib/shared_data.f'
      integer ifield,istep 
      real utrans,udiff,temp 
c     implementation

      udiff  = 0.0
      utrans = 0.0

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
      include 'lib/slf_data.f'  
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

