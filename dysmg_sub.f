c-----------------------------------------------------------------------
c>    initialize the turbulent statistics for the case that the initial
c>    velocity field is zero
      subroutine slf_init()
      implicit none
c     interface
      include 'SIZE'
      include 'lib/dysmg_data.f'
c     implementation
      integer n

      vxv = 0.0
      vyv = 0.0
      vzv = 0.0
      kin = 0.0
      eps = 0.0
      velrms = 0.0
      lambda = 0.0
      relambda = 0.0
      teddy = 0.0
      leddy = 0.0
      rel = 0.0
      teta = 0.0
      leta = 0.0

      n = ly1*lely
      call rzero(Ruu,n)
      call rzero(Rvv,n)
      call rzero(Rww,n)
      call rzero(Ruv,n)

      return
      end

c-----------------------------------------------------------------------
c>    compute volume average of velocity and turbulent kinetic energy
c>    "Turbulent Flows", Stephen B. Pope, 2011
      subroutine slf_vk_compute(istep,iostep,bm1,vx,vy,vz,volvm1)
      implicit none
c     interface
      include 'SIZE'
      include 'lib/shared_data.f'
      include 'lib/filtering_data.f'
      include 'lib/dysmg_data.f'
      integer istep,iostep
      real bm1(lx1,ly1,lz1,lelt),vx(lx1,ly1,lz1,lelv),
     + vy(lx1,ly1,lz1,lelv),vz(lx1,ly1,lz1,lelv),volvm1
      real glsc2, glsc3
c     implementation
      integer n ! local number of gridpoints
      real vxvx, vyvy, vzvz ! scalar products

c     compute volume average of velocity
      n = nx1*ny1*nz1*nelv
      if (iftles) then
         vxv = glsc2(bm1,uxdconv,n)/volvm1
         vyv = glsc2(bm1,uydconv,n)/volvm1
         vzv = glsc2(bm1,uzdconv,n)/volvm1
      else
         vxv = glsc2(bm1,vx,n)/volvm1
         vyv = glsc2(bm1,vy,n)/volvm1
         vzv = glsc2(bm1,vz,n)/volvm1
      endif

c     compute urms square and the kinetic energy
      if (mod(istep,statisticstep) .eq. 0 .or. (ifuf .and. iflfs)) then
         if (iftles) then
            vxvx = glsc3(uxdconv,uxdconv,bm1,n)/volvm1
            vyvy = glsc3(uydconv,uydconv,bm1,n)/volvm1
            vzvz = glsc3(uzdconv,uzdconv,bm1,n)/volvm1
         else
            vxvx = glsc3(vx,vx,bm1,n)/volvm1
            vyvy = glsc3(vy,vy,bm1,n)/volvm1
            vzvz = glsc3(vz,vz,bm1,n)/volvm1
         endif
         kin = 0.5*((vxvx-vxv**2)+(vyvy-vyv**2)+(vzvz-vzv**2))
      endif

c     update (scale) the forcing parameter Af
      if (ifuf .and. iflfs) then
         Af = epsr/(2.0*kin)
      endif

      return
      end

c-----------------------------------------------------------------------
c>    compute turbulent statistics
c>    "Turbulent Flows", Stephen B. Pope, 2011
      subroutine slf_trb_compute(bm1,vx,vy,vz,volvm1,nu)
      implicit none
c     interface
      include 'SIZE'
      include 'lib/shared_data.f'
      include 'lib/filtering_data.f'
      include 'lib/dysmg_data.f'
      real bm1(lx1,ly1,lz1,lelt),vx(lx1,ly1,lz1,lelv),
     + vy(lx1,ly1,lz1,lelv),vz(lx1,ly1,lz1,lelv),volvm1
      real glsc2
      real nu ! viscosity
c     implementation
      integer ntot ! local number of gridpoints
      integer i,j,k,l
      real vxx(lx1,ly1,lz1,lelv),vxy(lx1,ly1,lz1,lelv),
     + vxz(lx1,ly1,lz1,lelv), ! d_dx(vx), d_dy(vx), d_dz(vx)
     + vyx(lx1,ly1,lz1,lelv),vyy(lx1,ly1,lz1,lelv),
     + vyz(lx1,ly1,lz1,lelv), ! d_dx(vy), d_dy(vy), d_dz(vy)
     + vzx(lx1,ly1,lz1,lelv),vzy(lx1,ly1,lz1,lelv),
     + vzz(lx1,ly1,lz1,lelv), ! d_dx(vz), d_dy(vz), d_dz(vz)
     + sijsij(lx1,ly1,lz1,lelv)

c     initialize computation of sij
      ntot = nx1*ny1*nz1*nelv
      call rzero(vxx,ntot)
      call rzero(vxy,ntot)
      call rzero(vxz,ntot)
      call rzero(vyx,ntot)
      call rzero(vyy,ntot)
      call rzero(vyz,ntot)
      call rzero(vzx,ntot)
      call rzero(vzy,ntot)
      call rzero(vzz,ntot)
      call rzero(sijsij,ntot)

c     compute the velocity gradients
      if (iftles) then
         call gradm1(vxx,vxy,vxz,uxdconv)
         call gradm1(vyx,vyy,vyz,uydconv)
         call gradm1(vzx,vzy,vzz,uzdconv)
      else
         call gradm1(vxx,vxy,vxz,vx)
         call gradm1(vyx,vyy,vyz,vy)
         call gradm1(vzx,vzy,vzz,vz)
      endif

c     compute sijsij
      do l = 1,nelv ! spectral elements on this processor
      do k = 1,nz1
      do j = 1,ny1
      do i = 1,nx1
         sijsij(i,j,k,l) =
     + vxx(i,j,k,l)**2.0 + vyy(i,j,k,l)**2.0 + vzz(i,j,k,l)**2.0 +
     + 1.0/2.0*(vxy(i,j,k,l) + vyx(i,j,k,l))**2.0 +
     + 1.0/2.0*(vxz(i,j,k,l) + vzx(i,j,k,l))**2.0 +
     + 1.0/2.0*(vyz(i,j,k,l) + vzy(i,j,k,l))**2.0
      enddo
      enddo
      enddo
      enddo

      eps = 2.0*nu*glsc2(bm1,sijsij,ntot)/volvm1
      velrms = sqrt(2.0/3.0*kin)
      lambda = sqrt(15.0*nu*velrms**2.0/eps)
      relambda = velrms*lambda/nu
      teddy = (2.0/3.0*kin)/eps
      leddy = kin**(3.0/2.0)/eps
      teta = (nu/eps)**(1.0/2.0)
      leta = (nu**3.0/eps)**(1.0/4.0)

      return
      end

c-----------------------------------------------------------------------
c>    write out volume average of velocity and turbulent statistics
      subroutine slf_write(time)
      implicit none
c     interface
      include 'SIZE'
      include 'lib/shared_data.f'
      include 'lib/dysmg_data.f'
c     implementation
      integer ios
      real time
      character*200 fname

      integer iclld ! timestep
      save    iclld
      data    iclld  /0/

      write(fname,71)
      if (iclld .eq. 0) then
         open(unit=1,file=fname,status='unknown',iostat=ios)
         if (ios .eq. 0) close(1,status='delete')
         open(unit=1,file=fname,status="new",action="write")
      else
         open(unit=1,file=fname,status="old",position="append",
     + action="write")
      endif

      write(1,72) iclld, time
      write(1,73) vxv, vyv, vzv, Af, kin, eps, velrms,
     + lambda, relambda, teddy, leddy, rel, teta, leta
      close(1)

      iclld = iclld + 1

   71 format('out/fluid/slf_m.txt')
   72 format(I10,E17.9)
   73 format(14E17.9)

      return
      end

c-----------------------------------------------------------------------
c>    compute channel flow statistics
      subroutine slf_channel_compute(istep,bm1,vx,vy,vz,volvm1,rho,nu)
      implicit none
c     interface
      include 'SIZE'
      include 'lib/shared_data.f'
      include 'lib/dysmg_data.f'
      integer istep
      real bm1(lx1,ly1,lz1,lelt),vx(lx1,ly1,lz1,lelv),
     + vy(lx1,ly1,lz1,lelv),vz(lx1,ly1,lz1,lelv),volvm1
      real glsc2
      real rho,nu ! density,kinematic viscosity
c     implementation
      integer n ! local number of gridpoints
      integer i,j,k,l
      real vxx(lx1,ly1,lz1,lelv),vxy(lx1,ly1,lz1,lelv),
     + vxz(lx1,ly1,lz1,lelv), ! d_dx(vx), d_dy(vx), d_dz(vx)
     + vyx(lx1,ly1,lz1,lelv),vyy(lx1,ly1,lz1,lelv),
     + vyz(lx1,ly1,lz1,lelv), ! d_dx(vy), d_dy(vy), d_dz(vy)
     + vzx(lx1,ly1,lz1,lelv),vzy(lx1,ly1,lz1,lelv),
     + vzz(lx1,ly1,lz1,lelv), ! d_dx(vz), d_dy(vz), d_dz(vz)
     + sijsij(lx1,ly1,lz1,lelv)

c     compute the wall shear stress tauw
      call slf_channel_wss(tauw,iftles)
      tauw=tauw/(ldmx*ldmz) ! divide by area

c     average u(y)
      if (istep .eq. 0) then
         call slf_channel_profile(ycoords,iftles)
      else
         call slf_channel_profile(uprof,iftles)
      endif

c     initialize computation of sij
      call rzero(vxx,lx1*ly1*lz1*lelv)
      call rzero(vxy,lx1,ly1,lz1,lelv)
      call rzero(vxz,lx1,ly1,lz1,lelv)
      call rzero(vyx,lx1,ly1,lz1,lelv)
      call rzero(vyy,lx1,ly1,lz1,lelv)
      call rzero(vyz,lx1,ly1,lz1,lelv)
      call rzero(vzx,lx1,ly1,lz1,lelv)
      call rzero(vzy,lx1,ly1,lz1,lelv)
      call rzero(vzz,lx1,ly1,lz1,lelv)
      call rzero(sijsij,lx1,ly1,lz1,lelv)

c     compute the velocity gradients
      call gradm1(vxx,vxy,vxz,vx)
      call gradm1(vyx,vyy,vyz,vy)
      call gradm1(vzx,vzy,vzz,vz)

c     compute sijsij
      do l = 1,nelv ! spectral elements on this processor
      do k = 1,lz1
      do j = 1,ly1
      do i = 1,lx1
         sijsij(i,j,k,l) =
     + vxx(i,j,k,l)**2.0 + vyy(i,j,k,l)**2.0 + vzz(i,j,k,l)**2.0 +
     + 1.0/2.0*(vxy(i,j,k,l) + vyx(i,j,k,l))**2.0 +
     + 1.0/2.0*(vxz(i,j,k,l) + vzx(i,j,k,l))**2.0 +
     + 1.0/2.0*(vyz(i,j,k,l) + vzy(i,j,k,l))**2.0
      enddo
      enddo
      enddo
      enddo

      n = nx1*ny1*nz1*nelv
      eps = 2.0*nu*glsc2(bm1,sijsij,n)/volvm1
      velrms = sqrt(2.0/3.0*kin)
      lambda = sqrt(15.0*nu*velrms**2.0/eps)
      relambda = velrms*lambda/nu
      teddy = (2.0/3.0*kin)/eps
      leddy = kin**(3.0/2.0)/eps
      rel = sqrt(kin)*leddy/nu
      teta = (nu/eps)**(1.0/2.0)
      leta = (nu**3.0/eps)**(1.0/4.0)
      utau = sqrt(abs(tauw)/rho)
      Retau = utau*0.5*ldmy/nu

      return
      end

c-----------------------------------------------------------------------
c>    compute and the average wall shear stress
      subroutine slf_channel_wss(tauw,iftles)
c     interface
      include 'SIZE'
      include 'TOTAL'
      include 'ZPER'
      real tauw ! wall shear stress
      logical iftles
c     implementation
      real x0(3)
      save x0

      common /ctorqn/ dragxn(0:maxobj),dragpxn(0:maxobj),
     + dragvxn(0:maxobj), dragyn(0:maxobj),dragpyn(0:maxobj),
     + dragvyn(0:maxobj), dragzn(0:maxobj),dragpzn(0:maxobj),
     + dragvzn(0:maxobj), torqxn(0:maxobj),torqpxn(0:maxobj),
     + torqvxn(0:maxobj), torqyn(0:maxobj),torqpyn(0:maxobj),
     + torqvyn(0:maxobj), torqzn(0:maxobj),torqpzn(0:maxobj),
     + torqvzn(0:maxobj), dpdxn_mean,dpdyn_mean,dpdzn_mean, dgtqn(3,4)

c     initialization
      if (istep.eq.0) then
         call slf_channel_obj()  ! objects for surface integrals
         call rzero(x0,3)        ! torque w.r.t. x0
      endif

c     integrate the wall shear stress
      call slf_channel_torque(1.0,x0,iftles)

c     take the average of the upper and lower wall
      tauw=0.5*(dragxn(1)+dragxn(2))

      return
      end

c-----------------------------------------------------------------------
c>    define objects for surface integrals
      subroutine slf_channel_obj()
c     interface
      include 'SIZE'
      include 'TOTAL'
c     implementation
      integer e,f

c     define new objects
      nobj = 2			! for Periodic
      iobj = 0
      do ii=nhis+1,nhis+nobj
         iobj = iobj+1
         hcode(10,ii) = 'I'
         hcode( 1,ii) = 'F'
         hcode( 2,ii) = 'F'
         hcode( 3,ii) = 'F'
         lochis(1,ii) = iobj
      enddo
      nhis = nhis + nobj

      if (maxobj.lt.nobj) write(6,*) 'increase maxobj in SIZEu. rm *.o'
      if (maxobj.lt.nobj) call exitt

      nxyz = nx1*ny1*nz1
      do e=1,nelv
      do f=1,2*ndim
         if (cbc(f,e,1).eq.'W  ') then
            iobj = 0
            if (f.eq.1) iobj=1  ! lower wall
            if (f.eq.3) iobj=2  ! upper wall
            if (iobj.gt.0) then
               nmember(iobj) = nmember(iobj) + 1
               mem = nmember(iobj)
               ieg = lglel(e)
               object(iobj,mem,1) = ieg
               object(iobj,mem,2) = f
            endif
         endif
      enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
c>    compute torque about point x0
c>    scale is a user-supplied multiplier so that results may be
c>    scaled to any convenient non-dimensionalization.
      subroutine slf_channel_torque(scale,x0,iftles)
c     interface
      include 'SIZE'
      include 'TOTAL'
      include 'lib/filtering_data.f'
      logical iftles
c     implementation
      real x0(3),w1(0:maxobj)

      common /cvflow_r/ flow_rate,base_flow,domain_length,xsec
     $                , scale_vf(3)

      common /scrns/         sij (lx1*ly1*lz1*6*lelv)
      common /scrcg/         pm1 (lx1,ly1,lz1,lelv)
      common /scrsf/         xm0(lx1,ly1,lz1,lelt)
     $,                      ym0(lx1,ly1,lz1,lelt)
     $,                      zm0(lx1,ly1,lz1,lelt)

      parameter (lr=lx1*ly1*lz1)
      common /scruz/         ur(lr),us(lr),ut(lr)
     $                     , vr(lr),vs(lr),vt(lr)
     $                     , wr(lr),ws(lr),wt(lr)

      common /ctorqn/ dragxn(0:maxobj),dragpxn(0:maxobj),
     + dragvxn(0:maxobj), dragyn(0:maxobj),dragpyn(0:maxobj),
     + dragvyn(0:maxobj), dragzn(0:maxobj),dragpzn(0:maxobj),
     + dragvzn(0:maxobj), torqxn(0:maxobj),torqpxn(0:maxobj),
     + torqvxn(0:maxobj), torqyn(0:maxobj),torqpyn(0:maxobj),
     + torqvyn(0:maxobj), torqzn(0:maxobj),torqpzn(0:maxobj),
     + torqvzn(0:maxobj), dpdxn_mean,dpdyn_mean,dpdzn_mean, dgtqn(3,4)

      n = nx1*ny1*nz1*nelv

      call mappr(pm1,pr,xm0,ym0) ! map pressure onto Mesh 1

c    add mean_pressure_gradient.X to p:
      if (param(55).ne.0) then
         dpdxn_mean = -scale_vf(1)
         dpdyn_mean = -scale_vf(2)
         dpdzn_mean = -scale_vf(3)
      endif

      call add2s2(pm1,xm1,dpdxn_mean,n)  ! doesn't work if object is cut by
      call add2s2(pm1,ym1,dpdyn_mean,n)  ! periodic boundary. In this case,
      call add2s2(pm1,zm1,dpdzn_mean,n)  ! set ._mean=0 and compensate in

c    compute sij
      nij = 3
      if (if3d.or.ifaxis) nij=6
      if (iftles) then
         call comp_sij(sij,nij,uxdconv,uydconv,uzdconv,
     + ur,us,ut,vr,vs,vt,wr,ws,wt)
      else
         call comp_sij(sij,nij,vx,vy,vz,ur,us,ut,vr,vs,vt,wr,ws,wt)
      endif

c     fill up viscous array w/ default
      if (istep.lt.1) call cfill(vdiff,param(2),n)
      call cadd2(xm0,xm1,-x0(1),n)
      call cadd2(ym0,ym1,-x0(2),n)
      call cadd2(zm0,zm1,-x0(3),n)
      x1min=glmin(xm0(1,1,1,1),n)
      x2min=glmin(ym0(1,1,1,1),n)
      x3min=glmin(zm0(1,1,1,1),n)
      x1max=glmax(xm0(1,1,1,1),n)
      x2max=glmax(ym0(1,1,1,1),n)
      x3max=glmax(zm0(1,1,1,1),n)

      do i=0,maxobj
         dragpxn(i) = 0
         dragvxn(i) = 0
         dragxn (i) = 0
         dragpyn(i) = 0
         dragvyn(i) = 0
         dragyn (i) = 0
         dragpzn(i) = 0
         dragvzn(i) = 0
         dragzn (i) = 0
         torqpxn(i) = 0
         torqvxn(i) = 0
         torqxn (i) = 0
         torqpyn(i) = 0
         torqvyn(i) = 0
         torqyn (i) = 0
         torqpzn(i) = 0
         torqvzn(i) = 0
         torqzn (i) = 0
      enddo

      nobj = 0
      do ii=1,nhis
        if (hcode(10,ii).EQ.'I') then
          iobj   = lochis(1,ii)
          memtot = nmember(iobj)
          nobj   = max(iobj,nobj)
          if (hcode(1,ii).ne.' ' .or. hcode(2,ii).ne.' ' .or.
     $      hcode(3,ii).ne.' ' ) then
            ifield = 1
c           compute drag for this object
            do mem=1,memtot
               ieg   = object(iobj,mem,1)
               ifc   = object(iobj,mem,2)
               if (gllnid(ieg).eq.nid) then ! this processor has a contribution
                  ie = gllel(ieg)
                  call drgtrq(dgtq,xm0,ym0,zm0,sij,pm1,vdiff,ifc,ie)
                  call cmult(dgtq,scale,12)
                  dragpxn(iobj) = dragpxn(iobj) + dgtqn(1,1) ! pressure
                  dragpyn(iobj) = dragpyn(iobj) + dgtqn(2,1)
                  dragpzn(iobj) = dragpzn(iobj) + dgtqn(3,1)
                  dragvxn(iobj) = dragvxn(iobj) + dgtqn(1,2) ! viscous
                  dragvyn(iobj) = dragvyn(iobj) + dgtqn(2,2)
                  dragvzn(iobj) = dragvzn(iobj) + dgtqn(3,2)
                  torqpxn(iobj) = torqpxn(iobj) + dgtqn(1,3) ! pressure
                  torqpyn(iobj) = torqpyn(iobj) + dgtqn(2,3)
                  torqpzn(iobj) = torqpzn(iobj) + dgtqn(3,3)
                  torqvxn(iobj) = torqvxn(iobj) + dgtqn(1,4) ! viscous
                  torqvyn(iobj) = torqvyn(iobj) + dgtqn(2,4)
                  torqvzn(iobj) = torqvzn(iobj) + dgtqn(3,4)
               endif
            enddo
          endif
        endif
      enddo

c     sum contributions from all processors
      call gop(dragpxn,w1,'+  ',maxobj+1)
      call gop(dragpyn,w1,'+  ',maxobj+1)
      call gop(dragpzn,w1,'+  ',maxobj+1)
      call gop(dragvxn,w1,'+  ',maxobj+1)
      call gop(dragvyn,w1,'+  ',maxobj+1)
      call gop(dragvzn,w1,'+  ',maxobj+1)
      call gop(torqpxn,w1,'+  ',maxobj+1)
      call gop(torqpyn,w1,'+  ',maxobj+1)
      call gop(torqpzn,w1,'+  ',maxobj+1)
      call gop(torqvxn,w1,'+  ',maxobj+1)
      call gop(torqvyn,w1,'+  ',maxobj+1)
      call gop(torqvzn,w1,'+  ',maxobj+1)
      nobj = iglmax(nobj,1)
      do i=1,nobj
         dragxn(i) = dragpxn(i) + dragvxn(i)
         dragyn(i) = dragpyn(i) + dragvyn(i)
         dragzn(i) = dragpzn(i) + dragvzn(i)
         torqxn(i) = torqpxn(i) + torqvxn(i)
         torqyn(i) = torqpyn(i) + torqvyn(i)
         torqzn(i) = torqpzn(i) + torqvzn(i)
         dragpxn(0) = dragpxn (0) + dragpxn (i)
         dragvxn(0) = dragvxn (0) + dragvxn (i)
         dragxn (0) = dragxn  (0) + dragxn  (i)
         dragpyn(0) = dragpyn (0) + dragpyn (i)
         dragvyn(0) = dragvyn (0) + dragvyn (i)
         dragyn (0) = dragyn  (0) + dragyn  (i)
         dragpzn(0) = dragpzn (0) + dragpzn (i)
         dragvzn(0) = dragvzn (0) + dragvzn (i)
         dragzn (0) = dragzn  (0) + dragzn  (i)
         torqpxn(0) = torqpxn (0) + torqpxn (i)
         torqvxn(0) = torqvxn (0) + torqvxn (i)
         torqxn (0) = torqxn  (0) + torqxn  (i)
         torqpyn(0) = torqpyn (0) + torqpyn (i)
         torqvyn(0) = torqvyn (0) + torqvyn (i)
         torqyn (0) = torqyn  (0) + torqyn  (i)
         torqpzn(0) = torqpzn (0) + torqpzn (i)
         torqvzn(0) = torqvzn (0) + torqvzn (i)
         torqzn (0) = torqzn  (0) + torqzn  (i)
      enddo

      i0 = 0
      if (nobj.le.1) i0 = 1  ! one output for single-object case

      return
      end

c-----------------------------------------------------------------------
c>    compute the velocity profile u(y)
      subroutine slf_channel_profile(profile_dummy,iftles)
c     interface
      include 'SIZE'
      include 'TOTAL'
      include 'ZPER'  ! for nelx,nely,nelz
      include 'lib/filtering_data.f'
      real profile_dummy(ly1*lely)
      logical iftles
c     implementation

      common /plane/  uavg_pl(ly1*lely)
     $             ,  urms_pl(ly1*lely)
     $             ,  vrms_pl(ly1*lely)
     $             ,  wrms_pl(ly1*lely)
     $             ,  uvms_pl(ly1*lely)
     $             ,  yy(ly1*lely)
     $             ,  w1(ly1*lely),w2(ly1*lely)
     $             ,  ffx_avg, dragx_avg

      nelx  = NUMBER_ELEMENTS_X
      nely  = NUMBER_ELEMENTS_Y
      nelz  = NUMBER_ELEMENTS_Z

c     get the wall normal coordinates
      if (istep.eq.0) then
         call planar_average_s(profile_dummy,ym1,w1,w2)
      else
         if (iftles) then
            call planar_average_s(profile_dummy,uxdconv,w1,w2)
         else
            call planar_average_s(profile_dummy,vx,w1,w2)
         endif
      endif

      return
      end

c-----------------------------------------------------------------------
c>    compute r-t planar average of quantity u()
      subroutine planar_average_s(ua,u,w1,w2)
c     interface
      include 'SIZE'
      include 'GEOM'
      include 'PARALLEL'
      include 'WZ'
      include 'ZPER'
      real ua(ny1,nely),u(nx1,ny1,nx1,nelv),w1(ny1,nely),w2(ny1,nely)
      integer e,eg,ex,ey,ez

      ny = ny1*nely
      call rzero(ua,ny)
      call rzero(w1,ny)

      do e=1,nelt
         eg = lglel(e)
         call get_exyz(ex,ey,ez,eg,nelx,nely,nelz)

         do k=1,nz1
         do j=1,ny1
         do i=1,nx1
            zz = (1.-zgm1(j,2))/2.0 ! = 1 for i=1, = 0 for k=nx1
            aa = zz*area(i,k,1,e) + (1-zz)*area(i,k,3,e) ! wgtd jacobian
            w1(j,ey) = w1(j,ey) + aa
            ua(j,ey) = ua(j,ey) + aa*u(i,j,k,e)
         enddo
         enddo
         enddo
      enddo

      call gop(ua,w2,'+  ',ny)
      call gop(w1,w2,'+  ',ny)

      do i=1,ny
         ua(i,1) = ua(i,1) / w1(i,1) ! normalize
      enddo

      return
      end

c-----------------------------------------------------------------------
c>    compute the Reynolds stresses as a function of y
      subroutine slf_channel_compute_stresses(time,vx,vy,vz)
      implicit none
c     interface
      include 'SIZE'
      include 'lib/shared_data.f'
      include 'lib/dysmg_data.f'
      include 'lib/filtering_data.f'
      real time,vx(lx1,ly1,lz1,lelv),vy(lx1,ly1,lz1,lelv),
     + vz(lx1,ly1,lz1,lelv)
c     implementation
      integer i,n
      real atime,timel,dtime
      real alpha,beta

      integer iclld ! timestep
      save    iclld
      data    iclld  /0/
      save atime,timel

      n = nx1*ny1*nz1*nelv

      if (iclld .eq. 0) then
        call rzero(uavg,n)
        call rzero(urms,n)
        call rzero(vrms,n)
        call rzero(wrms,n)
        call rzero(uvms,n)
        call rzero(uavg_dconv,n)
        call rzero(urms_dconv,n)
        call rzero(vrms_dconv,n)
        call rzero(wrms_dconv,n)
        call rzero(uvms_dconv,n)
        call rzero(tijavg_xx,n)
        call rzero(tijavg_yy,n)
        call rzero(tijavg_zz,n)
        call rzero(tijavg_xy,n)
        atime = 0
        timel = time
        iclld = 1
      endif

      dtime = time - timel
      atime = atime + dtime

      if (atime .ne. 0.0 .and. dtime .ne. 0.0) then
         beta      = dtime/atime
         alpha     = 1.0-beta

c        compute the time averages
         call avg1(uavg,vx,alpha,beta,n,'uavg',.false.)
         call avg2(urms,vx,alpha,beta,n,'urms',.false.)
         call avg2(vrms,vy,alpha,beta,n,'vrms',.false.)
         call avg2(wrms,vz,alpha,beta,n,'wrms',.false.)
         call avg3(uvms,vx,vy,alpha,beta,n,'uvmm',.false.)

         if (iftles) then
            call avg1(uavg_dconv,uxdconv,alpha,beta,n,'uavg',.false.)
            call avg2(urms_dconv,uxdconv,alpha,beta,n,'urms',.false.)
            call avg2(vrms_dconv,uydconv,alpha,beta,n,'vrms',.false.)
            call avg2(wrms_dconv,uzdconv,alpha,beta,n,'wrms',.false.)
            call avg3(uvms_dconv,uxdconv,uydconv,alpha,beta,n,
     + 'uvmm',.false.)
            call avg1(tijavg_xx,tij(:,:,:,:,1,1),alpha,beta,n,
     + 'tijxx',.false.)
            call avg1(tijavg_yy,tij(:,:,:,:,2,2),alpha,beta,n,
     + 'tijyy',.false.)
            call avg1(tijavg_zz,tij(:,:,:,:,3,3),alpha,beta,n,
     + 'tijzz',.false.)
            call avg1(tijavg_xy,tij(:,:,:,:,1,2),alpha,beta,n,
     + 'tijxy',.false.)
         endif

c        average over statistical homogeneous directions (x-z plane)
         call slf_channel_compute_profiles(iftles)

c        compute the Reynolds stresses
         do i=1,ly1*lely
            Ruu(i) = urms_y(i) - uavg_y(i)**2
            Rvv(i) = vrms_y(i)
            Rww(i) = wrms_y(i)
            Ruv(i) = uvms_y(i)
         enddo
         if (iftles) then
            do i=1,ly1*lely
               Ruu_dconv(i) = urms_dconv_y(i) - uavg_dconv_y(i)**2
               Rvv_dconv(i) = vrms_dconv_y(i)
               Rww_dconv(i) = wrms_dconv_y(i)
               Ruv_dconv(i) = uvms_dconv_y(i)
            enddo
         endif

      endif

      return
      end

c-----------------------------------------------------------------------
c>    compute the planar averages of the time averages
      subroutine slf_channel_compute_profiles(iftles)
c     interface
      include 'SIZE'
      include 'TOTAL'
      include 'ZPER'  ! for nelx,nely,nelz
      include 'lib/dysmg_data.f'
      logical iftles
c     implementation
      integer m

      common /plane/  uavg_pl(ly1*lely)
     $             ,  urms_pl(ly1*lely)
     $             ,  vrms_pl(ly1*lely)
     $             ,  wrms_pl(ly1*lely)
     $             ,  uvms_pl(ly1*lely)
     $             ,  yy(ly1*lely)
     $             ,  w1(ly1*lely),w2(ly1*lely)
     $             ,  ffx_avg, dragx_avg

      nelx  = NUMBER_ELEMENTS_X
      nely  = NUMBER_ELEMENTS_Y
      nelz  = NUMBER_ELEMENTS_Z

c     compute the planar averages
      call planar_average_s(uavg_y,uavg,w1,w2)
      call planar_average_s(urms_y,urms,w1,w2)
      call planar_average_s(vrms_y,vrms,w1,w2)
      call planar_average_s(wrms_y,wrms,w1,w2)
      call planar_average_s(uvms_y,uvms,w1,w2)

      if (iftles) then
         call planar_average_s(uavg_dconv_y,uavg_dconv,w1,w2)
         call planar_average_s(urms_dconv_y,urms_dconv,w1,w2)
         call planar_average_s(vrms_dconv_y,vrms_dconv,w1,w2)
         call planar_average_s(wrms_dconv_y,wrms_dconv,w1,w2)
         call planar_average_s(uvms_dconv_y,uvms_dconv,w1,w2)

         call planar_average_s(tijavg_xx_y,tijavg_xx,w1,w2)
         call planar_average_s(tijavg_yy_y,tijavg_yy,w1,w2)
         call planar_average_s(tijavg_zz_y,tijavg_zz,w1,w2)
         call planar_average_s(tijavg_xy_y,tijavg_xy,w1,w2)
      endif

c     average over half the channel height
c      m = ny1*nely
c      do i=1,ny1*nely/2
c         uavg_pl(i) = 0.5 * (uavg_pl(i) + uavg_pl(m-i+1))
c         urms_pl(i) = 0.5 * (urms_pl(i) + urms_pl(m-i+1))
c         vrms_pl(i) = 0.5 * (vrms_pl(i) + vrms_pl(m-i+1))
c         wrms_pl(i) = 0.5 * (wrms_pl(i) + wrms_pl(m-i+1))
c      enddo

      return
      end

c-----------------------------------------------------------------------
c>    write out volume averaged velocity, and turbulent statistics
c>    and the average velocity profile u(y), as well as the profiles
c>    of the Reynolds stresses
      subroutine slf_channel_write(time)
      implicit none
c     interface
      include 'SIZE'
      include 'lib/shared_data.f'
      include 'lib/dysmg_data.f'
c     implementation
      integer i,ios
      real time
      character*200 fname

      integer iclld ! timestep
      save    iclld
      data    iclld  /0/

c     write out channel statistics
      write(fname,71)
      if (iclld .eq. 0) then
         open(unit=1,file=fname,status='unknown',iostat=ios)
         if (ios .eq. 0) close(1,status='delete')
         open(unit=1,file=fname,status="new",action="write")
      else
         open(unit=1,file=fname,status="old",position="append",
     + action="write")
         write(1,72) iclld, time
         write(1,73) vxv, vyv, vzv, kin, eps, velrms,
     + teddy, leddy, teta, leta, tauw, utau, Retau
         close(1)
      endif

c     write out the average velocity profile u(y)
      write(fname,74)
      if (iclld .eq. 0) then
         open(unit=2,file=fname,status='unknown',iostat=ios)
         if (ios .eq. 0) close(2,status='delete')
         open(unit=2,file=fname,status="new",action="write")
      else
         open(unit=2,file=fname,status="old",position="append",
     + action="write")
         write(2,72) iclld, time
         do i=1,lely*ly1
            write(2,75) ycoords(i),uprof(i)
         enddo
         close(2)
      endif

c     write out the Reynolds stresses as a function of y
      write(fname,76)
      if (iclld .eq. 0) then
         open(unit=3,file=fname,status='unknown',iostat=ios)
         if (ios .eq. 0) close(3,status='delete')
         open(unit=3,file=fname,status="new",action="write")
      else
         if (time .ge. 10.0) then
            open(unit=3,file=fname,status="old",position="append",
     + action="write")
            write(3,72) iclld, time
            if (iftles) then
               do i=1,lely*ly1
                  write(3,77) ycoords(i),Ruu(i),Rvv(i),Rww(i),Ruv(i),
     + Ruu_dconv(i),Rvv_dconv(i),Rww_dconv(i),Ruv_dconv(i),
     + tijavg_xx_y(i),tijavg_yy_y(i),tijavg_zz_y(i),tijavg_xy_y(i)
               enddo
            else
               do i=1,lely*ly1
                  write(3,78) ycoords(i),Ruu(i),Rvv(i),Rww(i),Ruv(i)
               enddo
            endif
            close(3)
         endif
      endif

      iclld = iclld + 1

   71 format('out/fluid/slf_m.txt')
   72 format(I10,E17.9)
   73 format(13E17.9)
   74 format('out/fluid/profile.txt')
   75 format(2E17.9)
   76 format('out/fluid/stresses.txt')
   77 format(13E17.9)
   78 format(5E17.9)

      return
      end

c-----------------------------------------------------------------------
c>    write out the time series u(t) at each time step after
c>    istep >= 0.5*nsteps
      subroutine slf_timeseries(time,xm1,ym1,zm1,vx,vy,vz)
      implicit none
c     interface
      include 'SIZE'
      include 'lib/shared_data.f'
c     implementation
      integer ios
      real time,xm1(lx1,ly1,lz1,lelt),ym1(lx1,ly1,lz1,lelt),
     + zm1(lx1,ly1,lz1,lelt),vx(lx1,ly1,lz1,lelv),vy(lx1,ly1,lz1,lelv),
     + vz(lx1,ly1,lz1,lelv)
      character*200 fname

      integer iclld ! timestep
      save    iclld
      data    iclld  /0/

c     write out channel statistics
      write(fname,76)
      if (iclld .eq. 0) then
         open(unit=1,file=fname,status='unknown',iostat=ios)
         if (ios .eq. 0) close(1,status='delete')
         open(unit=1,file=fname,status="new",action="write")
         write(1,77) time,xm1(4,4,4,4),ym1(4,4,4,4),zm1(4,4,4,4)
      else
         open(unit=1,file=fname,status="old",position="append",
     + action="write")
         write(1,77) time,vx(4,4,4,4),vy(4,4,4,4),vz(4,4,4,4)
         close(1)
      endif

      write(fname,78)
      if (iclld .eq. 0) then
         open(unit=1,file=fname,status='unknown',iostat=ios)
         if (ios .eq. 0) close(1,status='delete')
         open(unit=1,file=fname,status="new",action="write")
         write(1,77) time,xm1(5,4,4,4),ym1(5,4,4,4),zm1(5,4,4,4)
      else
         open(unit=1,file=fname,status="old",position="append",
     + action="write")
         write(1,77) time,vx(5,4,4,4),vy(5,4,4,4),vz(5,4,4,4)
         close(1)
      endif

      iclld = iclld + 1

   76 format('out/fluid/timeseries.txt')
   77 format(4E17.9)
   78 format('out/fluid/timeseries2.txt')

      return
      end

c-----------------------------------------------------------------------
c>    write out all velocity fields at each time step at one particular
c>    grid point until istep 0.05*nsteps to validate the approximate
c>    deconvolution regularization
      subroutine slf_ADRvalidation(time,xm1,ym1,zm1,vx,vy,vz)
      implicit none
c     interface
      include 'SIZE'
      include 'lib/shared_data.f'
      include 'lib/filtering_data.f'
c     implementation
      integer ios
      real time,xm1(lx1,ly1,lz1,lelt),ym1(lx1,ly1,lz1,lelt),
     + zm1(lx1,ly1,lz1,lelt),vx(lx1,ly1,lz1,lelv),vy(lx1,ly1,lz1,lelv),
     + vz(lx1,ly1,lz1,lelv)
      character*200 fname

      integer iclld ! timestep
      save    iclld
      data    iclld  /0/

c     write out channel statistics
      write(fname,76)
      if (iclld .eq. 0) then
         open(unit=1,file=fname,status='unknown',iostat=ios)
         if (ios .eq. 0) close(1,status='delete')
         open(unit=1,file=fname,status="new",action="write")
         write(1,77) time,xm1(4,4,4,4),ym1(4,4,4,4),zm1(4,4,4,4)
      endif
      open(unit=1,file=fname,status="old",position="append",
     + action="write")
         write(1,77) time,uxdconv(4,4,4,4),vx(4,4,4,4),u2(1,4,4,4,4),
     + u3(1,4,4,4,4),u4(1,4,4,4,4),u5(1,4,4,4,4),
     + wxuf(4,4,4,4),wxbar(4,4,4,4),dvdt(4,4,4,4,1)
         close(1)

      iclld = iclld + 1

   76 format('out/fluid/ADRval.txt')
   77 format(10E17.9)

      return
      end
