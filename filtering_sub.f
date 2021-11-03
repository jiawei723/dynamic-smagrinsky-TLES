c-----------------------------------------------------------------------
c>    initialize the filtering data
      subroutine filter_init(vx,vy,vz,dt)
      implicit none
c     interface
      include 'SIZE'
      include 'lib/shared_data.f'
      include 'lib/filtering_data.f'
      include 'lib/fluid_data.f'
      real vx(lx1,ly1,lz1,lelv),vy(lx1,ly1,lz1,lelv),
     + vz(lx1,ly1,lz1,lelv) 
      real dt
c     implementation
      integer h,i,j,k,l,m,n

c     initialize filter width
      if (iframp .and. (Tf .gt. 3.0)) then
         FTS=3.0*dt
      else
         FTS=Tf*dt
      endif
      FTS2=Tf2*dt

c     set the stress tensors and velocity derivatives to zero
      call rzero(tij_old,lx1*ly1*lz1*lelv*3*3)
      call rzero(tij,lx1*ly1*lz1*lelv*3*3)
      call rzero(dvdt,lx1*ly1*lz1*lelv*3)
      call rzero(dvdt_old,lx1*ly1*lz1*lelv*3)

c     set the regularization term to zero
      call rzero(freg,lx1*ly1*lz1*lelv*3)

c     set all filtered fields to the filtered velocity
      if (ifadreg) then
         do l = 1,lelv
         do k = 1,lz1
         do j = 1,ly1
         do i = 1,lx1 
            wxuf(i,j,k,l) = vx(i,j,k,l)
            wyuf(i,j,k,l) = vy(i,j,k,l)
            wzuf(i,j,k,l) = vz(i,j,k,l)
            wxbar(i,j,j,l) = vx(i,j,k,l)
            wybar(i,j,j,l) = vy(i,j,k,l)
            wzbar(i,j,j,l) = vz(i,j,k,l)
            u2(1,i,j,k,l) = vx(i,j,k,l)
            u2(2,i,j,k,l) = vy(i,j,k,l)
            u2(3,i,j,k,l) = vz(i,j,k,l)
            do h = 1,3
               u3(h,i,j,k,l) = u2(h,i,j,k,l)
               u4(h,i,j,k,l) = u2(h,i,j,k,l)
               u5(h,i,j,k,l) = u2(h,i,j,k,l)
            enddo
         enddo
         enddo
         enddo
         enddo
c     set qbar equal to the filtered velocity
      else
         do l = 1,lelv
         do k = 1,lz1
         do j = 1,ly1
         do i = 1,lx1 
            qbar(1,i,j,k,l) = vx(i,j,k,l)
            qbar(2,i,j,k,l) = vy(i,j,k,l)
            qbar(3,i,j,k,l) = vz(i,j,k,l)
         enddo
         enddo
         enddo
         enddo
      endif
      

      if (ifrestart) then
c        initialize velocity of with the first restart file
         call restart(1)
         do l = 1,lelv
         do k = 1,lz1
         do j = 1,ly1
         do i = 1,lx1 
            vxtemp(i,j,k,l) = vx(i,j,k,l)
            vytemp(i,j,k,l) = vy(i,j,k,l)
            vztemp(i,j,k,l) = vz(i,j,k,l)
         enddo
         enddo
         enddo
         enddo

c        compute the initial time derivative of the velocity using EE
         call restart(2)
         do l = 1,lelv
         do k = 1,lz1
         do j = 1,ly1
         do i = 1,lx1
            dvdt(i,j,k,l,1) = (vx(i,j,k,l) - vxtemp(i,j,k,l))/redt
            dvdt(i,j,k,l,2) = (vy(i,j,k,l) - vytemp(i,j,k,l))/redt 
            dvdt(i,j,k,l,3) = (vz(i,j,k,l) - vztemp(i,j,k,l))/redt
         enddo
         enddo
         enddo
         enddo 

c        set the initial velocity field back to the first restart file
         ! call restart(1)
     
         if (ifretles) then
c        approximate initial tij as stationary value / velocity is
c        already filtered         
            do n = 1,3 ! j
            do m = 1,3 ! i
            do l = 1,lelv
            do k = 1,lz1
            do j = 1,ly1
            do i = 1,lx1
               tij(i,j,k,l,m,n) = FTS*FTS*dvdt(i,j,k,l,m)*
     + dvdt(i,j,k,l,n)
            enddo
            enddo
            enddo
            enddo
            enddo
            enddo
         else
c        filter the velocity field (assuming d<u>/dt=du/dt)
            do l = 1,lelv
            do k = 1,lz1
            do j = 1,ly1
            do i = 1,lx1
               vx(i,j,k,l) = vx(i,j,k,l) - FTS*dvdt(i,j,k,l,1)
               vy(i,j,k,l) = vy(i,j,k,l) - FTS*dvdt(i,j,k,l,2)
               vz(i,j,k,l) = vz(i,j,k,l) - FTS*dvdt(i,j,k,l,3) 
            enddo
            enddo
            enddo
            enddo
         endif
      endif

      if (ifdudtbd3) then
         do l = 1,lelv
         do k = 1,lz1
         do j = 1,ly1
         do i = 1,lx1
            vxold(i,j,k,l)=vx(i,j,k,l)
            vyold(i,j,k,l)=vy(i,j,k,l)
            vzold(i,j,k,l)=vz(i,j,k,l)
         enddo
         enddo
         enddo
         enddo
      else
         call rzero(vxold,lx1*ly1*lz1*lelv)
         call rzero(vyold,lx1*ly1*lz1*lelv)
         call rzero(vzold,lx1*ly1*lz1*lelv)
      endif

      return
      end

c-----------------------------------------------------------------------
c>    restart TLES
      subroutine filter_restart(vx,vy,vz)
      implicit none
c     interface
      include 'SIZE'
      include 'lib/shared_data.f'
      include 'lib/filtering_data.f'
      include 'lib/fluid_data.f'
      real vx(lx1,ly1,lz1,lelv),vy(lx1,ly1,lz1,lelv),
     + vz(lx1,ly1,lz1,lelv)
c     implementation
      integer i,j,k,l
      character*200 fname

c     set the previous stress tensors to zero
      call rzero(tij_old,3*3*lx1*ly1*lz1*lelv)

c     read in the data
      write(fname,71) nid
      open(unit=1,file=fname)
      read(1,72) FTS ! filter width
      do l = 1,lelv
      do k = 1,lz1
      do j = 1,ly1
      do i = 1,lx1
         read (1,73)
     + dvdt(i,j,k,l,1),dvdt(i,j,k,l,2),dvdt(i,j,k,l,3),
     + tij(i,j,k,l,1,1),tij(i,j,k,l,2,1),tij(i,j,k,l,3,1),
     + tij(i,j,k,l,1,2),tij(i,j,k,l,2,2),tij(i,j,k,l,3,2),
     + tij(i,j,k,l,1,3),tij(i,j,k,l,2,3),tij(i,j,k,l,3,3)
      enddo
      enddo
      enddo
      enddo

      close(1)
   71 format('restart/nid',i3.3,'.txt')
   72 format(E17.9)
   73 format(12E17.9)

      if (ifdudtbd3) then
         do l = 1,lelv
         do k = 1,lz1
         do j = 1,ly1
         do i = 1,lx1
            vxold(i,j,k,l)=vx(i,j,k,l)
            vyold(i,j,k,l)=vy(i,j,k,l)
            vzold(i,j,k,l)=vz(i,j,k,l)
         enddo
         enddo
         enddo
         enddo
      else
         call rzero(vxold,lx1*ly1*lz1*lelv)
         call rzero(vyold,lx1*ly1*lz1*lelv)
         call rzero(vzold,lx1*ly1*lz1*lelv)
      endif

      return
      end

c-----------------------------------------------------------------------
c>    compute the stress components tij at the new time step
      subroutine filter_stress(istep,vx,vxlag,vy,vylag,vz,vzlag,
     + dt,torder)
      implicit none
c     interface
      include 'SIZE'
      include 'lib/shared_data.f'
      include 'lib/filtering_data.f'
      include 'lib/fluid_data.f'
      integer istep,torder
      real vx(lx1,ly1,lz1,lelv),vy(lx1,ly1,lz1,lelv),
     + vz(lx1,ly1,lz1,lelv) 
      real vxlag(lx1,ly1,lz1,lelv,torder), 
     + vylag(lx1,ly1,lz1,lelv,torder),vzlag(lx1,ly1,lz1,lelv,torder)
      real dt
c     implementation 
      real kzero,kone,ktwo,kthree
      real dvdt_halfm,dvdt_halfn
      integer i,j,k,l,m,n,o

c     adapt the filter width if necessary
      if (iframp .and. (istep .le. 50.0) .and. (Tf .gt. 3.0)) then
         FTS=3.0*dt+(Tf-3.0)*(istep/50.0)*dt
      endif

c     compute the temporal velocity gradients using 3rd order BD
      if (ifdudtbd3 .and. (istep .gt. 2)) then
         do l = 1,lelv
         do k = 1,lz1
         do j = 1,ly1
         do i = 1,lx1
            dvdt_old(i,j,k,l,1) = dvdt(i,j,k,l,1)
            dvdt_old(i,j,k,l,2) = dvdt(i,j,k,l,2)
            dvdt_old(i,j,k,l,3) = dvdt(i,j,k,l,3)
            dvdt(i,j,k,l,1) = (11.0*vx(i,j,k,l) - 
     + 18.0*vxlag(i,j,k,l,1) + 9.0*vxlag(i,j,k,l,2) - 
     + 2.0*vxold(i,j,k,l))/(6.0*dt)
            dvdt(i,j,k,l,2) = (11.0*vy(i,j,k,l) - 
     + 18.0*vylag(i,j,k,l,1) + 9.0*vylag(i,j,k,l,2) -
     + 2.0*vyold(i,j,k,l))/(6.0*dt)
            dvdt(i,j,k,l,3) = (11.0*vz(i,j,k,l) -
     + 18.0*vzlag(i,j,k,l,1) + 9.0*vzlag(i,j,k,l,2) -
     + 2.0*vzold(i,j,k,l))/(6.0*dt)
c           copy the velocity of t^(n-2) to vold
            vxold(i,j,k,l)=vxlag(i,j,k,l,2)
            vyold(i,j,k,l)=vylag(i,j,k,l,2)
            vzold(i,j,k,l)=vzlag(i,j,k,l,2)
         enddo
         enddo
         enddo
         enddo
c     compute the temporal velocity gradients using 2nd order BD
      elseif ((ifdudtbd2 .and. (istep .gt. 1)) .or. 
     + (ifdudtbd3 .and. (istep .eq. 2))) then
         do l = 1,lelv
         do k = 1,lz1
         do j = 1,ly1
         do i = 1,lx1
            dvdt_old(i,j,k,l,1) = dvdt(i,j,k,l,1)
            dvdt_old(i,j,k,l,2) = dvdt(i,j,k,l,2)
            dvdt_old(i,j,k,l,3) = dvdt(i,j,k,l,3)
            dvdt(i,j,k,l,1) = (3.0*vx(i,j,k,l) - 4.0*vxlag(i,j,k,l,1) + 
     + vxlag(i,j,k,l,2))/(2.0*dt)
            dvdt(i,j,k,l,2) = (3.0*vy(i,j,k,l) - 4.0*vylag(i,j,k,l,1) + 
     + vylag(i,j,k,l,2))/(2.0*dt)
            dvdt(i,j,k,l,3) = (3.0*vz(i,j,k,l) - 4.0*vzlag(i,j,k,l,1) +
     + vzlag(i,j,k,l,2))/(2.0*dt)
         enddo
         enddo
         enddo
         enddo
c     compute the temporal velocity gradients 1st order BD
      elseif (istep .gt. 0) then
         do l = 1,lelv
         do k = 1,lz1
         do j = 1,ly1
         do i = 1,lx1
            dvdt_old(i,j,k,l,1) = dvdt(i,j,k,l,1)
            dvdt_old(i,j,k,l,2) = dvdt(i,j,k,l,2)
            dvdt_old(i,j,k,l,3) = dvdt(i,j,k,l,3)
            dvdt(i,j,k,l,1) = (vx(i,j,k,l) - vxlag(i,j,k,l,1))/dt
            dvdt(i,j,k,l,2) = (vy(i,j,k,l) - vylag(i,j,k,l,1))/dt
            dvdt(i,j,k,l,3) = (vz(i,j,k,l) - vzlag(i,j,k,l,1))/dt
         enddo
         enddo
         enddo
         enddo
      endif

c     integrate the stress tensor tij with RK4
      if (ifrkf .and. istep .gt. 0) then
         do n = 1,3 ! j of the tensor tij
         do m = 1,3 ! i of the tensor tij
         do l = 1,lelv
         do k = 1,lz1
         do j = 1,ly1
         do i = 1,lx1
            tij_old(i,j,k,l,m,n) = tij(i,j,k,l,m,n)
            dvdt_halfm=0.5*(dvdt_old(i,j,k,l,m)+dvdt(i,j,k,l,m))
            dvdt_halfn=0.5*(dvdt_old(i,j,k,l,n)+dvdt(i,j,k,l,n))
            kzero = - tij_old(i,j,k,l,m,n)/FTS +
     + FTS*dvdt_old(i,j,k,l,m)*dvdt_old(i,j,k,l,n)
            kone = - (tij_old(i,j,k,l,m,n)+0.5*dt*kzero)/FTS +
     + FTS*dvdt_halfm*dvdt_halfn
            ktwo = - (tij_old(i,j,k,l,m,n)+0.5*dt*kone)/FTS +
     + FTS*dvdt_halfm*dvdt_halfn
            kthree = - (tij_old(i,j,k,l,m,n)+dt*ktwo)/FTS +
     + FTS*dvdt(i,j,k,l,m)*dvdt(i,j,k,l,n)
            tij(i,j,k,l,m,n) = tij_old(i,j,k,l,m,n) +
     + dt*(kzero+2.0*kone+2.0*ktwo+kthree)/6.0      
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
c     integrate tij using Euler explicit
      elseif (ifee .and. istep .gt. 0) then
         do n = 1,3 ! j
         do m = 1,3 ! i
         do l = 1,lelv
         do k = 1,lz1
         do j = 1,ly1
         do i = 1,lx1
            tij_old(i,j,k,l,m,n) = tij(i,j,k,l,m,n)
            kzero = - tij_old(i,j,k,l,m,n)/FTS +
     + FTS*dvdt(i,j,k,l,m)*dvdt(i,j,k,l,n) + 
     + diff*tij_lap(i,j,k,l,m,n)
            tij(i,j,k,l,m,n) = tij_old(i,j,k,l,m,n) +
     + dt*kzero
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
c     integrate tij using exact integration      
      elseif (ifei .and. istep .gt. 0) then
         do n = 1,3 ! j
         do m = 1,3 ! i
         do l = 1,lelv
         do k = 1,lz1
         do j = 1,ly1
         do i = 1,lx1
            tij_old(i,j,k,l,m,n) = tij(i,j,k,l,m,n)
            tij(i,j,k,l,m,n) = (tij_old(i,j,k,l,m,n) - FTS*FTS*
     + dvdt(i,j,k,l,m)*dvdt(i,j,k,l,n))*exp(-dt/FTS) + FTS*FTS*
     + dvdt(i,j,k,l,m)*dvdt(i,j,k,l,n)
         enddo
         enddo
         enddo
         enddo
         enddo
         enddo
      endif

      return
      end

c-----------------------------------------------------------------------
c>    overwrite the velocity to an analytical value for 
c>    the validation of tij
      subroutine filter_valid(xm1,ym1,zm1,vx,vy,vz,istep,dt,pi)
      implicit none
c     interface
      include 'SIZE'
      real xm1(lx1,ly1,lz1,lelt),ym1(lx1,ly1,lz1,lelt),
     + zm1(lx1,ly1,lz1,lelt),vx(lx1,ly1,lz1,lelv),vy(lx1,ly1,lz1,lelv),
     + vz(lx1,ly1,lz1,lelv),dt,pi
      integer istep
c     implementation
      real x,y,z
      integer i,j,k,l,m

      do m = 1,nelv
      do l = 1,lz1
      do k = 1,ly1
      do j = 1,lx1
         x = xm1(j,k,l,m)
         y = ym1(j,k,l,m)
         z = zm1(j,k,l,m)
         vx(j,k,l,m) = 10.0*sin(4.0*x/2.0)*
     + sin(2.0*pi*dt*istep/(25.0*dt))
         vy(j,k,l,m) = 0.0
         vz(j,k,l,m) = 0.0
         ! vx(j,k,l,m) = 1.0*sin(pi*y)*sin(x/2.0)*dt*istep
      enddo
      enddo
      enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
c>    compute divergence of tij
      subroutine filter_div(binvm1)
      implicit none
c     interface
      include 'SIZE'
      include 'lib/filtering_data.f'
      real binvm1(lx1,ly1,lz1,lelt)
c     implementation
      integer ntot

      ntot = nx1*ny1*nz1*nelv

c     compute d(tau_x,j)/(dxj)
      call opdiv(tij_xdiv,tij(:,:,:,:,1,1),tij(:,:,:,:,1,2),
     + tij(:,:,:,:,1,3))
      call dssum(tij_xdiv,nx1,ny1,nz1)
      call col2(tij_xdiv,binvm1,ntot)

c     compute d(tau_y,j)/(dxj)
      call opdiv(tij_ydiv,tij(:,:,:,:,2,1),tij(:,:,:,:,2,2),
     + tij(:,:,:,:,2,3))
      call dssum(tij_ydiv,nx1,ny1,nz1)
      call col2(tij_ydiv,binvm1,ntot)

c     compute d(tau_z,j)/(dxj)
      call opdiv(tij_zdiv,tij(:,:,:,:,3,1),tij(:,:,:,:,3,2),
     + tij(:,:,:,:,3,3))
      call dssum(tij_zdiv,nx1,ny1,nz1)
      call col2(tij_zdiv,binvm1,ntot)

      return
      end

c-----------------------------------------------------------------------
c>    compute the laplacian of tij for the diffusion term 
      subroutine filter_diffterm(binvm1)
      implicit none
c     interface
      include 'SIZE'
      include 'lib/filtering_data.f'
      real binvm1(lx1,ly1,lz1,lelt)
c     implementation
      real work1(lx1,ly1,lz1,lelv),work2(lx1,ly1,lz1,lelv),
     + work3(lx1,ly1,lz1,lelv)
      integer i,j,k,l,m,n

      do n = 1,3 ! j
      do m = 1,3 ! i
         call opgrad(work1,work2,work3,tij(:,:,:,:,n,m))
         call opdssum(work1,work2,work3)
         call opcolv(work1,work2,work3,binvm1)
         !call gradm1(work1,work2,work3,tij(:,:,:,:,1,1))
         call opdiv(work1,work1,work2,work3)
         call opdssum(work1,work2,work3)
         call opcolv(work1,work2,work3,binvm1)
c        copy to tij_lap
         do l = 1,lelv
         do k = 1,lz1
         do j = 1,ly1
         do i = 1,lx1
            tij_lap(i,j,k,l,m,n) = work1(i,j,k,l)
         enddo
         enddo
         enddo
         enddo
      enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
c>    deconvolute velocity field
      subroutine filter_deconv(vx,vy,vz)
      implicit none
c     interface
      include 'SIZE'
      include 'lib/filtering_data.f'
      real vx(lx1,ly1,lz1,lelv),vy(lx1,ly1,lz1,lelv),
     + vz(lx1,ly1,lz1,lelv) 
c     implementation
      integer i,j,k,l

      do l = 1,lelv
      do k = 1,lz1
      do j = 1,ly1
      do i = 1,lx1 
         uxdconv(i,j,k,l) = vx(i,j,k,l) + FTS*dvdt(i,j,k,l,1)
         uydconv(i,j,k,l) = vy(i,j,k,l) + FTS*dvdt(i,j,k,l,2)
         uzdconv(i,j,k,l) = vz(i,j,k,l) + FTS*dvdt(i,j,k,l,3)
      enddo
      enddo
      enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
c>    compute the regularization term based on Pruetts' Paper
c>    using his approx. deconvolution
      subroutine filter_ADregularization(istep,vx,vy,vz,dt)
      implicit none
c     interface
      include 'SIZE'
      include 'lib/shared_data.f'
      include 'lib/filtering_data.f'
      integer istep
      real vx(lx1,ly1,lz1,lelv),vy(lx1,ly1,lz1,lelv),
     + vz(lx1,ly1,lz1,lelv),dt
c     implementation
      integer h,i,j,k,l
      real chi_ramp

c     adapt the filter width if necessary
      if (iframp .and. (istep .le. 50.0)) then
         if (istep .lt. 20) then
            chi_ramp = 0.0
         else
            chi_ramp = (istep-20.0)/30.0*chi
         endif
      else
         chi_ramp = chi
      endif

      if (istep .gt. 0) then
         do l = 1,lelv
         do k = 1,lz1
         do j = 1,ly1
         do i = 1,lx1
c           integrate fields at new time step       
            u2(1,i,j,k,l) = u2(1,i,j,k,l) + 
     + dt*(vx(i,j,k,l) - u2(1,i,j,k,l))/FTS
            u2(2,i,j,k,l) = u2(2,i,j,k,l) + 
     + dt*(vy(i,j,k,l) - u2(2,i,j,k,l))/FTS
            u2(3,i,j,k,l) = u2(3,i,j,k,l) + 
     + dt*(vz(i,j,k,l) - u2(3,i,j,k,l))/FTS
            do h=1,3
               u3(h,i,j,k,l) = u3(h,i,j,k,l) + 
     + dt*(u2(h,i,j,k,l) - u3(h,i,j,k,l))/FTS
               u4(h,i,j,k,l) = u4(h,i,j,k,l) + 
     + dt*(u3(h,i,j,k,l) - u4(h,i,j,k,l))/FTS
               u5(h,i,j,k,l) = u5(h,i,j,k,l) + 
     + dt*(u4(h,i,j,k,l) - u5(h,i,j,k,l))/FTS
            enddo
c           approx. deconvolute
            if (regorder .eq. 3) then
               wxuf(i,j,k,l) = DC0*vx(i,j,k,l) + DC1*u2(1,i,j,k,l) +
     + DC2*u3(1,i,j,k,l) + DC3*u4(1,i,j,k,l)
               wyuf(i,j,k,l) = DC0*vy(i,j,k,l) + DC1*u2(2,i,j,k,l) +
     + DC2*u3(2,i,j,k,l) + DC3*u4(2,i,j,k,l)
               wzuf(i,j,k,l) = DC0*vz(i,j,k,l) + DC1*u2(3,i,j,k,l) +
     + DC2*u3(3,i,j,k,l) + DC3*u4(3,i,j,k,l)
            elseif (regorder .eq. 4) then
               wxuf(i,j,k,l) = DC0_q4*vx(i,j,k,l) + 
     + DC1_q4*u2(1,i,j,k,l) + DC2_q4*u3(1,i,j,k,l) + 
     + DC3_q4*u4(1,i,j,k,l) + DC4_q4*u5(1,i,j,k,l)
               wyuf(i,j,k,l) = DC0_q4*vy(i,j,k,l) + 
     + DC1_q4*u2(2,i,j,k,l) + DC2_q4*u3(2,i,j,k,l) + 
     + DC3_q4*u4(2,i,j,k,l) + DC4_q4*u5(2,i,j,k,l)
               wzuf(i,j,k,l) = DC0_q4*vz(i,j,k,l) +
     + DC1_q4*u2(3,i,j,k,l) + DC2_q4*u3(3,i,j,k,l) + 
     + DC3_q4*u4(3,i,j,k,l) + DC4_q4*u5(3,i,j,k,l)
            endif
c           filter w
            wxbar(i,j,k,l) = wxbar(i,j,k,l) + 
     + dt*(wxuf(i,j,k,l) - wxbar(i,j,k,l))/FTS
            wybar(i,j,k,l) = wybar(i,j,k,l) + 
     + dt*(wyuf(i,j,k,l) - wybar(i,j,k,l))/FTS
            wzbar(i,j,k,l) = wzbar(i,j,k,l) + 
     + dt*(wzuf(i,j,k,l) - wzbar(i,j,k,l))/FTS
c           compute the regularization term
            freg(1,i,j,k,l) = chi_ramp*(vx(i,j,k,l) - wxbar(i,j,k,l))
            freg(2,i,j,k,l) = chi_ramp*(vy(i,j,k,l) - wybar(i,j,k,l))
            freg(3,i,j,k,l) = chi_ramp*(vz(i,j,k,l) - wzbar(i,j,k,l))
         enddo
         enddo
         enddo
         enddo
      endif

      return
      end

c-----------------------------------------------------------------------
c>    compute the regularization term and integrate the filtered 
c>    velocity qbar
      subroutine filter_regularization(istep,vx,vy,vz,dt)
      implicit none
c     interface
      include 'SIZE'
      include 'lib/shared_data.f'
      include 'lib/filtering_data.f'
      integer istep
      real vx(lx1,ly1,lz1,lelv),vy(lx1,ly1,lz1,lelv),
     + vz(lx1,ly1,lz1,lelv),dt
c     implementation
      integer i,j,k,l
      real chi_ramp

c     adapt the filter width if necessary
      if (iframp .and. (istep .le. 50.0)) then
         if (istep .lt. 20) then
            chi_ramp = 0.0
         else
            chi_ramp = (istep-20.0)/30.0*chi
         endif
      else
         chi_ramp = chi
      endif

      if (istep .gt. 0) then
         do l = 1,lelv
         do k = 1,lz1
         do j = 1,ly1
         do i = 1,lx1
            qbar(1,i,j,k,l) = qbar(1,i,j,k,l) + 
     + dt*(uxdconv(i,j,k,l) - qbar(1,i,j,k,l))/FTS2
            qbar(2,i,j,k,l) = qbar(2,i,j,k,l) + 
     + dt*(uydconv(i,j,k,l) - qbar(2,i,j,k,l))/FTS2
            qbar(3,i,j,k,l) = qbar(3,i,j,k,l) + 
     + dt*(uzdconv(i,j,k,l) - qbar(3,i,j,k,l))/FTS2
            freg(1,i,j,k,l) = chi_ramp*
     + (uxdconv(i,j,k,l) - qbar(1,i,j,k,l))
            freg(2,i,j,k,l) = chi_ramp*
     + (uydconv(i,j,k,l) - qbar(2,i,j,k,l))
            freg(3,i,j,k,l) = chi_ramp*
     + (uzdconv(i,j,k,l) - qbar(3,i,j,k,l))
         enddo
         enddo
         enddo
         enddo
      endif

      return
      end

c-----------------------------------------------------------------------
c>    restart TLES
      subroutine filter_write_restart()
      implicit none
c     interface
      include 'SIZE'
      include 'lib/filtering_data.f'
c     implementation
      integer i,j,k,l,ios
      character*200 fname

c     read in the data
      write(fname,71) nid
      open(unit=1,file=fname,status='unknown',iostat=ios)
      if (ios .eq. 0) close(1,status='delete')
      open(unit=1,file=fname)
      write(1,72) FTS ! filter width
      do l = 1,lelv
      do k = 1,lz1
      do j = 1,ly1
      do i = 1,lx1
         write (1,73)
     + dvdt(i,j,k,l,1),dvdt(i,j,k,l,2),dvdt(i,j,k,l,3),
     + tij(i,j,k,l,1,1),tij(i,j,k,l,2,1),tij(i,j,k,l,3,1),
     + tij(i,j,k,l,1,2),tij(i,j,k,l,2,2),tij(i,j,k,l,3,2),
     + tij(i,j,k,l,1,3),tij(i,j,k,l,2,3),tij(i,j,k,l,3,3)
      enddo
      enddo
      enddo
      enddo

      close(1)
   71 format('out/restart/nid',i3.3,'.txt')
   72 format(E17.9)
   73 format(12E17.9)

      return
      end

c-----------------------------------------------------------------------
c>    diagnosis of the TLES closure - only x components of one node
      subroutine filter_diag(xm1,ym1,zm1,vx,vy,vz,dt)
      implicit none
c     interface
      include 'SIZE'
      include 'lib/shared_data.f'
      include 'lib/filtering_data.f'
      real xm1(lx1,ly1,lz1,lelt),ym1(lx1,ly1,lz1,lelt),
     + zm1(lx1,ly1,lz1,lelt),vx(lx1,ly1,lz1,lelv),vy(lx1,ly1,lz1,lelv),
     + vz(lx1,ly1,lz1,lelv),dt
c     implementation
      integer i,j,k,l,ios
      real x,y,z
      real tole
      parameter (tole = 1.0E-04)
      character*200 fname

      integer iclld ! timestep
      save    iclld
      data    iclld  /0/

      write(fname,71)
      if (iclld .eq. 0 .and. nid .eq. 0) then
         open(unit=2,file=fname,status='unknown',iostat=ios)
         if (ios .eq. 0) close(2,status='delete')
         open(unit=2,file=fname,status="new",action="write")
         write(2,72) FTS,FTS2,dt
         close(2)
      endif
   71 format('out/fluid/diag.txt')
   72 format(3E17.9)

      do l = 1,lelv
      do k = 1,lz1
      do j = 1,ly1
      do i = 1,lx1 
         x = xm1(i,j,k,l)
         y = ym1(i,j,k,l)
         z = zm1(i,j,k,l)
         if (x .gt. (ldmx/2.0-tole) .and.
     + x .lt. (ldmx/2.0+tole) .and.
     + y .gt. (ldmy/2.0-tole) .and.
     + y .lt. (ldmy/2.0+tole) .and. 
     + z .gt. (ldmz/2.0-tole) .and.
     + z .lt. (ldmz/2.0+tole)) then
            open(unit=2,file=fname,status="old",position="append", 
     + action="write")
            write(2,73) vx(i,j,k,l),dvdt(i,j,k,l,1),uxdconv(i,j,k,l),
     + tij(i,j,k,l,1,1),qbar(1,i,j,k,l),tij_xdiv(i,j,k,l)
            close(2)
            iclld = iclld + 1
   73       format(6E17.9)
            return
         endif
      enddo
      enddo
      enddo
      enddo

      return
      end

