c-----------------------------------------------------------------------
c>    initialize fluid velocity at grid point (x,y,z) with value from 
c>    kinematic simulation 
      subroutine fluid_init(xm1,ym1,zm1,vx,vy,vz,pi,pr)
      implicit none
c     interface
      include 'SIZE'
      include 'lib/fluid_data.f'
      include 'lib/shared_data.f'
      real xm1(lx1,ly1,lz1,lelt),ym1(lx1,ly1,lz1,lelt),
     + zm1(lx1,ly1,lz1,lelt),vx(lx1,ly1,lz1,lelv),vy(lx1,ly1,lz1,lelv),
     + vz(lx1,ly1,lz1,lelv),pi
      real pr(lx2,ly2,lz2,lelv)
c     implementation
      integer i,j,k,l,m  
      real arg,x,y,z,ux,uy,uz 
      real p1,p2,p3,q1,q2,q3,q4 ! fit parameters

      p1=52.17;
      p2=3.731;
      p3=0.0004528;
      q1=-2.006;
      q2=3.869;
      q3=0.1695;
      q4=0.0225;

      do m = 1,nelv
      do l = 1,lz1
      do k = 1,ly1
      do j = 1,lx1
         x = xm1(j,k,l,m)
         y = ym1(j,k,l,m)
         z = zm1(j,k,l,m)
         ux = 0.0
         uy = 0.0
         uz = 0.0
         do i = 1,nksim
            arg = kint(1,i)*x + kint(2,i)*y + kint(3,i)*z
            ux = ux + an(1,i)*cos(arg) + bn(1,i)*sin(arg)
            uy = uy + an(2,i)*cos(arg) + bn(2,i)*sin(arg)
            uz = uz + an(3,i)*cos(arg) + bn(3,i)*sin(arg)
         enddo
         if (ifchannel) then
            !vx(j,k,l,m) = 0.0
            !vy(j,k,l,m) = 0.0
            !vz(j,k,l,m) = 0.0
            !vx(j,k,l,m) = 23.58*(1.0-(y-1.0)**2.0)
            !vy(j,k,l,m) = -2.5*sin(1.0*2.0*pi*y/2.0)*sin(1.0*1.5*z)
            !vz(j,k,l,m) = -2.5*sin(1.0*2.0*pi*y/2.0)*sin(1.0*1.5*z)
c           vz(j,k,l,m) = 0.8*sin(1.0*2.0*pi*y/2.0)*sin(1.0*1.5*z)
c           if (y .gt. 1.0) y=abs(y-2.0)
c           vx(j,k,l,m) = (p1*y**2.0+p2*y+p3)/(y**4.0+q1*y**3.0+
c     + q2*y**2.0+q3*y+q4)
         else if (iftgv) then
            vx(j,k,l,m) = sin(x)*cos(y)*cos(z)
            vy(j,k,l,m) = -cos(x)*sin(y)*cos(z)
            vz(j,k,l,m) = 0.0
            pr(j,k,l,m) = pr(j,k,l,m) +
     + (cos(2.0*x) + cos(2.0*y)) * (cos(2.0*z) + 2.0)/16.0
         else ! 3d box
            vx(j,k,l,m) = 0.0 !ux
            vy(j,k,l,m) = 0.0 !uy
            vz(j,k,l,m) = 0.0 !uz
         endif
      enddo
      enddo
      enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
c>    perform kinematic simulation
c>    taken from P. Weiss
      subroutine fluid_ksim(pi)
      implicit none
c     interface
      include 'SIZE'
      include 'lib/shared_data.f' 
      include 'lib/fluid_data.f' 
      real pi
c     implementation
      integer n1,i,ios 
      real ldsim,l1,eta1,emdl1,kmin1,kdbl1(3),dkn1,kdblm1,rt1,vreal 
      real u1,theta1,gamma1(3),alpha1(3),beta1(3) 
      real abn1,crs1(3),acrs1,kdblt
      character*200 fname

c     set parameters 
      ldsim = ldmx
      l1 = 0.19/(2.0/3.0)**(3.0/2.0)*ldsim ! integral lengthscale 
      eta1 = l1*(3.0/20.0*rlsim**2.0)**(-3.0/4.0) ! Kolmogorov lengths.

      emdl1 = 0.0
      kmin1 = 2.0*pi/ldsim ! minimum wavenumber

c     compute wavenumber vectors 
      call srand(0) ! each processor starts with same random number seed
      do n1 = 1,nksim ! wavenumbers 
c        random direction    
         u1 = 2.0*rand(0) - 1.0
         theta1 = 2.0*pi*rand(0)
         gamma1(1) = sqrt(1.0 - u1**2.0)*cos(theta1) ! in x
         gamma1(2) = sqrt(1.0 - u1**2.0)*sin(theta1) ! in y
         gamma1(3) = u1  ! in z
         rt1 = ((real(n1)-1.0)/(real(nksim)-1.0))
         kdbl1(1) = gamma1(1)*(ldsim/eta1)**rt1 
         kdbl1(2) = gamma1(2)*(ldsim/eta1)**rt1 
         kdbl1(3) = gamma1(3)*(ldsim/eta1)**rt1 
c        discrete wavenumber vector
         kint(1,n1) = anint(kdbl1(1)) ! in x
         kint(2,n1) = anint(kdbl1(2)) ! in y
         kint(3,n1) = anint(kdbl1(3)) ! in z 
         if ((kint(1,n1) .eq. 0.0) .and. (kint(2,n1) .eq. 0.0)
     + .and. (kint(3,n1) .eq. 0.0)) then
            kdblm1 = max(abs(kdbl1(1)),abs(kdbl1(2)),abs(kdbl1(3)))
            do i = 1,3
               if (abs(kdbl1(i)) .eq. kdblm1) then
                  kint(i,n1) = 1
                  kint(i,n1) = sign(kint(i,n1),kdbl1(i))
               endif
            enddo
         endif
         kint(1,n1) = kmin1*kint(1,n1)
         kint(2,n1) = kmin1*kint(2,n1)
         kint(3,n1) = kmin1*kint(3,n1)
c        continuous wavenumber 
         kdbl(n1) = kmin1*sqrt(kdbl1(1)**2.0 + kdbl1(2)**2.0 + 
     + kdbl1(3)**2.0)
      enddo

c     compute Fourier vectors 
      do n1 = 1, nksim ! wavenumbers
         if (n1 .eq. 1) then
            dkn1 = (kdbl(2) - kdbl(1))/2.0
         elseif (n1 .eq. nksim) then
            dkn1 = (kdbl(nksim) - kdbl(nksim-1))/2.0
         else
            dkn1 = (kdbl(n1+1) - kdbl(n1-1))/2.0
         endif    
         kdblt = kdbl(n1)   
         call fluid_emdl(emdl1,kdblt,l1,eta1) 
         abn1 = sqrt(2.0*emdl1*dkn1)
         u1 = 2.0*rand(0) - 1.0
         theta1 = 2.0*pi*rand(0)
         alpha1(1) = sqrt(1.0 - u1**2.0)*cos(theta1)
         alpha1(2) = sqrt(1.0 - u1**2.0)*sin(theta1)
         alpha1(3) = u1
         crs1(1) = kint(2,n1)*alpha1(3) - kint(3,n1)*alpha1(2)      
         crs1(2) = kint(3,n1)*alpha1(1) - kint(1,n1)*alpha1(3)
         crs1(3) = kint(1,n1)*alpha1(2) - kint(2,n1)*alpha1(1)
         acrs1 = sqrt(crs1(1)**2.0 + crs1(2)**2.0 + crs1(3)**2.0)
         if (acrs1 .gt. 0.0) then
            an(1,n1) = crs1(1)*abn1/acrs1
            an(2,n1) = crs1(2)*abn1/acrs1
            an(3,n1) = crs1(3)*abn1/acrs1
         endif
         u1 = 2.0*rand(0) - 1.0
         theta1 = 2.0*pi*rand(0)
         beta1(1) = sqrt(1.0 - u1**2.0)*cos(theta1)
         beta1(2) = sqrt(1.0 - u1**2.0)*sin(theta1)
         beta1(3) = u1
         crs1(1) = kint(2,n1)*beta1(3) - kint(3,n1)*beta1(2)      
         crs1(2) = kint(3,n1)*beta1(1) - kint(1,n1)*beta1(3)
         crs1(3) = kint(1,n1)*beta1(2) - kint(2,n1)*beta1(1)
         acrs1 = sqrt(crs1(1)**2.0 + crs1(2)**2.0 + crs1(3)**2.0)
         if (acrs1 .gt. 0.0) then
            bn(1,n1) = crs1(1)*abn1/acrs1
            bn(2,n1) = crs1(2)*abn1/acrs1
            bn(3,n1) = crs1(3)*abn1/acrs1
         endif
      enddo

      return
      end

c-----------------------------------------------------------------------
c>    evaluate Pope's model spectrum 
c>    "Turbulent Flows", Stephen B. Pope, 2011
c>    taken from P. Weiss
      subroutine fluid_emdl(emdl2,k2,l2,eta2)
      implicit none
c     interface
      include 'SIZE'
      include 'lib/shared_data.f'  
      include 'lib/fluid_data.f'  
      real emdl2,k2,l2,eta2
c     implementation
      real c2,p02,beta2,cl2,ceta2,fl2,feta2

      c2 = 1.5
      p02 = 2.0
      beta2 = 5.2
      cl2 = 3.3659967589175017*(1.0 - tanh(2.104125639020907 - 
     + 0.5840061721497946*log(rlsim)))
      ceta2 = 0.40168569828192224 + 2.675542340609635/
     + (0.6144720501233607 + rlsim)**1.0623906384407826
      fl2 = (k2*l2/((k2*l2)**2.0 + cl2)**0.5)**(5.0/3.0 + p02)
      feta2 = exp(-beta2*(((k2*eta2)**4.0 + ceta2**4.0) 
     + **(1.0/4.0) - ceta2))
      emdl2 =  c2*(mugi/rhogi)**2.0/eta2*(k2*eta2)**(-5.0/3.0)*fl2*feta2

      return
      end

c-----------------------------------------------------------------------
c>    compute the forcing term to obtain a constant volume flow
c>    adapted from Nek examples "turbChannel"
      subroutine set_forcing(vx,dt,istep,bm1,volvm1)
      implicit none
c     interface
      include 'SIZE'
      include 'lib/shared_data.f'
      include 'lib/fluid_data.f'
      real vx(lx1,ly1,lz1,lelv),dt,bm1(lx1,ly1,lz1,lelt),volvm1
      integer istep
c     implementation
      real ubar,utilde,alpha,glsc2
      integer n

      n = nx1*ny1*nz1*nelv
      ubar = glsc2(vx,bm1,n)/volvm1
      utilde = ubar/utarget

      if (istep.eq.0) fvf = 0.32
      if (istep.gt.0) fvf = 0.5*(fvf + fvf/utilde)
      if (istep.gt.5) then
         alpha = 100*abs(utilde - utildeold)/dt
         alpha = min(alpha,0.05)
         if (utilde .gt. 1 .and. utilde .gt. utildeold) then
            fvf = (1.0 - alpha)*fvf
         elseif (utilde .lt. 1 .and. utilde .lt. utildeold) then
            fvf = (1.0 + alpha)*fvf
         endif
      endif
      utildeold = utilde

c     limit the forcing term with ad hoc limits
      fvf = max(fvf,0.00001)
      fvf = min(fvf,5.00000)

      return
      end

