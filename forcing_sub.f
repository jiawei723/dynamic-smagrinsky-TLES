c-----------------------------------------------------------------------
c>    initialization of the low wave number forcing
      subroutine forcing_init(pi)
      implicit none
c     interface
      include 'SIZE'
      include 'lib/shared_data.f' 
      include 'lib/forcing_data.f' 
      real pi
c     implementation
      integer i,j,k     ! loop variables
      integer ix,iy,iz  ! array indices

      call srand(1)     ! initialize random number generator
      call rzero(aFxr,NU*NU*NU)
      call rzero(aFyr,NU*NU*NU)
      call rzero(aFzr,NU*NU*NU)
      call rzero(flwx,lx1*ly1*lz1*lelv)
      call rzero(flwy,lx1*ly1*lz1*lelv)
      call rzero(flwz,lx1*ly1*lz1*lelv)
      call rzero(flwx_old,lx1*ly1*lz1*lelv)
      call rzero(flwy_old,lx1*ly1*lz1*lelv)
      call rzero(flwz_old,lx1*ly1*lz1*lelv)
      do k=1,NU
      do j=1,NU
      do i=1,NU
         bx(i,j,k)=cmplx(0.0,0.0)
         by(i,j,k)=cmplx(0.0,0.0)
         bz(i,j,k)=cmplx(0.0,0.0)
      enddo
      enddo
      enddo
      k0=2.0*pi/ldmx
      kF=sqrt(2.0)*k0
      udx=ldmx/NU

c     test if NU is an even number
      if (mod(NU,2) .eq. 1) then
         write (*,*) "NU has to be an even number!"
         STOP
      endif

c     build the wavenumber vector
      do k=-NU/2,NU/2-1,1
         do j=-NU/2,NU/2-1,1
            do i=-NU/2,NU/2-1,1
               ix=i+NU/2+1
               iy=j+NU/2+1
               iz=k+NU/2+1
               wnvx(ix,iy,iz)=i*k0
               wnvy(ix,iy,iz)=j*k0
               wnvz(ix,iy,iz)=k*k0
               abv(ix,iy,iz)=(i*k0)**2+(j*k0)**2+(k*k0)**2
            enddo
         enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
c>    set the forcing acceleration aF to zero after the interpolation
c>    to prevent a summation over time using gop 
      subroutine forcing_setzero()
      implicit none
c     interface
      include 'SIZE'
      include 'lib/forcing_data.f' 
c     implementation

      call rzero(aFxr,NU*NU*NU)
      call rzero(aFyr,NU*NU*NU)
      call rzero(aFzr,NU*NU*NU)

      return
      end

c-----------------------------------------------------------------------
c>    compute the forcing term in the physical domain
      subroutine forcing_computeterm(dt)
      implicit none
c     interface
      include 'SIZE'
      include 'lib/forcing_data.f' 
      real dt                    ! time step
c     implementation
      logical rflag(NU,NU,NU)    ! flag for complex pairs 
      integer i,j,k              ! loop variables
      integer ix,iy,iz           ! array indices
      integer ixr,iyr,izr        ! point reflection indices 
      real randn                 ! random even distr. number
      double precision rx,ry,rz  ! real path increments
      double precision cx,cy,cz  ! complex path increments
      double complex scp         ! scalar product wnv*b 
      double complex fft_in(NU,NU,NU)
      double complex fft_out(NU,NU,NU)
c     components of forcing acceleration aF (spectral and real)
      double complex aFx(NU,NU,NU),aFy(NU,NU,NU),aFz(NU,NU,NU)

c     initialization
      do k=1,NU
      do j=1,NU
      do i=1,NU
         rflag(i,j,k)=.false.
         aFx(i,j,k)=cmplx(0.0,0.0)
         aFy(i,j,k)=cmplx(0.0,0.0)
         aFz(i,j,k)=cmplx(0.0,0.0)
      enddo
      enddo
      enddo

c     build the acceleration term in the spectral space
      do k=-NU/2,NU/2-1,1
         do j=-NU/2,NU/2-1,1
            do i=-NU/2,NU/2-1,1
               ix=i+NU/2+1
               iy=j+NU/2+1
               iz=k+NU/2+1
c              point reflection indices
               ixr=-(i)+NU/2+1
               iyr=-(j)+NU/2+1
               izr=-(k)+NU/2+1
               if (sqrt(abv(ix,iy,iz)) .gt. 0.0 .and. 
     + sqrt(abv(ix,iy,iz)) .le. kF) then
c              compute the compl. random walk vector (6 processes)
               if (.not. rflag(ixr,iyr,izr)) then
c              build real components
               call forcing_evenrand(randn)
               rx=-1.0/TL*real(bx(ix,iy,iz))*dt +
     + sigma*sqrt(2.0*dt/TL)*randn
               call forcing_evenrand(randn)
               ry=-1.0/TL*real(by(ix,iy,iz))*dt +
     + sigma*sqrt(2.0*dt/TL)*randn
               call forcing_evenrand(randn)
               rz=-1.0/TL*real(bz(ix,iy,iz))*dt +
     + sigma*sqrt(2.0*dt/TL)*randn
c              build imaginary components
               call forcing_evenrand(randn)
               cx=-1.0/TL*imag(bx(ix,iy,iz))*dt +
     + sigma*sqrt(2.0*dt/TL)*randn
               call forcing_evenrand(randn)
               cy=-1.0/TL*imag(by(ix,iy,iz))*dt +
     + sigma*sqrt(2.0*dt/TL)*randn
               call forcing_evenrand(randn)
               cz=-1.0/TL*imag(bz(ix,iy,iz))*dt +
     + sigma*sqrt(2.0*dt/TL)*randn
               rflag(ix,iy,iz)=.true.
c              combine real and imag. parts to rand. walk vector
               bx(ix,iy,iz)=bx(ix,iy,iz)+dcmplx(rx,cx)
               by(ix,iy,iz)=by(ix,iy,iz)+dcmplx(ry,cy)
               bz(ix,iy,iz)=bz(ix,iy,iz)+dcmplx(rz,cz)
               else
               bx(ix,iy,iz)=conjg(bx(ixr,iyr,izr))
               by(ix,iy,iz)=conjg(by(ixr,iyr,izr))
               bz(ix,iy,iz)=conjg(bz(ixr,iyr,izr))
               endif
c              compute the forcing term by subtracting the  
c              projection of b onto k to fulfill the continuity eq.               
               scp=wnvx(ix,iy,iz)*bx(ix,iy,iz) +
     + wnvy(ix,iy,iz)*by(ix,iy,iz) +
     + wnvz(ix,iy,iz)*bz(ix,iy,iz)          
               aFx(ix,iy,iz)=bx(ix,iy,iz)-
     + scp*wnvx(ix,iy,iz)/abv(ix,iy,iz)
               aFy(ix,iy,iz)=by(ix,iy,iz)-
     + scp*wnvy(ix,iy,iz)/abv(ix,iy,iz)
               aFz(ix,iy,iz)=bz(ix,iy,iz)-
     + scp*wnvz(ix,iy,iz)/abv(ix,iy,iz)
               endif 
            enddo
         enddo
      enddo

c     transform the acceleration term into the phys. domain
      do i=1,3
         call forcing_cshift(aFx,NU/2,i)
         call forcing_cshift(aFy,NU/2,i)
         call forcing_cshift(aFz,NU/2,i)
      enddo
      call forcing_fft(aFx,fft_out)
      do k=1,NU
      do j=1,NU
      do i=1,NU
         aFxr(i,j,k)=real(fft_out(i,j,k))
      enddo
      enddo
      enddo
      call forcing_fft(aFy,fft_out)
      do k=1,NU
      do j=1,NU
      do i=1,NU
         aFyr(i,j,k)=real(fft_out(i,j,k))
      enddo
      enddo
      enddo
      call forcing_fft(aFz,fft_out)
      do k=1,NU
      do j=1,NU
      do i=1,NU
         aFzr(i,j,k)=real(fft_out(i,j,k))
      enddo
      enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
c>    even distributed random number generator
c>    based on "Marsaglia polar method"
      subroutine forcing_evenrand(rn)
      implicit none
c     interface
      real rn      ! random number
c     implementation
      logical ifspare
      real u,v,s   ! temp. variables
      real spare   ! second number of pair
      save spare
      save ifspare
      data ifspare  /.false./
 
      if (ifspare) then
         ifspare=.false.
         rn=spare
         return
      endif
      ifspare=.true.

      do 
         u=rand()*2.0-1.0
         v=rand()*2.0-1.0
         s=u*u+v*v
c        check if inside unit circle
         if (s .lt. 1.0 .and. s .gt. 0.0) then
            exit   ! leave the while loop
         endif
      enddo
      s=sqrt(-2.0*log(s)/s)
      spare=v*s
      rn=u*s

      return
      end

c-----------------------------------------------------------------------
c>    performs the cshift known from Fortran90
      subroutine forcing_cshift(A,shift,direction)
      implicit none
c     interface
      include 'SIZE'
      include 'lib/forcing_data.f' 
      double complex A(NU,NU,NU) 
      integer shift,direction
c     implementation 
      double complex B(NU,NU,NU)   
      integer i,j,k
      integer temp

      if (direction .eq. 1) then
         do k=1,NU
            do j=1,NU
               do i=1,NU
                  temp=mod(NU+i-shift,NU)
                  if (temp .eq. 0) temp=NU
                  B(i,j,k)=A(temp,j,k)
               enddo
            enddo
         enddo
      endif
      if (direction .eq. 2) then
         do k=1,NU
            do j=1,NU
               do i=1,NU
                  temp=mod(NU+j-shift,NU)
                  if (temp .eq. 0) temp=NU
                  B(i,j,k)=A(i,temp,k)
               enddo
            enddo
         enddo
      endif
      if (direction .eq. 3) then
         do k=1,NU
            do j=1,NU
               do i=1,NU
                  temp=mod(NU+k-shift,NU)
                  if (temp .eq. 0) temp=NU
                  B(i,j,k)=A(i,j,temp)
               enddo
            enddo
         enddo
      endif
c     return the shifted data (A=B)
      do k=1,NU
         do j=1,NU
            do i=1,NU
               A(i,j,k)=B(i,j,k)
            enddo
         enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
c>    applies backwards FFT transformation
      subroutine forcing_fft(fft_in,fft_out)
      implicit none
c     interface
      include 'SIZE'
      include 'lib/forcing_data.f'
      include 'lib/fftw3.f' 
      double complex fft_in(NU,NU,NU)
      double complex fft_out(NU,NU,NU)
c     implementation

      call dfftw_plan_dft_3D(plan,NU,NU,NU,fft_in,fft_out,
     + fftw_backward,fftw_estimate)
      call dfftw_execute (plan)
      call dfftw_destroy_plan (plan)	      

      return
      end

c-----------------------------------------------------------------------
c>    interpolate the real acceleration term from the uniform grid
c>    to the GLL grid
      subroutine forcing_interpolate(xm1,ym1,zm1,dt)
      implicit none
c     interface
      include 'SIZE'
      include 'lib/forcing_data.f'
      include 'lib/filtering_data.f'
      include 'lib/shared_data.f'
      real xm1(lx1,ly1,lz1,lelt),ym1(lx1,ly1,lz1,lelt),
     + zm1(lx1,ly1,lz1,lelt),dt
c     implementation
      integer i,j,k,l         ! loop indices
      integer iv,jv,kv        ! loop indices
      integer counter         ! sub volume index
      integer ind(3,2)        ! neighboring grid points (uniform)
      integer indx,indy,indz  ! periodic indices
      real vtot               ! sum of sub volumes
      real vol(8)             ! 8 sub volumes
      real aFtemp(NU,NU,NU)   ! temporary storage for parallel comm.

c     initialization
      vtot=udx**3
      call rzero(flwx,lx1*ly1*lz1*lelv)
      call rzero(flwy,lx1*ly1*lz1*lelv)
      call rzero(flwz,lx1*ly1*lz1*lelv)

c     parallel communication
      call rzero(aFtemp,NU*NU*NU)
      call gop(aFxr,aFtemp,'+  ',NU*NU*NU)
      call rzero(aFtemp,NU*NU*NU)
      call gop(aFyr,aFtemp,'+  ',NU*NU*NU)
      call rzero(aFtemp,NU*NU*NU)
      call gop(aFzr,aFtemp,'+  ',NU*NU*NU)

c     trilinear interpolation
      do l=1,lelv
      do k=1,lz1
      do j=1,ly1
      do i=1,lx1
         ind(1,1)=floor(xm1(i,j,k,l)/udx)+1
         ind(1,2)=ind(1,1)+1
         ind(2,1)=floor(ym1(i,j,k,l)/udx)+1
         ind(2,2)=ind(2,1)+1
         ind(3,1)=floor(zm1(i,j,k,l)/udx)+1
         ind(3,2)=ind(3,1)+1
c        compute the 8 sub volumes
         counter=1
         do kv=1,2
         do jv=1,2
         do iv=1,2
            vol(counter)=abs((xm1(i,j,k,l)-(ind(1,iv)-1)*udx)*
     + (ym1(i,j,k,l)-(ind(2,jv)-1)*udx)*
     + (zm1(i,j,k,l)-(ind(3,kv)-1)*udx))       
            counter=counter+1
         enddo
         enddo
         enddo
c        interpolate
         counter=8
         do kv=1,2
         do jv=1,2
         do iv=1,2
            indx=mod(ind(1,iv)-1,NU)+1
            indy=mod(ind(2,jv)-1,NU)+1
            indz=mod(ind(3,kv)-1,NU)+1
            flwx(i,j,k,l)=flwx(i,j,k,l)+vol(counter)/vtot*
     + aFxr(indx,indy,indz)
            flwy(i,j,k,l)=flwy(i,j,k,l)+vol(counter)/vtot*
     + aFyr(indx,indy,indz) 
            flwz(i,j,k,l)=flwz(i,j,k,l)+vol(counter)/vtot*
     + aFzr(indx,indy,indz)   
            counter=counter-1
         enddo
         enddo
         enddo
c        filter the forcing term using Euler expl. (approximate)
         if (iftles) then
            flwx(i,j,k,l)=(dt*flwx(i,j,k,l)+
     + FTS*flwx_old(i,j,k,l))/(dt+FTS)
            flwy(i,j,k,l)=(dt*flwy(i,j,k,l)+
     + FTS*flwy_old(i,j,k,l))/(dt+FTS)
            flwz(i,j,k,l)=(dt*flwz(i,j,k,l)+
     + FTS*flwz_old(i,j,k,l))/(dt+FTS)
            flwx_old(i,j,k,l)=flwx(i,j,k,l)
            flwy_old(i,j,k,l)=flwy(i,j,k,l)
            flwz_old(i,j,k,l)=flwz(i,j,k,l)
         endif
      enddo
      enddo
      enddo
      enddo

      return
      end

c-----------------------------------------------------------------------
c>    create the point source field
      subroutine forcing_ps(xm1,ym1,zm1)
      implicit none
c     interface
      include 'SIZE'
      include 'lib/forcing_data.f'
      include 'lib/shared_data.f'
      real xm1(lx1,ly1,lz1,lelt),ym1(lx1,ly1,lz1,lelt),
     + zm1(lx1,ly1,lz1,lelt)
c     implementation
      integer i,j,k,l         ! loop indices
      real tole
      parameter (tole = 1.0E-03)

      call rzero(fpsx,lx1*ly1*lz1*lelv)
      call rzero(fpsy,lx1*ly1*lz1*lelv)
      call rzero(fpsz,lx1*ly1*lz1*lelv)

c     create the forcing field
      do l=1,lelv
      do k=1,lz1
      do j=1,ly1
      do i=1,lx1
         if (xm1(i,j,k,l) .gt. (ldmx/2.0-tole) .and.
     + xm1(i,j,k,l) .lt. (ldmx/2.0+tole) .and.
     + ym1(i,j,k,l) .gt. (ldmy/2.0-tole) .and.
     + ym1(i,j,k,l) .lt. (ldmy/2.0+tole)) then
            fpsx(i,j,k,l)=1.0
            fpsy(i,j,k,l)=1.0
            fpsz(i,j,k,l)=1.0
         endif
      enddo
      enddo
      enddo
      enddo

      return
      end

c         if (xm1(i,j,k,l) .gt. (ldmx/2.0-tole) .and.
c     + xm1(i,j,k,l) .lt. (ldmx/2.0+tole) .and.
c     + ym1(i,j,k,l) .gt. (ldmy/2.0-tole) .and.
c     + ym1(i,j,k,l) .lt. (ldmy/2.0+tole) .and.
c     + zm1(i,j,k,l) .gt. (ldmz/2.0-tole) .and.
c     + zm1(i,j,k,l) .lt. (ldmz/2.0+tole)) then
c            fpsx(i,j,k,l)=100.0
c            fpsy(i,j,k,l)=100.0
c            fpsz(i,j,k,l)=100.0
c         endif


