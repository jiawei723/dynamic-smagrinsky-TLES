c-----------------------------------------------------------------------
c>    initialize the statistics
      subroutine statistics_init()
      implicit none
c     interface
      include 'SIZE'
      include 'lib/shared_data.f'
      include 'lib/statistics_data.f'
c     implementation
      integer n,i,ios
      character*200 fname

      n = lx1*ly1*lz1*nelv
      call rzero(ux_avgz,n)
      call rzero(uy_avgz,n)
      call rzero(uz_avgz,n)
      call gtpp_gs_setup(igs_z,resx,resy,resz,3)

c     read in the sample points 'phill.his'
      write(fname,71)
      open(unit=1,file=fname,status='old',iostat=ios)
      read(1,72) nprobes
      do i=1,nprobes
         read(1,73) probe_cords(1,i),probe_cords(2,i),probe_cords(3,i)
      enddo
      close(1)

      write(*,*) 'nprobes = ',nprobes

   71 format('phill.his')
   72 format(I5,A)
   73 format(3E17.9)

      return
      end

c-----------------------------------------------------------------------
c>    compute the average in the z-direction
      subroutine statistics_computeaverage_z()
      implicit none
c     interface
      include 'SIZE'
      include 'TOTAL'
      include 'lib/shared_data.f'
      include 'lib/statistics_data.f'
      include 'lib/filtering_data.f'
c     implementation

c     contract in z: u_avgz=f(x,y)
      call planar_avg(ux_avgz,vx,igs_z)
      call planar_avg(uy_avgz,vy,igs_z)
      call planar_avg(uz_avgz,vz,igs_z)

c     TLES
      if (iftles) then
         call planar_avg(uxdconv_avgz,uxdconv,igs_z)
         call planar_avg(uydconv_avgz,uydconv,igs_z)
         call planar_avg(uzdconv_avgz,uzdconv,igs_z)
      endif

      return
      end

c----------------------------------------------------------------------
c>    interpolates the 'u_avgz' field to probe points 
      subroutine statistics_interpolate_samples(probe_length)
      implicit none
c     interface
      include 'SIZE'
      include 'lib/shared_data.f'
      include 'lib/statistics_data.f'
      integer probe_length
c     implementation    
      integer i, ih, iwk(maxsamples,3),ntot
      real tolin,rwk(maxsamples,3+1),pts(3*maxsamples)
      real fwrk(lx1*ly1*lz1*lelv,3),fpts(3*maxsamples)
      
      common /rwk_stat_intp/ fwrk, rwk, fpts, pts
      common /iwk_stat_intp/ iwk
      
      integer icalld,e
      save    icalld
      data    icalld /0/
      save    ih
      data    ih /0/      

      ntot  = nx1*ny1*nz1*nelt ! total number of elements

c     initialize the interpolation routine
      if (icalld .eq. 0) then
        icalld = 1
        tolin  = 1.e-8
        call interp_setup(ih,tolin,lx1,nelt)
      endif

c     check if probe_length < maxsamples
      if (probe_length .gt. maxsamples) then
         call exitti ('probe_length > maxsamples!$',probe_length)
      endif

c     assemble the coordinates of the sample points
      do i=1,probe_length
         pts(i)                  = probe_cords(1,i)
         pts(i + probe_length)   = probe_cords(2,i)
         pts(i + probe_length*2) = probe_cords(3,i)
      enddo
  
c     pack the working array 'fwrk'
      call opcopy(fwrk(1,1),fwrk(1,2),fwrk(1,3),ux_avgz,uy_avgz,uz_avgz)

c     interpolate u_avgz to the sample points array 'pts'
      call interp_nfld(fpts,fwrk,ndim,pts(1),pts(1+probe_length),
     + pts(2*probe_length+1),probe_length,iwk,rwk,maxsamples,.true.,ih)

c     save the values as u_profile
      do i=1,probe_length
         u_profile(1,i) = fpts(i)
         u_profile(2,i) = fpts(i + probe_length)
         u_profile(3,i) = fpts(i + probe_length*2)
      enddo

c     TLES
      if (iftles) then
         call opcopy(fwrk(1,1),fwrk(1,2),fwrk(1,3),uxdconv_avgz,
     + uydconv_avgz,uzdconv_avgz)
         call interp_nfld(fpts,fwrk,ndim,pts(1),pts(1+probe_length),
     + pts(2*probe_length+1),probe_length,iwk,rwk,maxsamples,.true.,ih)
         do i=1,probe_length
            udconv_profile(1,i) = fpts(i)
            udconv_profile(2,i) = fpts(i + probe_length)
            udconv_profile(3,i) = fpts(i + probe_length*2)
         enddo
      endif
      
      return
      end

c-----------------------------------------------------------------------
c>    write out the velocity slice averaged in the z-direction
      subroutine statistics_write_slice(xm1,ym1,zm1,time)
      implicit none
c     interface
      include 'SIZE'
      include 'lib/shared_data.f'
      include 'lib/statistics_data.f'
      real xm1(lx1,ly1,lz1,lelt),ym1(lx1,ly1,lz1,lelt)
      real zm1(lx1,ly1,lz1,lelt)
      real time
c     implementation
      integer ios,m,k,j 
      character*200 fname
      character*200 fname_dconv
      character*5 fnid
      character*5 ficlld

      integer iclld ! timestep
      save    iclld
      data    iclld  /0/

      write(fnid,74) nid
      write(ficlld,74) iclld
      write(fname,71)
      fname = trim(fname)//trim(fnid)//"_t_"//trim(ficlld)//".txt"

c     TLES
      write(fname_dconv,75)
      fname_dconv = trim(fname_dconv)//trim(fnid)//
     + "_t_"//trim(ficlld)//".txt"

c     slice
      open(unit=1,file=fname,status='unknown',iostat=ios)
      if (ios .eq. 0) close(1,status='delete')
      open(unit=1,file=fname,status="new",action="write")

c     TLES
      if (iftles) then
         open(unit=2,file=fname_dconv,status='unknown',iostat=ios)
         if (ios .eq. 0) close(2,status='delete')
         open(unit=2,file=fname_dconv,status="new",action="write")
         write(2,72) iclld, time, lx1, ly1, nelv
      endif

      write(1,72) iclld, time, lx1, ly1, nelv
      do m = 1,nelv
      do k = 1,ly1
      do j = 1,lx1
         if (zm1(j,k,1,m) .eq. 0.0) then
            write(1,73) xm1(j,k,1,m),ym1(j,k,1,m),ux_avgz(j,k,1,m)

c           TLES
            if (iftles) then
               write(2,73) xm1(j,k,1,m),ym1(j,k,1,m),
     + uxdconv_avgz(j,k,1,m)
            endif
         endif
      enddo
      enddo
      enddo
      close(1)

      iclld = iclld + 1

   71 format('out/stat/slice_n_')
   72 format(I10,E17.9,3I10)
   73 format(3E17.9)
   74 format(I0.3)
   75 format('out/stat_dconv/slice_n_')

      return
      end

c-----------------------------------------------------------------------
c>    write out the velocity profiles
      subroutine statistics_write_profiles(time)
      implicit none
c     interface
      include 'SIZE'
      include 'lib/shared_data.f'
      include 'lib/statistics_data.f'
      real time
c     implementation
      integer i,ios 
      character*200 fname
      character*200 fname_dconv
      character*5 ficlld

      integer iclld ! timestep
      save    iclld
      data    iclld  /0/

      write(fname,71)
      write(ficlld,74) iclld
      fname = trim(fname)//trim(ficlld)//".txt"
      write(fname_dconv,75)
      fname_dconv = trim(fname_dconv)//trim(ficlld)//".txt"

c     slice
      open(unit=1,file=fname,status='unknown',iostat=ios)
      if (ios .eq. 0) close(1,status='delete')
      open(unit=1,file=fname,status="new",action="write")

c     TLES
      if (iftles) then
         open(unit=2,file=fname_dconv,status='unknown',iostat=ios)
         if (ios .eq. 0) close(2,status='delete')
         open(unit=2,file=fname_dconv,status="new",action="write")
         write(2,72) iclld, time
      endif

      write(1,72) iclld, time
      do i = 1,nprobes
         write(1,73) u_profile(1,i),u_profile(2,i),u_profile(3,i)
   
c        TLES
         if (iftles) then
            write(2,73) udconv_profile(1,i),udconv_profile(2,i),
     + udconv_profile(3,i)
         endif

      enddo
      close(1)

      iclld = iclld + 1

   71 format('out/profiles/t_')
   72 format(I10,E17.9)
   73 format(3E17.9)
   74 format(I0.3)
   75 format('out/profiles_dconv/t_')

      return
      end

