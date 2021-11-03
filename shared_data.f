c     domain size
      real ldmx,ldmy,ldmz

c     number of elements - needs to match value for nelx and co in .box file
      integer resx
      parameter (resx = 28)
      integer resy
      parameter (resy = 18)
      integer resz
      parameter (resz = 26)

c     fluid properties
      real rhogi ! initial density
      parameter (rhogi = 1.0)
      real mugi ! initial dynamic viscosity
      parameter (mugi = 4.491E-03)

c     choose simulation case
      logical ifchannel ! channel flow simulation
      parameter (ifchannel = .false.)
      logical iftgv ! Taylor-Green vortex simulation
      parameter (iftgv = .false.)

c     forcing schemes and parameters
      logical ifcvf  ! constant volume flow
      parameter (ifcvf = .false.)
      real utarget ! target velocity
      parameter (utarget = 15.71808)
      logical iflfs ! linear forcing scheme
      parameter (iflfs = .false.)
      logical ifuf ! update forcing parameter
      parameter (ifuf = .false.)
      logical iflwf ! low wave number forcing
      parameter (iflwf = .false.)
      logical ifps ! point source field to initiate transition
      parameter (ifps = .false.)
      integer ips ! time step at which the point source is applied
      parameter (ips = 500)
      real Afr ! linear forcing parameter
      parameter (Afr = 0.20)
      real epsr ! reference dissipation (DNS sim.)
      parameter (epsr = 0.254384577)

c     misc parameters and data
      logical ifksm ! initialize fluid velocity randomly
      parameter (ifksm = .false.)
      logical ifvalid ! overwrite velocity for validation study
      parameter (ifvalid = .false.)
      logical iftimeseries ! write out u(t) of a grid point
      parameter (iftimeseries = .false.)
      integer statisticstep ! iostep for turb. statistics
      parameter (statisticstep = 10)

c     restart parameters
      logical ifrestart ! approx. the vel. gradients with 2 prev. sim.
      parameter (ifrestart = .false.)
      logical ifretles ! input velocity fields are already filtered
      parameter (ifretles = .false.)
      logical iffullrestart ! read in everything (only TLES)
      parameter (iffullrestart = .false.)
      real redt ! time step size of restart files
      parameter (redt = 5.0E-04)

c     TLES parameters
      logical iftles ! use TLES
      parameter (iftles = .true.)
      real Tf ! filter time scale factor FTS=Tf*dt
      parameter (Tf = 20.0)
      logical if2ndreg ! secondary regularization
      parameter (if2ndreg = .false.)
      logical ifadreg ! approx. deconv. regularization
      parameter (ifadreg = .true.)
      integer regorder! order q of the reg. term (3 or 4)
      parameter (regorder = 4)
      real chi ! regularization coefficient
      parameter (chi = 2.0)
      real Tf2 ! filter width for regularization
      parameter (Tf2 = 100.0)
      real diff ! diffusion coefficient for tij
      parameter (diff=0.0)
      logical iframp ! increase the filter width slowly (num. stab.)
      parameter (iframp = .false.)
      logical ifdudtbd2 ! compute du/dt using 2nd order backw difference
      parameter (ifdudtbd2 = .false.)
      logical ifdudtbd3 ! compute du/dt using 3rd order backw difference
      parameter (ifdudtbd3 = .true.)
      logical ifrkf ! integrate the stress tensor using 4th order RK
      parameter (ifrkf = .false.)
      logical ifee ! integrate the stress tensor using Euler explicit
      parameter (ifee = .true.)
      logical ifei ! integrate the stress tensor using exact integration
      parameter (ifei = .false.)
      logical ifdiag ! write out tij, ui and more for diagnostics
      parameter (ifdiag = .false.)

c     options for statistics and velocity profile generation
      logical ifstat ! generate statistics
      parameter(ifstat = .true.)
      integer numProfiles ! number of profiles generated at equal distances
      parameter(numProfiles = 9)
      real profileDist ! distance between two profiles, first profile locatet at 0
      parameter(profileDist = 1)
      integer profileStart ! first time step at which profile is generated
      parameter(profileStart = 10000)
      integer profileStop ! last time step at which profile is generated
      parameter(profileStop = 1000000)
      integer profileStep ! time step interval between two sets of profiles
      parameter(profileStep = 500)

c     common blocks
      common /shared_real/ ldmx,ldmy,ldmz
