c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     GLOBAL CONSTANTS
c     ================
c
c     NRECMAX = max. number of observation positions
c     NZSMAX = max. number of the discrete source depths
c     NRMAX = max. number of the discrete radial diatances
c     NSMAX = max. number of the source rectangles
c     NPSMAX = max. number of discrete point sources per source depth
c     NTMAX = max. number of time samples used for Green's functions
c     NTRMAX = max. number of time samples of the outputs
c     NSCMAX = max. number of scenario outputs
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer NZSMAX,NRMAX,NSMAX,NPSMAX,NRECMAX,NTMAX,NTRMAX,NSCMAX
C code modified version 4.4
c      parameter(NZSMAX=100,NRMAX=151)
c      parameter(NZSMAX=100,NRMAX=501)
      parameter(NZSMAX=150,NRMAX=1501)
      parameter(NSMAX=1000,NPSMAX=50000)
      parameter(NRECMAX=3000)
C code modified version 4.4
c     parameter(NRECMAX=22801)
c     parameter(NRECMAX=15000)
      parameter(NTMAX=256,NTRMAX=NTMAX)
      parameter(NSCMAX=100)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     RECTANGULAR SOURCE PLANES
c     =========================
c
c     (xs0,ys0,zs0) = coordinates of the reference point
c     with x = north, y = east, z = downward.
c     all angles in degree.
c     NSMAX = the max. number of source rectangles
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      double precision slip(NSMAX),openning(NSMAX)
      double precision xref(NSMAX),yref(NSMAX),zref(NSMAX)
      double precision length(NSMAX),width(NSMAX)
      double precision strike(NSMAX),dip(NSMAX),rake(NSMAX)
      double precision magfac(NSMAX)
      double precision taplt(NSMAX),taprt(NSMAX)
      double precision tapup(NSMAX),taplw(NSMAX)
c
      common/rects/slip,openning,xref,yref,zref,length,width,
     &             strike,dip,rake,magfac,taplt,taprt,tapup,taplw
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     DISTRETE POINT SOURCES
c     ======================
c
c     (xs,ys,zs) = coordinates of the discrete point sources
c     with x = north, y = east, z = downward
c     angles in degree.
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer nps(NZSMAX)
      double precision px(NPSMAX,NZSMAX),py(NPSMAX,NZSMAX)
      double precision pz(NPSMAX,NZSMAX),pmwei(6,NPSMAX,NZSMAX)
c
      common/points/px,py,pz,pmwei
      common/ipoints/nps
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     GREEN'S FUNCTION INFO
c     =====================
c
c     nzs,zs1,zs2 = number of depth samples, start and end depths used
c           in Green's functions
c     nr,r1,r2 = number of distance samples, start and end distances used
c           in Green's functions
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer nr,nzs,nt
      double precision r1,r2,zs1,zs2
      double precision zrec,larec,murec,rhorec,etarec,relaxrec,
     +                 twindow
      common /greeninfo/r1,r2,zs1,zs2,zrec,larec,murec,rhorec,
     +                  etarec,relaxrec,twindow
      common /igreeninfo/nr,nzs,nt
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     OBSERVATION POSITIONS AND OBSERVABLES
c     =====================================
c
c     (xrec(i),yrec(i))=coordinates of the observation positions
c     the 3 displcement/velocity/acceleration components: ux,uy,uz
c     NRECMAX = the max. number of observation positions
c     latlon: 1/0 = geographic/Cartesian
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer ntrec,latlon
      double precision xrec0(NRECMAX),yrec0(NRECMAX)
      double precision xrec(NRECMAX),yrec(NRECMAX)
      double precision obs(NTRMAX,NRECMAX,13)
c
      common/obsarray/xrec0,yrec0,xrec,yrec,obs
      common/iobsarray/ntrec,latlon
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     EARTHQUAKE'S GLOBAL PARAMETERS
c     ==============================
c
c     lat0 = epicentre's latitude
c     lon0 = epicentre's east longitude
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      double precision lat0,lon0
c
      common/coororigin/lat0,lon0
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     OUTPUT SWITCH
c     =============
c
c     icopost = 1: co and postseismic changes
c               0: only post-seismic changes
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer icopost
c
      common/cpseismic/icopost
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     COSINES OF LOS TO INSAR ORBIT
c     =============================
c
c     insar = 1: output los displacements
c             0: not output los displacements
c     xlos, ylos, zlos = cosines of the los
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer insar
      double precision xlos,ylos,zlos
c
      common/coororigin/xlos,ylos,zlos,insar
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     WARNING STATISTICS
c     ==================
c
c     nwarn = total number of warnings
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer nwarn
c
      common/warnings/nwarn
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     LOCAL CONSTANTS
c     ==============
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      double precision DEG2RAD,KM2M,DAY2SEC,REARTH,G0
      parameter(DEG2RAD=1.745329252d-02,KM2M=1.0d+03)
      parameter(DAY2SEC=8.64d+04,REARTH=6.371d+06,G0=9.82d+00)
