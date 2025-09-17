C Fortran routine.
      subroutine pscmain(disgeom,ns,cxyrec,nrec,Ux,Uy,Uz,Sxx,Syy,Szz,
     &  Sxy,Syz,Szx,Tx,Ty,Gd,Gr,X,Y,grndir)
      implicit none
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     this program synthesizes seismograms due to a number of          c
c     rectanglar rupture planes using the Green's function approach.   c
c                                                                      c
c     The input data will be read from an input file                   c
c                                                                      c
c     Last modified: Potsdam, Feb, 2006, by R. Wang                    c
c                                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     BEGIN DECLARATIONS
c     ==================
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     GLOBAL CONSTANTS
c     ================
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      include 'pscglob.h'
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     MEMORIES FOR INPUTS
c     ====================
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      character*80 grndir
      character*80 green(13)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     MEMORIES FOR OUTPUTS
c     ====================
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer itout(13)
      character*80 toutfile(13)
      integer nsc,itsc(NSCMAX)
      character*80 scoutfile(NSCMAX)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     LOCAL WORK SPACES
c     =================
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer i,is,isc,irec,ixyrec,ixrec,iyrec,nxrec,nyrec
      integer nrec,ns
      double precision xrec1,xrec2,yrec1,yrec2,dxrec,dyrec
      double precision lats,lons
      double precision disgeom(14,3000)
      double precision cxyrec(2,3000)
      character*180 dataline
      logical onlysc
      character*3 nwarnc
      double precision Ux(250,3000),Uy(250,3000)
      double precision Uz(250,3000)
      double precision Sxx(250,3000),Syy(250,3000)
      double precision Szz(250,3000)
      double precision Sxy(250,3000),Syz(250,3000)
      double precision Szx(250,3000)
      double precision Tx(250,3000),Ty(250,3000)
      double precision Gd(250,3000),Gr(250,3000)
      double precision X(1,3000),Y(1,3000)

      integer iout,jout,kout
      character*6 nrecch
      double precision PI, PI2
      data PI,PI2/3.14159265358979d0,6.28318530717959d0/

c added by zhao bin
c make signature for input and ouput variables
!f2py intent(in) :: disgeom,ns,cxyrec,nrec,grndir
!f2py intent(out) :: Ux,Uy,Uz,Sxx,Syy,Szz,Sxy,Syz,Szx,Tx,Ty,Gd,Gr,X,Y
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     END DECLARATIONS
c     ================
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c00000000000000000000000000000000000000000000000000000000000000000000000
c     BEGIN INPUT PARAMETERS
c     ======================
c00000000000000000000000000000000000000000000000000000000000000000000000
c
      nwarn=0
c
c     call mexPrintf('##############################################\n')
c     call mexPrintf('#                                            #\n')
c     call mexPrintf('#                Welcome to                  #\n')
c     call mexPrintf('#                                            #\n')
c     call mexPrintf('#                                            #\n')
c     call mexPrintf('#   PPPP     SSSS    CCCC   M   M    PPPP    #\n')
c     call mexPrintf('#   P   P   S       C       MM MM    P   P   #\n')
c     call mexPrintf('#   PPPP     SSS    C       M M M    PPPP    #\n')
c     call mexPrintf('#   P           S   C       M   M    P       #\n')
c     call mexPrintf('#   P       SSSS     CCCC   M   M    P       #\n')
c     call mexPrintf('#                                            #\n')
c     call mexPrintf('#                Version 4.4                 #\n')
c     call mexPrintf('#                                            #\n')
c     call mexPrintf('#                    by                      #\n')
c     call mexPrintf('#                                            #\n')
c     call mexPrintf('#               Rongjiang Wang               #\n')
c     call mexPrintf('#            (wang@gfz-potsdam.de)           #\n')
c     call mexPrintf('#                                            #\n')
c     call mexPrintf('#         GeoForschungsZentrum Potsdam       #\n')
c     call mexPrintf('#                   July 2006                #\n')
c     call mexPrintf('##############################################\n')
c     call mexPrintf('                                              \n')
c00000000000000000000000000000000000000000000000000000000000000000000000
c     PARAMETERS FOR OBSERVATION ARRAY
c     ================================
c00000000000000000000000000000000000000000000000000000000000000000000000
      latlon=0
      ixyrec=0
        do irec=1,nrec
          xrec0(irec)=cxyrec(1,irec)
          yrec0(irec)=cxyrec(2,irec)
        enddo
c00000000000000000000000000000000000000000000000000000000000000000000000
c      OUTPUT PARAMETERS
c      =================
c00000000000000000000000000000000000000000000000000000000000000000000000
      icopost = 1
      insar = 0
      xlos = -0.072
      ylos = 0.408
      zlos = -0.910
c
      itout(1) = 0
      itout(2) = 0
      itout(3) = 0
      toutfile(1) = 'ux.dat'
      toutfile(2) = 'uy.dat'
      toutfile(3) = 'uz.dat'
      itout(4) = 0
      itout(5) = 0
      itout(6) = 0
      itout(7) = 0
      itout(8) = 0
      itout(9) = 0
      toutfile(4) = 'sxx.dat'
      toutfile(5) = 'syy.dat'
      toutfile(6) = 'szz.dat'
      toutfile(7) = 'sxy.dat'
      toutfile(8) = 'syz.dat'
      toutfile(9) = 'szx.dat'
      itout(10) = 0
      itout(11) = 0
      itout(12) = 0
      itout(13) = 0
      toutfile(10) = 'tx.dat'
      toutfile(11) = 'ty.dat'
      toutfile(12) = 'gd.dat'
      toutfile(13) = 'gr.dat'
      nsc = 2
        itsc(1) = 1
        scoutfile(1) = 'all001.dat'
        itsc(2) = 2
        scoutfile(2) = 'all002.dat'
      onlysc=.true.
c00000000000000000000000000000000000000000000000000000000000000000000000
c     PARAMETERS FOR EARTH MODEL CHOICE
c     =================================
c00000000000000000000000000000000000000000000000000000000000000000000000
      green(1) = 'uz'
      green(2) = 'ur'
      green(3) = 'ut'
      green(4) = 'szz'
      green(5) = 'srr'
      green(6) = 'stt'
      green(7) = 'szr'
      green(8) = 'srt'
      green(9) = 'stz'
      green(10) = 'tr'
      green(11) = 'tt'
      green(12) = 'gd'
      green(13) = 'gr'
c00000000000000000000000000000000000000000000000000000000000000000000000
c     PARAMETERS FOR RECTANGULAR SOURCES
c     ==================================
c00000000000000000000000000000000000000000000000000000000000000000000000
      lat0 = 0.0
      lon0 = 0.0
      do is=1,ns
      slip(is) = sqrt(disgeom(9,is)**2 + disgeom(8,is)**2)
      openning(is) = disgeom(10,is)
c convert disgeom km to lat/lon
      lats = (KM2M*disgeom(7,is)/REARTH)/DEG2RAD
      lons = (KM2M*disgeom(6,is)/REARTH)/DEG2RAD
      zref(is) = disgeom(3,is)
      length(is) = disgeom(1,is)
      width(is) = disgeom(2,is)
      strike(is) = disgeom(5,is)
      dip(is) = disgeom(4,is)
      rake(is) = atan2(disgeom(9,is),disgeom(8,is))/DEG2RAD
      taplt(is) = disgeom(11,is)
      taprt(is) = disgeom(12,is)
      tapup(is) = disgeom(13,is)
      taplw(is) = disgeom(14,is)
        zref(is)=KM2M*zref(is)
        length(is)=KM2M*length(is)
        width(is)=KM2M*width(is)
        call disazi(REARTH,lat0,lon0,lats,lons,xref(is),yref(is))
c
c       Note unit of xref(is) and yref(is) returned from disazi is meter!
c
        if(zref(is).lt.0.d0)then
c         call mexErrMsgTxt(' Error: source depth zs0 < 0!\n')
        endif
        if(dabs(strike(is)).gt.360.d0)then
c         call mexErrMsgTxt(' Error: wrong strike angle!\n')
        endif
        if(strike(is).lt.0.d0)then
          strike(is)=strike(is)+360.d0
        endif
        if(dip(is).gt.90.d0.or.dip(is).lt.0.d0)then
c         call mexErrMsgTxt(' Error: wrong dip angle!\n')
        endif
        if(dabs(rake(is)).gt.360.d0)then
c         call mexErrMsgTxt(' Error: wrong rake angle!\n')
        else if(rake(is).gt.180.d0)then
          rake(is)=rake(is)-360.d0
        else if(rake(is).le.-180.d0)then
          rake(is)=rake(is)+360.d0
        endif
      enddo
c
c00000000000000000000000000000000000000000000000000000000000000000000000
c      END INPUT PARAMETERS
c      ====================
c00000000000000000000000000000000000000000000000000000000000000000000000
c     call mexPrintf('... input data successful ...\n')
c00000000000000000000000000000000000000000000000000000000000000000000000
c      BEGIN PROCESSING
c      ================
c00000000000000000000000000000000000000000000000000000000000000000000000
      nwarn=0
      if(latlon.eq.1)then
c       call mexPrintf('... transforming geographic to cartesian ...\n')
        do irec=1,nrec
          call disazi(REARTH,lat0,lon0,xrec0(irec),yrec0(irec),
     &                  xrec(irec),yrec(irec))
c       call mexPrintf('transforming geographic to cartesian\n')
        enddo
      else
        do irec=1,nrec
          xrec(irec)=xrec0(irec)*KM2M
          yrec(irec)=yrec0(irec)*KM2M
        enddo
      endif
c     call mexPrintf('... using the Green function approach ...\n')
      call pscgrn(ns,nrec,grndir,green,onlysc,nsc,itsc)
c00000000000000000000000000000000000000000000000000000000000000000000000
c      END OF STANDARD PROCESSING
c      ==========================
c00000000000000000000000000000000000000000000000000000000000000000000000
      if(nwarn.eq.0)then
c       call mexPrintf('############################################\n')
c       call mexPrintf('#                                          #\n')
c       call mexPrintf('#      End of computations with PSCMP      #\n')
c       call mexPrintf('#                                          #\n')
c       call mexPrintf('############################################\n')
      else
        write(nwarnc,888)nwarn
  888   format(i3)
c       call mexPrintf('###########################################\n')
c       call mexPrintf(' Sorry, there have been'//nwarnc
c    &     //' warnings. \n')
c       call mexPrintf('         Results may be inaccurate!        \n')
c       call mexPrintf('###########################################\n')
      endif
c
      do iout=1,250
        do jout=1,3000
          Ux(iout,jout)=obs(iout,jout,1)
          Uy(iout,jout)=obs(iout,jout,2)
          Uz(iout,jout)=obs(iout,jout,3)
          Sxx(iout,jout)=obs(iout,jout,4)
          Syy(iout,jout)=obs(iout,jout,5)
          Szz(iout,jout)=obs(iout,jout,6)
          Sxy(iout,jout)=obs(iout,jout,7)
          Syz(iout,jout)=obs(iout,jout,8)
          Szx(iout,jout)=obs(iout,jout,9)
          Tx(iout,jout)=obs(iout,jout,10)
          Ty(iout,jout)=obs(iout,jout,11)
          Gd(iout,jout)=obs(iout,jout,12)
          Gr(iout,jout)=obs(iout,jout,13)
          if(iout.eq.1)then
            X(iout,jout)=xrec(jout)
            Y(iout,jout)=yrec(jout)
          endif
        enddo
      enddo
      write(nrecch,889)nrec
 889  format(i6)
c     call mexprintf(' 15 arrays of'//nrecch//' values returned.\n')
      print *, Ux(1,1), Uy(1,1), Uz(1,1)
      end
