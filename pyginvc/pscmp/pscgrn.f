      subroutine pscgrn(ns,nrec,grndir,green,onlysc,nsc,itsc)
      implicit none
c
c     Last modified: Potsdam, Feb, 2002, by R. Wang
c
      integer ns,nrec,nsc,itsc(nsc)
      character*80 grndir,green(13)
      logical onlysc
c
      include 'pscglob.h'
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     GREEN'S FUNNCTIONN PARAMETERS
c     =============================
c
c     Green's function source types:
c       1 = explosion (m11=m22=m33=kappa)
c       2 = strike-slip (m12=m21=mue)
c       3 = dip-slip (m13=m31=mue)
c       4 = compensated linear vector dipole (m33=mue, m11=m22=-m33/2)
c     Green's function coordinate system:
c       (z,r,t) = cylindrical with z being downward(!)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer unit(13,4)
      integer idec(NRMAX,13,4)
      integer igrns(NTMAX,NRMAX,13,4)
      character*163 greens(13,4)
      logical select(13,4)
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     LOCAL WORK SPACES
c     =================
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer i,n,irec,ips,isc,ir,izs,nzs1,nzs2,ntr,nsmall,nlarge
      integer it,itr,lend,lenf,id1,id2,istp,npsum
      double precision si,co,si2,co2,dis
      double precision azi,ur,ut,uz,gr
      double precision ezz,err,ett,ezr,ert,etz
      double precision tr,tt,gd
      double precision dr,dt
      double precision zs,las,mus,rhos,etas,relaxs
      double precision psss,shss,psds,shds,pscl,psep
      double precision d1,d2,d3,d4,d5,d6
      character*180 dataline
      character npsc*4,zsc*12,nsmallc*5,nlargec*5
c     OPEN GREEN'S FUNCTIONS FILES
c     ============================
c
      do lend=80,1,-1
        if(grndir(lend:lend).ne.' ')goto 100
      enddo
100   continue
      do i=1,13
        do lenf=80,1,-1
          if(green(i)(lenf:lenf).ne.' ')goto 110
        enddo
110     continue
        greens(i,1)=grndir(1:lend)//green(i)(1:lenf)//'.ep'
        greens(i,2)=grndir(1:lend)//green(i)(1:lenf)//'.ss'
        greens(i,3)=grndir(1:lend)//green(i)(1:lenf)//'.ds'
        greens(i,4)=grndir(1:lend)//green(i)(1:lenf)//'.cl'
      enddo
c
      do istp=1,4
        do i=1,13
          select(i,istp)=.true.
        enddo
      enddo
      select(3,1)=.false.
      select(8,1)=.false.
      select(9,1)=.false.
      select(11,1)=.false.
      select(3,4)=.false.
      select(8,4)=.false.
      select(9,4)=.false.
      select(11,4)=.false.
c
      do istp=1,4
        do i=1,13
          if(.not.select(i,istp))goto 300
          unit(i,istp)=10+13*(istp-1)+i
          open(unit(i,istp),file=greens(i,istp),status='old')
          if(i*istp.eq.1)then
            call getdata(unit(i,istp),dataline)
            read(dataline,*)nr,r1,r2
            if(nr.gt.NRMAX)then
c             call mexErrMsgTxt('Error: NRMAX too small defined!')
            endif
            call getdata(unit(i,istp),dataline)
            read(dataline,*)zrec,larec,murec,rhorec,etarec,
     &                      relaxrec
            call getdata(unit(i,istp),dataline)
            read(dataline,*)nzs,zs1,zs2
            if(nzs.gt.NZSMAX)then
c             call mexErrMsgTxt('Error: NZSMAX too small defined!')
            endif
            call getdata(unit(i,istp),dataline)
            read(dataline,*)nt,twindow
            if(nt.gt.NTMAX)then
c             call mexErrMsgTxt('Error: NTMAX too small defined!')
            endif
          else
            call getdata(unit(i,istp),dataline)
            read(dataline,*)n,d1,d2
            if(n.ne.nr.or.d1.ne.r1.or.d2.ne.r2)then
c             call mexErrMsgTxt('Error: different observation
c    &           sampling in Greens!')
            endif
            call getdata(unit(i,istp),dataline)
            read(dataline,*)d1,d2,d3,d4,d5,d6
            if(d1.ne.zrec.or.d2.ne.larec.or.
     &         d3.ne.murec.or.d4.ne.rhorec.or.
     &         d5.ne.etarec.or.d6.ne.relaxrec)then
c             call mexErrMsgTxt('Error: diff. observation site
c    &           parameters in Greens!')
            endif
            call getdata(unit(i,istp),dataline)
            read(dataline,*)n,d1,d2
            if(n.ne.nzs.or.d1.ne.zs1.or.d2.ne.zs2)then
c             call mexErrMsgTxt('Error: different source sampling in
c    &           Greens!')
            endif
            call getdata(unit(i,istp),dataline)
            read(dataline,*)n,d1
            if(n.ne.nt.or.d1.ne.twindow)then
c             call mexErrMsgTxt('Error: different time sampling in
c    &           Greens!')
            endif
          endif
300       continue
        enddo
      enddo
c
c     all Green's function files have been opened
c     ===========================================
c
      if(nt.gt.1)then
        dt=twindow/dble(nt-1)
      else
        dt=twindow
      endif
      dr=(r2-r1)/dble(nr-1)
c
c     INITIALIZATION
c     ==============
c
      if(onlysc)then
        ntr=nsc
      else
        ntr=nt
      endif
      if(ntr.gt.NTRMAX)then
c       call mexErrMsgTxt(' Error: NTRMAX defined too small!')
      endif
      do isc=1,nsc
        if(itsc(isc).lt.1)then
          nwarn=nwarn+1
c         call mexPrintf(' Warning in pscgrn: time point')
c         call mexPrintf(' outside time window,\n')
          itsc(isc)=1
        else if(itsc(isc).gt.nt)then
          nwarn=nwarn+1
c         call mexPrintf(' Warning in pscgrn: time point')
c         call mexPrintf(' outside time window,\n')
          itsc(isc)=nt
        endif
      enddo
c
c     DISCRETISATION OF RECTANGULAR PLANE SOURCES
c     ===========================================
c
      call pscdisc(ns,npsum)
c
c     SUPERPOSITION OF ALL DISCRETE POINT SOURCES
c     ===========================================
c
c     call mexPrintf('... superposition of all discrete point')
c     call mexPrintf(' sources ...\n')
      nzs1=1
      nzs2=0
      do izs=1,nzs
        if(nps(izs).gt.0)then
          nzs2=izs
        endif
      enddo
c
      if(icopost.eq.0)then
        do i=1,13
          do irec=1,nrec
            do itr=1,ntr
              obs(itr,irec,i)=0.d0
            enddo
          enddo
        enddo
      else
        call pscokada(ns,nrec,nzs1,nzs2)
        do i=1,13
          do irec=1,nrec
            do itr=2,ntr
              obs(itr,irec,i)=obs(1,irec,i)
            enddo
          enddo
        enddo
      endif
c
      do izs=nzs1,nzs2
        nsmall=0
        nlarge=0
c
c       read in Green's functions
c
        do istp=1,4
          do i=1,13
            if(.not.select(i,istp))goto 400
            if(i.eq.1)then
              call getdata(unit(i,istp),dataline)
              read(dataline,*)zs,las,mus,rhos,etas,relaxs
            else
              call getdata(unit(i,istp),dataline)
              read(dataline,*)d1,d2,d3,d4,d5,d6
              if(d1.ne.zs.or.d2.ne.las.or.
     &           d3.ne.mus.or.d4.ne.rhos.or.
     &           d5.ne.etas.or.d6.ne.relaxs)then
c               call mexErrMsgTxt('Error: different s. layer
c    &             parameters in greens!')
              endif
            endif
            read(unit(i,istp),*)(idec(ir,i,istp),ir=1,nr)
            read(unit(i,istp),*)((igrns(it,ir,i,istp),
     &                           ir=1,nr),it=1,nt)
            if(icopost.eq.0)then
              do ir=1,nr
                do it=2,nt
                  igrns(it,ir,i,istp)=igrns(it,ir,i,istp)
     &                               -igrns(1,ir,i,istp)
                enddo
                igrns(1,ir,i,istp)=0
              enddo
            endif
400         continue
          enddo
        enddo
        if(nps(izs).ge.1)then
c
          write(npsc,886)nps(izs)
 886      format(i4)
          write(zsc,887)zs/KM2M
 887      format(f12.4)
c         call mexPrintf(' processing '//npsc//
c    &       ' point sources at depth '//zsc//' km\n')
        endif
c
        do ips=1,nps(izs)
c
          do irec=1,nrec
            dis=dsqrt((xrec(irec)-px(ips,izs))**2
     &               +(yrec(irec)-py(ips,izs))**2)
            id1=idint((dis-r1)/dr)+1
            if(id1.gt.nr)then
C             call mexPrintf(' Warning: too large distances ignored!\n')
              nlarge=nlarge+1
            else
              if(dis.gt.0.d0)then
                azi=datan2(yrec(irec)-py(ips,izs),
     &                     xrec(irec)-px(ips,izs))
              else
                azi=0.d0
              endif
              if(id1.le.0)then
                id1=1
                id2=1
                d1=1.d0
                d2=0.d0
              else if(id1.lt.nr)then
                d2=dmod(dis/dr,1.d0)
                d1=1.d0-d2
                id2=id1+1
              else
                d1=1.d0
                d2=0.d0
                id2=id1
              endif
c
              co=dcos(azi)
              si=dsin(azi)
              co2=dcos(2.d0*azi)
              si2=dsin(2.d0*azi)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c             pmwei(1-6):
c             1 = weight for strike-slip: m12=m21=1;
c             poloidal*sin(2 * theta), toroidal*cos(2 * theta)
c
c             2 = weight for dip-slip: m13=m31=1
c             poloidal * cos(theta), toroidal * sin(theta)
c
c             3 = weight for clvd: m33=-m11=-m22=1
c             axisymmetric
c
c             4 = weight for 45 deg strike-slip: m11=-m22=1
c             greenfct4(theta) = green1(theta + 45 deg)
c
c             5 = weight for 45 deg dip-slip: m23=m32=1
c             greenfct5(theta) = green2(theta - 90 deg)
c
c             6 = weight for explosion
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
              psep=pmwei(6,ips,izs)
              psss=pmwei(1,ips,izs)*si2+pmwei(4,ips,izs)*co2
              shss=pmwei(1,ips,izs)*co2-pmwei(4,ips,izs)*si2
              psds=pmwei(2,ips,izs)*co+pmwei(5,ips,izs)*si
              shds=pmwei(2,ips,izs)*si-pmwei(5,ips,izs)*co
              pscl=pmwei(3,ips,izs)
c
              do itr=1,ntr
                if(onlysc)then
                  it=itsc(itr)
                else
                  it=itr
                endif
                uz=0.d0
                ur=0.d0
                ut=0.d0
                ezz=0.d0
                err=0.d0
                ett=0.d0
                ezr=0.d0
                ert=0.d0
                etz=0.d0
                tr=0.d0
                tt=0.d0
                gd=0.d0
                gr=0.d0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c               contributions from the explosion components
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                uz=uz+psep*(d1*dble(igrns(it,id1,1,1))
     +                        *10.d0**idec(id1,1,1)
     +                     +d2*dble(igrns(it,id2,1,1))
     +                        *10.d0**idec(id2,1,1))
                ur=ur+psep*(d1*dble(igrns(it,id1,2,1))
     +                        *10.d0**idec(id1,2,1)
     +                     +d2*dble(igrns(it,id2,2,1))
     +                        *10.d0**idec(id2,2,1))
c
                ezz=ezz+psep*(d1*dble(igrns(it,id1,4,1))
     +                          *10.d0**idec(id1,4,1)
     +                       +d2*dble(igrns(it,id2,4,1))
     +                          *10.d0**idec(id2,4,1))
                err=err+psep*(d1*dble(igrns(it,id1,5,1))
     +                          *10.d0**idec(id1,5,1)
     +                       +d2*dble(igrns(it,id2,5,1))
     +                          *10.d0**idec(id2,5,1))
                ett=ett+psep*(d1*dble(igrns(it,id1,6,1))
     +                          *10.d0**idec(id1,6,1)
     +                       +d2*dble(igrns(it,id2,6,1))
     +                          *10.d0**idec(id2,6,1))
                ezr=ezr+psep*(d1*dble(igrns(it,id1,7,1))
     +                          *10.d0**idec(id1,7,1)
     +                       +d2*dble(igrns(it,id2,7,1))
     +                          *10.d0**idec(id2,7,1))
c
                tr=tr+psep*(d1*dble(igrns(it,id1,10,1))
     +                        *10.d0**idec(id1,10,1)
     +                     +d2*dble(igrns(it,id2,10,1))
     +                        *10.d0**idec(id2,10,1))
                gd=gd+psep*(d1*dble(igrns(it,id1,12,1))
     +                        *10.d0**idec(id1,12,1)
     +                     +d2*dble(igrns(it,id2,12,1))
     +                        *10.d0**idec(id2,12,1))
c
                gr=gr+psep*(d1*dble(igrns(it,id1,13,1))
     +                        *10.d0**idec(id1,13,1)
     +                     +d2*dble(igrns(it,id2,13,1))
     +                        *10.d0**idec(id2,13,1))
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c               contributions from the strike-slip components
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                uz=uz+psss*(d1*dble(igrns(it,id1,1,2))
     +                     *10.d0**idec(id1,1,2)
     +                  +d2*dble(igrns(it,id2,1,2))
     +                     *10.d0**idec(id2,1,2))
                ur=ur+psss*(d1*dble(igrns(it,id1,2,2))
     +                     *10.d0**idec(id1,2,2)
     +                  +d2*dble(igrns(it,id2,2,2))
     +                     *10.d0**idec(id2,2,2))
                ut=ut+shss*(d1*dble(igrns(it,id1,3,2))
     +                     *10.d0**idec(id1,3,2)
     +                  +d2*dble(igrns(it,id2,3,2))
     +                     *10.d0**idec(id2,3,2))
c
                ezz=ezz+psss*(d1*dble(igrns(it,id1,4,2))
     +                      *10.d0**idec(id1,4,2)
     +                   +d2*dble(igrns(it,id2,4,2))
     +                      *10.d0**idec(id2,4,2))
                err=err+psss*(d1*dble(igrns(it,id1,5,2))
     +                      *10.d0**idec(id1,5,2)
     +                   +d2*dble(igrns(it,id2,5,2))
     +                      *10.d0**idec(id2,5,2))
                ett=ett+psss*(d1*dble(igrns(it,id1,6,2))
     +                      *10.d0**idec(id1,6,2)
     +                   +d2*dble(igrns(it,id2,6,2))
     +                      *10.d0**idec(id2,6,2))
                ezr=ezr+psss*(d1*dble(igrns(it,id1,7,2))
     +                      *10.d0**idec(id1,7,2)
     +                   +d2*dble(igrns(it,id2,7,2))
     +                      *10.d0**idec(id2,7,2))
                ert=ert+shss*(d1*dble(igrns(it,id1,8,2))
     +                      *10.d0**idec(id1,8,2)
     +                   +d2*dble(igrns(it,id2,8,2))
     +                      *10.d0**idec(id2,8,2))
                etz=etz+shss*(d1*dble(igrns(it,id1,9,2))
     +                      *10.d0**idec(id1,9,2)
     +                   +d2*dble(igrns(it,id2,9,2))
     +                      *10.d0**idec(id2,9,2))
c
                tr=tr+psss*(d1*dble(igrns(it,id1,10,2))
     +                     *10.d0**idec(id1,10,2)
     +                  +d2*dble(igrns(it,id2,10,2))
     +                     *10.d0**idec(id2,10,2))
                tt=tt+shss*(d1*dble(igrns(it,id1,11,2))
     +                     *10.d0**idec(id1,11,2)
     +                  +d2*dble(igrns(it,id2,11,2))
     +                     *10.d0**idec(id2,11,2))
c
                gd=gd+psss*(d1*dble(igrns(it,id1,12,2))
     +                     *10.d0**idec(id1,12,2)
     +                  +d2*dble(igrns(it,id2,12,2))
     +                     *10.d0**idec(id2,12,2))
c
                gr=gr+psss*(d1*dble(igrns(it,id1,13,2))
     +                     *10.d0**idec(id1,13,2)
     +                  +d2*dble(igrns(it,id2,13,2))
     +                     *10.d0**idec(id2,13,2))
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c               contributions from the dip-slip components
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                uz=uz+psds*(d1*dble(igrns(it,id1,1,3))
     +                        *10.d0**idec(id1,1,3)
     +                     +d2*dble(igrns(it,id2,1,3))
     +                        *10.d0**idec(id2,1,3))
                ur=ur+psds*(d1*dble(igrns(it,id1,2,3))
     +                        *10.d0**idec(id1,2,3)
     +                     +d2*dble(igrns(it,id2,2,3))
     +                        *10.d0**idec(id2,2,3))
                ut=ut+shds*(d1*dble(igrns(it,id1,3,3))
     +                        *10.d0**idec(id1,3,3)
     +                     +d2*dble(igrns(it,id2,3,3))
     +                        *10.d0**idec(id2,3,3))
c
                ezz=ezz+psds*(d1*dble(igrns(it,id1,4,3))
     +                          *10.d0**idec(id1,4,3)
     +                       +d2*dble(igrns(it,id2,4,3))
     +                          *10.d0**idec(id2,4,3))
                err=err+psds*(d1*dble(igrns(it,id1,5,3))
     +                          *10.d0**idec(id1,5,3)
     +                       +d2*dble(igrns(it,id2,5,3))
     +                          *10.d0**idec(id2,5,3))
                ett=ett+psds*(d1*dble(igrns(it,id1,6,3))
     +                          *10.d0**idec(id1,6,3)
     +                       +d2*dble(igrns(it,id2,6,3))
     +                          *10.d0**idec(id2,6,3))
                ezr=ezr+psds*(d1*dble(igrns(it,id1,7,3))
     +                          *10.d0**idec(id1,7,3)
     +                       +d2*dble(igrns(it,id2,7,3))
     +                          *10.d0**idec(id2,7,3))
                ert=ert+shds*(d1*dble(igrns(it,id1,8,3))
     +                          *10.d0**idec(id1,8,3)
     +                       +d2*dble(igrns(it,id2,8,3))
     +                          *10.d0**idec(id2,8,3))
                etz=etz+shds*(d1*dble(igrns(it,id1,9,3))
     +                          *10.d0**idec(id1,9,3)
     +                       +d2*dble(igrns(it,id2,9,3))
     +                          *10.d0**idec(id2,9,3))
c
                tr=tr+psds*(d1*dble(igrns(it,id1,10,3))
     +                        *10.d0**idec(id1,10,3)
     +                     +d2*dble(igrns(it,id2,10,3))
     +                        *10.d0**idec(id2,10,3))
                tt=tt+shds*(d1*dble(igrns(it,id1,11,3))
     +                        *10.d0**idec(id1,11,3)
     +                     +d2*dble(igrns(it,id2,11,3))
     +                        *10.d0**idec(id2,11,3))
c
                gd=gd+psds*(d1*dble(igrns(it,id1,12,3))
     +                        *10.d0**idec(id1,12,3)
     +                     +d2*dble(igrns(it,id2,12,3))
     +                        *10.d0**idec(id2,12,3))
c
                gr=gr+psds*(d1*dble(igrns(it,id1,13,3))
     +                        *10.d0**idec(id1,13,3)
     +                     +d2*dble(igrns(it,id2,13,3))
     +                        *10.d0**idec(id2,13,3))
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c               contributions from the clvd components
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                uz=uz+pscl*(d1*dble(igrns(it,id1,1,4))
     +                        *10.d0**idec(id1,1,4)
     +                     +d2*dble(igrns(it,id2,1,4))
     +                        *10.d0**idec(id2,1,4))
                ur=ur+pscl*(d1*dble(igrns(it,id1,2,4))
     +                        *10.d0**idec(id1,2,4)
     +                     +d2*dble(igrns(it,id2,2,4))
     +                        *10.d0**idec(id2,2,4))
c
                ezz=ezz+pscl*(d1*dble(igrns(it,id1,4,4))
     +                          *10.d0**idec(id1,4,4)
     +                       +d2*dble(igrns(it,id2,4,4))
     +                          *10.d0**idec(id2,4,4))
                err=err+pscl*(d1*dble(igrns(it,id1,5,4))
     +                          *10.d0**idec(id1,5,4)
     +                       +d2*dble(igrns(it,id2,5,4))
     +                          *10.d0**idec(id2,5,4))
                ett=ett+pscl*(d1*dble(igrns(it,id1,6,4))
     +                          *10.d0**idec(id1,6,4)
     +                       +d2*dble(igrns(it,id2,6,4))
     +                          *10.d0**idec(id2,6,4))
                ezr=ezr+pscl*(d1*dble(igrns(it,id1,7,4))
     +                          *10.d0**idec(id1,7,4)
     +                       +d2*dble(igrns(it,id2,7,4))
     +                          *10.d0**idec(id2,7,4))
c
                tr=tr+pscl*(d1*dble(igrns(it,id1,10,4))
     +                        *10.d0**idec(id1,10,4)
     +                     +d2*dble(igrns(it,id2,10,4))
     +                        *10.d0**idec(id2,10,4))
                gd=gd+pscl*(d1*dble(igrns(it,id1,12,4))
     +                        *10.d0**idec(id1,12,4)
     +                     +d2*dble(igrns(it,id2,12,4))
     +                        *10.d0**idec(id2,12,4))
c
                gr=gr+pscl*(d1*dble(igrns(it,id1,13,4))
     +                        *10.d0**idec(id1,13,4)
     +                     +d2*dble(igrns(it,id2,13,4))
     +                        *10.d0**idec(id2,13,4))
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
                obs(itr,irec,1)=obs(itr,irec,1)+ur*co-ut*si
                obs(itr,irec,2)=obs(itr,irec,2)+ur*si+ut*co
                obs(itr,irec,3)=obs(itr,irec,3)+uz
c
                obs(itr,irec,4)=obs(itr,irec,4)+err*co*co
     +                                       +ett*si*si-ert*si2
                obs(itr,irec,5)=obs(itr,irec,5)+err*si*si
     +                                       +ett*co*co+ert*si2
                obs(itr,irec,6)=obs(itr,irec,6)+ezz
                obs(itr,irec,7)=obs(itr,irec,7)+0.5d0*(err-ett)*si2
     +                                       +ert*co2
                obs(itr,irec,8)=obs(itr,irec,8)+ezr*si+etz*co
                obs(itr,irec,9)=obs(itr,irec,9)+ezr*co-etz*si
c
                obs(itr,irec,10)=obs(itr,irec,10)+tr*co-tt*si
                obs(itr,irec,11)=obs(itr,irec,11)+tr*si+tt*co
                obs(itr,irec,12)=obs(itr,irec,12)+gd
                obs(itr,irec,13)=obs(itr,irec,13)+gr
              enddo
            endif
          enddo
        enddo
        if(nsmall.gt.0)then
          nwarn=nwarn+nsmall
          write(nsmallc,888)nsmall
 888      format(i5)
c         call mexPrintf(' Warning: '//nsmallc//'too small distances'
c    &                      //' exceed the Green-function coverage!\n')
        endif
        if(nlarge.gt.0)then
          nwarn=nwarn+nlarge
          write(nlargec,888)nlarge
c         call mexPrintf(' Warning: '//nlargec//' too large distances'
c    &                     //' exceed the Green-function coverage!\n')
        endif
      enddo
c
      do istp=1,4
        do i=1,13
          if(select(i,istp))close(unit(i,istp))
        enddo
      enddo
      return
      end
