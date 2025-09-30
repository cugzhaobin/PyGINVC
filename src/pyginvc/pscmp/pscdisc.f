      subroutine pscdisc(ns,npsum)
      implicit none
c
c     Last modified: Potsdam, April, 2003, by R. Wang
c
      integer ns,npsum
c
c     returned outputs:
c     outputs through common blocks
c
      include 'pscglob.h'
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c     LOCAL WORK SPACES
c     =================
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer is,izs,ix,iy,ixy,nx,ny,ips
      double precision dr,x,y,x1,x2,y1,y2,dx,dy,dzs,st,di,ra
      double precision slparea,opnarea
      double precision pz0,dmslp,dmopn_cl,dmopn_ep,dmsum
      double precision sss,sss2,ss2s,ssd,ssd2,ss2d,ssr,ssr2,ss2r
      double precision css,css2,cs2s,csd,csd2,cs2d,csr,csr2,cs2r
      double precision sm(3,3),om(3,3),dm(NPSMAX)
c
      double precision PI
      data PI/3.14159265358979d0/
      character isc*2,ixyc*6,npsumc*7
c
c     call mexPrintf(' ... discretise rectangular plane sources:\n')
      dr=(r2-r1)/dble(nr-1)
      do izs=1,nzs
        nps(izs)=0
      enddo
      if(nzs.gt.1)then
        dzs=(zs2-zs1)/dble(nzs-1)
      else
        dzs=dr
      endif
c
      npsum=0
c
      do is=1,ns
c
        st=strike(is)*DEG2RAD
        di=dip(is)*DEG2RAD
        ra=rake(is)*DEG2RAD
c
        sss=dsin(st)
        sss2=sss*sss
        ss2s=dsin(2.d0*st)
        ssd=dsin(di)
        ssd2=ssd*ssd
        ss2d=dsin(2.d0*di)
        ssr=dsin(ra)
        ssr2=ssr*ssr
        ss2r=dsin(2.d0*ra)
c
        css=dcos(st)
        css2=css*css
        cs2s=dcos(2.d0*st)
        csd=dcos(di)
        csd2=csd*csd
        cs2d=dcos(2.d0*di)
        csr=dcos(ra)
        csr2=csr*csr
        cs2r=dcos(2.d0*ra)
c
        sm(1,1)=-ssd*csr*ss2s-ss2d*ssr*sss2
        sm(2,2)= ssd*csr*ss2s-ss2d*ssr*css2
        sm(3,3)=-(sm(1,1)+sm(2,2))
        sm(1,2)= ssd*csr*cs2s+0.5d0*ss2d*ssr*ss2s
        sm(2,1)=sm(1,2)
        sm(2,3)=-csd*csr*sss+cs2d*ssr*css
        sm(3,2)=sm(2,3)
        sm(3,1)=-csd*csr*css-cs2d*ssr*sss
        sm(1,3)=sm(3,1)
c
c       openning => explosion (ep) + clvd (om)
c
        om(1,1)=-0.5d0+1.5d0*sss2*ssd2
        om(2,2)=-0.5d0+1.5d0*css2*ssd2
        om(3,3)=-(om(1,1)+om(2,2))
        om(1,2)=-1.5d0*sss*css*ssd2
        om(2,1)=om(1,2)
        om(2,3)=-1.5d0*css*ssd*csd
        om(3,2)=om(2,3)
        om(3,1)= 1.5d0*sss*ssd*csd
        om(1,3)=om(3,1)
c
        slparea=slip(is)
        opnarea=openning(is)
c
        dx=dr
        nx=max0(1,idnint(length(is)/dx))
        if(length(is).gt.0.d0)then
          dx=length(is)/dble(nx)
          slparea=slparea*length(is)
          opnarea=opnarea*length(is)
        else
          dx=0.d0
        endif
c
C code modified version 4.4
c
        if(dabs(dzs*ssd).gt.0.d0)then
          dy=dmin1(dr,dzs/ssd)
        else
          dy=dr
        endif
c       dy=dmax1(dr,dzs)
        ny=max0(1,idnint(width(is)/dy))
        if(width(is).gt.0.d0)then
          dy=width(is)/dble(ny)
          slparea=slparea*width(is)
          opnarea=opnarea*width(is)
        else
          dy=0.d0
        endif
c
        if(nx*ny.gt.NPSMAX)then
c         call mexErrMsgTxt('Error: NPSMAX defined too small!')
        endif
c
C code modified version 4.4
c
        if(zref(is)+1.5d0*dy*ssd.lt.zs1.or.
     &     zref(is)+(width(is)-1.5d0*dy)*ssd.gt.zs2)then
          nwarn=nwarn+1
c         call mexPrintf('Warning in pscdisc: the fault')
c         call mexPrintf(' plane in the depth range not')
c         call mexPrintf(' covered by Green funcions!\n')
        endif
c
c       for taper
c
        x1=taplt(is)
        x2=taprt(is)
        y1=tapup(is)
        y2=taplw(is)
c
        dmsum=0.d0
c
        do iy=1,ny
          y=dy*(dble(iy)-0.5d0)
          do ix=1,nx
            x=dx*(dble(ix)-0.5d0)
            ixy=(iy-1)*nx+ix
            npsum=npsum+1
c
            dm(ixy)=1.d0
c
            if(x1.gt.0.d0.and.x.lt.x1)then
              dm(ixy)=dm(ixy)*dcos(0.5d0*PI*(x1-x)/x1)
            endif
            if(x2.gt.0.d0.and.x.gt.length(is)-x2)then
              dm(ixy)=dm(ixy)*dcos(0.5d0*PI*(x-length(is)+x2)/x2)
            endif
            if(y1.gt.0.d0.and.y.lt.y1)then
              dm(ixy)=dm(ixy)*dcos(0.5d0*PI*(y1-y)/y1)
            endif
            if(y2.gt.0.d0.and.y.gt.width(is)-y2)then
              dm(ixy)=dm(ixy)*dcos(0.5d0*PI*(y-width(is)+y2)/y2)
            endif
c
            dmsum=dmsum+dm(ixy)
c
          enddo
        enddo
c
        do iy=1,ny
          y=dy*(dble(iy)-0.5d0)
          pz0=zref(is)+y*ssd
          if(dzs.gt.0.d0)then
            izs=idnint((pz0-zs1)/dzs)+1
            izs=max0(min0(izs,nzs),1)
          else
            izs=1
          endif
          do ix=1,nx
            x=dx*(dble(ix)-0.5d0)
            nps(izs)=nps(izs)+1
            ips=nps(izs)
            if(nps(izs).gt.NPSMAX)then
c             call mexErrMsgTxt('Error: NPSMAX too small defined!')
            endif
            px(ips,izs)=xref(is)+x*css-y*csd*sss
            py(ips,izs)=yref(is)+x*sss+y*csd*css
            pz(ips,izs)=pz0
c
c           1 = weight for strike-slip: m12=m21=1;
c           2 = weight for dip-slip: m13=m31=1
c           3 = weight for clvd: m33=-m11=-m22=1
c           4 = weight for 45 deg strike-slip: m11=-m22=1
c           5 = weight for 45 deg dip-slip: m23=m32=1
c           6 = weight for explosion: m11=m22=m33=1
c
            ixy=(iy-1)*nx+ix
            dmslp=dm(ixy)*slparea/dmsum
            dmopn_ep=dm(ixy)*opnarea/dmsum
            dmopn_cl=dmopn_ep*4.d0/3.d0
c
            pmwei(1,ips,izs)=sm(1,2)*dmslp+om(1,2)*dmopn_cl
            pmwei(2,ips,izs)=sm(3,1)*dmslp+om(3,1)*dmopn_cl
            pmwei(3,ips,izs)=sm(3,3)*dmslp+om(3,3)*dmopn_cl
            pmwei(4,ips,izs)=0.5d0*(sm(1,1)-sm(2,2))*dmslp
     &                      +0.5d0*(om(1,1)-om(2,2))*dmopn_cl
            pmwei(5,ips,izs)=sm(2,3)*dmslp+om(2,3)*dmopn_cl
c
            pmwei(6,ips,izs)=dmopn_ep
          enddo
        enddo
        write(isc,886)is
 886    format(i2)
        write(ixyc,887)ixy
 887    format(i6)
c       call mexPrintf(' the '//isc//'. rectangle => '//ixyc//
c    &     ' point sources.\n')
      enddo
c     call mexPrintf('----------------------------------------------\n')
      write(npsumc,888)npsum
 888  format(i7)
c     call mexPrintf(' total number of point sources: '//npsumc//'\n')
      return
      end
