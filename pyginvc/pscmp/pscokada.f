      subroutine pscokada(ns,nrec,nzs1,nzs2)
      implicit none
c
      integer ns,nrec,nzs1,nzs2
c
      include 'pscglob.h'
c
c     Last modified: Potsdam, Nov, 2003, by R. Wang
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c     from Okada's subroutine DC3D0:
c
      INTEGER IRET
      REAL*4 ALPHA,X,Y,Z,DEPTH,DIPS,
     &         UX,UY,UZ,UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ
c
c     more from Okada's subroutine DC3D:
c
      REAL*4 AL1,AL2,AW1,AW2,DISL1,DISL2,DISL3
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c     LOCAL CONSTANTS
c     ===============
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      double precision eps,eps2
      data eps,eps2/1.0d-03,1.0d-06/
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c     LOCAL WORK SPACES
c     =================
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      integer i,is,irec,ixs,iys,nxs,nys,ixsys
      double precision st,di,ra
      double precision dr,xs,ys,xs1,xs2,ys1,ys2,dxs,dys
      double precision slparea,opnarea,dmsum
      double precision csst,ssst,csra,ssra,csdi,ssdi
      double precision cs2st,ss2st,eii
      double precision dm(NPSMAX),disp0(3),tilt0(2),strain0(6)
c
      double precision PI
      data PI/3.14159265358979d0/
c
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c     PROCESSING
c     ==========
c~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c
c     receiver and source independent variables
c
      ALPHA=sngl((larec+murec)/(larec+2.d0*murec))
      Z=-sngl(zrec)
      dr=(r2-r1)/dble(nr-1)
c
      do irec=1,nrec
c
c	initialization
c
c       obs(1)=' Ux'
c       obs(2)=' Uy'
c       obs(3)=' Uz'
c       obs(4)='Sxx'
c       obs(5)='Syy'
c       obs(6)='SZZ'
c       obs(7)='Sxy'
c       obs(8)='Syz'
c       obs(9)='Szx'
c       obs(10)=' Tx'
c       obs(11)=' Ty'
c       obs(12)=' Gd'
c       obs(13)=' Gr'
c
        do i=1,13
          obs(1,irec,i)=0.d0
        enddo
      enddo
      do is=1,ns
        st=strike(is)*DEG2RAD
        csst=dcos(st)
        ssst=dsin(st)
        cs2st=dcos(2.d0*st)
        ss2st=dsin(2.d0*st)
c
	di=dip(is)*DEG2RAD
        csdi=dcos(di)
        ssdi=dsin(di)
c
        ra=rake(is)*DEG2RAD
        csra=dcos(ra)
        ssra=dsin(ra)
c
        slparea=slip(is)
        opnarea=openning(is)
c
        if(length(is).eq.0.d0)then
C code modified version 4.4
          dxs=1.d-02*dr
c         dxs=1.d-06*dr
        else if(taplt(is).gt.0.d0.or.taprt(is).gt.0.d0)then
          dxs=dr
        else
          dxs=dmax1(dr,length(is))
        endif
        nxs=max0(1,idnint(length(is)/dxs))
        if(length(is).eq.0.d0)then
          slparea=slparea/dxs
          opnarea=opnarea/dxs
        else
          dxs=length(is)/dble(nxs)
        endif
c
        if(width(is).eq.0.d0)then
C code modified version 4.4
          dys=1.d-02*dr
c         dys=1.d-06*dr
        else if(tapup(is).gt.0.d0.or.taplw(is).gt.0.d0)then
          dys=dr
        else
          dys=dmax1(dr,width(is))
        endif
        nys=max0(1,idnint(width(is)/dys))
        if(width(is).eq.0.d0)then
          slparea=slparea/dys
          opnarea=opnarea/dys
        else
          dys=width(is)/dble(nys)
        endif
c
        xs1=taplt(is)
        xs2=taprt(is)
        ys1=tapup(is)
        ys2=taplw(is)
c
        dmsum=0.d0
        do iys=1,nys
          ys=dys*(dble(iys)-0.5d0)
          do ixs=1,nxs
            xs=dxs*(dble(ixs)-0.5d0)
            ixsys=(iys-1)*nxs+ixs
c
            dm(ixsys)=1.d0
c
            if(xs1.gt.0.d0.and.xs.lt.xs1)then
              dm(ixsys)=dm(ixsys)*dcos(0.5d0*PI*(xs1-xs)/xs1)
            endif
            if(xs2.gt.0.d0.and.xs.gt.length(is)-xs2)then
              dm(ixsys)=dm(ixsys)
     &                 *dcos(0.5d0*PI*(xs-length(is)+xs2)/xs2)
            endif
            if(ys1.gt.0.d0.and.ys.lt.ys1)then
              dm(ixsys)=dm(ixsys)*dcos(0.5d0*PI*(ys1-ys)/ys1)
            endif
            if(ys2.gt.0.d0.and.ys.gt.width(is)-ys2)then
              dm(ixsys)=dm(ixsys)
     &                 *dcos(0.5d0*PI*(ys-width(is)+ys2)/ys2)
            endif
c
            dmsum=dmsum+dm(ixsys)
c
          enddo
        enddo
        if(dmsum.ne.dble(nxs*nys))then
          do ixsys=1,nxs*nys
            dm(ixsys)=dm(ixsys)*dble(nxs*nys)/dmsum
          enddo
        endif
c
        DEPTH=sngl(zref(is))
        DIPS=sngl(dip(is))
c
        do iys=1,nys
          ys=dys*(dble(iys)-0.5d0)
          do ixs=1,nxs
            xs=dxs*(dble(ixs)-0.5d0)
            ixsys=(iys-1)*nxs+ixs
c
            DISL1=sngl(slparea*csra*dm(ixsys))
            DISL2=sngl(slparea*ssra*dm(ixsys))
            DISL3=sngl(opnarea*dm(ixsys))
c
c           for extended source
c
            AL1=sngl(xs-0.5d0*dxs)
            AL2=sngl(xs+0.5d0*dxs)
            AW2=-sngl(ys-0.5d0*dys)
            AW1=-sngl(ys+0.5d0*dys)
            do irec=1,nrec
c
c             transform from Aki's to Okada's system
c
              X=sngl((xrec(irec)-xref(is))*csst
     &         +(yrec(irec)-yref(is))*ssst)
              Y=sngl((xrec(irec)-xref(is))*ssst
     &         -(yrec(irec)-yref(is))*csst)
              IRET=1
              call DC3D(ALPHA,X,Y,Z,DEPTH,DIPS,AL1,AL2,AW1,AW2,
     &               DISL1,DISL2,DISL3,UX,UY,UZ,
     &               UXX,UYX,UZX,UXY,UYY,UZY,UXZ,UYZ,UZZ,IRET)
c              if(IRET.eq.1)then
c                stop ' There is a problem in Okada subroutine!'
c              endif
c
c             transform from Okada's to Aki's system
c
              disp0(1)=dble(UX)*csst+dble(UY)*ssst
              disp0(2)=dble(UX)*ssst-dble(UY)*csst
              disp0(3)=-dble(UZ)
c
              tilt0(1)=-(dble(UXZ)*csst+dble(UYZ)*ssst)
              tilt0(2)=-(dble(UXZ)*ssst-dble(UYZ)*csst)
c
              strain0(1)=dble(UXX)*csst*csst+dble(UYY)*ssst*ssst
     &                  +0.5d0*dble(UXY+UYX)*ss2st
              strain0(2)=dble(UXX)*ssst*ssst+dble(UYY)*csst*csst
     &                  -0.5d0*dble(UXY+UYX)*ss2st
              strain0(3)=dble(UZZ)
              strain0(4)=0.5d0*dble(UXX-UYY)*ss2st
     &                  -0.5d0*dble(UXY+UYX)*cs2st
              strain0(5)=-0.5d0*dble(UZX+UXZ)*ssst
     &                   +0.5d0*dble(UYZ+UZY)*csst
              strain0(6)=-0.5d0*dble(UZX+UXZ)*csst
     &                   -0.5d0*dble(UYZ+UZY)*ssst
c
              do i=1,3
                obs(1,irec,i)=obs(1,irec,i)+disp0(i)
              enddo
              do i=4,9
                obs(1,irec,i)=obs(1,irec,i)+strain0(i-3)
              enddo
              do i=10,11
                obs(1,irec,i)=obs(1,irec,i)+tilt0(i-9)
              enddo
            enddo
          enddo
        enddo
C code modified version 4.4
      enddo
      do irec=1,nrec
        do i=1,6
          strain0(i)=obs(1,irec,i+3)
        enddo
        eii=strain0(1)+strain0(2)+strain0(3)
        obs(1,irec,4)=larec*eii+2.d0*murec*strain0(1)
        obs(1,irec,5)=larec*eii+2.d0*murec*strain0(2)
        obs(1,irec,6)=larec*eii+2.d0*murec*strain0(3)
        obs(1,irec,7)=2.d0*murec*strain0(4)
        obs(1,irec,8)=2.d0*murec*strain0(5)
        obs(1,irec,9)=2.d0*murec*strain0(6)
        if(zrec.eq.0.d0)then
          obs(1,irec,6)=0.d0
          obs(1,irec,8)=0.d0
          obs(1,irec,9)=0.d0
        endif
        obs(1,irec,13)=(2.d0*G0/REARTH)*obs(1,irec,3)
      enddo
c
      return
      end
