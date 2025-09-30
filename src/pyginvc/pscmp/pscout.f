      subroutine pscout(nrec,outdir,itout,toutfile,
     &                  onlysc,nsc,itsc,scoutfile)
      implicit none
c
c     Last modified: Potsdam, April, 2003, by R. Wang
c
      integer nrec,nsc,itout(13),itsc(nsc)
      character*80 outdir,toutfile(13),scoutfile(nsc)
      logical onlysc
c
      include 'pscglob.h'
c
      integer i,it,isc,itr,j,inx,irec,lend,m
      double precision dt,t,los,losmax,xmax,ymax,losmin,xmin,ymin
      character*3 cmptxt(13)
      character*7 rtxt(NRECMAX)
      character*160 outfile
      character*3 iscc
c
c     DATA OUTPUT
c     ===========
c
      do lend=80,1,-1
        if(outdir(lend:lend).ne.' ')goto 100
      enddo
100   continue
c
      if(lend.lt.1)then
        call mexErrMsgTxt('Error in edcmain: wrong for output dirtory!')
      endif
c
      cmptxt(1)=' Ux'
      cmptxt(2)=' Uy'
      cmptxt(3)=' Uz'
      cmptxt(4)='Sxx'
      cmptxt(5)='Syy'
      cmptxt(6)='Szz'
      cmptxt(7)='Sxy'
      cmptxt(8)='Syz'
      cmptxt(9)='Szx'
      cmptxt(10)=' Tx'
      cmptxt(11)=' Ty'
      cmptxt(12)=' Gd'
      cmptxt(13)=' Gr'
      inx=int(alog10(0.1+float(nrec)))+1
      do irec=1,nrec
        m=irec
        do j=1,inx
          i=m/10**(inx-j)
          rtxt(irec)(j:j)=char(ichar('0')+i)
          m=m-i*10**(inx-j)
        enddo
      enddo
c
      if(nt.gt.1)then
        dt=twindow/dble(nt-1)
      else
        dt=twindow
      endif
c
      do i=1,13
        if(itout(i).eq.1)then
          outfile=outdir(1:lend)//toutfile(i)
          open(30,file=outfile,status='unknown')
          write(30,'(a12,$)')'    Time_day'
          do irec=1,nrec-1
            write(30,'(a12,$)')cmptxt(i)//rtxt(irec)(1:inx)
          enddo
          write(30,'(a12)')cmptxt(i)//rtxt(nrec)(1:inx)
          do it=1,nt
            t=dble(it-1)*dt
            write(30,'(E12.4,$)')t/DAY2SEC
            do irec=1,nrec-1
              write(30,'(E12.4,$)')obs(it,irec,i)
            enddo
            write(30,'(E12.4)')obs(it,nrec,i)
          enddo
          close(30)
        endif
      enddo
c
c     scenario outputs
c
      do isc=1,nsc
        if(itsc(isc).lt.1.or.itsc(isc).gt.nt)then
          nwarn=nwarn+1
          write(iscc,888)isc
  888     format(i3)
          call mexPrintf(' Warning in pecout: time point')
          call mexPrintf(' outside time window, no output')
          call mexPrintf(' for the'//iscc//'. scenario!\n')
          goto 500
        endif
        outfile=outdir(1:lend)//scoutfile(isc)
        open(30,file=outfile,status='unknown')
        if(latlon.eq.1)then
          write(30,'(a24,$)')'    Lat[deg]    Lon[deg]'
        else
          write(30,'(a24,$)')'       X[m]         Y[m]'
        endif
        do i=1,13
          write(30,'(a12,$)')cmptxt(i)
        enddo
        if(insar.eq.1)then
          write(30,'(a12)')'LOS[cm]'
        else
          write(30,'(a1)')' '
        endif
        if(onlysc)then
          itr=isc
        else
          itr=itsc(isc)
        endif
        losmax=100.d0*(xlos*obs(itr,1,1)
     &        +ylos*obs(itr,1,2)+zlos*obs(itr,1,3))
        losmin=losmax
        do irec=1,nrec
          los=100.d0*(xlos*obs(itr,irec,1)
     &       +ylos*obs(itr,irec,2)+zlos*obs(itr,irec,3))
          if(losmax.lt.los)then
            losmax=los
            if(latlon.eq.1)then
              xmax=xrec0(irec)
              ymax=yrec0(irec)
            else
              xmax=xrec(irec)
              ymax=yrec(irec)
            endif
          endif
          if(losmin.gt.los)then
            losmin=los
            if(latlon.eq.1)then
              xmin=xrec0(irec)
              ymin=yrec0(irec)
            else
              xmin=xrec(irec)
              ymin=yrec(irec)
            endif
          endif
          if(latlon.eq.1)then
	      if(insar.eq.1)then
              write(30,1000)xrec0(irec),yrec0(irec),
     &                      (obs(itr,irec,i),i=1,13),los
	      else
	        write(30,1000)xrec0(irec),yrec0(irec),
     &                      (obs(itr,irec,i),i=1,13)
	      endif
          else
	      if(insar.eq.1)then
              write(30,1000)xrec(irec),yrec(irec),
     &                      (obs(itr,irec,i),i=1,13),los
	      else
	        write(30,1000)xrec(irec),yrec(irec),
     &                      (obs(itr,irec,i),i=1,13)
	      endif
          endif
        enddo
        close(30)
500     continue
      enddo
c
1000  format(2f12.2,14E12.4)
      return
      end
