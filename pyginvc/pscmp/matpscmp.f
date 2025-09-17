C The gateway routine
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
      implicit none
      include 'pscglob.h'

      integer plhs(*), prhs(*)
      integer nlhs, nrhs
      integer strlen, status
      character*80 grndir

      integer nmax, mmax
      parameter (mmax = 250)
      parameter (nmax = 3000)

c     Declare mex functions
      integer mxisnumeric, mxgetm, mxgetn, mxgetpr, mxcreatedoublematrix
      integer mxIsChar,mxGetString

      integer disgeom_pr
      integer cxyrec_pr
      integer Ux_pr,Uy_pr,Uz_pr
      integer Sxx_pr,Syy_pr,Szz_pr
      integer Sxy_pr,Syz_pr,Szx_pr
      integer Tx_pr,Ty_pr
      integer Gd_pr,Gr_pr
      integer X_pr,Y_pr
      integer m, n, size, d_m, d_n, i, j, c_m, c_n
      real*8 disgeom(14,10), dgeom10(10), dgeom14(14)
      real*8 cxyrec(2,nmax), crec(2)
      real*8 Ux(mmax,nmax), Uy(mmax,nmax), Uz(mmax,nmax)
      real*8 Sxx(mmax,nmax), Syy(mmax,nmax), Szz(mmax,nmax)
      real*8 Sxy(mmax,nmax), Syz(mmax,nmax), Szx(mmax,nmax)
      real*8 Tx(mmax,nmax), Ty(mmax,nmax)
      real*8 Gd(mmax,nmax), Gr(mmax,nmax)
      real*8 X(1,nmax), Y(1,nmax)
      real*8 midpoint_x, midpoint_y, depth_shift, horiz_shift

      character dischar*12


C Check for proper number of arguments.
      if (nrhs .ne. 2 .and. nrhs .ne. 3) then
         call mexErrMsgTxt('Two or three inputs required.')
      endif
      if (nlhs .ne. 15) then
         call mexErrMsgTxt('Fifteen outputs required.')
      endif

C Optional Input #3 must be a string.
      if (nrhs .eq. 3) then
      if (mxIsChar(prhs(3)) .ne. 1) then
         call mexErrMsgTxt('Input #3 must be a directory name.')
      endif

C Input #3 must be a row vector.
      if (mxGetM(prhs(3)) .ne. 1) then
         call mexErrMsgTxt('Input must be a row vector.')
      endif

C Get the length of input #3 string.
      strlen = mxGetM(prhs(3))*mxGetN(prhs(3))

C Get the string contents.
      status = mxGetString(prhs(3), grndir, strlen)

C Check if mxGetString is successful.
      if (status .ne. 0) then
         call mexErrMsgTxt('String length must be less than 80.')
      endif
      else
C Default Green's Function directory to the current psgrnfcts directory.
         grndir = './psgrnfcts/'
      endif

C Input #1 must be numeric.
      if (mxIsNumeric(prhs(1)) .ne. 1) then
         call mexErrMsgTxt('Input #1 is not a numeric array.')
      endif

C Get the size of the input matrix.
      d_m = mxGetM(prhs(1))
      d_n = mxGetN(prhs(1))
      if (d_m .ne. 10 .and. d_m .ne. 14) then
      call mexErrMsgTxt('disgeom array must be 10xn or 14xn array')
      endif

C Load the input matrix into a Fortran array.
C disgeom = [850.0;130.0;4.0;20.0;8.0;-75.89;-44.66;-4.4;16.42;0.0;0.0;0.0;0.0;0.0]
      disgeom_pr = mxGetPr(prhs(1))
      if(d_m .eq.14) then
C disgeom = [850.0;130.0;4.0;20.0;8.0;-75.89;-44.66;-4.4;16.42;0.0;0.0;0.0;0.0;0.0]
      do i = 1, d_n
        call mxCopyPtrToReal8(disgeom_pr+(i-1)*112, dgeom14, 14)
c
c         Adjust X and Y coordinates of fault to reflect midpoint of fault
c         so that this code matches the standard disloc code
c
        midpoint_x = (dgeom14(1)/2)*cos( (90-dgeom14(5))*DEG2RAD )
        midpoint_y = (dgeom14(1)/2)*sin( (90-dgeom14(5))*DEG2RAD )
        dgeom14(6) = dgeom14(6) - midpoint_x
        dgeom14(7) = dgeom14(7) - midpoint_y

        disgeom(1,i) = dgeom14(1)
        disgeom(2,i) = dgeom14(2)
c
c         Adjust depth to refer to the TOP of the fault (Wang et al. convention)
c         instead of the BOTTOM of the fault (Okada convention). For a fault with
c         a dip not equal to 90 degrees, this also shifts the horizontal position
c         of the fault.
c
        depth_shift = dgeom14(2)*sin(dgeom14(4)*DEG2RAD)
        horiz_shift = dgeom14(2)*cos(dgeom14(4)*DEG2RAD)
        disgeom(3,i) = dgeom14(3) - depth_shift
        dgeom14(6) = dgeom14(6)+horiz_shift*sin(DEG2RAD*(90-dgeom14(5)))
        dgeom14(7) = dgeom14(7)-horiz_shift*cos(DEG2RAD*(90-dgeom14(5)))

        disgeom(4,i) = dgeom14(4)
        disgeom(5,i) = dgeom14(5)
        disgeom(6,i) = dgeom14(6)
        disgeom(7,i) = dgeom14(7)
        disgeom(8,i) = dgeom14(8)
        disgeom(9,i) = dgeom14(9)
        disgeom(10,i) = dgeom14(10)
        disgeom(11,i) = dgeom14(11)
        disgeom(12,i) = dgeom14(12)
        disgeom(13,i) = dgeom14(13)
        disgeom(14,i) = dgeom14(14)
      enddo

      else if (d_m .eq. 10) then
C disgeom = [850.0;130.0;4.0;20.0;8.0;-75.89;-44.66;-4.4;16.42;0.0]
      do i = 1, d_n
        call mxCopyPtrToReal8(disgeom_pr+(i-1)*80, dgeom10, 10)

c
c         Adjust X and Y coordinates of fault to reflect midpoint of fault
c         so that this code matches the standard disloc code
c
        midpoint_x = (dgeom10(1)/2)*cos( (90-dgeom10(5))*DEG2RAD )
        midpoint_y = (dgeom10(1)/2)*sin( (90-dgeom10(5))*DEG2RAD )
        dgeom10(6) = dgeom10(6) - midpoint_x
        dgeom10(7) = dgeom10(7) - midpoint_y

        disgeom(1,i) = dgeom10(1)
        disgeom(2,i) = dgeom10(2)
c
c         Adjust depth to refer to the TOP of the fault (Wang et al. convention)
c         instead of the BOTTOM of the fault (Okada convention). For a fault with
c         a dip not equal to 90 degrees, this also shifts the horizontal position
c         of the fault.
c
        depth_shift = dgeom10(2)*sin(dgeom10(4)*DEG2RAD)
        horiz_shift = dgeom10(2)*cos(dgeom10(4)*DEG2RAD)
        disgeom(3,i) = dgeom10(3) - depth_shift
        dgeom10(6) = dgeom10(6)-horiz_shift*sin(DEG2RAD*(90-dgeom10(5)))
        dgeom10(7) = dgeom10(7)+horiz_shift*cos(DEG2RAD*(90-dgeom10(5)))

        disgeom(4,i) = dgeom10(4)
        disgeom(5,i) = dgeom10(5)
        disgeom(6,i) = dgeom10(6)
        disgeom(7,i) = dgeom10(7)
        disgeom(8,i) = dgeom10(8)
        disgeom(9,i) = dgeom10(9)
        disgeom(10,i) = dgeom10(10)

C Set the tapers to zero by default
        disgeom(11,i) = 0.0
        disgeom(12,i) = 0.0
        disgeom(13,i) = 0.0
        disgeom(14,i) = 0.0
      enddo
      endif

C Input #2 must be numeric.
      if (mxIsNumeric(prhs(2)) .ne. 1) then
         call mexErrMsgTxt('Input #2 is not a numeric array.')
      endif

C Get the size of the input matrix.
      c_m = mxGetM(prhs(2))
      c_n = mxGetN(prhs(2))
      if (c_m .ne. 2) then
      call mexErrMsgTxt('enter matrix as two dimensional array x x;y y')
      endif

C Load the input matrix into a Fortran array.
C c = [5.0 25.0 100.0 0.0 0.0 0.0;50.0 50.0 50.0 1.5 3.0 10.0]
c     First load the X,Y coordinates, correcting for coordinate system
c     Must swap the X and Y coordinates because Wang et al use Aki's
c     coordinate system convention (X, Y, Z) = (North, East, Down)
c
      cxyrec_pr = mxGetPr(prhs(2))
      do i = 1, c_n
        call mxCopyPtrToReal8(cxyrec_pr+(i-1)*16,crec,2)
c            swap X and Y (see above)
        cxyrec(1,i) = crec(2)
        cxyrec(2,i) = crec(1)
      enddo

C Create 13 matrices for return data and an X and Y array.
      m = mmax
      n = nmax
      size = m*n
      plhs(1) = mxCreateDoubleMatrix(m,n,0)
      plhs(2) = mxCreateDoubleMatrix(m,n,0)
      plhs(3) = mxCreateDoubleMatrix(m,n,0)
      plhs(4) = mxCreateDoubleMatrix(m,n,0)
      plhs(5) = mxCreateDoubleMatrix(m,n,0)
      plhs(6) = mxCreateDoubleMatrix(m,n,0)
      plhs(7) = mxCreateDoubleMatrix(m,n,0)
      plhs(8) = mxCreateDoubleMatrix(m,n,0)
      plhs(9) = mxCreateDoubleMatrix(m,n,0)
      plhs(10) = mxCreateDoubleMatrix(m,n,0)
      plhs(11) = mxCreateDoubleMatrix(m,n,0)
      plhs(12) = mxCreateDoubleMatrix(m,n,0)
      plhs(13) = mxCreateDoubleMatrix(m,n,0)
      plhs(14) = mxCreateDoubleMatrix(1,n,0)
      plhs(15) = mxCreateDoubleMatrix(1,n,0)

C Get output array pointer
      Ux_pr = mxGetPr(plhs(1))
      Uy_pr = mxGetPr(plhs(2))
      Uz_pr = mxGetPr(plhs(3))
      Sxx_pr = mxGetPr(plhs(4))
      Syy_pr = mxGetPr(plhs(5))
      Szz_pr = mxGetPr(plhs(6))
      Sxy_pr = mxGetPr(plhs(7))
      Syz_pr = mxGetPr(plhs(8))
      Szx_pr = mxGetPr(plhs(9))
      Tx_pr = mxGetPr(plhs(10))
      Ty_pr = mxGetPr(plhs(11))
      Gd_pr = mxGetPr(plhs(12))
      Gr_pr = mxGetPr(plhs(13))
      X_pr = mxGetPr(plhs(14))
      Y_pr = mxGetPr(plhs(15))

C Call the computational subroutine
      call pscmain(disgeom,d_n,cxyrec,c_n,
     * Ux,Uy,Uz,Sxx,Syy,Szz,Sxy,Syz,Szx,Tx,Ty,Gd,Gr,X,Y,grndir)

C Load the data into the output to MATLAB.
c     Need to swap x and y axes
c     Need to account for the vertical sign convention
c          multiply z components by -1 (except Szz)
c
c        Displacements: multiply z by -1
      call mxCopyReal8ToPtr(Ux, Uy_pr, size)
      call mxCopyReal8ToPtr(Uy, Ux_pr, size)
      do i = 1, c_m
         do j = 1, c_n
            Uz(i,j) = -1.d0*Uz(i,j)
         enddo
      enddo
      call mxCopyReal8ToPtr(Uz, Uz_pr, size)
     
c        Strains: multiply by -1 except Szz (by -1*-1 = 1)     
c     Need to swap x and y axes
      call mxCopyReal8ToPtr(Sxx, Syy_pr, size)
      call mxCopyReal8ToPtr(Syy, Sxx_pr, size)
      call mxCopyReal8ToPtr(Szz, Szz_pr, size)
      call mxCopyReal8ToPtr(Sxy, Sxy_pr, size)
      do i = 1, c_m
         do j = 1, c_n
            Syz(i,j) = -1.d0*Syz(i,j)
            Szx(i,j) = -1.d0*Szx(i,j)
         enddo
      enddo
      call mxCopyReal8ToPtr(Syz, Szx_pr, size)
      call mxCopyReal8ToPtr(Szx, Syz_pr, size)
      
c         Tilts: multiply by -1
c     Need to swap x and y axes
      do i = 1, c_m
         do j = 1, c_n
            Tx(i,j) = -1.d0*Tx(i,j)
            Ty(i,j) = -1.d0*Ty(i,j)
         enddo
      enddo
      call mxCopyReal8ToPtr(Tx, Ty_pr, size)
      call mxCopyReal8ToPtr(Ty, Tx_pr, size)
      call mxCopyReal8ToPtr(Gd, Gd_pr, size)
      call mxCopyReal8ToPtr(Gr, Gr_pr, size)

c         X and Y: swap axes because these variables reflect the 
c          Aki coordinate convention used internally
      call mxCopyReal8ToPtr(X, Y_pr, n)
      call mxCopyReal8ToPtr(Y, X_pr, n)

      return
      end
