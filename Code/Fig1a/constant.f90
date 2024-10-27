


      Module global

       integer,parameter :: Ne=8 !Number the electrons in the system

       integer,parameter :: Nx=16 !System lattice size
       integer,parameter :: Ny=16 !4


       double precision,parameter :: Pi= 4.0d0*datan(1.0d0)
       double complex,parameter :: ii=(0.0d0,1.0d0)

      !double precision,parameter :: t1=1.0d0
      !double precision,parameter :: t2=1.0d0/(2.+sqrt(2.0))
      !double precision,parameter :: t3=1.0d0/(2.+2.*sqrt(2.0))
      !double precision,parameter :: phi=Pi/4.0
      !double precision,parameter :: M=0.0d0
      !double precision,parameter :: Hubb_U=1.0d0


      integer :: Total_Kx !=0
      integer :: Total_Ky !=3

      integer :: dimen !1044

      double precision,parameter :: EX=1.0E-8

      integer(kind=8),parameter :: MAXmatrix=750000000

      integer,parameter :: MAXbasis=350000

      integer,parameter :: MAXnonzero=5000

      integer,parameter :: OMPThreads=8


         integer :: countnum
         double complex,allocatable :: MM(:)
         !integer,allocatable :: MM(:)
         integer(kind=8),allocatable :: Mi(:)
         integer,allocatable :: Mj(:)


        integer,parameter :: ns0=Nx*Ny
        integer,parameter :: kx0=Nx,ky0=Ny,nst=ns0
        integer, parameter::nq1=21,nq2=21, nl=2, nk1=kx0,nk2=ky0


        integer,parameter :: mb = nq1*nq2*nl !98
        double complex,allocatable :: wf(:,:) !wf(1:mb,0:Nx*Ny-1)
        double precision,allocatable :: Ek(:,:) !Ek(1:mb,0:Nx*Ny-1)

        !complex*16 :: fq11(ns0,ns0,-nq1:nq1,-nq2:nq2)
        !complex*16 :: fq00(ns0,ns0,-nq1:nq1,-nq2:nq2)
        !complex*16 :: fq01(ns0,ns0,-nq1:nq1,-nq2:nq2)
        !complex*16 :: fq10(ns0,ns0,-nq1:nq1,-nq2:nq2)
        !complex*16 :: fq22(ns0,ns0,-nq1:nq1,-nq2:nq2)
        !complex*16 :: fq02(ns0,ns0,-nq1:nq1,-nq2:nq2)
        !complex*16 :: fq20(ns0,ns0,-nq1:nq1,-nq2:nq2)
        !complex*16 :: fq21(ns0,ns0,-nq1:nq1,-nq2:nq2)
        !complex*16 :: fq12(ns0,ns0,-nq1:nq1,-nq2:nq2)

        integer  mq0(2,2)
        integer  mg0(2,2)
        integer mnt_d(ns0, ns0)
        integer mnt(kx0*ky0*2,2), ind_mnt(0:kx0,0:ky0)

        double complex  glist(7), g1,g2

        double precision :: val_hf(ns0),valf(ns0,1:10),valf0(ns0)

        complex*16 :: vecf0(nq1*nq2*nl,ns0),vecf3(nq1*nq2*nl,ns0),vecf(nq1*nq2*nl,ns0)

        double precision :: area, ee1, a0, aM, tw_angle, tw1, factor, prefactor, epsilon
        double precision :: mass, m00, me, hbar, rj
         double precision phi,V1,V2,W1,W2

         integer,parameter :: Nb=3
         double precision :: bandk(1:Nb,1:Nx*Ny)
         double complex :: eigveck(1:Nb,1:Nb,1:Nx*Ny)
         double complex :: rho(1:Nb,1:Nb,1:Nx*Ny)

  


      End module global




