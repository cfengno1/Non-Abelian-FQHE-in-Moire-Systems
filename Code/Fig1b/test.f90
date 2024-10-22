      !*************************************************************
      !The program for eigensystem of Flat band model on checkboard
      !lattice. The results are the same with Phys. Rev. X 1, 021014
      !(2011). 
      !
      !*************************************************************



      Program main
       use global
       implicit none

         !integer :: int_i,int_j,int_k,int_l,int_m,int_n,&
         !           int_p,int_q,int_r,int_s,int_t,int_u

         !integer,allocatable :: basis(:,:)
         integer :: countbasis
         integer(kind=8) :: two(0:Nx*Ny-1)
         integer(kind=8),allocatable :: labelbasis(:)
 
         !double complex,allocatable :: MM(:)
         !integer,allocatable :: MM(:)
         !integer(kind=8),allocatable :: Mi(:)
         !integer,allocatable :: Mj(:)

         !integer :: countnum

         !double complex :: Vn(0:Nx*Ny-1,0:Nx*Ny-1,0:Nx*Ny-1,0:Nx*Ny-1)
         !*********************************************************
         !character *2 ci

         double precision :: real,imag !,qq,tmp
         character*1 ci,cj
         integer :: kx,ky
         integer :: int_l
         !********************************************************
         double precision :: starttime,finishtime
         !********************************************************

          call cpu_time(starttime)

          call random_seed()

          !Allocate(basis(0:MAXbasis,0:Ne-1))
           allocate(labelbasis(1:MAXbasis))
 
            Allocate(MM(1:MAXmatrix))
            !MM(:)=0.0D-0
            !MM(:)=0
            Allocate(Mi(1:MAXbasis))
            Allocate(Mj(1:MAXmatrix))
            !Mi(:)=0
            !Mj(:)=0

            allocate(wf(1:mb,0:Nx*Ny-1))
            allocate(Ek(1:mb,0:Nx*Ny-1))


            print *, 'Allocate MM matrix successfully.'


            call make_band()
            print*, 'Single-particle Band structure.'
            stop

            !call make_harteefock()
            !stop

            !!call makeVn_TMD(Vn)
            !call make_Vn(Vn)
            print*, 'Vn elements.'
            !stop


           !Total_Kx=0
           !Total_Ky=0
           !do Total_Kx=0,Nx-1
           !do Total_Ky=Total_Kx,Ny-1

           ! call constructbasis(labelbasis,countbasis,two)

           !**********************************************************************

            !!call makeVn_TMD(Vn)
           ! call makehamilton(labelbasis,countbasis,two, Vn)
            
           !**********************************************************************
           ! print *, 'Number of non-zero elemets of Hamilton:', countnum
           !*************************************************************Diagonization     
            !stop

           ! call zndrv1(labelbasis, two)          

           !***************************************************************
           !end do
           !end do

            call cpu_time(finishtime)

            write(*,*) 'Total running time is', finishtime-starttime
     


        End program main






