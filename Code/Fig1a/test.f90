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

            !do kx=0,Nx-1
            !do ky=0,Ny-1

            ! write(ci,'(I1.1)') Kx
            ! write(cj,'(I1.1)') Ky

            ! open(11,file='N12/TMD_band-1_spstate_Nx3_Ny4_Kx'//ci//'_Ky'//cj//'.txt',status='old')
            ! do int_l=1,mb
            !  read(11,*) real,imag !int_l,labelbasis(int_l)
            !  wf(int_l,kx*Ny+ky)=real+ii*imag
            ! end do
            ! close(11)

            !end do
            !end do

            print *, 'Allocate MM matrix successfully.'


            V1=17.50d0
            W1=-6.50d0
            V2=-0.0d0 !-10.5 !-11.00d0
            W2=11.75 !12.0d0
            phi= -1.020843d0 !-0.97  !-56.49d0/180.0d0*pi !-91.d0/180.0d0*pi !-0.9859365d0

            tw_angle=1.6d0/180.0d0*pi

            do while(tw_angle<2.4/180.d0*pi)

            do while(V2>=-16.0d0)


            !call make_band()
            !print*, 'Single-particle Band structure.'

            call berryphase()

            V2=V2-1.d0
            end do

            tw_angle=tw_angle+0.5d0/180.0*pi
            end do
            !call make_harteefock()
            stop

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






