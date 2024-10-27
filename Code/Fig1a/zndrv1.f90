

       subroutine zndrv1(labelbasis, two)
        use global
        implicit none 

        
        integer(kind=8) :: two(0:Nx*Ny-1)
        integer(kind=8) :: labelbasis(1:MAXbasis)
        !integer :: countbasis


        !double complex :: Vn(0:Nx*Ny-1,0:Nx*Ny-1,0:Nx*Ny-1,0:Nx*Ny-1)

        character*2 ci

        double complex :: checkz,check0(1:dimen),check(1:dimen)

        !****************************************************

        integer :: maxn != dimen !MAXbasis !n is the Maximum dimension of the matrix
        integer,parameter :: maxnev = 10 !nev is the number of eigenvalues requested
        integer,parameter :: maxncv = 48 !
        integer :: ldv != maxn    !generally, ldv=n

        double complex,allocatable :: resid(:) !resid(1:maxn)

        double complex,allocatable :: v(:,:) !v(1:ldv,1:maxncv)

        integer :: iparam(1:11)
        integer :: ipntr(1:14)

        double complex,allocatable :: workd(:) !workd(1:3*maxn)
        double complex :: workev(1:3*maxncv)
        double complex :: workl(1:3*maxncv*maxncv+5*maxncv)

        double precision :: rwork(1:maxncv)


        character*1 bmat
        character*2  which
        integer :: ido,info,nev,ncv,nd,lworkl
        double precision :: tol

        !***********************************************************

        logical :: rvec

        integer :: ldz,nconv
        logical :: select0(1:maxncv)

        double complex :: sigma

        double complex :: d(1:maxncv)
        double complex,allocatable :: z(:,:) !z(1:maxn,1:maxnev)

        character*1 :: howmny

        integer :: dotimes,int_i,int_j,int_n


        !***********************************************************




        nd=dimen
        maxn=dimen
        ldv=dimen
        allocate(z(1:maxn,1:maxnev))
        allocate(workd(1:3*maxn))
        allocate(v(1:ldv,1:maxncv))
        allocate(resid(1:maxn))

        nev=10
        ncv=24
        if(nd>maxn) print *, 'Nd is greater than MAXN.'
        if(nev>maxnev) print *, 'nev is great than MAXNEV.'
        if(ncv>maxncv) print *, 'ncv is great than MAXNCV.'

        bmat = "I"   !   local scalars
        which = "SR" !   Smallest real part

        lworkl=3*ncv*ncv+5*ncv
        tol = 1.0D-6 !0.0      !determines the stopping criterion. If tol=0, machine precision is used
        info = 0       !info=0 means a random starting vector is requested
        ido = 0        !must initially be set to 0

        iparam(1) = 1
        iparam(3) = 300
        iparam(7) = 1     !mode 1


        rvec=.true.
        howmny= "A"
        ldz=nd


        !call makeVn(Vn)


        call znaupd(ido,bmat,nd,which,nev,tol,resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,rwork,info)

        dotimes=0
        do while(ido==1.or.ido==-1)

         call Mv(workd(ipntr(1)), workd(ipntr(2)) )

         call znaupd(ido,bmat,nd,which,nev,tol,resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,rwork,info)

         dotimes=dotimes+1
         if(mod(dotimes,100)==0) print *, 'ARPACK process =', dotimes
        end do

        if(info<0.or.ido/=99) print *, 'Errors in znaupd, info=', info


         call zneupd(rvec,howmny,select0,d,z,ldz,sigma,workev,bmat,nd,which,&
         nev,tol,resid,ncv,v,ldv,iparam,ipntr,workd,workl,lworkl,rwork,info)


         if(info==0) then

           print *, 'Diagonization is successfully in ARPACK.'

           !write(ci,'(I2)') Total_Kx+Total_Ky*Nx
           open(20,file='energyspec.dat',position='append')

            do int_i=1,nev
              write(20, '(I8,2X,F16.8)') Total_Kx*Ny+Total_Ky, dble(d(int_i))
              print*, Total_Kx*Ny+Total_Ky, d(int_i)
            end do
            write(20,*) ' '
           close(20)

          else
          
             print *, ' Error with _neupd, info = ', info

          end if



          if(info==0.and.rvec==.true.) then

           nconv=iparam(5)
           print*, 'Number of converge vectors =', nconv


           do int_n=1,nconv

            do int_i=1,dimen
             check0(int_i)=z(int_i,int_n)
            end do

            call Mv(check0, check )
            
            checkz=0.0D-0
            do int_i=1,dimen
             checkz=checkz+check(int_i)-d(int_n)*check0(int_i)
            end do

            !print *, '****************************************'
            print *, 'Errors of',int_n,'vector =',abs(checkz)
            !print *, '****************************************'

           end do

           !do int_n=1,3
           ! do int_i=1,dimen
           !  wf(int_i,wfnum)=z(int_i,wfnum)  !Needs mofifications here
           ! end do
           !end do

           write(ci,'(I2.2)') Total_Ky+Total_Kx*Ny
           open(30,file='vector1'//ci//'.dat')
           write(30,*) d(1)
           write(30,*) dimen 
            do int_i=1,dimen
              write(30,*) int_i,z(int_i,1)
            end do
           close(30)

           write(ci,'(I2.2)') Total_Ky+Total_Kx*Ny
           open(31,file='vector2'//ci//'.dat')
           write(31,*) d(2)
           write(31,*) dimen
            do int_i=1,dimen
              write(31,*) int_i,z(int_i,2)
            end do
           close(31)

           write(ci,'(I2.2)') Total_Ky+Total_Kx*Ny
           open(32,file='vector3'//ci//'.dat')
           write(32,*) d(3)
           write(32,*) dimen
            do int_i=1,dimen
              write(32,*) int_i,z(int_i,3)
            end do
           close(32)


          end if 




          return
         end subroutine zndrv1




         subroutine Mv(vec1, vec2)
          use global
          implicit none

          
          !integer(kind=8) :: two(0:Nx*Ny-1)
          !integer(kind=8) :: labelbasis(1:MAXbasis)
          !integer :: countbasis

          double complex :: element
                  
          integer :: exchange
          integer :: j1,j2,j3,j4
          integer :: countnonzero

          !double complex :: Vn(0:Nx*Ny-1,0:Nx*Ny-1,0:Nx*Ny-1,0:Nx*Ny-1)

          !integer :: tempMM(0:MAXnonzero),tempposi(0:MAXnonzero)

         
          integer :: int_i,int_j,int_n
          double complex :: vec1(1:dimen),vec2(1:dimen)

          integer :: Tid,Omp_Get_Thread_Num,NTHREADS,OMP_GET_NUM_THREADS
          integer,parameter :: CHUNK = 100

          !*********************************************************            


          !$Omp parallel shared(vec1,vec2)&
          !$Omp private(Tid,int_n,int_i,int_j)
           Tid= OMP_GET_THREAD_NUM()
          !$Omp do schedule(dynamic, chunk)

            do int_i=1,dimen
             vec2(int_i)=0.0

             int_n=Mi(int_i)
             do while(int_n<Mi(int_i+1))
               vec2(int_i)=vec2(int_i)+  MM(int_n)*vec1(Mj(int_n))
               int_n=int_n+1
             end do
            end do

          !$Omp end parallel  


          !!$Omp parallel shared(vec1,vec2,labelbasis,two)&
          !!$Omp private(Tid,int_i,int_j,int_h,exchange,element)

          ! Tid= OMP_GET_THREAD_NUM()

          !!$Omp do schedule(dynamic, chunk)



            !do int_i=1,dimen

             !vec2(int_i)=0.0d0
             

             !call makenonzero(labelbasis,two,int_i,tempMM,tempposi,countnonzero,element)
             
             !vec2(int_i)=element*vec1(int_i)  !diagonal Ek
             !--------------------------------------------------                          
                                            !the other non-zero part is V(j1,j2,j3,j4)
          

             !int_h=1
             !do while(int_h<=countnonzero)
              
              !int_j=tempMM(int_h)
              
              !if(int_j<0) then
              ! exchange=-1
              ! int_j=-int_j
              !else
              ! exchange=1
              !end if
              !!print*, exchange,int_l              

              !j1=mod(int_j,Nx*Ny)
              !int_j=(int_j-j1)/(Nx*Ny)
              !j2=mod(int_j,Nx*Ny)
              !int_j=(int_j-j2)/(Nx*Ny)
              !j3=mod(int_j,Nx*Ny)
              !int_j=(int_j-j3)/(Nx*Ny)
              !j4=mod(int_j,Nx*Ny)
              
               
               !element=exchange*Vn(j1,j2,j3,j4)

               !vec2(int_i)=vec2(int_i)+conjg(element)*vec1(tempposi(int_h))
               
               !int_h=int_h+1
               
             !end do !int_h


            !end do !int_i



           !!$Omp end parallel


           return
          end subroutine Mv





          subroutine mkl_zheev(nd,Ma,eigvalue)
           !use global
           implicit none



         !*********************************************************

            integer :: nd

            double complex :: Ma(1:nd,1:nd)

            double precision :: eigvalue(1:nd)



            double precision :: w(1:nd)

            integer,parameter :: lwork=4000000 !3*nd+1

            integer,parameter :: lrwork=4000000 !3*nd+1

            integer,parameter :: liwork=500000 !3*nd+1

            double complex :: work(1:lwork)

            double precision :: rwork(1:lrwork)

            integer :: iwork(1:liwork)

            integer :: info

          !*******************************************************
            integer :: int_l
            character*2 ci


          !if(dimen/=0) then
          call ZHEEVD('V','U',nd,Ma,nd,W,WORK,LWORK,RWORK,LRWORK,IWORK,LIWORK,INFO)
          !call zheev('V','U',nd,Ma,nd,w,work,lwork,rwork,info)
          !end if


          if(info==0) then

           !write(ci,'(I2)') Total_Kx+Total_Ky*Nx
           !open(20,file='eigen'//ci//'.dat')
           !open(20,file='eigen.dat')
           ! do int_l=1,dimen
           !   write(20,*) Nx+Total_Ky*Nx, w(int_l)
           ! end do
           !close(20)

           eigvalue(:)=w(:)
          end if



           return
          end subroutine mkl_zheev



       double precision function Converttheta(key)
        use global
        implicit none

        double complex :: key
        double precision :: xx,yy

        xx=dble(key)
        yy=aimag(key)

       if(abs(xx-0.d0)<EX) then
         if(yy>0.0d0) then
           Converttheta=PI/2.d0
         else
           Converttheta=-PI/2.d0      
          !if(yy<0.0d0) then
          ! Converttheta=-PI/2.d0
          !else
          ! Converttheta=0.d0
          !end if
         end if
       else

        Converttheta=datan(yy/xx)
        if(xx<0.0d0) then
         if(yy>=0.0d0) then
          Converttheta=Converttheta+PI
         else
          Converttheta=Converttheta-PI
         end if
        end if

       end if

       !if(Converttheta<0.0d0) Converttheta=Converttheta+2*pi

      return
      END function Converttheta


