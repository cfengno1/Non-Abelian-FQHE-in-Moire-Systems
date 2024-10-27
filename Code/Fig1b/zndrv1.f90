



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


