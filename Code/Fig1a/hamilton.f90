          subroutine makehamilton(labelbasis,countbasis,two, Vn)
           use global
           implicit none


          integer(kind=8) :: two(0:Nx*Ny-1)
          integer(kind=8) :: labelbasis(1:MAXbasis)
          integer :: countbasis

          double complex :: Vn(0:Nx*Ny-1,0:Nx*Ny-1,0:Nx*Ny-1,0:Nx*Ny-1)


          integer :: int_l,int_k,int_m,int_n,int_r,int_s,int_t,int_u

          integer(kind=8) :: int_o

          integer :: j1,j2,j3,j4

          integer :: px1,px2,px3,px4,py1,py2,py3,py4

          integer :: kix,kiy


          integer :: tempbasis(0:Ne-1),step0(0:Ne-1),step1(0:Ne-1),&
                     step2(0:Ne-1),step3(0:Ne-1),step4(0:Ne-1),temp

          integer :: exchange,jump1,jump2,jump3,jump4

          integer :: flag1,flag2,flag3,flag4

          !double complex :: func

          double complex :: element,tempz,diagEk


         !*********************************************************

          !double complex :: tempHam0(0:MAXnonzero)
          integer :: tempposi(0:MAXnonzero)
          double complex :: tempMM(0:MAXnonzero)

         !********************************************************


           ! do kx=0,3-1
           ! do ky=0,4-1

           !  write(ci,'(I1.1)') Kx
           !  write(cj,'(I1.1)') Ky

           !  open(11,file='N12/TMD_band-1_spstate_Nx3_Ny4_Kx'//ci//'_Ky'//cj//'.txt',status='old')
           !  do int_l=1,mb
           !   read(11,*) real,imag !int_l,labelbasis(int_l)
           !   wf(int_l,kx*Ny+ky)=real+ii*imag
           !  end do
           !  close(11)

           ! end do
           ! end do



            countnum=0
            int_l=0
            do int_l=1,dimen  !int_l related to the basis(int_l)=|Phi>

             !flag3=0     
             int_u=-1
             !tempHam0(:)=0
             tempposi(:)=-1
             tempMM(:)=0.0d0

             int_m=0
             int_o=labelbasis(int_l)
             do int_n=0,Nx*Ny-1
              if(int_o>=two(int_n)) then
               step0(int_m)=int_n
               int_m=int_m+1
               int_o=int_o-two(int_n)
              end if
             end do
             !print*, 'int_l=',int_l

             !!do int_m=0,Ne-1
             !! step0(int_m)=basis(int_l,int_m)
             !!end do

             int_u=int_u+1
             tempposi(int_u)=int_l    !dynamic Ek term    

             diagEk=0.0d0
             do int_n=0,Ne-1
              !kix=int(step0(int_n)/Ny)
              !kiy=mod(step0(int_n),Ny)
              !int_k= kiy+1+kix*ky0
              int_k=step0(int_n)+1
              !print*, int_k
              diagEk=diagEk+ bandk(2,int_k) !valf(int_k,2) + val_hf(int_k)    !func(kix,kiy,0)
             end do
             tempMM(int_u)=diagEk




             do int_s=0,Ne-1    !int_s related to the int_s th electron in basis(int_l)

              j4=step0(int_s)
              px4=int(j4/Ny)
              py4=mod(j4,Ny)


              do int_t=0,Ne-1 !int_s+1,Ne-1

               j3=step0(int_t)
               px3=int(j3/Ny)
               py3=mod(j3,Ny)


               do j2=0,Nx*Ny-1 !j2 from 1 to Nx*Ny-1

                 flag2=1;int_m=0
                 do while(flag2==1.and.int_m<=Ne-1)
                  if(int_m/=int_s.and.int_m/=int_t.and.j2==step0(int_m)) flag2=0
                  int_m=int_m+1
                 end do  !j2 cannot be the existed electrons

                 px2=int(j2/Ny)
                 py2=mod(j2,Ny)


                 do j1=0,Nx*Ny-1

                  flag1=1;int_m=0
                  do while(flag1==1.and.int_m<=Ne-1)
                    if(int_m/=int_s.and.int_m/=int_t.and.j1==step0(int_m)) flag1=0
                    int_m=int_m+1
                  end do

                  px1=int(j1/Ny)
                  py1=mod(j1,Ny)


                  if(j4/=j3.and.j2/=j1.and.&
                     flag1==1.and.flag2==1) then !.and.j1<j3) then

                  !print*, j1,j2,j3,j4
                  !call Vn_TMD(j1,j2,j3,j4,element)
                  !print*, element
                  element=Vn(j1,j2,j3,j4)

                  !if(px1+px2==px3+px4.and.py1+py2==py3+py4) then
                  if(abs(element)>0.000001d0) then

                    do int_m=0,Ne-1
                     tempbasis(int_m)=step0(int_m)
                    end do
                    tempbasis(int_s)=j1
                    tempbasis(int_t)=j2


                     exchange=1
                     do int_m=Ne-2,0,-1  !make the sequence 
                      do int_n=0,int_m  ! 
                       if(tempbasis(int_n)>tempbasis(int_n+1)) then
                        temp=tempbasis(int_n)
                        tempbasis(int_n)=tempbasis(int_n+1)
                        tempbasis(int_n+1)=temp
                        exchange=exchange*(-1)
                       end if
                      end do
                     end do

                     !***************************************

                     int_o=0
                     do int_r=0,Ne-1
                      int_o=int_o+two(tempbasis(int_r))
                     end do

                     call search(int_o,labelbasis,dimen,flag4,int_m)


                     if(flag4==1.and.int_m<=dimen) then

                       int_u=int_u+1
                       if(int_u>MAXnonzero) &
                       print*, 'MAXnonzero elements overflow.'
                       tempposi(int_u)=int_m
                       tempMM(int_u)=exchange*element !(j1+j2*(Nx*Ny)+j3*(Nx*Ny)**2+j4*(Nx*Ny)**3)


                     end if !flag4


                    end if 
                    end if !flag2&flag1&p1+p3-p2-p4=0  

                   end do !j1 

                  end do !j2 


               end do !int_t

              end do !int_s 


              !!countnum=int_u
              !do int_r=0,int_u                      

              !  countnum=countnum+1
              !  MM(countnum)=conjg(tempMM(int_r))
              !  !Mi(countnum)=int_l+1
              !  Mj(countnum)=tempposi(int_r)+1
              ! if(countnum>MAXmatrix) then
              !   print *, 'non-zero elements overflow.'
              !   stop
              ! end if

              !end do


              do int_r=0,int_u

               countnum=countnum+1
               if(countnum>MAXmatrix) then
                 print *, 'MAXmatrix overflow.'
                 stop
               end if

               MM(countnum)=tempMM(int_r)
               Mj(countnum)=tempposi(int_r)

              end do

              Mi(int_l+1)=countnum+1

              if(mod(int_l,10000)==0) print*,'Hamiltonian contructed :', int_l,countnum



            end do !int_l  !choose one basis

            Mi(1)=1

            !print *, 'Hamiltonian has been contructed.'

            !do int_l=1,dimen
            ! print*, int_l,Mi(int_l),Mj(Mi(int_l))
            !end do

           return
           end subroutine makehamilton











          subroutine makenonzero(labelbasis,two,int_l,tempMM,tempposi,countnum0,diagEk)
           use global
           implicit none
          

          integer(kind=8) :: two(0:Nx*Ny-1)
          integer(kind=8) :: labelbasis(1:MAXbasis)
          !integer :: countbasis

          integer :: countnum0
          double complex :: diagEk

          integer :: int_l

          integer :: int_k,int_m,int_n,int_r,int_s,int_t,int_u

          integer(kind=8) :: int_o

          integer :: j1,j2,j3,j4
                     
          integer :: px1,px2,px3,px4,py1,py2,py3,py4
          
          integer :: kix,kiy


          integer :: tempbasis(0:Ne-1),step0(0:Ne-1),step1(0:Ne-1),&
                     step2(0:Ne-1),step3(0:Ne-1),step4(0:Ne-1),temp

          integer :: exchange,jump1,jump2,jump3,jump4

          integer :: flag1,flag2,flag3,flag4

          double complex :: func

          !double complex :: element,tempz

          character *2 ci

         !*********************************************************

          !double complex :: tempHam0(0:MAXnonzero)
          integer :: tempposi(0:MAXnonzero)
          integer :: tempMM(0:MAXnonzero)
          !integer :: tempsign(0:MAXnonzero)

          !double complex :: MM(1:MAXmatrix)
          !integer :: MM(1:MAXmatrix)
          !integer(kind=8) :: Mi(1:MAXbasis)
          !integer :: Mj(1:MAXmatrix)

          !integer :: countnum
         !*********************************************************




            !countnum=0
            !int_l=0
            !do int_l=0,dimen-1  !int_l related to the basis(int_l)=|Phi>

             !flag3=0     
             int_u=-1
             !tempHam0(:)=0
             tempposi(:)=-1
             tempMM(:)=-1

             int_m=0
             int_o=labelbasis(int_l)
             do int_n=0,Nx*Ny-1
              if(int_o>=two(int_n)) then
               step0(int_m)=int_n
               int_m=int_m+1
               int_o=int_o-two(int_n)
              end if
             end do

             !!do int_m=0,Ne-1
             !! step0(int_m)=basis(int_l,int_m)
             !!end do
                                    
             int_u=int_u+1
             tempposi(int_u)=int_l    !dynamic Ek term    
             
             diagEk=0.0d0     
             !do int_n=0,Ne-1
             ! kix=int(step0(int_n)/Ny)
             ! kiy=mod(step0(int_n),Ny)
             ! diagEk=diagEk+func(kix,kiy,0)
             !end do


               
             do int_s=0,Ne-1    !int_s related to the int_s th electron in basis(int_l)
             
              j3=step0(int_s)
              px3=int(j3/Ny)
              py3=mod(j3,Ny)
             
              
              do int_t=0,Ne-1 !int_s+1,Ne-1

               j4=step0(int_t)
               px4=int(j4/Ny)
               py4=mod(j4,Ny)        

               
               do j2=0,Nx*Ny-1 !j2 from 1 to Nx*Ny-1


                 flag2=1;int_m=0
                 do while(flag2==1.and.int_m<=Ne-1)
                  if(int_m/=int_s.and.int_m/=int_t.and.j2==step0(int_m)) flag2=0
                  int_m=int_m+1
                 end do  !j2 cannot be the existed electrons

                 px2=int(j2/Ny)
                 py2=mod(j2,Ny)


                  !px1=mod(px2+px4-px3+Nx,Nx)
                  !py1=mod(py2+py4-py3+Ny,Ny)
                  !j1=px1*Ny+py1

                 do j1=0,Nx*Ny-1

                  flag1=1;int_m=0
                  do while(flag1==1.and.int_m<=Ne-1)
                    if(int_m/=int_s.and.int_m/=int_t.and.j1==step0(int_m)) flag1=0
                    int_m=int_m+1
                  end do

                  px1=int(j1/Ny)
                  py1=mod(j1,Ny)



                  if(j4/=j3.and.j2/=j1.and.&
                     flag1==1.and.flag2==1) then !.and.j1<j3) then

                  !if(px1+px2==px3+px4.and.py1+py2==py3+py4) then


                    do int_m=0,Ne-1
                     tempbasis(int_m)=step0(int_m)
                    end do
                    tempbasis(int_s)=j1
                    tempbasis(int_t)=j2


                     exchange=1
                     do int_m=Ne-2,0,-1  !make the sequence 
                      do int_n=0,int_m  ! 
                       if(tempbasis(int_n)>tempbasis(int_n+1)) then
                        temp=tempbasis(int_n)
                        tempbasis(int_n)=tempbasis(int_n+1)
                        tempbasis(int_n+1)=temp
                        exchange=exchange*(-1)
                       end if
                      end do
                     end do

                     !*****************************************
                                                         !calculate the phase
                     
                     !do int_n=0,Ne-1
                     ! step3(int_n)=step0(int_n)
                     !end do
                     
                     !jump4=1
                     !int_r=0
                     !do while(step0(int_s)/=step3(int_r))
                     ! int_r=int_r+1
                     ! jump4=jump4*(-1)
                     !end do

                     !step3(int_s)=Nx*Ny

                     !int_r=0
                     !jump2=1
                     !do while(step0(int_t)/=step3(int_r))
                     ! if(step3(int_r)/=Nx*Ny) jump2=jump2*(-1)
                     ! int_r=int_r+1
                     !end do

                     !step3(int_t)=Nx*Ny

                     !int_n=2
                     !do int_r=0,Ne-1
                     ! if(step3(int_r)/=Nx*Ny) then
                     !  step4(int_n)=step3(int_r)
                     !  int_n=int_n+1
                     ! end if
                     !end do
                     !step4(1)=j3
                     !step4(0)=j1

                     !int_r=0                 
                     !jump3=1
                     !do int_r=Ne-2,0,-1  !make the sequence 
                     ! do int_n=0,int_r  ! 
                     !  if(step4(int_n)>step4(int_n+1)) then
                     !   temp=step4(int_n)
                     !   step4(int_n)=step4(int_n+1)
                     !   step4(int_n+1)=temp
                     !   jump3=jump3*(-1)
                     !  end if
                     ! end do
                     !end do


                     !***************************************
                           

                     int_o=0
                     do int_r=0,Ne-1
                      int_o=int_o+two(tempbasis(int_r))
                     end do
                     
                     call search(int_o,labelbasis,dimen,flag4,int_m)
          
                  

         !if(int_i==5.and.int_l==10.and.int_j==1) &
         !write(*,'(A2,I3,I3,I3,I3,I3,I3,I3,I3,I3,I3,A2,I4,I4,I4,I4,I16)') &
         !        '(',tempbasis(0),tempbasis(1),tempbasis(2),tempbasis(3),&
         !            tempbasis(4),tempbasis(5),tempbasis(6),tempbasis(7),&
         !            tempbasis(8),tempbasis(9),')',j1,j3,&
         !            exchange,jump2*jump3*jump4,&
         !            int_m


   
     
                     if(flag4==1.and.int_m<=dimen) then
                                                  
                        
                       !element=exchange* & !jump2*jump3*jump4*&
                       !(   (Conjg(func(px1,py1,1))*func(px2,py2,1)*&
                       !     Conjg(func(px3,py3,2))*func(px4,py4,2))*&
                       !     4*cos((px4-px3)*2.*Pi/(2.*Nx))*&
                       !       cos((py4-py3)*2.*Pi/(2.*Ny)) &

                       !    +(Conjg(func(px3,py3,1))*func(px4,py4,1)*&
                       !     Conjg(func(px1,py1,2))*func(px2,py2,2))*&
                       !     4*cos((px2-px1)*2.*Pi/(2.*Nx))*&
                       !       cos((py2-py1)*2.*Pi/(2.*Ny)) &

                       !    -(Conjg(func(px1,py1,1))*func(px4,py4,1)*&
                       !      Conjg(func(px3,py3,2))*func(px2,py2,2))*&
                       !     4*cos((px2-px3)*2.*Pi/(2.*Nx))*&
                       !       cos((py2-py3)*2.*Pi/(2.*Ny)) +&

                       !    -(Conjg(func(px3,py3,1))*func(px2,py2,1)*&
                       !     Conjg(func(px1,py1,2))*func(px4,py4,2))*&
                       !     4*cos((px4-px1)*2.*Pi/(2.*Nx))*&
                       !       cos((py4-py1)*2.*Pi/(2.*Ny)) & 
                       !)*Hubb_U/(1.0*Nx*Ny)

                       !element=exchange* & !jump2*jump3*jump4*&
                       !(   (Conjg(func(px1,py1,1))*func(px2,py2,1)*&
                       !     Conjg(func(px3,py3,2))*func(px4,py4,2))*&
                       !     (1.+exp(ii*(py4-py3)*2.*pi/(1.0*Ny)))*&
                       !     (1.+exp(ii*(px3-px4)*2.*pi/(1.0*Nx)))   &

                       !    +(Conjg(func(px3,py3,1))*func(px4,py4,1)*&
                       !     Conjg(func(px1,py1,2))*func(px2,py2,2))*&
                       !     (1.+exp(ii*(py2-py1)*2.*pi/(1.0*Ny)))*&
                       !     (1.+exp(ii*(px1-px2)*2.*pi/(1.0*Nx)))   &

                       !    -(Conjg(func(px1,py1,1))*func(px4,py4,1)*&
                       !      Conjg(func(px3,py3,2))*func(px2,py2,2))*&
                       !     (1.+exp(ii*(py2-py3)*2.*pi/(1.0*Ny)))*&
                       !     (1.+exp(ii*(px3-px2)*2.*pi/(1.0*Nx)))   &

                       !    -(Conjg(func(px3,py3,1))*func(px2,py2,1)*&
                       !     Conjg(func(px1,py1,2))*func(px4,py4,2))*&
                       !     (1.+exp(ii*(py4-py1)*2.*pi/(1.0*Ny)))*&
                       !     (1.+exp(ii*(px1-px4)*2.*pi/(1.0*Nx)))   & 
                       !)*Hubb_U/(1.0*Nx*Ny)



                      !if(int_m==int_l.and.flag3==0) then
                      ! do int_n=0,Ne-1
                      !  kix=int(step0(int_n)/Ny)
                      !  kiy=mod(step0(int_n),Ny)
                      !  element=element+func(kix,kiy,0)
                      ! end do
                      ! flag3=1
                      !end if


                      !int_r=0
                      !int_t=0
                      !do while(int_t==0.and.int_r<=int_u)
                      ! if(tempposi(int_r)==int_m) then
                      !  int_t=1
                      ! else
                      !  int_r=int_r+1
                      ! end if
                      !end do

                      !if(int_t==1) then
                      !  tempHam0(int_r)=tempHam0(int_r)+element
                      !else if(int_t==0.and.int_r>int_u) then
                      !  int_u=int_u+1
                      !  tempHam0(int_u)=element
                      !  tempposi(int_u)=int_m
                      !  if(int_u>MAXnonzero) & 
                      !    print*, 'non-zero elements overflow.'
                      !end if
                      
                      
                       int_u=int_u+1
                       if(int_u>MAXnonzero) & 
                       print*, 'MAXnonzero elements overflow.'
                       tempposi(int_u)=int_m
                       tempMM(int_u)=exchange*(j1+j2*(Nx*Ny)+j3*(Nx*Ny)**2+j4*(Nx*Ny)**3)
                       
 
                     end if !flag4


                    !end if 
                    end if !flag2&flag1&p1+p3-p2-p4=0  

                   end do !j1 

                  end do !j2 
                           
                
               end do !int_t

              end do !int_s 

              countnum0=int_u       
              


               !do int_r=int_u-1,0,-1  !make the sequence
               ! do int_n=0,int_r  !
               !  if(tempposi(int_n)>tempposi(int_n+1)) then
               !   int_s=tempposi(int_n)
               !   tempposi(int_n)=tempposi(int_n+1)
               !   tempposi(int_n+1)=int_s
               !   tempz=tempHam0(int_n)
               !   tempHam0(int_n)=tempHam0(int_n+1)
               !   tempHam0(int_n+1)=tempz
               !  end if
               ! end do
               !end do



              !do int_r=0,int_u                      

              !  countnum=countnum+1
              !  MM(countnum)=conjg(tempHam0(int_r))
              !  !Mi(countnum)=int_l+1
              !  Mj(countnum)=tempposi(int_r)+1
              ! if(countnum>MAXmatrix) then
              !   print *, 'non-zero elements overflow.'
              !   stop
              ! end if

              !end do
              

              !do int_r=0,int_u
              
              ! countnum=countnum+1
              ! if(countnum>MAXmatrix) then
              !   print *, 'MAXmatrix overflow.'
              !   stop
              ! end if
               
              ! MM(countnum)=tempMM(int_r)
              ! Mj(countnum)=tempposi(int_r)+1
                             
              !end do

              !Mi(int_l+2)=countnum+1

              !if(mod(int_l,10000)==0) print*,'Hamiltonian contructed :', int_l,countnum

                         

            !end do !int_l  !choose one basis

            !Mi(1)=1

            !print *, 'Hamiltonian has been contructed.'
            
            !do int_l=1,dimen
            ! print*, int_l,Mi(int_l),Mj(Mi(int_l))
            !end do
            
           return
           end subroutine makenonzero
           
           
           
           
           
           subroutine search(key,labelbasis,countbasis,flag4,location)
            use global
            implicit none

             integer :: countbasis,flag4,location
             integer(kind=8) :: labelbasis(1:MAXbasis)


             integer(kind=8) :: key     ! 所要寻找的值
             integer :: up      ! 记录每一个小组的类型起始位置
             integer :: down    ! 记录每一个小组的类型结束位置
             integer :: middle  ! 记录每一个小组的类型中间位置

                                ! 一开始的小组范围就是整个类型
             down=1
             up=countbasis
             middle=int((down+up)/2)

                                ! 如果key值超出范围, 铁定不存在类型中
             if((key>labelbasis(down)).or.(key<labelbasis(up))) then
               flag4=0
               return
             end if


             do while(down<=up)
              if(key>labelbasis(middle)) then
                                ! 如果key>中间值，那数据就落在上半部
                up=middle-1
                middle=int((down+up)/2)
              else if(key<labelbasis(middle)) then
                                ! 如果 key < 中间值，那数据就落在下半部
                down=middle+1
                middle=int((down+up)/2)
              else if(key==labelbasis(middle)) then
                location=middle
                flag4=1
                return
              end if
             end do

             flag4=0
            return
           end subroutine search





           subroutine makeVn_TMD(Vn)
            use global
            implicit none


            double complex :: element
            integer :: j1,j2,j3,j4

            integer :: px1,px2,px3,px4,py1,py2,py3,py4

            integer :: gg1,gg2,gg3,gg4,m1,m2,m3,m4,n1,n2,n3,n4

            double complex :: Vn(0:Nx*Ny-1,0:Nx*Ny-1,0:Nx*Ny-1,0:Nx*Ny-1)

            double precision :: qq
            double complex :: tmp

            character*1 ci,cj

            integer :: kx,ky,int_l,int_m,int_n

            !integer,parameter :: aM=1.0d0 
            aM=1.0d0

            !do kx=0,3-1
            !do ky=0,4-1

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



          Vn=0.0d0

          do j1=0,Nx*Ny-1

           px1=int(j1/Ny)
           py1=mod(j1,Ny)


           do j2=0,Nx*Ny-1

            px2=int(j2/Ny)
            py2=mod(j2,Ny)

            do j3=0,Nx*Ny-1

             px3=int(j3/Ny)
             py3=mod(j3,Ny)

             do j4=0,Nx*Ny-1

              px4=int(j4/Ny)
              py4=mod(j4,Ny)

              if(mod(px1+px2-px3-px4+7*Nx,Nx)==0.and.mod(py1+py2-py3-py4+7*Ny,Ny)==0) then
              !---------------------------------
              tmp=0.0d0

              gg1=25-1; gg2=25-1; gg3=25-1; gg4=25-1
              do gg1=0,49-1
              m1= int(gg1/7); n1= mod(gg1,7)
              if((m1-3)*(m1-3)+(n1-3)*(n1-3)<12) then
              do gg2=0,49-1
              m2= int(gg2/7); n2= mod(gg2,7)
              if((m2-3)*(m2-3)+(n2-3)*(n2-3)<12)then
              do gg3=0,49-1
              m3= int(gg3/7); n3= mod(gg3,7)
              if((m3-3)*(m3-3)+(n3-3)*(n3-3)<12)then
              do gg4=0,49-1
              m4= int(gg4/7); n4= mod(gg4,7)
              if((m4-3)*(m4-3)+(n4-3)*(n4-3)<12) then

                !m1= int(g1/7); n1= mod(g1,7)
                !m2= int(g2/7); n2= mod(g2,7)
                !m3= int(g3/7); n3= mod(g3,7)
                !m4= int(g4/7); n4= mod(g4,7)
                !if(abs((m1-3)*(m1-3)+(n1-3)*(n1-3))<1.and.&
                !   abs((m2-3)*(m2-3)+(n2-3)*(n2-3))<1.and.&
                !   abs((m3-3)*(m3-3)+(n3-3)*(n3-3))<1.and.&
                !   abs((m4-3)*(m4-3)+(n4-3)*(n4-3))<1    ) then


                kx=px1-px4  !+ (m1-m4) 
                ky=py1-py4  !+ (n1-n4)              

                if(px1+px2+(m1+m2)*Nx==px3+px4+(m3+m4)*Nx.and.&
                   py1+py2+(n1+n2)*Ny==py3+py4+(n3+n4)*Ny ) then

                qq=sqrt( (1.0*kx/Nx +1.0*(m1-m4))**2 + (1.0*ky/Ny +1.0*(n1-n4))**2 )
                if(abs(qq)>0.0001d0) then
                   !qq= sqrt(1.0*(kx*kx+ky*ky))      

                   tmp=tmp+2*pi/(Nx*Ny*aM*qq)*&
                               (conjg(wf(gg1*2+1,j1)*wf(gg2*2+1,j2))*wf(gg3*2+1,j3)*wf(gg4*2+1,j4)+&
                                conjg(wf(gg1*2+2,j1)*wf(gg2*2+2,j2))*wf(gg3*2+2,j3)*wf(gg4*2+2,j4)+&
                                conjg(wf(gg1*2+1,j1)*wf(gg2*2+2,j2))*wf(gg3*2+2,j3)*wf(gg4*2+1,j4)+&
                                conjg(wf(gg1*2+2,j1)*wf(gg2*2+1,j2))*wf(gg3*2+1,j3)*wf(gg4*2+2,j4)  )

                end if
                end if

               end if
              end do
               end if
              end do
               end if
              end do
               end if
              end do


              Vn(j1,j2,j3,j4)= tmp
              !element=tmp
              if(abs(tmp)>0.001) &
                write(1100,'(I4,I4,I4,I4,4X,F12.6,F12.6)') j1,j2,j3,j4, dble(tmp),aimag(tmp)
              end if

             end do

            end do

           end do

          end do


          return
          end subroutine makeVn_TMD






           subroutine makeVncb(Vn)
            use global
            implicit none
   
            
            integer :: j1,j2,j3,j4

            integer :: px1,px2,px3,px4,py1,py2,py3,py4

            double complex :: func

            double complex :: Vn(0:Nx*Ny-1,0:Nx*Ny-1,0:Nx*Ny-1,0:Nx*Ny-1)

            double precision,parameter :: Hubb_U=1.0d0


          do j1=0,Nx*Ny-2

           px1=int(j1/Ny)
           py1=mod(j1,Ny)

          
           do j2=0,Nx*Ny-2

            px2=int(j2/Ny)
            py2=mod(j2,Ny)

            do j3=j1+1,Nx*Ny-1

             px3=int(j3/Ny)
             py3=mod(j3,Ny)

             do j4=j2+1,Nx*Ny-1

              px4=int(j4/Ny)
              py4=mod(j4,Ny)


              Vn(j1,j2,j3,j4)= & 
                       (   (Conjg(func(px1,py1,1))*func(px2,py2,1)*&
                            Conjg(func(px3,py3,2))*func(px4,py4,2))*&
                            (1.+exp(ii*(py4-py3)*2.*pi/(1.0*Ny)))*&
                            (1.+exp(ii*(px3-px4)*2.*pi/(1.0*Nx)))   &

                           +(Conjg(func(px3,py3,1))*func(px4,py4,1)*&
                            Conjg(func(px1,py1,2))*func(px2,py2,2))*&
                            (1.+exp(ii*(py2-py1)*2.*pi/(1.0*Ny)))*&
                            (1.+exp(ii*(px1-px2)*2.*pi/(1.0*Nx)))   &

                           -(Conjg(func(px1,py1,1))*func(px4,py4,1)*&
                             Conjg(func(px3,py3,2))*func(px2,py2,2))*&
                            (1.+exp(ii*(py2-py3)*2.*pi/(1.0*Ny)))*&
                            (1.+exp(ii*(px3-px2)*2.*pi/(1.0*Nx)))   &

                           -(Conjg(func(px3,py3,1))*func(px2,py2,1)*&
                            Conjg(func(px1,py1,2))*func(px4,py4,2))*&
                            (1.+exp(ii*(py4-py1)*2.*pi/(1.0*Ny)))*&
                            (1.+exp(ii*(px1-px4)*2.*pi/(1.0*Nx)))   & 
                       )*Hubb_U/(1.0*Nx*Ny)


             end do

            end do

           end do

          end do


          return
          end subroutine makeVncb



           double complex Function func(int_kx,int_ky,int_x)
            use global
            implicit none

            integer :: int_kx,int_ky,int_x

            double complex :: h12,h21,Apart1,Bpart1,h11,h22

            double precision :: kx,ky

      double precision,parameter :: t1=1.0d0
      double precision,parameter :: t2=1.0d0/(2.+sqrt(2.0))
      double precision,parameter :: t3=1.0d0/(2.+2.*sqrt(2.0))
      double precision,parameter :: xphi=Pi/4.0
      !double precision,parameter :: M=0.0d0
      !double precision,parameter :: Hubb_U=1.0d0


             kx=int_kx*2.*Pi/(1.*Nx)
             ky=int_ky*2.*Pi/(1.*Ny)


             h11=4*t3*cos(kx)*cos(ky)+2*t2*(cos(kx)-cos(ky))
             h22=4*t3*cos(kx)*cos(ky)-2*t2*(cos(kx)-cos(ky))
             !h12=4*t1*cos(phi)*cos(kx/2.)*cos(ky/2.)-ii*4*t1*sin(phi)*sin(kx/2.)*sin(ky/2.)
             !h21=4*t1*cos(phi)*cos(kx/2.)*cos(ky/2.)+ii*4*t1*sin(phi)*sin(kx/2.)*sin(ky/2.)
       
             h12=t1*exp(ii*phi)*(1+exp(ii*(ky-kx)))+t1*exp(-ii*xphi)*(exp(ii*ky)+exp(-ii*kx))
             h21=conjg(h12)

             if(int_x==0) then

              func=(h11+h22)/2-sqrt(h12*h21-h11*h22+(h11+h22)*(h11+h22)/4)

             else

      Apart1=-(-h11+h22+sqrt(h11*h11+4*h12*h21-2*h11*h22+h22*h22))/(2*h21)
               Bpart1=1


               if(int_x==1) then
                func=Apart1/sqrt(abs(Apart1*conjg(Apart1)+Bpart1*conjg(Bpart1)))
               else if(int_x==2) then
                func=Bpart1/sqrt(abs(Apart1*conjg(Apart1)+Bpart1*conjg(Bpart1)))
               end if

              end if

            return
           end function func

  
  
  
           subroutine checkhamilton()
            use global
            implicit none
            
            !double complex :: MM(1:MAXmatrix)
            !integer :: Mi(1:MAXbasis),Mj(1:MAXmatrix)
  
  
            integer :: int_i,int_j,int_k,int_l,int_m,int_n,int_r
            
            integer :: flag3
                   
            !integer :: countnum
            
            
         !********************************************************
         integer :: Tid,Omp_Get_Thread_Num,NTHREADS,OMP_GET_NUM_THREADS
         integer,parameter :: CHUNK = 100
            
            
  
            !!$Omp parallel shared(Mi,Mj,MM,countnum) &
            !!$Omp private(Tid,int_n,int_l,int_m,int_k)

            !!$ Tid= OMP_GET_THREAD_NUM()

            int_k=1
            int_l=0
            !!$Omp do schedule(dynamic, 1000)

              !int_k=1
              do int_n=1,countnum,1
               if(int_n==Mi(int_k)) then
                int_l=int_l+1
                int_k=int_k+1
               end if

               int_m=Mj(int_n)
         
               if(int_l==int_m.and.abs(aimag(MM(int_n)))>EX) then
                write(*,*) int_l,int_m,MM(int_n),'non-Hermin'
               end if
           
              end do

            !!$Omp end parallel



            !$Omp parallel shared(Mi,Mj,MM,countnum)&
            !$Omp private(Tid,int_n,int_l,int_m,int_k,int_r,flag3)

            !$ Tid= OMP_GET_THREAD_NUM()

             int_r=1
             int_l=1
            !$Omp do schedule(dynamic, 100000)

             !int_r=1
             do int_n=10,countnum,1000
               if(int_n==Mi(int_r)) then
                int_l=int_l+1
                int_r=int_r+1
               end if

              int_m=Mj(int_n)

              flag3=0;int_k=1
              do while(flag3==0.and.int_k<=countnum)
               if(int_l==Mj(int_k).and.int_m==Mi(int_k)) then
                if(abs(aimag(MM(int_n)+MM(int_k)))>EX.or.&
                   abs(dble(MM(int_n)-MM(int_k)))>EX) &
                   write(*,*) int_l,int_m,MM(int_n),'non-Hermin'
                flag3=1
               end if
               int_k=int_k+1
              end do   
             
             end do

            !$Omp end parallel


            print *, 'Hermin has been checked.'
            
           return
           end subroutine checkhamilton






        subroutine make_band(k1,k2,vec)
        use global
        implicit none


        !integer  kp1, kp2

        integer,parameter :: nsp=nq1*nq2*nl, nsk=nk1*nk2
        integer, parameter :: nnp=nsp, np1=nq1, np2=nq2


        integer, parameter :: neibt=12
        integer :: neib(1:nq1*nq2,12)

        double precision :: VD1
        double precision :: angle, x,y, g, prefac
        double complex :: tensor(0:30,2,2), ck, ck0, ckshift

        !double precision phi,V1,V2,W1,W2

        real*8  valp(nnp)
        complex*16 hamp(nnp,nnp),vecp(nnp,nnp) !,vecf(nnp,nsk)
        !complex*16 vecf0(nnp,nsk),vecf3(nnp,nsk)

        integer  kspin, ks1

        double complex ci, k_plus,k_minus, kg1,kg2
        integer m1(2), m2(2)

        double precision :: mom_line(0:100,2)


        integer :: i,j,j1,j2,k1,k2,ki,ik

        integer :: ix,iy,ix0,iy0,ix1,iy1,ix2,iy2

        integer :: jj,j11,j12,j13,jy,jy1

        integer :: k11,k12,q1,q2,q11,q12,iq1,iq2,ia1,ib1

        character *2  chx,chy

        double complex :: curvature,tmp1,tmp2,tmp3,tmp4
        double precision :: C1,C2,C3
        double precision, external :: Converttheta

        double complex :: vec(1:nnp,1:Nb)

        !-----------------------------------------
        me= 9.1093837e-31
        m00 = 0.620d0*me
        hbar=1.054571817e-34
        rJ=6.582119569509066e-13/hbar

        prefactor = hbar**2/(2.0d0*m00)*rJ*1.0e20
        ci=cmplx(0.0d0, 1.0d0)
        VD1= 0.0d0

        a0=3.52d0
        factor=prefactor
        !!theta=pi/3.0d0
        !!antip=0
        !tw_angle=1.8d0/180.0d0*pi !1.800d0/180.0d0*pi
        tw1=tw_angle*180.0d0/pi
        aM=a0/(2.0d0*sin(tw_angle/2.0d0))

        do j=1, 7
        angle=(j-1)*pi/3.0d0
        prefac=4*pi/(dsqrt(3.0d0)*aM)
        glist(j)=prefac*exp(ci*angle)
        enddo



        g1=glist(1)
        g2=glist(3)  !! g(2)-g(1)
        m1(1)=1
        m1(2)=0    !!! for g1
        m2(1)=-1
        m2(2)=1

        !mg0=0
        mq0=0
        mq0(1,1)=m1(1)   !1
        mq0(1,2)=m1(2)   !0
        mq0(2,1)=m2(1)   !-1
        mq0(2,2)=m2(2)   ! 1
        !mg0(1,1)=dfloat(m2(2))/(m2(2)*m1(1)-m2(1)*m1(2))
        !mg0(1,2)=dfloat(-m1(2))/(m2(2)*m1(1)-m2(1)*m1(2))
        !mg0(2,1)=dfloat(m2(1))/(m2(1)*m1(2)-m1(1)*m2(2))
        !mg0(2,2)=dfloat(-m1(1))/(m2(1)*m1(2)-m1(1)*m2(2))

        k_minus= (glist(1)+glist(6))/3.0d0
        k_plus=(glist(1)+glist(2))/3.0d0
        kg1=(glist(2)+glist(3))/3.0d0



        neib=0
        do ix=1, np1
        do jy=1, np2
        i= jy+(ix-1)*np2

        ix1=ix+1
        if(ix==np1)ix1=1
        jy1=jy+1
        if(jy==np2)jy1=1
        ix2=ix-1
        if(ix==1)ix2=np1

        if(ix1.ne.0)neib(i,1)=jy+(ix1-1)*np2
        if(jy1.ne.0)neib(i,2)=jy1+(ix-1)*np2
        if(jy1*ix2.ne.0)neib(i,3)=jy1+(ix2-1)*np2
        j11=jy+(ix1-1)*np2
        j12=jy1+(ix-1)*np2
        j13=jy1+(ix2-1)*np2
        neib(j11,neibt/2+1)=i
        neib(j12,neibt/2+2)=i
        neib(j13,neibt/2+3)=i
        enddo
        enddo

        do i=1, np1*np2
                j1=neib(i,1)
                if(j1.ne.0)then
                j2=neib(j1,2)
                neib(i,4)=j2
                neib(j2, 4+neibt/2)=i
                endif

                j1=neib(i,2)
                j2=neib(j1,3)
                        if(j1*j2.ne.0)then
                        neib(i,5)=j2
                        neib(j2, 5+neibt/2)=i
                        endif

                j1=neib(i,3)
                j2=neib(j1,1+neibt/2)
                if(j1*j2.ne.0)then
                neib(i,6)=j2
                neib(j2, 6+neibt/2)=i
                endif
        enddo

        !do i=1,np1*np2
        !do ik=1,6
        !if(abs(neib(i,ik))>1000.or.neib(i,ik)==0) print*,i,ik,neib(i,ik)
        !end do
        !end do
        !print*,neib(6,1)
        !stop

        !V1=17.50d0
        !W1=-6.50d0
        !V2=-10.0d0 !-10.5 !-11.00d0
        !W2=11.75 !12.0d0
        !phi= -0.97  !-56.49d0/180.0d0*pi !-91.d0/180.0d0*pi !-0.9859365d0





        tensor=0.0d0
        tensor(0,1,2)=W1 !+W2
        tensor(0,2,1)=W1 !+W2
        do j1=1,3
        tensor(j1,1,1)=tensor(j1,1,1)-V1*cdexp(ci*(-1)**(j1-1)*phi)
        tensor(j1,2,2)=tensor(j1,2,2)-V1*cdexp(ci*(-1)**(j1-1)*(-phi))
        enddo

        do j1=4,6
        tensor(j1,1,1)=tensor(j1,1,1)-V2 !*cdexp(ci*(-1)**(j1-1)*phi)
        tensor(j1,2,2)=tensor(j1,2,2)-V2  !*cdexp(ci*(-1)**(j1-1)*(-phi))
        enddo

        do j1=2,3
        tensor(j1,2,1)=tensor(j1,2,1)+W1
        enddo
        j1=1
        tensor(j1,2,1)=tensor(j1,2,1)+W2   !! *2
        tensor(j1,1,2)=tensor(j1,1,2)+W2   !! *2
        do j1=5, 5
        tensor(j1,2,1)=tensor(j1,2,1)+W2
        enddo


        !do k1=0,40
        !  mom_line(k1,1)=0.0
        !  mom_line(k1,2)=k1*0.5/40
        !end do
        !do k1=41,65
        !  mom_line(k1,1)=(k1-40)*0.333333/25
        !  mom_line(k1,2)=0.5+(k1-40)*(0.666666-0.5)/25
        !end do
        !do k1=66,100
        !  mom_line(k1,1)=0.333333-(k1-65)*0.333333/35
        !  mom_line(k1,2)=0.666666-(k1-65)*0.666666/35
        !end do



        kspin=1
        ks1=1
        !!if(kspin==2)ks1=-1
        !!if(kspin.eq.1)then
        !k1=1
        !k2=1
        !do k1=1,kx0
        !do k2=1,ky0 !!nk2
        !ki=k2+(k1-1)*ky0

         ck0=(k1-1)/dble(kx0)*g1+(k2-1)/dble(ky0)*g2

        !do ki=100,0,-1

         !ck0=mom_line(ki,1)*g1 + mom_line(ki,2)*g2


         do i=1,nnp
         do j=1,nnp
           hamp(i,j)=cmplx(0.d0,0.d0)
         enddo
         enddo


        ix0=(np1-1)/2+1
        iy0=(np2-1)/2+1 

        ix=ix0
        iy=iy0
        do ix=1, np1
        do iy=1, np2

        i=iy+(ix-1)*np2          !!iy+(ix-1)*ny
        !ix0=(np1-1)/2+1
        !iy0=(np2-1)/2+1

        ck=ck0+glist(1)*(ix-ix0)+glist(2)*(iy-iy0)  !+ckshift

        hamp(2*i-1,2*i-1)=hamp(2*i-1,2*i-1)+factor*abs(ck-ks1*k_plus)**2-VD1/2.0d0
        hamp(2*i,2*i)    =hamp(2*i,2*i)    +factor*abs(ck-ks1*k_minus)**2+VD1/2.0d0
        hamp(2*i-1,2*i)  =hamp(2*i-1,2*i)  +W1  !tensor(0,1,2)
        hamp(2*i,2*i-1)  =hamp(2*i,2*i-1)  +W1  !tensor(0,2,1)

        go to 14
        do ik=1,6
        jj=neib(i,ik)
        !if(abs(jj)>1000) print*,i,ik,jj
        !!if(j.ne.0)then
        !print*, i,ik,jj
        !if(ik==5) print*, ix,iy, int(jj/np2)+1,mod(jj,np2)

        hamp(2*i-1,2*jj-1)=hamp(2*i-1,2*jj-1) +tensor(ik,1,1)
        hamp(2*i,2*jj)=hamp(2*i,2*jj) +tensor(ik,2,2)
        hamp(2*i-1,2*jj)=hamp(2*i-1,2*jj) +tensor(ik,1,2)
        hamp(2*i,2*jj-1)=hamp(2*i,2*jj-1) +tensor(ik,2,1)

        hamp(2*jj-1,2*i-1)=hamp(2*jj-1,2*i-1) +conjg(tensor(ik,1,1))
        hamp(2*jj,2*i)=hamp(2*jj,2*i) +conjg(tensor(ik,2,2))
        hamp(2*jj,2*i-1)=hamp(2*jj,2*i-1) +conjg(tensor(ik,1,2))
        hamp(2*jj-1,2*i)=hamp(2*jj-1,2*i) +conjg(tensor(ik,2,1))
        !!endif
        enddo!ik
        14 continue

        !go to 16
        !-------------------------------------------------------
        !neighbor 1
        ix1=ix+1
        iy1=iy
        if(ix1<=np1) then
        jj=iy1+(ix1-1)*np2

        !print*, ix1,iy1
        hamp(2*i-1,2*jj-1)=hamp(2*i-1,2*jj-1) -V1*cdexp(ii*phi) !
        hamp(2*i,2*jj)=hamp(2*i,2*jj) -V1*cdexp(-ii*phi)

        hamp(2*jj-1,2*i-1)=hamp(2*jj-1,2*i-1) -V1*cdexp(-ii*phi) 
        hamp(2*jj,2*i)=hamp(2*jj,2*i) -V1*cdexp(ii*phi)

        hamp(2*i,2*jj-1)=hamp(2*i,2*jj-1) + W2
        hamp(2*jj-1,2*i)=hamp(2*jj-1,2*i) + W2

        hamp(2*i-1,2*jj)=hamp(2*i-1,2*jj) + W2
        hamp(2*jj,2*i-1)=hamp(2*jj,2*i-1) + W2
        end if

        !neighbor 2
        ix1=ix
        iy1=iy+1
        if(iy1<=np2) then
        jj=iy1+(ix1-1)*np2

        !print*, ix1,iy1
        hamp(2*i-1,2*jj-1)=hamp(2*i-1,2*jj-1) -V1*cdexp(-ii*phi) !-V1*cdexp(ii*phi)
        hamp(2*i,2*jj)=hamp(2*i,2*jj) -V1*cdexp(ii*phi)

        hamp(2*jj-1,2*i-1)=hamp(2*jj-1,2*i-1) -V1*cdexp(ii*phi) !-V1*cdexp(ii*phi)
        hamp(2*jj,2*i)=hamp(2*jj,2*i) -V1*cdexp(-ii*phi)

        hamp(2*i,2*jj-1)=hamp(2*i,2*jj-1) + W1
        hamp(2*jj-1,2*i)=hamp(2*jj-1,2*i) + W1
        end if

        !neighbor 3
        ix1=ix-1
        iy1=iy+1
        if(ix1>=1.and.iy1<=np2) then
        jj=iy1+(ix1-1)*np2

        !print*, ix1,iy1
        hamp(2*i-1,2*jj-1)=hamp(2*i-1,2*jj-1) -V1*cdexp(ii*phi) !-V1*cdexp(-ii*phi)
        hamp(2*i,2*jj)=hamp(2*i,2*jj) -V1*cdexp(-ii*phi)

        hamp(2*jj-1,2*i-1)=hamp(2*jj-1,2*i-1) -V1*cdexp(-ii*phi) !-V1*cdexp(-ii*phi)
        hamp(2*jj,2*i)=hamp(2*jj,2*i) -V1*cdexp(ii*phi)

        hamp(2*i,2*jj-1)=hamp(2*i,2*jj-1) + W1
        hamp(2*jj-1,2*i)=hamp(2*jj-1,2*i) + W1
        end if


        !go to 16
        !next-neighbor 1
        ix1=ix+1
        iy1=iy+1
        if(ix1<=np1.and.iy1<=np2) then
        jj=iy1+(ix1-1)*np2

        !print*, ix1,iy1
        hamp(2*i-1,2*jj-1)=hamp(2*i-1,2*jj-1) -V2  !-V1*cdexp(ii*phi)
        hamp(2*i,2*jj)=hamp(2*i,2*jj) -V2 

        hamp(2*jj-1,2*i-1)=hamp(2*jj-1,2*i-1) -V2  !-V1*cdexp(ii*phi)
        hamp(2*jj,2*i)=hamp(2*jj,2*i) -V2 
        end if

        !next-neighbor 2
        ix1=ix-1
        iy1=iy+2
        if(ix1>=1.and.iy1<=np2) then
        jj=iy1+(ix1-1)*np2

        !print*, ix1,iy1
        hamp(2*i-1,2*jj-1)=hamp(2*i-1,2*jj-1) -V2  !-V1*cdexp(ii*phi)
        hamp(2*i,2*jj)=hamp(2*i,2*jj) -V2 

        hamp(2*jj-1,2*i-1)=hamp(2*jj-1,2*i-1) -V2  !-V1*cdexp(ii*phi)
        hamp(2*jj,2*i)=hamp(2*jj,2*i) -V2 

        hamp(2*i,2*jj-1)=hamp(2*i,2*jj-1) + W2
        hamp(2*jj-1,2*i)=hamp(2*jj-1,2*i) + W2
        end if

        !next-neighbor 3
        ix1=ix-2
        iy1=iy+1
        if(ix1>=1.and.iy1<=np2) then
        jj=iy1+(ix1-1)*np2

        hamp(2*i-1,2*jj-1)=hamp(2*i-1,2*jj-1) -V2  !-V1*cdexp(-ii*phi)
        hamp(2*i,2*jj)=hamp(2*i,2*jj) -V2 

        hamp(2*jj-1,2*i-1)=hamp(2*jj-1,2*i-1) -V2 !-V1*cdexp(-ii*phi)
        hamp(2*jj,2*i)=hamp(2*jj,2*i) -V2         
        end if
        16 continue

        end do!iy
        end do!ix
        !-------------------------------------



        do i=1, nnp-1
        do j=i+1, nnp
        if(hamp(i,j).ne.conjg(hamp(j,i)))then
        write(*,*)hamp(i,j),hamp(j,i), 'non-hermitian', i,j
        stop
        endif
        ! if(abs(hamp(i,j))>0.1d0) write(150,*) i,j,hamp(i,j)
        enddo
        enddo
        !stop
        !call eigen1(nnp, hamp, valp, vecp)
        call mkl_zheev(nnp, hamp, valp)


             ck=ck0
             open(11,file="val.dat", access='append')
             do i=1,10
               write(11,'(I4,F12.6,F12.6,4X,I4, F24.16)') ki,real(ck), aimag(ck), i,valp(i)
             end do
             write(11,*) ' '
             close(11)

         !vecf(1:nnp, ki)= hamp(1:nnp,2)
         !vecf0(1:nnp, ki)= hamp(1:nnp,1)
         !vecf3(1:nnp, ki)= hamp(1:nnp,3)
         !valf(ki,1:10) = valp(1:10)
         vec(:,1:Nb)=hamp(:,1:Nb)

         !write(800,'(I4,2F8.4,2X,10F12.4)') k1,mom_line(k1,1),mom_line(k1,2),valp(1:10)


        !enddo !k2                !!!for maxnum loop
        !enddo !k1
        !enddo

        return
        end subroutine make_band

       


        subroutine berryphase()
        use global 


        integer :: k1,k2
        integer :: ix0,iy0,ix,iy,ix1,iy1
        integer :: j11,j12,j13,jj

        double complex :: curvature,tmp1,tmp2,tmp3,tmp4
        double precision :: C(1:Nb)
        double precision, external :: Converttheta

        integer,parameter :: nnp=nq1*nq2*nl
        integer :: int_n

        double complex :: vec1(1:nnp,1:Nb),vec2(1:nnp,1:Nb),vec3(1:nnp,1:Nb),vec4(1:nnp,1:Nb)







        !-------------------------------------------
        !Berry phase

        

        C(:)=0.0d0
        do k1=1,kx0
        do k2=1,ky0 


          !ki=k2+(k1-1)*ky0

          ix0=k1
          iy0=k2
          jj=iy0+(ix0-1)*ky0
          call make_band(ix0,iy0,vec1)

          ix=ix0+1
          !if(ix0==kx0) ix=1
          iy=iy0
          j11=iy+(ix-1)*ky0
          call make_band(ix,iy,vec2)

          ix1=ix0+1
          !if(ix0==kx0) ix1=1
          iy1=iy0+1
          !if(iy0==ky0) iy1=1
          j12=iy1+(ix1-1)*ky0
          call make_band(ix1,iy1,vec3)

          ix2=ix0
          iy2=iy0+1
          !if(iy0==ky0) iy2=1
          j13=iy2+(ix2-1)*ky0
          call make_band(ix2,iy2,vec4)


          do int_n=1,Nb
          !---------------------------------- 
          tmp1=0.0d0 !     
          do i=1,nnp
            tmp1=tmp1+ conjg(vec1(i,int_n))*vec2(i,int_n)
          end do     
          tmp2=0.0d0 !
          do i=1,nnp
            tmp2=tmp2+ conjg(vec2(i,int_n))*vec3(i,int_n)
          end do          
          tmp3=0.0d0 !
          do i=1,nnp
            tmp3=tmp3+ conjg(vec3(i,int_n))*vec4(i,int_n)
          end do          
          tmp4=0.0d0 !
          do i=1,nnp
            tmp4=tmp4+ conjg(vec4(i,int_n))*vec1(i,int_n)
          end do

          curvature=tmp1*tmp2*tmp3*tmp4

          !print*, jj,j11,j12,j13,abs(curvature),aimag(curvature)
          !write(*,'(I4,I4,I4,I4,2X,F8.4,F8.4)') jj,j11,j12,j13,abs(curvature),Converttheta(curvature)
          if(int_n==1) write(21,'(I4,I4,2X,F12.6,F12.6)') k1,k2,abs(curvature),Converttheta(curvature)
          if(int_n==2) write(22,'(I4,I4,2X,F12.6,F12.6)') k1,k2,abs(curvature),Converttheta(curvature)
          if(int_n==3) write(23,'(I4,I4,2X,F12.6,F12.6)') k1,k2,abs(curvature),Converttheta(curvature)          

          C(int_n)=C(int_n)+aimag(log(curvature)) !Converttheta(curvature)
          !---------------------------------- 
          end do !int_n


        end do
        !write(21,*) ' '
        !write(22,*) ' '
        end do

        print*, C(1)/2/pi, C(2)/2/pi, C(3)/2/pi
        write(50,'(F8.4,F12.6,2X,F12.6,F12.6,F12.6)') tw_angle/pi*180,V2,C(1)/2/pi, C(2)/2/pi, C(3)/2/pi 
        !stop


        return
        end subroutine berryphase












        subroutine make_Vn(Vn)                !  ns ~ ab/2*pi*ell^2
        !use paramsizew
        !use tmd
        use global
        implicit none

        double complex :: Vn(0:Nx*Ny-1,0:Nx*Ny-1,0:Nx*Ny-1,0:Nx*Ny-1)


        double precision :: vq0

        integer :: j1,j2,j3,j4,jj,jx1,jx4,jy1,jy4

        integer :: kx,ky,kx1,ky1

        integer :: ni,nj,n,i,j,j41,j12,iaa,ibb,ij

        double complex :: cqx ,cqx1

        integer :: q1,q2,q3,q4, qd1,qd2



        Vn=0.0d0

        ee1=1.60217663*8.9875517923d0*1e+03
        !dgate=0
        !if(dgate.le.0.1)epsilon=5.0d0
        epsilon=5.0d0

        !a0=3.52d0
        !factor=prefactor
        !!theta=pi/3.0d0
        !!antip=0
        !tw_angle=1.800d0/180.0d0*pi
        !tw1=tw_angle*180.0d0/pi
        !aM=a0/(2.0d0*sin(tw_angle/2.0d0))

        area=1.0d0
        ee1=ee1*2.0d0*pi/area/epsilon
        area=aM*aM*dsqrt(3.0d0)/2.0d0*kx0*ky0
        ee1=ee1/area



        mnt=0
        ind_mnt=0
        do kx=0,kx0-1
        do ky=0,ky0-1
         j=ky+1+kx*ky0

         mnt(j,1)=kx
         mnt(j,2)=ky
         mnt(j+ns0,1)=kx
         mnt(j+ns0,2)=ky

         ind_mnt(kx,ky)=j
        enddo
        enddo

        !vc=0.0d0


        do kx=0,kx0-1
        do ky=0,ky0-1

        j1=ind_mnt(kx,ky)

        do kx1=0,kx0-1
        do ky1=0,ky0-1
        j2=ind_mnt(kx1,ky1)

        jx1=mod(kx-kx1+kx0,kx0)
        jy1=mod(ky-ky1+ky0,ky0)

        mnt_d(j1,j2)=ind_mnt(jx1,jy1)
        enddo
        enddo
        enddo
        enddo



        !vsq=0
        !vsq1=0.0d0
        !do j3=1,ns0
        !do j2=1,ns0
        !do j1=1,ns0
        !jx4=mod(mnt(j1,1)+mnt(j2,1)-mnt(j3,1)+kx0,kx0)
        !jy4=mod(mnt(j1,2)+mnt(j2,2)-mnt(j3,2)+ky0,ky0)
        !j41=ind_mnt(jx4, jy4)
        !ni=j1+(j2-1)*ns0!!!!!+(j3-1)*ns*ns
        !nj=j3+(j41-1)*ns0!!!!!+(j3-1)*ns*ns
        !!!vcoe(ni,nj)=1
        !!!!!!jq1=mnt_d(j41,j1)
        !enddo
        !enddo
        !enddo
        !v=dcmplx(0.d0,0.d0)
        !v10=dcmplx(0.0d0,0.0d0)
        !v00=dcmplx(0.0d0,0.0d0)
        !u=dcmplx(0.0d0, 1.0d0)

        do j1=1,ns0
        do j2=1,ns0
        j12=j1+(j2-1)*ns0
        do j3=1,ns0
        !do j41=1,1!!!ns0
        !!!!!!if(vcoe(j12,j34).eq.1)then
        jx4=mod(mnt(j1,1)+mnt(j2,1)-mnt(j3,1)+kx0,kx0)
        jy4=mod(mnt(j1,2)+mnt(j2,2)-mnt(j3,2)+ky0,ky0)
        j4=ind_mnt(jx4, jy4)
!       jq=mnt_d(j4,j1) !!! momentum transfer?
        !!!!do j4=1,ns
!!!! this part,  jq has to be determined from
        jj=j1+(j2-1)*ns0+(j3-1)*ns0**2
        !!jj=jj+(j4-1+ji4*ns0)*ns*ns*ns
        n=jj

        vq0=0.0d0
        do q1=-nq1,nq1
        do q2=-nq2,nq2
        qd1=mnt(j1,1)+mnt(j2,1)-mnt(j3,1)-mnt(j4,1)
        qd2=mnt(j1,2)+mnt(j2,2)-mnt(j3,2)-mnt(j4,2)
        qd1=qd1/kx0
        qd2=qd2/ky0
        q3=-q1-qd1*mq0(1,1)-qd2*mq0(2,1)
        q4=-q2-qd1*mq0(1,2)-qd2*mq0(2,2)

        cqx=(mnt(j1,1)-mnt(j4,1))/dble(kx0)*g1+q1*glist(1)+&
            (mnt(j1,2)-mnt(j4,2))/dble(ky0)*g2+q2*glist(2)


        vq0=abs(cqx)
        if(vq0.ge.1e-15.and.abs(q3)<=nq1.and.abs(q4)<=nq2)then
        !if(dgate.gt.0.1)then
        ! vq0=tanh(dgate*vq0)/vq0
        !else
         vq0=1.0d0/vq0
        !endif
        !v(1,n)=v(1,n)+vq0*fq(j1,j4,q1,q2)*fq(j2,j3,q3,q4)
!!!!!        Vn(j1-1,j2-1,j3-1,j4-1)=Vn(j1-1,j2-1,j3-1,j4-1)+vq0*fq11(j1,j4,q1,q2)*fq11(j2,j3,q3,q4)


        !iaa=mnt(j1,1)-mnt(j4,1) !!!!+ns0+(q1+nq1)*ns0*2
        !ibb=mnt(j1,2)-mnt(j4,2) !!!!+ns0+(q1+nq1)*ns0*2
        !iaa=iaa+(q1*mg0(1,1)+q2*mg0(2,1))*kx0
        !!!iqb=q2-nb1b
        !ibb=ibb+(q1*mg0(1,2)+q2*mg0(2,2))*ky0
        !cqx1=iaa/dble(kx0)*g1+ibb/dble(ky0)*g2

        !if(abs(iaa).le.ms0.and.abs(ibb).le.ms0)then
        !vsq(n,iaa,ibb)=vsq(n,iaa,ibb)+fq(j1,j4,q1,q2)*fq(j2,j3,q3,q4)
        !vsq1(n,iaa,ibb)=vsq1(n,iaa,ibb)+vq0*fq(j1,j4,q1,q2)*fq(j2,j3,q3,q4)
        !endif
        endif

        !enddo
        enddo
        enddo
        enddo
        enddo
        enddo
        !v=v*ee1/2.0d0
        !v10=v10*ee1/2.0d0
        !v00=v00*ee1/2.0d0
        !vsq1=vsq1*ee1/2.0d0

        Vn=Vn*ee1/2.0d0


        return
        !Hartee-Fock  
        val_hf=0.0d0
        do j1=1,ns0
        do j2=1,ns0

        do ij=1,2
        if(ij==1)then
        j3=j2
        j4=j1
        else
        j3=j1
        j4=j2
        endif

        vq0=0.0d0
        do q1=-nq1,nq1
        do q2=-nq2,nq2
        q3=-q1
        q4=-q2
        cqx=(mnt(j1,1)-mnt(j4,1))/dble(kx0)*g1+q1*glist(1)
        cqx=cqx+(mnt(j1,2)-mnt(j4,2))/dble(ky0)*g2
        cqx=cqx+q2*glist(2)
        vq0=cdabs(cqx)


        if(vq0.ge.1e-15)then
        !if(dgate.gt.0.1)then
        ! vq0=tanh(dgate*vq0)/vq0
        !else
         vq0=1.0d0/vq0
        !endif
        if(ij==1)then
!!!!!        val_hf(j2)=val_hf(j2)+2.0d0*vq0*fq0(j1,j4,q1,q2)*fq(j2,j3,q3,q4)
        !val_hf(j2)=val_hf(j2)+2.0d0*vq0*fq00(j1,j4,q1,q2)*fq(j2,j3,q3,q4)
        else
!!!!!        val_hf(j2)=val_hf(j2)-2.0d0*vq0*abs(fq10(j1,j4,q1,q2)**2)
        endif
        endif
        enddo
        enddo
        enddo
!!!switching
        enddo
        enddo
        val_hf=val_hf*ee1/2.0d0

        !do j1=1,ns0
        !write(23,*)j1, val_hf(j1)
        !enddo
        val_hf=0.0d0

        return
        end subroutine make_Vn
