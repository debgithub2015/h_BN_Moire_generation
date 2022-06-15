!! This program calculates the order parameter in h-BN moire structure
!!updates
!! Debajit Chakraborty : 02/16/2022

!! Declaring the global variable
!        MODULE kinds
!
!        IMPLICIT NONE
!
!        integer, parameter :: sp = selected_real_kind(6, 37)
!        integer, parameter :: dp = selected_real_kind(15, 307)
!        integer, parameter :: qp = selected_real_kind(33, 4931)
!
!        END MODULE kinds      


        MODULE global_variable


        IMPLICIT NONE
        INTEGER*4:: natoms_init, natoms_fin, timestep_init,&
                    timestep_fin,x_cell, y_cell,tot_hex,tot_sort_l2,&
                    unit_output_1,unit_output_2,unit_output_3,&
                    unit_output_4,layer1_output,layer2_output,&
                    rcom_output
                    
        REAL*8 ::xlo_init,xhi_init,ylo_init,yhi_init,zlo_init,&
                 zhi_init, xlo_fin,xhi_fin,ylo_fin,yhi_fin,zlo_fin,&
                 zhi_fin,box_orig(3),box_init(3),box_fin(3),rBN_l1,&
                 rBN_l2,lay1_fin_inn1(3),lay1_fin_inn2(3),&
                 lay1_fin_out1(3),lay1_fin_out2(3)
        REAL*8, DIMENSION(:), ALLOCATABLE :: rx_init, ry_init,rz_init, &
                                           rx_fin, ry_fin, rz_fin, &
                                           rx_sorted_l1,ry_sorted_l1, &
                                           rz_sorted_l1,rx_sorted_l2, &
                                           ry_sorted_l2,rz_sorted_l2
        INTEGER*4, DIMENSION(:), ALLOCATABLE ::atom_index_init, &
                                             atom_index_fin,&
                                             layer1_init,layer1_fin, &
                                             layer2_init,layer2_fin, &
                                             sorted_index_l1,&
                                             sorted_index_l2
        INTEGER*4,DIMENSION(:,:),ALLOCATABLE::nn4_sorted_l1,&
                                            nn4_sorted_l2, hexa_l1,&
                                            cnt_atom_mb_l1_init,&
                                            cnt_atom_mb_l1_fin,&
                                            cnt_atom_mb_l2_init,&
                                            cnt_atom_mb_l2_fin,&
                                            sum_atom_mb_l1_init,&
                                            sum_atom_mb_l1_fin,&
                                            sum_atom_mb_l2_init,&
                                            sum_atom_mb_l2_fin
  
        INTEGER*4,DIMENSION(:,:,:),ALLOCATABLE::layer1_block_init,&
                                              layer1_block_fin,&
                                              layer2_block_init,&
                                              layer2_block_fin
        CHARACTER, DIMENSION(:), ALLOCATABLE ::b_init,b_fin,&
                                               b_sorted_l1, b_sorted_l2

        CHARACTER, DIMENSION(:,:,:),ALLOCATABLE:: layer1_block_init_b,&
                                                  layer1_block_fin_b,&
                                                  layer2_block_init_b,&
                                                  layer2_block_fin_b
        REAL*8,DIMENSION(:,:),ALLOCATABLE :: hexa_rcom_l1,hexa_rcom_l2
        REAL*8,DIMENSION(:,:),ALLOCATABLE :: hexa_rgeocp_l1,hexa_rgeocp_l2
        REAL*8,DIMENSION(:,:,:),ALLOCATABLE :: hexa_norm_l1
        REAL*8,DIMENSION(:,:,:,:),ALLOCATABLE:: layer1_block_init_r,&
                                              layer1_block_fin_r, &
                                              layer2_block_init_r,&
                                              layer2_block_fin_r

        LOGICAL,DIMENSION(:,:),ALLOCATABLE:: mesh_block_l1_init,&
                                             mesh_block_l1_fin,&
                                             mesh_block_l2_init,&
                                             mesh_block_l2_fin    

        CHARACTER(len=1000)::order_param_file,hexagon_data_dist_file, &
                            hexagon_data_up_xyz_file,file_layer1,&
                            file_layer2,hexagon_data_down_xyz_file,&
                            hexagon_rcom_data_file

        END MODULE global_variable

        MODULE my_functions

        IMPLICIT NONE
        
        CONTAINS
 
        REAL*8 FUNCTION dist_adjusted(rx2,rx1,ry2,ry1,rz2,rz1,xhi2,&
                         xlo2,ylo2,yhi2,zlo2,zhi2)
!        USE kinds
        USE global_variable

        IMPLICIT NONE
        REAL*8,INTENT(IN) :: rx2,rx1,ry2,ry1,rz2,rz1,xhi2,xlo2,&
                           ylo2,yhi2,zlo2,zhi2
        REAL*8::dx,dy,dz,hx,hy,hz
       
        hx=xhi2-xlo2
        hy=yhi2-ylo2
        hz=zhi2-zlo2

        dx=rx2-rx1
        dy=ry2-ry1
        dz=rz2-rz1

           IF  (dx.GE.hx/2.d0)      dx = dx - hx
           IF  (dx.LT.-hx/2.d0)     dx = dx + hx
           IF  (dy.GE.hy/2.d0)      dy = dy - hy
           IF  (dy.LT.-hy/2.d0)     dy = dy + hy
           IF  (dz.GE.hz/2.d0)      dz = dz - hz
           IF  (dz.LT.-hz/2.d0)     dz = dz + hz

        dist_adjusted = SQRT(dx*dx+dy*dy+dz*dz)
       
        RETURN
        END FUNCTION dist_adjusted

        REAL*8 FUNCTION dist(rx2,rx1,ry2,ry1,rz2,rz1,xhi2,xlo2,&
                           ylo2,yhi2,zlo2,zhi2)
        
        USE global_variable

        IMPLICIT NONE
        REAL*8,INTENT(IN) :: rx2,rx1,ry2,ry1,rz2,rz1,xhi2,xlo2,&
                           ylo2,yhi2,zlo2,zhi2
        REAL*8::dx,dy,dz,hx,hy,hz
       
        hx=xhi2-xlo2
        hy=yhi2-ylo2
        hz=zhi2-zlo2

        dx=rx2-rx1
        dy=ry2-ry1
        dz=rz2-rz1

        dist = SQRT(dx*dx+dy*dy+dz*dz)
        
        RETURN
        END FUNCTION dist
        
        FUNCTION rgeocp(xhex,yhex,zhex)

        USE global_variable

        IMPLICIT NONE
        
        REAL*8,DIMENSION(6),INTENT(IN)::xhex,yhex,zhex
        
        REAL*8,DIMENSION(3):: rgeocp

        INTEGER*4:: i

        REAL*8:: sum_rx,sum_ry,sum_rz

        sum_rx=0.d0
        sum_ry=0.d0
        sum_rz=0.d0

        DO i=1,6

          sum_rx=sum_rx+xhex(i)
          sum_ry=sum_ry+yhex(i)
          sum_rz=sum_rz+zhex(i) 
        
          !WRITE(6,*) i,xhex(i),yhex(i),zhex(i),bhex(i),&
          !           sum_rx,sum_ry,sum_rz,mass(i) 
        
        END DO

        rgeocp(1)=sum_rx/6
        rgeocp(2)=sum_ry/6
        rgeocp(3)=sum_rz/6

        RETURN
        END FUNCTION rgeocp

        FUNCTION rcom(xhex,yhex,zhex,bhex)

        USE global_variable

        IMPLICIT NONE
        
        REAL*8,DIMENSION(6),INTENT(IN)::xhex,yhex,zhex

        CHARACTER,DIMENSION(6),INTENT(IN)::bhex  
        
        REAL*8,DIMENSION(3):: rcom

        INTEGER*4:: i

        REAL*8:: sum_mass,sum_rx,sum_ry,sum_rz, mass(6)

        sum_mass=0.d0
        sum_rx=0.d0
        sum_ry=0.d0
        sum_rz=0.d0

        DO i=1,6

          IF(bhex(i)=='B') mass(i)= 10.811d0
          IF(bhex(i)=='N') mass(i)= 14.0067d0

          sum_mass=sum_mass+mass(i)
          sum_rx=sum_rx+xhex(i)*mass(i)
          sum_ry=sum_ry+yhex(i)*mass(i)
          sum_rz=sum_rz+zhex(i)*mass(i)       
        
          !WRITE(6,*) i,xhex(i),yhex(i),zhex(i),bhex(i),&
          !           sum_rx,sum_ry,sum_rz,mass(i) 
        
        END DO

        rcom(1)=sum_rx/sum_mass
        rcom(2)=sum_ry/sum_mass
        rcom(3)=sum_rz/sum_mass

        RETURN
        END FUNCTION rcom

        FUNCTION cross(x1, x2, x3, y1, y2, y3, z1, z2, z3)
        IMPLICIT NONE
         
        REAL*8, INTENT(IN) :: x1, x2, x3, y1, y2, y3, z1, z2, z3
        REAL*8, DIMENSION(3) :: a, b, cross
      

        a(1) = x2 - x1
        a(2) = y2 - y1
        a(3) = z2 - z1
     
        b(1) = x3 - x1
        b(2) = y3 - y1
        b(3) = z3 - z1

        cross(1) = a(2) * b(3) - b(2) * a(3)
        cross(2) = a(3) * b(1) - a(1) * b(3)
        cross(3) = a(1) * b(2) - a(2) * b(1)
     
        RETURN
        END FUNCTION cross

        FUNCTION ROT_TRANSFORM(nx,ny,nz)

! This transform the unit normal vector to unit z-axis in
! unitary sphere: Rot_y(-beta)*Rot_z*(alpha)<nx,ny,nz>=<0,0,1>
        IMPLICIT NONE
        
        REAL*8, INTENT(IN) :: nx,ny,nz
        REAL*8, DIMENSION(3,3) :: ROT_TRANSFORM
         
        REAL*8 :: a1, b1

        a1=SQRT(nx**2+ny**2+nz**2)
        b1=SQRT(nx**2+ny**2)

        ROT_TRANSFORM(1,1)=nz*nx/(a1*b1)
        ROT_TRANSFORM(1,2)=nz*ny/(a1*b1)
        ROT_TRANSFORM(1,3)=-b1/a1
        ROT_TRANSFORM(2,1)=-ny/b1
        ROT_TRANSFORM(2,2)=nx/b1
        ROT_TRANSFORM(2,3)=0.0d0
        ROT_TRANSFORM(3,1)=nx/a1
        ROT_TRANSFORM(3,2)=ny/a1
        ROT_TRANSFORM(3,3)=nz/a1
      
        END FUNCTION ROT_TRANSFORM

        FUNCTION ROT_TRANSFORM_INV(nx,ny,nz)

! This transform back the unit normal vector at unit z-axis to the
! original unit normal
        IMPLICIT NONE
        
        REAL*8, INTENT(IN) :: nx,ny,nz
        REAL*8, DIMENSION(3,3) :: ROT_TRANSFORM_INV
         
        REAL*8 :: a1, b1

        a1=SQRT(nx**2+ny**2+nz**2)
        b1=SQRT(nx**2+ny**2)

        ROT_TRANSFORM_INV(1,1)=nz*nx/(a1*b1)
        ROT_TRANSFORM_INV(1,2)=-ny/b1
        ROT_TRANSFORM_INV(1,3)=nx/a1
        ROT_TRANSFORM_INV(2,1)=ny*nz/(a1*b1)
        ROT_TRANSFORM_INV(2,2)=nx/b1
        ROT_TRANSFORM_INV(2,3)=ny/a1
        ROT_TRANSFORM_INV(3,1)=-(b1/a1)
        ROT_TRANSFORM_INV(3,2)=0.d0
        ROT_TRANSFORM_INV(3,3)=nz/a1
      
        END FUNCTION ROT_TRANSFORM_INV

        FUNCTION ROT_MAT_Z(alpha) 
        IMPLICIT NONE
        
        REAL*8, INTENT(IN) :: alpha
        REAL*8, DIMENSION(3,3) :: ROT_MAT_Z

        REAL :: costhet, sinthet

        costhet=cos(alpha)
        sinthet=sin(alpha)

        ROT_MAT_Z(1,1)=costhet
        ROT_MAT_Z(1,2)=-sinthet
        ROT_MAT_Z(1,3)=0.d0
        ROT_MAT_Z(2,1)=sinthet
        ROT_MAT_Z(2,2)=costhet
        ROT_MAT_Z(2,3)=0.d0
        ROT_MAT_Z(3,1)=0.d0
        ROT_MAT_Z(3,2)=0.d0
        ROT_MAT_Z(3,3)=1.d0
        
        RETURN
        END FUNCTION ROT_MAT_Z

        REAL*8 FUNCTION theta(u,v)
 
        IMPLICIT NONE

        REAL*8,INTENT(IN)::u(3),v(3)

        REAL*8::udotu,udotv,vdotv,pi

        pi=4.d0*datan(1.d0)

        udotu=DOT_PRODUCT(u,u)
        vdotv=DOT_PRODUCT(v,v)
        udotv=DOT_PRODUCT(u,v)

        theta=(180.d0/pi)*dacos((udotv) &
                   /(SQRT(udotu)*SQRT(vdotv)))

        END FUNCTION theta

        END MODULE my_functions
        
        module DynamicalArrays

          contains

          subroutine AddToList_real(list, element)

          IMPLICIT NONE

          integer*4 :: i, isize
          real*8, intent(in) :: element
          real*8, dimension(:), allocatable, intent(inout) ::list
          real*8, dimension(:), allocatable :: clist


          if(allocated(list)) then
              isize = size(list)
              allocate(clist(isize+1))
              do i=1,isize          
              clist(i) = list(i)
              end do
              clist(isize+1) = element

              deallocate(list)
              call move_alloc(clist, list)

          else
              allocate(list(1))
              list(1) = element
          end if


         end subroutine AddToList_real
          
          subroutine AddToList_char(list, element)

          IMPLICIT NONE

          integer*4 :: i, isize
          character, intent(in) :: element
          character, dimension(:), allocatable, intent(inout) ::list
          character, dimension(:), allocatable :: clist


          if(allocated(list)) then
              isize = size(list)
              allocate(clist(isize+1))
              do i=1,isize          
              clist(i) = list(i)
              end do
              clist(isize+1) = element

              deallocate(list)
              call move_alloc(clist, list)

          else
              allocate(list(1))
              list(1) = element
          end if


         end subroutine AddToList_char

         subroutine AddToList_int(list, element)

          IMPLICIT NONE

          integer*4 :: i, isize
          integer*4, intent(in) :: element
          integer*4, dimension(:), allocatable, intent(inout) ::list
          integer*4, dimension(:), allocatable :: clist


          if(allocated(list)) then
              isize = size(list)
              allocate(clist(isize+1))
              do i=1,isize          
              clist(i) = list(i)
              end do
              clist(isize+1) = element

              deallocate(list)
              call move_alloc(clist, list)

          else
              allocate(list(1))
              list(1) = element
          end if


      end subroutine AddToList_int

      end module DynamicalArrays

        PROGRAM main
        USE global_variable
        USE my_functions        
        USE DynamicalArrays

        IMPLICIT NONE
        
        INTEGER cnt0,cnt1
        INTEGER cnt_rate,cnt_max,cnt_diff

        REAL*8 :: start,finish,wall_time

        CALL system_clock(COUNT_RATE=cnt_rate,COUNT_MAX=cnt_max)

        CALL system_clock(COUNT=cnt0)

        CALL cpu_time(start)

! Reading the structures from the MD dump file
        
        CALL ReadStructure_xyz

! Separate out the first layers
  
        CALL Find_layer_l1

! Separate out the second layers

        CALL Find_layer_l2

!! Find the nearest neighbors in the sorted layer 1

        CALL nearest_neighbor_l1

! Find the nearest neighbors in the sorted layer 2

        CALL nearest_neighbor_l2

!! Find the hexagons in the layer 1

        CALL find_hexagon_lower_layer

! Find the order parameter

        CALL order_parameter

        DEALLOCATE(rx_init, ry_init, rz_init, rx_fin, ry_fin, rz_fin, &
                atom_index_init,atom_index_fin, layer1_init,layer1_fin,&
                layer2_init,layer2_fin,nn4_sorted_l1,nn4_sorted_l2,&
                hexa_l1,b_fin,hexa_rcom_l1,hexa_rcom_l2,hexa_rgeocp_l1,&
                hexa_rgeocp_l2,hexa_norm_l1) 
        
        CALL cpu_time(finish)
    
        WRITE(6,*) '("Time taken for all finished = ",f6.3,&
                                    " seconds.")',finish-start

        CALL system_clock(COUNT=cnt1)
        cnt_diff=cnt1-cnt0
        IF(cnt1<cnt0) cnt_diff=cnt_diff+cnt_max

        wall_time=REAL(cnt_diff)/cnt_rate

        WRITE(6,*) '("Walltime = ",f6.3,"seconds")',wall_time

        END PROGRAM

        SUBROUTINE order_parameter

        USE global_variable
        USE my_functions

        IMPLICIT NONE

        INTEGER*4::i,j,ii,jj,k,cnt1,cnt2,cnt3,cnt4,I_dim,J_dim,item,&
                   cnt5,cnt6,cnt7,cnt8,cnt9,cnt10,cnt11,l2_atom_num,&
                   i1,i2,i3,i4,i5,i6,k1,k2

        REAL*8::xmin,xmax,ymin,ymax,a0,xmin_adj,xmax_adj,ymin_adj,&
              ymax_adj,x_adjust,y_adjust,r_hexa_l2(3),pi, &
              sum_z_l1,sum_z_l2,ave_z_l1,ave_z_l2,cone_height,&
              r_hexa_l1(3),norm_hexa_l1(3),rcom_hexa_l1(3),&
              rcom_hexa_l2(3),rgeocp_hexa_l2(3),rgeocp_hexa_l1(3),&
              angle_two_lay,dist_l2,ord_pm_AA,a0_limit,a0_norm,&
              ord_pm_AA_1,ord_pm_AB,ord_pm_AB_1,ord_pm_AB_op2,&
              ord_pm_AB_inter,ord_pm_AB_no_up_atom,l1_sq,l2_sq,&
              cone_radius,temp_x,temp_y,temp_z,d_rcom,length_scale
        
        LOGICAL::latoms,lsameatom,langlezero,ldisthalf,l1(6),l2(6),&
                 l3(6),l4(6),lhexagon=.TRUE.
 
        CHARACTER::b_hexa_l1,b_hexa_l2,l2_atom_b

        INTEGER*4,DIMENSION(:,:),ALLOCATABLE::item_upper_match

        REAL*8,DIMENSION(:),ALLOCATABLE::color_index,cosphi,sinphi,phi

        REAL*8,DIMENSION(:,:,:),ALLOCATABLE::r_upper_match
                                                
        REAL*8,DIMENSION(:,:),ALLOCATABLE::theta_up_down, &
                                           dist_upper_match

        CHARACTER(len=20),DIMENSION(:),ALLOCATABLE::hexa_op
        CHARACTER,DIMENSION(:,:),ALLOCATABLE::b_upper_match
              
        LOGICAL,DIMENSION(:,:),ALLOCATABLE::mesh_atoms,no_atom_found
        LOGICAL,DIMENSION(:,:,:),ALLOCATABLE::hexa_atom_match

        REAL*8 :: start, finish

        CALL cpu_time(start)

! Finding the cone height (i.e. average distance between the two layers)
        
         sum_z_l1=0.d0
         DO i=1,natoms_fin/2
               sum_z_l1=sum_z_l1+rz_fin(layer1_fin(i))
         END DO
         ave_z_l1=sum_z_l1/FLOAT(natoms_fin/2)
         !write(6,*) sum_z_l1, ave_z_l1
        
         sum_z_l2=0.d0
         DO i=1,natoms_fin/2
               sum_z_l2=sum_z_l2+rz_fin(layer2_fin(i))
         END DO
         ave_z_l2=sum_z_l2/FLOAT(natoms_fin/2)
         !write(6,*) sum_z_l2, ave_z_l2

         cone_height=ABS(ave_z_l2-ave_z_l1)

         write(6,*) 'c_h',cone_height
         pi=4.d0*datan(1.d0)

!To find the order parameter divide the upper layer in M by M cells and
!find the atoms(atomtypes) in each cell.
! Now search the atoms on upper layer for each hexagon at the lower
! layer
         
         I_dim=x_cell
         J_dim=y_cell
         
         a0=2.54d0   ! B-N hex lattice constant

         ALLOCATE(b_upper_match(tot_hex,6),color_index(tot_hex),&
                 phi(tot_hex),sinphi(tot_hex),cosphi(tot_hex),&
                 mesh_atoms(I_dim,J_dim),item_upper_match(tot_hex,6),&
                 r_upper_match(tot_hex,6,3),theta_up_down(tot_hex,6),&
                 dist_upper_match(tot_hex,6),no_atom_found(tot_hex,6),&
                 hexa_atom_match(tot_hex,6,3),hexa_op(tot_hex),&
                 hexa_rcom_l2(tot_hex,3),hexa_rgeocp_l2(tot_hex,3))

        l1=(/.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE./)
        l2=(/.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE./)
        l3=(/.FALSE.,.TRUE.,.FALSE.,.TRUE.,.FALSE.,.TRUE./)
        l4=(/.TRUE.,.FALSE.,.TRUE.,.FALSE.,.TRUE.,.FALSE./)
        
        OPEN(UNIT=unit_output_1,FILE=hexagon_data_dist_file)
        OPEN(UNIT=unit_output_2,FILE=order_param_file)
        OPEN(UNIT=unit_output_3,FILE=hexagon_data_up_xyz_file)
        OPEN(UNIT=unit_output_4,FILE=hexagon_data_down_xyz_file)
        OPEN(UNIT=rcom_output,FILE=hexagon_rcom_data_file)

        cnt5=0
        cnt6=0
        cnt7=0
        cnt8=0
        cnt9=0        
        cnt10=0
        cnt11=0

        a0_limit=a0/sqrt(3.d0)

        a0_norm=sqrt(3.d0)/2.d0

        DO ii=1,tot_hex
              rcom_hexa_l1(:)=hexa_rcom_l1(ii,:) 
              rgeocp_hexa_l1(:)=hexa_rgeocp_l1(ii,:)
              no_atom_found(ii,:)=l1
              !WRITE(6,*) rcom_hexa_l1
       inner: DO jj=1,6
! Lower layer hexagon points
              item=hexa_l1(ii,jj)
              r_hexa_l1(1)=rx_fin(item)
              r_hexa_l1(2)=ry_fin(item)
              r_hexa_l1(3)=rz_fin(item)
              b_hexa_l1=b_fin(item)  
       
              norm_hexa_l1(:)=hexa_norm_l1(ii,jj,:)     

             ! WRITE(6,*) ii, jj, r_hexa_l1, norm_hexa_l1, rcom_hexa_l1,&
             !            b_hexa_l1,hexa_l1(ii,jj)  
             ! WRITE(6,*) ii, jj, item
                      
              WRITE(6,*) 'Lower Hexagon', ii, jj, item, r_hexa_l1, &
                            b_hexa_l1             
              WRITE(unit_output_4,*) b_hexa_l1, r_hexa_l1           
 
             DO i=1,I_dim
              DO j=1,J_dim
                
                k1=4*J_dim*(i-1)+4*(j-1)
                k2=4*J_dim*(i-1)+4*J_dim
            
                IF(item.gt.k1.and.item.lt.k2) THEN
                
                !WRITE(6,*) i,j,k1,k2,item

              CALL search_atom_up_pbc(item,i,j,r_hexa_l1,norm_hexa_l1,&
                      cone_height,l2_atom_num,r_hexa_l2,l2_atom_b,&
                      angle_two_lay,dist_l2,lhexagon) 
 
              !WRITE(6,*) 'Up_atom',item,i,j,r_hexa_l1,norm_hexa_l1,&
              !        cone_height,l2_atom_num,r_hexa_l2,l2_atom_b,&
              !        angle_two_lay,dist_l2,lhexagon          

 
                IF(lhexagon)THEN
                WRITE(6,*) ii,jj, item,&
                           'No atom found at the upper-layer',lhexagon
                            no_atom_found(ii,jj)=lhexagon
                WRITE(6,*) 'Upper Hexagon',ii,jj,item,'No atom found'

               ! item_upper_match(ii,jj)=0
               ! r_upper_match(ii,jj,1)=0.d0
               ! r_upper_match(ii,jj,2)=0.d0
               ! r_upper_match(ii,jj,3)=0.d0
               ! b_upper_match(ii,jj)='Ud'

               ! WRITE(6,*) 'Upper Hexagon',ii,jj,item,&
               ! item_upper_match(ii,jj),r_upper_match(ii,jj,:),&
               ! b_upper_match(ii,jj)

               ! WRITE(unit_output_3,*) b_upper_match(ii,jj), &
               !                        r_upper_match(ii,jj,:)
               ! 
                ELSE

                item_upper_match(ii,jj)=l2_atom_num
                r_upper_match(ii,jj,1)=r_hexa_l2(1)
                r_upper_match(ii,jj,2)=r_hexa_l2(2)
                r_upper_match(ii,jj,3)=r_hexa_l2(3)
                b_upper_match(ii,jj)=l2_atom_b        
                
                theta_up_down(ii,jj)=angle_two_lay
                dist_upper_match(ii,jj)=dist_l2
                
                no_atom_found(ii,jj)=lhexagon
                
                WRITE(6,*) 'Upper Hexagon',ii,jj,item,&
                item_upper_match(ii,jj),r_upper_match(ii,jj,:),&
                b_upper_match(ii,jj),theta_up_down(ii,jj),&
                dist_upper_match(ii,jj)
               
                WRITE(unit_output_3,*) b_upper_match(ii,jj),&
                                       r_upper_match(ii,jj,:)

               END IF 
                 CYCLE inner
              END IF
             END DO
            END DO
                
                
       END DO inner

      IF(ALL([(ALL(no_atom_found(ii,jj).eqv.l2),jj=1,6)]))THEN
       hexa_rcom_l2(ii,:)=rcom(r_upper_match(ii,:,1),&
                          r_upper_match(ii,:,2),r_upper_match(ii,:,3),&
                          b_upper_match(ii,:))
       hexa_rgeocp_l2(ii,:)=rgeocp(r_upper_match(ii,:,1),&
                          r_upper_match(ii,:,2),r_upper_match(ii,:,3))

       rcom_hexa_l2(:)=hexa_rcom_l2(ii,:) 
       rgeocp_hexa_l2(:)=hexa_rgeocp_l2(ii,:)
       
       WRITE(6,*) 'R_COM:',ii,rcom_hexa_l1,rcom_hexa_l2 
       WRITE(6,*) 'R_COM:',ii,rgeocp_hexa_l1,rgeocp_hexa_l2 

       temp_x=rcom_hexa_l2(1)-rcom_hexa_l1(1)
       temp_y=rcom_hexa_l2(2)-rcom_hexa_l1(2)
       temp_z=rcom_hexa_l2(3)-rcom_hexa_l1(3)-cone_height
       d_rcom=SQRT(temp_x**2+temp_y**2)
       WRITE(rcom_output,*) ii,rcom_hexa_l1,rcom_hexa_l2,d_rcom

       WRITE(6,*) ii,rcom_hexa_l2(:),no_atom_found(ii,:),d_rcom
       IF (d_rcom==0.d0) THEN
           color_index(ii)=0.d0
       ELSE
           cosphi(ii)=temp_x/d_rcom 
           sinphi(ii)=temp_y/d_rcom
           WRITE(6,*) ii,cosphi(ii),sinphi(ii) 
           IF (cosphi(ii).ge.0.d0) THEN
                IF (sinphi(ii).ge.0.d0)THEN
!1st Quadrant
                        phi(ii)=dacos((temp_x/d_rcom)*pi/180.d0)
                    WRITE(6,*)temp_x/d_rcom,180.d0*(temp_x/d_rcom)/pi
                ELSEIF(sinphi(ii).lt.0.d0) THEN
!4th Quadrant
                        phi(ii)=dacos(-(temp_x/d_rcom)*pi/180.d0)
                    WRITE(6,*) temp_x/d_rcom,180.d0*(temp_x/d_rcom)/pi
                END IF
           ELSEIF (cosphi(ii).lt.0.d0) THEN
                IF (sinphi(ii).ge.0.d0)THEN
!2nd Quadrant
                        phi(ii)=dacos((temp_x/d_rcom)*pi/180.d0)     
                    WRITE(6,*) temp_x/d_rcom,180.d0*(temp_x/d_rcom)/pi
                ELSEIF(sinphi(ii).lt.0.d0) THEN
!3rd Quadrant
                        phi(ii)=dacos(-(temp_x/d_rcom)*pi/180.d0)
                    WRITE(6,*) temp_x/d_rcom,180.d0*(temp_x/d_rcom)/pi
                END IF
           END IF
           WRITE(6,*) ii, cosphi(ii), sinphi(ii), phi(ii)

!#Case I:
           IF((-15.d0.lt.phi(ii).and.phi(ii).le.0.d0).OR.&
               (0.d0.lt.phi(ii).and.phi(ii).le.15.d0).OR. &
              (105.d0.lt.phi(ii).and.phi(ii).le.135.d0).OR.&
              (-75.d0.lt.phi(ii).and.phi(ii).lt.-35.d0)) THEN
               IF(0.d0.lt.d_rcom.and.d_rcom.le.a0_limit) THEN
                    color_index(ii)=d_rcom/a0_limit
               END IF
!#Case II: 
           ELSEIF((45.d0.lt.phi(ii).and.phi(ii).le.75.d0).OR.&
              (165.d0.lt.phi(ii).and.phi(ii).le.180.d0).OR.&
             (-180.d0.lt.phi(ii).and.phi(ii).le.-165.d0).OR.&
             (-135.d0.lt.phi(ii).and.phi(ii).le.-105.d0)) THEN
               IF(0.d0.lt.d_rcom.and.d_rcom.le.a0_limit)THEN
                    color_index(ii)=-d_rcom/a0_limit
               END IF
!#Case III:
            ELSEIF((15.d0.lt.phi(ii).and.phi(ii).le.45.d0).OR.&
               (135.d0.lt.phi(ii).and.phi(ii).le.165.d0).OR.&
              (-105.d0.lt.phi(ii).and.phi(ii).le.-75.d0))THEN
                IF(0.d0.lt.(d_rcom*a0_norm).and. &
                           (d_rcom*a0_norm).le.(a0/2.d0))THEN
                    color_index(ii)=((d_rcom*a0_norm)/a0_limit)
                END IF
!#Case IV:
            ELSEIF(75.d0.lt.phi(ii).and.phi(ii).le.105.d0.OR.&
              -45.d0.lt.phi(ii).and.phi(ii).le.-15.d0.OR.&
              -165.d0.lt.phi(ii).and.phi(ii).le.-135.d0)THEN
                IF(0.d0.lt.(d_rcom*a0_norm).and. &
                           (d_rcom*a0_norm).le.(a0/2.d0))THEN
                     color_index(ii)=(-(d_rcom*a0_norm)/a0_limit)
                END IF
            END IF
       END IF
       WRITE(6,*)ii,color_index(ii)

       WRITE(unit_output_2,*) ii, hexa_rcom_l1(ii,:), &
                                  hexa_rcom_l2(ii,:), &
                                  d_rcom,color_index(ii)

       WRITE(unit_output_3,*) b_upper_match(ii,:),&
                                       r_upper_match(ii,:,:)

      END IF

      END DO

      CLOSE(unit_output_1)
      CLOSE(unit_output_2)
 
      CLOSE(unit_output_3)
      CLOSE(unit_output_4)
      CLOSE(rcom_output)
        

         DEALLOCATE(b_upper_match,mesh_atoms,item_upper_match,&
                 r_upper_match,theta_up_down,hexa_atom_match,hexa_op)

        
        CALL cpu_time(finish)
    
        WRITE(6,*) '("Time taken for order_parameter = ",f6.3,&
                                    " seconds.")',finish-start
    
         END SUBROUTINE order_parameter

        SUBROUTINE find_hexagon_lower_layer

!Here the idea is to find the hexagons from the unit cell. Here unit
!cell dimensions needs to be mentioned

        USE global_variable
        USE my_functions

        IMPLICIT NONE

        INTEGER*4::i,j,J_dim,I_dim,cnt,J_hex,J_hex_prim,item,item1,&
                   item2,item3,item4
      
        REAL*8::v1(3),v2(3),v3(3),v4(3),v_norm(3)

        REAL*8,DIMENSION(:),ALLOCATABLE::x_hexa,y_hexa,z_hexa,rcom_hexa
        REAL*8,DIMENSION(:),ALLOCATABLE::x_hexa_l2,y_hexa_l2,z_hexa_l2
        REAL*8,DIMENSION(:,:),ALLOCATABLE::hexa_l2
                                         

        CHARACTER,DIMENSION(:),ALLOCATABLE::b_hexa,b_hexa_l2
 
        REAL*8 :: start, finish

        CALL cpu_time(start)
        
        J_dim=y_cell

        I_dim=x_cell

        J_hex=INT(FLOOR((J_dim)/3.d0))+INT(FLOOR((J_dim-1)/3.d0))
 
        tot_hex=J_hex*(I_dim-1)
         
        WRITE(6,*) J_hex,I_dim-1,tot_hex

        cnt=0

        ALLOCATE(hexa_l1(tot_hex,6),hexa_rcom_l1(tot_hex,3),&
                 hexa_rgeocp_l1(tot_hex,3),hexa_l2(tot_hex,6))

!In J direction we look for the hexagons. For I=1 and I=M-1, consider
!every isolated hexagon at lower layer after 3 unit cells each 

         DO i=1,I_dim-1

                IF(MOD(J_hex,2).EQ.0) THEN
 
                DO j=0,J_hex-1,2

                  cnt=cnt+1
         
                  hexa_l1(cnt,1)=4*J_dim*(i-1)+3+12*(j/2)
                  hexa_l1(cnt,2)=4*J_dim*(i-1)+4+12*(j/2)
                  hexa_l1(cnt,3)=4*J_dim*i+5+12*(j/2)
                  hexa_l1(cnt,4)=4*J_dim*(i-1)+8+12*(j/2)
                  hexa_l1(cnt,5)=4*J_dim*(i-1)+7+12*(j/2)
                  hexa_l1(cnt,6)=4*J_dim*(i-1)+6+12*(j/2)

                  WRITE(6,*)cnt,i,j+1,hexa_l1(cnt,1:6)

                  cnt=cnt+1

                  hexa_l1(cnt,1)=4*J_dim*i+9+12*(j/2)
                  hexa_l1(cnt,2)=4*J_dim*i+10+12*(j/2)
                  hexa_l1(cnt,3)=4*J_dim*i+11+12*(j/2)
                  hexa_l1(cnt,4)=4*J_dim*i+14+12*(j/2)
                  hexa_l1(cnt,5)=4*J_dim*i+13+12*(j/2)
                  hexa_l1(cnt,6)=4*J_dim*(i-1)+12+12*(j/2)

                       
                  WRITE(6,*)cnt,i,j+2,hexa_l1(cnt,1:6)     
               
                END DO 

                ELSE
 
                DO j=0,J_hex-2,2

                  cnt=cnt+1
         
                  hexa_l1(cnt,1)=4*J_dim*(i-1)+3+12*(j/2)
                  hexa_l1(cnt,2)=4*J_dim*(i-1)+4+12*(j/2)
                  hexa_l1(cnt,3)=4*J_dim*i+5+12*(j/2)
                  hexa_l1(cnt,4)=4*J_dim*(i-1)+8+12*(j/2)
                  hexa_l1(cnt,5)=4*J_dim*(i-1)+7+12*(j/2)
                  hexa_l1(cnt,6)=4*J_dim*(i-1)+6+12*(j/2)

                  WRITE(6,*)cnt,i,j+1,hexa_l1(cnt,1:6)

                  cnt=cnt+1

                  hexa_l1(cnt,1)=4*J_dim*i+9+12*(j/2)
                  hexa_l1(cnt,2)=4*J_dim*i+10+12*(j/2)
                  hexa_l1(cnt,3)=4*J_dim*i+11+12*(j/2)
                  hexa_l1(cnt,4)=4*J_dim*i+14+12*(j/2)
                  hexa_l1(cnt,5)=4*J_dim*i+13+12*(j/2)
                  hexa_l1(cnt,6)=4*J_dim*(i-1)+12+12*(j/2)

                       
                  WRITE(6,*)cnt,i,j+2,hexa_l1(cnt,1:6)     
               
                END DO 
       
                cnt=cnt+1
                
                j=J_hex-1
         
                hexa_l1(cnt,1)=4*J_dim*(i-1)+3+12*(j/2)
                hexa_l1(cnt,2)=4*J_dim*(i-1)+4+12*(j/2)
                hexa_l1(cnt,3)=4*J_dim*i+5+12*(j/2)
                hexa_l1(cnt,4)=4*J_dim*(i-1)+8+12*(j/2)
                hexa_l1(cnt,5)=4*J_dim*(i-1)+7+12*(j/2)
                hexa_l1(cnt,6)=4*J_dim*(i-1)+6+12*(j/2)

                WRITE(6,*)cnt,i,j+1,hexa_l1(cnt,1:6)

                END IF 

         END DO  
      

! Find the COM for each hexagon
        open(unit=111,file="initial_hexagon_test.xyz")
        open(unit=112,file="initial_hexagon_com_test.xyz")


         DO i=1,tot_hex
           
           ALLOCATE(x_hexa(6),y_hexa(6),z_hexa(6),b_hexa(6),&
                    x_hexa_l2(6),y_hexa_l2(6),z_hexa_l2(6),&
                    b_hexa_l2(6))    

           DO j=1,6
        
                item=hexa_l1(i,j)
              
                item1=hexa_l1(i,j)+natoms_fin/2

                
                x_hexa(j)=rx_fin(item)
                y_hexa(j)=ry_fin(item)
                z_hexa(j)=rz_fin(item)
                b_hexa(j)=b_fin(item)

                x_hexa_l2(j)=rx_fin(item1)
                y_hexa_l2(j)=ry_fin(item1)
                z_hexa_l2(j)=rz_fin(item1)
                b_hexa_l2(j)=b_fin(item1)
                
               WRITE(6,*) i,j,hexa_l1(i,j),item,x_hexa(j),&
                          y_hexa(j),z_hexa(j),b_hexa(j)
           
               WRITE(111,*) i,j, b_hexa(j),x_hexa(j),y_hexa(j),z_hexa(j)
               WRITE(111,*) i,j, b_hexa_l2(j),x_hexa_l2(j),&
                                 y_hexa_l2(j),z_hexa_l2(j)

           END DO

           !WRITE(111,*) i, x_hexa(:),y_hexa(:),z_hexa(:),&
           !                x_hexa_l2(:),y_hexa_l2(:),z_hexa_l2(:)


           hexa_rcom_l1(i,:)=rcom(x_hexa(:),y_hexa(:),z_hexa(:),&
                                                         b_hexa(:))

!           hexa_rcom_l2(i,:)=rcom(x_hexa_l2(:),y_hexa_l2(:),&
!                                     z_hexa_l2(:),b_hexa_l2(:))
!
!           WRITE(112,*) i,hexa_rcom_l1(i,:),hexa_rcom_l2(i,:)
!
!           hexa_rgeocp_l1(i,:)=rgeocp(x_hexa(:),y_hexa(:),z_hexa(:))
!
!           WRITE(6,*) 'hexagon_com',i, hexa_l1(i,:), hexa_rcom_l1(i,:)
!
!           WRITE(6,*) 'hexagon_geocp',i,hexa_l1(i,:),hexa_rgeocp_l1(i,:)
             
           DEALLOCATE(x_hexa,y_hexa,z_hexa,b_hexa,x_hexa_l2,y_hexa_l2,&
                      z_hexa_l2,b_hexa_l2)                

         END DO
        
         CLOSE(111)
         CLOSE(112)

        
!Find normal for each of the points of the hexagon

          ALLOCATE(hexa_norm_l1(tot_hex,6,3))

          DO i=1,tot_hex

             DO j=1,6

                item=hexa_l1(i,j)
                item1=nn4_sorted_l1(item,1)
                item2=nn4_sorted_l1(item,2)
                item3=nn4_sorted_l1(item,3)
                item4=nn4_sorted_l1(item,4)

                v1(1)=rx_fin(item1)
                v1(2)=ry_fin(item1)
                v1(3)=rz_fin(item1)
                
                v2(1)=rx_fin(item2)
                v2(2)=ry_fin(item2)
                v2(3)=rz_fin(item2)

                v3(1)=rx_fin(item3)
                v3(2)=ry_fin(item3)
                v3(3)=rz_fin(item3)
                
                v4(1)=rx_fin(item4)
                v4(2)=ry_fin(item4)
                v4(3)=rz_fin(item4)

                CALL Normal(v1(:),v2(:),v3(:),v4(:),v_norm(:))

                hexa_norm_l1(i,j,:)=v_norm(:)

                !WRITE(6,*) i,j,item,item1,item2,item3,item4
                !WRITE(6,*) v1(:),v2(:),v3(:),v4(:),v_norm(:)
                WRITE(6,*) 'Hexa_norm',i,j,hexa_norm_l1(i,j,:)

            END DO

          END DO
        
        CALL cpu_time(finish)
    
        WRITE(6,*) '("Time taken for find_hexagon_lower_layer = ",f6.3,&
                                    " seconds.")',finish-start
    
          
       END SUBROUTINE find_hexagon_lower_layer 

        SUBROUTINE nearest_neighbor_l2

        USE global_variable
        USE my_functions
        
        IMPLICIT NONE

        INTEGER*4::i,j,k,i1,i2,i3,i4,i5,i6,ii,jj,atom_num(12),item1, &
                   item2,item3,cnt,cnt1,I_dim,J_dim

        REAL*8:: rBN

        REAL*8 :: start, finish

        INTEGER,DIMENSION(:),ALLOCATABLE:: atom_num_radius

        CALL cpu_time(start)

        ALLOCATE(nn4_sorted_l2(size(layer2_fin),12))
       
        !rBN=rBN_l1+0.5d0
        
        rBN=1.60d0

        WRITE(6,*) rBN

        I_dim=x_cell
        J_dim=y_cell

        nn4_sorted_l2=0

        k=natoms_fin/2

!For j.eq 1 and i.eq.1

         j=1
         i=1
                          
         ALLOCATE(atom_num_radius(16))
         i1=4*J_dim*(i-1)+4*(j-1)+1
         i2=4*J_dim*(i-1)+4*(j+1)
         i3=4*J_dim*i+4*(j-1)+1
         i4=4*J_dim*i+4*(j+1)
         
        !WRITE(6,*) i1,i2,i3,i4

         cnt=0
         DO ii=i1,i2
               cnt=cnt+1
               atom_num_radius(cnt)=ii
        !       WRITE(6,*) cnt, atom_num_radius(cnt),ii
         END DO
         DO ii=i3,i4
               cnt=cnt+1
               atom_num_radius(cnt)=ii
         !      WRITE(6,*) cnt, atom_num_radius(cnt),ii
         END DO

        !WRITE(6,*) i,j, atom_num_radius(:)
              
        DO ii=1,4
           item1=atom_num_radius(ii)
           nn4_sorted_l2(atom_num_radius(ii),1)=item1+k
          ! WRITE(6,*) item1,nn4_sorted_l1(atom_num_radius(ii),1)
           cnt1=1
          DO jj=1,size(atom_num_radius)
             item2=atom_num_radius(jj)
             IF(item1.ne.item2) THEN
              IF(dist(rx_fin(item1+k),rx_fin(item2+k),ry_fin(item1+k), &
                      ry_fin(item2+k),rz_fin(item1+k),rz_fin(item2+k),&
                      xlo_fin,xhi_fin,ylo_fin,yhi_fin,zlo_fin,&
                      zhi_fin).LT.rBN) THEN
                      cnt1=cnt1+1 
                      nn4_sorted_l2(atom_num_radius(ii),cnt1)=item2+k
         !  WRITE(6,*) cnt1,item1,nn4_sorted_l1(atom_num_radius(ii),cnt1)
         !       ! WRITE(6,*)'inner_inner_n4','ii',ii,'jj',jj,'kk',&
         !       !            kk,cnt,i1,i2,item3
              END IF
             END IF
              ! WRITE(6,*) 'inner_n4','i',ii,'j',jj,'k',kk,&
              !             nn4_sorted_l2(atom_num(k),:)
           END DO
         !    WRITE(6,*) item1, nn4_sorted_l1(atom_num_radius(ii),:)
         END DO
       
       DEALLOCATE(atom_num_radius)

!For i.eq 1 and 2<=j<=J_dim-1

         i=1
         DO j=2,J_dim-1

         ALLOCATE(atom_num_radius(24))
         i1=4*J_dim*(i-1)+4*(j-2)+1
         i2=4*J_dim*(i-1)+4*(j+1)
         i3=4*J_dim*i+4*(j-2)+1
         i4=4*J_dim*i+4*(j+1)
        
         !WRITE(6,*) i1,i2,i3,i4
         
         cnt=0
         DO ii=i1,i2
               cnt=cnt+1
               atom_num_radius(cnt)=ii
         END DO
         DO ii=i3,i4
               cnt=cnt+1
               atom_num_radius(cnt)=ii
         END DO

        !WRITE(6,*) i,j, atom_num_radius
              
        DO ii=5,8
           item1=atom_num_radius(ii)
           nn4_sorted_l2(atom_num_radius(ii),1)=item1+k
          ! WRITE(6,*) item1,nn4_sorted_l1(atom_num_radius(ii),1)
           cnt1=1
          DO jj=1,size(atom_num_radius)
             item2=atom_num_radius(jj)
             IF(item1.ne.item2) THEN
              IF(dist(rx_fin(item1+k),rx_fin(item2+k),ry_fin(item1+k), &
                      ry_fin(item2+k),rz_fin(item1+k),rz_fin(item2+k),&
                      xlo_fin,xhi_fin,ylo_fin,yhi_fin,zlo_fin,&
                      zhi_fin).LT.rBN) THEN
                      cnt1=cnt1+1 
                      nn4_sorted_l2(atom_num_radius(ii),cnt1)=item2+k
         !  WRITE(6,*) cnt1,item1,nn4_sorted_l1(atom_num_radius(ii),cnt1)
         !       ! WRITE(6,*)'inner_inner_n4','ii',ii,'jj',jj,'kk',&
         !       !            kk,cnt,i1,i2,item3
              END IF
             END IF
              ! WRITE(6,*) 'inner_n4','i',ii,'j',jj,'k',kk,&
              !             nn4_sorted_l2(atom_num(k),:)
           END DO
         !    WRITE(6,*) item1, nn4_sorted_l1(atom_num_radius(ii),:)
         END DO

         DEALLOCATE(atom_num_radius)

         END DO

!For j.eq 1 and 2<=i<=I_dim-1

         j=1
         DO i=2,I_dim-1

         ALLOCATE(atom_num_radius(24))
         i1=4*J_dim*(i-2)+4*(j-1)+1
         i2=4*J_dim*(i-2)+4*(j+1)
         i3=4*J_dim*(i-1)+4*(j-1)+1
         i4=4*J_dim*(i-1)+4*(j+1)
         i5=4*J_dim*i+4*(j-1)+1
         i6=4*J_dim*i+4*(j+1)
         
         !WRITE(6,*) i1,i2,i3,i4,i5,i6
         
         cnt=0
         DO ii=i1,i2
               cnt=cnt+1
               atom_num_radius(cnt)=ii
         END DO
         DO ii=i3,i4
               cnt=cnt+1
               atom_num_radius(cnt)=ii
         END DO
         DO ii=i5,i6
               cnt=cnt+1
               atom_num_radius(cnt)=ii
         END DO

        !WRITE(6,*) i,j, atom_num_radius
              
        DO ii=9,12
           item1=atom_num_radius(ii)
           nn4_sorted_l2(atom_num_radius(ii),1)=item1+k
          ! WRITE(6,*) item1,nn4_sorted_l1(atom_num_radius(ii),1)
           cnt1=1
          DO jj=1,size(atom_num_radius)
             item2=atom_num_radius(jj)
             IF(item1.ne.item2) THEN
              IF(dist(rx_fin(item1+k),rx_fin(item2+k),ry_fin(item1+k), &
                      ry_fin(item2+k),rz_fin(item1+k),rz_fin(item2+k),&
                      xlo_fin,xhi_fin,ylo_fin,yhi_fin,zlo_fin,&
                      zhi_fin).LT.rBN) THEN
                      cnt1=cnt1+1 
                      nn4_sorted_l2(atom_num_radius(ii),cnt1)=item2+k
         !  WRITE(6,*) cnt1,item1,nn4_sorted_l1(atom_num_radius(ii),cnt1)
         !       ! WRITE(6,*)'inner_inner_n4','ii',ii,'jj',jj,'kk',&
         !       !            kk,cnt,i1,i2,item3
              END IF
             END IF
              ! WRITE(6,*) 'inner_n4','i',ii,'j',jj,'k',kk,&
              !             nn4_sorted_l2(atom_num(k),:)
           END DO
         !    WRITE(6,*) item1, nn4_sorted_l1(atom_num_radius(ii),:)
         END DO

         DEALLOCATE(atom_num_radius)

         END DO


! For 2 <= i <= I_dim-1 and 2 <= j <= J_dim-1

        DO i=2,I_dim-1
         DO j=2,J_dim-1
             
         ALLOCATE(atom_num_radius(36))
         i1=4*J_dim*(i-2)+4*(j-2)+1
         i2=4*J_dim*(i-2)+4*(j+1)
         i3=4*J_dim*(i-1)+4*(j-2)+1
         i4=4*J_dim*(i-1)+4*(j+1)
         i5=4*J_dim*i+4*(j-2)+1
         i6=4*J_dim*i+4*(j+1)
         
         !WRITE(6,*) i1,i2,i3,i4,i5,i6
        
         cnt=0
         DO ii=i1,i2
               cnt=cnt+1
               atom_num_radius(cnt)=ii
         END DO
         DO ii=i3,i4
               cnt=cnt+1
               atom_num_radius(cnt)=ii
         END DO
         DO ii=i5,i6
               cnt=cnt+1
               atom_num_radius(cnt)=ii
         END DO

        !WRITE(6,*) i,j, atom_num_radius
              
        DO ii=17,20
           item1=atom_num_radius(ii)
           nn4_sorted_l2(atom_num_radius(ii),1)=item1+k
          ! WRITE(6,*) item1,nn4_sorted_l1(atom_num_radius(ii),1)
           cnt1=1
          DO jj=1,size(atom_num_radius)
             item2=atom_num_radius(jj)
             IF(item1.ne.item2) THEN
              IF(dist(rx_fin(item1+k),rx_fin(item2+k),ry_fin(item1+k), &
                      ry_fin(item2+k),rz_fin(item1+k),rz_fin(item2+k),&
                      xlo_fin,xhi_fin,ylo_fin,yhi_fin,zlo_fin,&
                      zhi_fin).LT.rBN) THEN
                      cnt1=cnt1+1 
                      nn4_sorted_l2(atom_num_radius(ii),cnt1)=item2+k
         !  WRITE(6,*) cnt1,item1,nn4_sorted_l1(atom_num_radius(ii),cnt1)
         !       ! WRITE(6,*)'inner_inner_n4','ii',ii,'jj',jj,'kk',&
         !       !            kk,cnt,i1,i2,item3
              END IF
             END IF
              ! WRITE(6,*) 'inner_n4','i',ii,'j',jj,'k',kk,&
              !             nn4_sorted_l2(atom_num(k),:)
           END DO
         !    WRITE(6,*) item1, nn4_sorted_l1(atom_num_radius(ii),:)
         END DO

          DEALLOCATE(atom_num_radius)

         END DO
        END DO
                                
!For i=I_dim and j=1

         i=I_dim
         j=1
          ALLOCATE(atom_num_radius(16))
          i1=4*J_dim*(i-2)+4*(j-1)+1
          i2=4*J_dim*(i-2)+4*(j+1)
          i3=4*J_dim*(i-1)+4*(j-1)+1
          i4=4*J_dim*(i-1)+4*(j+1)

        ! WRITE(6,*) i1,i2,i3,i4
        
          cnt=0
          DO ii=i1,i2
                cnt=cnt+1
                atom_num_radius(cnt)=ii
          END DO
          DO ii=i3,i4
                cnt=cnt+1
                atom_num_radius(cnt)=ii
          END DO

        !WRITE(6,*) i,j, atom_num_radius
              
        DO ii=9,12
           item1=atom_num_radius(ii)
           nn4_sorted_l2(atom_num_radius(ii),1)=item1+k
          ! WRITE(6,*) item1,nn4_sorted_l1(atom_num_radius(ii),1)
           cnt1=1
          DO jj=1,size(atom_num_radius)
             item2=atom_num_radius(jj)
             IF(item1.ne.item2) THEN
              IF(dist(rx_fin(item1+k),rx_fin(item2+k),ry_fin(item1+k), &
                      ry_fin(item2+k),rz_fin(item1+k),rz_fin(item2+k),&
                      xlo_fin,xhi_fin,ylo_fin,yhi_fin,zlo_fin,&
                      zhi_fin).LT.rBN) THEN
                      cnt1=cnt1+1 
                      nn4_sorted_l2(atom_num_radius(ii),cnt1)=item2+k
         !  WRITE(6,*) cnt1,item1,nn4_sorted_l1(atom_num_radius(ii),cnt1)
         !       ! WRITE(6,*)'inner_inner_n4','ii',ii,'jj',jj,'kk',&
         !       !            kk,cnt,i1,i2,item3
              END IF
             END IF
              ! WRITE(6,*) 'inner_n4','i',ii,'j',jj,'k',kk,&
              !             nn4_sorted_l2(atom_num(k),:)
           END DO
         !    WRITE(6,*) item1, nn4_sorted_l1(atom_num_radius(ii),:)
         END DO

          DEALLOCATE(atom_num_radius)


!For i=I_dim and 2<=j<=J_dim-1

         i=I_dim
         DO j=2,J_dim-1
          ALLOCATE(atom_num_radius(24))
          i1=4*J_dim*(i-2)+4*(j-2)+1
          i2=4*J_dim*(i-2)+4*(j+1)
          i3=4*J_dim*(i-1)+4*(j-2)+1
          i4=4*J_dim*(i-1)+4*(j+1)

         !WRITE(6,*) i1,i2,i3,i4
        
          cnt=0
          DO ii=i1,i2
                cnt=cnt+1
                atom_num_radius(cnt)=ii
          END DO
          DO ii=i3,i4
                cnt=cnt+1
                atom_num_radius(cnt)=ii
          END DO

        !WRITE(6,*) i,j, atom_num_radius
              
        DO ii=17,20
           item1=atom_num_radius(ii)
           nn4_sorted_l2(atom_num_radius(ii),1)=item1+k
          ! WRITE(6,*) item1,nn4_sorted_l1(atom_num_radius(ii),1)
           cnt1=1
          DO jj=1,size(atom_num_radius)
             item2=atom_num_radius(jj)
             IF(item1.ne.item2) THEN
              IF(dist(rx_fin(item1+k),rx_fin(item2+k),ry_fin(item1+k), &
                      ry_fin(item2+k),rz_fin(item1+k),rz_fin(item2+k),&
                      xlo_fin,xhi_fin,ylo_fin,yhi_fin,zlo_fin,&
                      zhi_fin).LT.rBN) THEN
                      cnt1=cnt1+1 
                      nn4_sorted_l2(atom_num_radius(ii),cnt1)=item2+k
         !  WRITE(6,*) cnt1,item1,nn4_sorted_l1(atom_num_radius(ii),cnt1)
         !       ! WRITE(6,*)'inner_inner_n4','ii',ii,'jj',jj,'kk',&
         !       !            kk,cnt,i1,i2,item3
              END IF
             END IF
              ! WRITE(6,*) 'inner_n4','i',ii,'j',jj,'k',kk,&
              !             nn4_sorted_l2(atom_num(k),:)
           END DO
         !    WRITE(6,*) item1, nn4_sorted_l1(atom_num_radius(ii),:)
         END DO

          DEALLOCATE(atom_num_radius)
         END DO

!For i=1 and j=J_dim

         j=J_dim
         i=1

          ALLOCATE(atom_num_radius(16))
          i1=4*J_dim*(i-1)+4*(j-2)+1
          i2=4*J_dim*(i-1)+4*j
          i3=4*J_dim*i+4*(j-2)+1
          i4=4*J_dim*i+4*j
          
          !WRITE(6,*) i1,i2,i3,i4

          cnt=0
          DO ii=i1,i2
                cnt=cnt+1
                atom_num_radius(cnt)=ii
          END DO
          DO ii=i3,i4
                cnt=cnt+1
                atom_num_radius(cnt)=ii
          END DO

        !WRITE(6,*) i,j, atom_num_radius
              
        DO ii=5,8
           item1=atom_num_radius(ii)
           nn4_sorted_l2(atom_num_radius(ii),1)=item1+k
          ! WRITE(6,*) item1,nn4_sorted_l1(atom_num_radius(ii),1)
           cnt1=1
          DO jj=1,size(atom_num_radius)
             item2=atom_num_radius(jj)
             IF(item1.ne.item2) THEN
              IF(dist(rx_fin(item1+k),rx_fin(item2+k),ry_fin(item1+k), &
                      ry_fin(item2+k),rz_fin(item1+k),rz_fin(item2+k),&
                      xlo_fin,xhi_fin,ylo_fin,yhi_fin,zlo_fin,&
                      zhi_fin).LT.rBN) THEN
                      cnt1=cnt1+1 
                      nn4_sorted_l2(atom_num_radius(ii),cnt1)=item2+k
         !  WRITE(6,*) cnt1,item1,nn4_sorted_l1(atom_num_radius(ii),cnt1)
         !       ! WRITE(6,*)'inner_inner_n4','ii',ii,'jj',jj,'kk',&
         !       !            kk,cnt,i1,i2,item3
              END IF
             END IF
              ! WRITE(6,*) 'inner_n4','i',ii,'j',jj,'k',kk,&
              !             nn4_sorted_l2(atom_num(k),:)
           END DO
         !    WRITE(6,*) item1, nn4_sorted_l1(atom_num_radius(ii),:)
         END DO

          DEALLOCATE(atom_num_radius)

!For 2<=i<=I_dim-1 and j=J_dim

         j=J_dim
         DO i=2,I_dim-1

          ALLOCATE(atom_num_radius(24))
          i1=4*J_dim*(i-2)+4*(j-2)+1
          i2=4*J_dim*(i-2)+4*j
          i3=4*J_dim*(i-1)+4*(j-2)+1
          i4=4*J_dim*(i-1)+4*j
          i5=4*J_dim*i+4*(j-2)+1
          i6=4*J_dim*i+4*j
          
          !WRITE(6,*) i1,i2,i3,i4,i5,i6

          cnt=0
          DO ii=i1,i2
                cnt=cnt+1
                atom_num_radius(cnt)=ii
          END DO
          DO ii=i3,i4
                cnt=cnt+1
                atom_num_radius(cnt)=ii
          END DO
          DO ii=i5,i6
                cnt=cnt+1
                atom_num_radius(cnt)=ii
          END DO

        !WRITE(6,*) i,j, atom_num_radius
              
        DO ii=13,16
           item1=atom_num_radius(ii)
           nn4_sorted_l2(atom_num_radius(ii),1)=item1+k
          ! WRITE(6,*) item1,nn4_sorted_l1(atom_num_radius(ii),1)
           cnt1=1
          DO jj=1,size(atom_num_radius)
             item2=atom_num_radius(jj)
             IF(item1.ne.item2) THEN
              IF(dist(rx_fin(item1+k),rx_fin(item2+k),ry_fin(item1+k), &
                      ry_fin(item2+k),rz_fin(item1+k),rz_fin(item2+k),&
                      xlo_fin,xhi_fin,ylo_fin,yhi_fin,zlo_fin,&
                      zhi_fin).LT.rBN) THEN
                      cnt1=cnt1+1 
                      nn4_sorted_l2(atom_num_radius(ii),cnt1)=item2+k
         !  WRITE(6,*) cnt1,item1,nn4_sorted_l1(atom_num_radius(ii),cnt1)
         !       ! WRITE(6,*)'inner_inner_n4','ii',ii,'jj',jj,'kk',&
         !       !            kk,cnt,i1,i2,item3
              END IF
             END IF
              ! WRITE(6,*) 'inner_n4','i',ii,'j',jj,'k',kk,&
              !             nn4_sorted_l2(atom_num(k),:)
           END DO
         !    WRITE(6,*) item1, nn4_sorted_l1(atom_num_radius(ii),:)
         END DO

          DEALLOCATE(atom_num_radius)
         END DO


!For i=I_dim and j=J_dim
           
          i=I_dim
          j=J_dim
      
          ALLOCATE(atom_num_radius(16))
          i1=4*J_dim*(i-2)+4*(j-2)+1
          i2=4*J_dim*(i-2)+4*j
          i3=4*J_dim*(i-1)+4*(j-2)+1
          i4=4*J_dim*(i-1)+4*j
         
          !WRITE(6,*) i1,i2,i3,i4

          cnt=0
          DO ii=i1,i2
                cnt=cnt+1
                atom_num_radius(cnt)=ii
          END DO
          DO ii=i3,i4
                cnt=cnt+1
                atom_num_radius(cnt)=ii
          END DO
      
        !WRITE(6,*) i,j, atom_num_radius
              
        DO ii=13,16
           item1=atom_num_radius(ii)
           nn4_sorted_l2(atom_num_radius(ii),1)=item1+k
          ! WRITE(6,*) item1,nn4_sorted_l1(atom_num_radius(ii),1)
           cnt1=1
          DO jj=1,size(atom_num_radius)
             item2=atom_num_radius(jj)
             IF(item1.ne.item2) THEN
              IF(dist(rx_fin(item1+k),rx_fin(item2+k),ry_fin(item1+k), &
                      ry_fin(item2+k),rz_fin(item1+k),rz_fin(item2+k),&
                      xlo_fin,xhi_fin,ylo_fin,yhi_fin,zlo_fin,&
                      zhi_fin).LT.rBN) THEN
                      cnt1=cnt1+1 
                      nn4_sorted_l2(atom_num_radius(ii),cnt1)=item2+k
         !  WRITE(6,*) cnt1,item1,nn4_sorted_l1(atom_num_radius(ii),cnt1)
         !       ! WRITE(6,*)'inner_inner_n4','ii',ii,'jj',jj,'kk',&
         !       !            kk,cnt,i1,i2,item3
              END IF
             END IF
              ! WRITE(6,*) 'inner_n4','i',ii,'j',jj,'k',kk,&
              !             nn4_sorted_l2(atom_num(k),:)
           END DO
         !    WRITE(6,*) item1, nn4_sorted_l1(atom_num_radius(ii),:)
         END DO

          DEALLOCATE(atom_num_radius)

       DO i=1,size(layer2_fin)  

     WRITE(6,*)'nearest_neighbor_l2', i+natoms_fin/2, nn4_sorted_l2(i,:)
             
       END DO 
        
        CALL cpu_time(finish)
    
        WRITE(6,*) '("Time taken for nearest_neighbor_l2 = ",f6.3,&
                                    " seconds.")',finish-start
    
        
        END SUBROUTINE nearest_neighbor_l2

        SUBROUTINE nearest_neighbor_l1

        USE global_variable
        USE my_functions
        
        IMPLICIT NONE

        INTEGER*4::i,j,k,i1,i2,i3,i4,i5,i6,ii,jj,atom_num(12),item1, &
                   item2,item3,cnt,cnt1,I_dim,J_dim

        REAL*8:: rBN

        REAL*8 :: start, finish

        INTEGER,DIMENSION(:),ALLOCATABLE:: atom_num_radius

        CALL cpu_time(start)

        ALLOCATE(nn4_sorted_l1(natoms_fin/2,12))
       
        !rBN=rBN_l1+0.5d0
        
        rBN=1.60d0

        WRITE(6,*) rBN

        I_dim=x_cell
        J_dim=y_cell

        nn4_sorted_l1=0

        

!For j.eq 1 and i.eq.1

         j=1
         i=1
                          
         ALLOCATE(atom_num_radius(16))
         i1=4*J_dim*(i-1)+4*(j-1)+1
         i2=4*J_dim*(i-1)+4*(j+1)
         i3=4*J_dim*i+4*(j-1)+1
         i4=4*J_dim*i+4*(j+1)
         
       ! WRITE(6,*) i1,i2,i3,i4

         cnt=0
         DO ii=i1,i2
               cnt=cnt+1
               atom_num_radius(cnt)=ii
        !       WRITE(6,*) cnt, atom_num_radius(cnt),ii
         END DO
         DO ii=i3,i4
               cnt=cnt+1
               atom_num_radius(cnt)=ii
         !      WRITE(6,*) cnt, atom_num_radius(cnt),ii
         END DO

        !WRITE(6,*) i,j, atom_num_radius(:)
              
        DO ii=1,4
           item1=atom_num_radius(ii)
           nn4_sorted_l1(atom_num_radius(ii),1)=item1
          ! WRITE(6,*) item1,nn4_sorted_l1(atom_num_radius(ii),1)
           cnt1=1
          DO jj=1,size(atom_num_radius)
             item2=atom_num_radius(jj)
             IF(item1.ne.item2) THEN
              IF(dist(rx_fin(item1),rx_fin(item2),ry_fin(item1), &
                      ry_fin(item2),rz_fin(item1),rz_fin(item2),&
                      xlo_fin,xhi_fin,ylo_fin,yhi_fin,zlo_fin,&
                      zhi_fin).LT.rBN) THEN
                      cnt1=cnt1+1 
                      nn4_sorted_l1(atom_num_radius(ii),cnt1)=item2
         !  WRITE(6,*) cnt1,item1,nn4_sorted_l1(atom_num_radius(ii),cnt1)
         !       ! WRITE(6,*)'inner_inner_n4','ii',ii,'jj',jj,'kk',&
         !       !            kk,cnt,i1,i2,item3
              END IF
             END IF
              ! WRITE(6,*) 'inner_n4','i',ii,'j',jj,'k',kk,&
              !             nn4_sorted_l2(atom_num(k),:)
           END DO
         !    WRITE(6,*) item1, nn4_sorted_l1(atom_num_radius(ii),:)
         END DO
       
       DEALLOCATE(atom_num_radius)

!For i.eq 1 and 2<=j<=J_dim-1

         i=1
         DO j=2,J_dim-1
         
         ALLOCATE(atom_num_radius(24))
         i1=4*J_dim*(i-1)+4*(j-2)+1
         i2=4*J_dim*(i-1)+4*(j+1)
         i3=4*J_dim*i+4*(j-2)+1
         i4=4*J_dim*i+4*(j+1)
        
         !WRITE(6,*) i1,i2,i3,i4
         
         cnt=0
         DO ii=i1,i2
               cnt=cnt+1
               atom_num_radius(cnt)=ii
         END DO
         DO ii=i3,i4
               cnt=cnt+1
               atom_num_radius(cnt)=ii
         END DO

        !WRITE(6,*) i,j, atom_num_radius
              
        DO ii=5,8
           item1=atom_num_radius(ii)
           nn4_sorted_l1(atom_num_radius(ii),1)=item1
          ! WRITE(6,*) item1,nn4_sorted_l1(atom_num_radius(ii),1)
           cnt1=1
          DO jj=1,size(atom_num_radius)
             item2=atom_num_radius(jj)
             IF(item1.ne.item2) THEN
              IF(dist(rx_fin(item1),rx_fin(item2),ry_fin(item1), &
                      ry_fin(item2),rz_fin(item1),rz_fin(item2),&
                      xlo_fin,xhi_fin,ylo_fin,yhi_fin,zlo_fin,&
                      zhi_fin).LT.rBN) THEN
                      cnt1=cnt1+1 
                      nn4_sorted_l1(atom_num_radius(ii),cnt1)=item2
         !  WRITE(6,*) cnt1,item1,nn4_sorted_l1(atom_num_radius(ii),cnt1)
         !       ! WRITE(6,*)'inner_inner_n4','ii',ii,'jj',jj,'kk',&
         !       !            kk,cnt,i1,i2,item3
              END IF
             END IF
              ! WRITE(6,*) 'inner_n4','i',ii,'j',jj,'k',kk,&
              !             nn4_sorted_l2(atom_num(k),:)
           END DO
         !    WRITE(6,*) item1, nn4_sorted_l1(atom_num_radius(ii),:)
         END DO

         DEALLOCATE(atom_num_radius)

         END DO

!For j.eq 1 and 2<=i<=I_dim-1

         j=1
         DO i=2,I_dim-1

         ALLOCATE(atom_num_radius(24))
         i1=4*J_dim*(i-2)+4*(j-1)+1
         i2=4*J_dim*(i-2)+4*(j+1)
         i3=4*J_dim*(i-1)+4*(j-1)+1
         i4=4*J_dim*(i-1)+4*(j+1)
         i5=4*J_dim*i+4*(j-1)+1
         i6=4*J_dim*i+4*(j+1)
         
         !WRITE(6,*) i1,i2,i3,i4,i5,i6
         
         cnt=0
         DO ii=i1,i2
               cnt=cnt+1
               atom_num_radius(cnt)=ii
         END DO
         DO ii=i3,i4
               cnt=cnt+1
               atom_num_radius(cnt)=ii
         END DO
         DO ii=i5,i6
               cnt=cnt+1
               atom_num_radius(cnt)=ii
         END DO

        !WRITE(6,*) i,j, atom_num_radius
              
        DO ii=9,12
           item1=atom_num_radius(ii)
           nn4_sorted_l1(atom_num_radius(ii),1)=item1
          ! WRITE(6,*) item1,nn4_sorted_l1(atom_num_radius(ii),1)
           cnt1=1
          DO jj=1,size(atom_num_radius)
             item2=atom_num_radius(jj)
             IF(item1.ne.item2) THEN
              IF(dist(rx_fin(item1),rx_fin(item2),ry_fin(item1), &
                      ry_fin(item2),rz_fin(item1),rz_fin(item2),&
                      xlo_fin,xhi_fin,ylo_fin,yhi_fin,zlo_fin,&
                      zhi_fin).LT.rBN) THEN
                      cnt1=cnt1+1 
                      nn4_sorted_l1(atom_num_radius(ii),cnt1)=item2
         !  WRITE(6,*) cnt1,item1,nn4_sorted_l1(atom_num_radius(ii),cnt1)
         !       ! WRITE(6,*)'inner_inner_n4','ii',ii,'jj',jj,'kk',&
         !       !            kk,cnt,i1,i2,item3
              END IF
             END IF
              ! WRITE(6,*) 'inner_n4','i',ii,'j',jj,'k',kk,&
              !             nn4_sorted_l2(atom_num(k),:)
           END DO
         !    WRITE(6,*) item1, nn4_sorted_l1(atom_num_radius(ii),:)
         END DO

         DEALLOCATE(atom_num_radius)

         END DO


! For 2 <= i <= I_dim-1 and 2 <= j <= J_dim-1

        DO i=2,I_dim-1
         DO j=2,J_dim-1
             
         ALLOCATE(atom_num_radius(36))
         i1=4*J_dim*(i-2)+4*(j-2)+1
         i2=4*J_dim*(i-2)+4*(j+1)
         i3=4*J_dim*(i-1)+4*(j-2)+1
         i4=4*J_dim*(i-1)+4*(j+1)
         i5=4*J_dim*i+4*(j-2)+1
         i6=4*J_dim*i+4*(j+1)
         
       ! WRITE(6,*) i1,i2,i3,i4,i5,i6
        
         cnt=0
         DO ii=i1,i2
               cnt=cnt+1
               atom_num_radius(cnt)=ii
         END DO
         DO ii=i3,i4
               cnt=cnt+1
               atom_num_radius(cnt)=ii
         END DO
         DO ii=i5,i6
               cnt=cnt+1
               atom_num_radius(cnt)=ii
         END DO

        !WRITE(6,*) i,j, atom_num_radius
              
        DO ii=17,20
           item1=atom_num_radius(ii)
           nn4_sorted_l1(atom_num_radius(ii),1)=item1
          ! WRITE(6,*) item1,nn4_sorted_l1(atom_num_radius(ii),1)
           cnt1=1
          DO jj=1,size(atom_num_radius)
             item2=atom_num_radius(jj)
             IF(item1.ne.item2) THEN
              IF(dist(rx_fin(item1),rx_fin(item2),ry_fin(item1), &
                      ry_fin(item2),rz_fin(item1),rz_fin(item2),&
                      xlo_fin,xhi_fin,ylo_fin,yhi_fin,zlo_fin,&
                      zhi_fin).LT.rBN) THEN
                      cnt1=cnt1+1 
                      nn4_sorted_l1(atom_num_radius(ii),cnt1)=item2
         !  WRITE(6,*) cnt1,item1,nn4_sorted_l1(atom_num_radius(ii),cnt1)
         !       ! WRITE(6,*)'inner_inner_n4','ii',ii,'jj',jj,'kk',&
         !       !            kk,cnt,i1,i2,item3
              END IF
             END IF
              ! WRITE(6,*) 'inner_n4','i',ii,'j',jj,'k',kk,&
              !             nn4_sorted_l2(atom_num(k),:)
           END DO
         !    WRITE(6,*) item1, nn4_sorted_l1(atom_num_radius(ii),:)
         END DO

          DEALLOCATE(atom_num_radius)

         END DO
        END DO
                                
!For i=I_dim and j=1

         i=I_dim
         j=1
          ALLOCATE(atom_num_radius(16))
          i1=4*J_dim*(i-2)+4*(j-1)+1
          i2=4*J_dim*(i-2)+4*(j+1)
          i3=4*J_dim*(i-1)+4*(j-1)+1
          i4=4*J_dim*(i-1)+4*(j+1)

         !WRITE(6,*) i1,i2,i3,i4
        
          cnt=0
          DO ii=i1,i2
                cnt=cnt+1
                atom_num_radius(cnt)=ii
          END DO
          DO ii=i3,i4
                cnt=cnt+1
                atom_num_radius(cnt)=ii
          END DO

        !WRITE(6,*) i,j, atom_num_radius
              
        DO ii=9,12
           item1=atom_num_radius(ii)
           nn4_sorted_l1(atom_num_radius(ii),1)=item1
          ! WRITE(6,*) item1,nn4_sorted_l1(atom_num_radius(ii),1)
           cnt1=1
          DO jj=1,size(atom_num_radius)
             item2=atom_num_radius(jj)
             IF(item1.ne.item2) THEN
              IF(dist(rx_fin(item1),rx_fin(item2),ry_fin(item1), &
                      ry_fin(item2),rz_fin(item1),rz_fin(item2),&
                      xlo_fin,xhi_fin,ylo_fin,yhi_fin,zlo_fin,&
                      zhi_fin).LT.rBN) THEN
                      cnt1=cnt1+1 
                      nn4_sorted_l1(atom_num_radius(ii),cnt1)=item2
         !  WRITE(6,*) cnt1,item1,nn4_sorted_l1(atom_num_radius(ii),cnt1)
         !       ! WRITE(6,*)'inner_inner_n4','ii',ii,'jj',jj,'kk',&
         !       !            kk,cnt,i1,i2,item3
              END IF
             END IF
              ! WRITE(6,*) 'inner_n4','i',ii,'j',jj,'k',kk,&
              !             nn4_sorted_l2(atom_num(k),:)
           END DO
         !    WRITE(6,*) item1, nn4_sorted_l1(atom_num_radius(ii),:)
         END DO

          DEALLOCATE(atom_num_radius)


!For i=I_dim and 2<=j<=J_dim-1

         i=I_dim
         DO j=2,J_dim-1
          ALLOCATE(atom_num_radius(24))
          i1=4*J_dim*(i-2)+4*(j-2)+1
          i2=4*J_dim*(i-2)+4*(j+1)
          i3=4*J_dim*(i-1)+4*(j-2)+1
          i4=4*J_dim*(i-1)+4*(j+1)

         !WRITE(6,*) i1,i2,i3,i4
        
          cnt=0
          DO ii=i1,i2
                cnt=cnt+1
                atom_num_radius(cnt)=ii
          END DO
          DO ii=i3,i4
                cnt=cnt+1
                atom_num_radius(cnt)=ii
          END DO

        !WRITE(6,*) i,j, atom_num_radius
              
        DO ii=17,20
           item1=atom_num_radius(ii)
           nn4_sorted_l1(atom_num_radius(ii),1)=item1
          ! WRITE(6,*) item1,nn4_sorted_l1(atom_num_radius(ii),1)
           cnt1=1
          DO jj=1,size(atom_num_radius)
             item2=atom_num_radius(jj)
             IF(item1.ne.item2) THEN
              IF(dist(rx_fin(item1),rx_fin(item2),ry_fin(item1), &
                      ry_fin(item2),rz_fin(item1),rz_fin(item2),&
                      xlo_fin,xhi_fin,ylo_fin,yhi_fin,zlo_fin,&
                      zhi_fin).LT.rBN) THEN
                      cnt1=cnt1+1 
                      nn4_sorted_l1(atom_num_radius(ii),cnt1)=item2
         !  WRITE(6,*) cnt1,item1,nn4_sorted_l1(atom_num_radius(ii),cnt1)
         !       ! WRITE(6,*)'inner_inner_n4','ii',ii,'jj',jj,'kk',&
         !       !            kk,cnt,i1,i2,item3
              END IF
             END IF
              ! WRITE(6,*) 'inner_n4','i',ii,'j',jj,'k',kk,&
              !             nn4_sorted_l2(atom_num(k),:)
           END DO
         !    WRITE(6,*) item1, nn4_sorted_l1(atom_num_radius(ii),:)
         END DO

          DEALLOCATE(atom_num_radius)
         END DO

!For i=1 and j=J_dim

         j=J_dim
         i=1

          ALLOCATE(atom_num_radius(16))
          i1=4*J_dim*(i-1)+4*(j-2)+1
          i2=4*J_dim*(i-1)+4*j
          i3=4*J_dim*i+4*(j-2)+1
          i4=4*J_dim*i+4*j
          
          !WRITE(6,*) i1,i2,i3,i4

          cnt=0
          DO ii=i1,i2
                cnt=cnt+1
                atom_num_radius(cnt)=ii
          END DO
          DO ii=i3,i4
                cnt=cnt+1
                atom_num_radius(cnt)=ii
          END DO

        !WRITE(6,*) i,j, atom_num_radius
              
        DO ii=5,8
           item1=atom_num_radius(ii)
           nn4_sorted_l1(atom_num_radius(ii),1)=item1
          ! WRITE(6,*) item1,nn4_sorted_l1(atom_num_radius(ii),1)
           cnt1=1
          DO jj=1,size(atom_num_radius)
             item2=atom_num_radius(jj)
             IF(item1.ne.item2) THEN
              IF(dist(rx_fin(item1),rx_fin(item2),ry_fin(item1), &
                      ry_fin(item2),rz_fin(item1),rz_fin(item2),&
                      xlo_fin,xhi_fin,ylo_fin,yhi_fin,zlo_fin,&
                      zhi_fin).LT.rBN) THEN
                      cnt1=cnt1+1 
                      nn4_sorted_l1(atom_num_radius(ii),cnt1)=item2
         !  WRITE(6,*) cnt1,item1,nn4_sorted_l1(atom_num_radius(ii),cnt1)
         !       ! WRITE(6,*)'inner_inner_n4','ii',ii,'jj',jj,'kk',&
         !       !            kk,cnt,i1,i2,item3
              END IF
             END IF
              ! WRITE(6,*) 'inner_n4','i',ii,'j',jj,'k',kk,&
              !             nn4_sorted_l2(atom_num(k),:)
           END DO
         !    WRITE(6,*) item1, nn4_sorted_l1(atom_num_radius(ii),:)
         END DO

          DEALLOCATE(atom_num_radius)

!For 2<=i<=I_dim-1 and j=J_dim

         j=J_dim
         DO i=2,I_dim-1

          ALLOCATE(atom_num_radius(24))
          i1=4*J_dim*(i-2)+4*(j-2)+1
          i2=4*J_dim*(i-2)+4*j
          i3=4*J_dim*(i-1)+4*(j-2)+1
          i4=4*J_dim*(i-1)+4*j
          i5=4*J_dim*i+4*(j-2)+1
          i6=4*J_dim*i+4*j
          
          !WRITE(6,*) i1,i2,i3,i4,i5,i6

          cnt=0
          DO ii=i1,i2
                cnt=cnt+1
                atom_num_radius(cnt)=ii
          END DO
          DO ii=i3,i4
                cnt=cnt+1
                atom_num_radius(cnt)=ii
          END DO
          DO ii=i5,i6
                cnt=cnt+1
                atom_num_radius(cnt)=ii
          END DO

        !WRITE(6,*) i,j, atom_num_radius
              
        DO ii=13,16
           item1=atom_num_radius(ii)
           nn4_sorted_l1(atom_num_radius(ii),1)=item1
          ! WRITE(6,*) item1,nn4_sorted_l1(atom_num_radius(ii),1)
           cnt1=1
          DO jj=1,size(atom_num_radius)
             item2=atom_num_radius(jj)
             IF(item1.ne.item2) THEN
              IF(dist(rx_fin(item1),rx_fin(item2),ry_fin(item1), &
                      ry_fin(item2),rz_fin(item1),rz_fin(item2),&
                      xlo_fin,xhi_fin,ylo_fin,yhi_fin,zlo_fin,&
                      zhi_fin).LT.rBN) THEN
                      cnt1=cnt1+1 
                      nn4_sorted_l1(atom_num_radius(ii),cnt1)=item2
         !  WRITE(6,*) cnt1,item1,nn4_sorted_l1(atom_num_radius(ii),cnt1)
         !       ! WRITE(6,*)'inner_inner_n4','ii',ii,'jj',jj,'kk',&
         !       !            kk,cnt,i1,i2,item3
              END IF
             END IF
              ! WRITE(6,*) 'inner_n4','i',ii,'j',jj,'k',kk,&
              !             nn4_sorted_l2(atom_num(k),:)
           END DO
         !    WRITE(6,*) item1, nn4_sorted_l1(atom_num_radius(ii),:)
         END DO

          DEALLOCATE(atom_num_radius)
         END DO


!For i=I_dim and j=J_dim
           
          i=I_dim
          j=J_dim
      
          ALLOCATE(atom_num_radius(16))
          i1=4*J_dim*(i-2)+4*(j-2)+1
          i2=4*J_dim*(i-2)+4*j
          i3=4*J_dim*(i-1)+4*(j-2)+1
          i4=4*J_dim*(i-1)+4*j
         
          !WRITE(6,*) i1,i2,i3,i4

          cnt=0
          DO ii=i1,i2
                cnt=cnt+1
                atom_num_radius(cnt)=ii
          END DO
          DO ii=i3,i4
                cnt=cnt+1
                atom_num_radius(cnt)=ii
          END DO
      
        !WRITE(6,*) i,j, atom_num_radius
              
        DO ii=13,16
           item1=atom_num_radius(ii)
           nn4_sorted_l1(atom_num_radius(ii),1)=item1
          ! WRITE(6,*) item1,nn4_sorted_l1(atom_num_radius(ii),1)
           cnt1=1
          DO jj=1,size(atom_num_radius)
             item2=atom_num_radius(jj)
             IF(item1.ne.item2) THEN
              IF(dist(rx_fin(item1),rx_fin(item2),ry_fin(item1), &
                      ry_fin(item2),rz_fin(item1),rz_fin(item2),&
                      xlo_fin,xhi_fin,ylo_fin,yhi_fin,zlo_fin,&
                      zhi_fin).LT.rBN) THEN
                      cnt1=cnt1+1 
                      nn4_sorted_l1(atom_num_radius(ii),cnt1)=item2
         !  WRITE(6,*) cnt1,item1,nn4_sorted_l1(atom_num_radius(ii),cnt1)
         !       ! WRITE(6,*)'inner_inner_n4','ii',ii,'jj',jj,'kk',&
         !       !            kk,cnt,i1,i2,item3
              END IF
             END IF
              ! WRITE(6,*) 'inner_n4','i',ii,'j',jj,'k',kk,&
              !             nn4_sorted_l2(atom_num(k),:)
           END DO
         !    WRITE(6,*) item1, nn4_sorted_l1(atom_num_radius(ii),:)
         END DO

          DEALLOCATE(atom_num_radius)

           DO i=1,natoms_fin/2  

               WRITE(6,*)'nearest_neighbor_l1', i, nn4_sorted_l1(i,:)
             
           END DO 


!        DO i=1,natoms_fin/2
!           item1=sorted_index_l1(i)
!           nn4_sorted_l1(i,1)=item1
!           
!           cnt=0
!           DO j=1,natoms_fin/2
!              IF(i.ne.j) THEN
!                item2=sorted_index_l1(j)
!                IF(dist(rx_sorted_l1(item2),rx_sorted_l1(item1),&
!                   ry_sorted_l1(item2),ry_sorted_l1(item1),&
!                   rz_sorted_l1(item2),rz_sorted_l1(item1),&
!                   xlo_fin,xhi_fin,ylo_fin,yhi_fin,zlo_fin,zhi_fin)&
!                   .LE.rBN)THEN
!                           cnt=cnt+1
!                           nn4_sorted_l1(i,cnt+1)=item2
!                          ! write(6,*)'nn4_inner',cnt,item1,item2
!                END IF
!              END IF
!           END DO
!           WRITE(6,*) 'nn4_l1',i,nn4_sorted_l1(i,:)
!        END DO
!        
        
        CALL cpu_time(finish)
    
        WRITE(6,*) '("Time taken for nearest_neighbor_l1 = ",f6.3,&
                                    " seconds.")',finish-start
    
        
        END SUBROUTINE nearest_neighbor_l1


        SUBROUTINE Find_layer_l2

        USE global_variable
        USE my_functions
        USE DynamicalArrays

        IMPLICIT NONE
      
        INTEGER*4::i,j,k
        REAL*8 :: dx,dy,dz,dr 
        REAL*8 :: start, finish

        INTEGER*4::nlines,io,cnt,cnt1,I_dim,J_dim,dummy,dummy1,dummy2,&
             dummy3,dummy4,dummy5,dummy6,dummy7,dummy8,item1,item2,item
        REAL*8::alpha,item_r(3),item_period_1(3),item_period_2(3),&
                item_period_3(3),item_period_4(3),item_period_5(3),&
                item_period_6(3),item_period_7(3),item_period_8(3)
        REAL*8,DIMENSION(:),ALLOCATABLE:: rx_dummy,ry_dummy,rz_dummy
        CHARACTER::item_b        


        INTEGER*4,DIMENSION(:),ALLOCATABLE::indices
        CHARACTER(len=9),DIMENSION(:),ALLOCATABLE::box_pos
        CHARACTER,DIMENSION(:),ALLOCATABLE::b_dummy

        CALL cpu_time(start)

!        nlines = 0
!        OPEN (1, file = 'lay2_info_1p0_AA.xyz')
!        DO
!        READ(1,*,iostat=io)
!        IF (io/=0) EXIT
!        nlines = nlines + 1
!        END DO
!
!        WRITE(6,*) nlines
!
        ALLOCATE(layer2_init(natoms_init/2),layer2_fin(natoms_fin/2))!, &
!                indices(nlines),box_pos(nlines)) 
!
!        REWIND(1)
!
!        DO i=1,nlines
!         READ(1,*) indices(i)
!         WRITE(6,*) indices(i)
!        END DO  
!
!        REWIND(1)
!
!        DO i=1,nlines
!         READ(1,110) box_pos(i)
!     110 FORMAT(13X, A9)  
!         WRITE(6,*) box_pos(i)
!        END DO  
!
!        CLOSE(1)
!
!
!        DO i=1,natoms_init/2
!              
!             k=i+natoms_init/2             
!  
!             layer2_init(i)=atom_index_init(k)
!        
!                DO j=1,nlines
!
!                  IF(layer2_init(i)==indices(j)) THEN
!                     IF(box_pos(j)=='neg_box_x') THEN
!                       rx_init(layer2_init(i))=rx_init(layer2_init(i))-&
!                                               box_init(1)*cos(1.0d0)
!                       write(6,*) 'i-term, -x',j 
!                     ELSEIF(box_pos(j)=='pos_box_x') THEN
!                       rx_init(layer2_init(i))=rx_init(layer2_init(i))+&
!                                               box_init(1) 
!                       write(6,*) 'i-term, +x',j 
!                     ELSEIF(box_pos(j)=='neg_box_y') THEN
!                       ry_init(layer2_init(i))=ry_init(layer2_init(i))-&
!                                               box_init(2) 
!                       write(6,*) 'i-term, -y',j 
!                     ELSEIF(box_pos(j)=='pos_box_y') THEN
!                       ry_init(layer2_init(i))=ry_init(layer2_init(i))+&
!                                               box_init(2) 
!                       write(6,*) 'i-term, +y',j 
!                     END IF    
!                  END IF            
!
!                END DO
!
!         
!             IF(i.ge.1.and.i.le.x_cell*4/2) THEN
!                IF(rx_init(layer2_init(i))>(xhi_init-1.d0))THEN
!                   rx_init(layer2_init(i))=ABS(rx_init(layer2_init(i)) &
!                                      - box_init(1))
!                   WRITE(6,*)'i-term',i
!                END IF
!             END IF
!         
!             IF(i.ge.1.and.i.le.y_cell*4/2) THEN
!                IF(ry_init(layer2_init(i))>(yhi_init-1.d0))THEN
!                   ry_init(layer2_init(i))=ABS(ry_init(layer2_init(i)) &
!                                      - box_init(2))
!                   WRITE(6,*)'i-term',i
!                END IF
!             END IF
!              
!             WRITE(6,*) 'Layer2_initial',layer2_init(i),&
!             rx_init(layer2_init(i)),ry_init(layer2_init(i)),&
!             rz_init(layer2_init(i)),xhi_init,yhi_init,zhi_init
!
!        END DO
!           
!        OPEN(UNIT=layer2_output,FILE=file_layer2)
!
!        DO i=1,natoms_fin/2
!              
!             k=i+natoms_fin/2             
! 
!             layer2_fin(i)=atom_index_fin(k)
!        
!                DO j=1,nlines
!
!                  IF(layer2_fin(i)==indices(j)) THEN
!                     IF(box_pos(j)=='neg_box_x') THEN
!                       rx_fin(layer2_fin(i))=rx_fin(layer2_fin(i))-&
!                                               box_fin(1) 
!                       write(6,*) 'i-term, -x',j 
!                     ELSEIF(box_pos(j)=='pos_box_x') THEN
!                       rx_fin(layer2_fin(i))=rx_fin(layer2_fin(i))+&
!                                               box_fin(1) 
!                       write(6,*) 'i-term, +x',j 
!                     ELSEIF(box_pos(j)=='neg_box_y') THEN
!                       ry_fin(layer2_fin(i))=ry_fin(layer2_fin(i))-&
!                                               box_fin(2) 
!                       write(6,*) 'i-term, -y',j 
!                     ELSEIF(box_pos(j)=='pos_box_y') THEN
!                       ry_fin(layer2_fin(i))=ry_fin(layer2_fin(i))+&
!                                               box_fin(2) 
!                       write(6,*) 'i-term, +y',j 
!                     END IF    
!                  END IF            
!
!                END DO
!
!             IF(i.ge.1.and.i.le.x_cell*4/2) THEN
!                IF(rx_fin(layer2_fin(i))>(xhi_fin-1.d0))THEN
!                   rx_fin(layer2_fin(i))=ABS(rx_fin(layer2_fin(i)) &
!                                      - box_fin(1))
!                   WRITE(6,*)'i-term',i
!                END IF
!             END IF
!         
!             IF(i.ge.1.and.i.le.y_cell*4/2) THEN
!                IF(ry_fin(layer2_fin(i))>(yhi_fin-1.d0))THEN
!                   ry_fin(layer2_fin(i))=ABS(ry_fin(layer2_fin(i)) &
!                                      - box_fin(2))
!                   WRITE(6,*)'i-term',i
!                END IF
!             END IF

        DO i=1,natoms_init/2

             k=i+natoms_init/2             
  
             layer2_init(i)=atom_index_init(k)
 
       !      WRITE(6,*) i, k, atom_index_init(k), layer2_init(i)
               
             IF(rx_init(layer2_init(i))<xlo_init) THEN
                rx_init(layer2_init(i))=rx_init(layer2_init(i)) + &     
                                       box_init(1)
             ELSEIF(rx_init(layer2_init(i))>xhi_init) THEN
                rx_init(layer2_init(i))=rx_init(layer2_init(i)) - &
                                       box_init(1)
             END IF
         
                
             IF(ry_init(layer2_init(i))<ylo_init) THEN
                ry_init(layer2_init(i))=ry_init(layer2_init(i)) + &     
                                       box_init(2)
             ELSEIF(ry_init(layer2_init(i))>yhi_init) THEN
                ry_init(layer2_init(i))=ry_init(layer2_init(i)) - &
                                       box_init(2)
             END IF
         
                
             IF(rz_init(layer2_init(i))<zlo_init) THEN
                rz_init(layer2_init(i))=rz_init(layer2_init(i)) + &     
                                       box_init(3)
             ELSEIF(rz_init(layer2_init(i))>zhi_init) THEN
                rz_init(layer2_init(i))=rz_init(layer2_init(i)) - &
                                       box_init(3)
             END IF
         
             WRITE(6,*) 'Layer2_initial',layer2_init(i),&
             rx_init(layer2_init(i)),ry_init(layer2_init(i)),&
             rz_init(layer2_init(i)),xhi_init,yhi_init,zhi_init
            
        END DO
 
        
        OPEN(UNIT=layer1_output,FILE=file_layer1)

        DO i=1,natoms_fin/2
        
             k=i+natoms_fin/2             
 
             layer2_fin(i)=atom_index_fin(k)
 
       !      WRITE(6,*) i, k, atom_index_fin(k), layer2_fin(i)
             
             IF(rx_fin(layer2_fin(i))<xlo_fin) THEN
                rx_fin(layer2_fin(i))=rx_fin(layer2_fin(i)) + &     
                                       box_fin(1)
             ELSEIF(rx_fin(layer2_fin(i))>xhi_fin) THEN
                rx_fin(layer2_fin(i))=rx_fin(layer2_fin(i)) - &
                                       box_fin(1)
             END IF
         
                
             IF(ry_fin(layer2_fin(i))<ylo_fin) THEN
                ry_fin(layer2_fin(i))=ry_fin(layer2_fin(i)) + &     
                                       box_fin(2)
             ELSEIF(ry_fin(layer2_fin(i))>yhi_fin) THEN
                ry_fin(layer2_fin(i))=ry_fin(layer2_fin(i)) - &
                                       box_fin(2)
             END IF
         
                
             IF(rz_fin(layer2_fin(i))<zlo_fin) THEN
                rz_fin(layer2_fin(i))=rz_fin(layer2_fin(i)) + &     
                                       box_fin(3)
             ELSEIF(rz_fin(layer2_fin(i))>zhi_fin) THEN
                rz_fin(layer2_fin(i))=rz_fin(layer2_fin(i)) - &
                                       box_fin(3)
             END IF
                
             WRITE(6,*) 'Layer2_final',layer2_fin(i),&
             rx_fin(layer2_fin(i)),ry_fin(layer2_fin(i)),&
             rz_fin(layer2_fin(i)),xhi_fin,yhi_fin,zhi_fin
          
             WRITE(layer2_output,*)layer2_fin(i),&
             rx_fin(layer2_fin(i)),ry_fin(layer2_fin(i)),&
             rz_fin(layer2_fin(i))


        END DO

        CLOSE(layer2_output)

!# Similarity check

        DO i=1,natoms_fin/2
         
             dx=rx_fin(layer2_fin(i))-rx_init(layer2_init(i))
             dy=ry_fin(layer2_fin(i))-ry_init(layer2_init(i))
             dz=rz_fin(layer2_fin(i))-rz_init(layer2_init(i))
             dr=SQRT(dx*dx+dy*dy+dz*dz)

             WRITE(6,*) 'Lay2_dist_comp',layer2_fin(i),layer2_init(i), &
                                      dx,dy,dz,dr
        
             IF (dr > 1.0d0) THEN

            WRITE(6,*) 'The atom',i,'at the 2nd layer for initial',&
                          'frame is not associated with final frame '

            WRITE(layer2_output,*)'The atom',i,'at the 2nd layer for',&
                     'initial frame is not associated with final frame'
             END IF
             
                 
         END DO 

        DO i=1,natoms_fin/2
         
             dx=rx_fin(layer2_fin(i))-rx_init(layer1_fin(i))
             dy=ry_fin(layer2_fin(i))-ry_init(layer1_fin(i))
             dz=rz_fin(layer2_fin(i))-rz_init(layer1_fin(i))
             dr=SQRT(dx*dx+dy*dy)

             WRITE(6,*) 'Lay1_Lay2_dist_comp',layer2_fin(i), &
                                              layer1_fin(i), &
                                              dx,dy,dz,dr
        
         END DO 

        CALL cpu_time(finish)
    
        WRITE(6,*) '("Time taken for find_layer_l2 = ",f6.3,&
                                    " seconds.")',finish-start
    
        END SUBROUTINE Find_layer_l2

        SUBROUTINE Find_layer_l1

        USE global_variable
        USE my_functions

      
        INTEGER*4 :: i,k,sum_r_l1
        REAL*8 :: dx,dy,dz,dr
        REAL*8 :: start, finish

        CALL cpu_time(start)

        ALLOCATE(layer1_init(natoms_init/2),layer1_fin(natoms_fin/2))

        DO i=1,natoms_init/2

             layer1_init(i)=atom_index_init(i)
                
             IF(rx_init(layer1_init(i))<xlo_init) THEN
                rx_init(layer1_init(i))=rx_init(layer1_init(i)) + &     
                                       box_init(1)
             ELSEIF(rx_init(layer1_init(i))>xhi_init) THEN
                rx_init(layer1_init(i))=rx_init(layer1_init(i)) - &
                                       box_init(1)
             END IF
         
                
             IF(ry_init(layer1_init(i))<ylo_init) THEN
                ry_init(layer1_init(i))=ry_init(layer1_init(i)) + &     
                                       box_init(2)
             ELSEIF(ry_init(layer1_init(i))>yhi_init) THEN
                ry_init(layer1_init(i))=ry_init(layer1_init(i)) - &
                                       box_init(2)
             END IF
         
                
             IF(rz_init(layer1_init(i))<zlo_init) THEN
                rz_init(layer1_init(i))=rz_init(layer1_init(i)) + &     
                                       box_init(3)
             ELSEIF(rz_init(layer1_init(i))>zhi_init) THEN
                rz_init(layer1_init(i))=rz_init(layer1_init(i)) - &
                                       box_init(3)
             END IF
         
             WRITE(6,*) 'Layer1_initial',layer1_init(i),&
             rx_init(layer1_init(i)),ry_init(layer1_init(i)),&
             rz_init(layer1_init(i)),xhi_init,yhi_init,zhi_init
            
        END DO
 
        
        OPEN(UNIT=layer1_output,FILE=file_layer1)

        DO i=1,natoms_fin/2
        
             layer1_fin(i)=atom_index_fin(i)
             
             IF(rx_fin(layer1_fin(i))<xlo_fin) THEN
                rx_fin(layer1_fin(i))=rx_fin(layer1_fin(i)) + &     
                                       box_fin(1)
             ELSEIF(rx_fin(layer1_fin(i))>xhi_fin) THEN
                rx_fin(layer1_fin(i))=rx_fin(layer1_fin(i)) - &
                                       box_fin(1)
             END IF
         
                
             IF(ry_fin(layer1_fin(i))<ylo_fin) THEN
                ry_fin(layer1_fin(i))=ry_fin(layer1_fin(i)) + &     
                                       box_fin(2)
             ELSEIF(ry_fin(layer1_fin(i))>yhi_fin) THEN
                ry_fin(layer1_fin(i))=ry_fin(layer1_fin(i)) - &
                                       box_fin(2)
             END IF
         
                
             IF(rz_fin(layer1_fin(i))<zlo_fin) THEN
                rz_fin(layer1_fin(i))=rz_fin(layer1_fin(i)) + &     
                                       box_fin(3)
             ELSEIF(rz_fin(layer1_fin(i))>zhi_fin) THEN
                rz_fin(layer1_fin(i))=rz_fin(layer1_fin(i)) - &
                                       box_fin(3)
             END IF
         
             WRITE(6,*) 'Layer1_final',layer1_fin(i),&
             rx_fin(layer1_fin(i)),ry_fin(layer1_fin(i)),&
             rz_fin(layer1_fin(i)),xhi_fin,yhi_fin,zhi_fin
             
             WRITE(layer1_output,*)layer1_fin(i),&
                     rx_fin(layer1_fin(i)),ry_fin(layer1_fin(i)),&
                     rz_fin(layer1_fin(i))
                 
        END DO

        CLOSE(layer1_output)

         sum_r_l1=0.d0
         DO i=1,natoms_fin/2-1
               dr=dist(rx_fin(layer1_fin(i+1)),rx_fin(layer1_fin(i)),&
                       ry_fin(layer1_fin(i+1)),ry_fin(layer1_fin(i)),&
                       rz_fin(layer1_fin(i+1)),rz_fin(layer1_fin(i)),&
                       xlo_fin,xhi_fin,ylo_fin,yhi_fin,zlo_fin,zhi_fin)
               sum_r_l1=sum_r_l1+SQRT(dr)
         END DO
         rBN_l1=sum_r_l1/FLOAT((natoms_fin/2)-1)
         write(6,*) sum_r_l1, rBN_l1

! Similarity check

!        DO i=1,natoms_fin/2
!         
!             dx=rx_fin(layer1_fin(i))-rx_init(layer1_init(i))
!             dy=ry_fin(layer1_fin(i))-ry_init(layer1_init(i))
!             dz=rz_fin(layer1_fin(i))-rz_init(layer1_init(i))
!             dr=SQRT(dx*dx+dy*dy+dz*dz)
!
!             WRITE(6,*) 'Lay1_dist_comp',layer1_fin(i),layer1_init(i), &
!                                      dx,dy,dz,dr
!        
!             IF (dr > 1.0d0) THEN
!
!               WRITE(6,*) 'The atom',i,'at the 1st layer for initial',&
!                          'frame is not associated with final frame '
!             END IF
!             
!                 
!         END DO 


        CALL cpu_time(finish)
    
        WRITE(6,*) '("Time taken for find_layer_l1 = ",f6.3,&
                                    " seconds.")',finish-start
    

        END SUBROUTINE Find_layer_l1

        SUBROUTINE ReadStructure_xyz

        USE global_variable
       
        IMPLICIT NONE
        INTEGER*4::i,lskip,unit_input_1,unit_input_2
        CHARACTER:: a
        CHARACTER(len=100)::file_init,file_fin
        REAL*8 :: start, finish

        CALL cpu_time(start)
        
        WRITE(6,*) 'Tell us about xmin,xmax,ymin,ymax,zmin,zmax'
        WRITE(6,*) 'Initial File'
        READ(5,*) xlo_init,xhi_init,ylo_init,yhi_init,zlo_init,zhi_init
        WRITE(6,*) xlo_init,xhi_init,ylo_init,yhi_init,zlo_init,zhi_init
        WRITE(6,*) 'Final File'
        READ(5,*) xlo_fin,xhi_fin,ylo_fin,yhi_fin,zlo_fin,zhi_fin
        WRITE(6,*)xlo_fin,xhi_fin,ylo_fin,yhi_fin,zlo_fin,zhi_fin
        WRITE(6,*) 'Dimension of the system'
        READ(5,*) x_cell       
        READ(5,*) y_cell      
        WRITE(6,*) 'x_cell', x_cell, 'y_cell', y_cell 
        READ(5,*) file_init
        READ(5,*) file_fin
        WRITE(6,*) 'initial file',"\t", file_init,"\n", &
                       'final file',"\t",file_fin
        
        READ(5,*) order_param_file
        READ(5,*) hexagon_data_dist_file
        READ(5,*) hexagon_data_up_xyz_file
        READ(5,*) hexagon_data_down_xyz_file
        READ(5,*) file_layer1 
        READ(5,*) file_layer2 
        READ(5,*) hexagon_rcom_data_file
 
        WRITE(6,*)"Look for output file 1","\t", order_param_file,"\n",&
                 'output file2',"\t",hexagon_data_dist_file,"\n",&
                 'output file3',"\t",hexagon_data_up_xyz_file,"\n",&
                 'output file4',"\t",hexagon_data_down_xyz_file,"\n",&
                 'output_file5',"\t",file_layer1,"\n",'output_file6',&
                  file_layer2,"\n","output_file6",hexagon_rcom_data_file
        READ(5,*)unit_input_1,unit_input_2,unit_output_1,unit_output_2,&
                 unit_output_3,unit_output_4,layer1_output,&
                 layer2_output,rcom_output
        WRITE(6,*) unit_input_1,unit_input_2,unit_output_1,&
                   unit_output_2,unit_output_3,unit_output_4,&
                   layer1_output,layer2_output,rcom_output
       
        box_orig(1)=429.462d0
        box_orig(2)=247.95d0
        box_orig(3)=100.0d0
 
        OPEN(UNIT=unit_input_1,FILE=file_init)
        OPEN(UNIT=unit_input_2,FILE=file_fin)
        
        READ(unit_input_1,*)natoms_init
        READ(unit_input_1,10)timestep_init
     10 FORMAT(17X, I10)  
        WRITE(6,*)'Initial Step',natoms_init,timestep_init
        ALLOCATE(b_init(natoms_init), rx_init(natoms_init), &
                 ry_init(natoms_init), rz_init(natoms_init), &
                 atom_index_init(natoms_init))
       
        box_init(1)=ABS(xhi_init-xlo_init)
        box_init(2)=ABS(yhi_init-ylo_init)
        box_init(3)=ABS(zhi_init-zlo_init)
 
        DO i=1,natoms_init
          READ(unit_input_1,*)b_init(i),rx_init(i),ry_init(i),rz_init(i)
        ENDDO

        DO i=1,natoms_init
          atom_index_init(i)=i
          WRITE(6,*)b_init(i),atom_index_init(i),rx_init(i),ry_init(i),&
                    rz_init(i)
        END DO

        CLOSE(unit_input_1)

        READ(unit_input_2,*)natoms_fin
        READ(unit_input_2,20)timestep_fin
     20 FORMAT(17X, I10)  
        WRITE(6,*)'Final Step', natoms_fin,timestep_fin
        ALLOCATE(b_fin(natoms_fin), rx_fin(natoms_fin),&
                 ry_fin(natoms_fin), rz_fin(natoms_fin), &
                 atom_index_fin(natoms_fin))
        
        box_fin(1)=ABS(xhi_fin-xlo_fin)
        box_fin(2)=ABS(yhi_fin-ylo_fin)
        box_fin(3)=ABS(zhi_fin-zlo_fin)

        DO i=1,natoms_fin
          READ(unit_input_2,*)b_fin(i),rx_fin(i),ry_fin(i),rz_fin(i)
        ENDDO

        DO i=1,natoms_fin
          atom_index_fin(i)=i
          WRITE(6,*)b_fin(i),atom_index_fin(i),rx_fin(i),ry_fin(i),&
                    rz_fin(i)
        END DO
    
        CLOSE(unit_input_2)
        
        CALL cpu_time(finish)
    
        WRITE(6,*) '("Time taken for read_structure = ",f6.3,&
                                    " seconds.")',finish-start
    
        END SUBROUTINE ReadStructure_xyz
 

        SUBROUTINE find_min(array,array_dim,min_val,min_pos)

        IMPLICIT NONE

        INTEGER*4::i
        
        INTEGER*4,INTENT(IN)::array_dim
        INTEGER*4,INTENT(OUT)::min_pos  
             
        REAL*8,INTENT(IN)::array(array_dim)
        REAL*8,INTENT(OUT)::min_val

        min_val=array(1)
        min_pos=1

         DO i=2,array_dim
            IF(array(i)<min_val) THEN
                min_val=array(i)
                min_pos=i
            END IF
         END DO

        END SUBROUTINE find_min

        SUBROUTINE heap_sort_nr(nterms,array)

        IMPLICIT NONE        

        INTEGER*4,INTENT(IN)::nterms
        INTEGER*4,INTENT(INOUT)::array(nterms)
       
        INTEGER*4::i,last,temp
        
        last=SIZE(array)

        !Building the heap

        DO i=last/2,1,-1
            CALL heapify_nr(nterms,array,i)
        END DO
      
        ! Unpick the heap

        DO i=last,2,-1
             temp=array(1)
             array(1)=array(i)
             array(i)=temp
             CALL heapify_nr(nterms,array(1:i-1),1)
        END DO

        END SUBROUTINE

        SUBROUTINE heapify_nr(nterms,array,i)
   
        IMPLICIT NONE
        
        INTEGER*4,INTENT(IN)::nterms,i
        INTEGER*4,INTENT(INOUT)::array(nterms)

        INTEGER*4:: left,right,root,last,largest,root_val

        last=size(array)
        root=i
        left=2*root
        root_val=array(root)
        largest=root

        DO WHILE(left.LE.last)
             right=left+1
            
             IF(left.LE.last)THEN
                IF(array(left).gt.array(largest)) largest=left
             END IF
  
             IF (right.LE.last) THEN
                IF(array(right).gt.array(largest)) largest=right
             END IF

             IF (largest.EQ.root) EXIT

             array(root)=array(largest)
             array(largest)=root_val
             root=largest
             left=2*root
        END DO

        !    array(largest)=root_val

        END SUBROUTINE heapify_nr
             

        SUBROUTINE Normal(v1,v2,v3,v4,vnorm)

        USE global_variable
        USE my_functions

        IMPLICIT NONE

        INTEGER*4:: i,j
       
        REAL*8,INTENT(IN)::v1(3),v2(3),v3(3),v4(3)

        REAL*8,INTENT(OUT)::vnorm(3)
   
        REAL*8::norm_dist

        vnorm(:)=cross(v2(1),v3(1),v4(1),v2(2),v3(2),v4(2),&
                        v2(3),v3(3),v4(3))
        !WRITE(6,*) 'norm', vnorm

        IF(vnorm(3).lt.0.d0) THEN 
              vnorm(:)=cross(v2(1),v4(1),v3(1),v2(2),v4(2),v3(2),&
                        v2(3),v4(3),v3(3))        
        !        WRITE(6,*) 'inverted_norm', vnorm
        END IF

        !norm_dist=dist(vnorm(1),0.0,vnorm(2),0.0,vnorm(3),0.0,&
        !          xhi_fin,xlo_fin,yhi_fin,ylo_fin,zhi_fin,zlo_fin)

        norm_dist=SQRT(DOT_PRODUCT(vnorm,vnorm))

        vnorm(:)=vnorm(:)/norm_dist

        !WRITE(6,*) norm_dist,vnorm

        END SUBROUTINE Normal

        SUBROUTINE inside_cone(x,norm,h,r,p,angle,linsidecone)

        USE global_variable
        USE my_functions

        IMPLICIT NONE
        
        REAL*8,INTENT(IN)::x(3),norm(3),h,r,p(3)

        REAL*8,INTENT(OUT)::angle
        LOGICAL,INTENT(OUT)::linsidecone

        REAL*8::point_z_dist,dist_vec(3),cone_slope,cone_limit,&
              ortho_dist_vec(3),ortho_dist2

        dist_vec(:)=p(:)-x(:)

        point_z_dist=SQRT(DOT_PRODUCT(dist_vec,norm))
     
        angle=theta(dist_vec,norm)
        
        IF(point_z_dist.lt.0.d0) THEN
                WRITE(6,*) 'The point p is below the apex'
        ELSEIF(point_z_dist.gt.h) THEN
                WRITE(6,*) 'The point p is above the cone, increase h'
        !ELSE
        !        WRITE(6,*) 'The point p is within the cone-height'
        END IF
 
        !angle=2.d0*datan(r/h)        

        cone_slope=r*r/h

!        cone_slope=r*dtan(angle/2.d0)

        cone_limit=cone_slope*point_z_dist
         
        ortho_dist_vec(:)=dist_vec(:)-point_z_dist*norm(:)
        
        ortho_dist2=DOT_PRODUCT(ortho_dist_vec,ortho_dist_vec)

        !WRITE(6,*) 'cone_info', dist_vec, point_z_dist, cone_slope, &
        !           cone_limit, ortho_dist_vec, ortho_dist2

        IF (ortho_dist2.le.cone_limit) THEN
           linsidecone=.TRUE.
         !  WRITE(6,*) 'The point is inside the cone',linsidecone,& 
         !             cone_limit,ortho_dist2
        ELSE
           linsidecone=.FALSE.
         !  WRITE(6,*) 'The point is outside the cone',linsidecone,&
         !             cone_limit,ortho_dist2
        END IF 

        END SUBROUTINE inside_cone


      SUBROUTINE find_upper_layer_atom(block,r_hexa_l1,&
                         norm_hexa_l1,cone_height,cone_radius,atom_num,&
                         atom_pos,atom_type,angle,dist_plane,lnohexagon)

      USE global_variable
      USE my_functions

      IMPLICIT NONE
 
      INTEGER*4,INTENT(IN)::block(36) 
      REAL*8,INTENT(IN)::r_hexa_l1(3),norm_hexa_l1(3),cone_height,&
                         cone_radius
      INTEGER*4,INTENT(OUT)::atom_num
      REAL*8,INTENT(OUT)::angle,dist_plane,atom_pos(3)
      CHARACTER,INTENT(OUT)::atom_type
      LOGICAL,INTENT(OUT)::lnohexagon
 
      INTEGER*4:: i,item,k
      REAL*8::dist_pt(3),r_hexa_l2(3),norm_hexa_tr_l1(3),angle_cone,&
              r_hexa_tr_l1(3),r_hexa_tr_l2(3),atom_tr_pos(3),nx,ny,nz
      LOGICAL::ptinsidecone
 
      LOGICAL,DIMENSION(:),ALLOCATABLE::lupperatom(:),l1(:)

        k=natoms_fin/2 

        ALLOCATE(lupperatom(36),l1(36))

        !WRITE(6,*) block

        DO i=1,36
        
           item=abs(block(i))+k           

           l1(i)=.FALSE.        

           lupperatom(i)=.FALSE.
           
        IF((norm_hexa_l1(1).eq.0.d0).AND.&
           (norm_hexa_l1(2).eq.0.d0).AND.&
           (norm_hexa_l1(3).eq.1.d0))THEN
           
  
            r_hexa_l2(1)=rx_fin(item)
            r_hexa_l2(2)=ry_fin(item)
            r_hexa_l2(3)=rz_fin(item)

            CALL inside_cone(r_hexa_l1,norm_hexa_l1,cone_height,&
                        cone_radius,r_hexa_l2,angle_cone,ptinsidecone)

            lupperatom(i)=ptinsidecone
            
           IF(ptinsidecone) THEN
                        atom_num=item
                        atom_type=b_fin(item)
                        atom_pos(1)=r_hexa_l2(1) 
                        atom_pos(2)=r_hexa_l2(2) 
                        atom_pos(3)=r_hexa_l2(3) 
                
             
                        dist_pt(1)=r_hexa_l2(1)-r_hexa_l1(1)
                        dist_pt(2)=r_hexa_l2(2)-r_hexa_l1(2) 
                        dist_pt(3)=r_hexa_l2(3)-r_hexa_l1(3) &
                                  -cone_height*norm_hexa_l1(3)
             
                        angle=angle_cone
                        
                        dist_plane=SQRT(DOT_PRODUCT(dist_pt,dist_pt))
                                             
        !                WRITE(6,*) 'Inside the cone', ptinsidecone,&
        !                i,atom_type,atom_num,atom_pos,angle,dist_plane
            END IF


        ELSE

               nx=norm_hexa_l1(1)
               ny=norm_hexa_l1(2)
               nz=norm_hexa_l1(3)
 
               norm_hexa_tr_l1=MATMUL(ROT_TRANSFORM(nx,ny,nz),&
                                                    norm_hexa_l1)

               
               !WRITE(6,*)nx,ny,nz,norm_hexa_tr_l1(:)
                   
               r_hexa_tr_l1=MATMUL(ROT_TRANSFORM(nx,ny,nz),& 
                                                    r_hexa_l1)

               !WRITE(6,*)r_hexa_tr_l1(:)

               r_hexa_l2(1)=rx_fin(item)
               r_hexa_l2(2)=ry_fin(item)
               r_hexa_l2(3)=rz_fin(item)

               r_hexa_tr_l2=MATMUL(ROT_TRANSFORM(nx,ny,nz),& 
                                                  r_hexa_l2)
               
                !WRITE(6,*)r_hexa_l2(:),r_hexa_tr_l2(:)
                 

             CALL inside_cone(r_hexa_tr_l1,norm_hexa_tr_l1,cone_height,&
                              cone_radius,r_hexa_tr_l2,angle_cone,&
                              ptinsidecone)

               lupperatom(i)=ptinsidecone

              IF(ptinsidecone) THEN
                         atom_num=item
                         atom_type=b_fin(item)
                         atom_tr_pos(1)=r_hexa_tr_l2(1) 
                         atom_tr_pos(2)=r_hexa_tr_l2(2) 
                         atom_tr_pos(3)=r_hexa_tr_l2(3) 
                 
              
                         atom_pos=MATMUL(ROT_TRANSFORM_INV(nx,ny,nz),&
                                                        atom_tr_pos)
                        
                         dist_pt(1)=r_hexa_tr_l2(1)-r_hexa_tr_l1(1)
                         dist_pt(2)=r_hexa_tr_l2(2)-r_hexa_tr_l1(2) 
                         dist_pt(3)=r_hexa_tr_l2(3)-r_hexa_tr_l1(3) &
                                   -cone_height*norm_hexa_tr_l1(3)
              
                         angle=angle_cone
                         
                         dist_plane=SQRT(DOT_PRODUCT(dist_pt,dist_pt))
                                              
         !                WRITE(6,*) 'Inside the cone', ptinsidecone,&
         !                i,atom_type,atom_num,atom_pos,angle,dist_plane
                 END IF
           END IF
                
         END DO

              IF(ALL([(ALL(lupperatom(i).eqv.l1),i=1,36)])) THEN
                     lnohexagon=.TRUE.
              ELSE
                     lnohexagon=.FALSE.
              END IF

        !WRITE(6,*) lnohexagon

        END SUBROUTINE find_upper_layer_atom
        
        SUBROUTINE search_atom_up_pbc(hexa_atom,i,j,r_low,norm_low,&
                    cone_h,item_up,r_up,b_up,angle_up,dist_up,latom)

        USE global_variable
        USE my_functions

        IMPLICIT NONE

        INTEGER*4,INTENT(IN)::i,j,hexa_atom
        REAL*8,INTENT(IN)::cone_h,r_low(3),norm_low(3)

        INTEGER*4,INTENT(OUT)::item_up
        REAL*8,INTENT(OUT)::r_up(3),angle_up,dist_up
        CHARACTER,INTENT(OUT)::b_up
        LOGICAL,INTENT(OUT)::latom


        INTEGER*4::i1,i2,i3,i4,i5,i6,ii,i_mod,j_mod,I_dim,J_dim,item,&
                   pos_min_dist,cnt,j1
        INTEGER*4,DIMENSION(:),ALLOCATABLE::block_radius

        REAL*8::cone_r,r_item(3),angle_item,dist_item,&
                min_dist_up

        CHARACTER::b_item

        LOGICAL::noatom

        INTEGER*4,DIMENSION(:),ALLOCATABLE::item_to_find
        REAL*8,DIMENSION(:),ALLOCATABLE::angle_to_find,dist_to_find
        REAL*8,DIMENSION(:,:),ALLOCATABLE::r_to_find
        CHARACTER,DIMENSION(:),ALLOCATABLE::b_to_find
        LOGICAL,DIMENSION(:),ALLOCATABLE::not_found
           
             I_dim=x_cell
             J_dim=y_cell

             ALLOCATE(item_to_find(5),r_to_find(5,3),b_to_find(5),&
                      angle_to_find(5),dist_to_find(5),not_found(5))
       
             j_mod= (hexa_atom-4*J_dim*(i-1)-4*(j-1))/4

             j1=j+j_mod
                
             !WRITE(6,*) hexa_atom,'Modified_j', j_mod,j,j1           
! Block 1

             ALLOCATE(block_radius(36))

             i1=4*J_dim*(i-3)+4*(j1-3)+1         
             i2=4*J_dim*(i-3)+4*j1
             i3=4*J_dim*(i-2)+4*(j1-3)+1
             i4=4*J_dim*(i-2)+4*j1
             i5=4*J_dim*(i-1)+4*(j1-3)+1
             i6=4*J_dim*(i-1)+4*j1
             
           !  WRITE(6,*) 'Block1',i,j,i1,i2,i3,i4,i5,i6
           
             cnt=0
             DO ii=i1,i2
                  cnt=cnt+1
                  block_radius(cnt)=ii
             END DO
             DO ii=i3,i4
                  cnt=cnt+1
                  block_radius(cnt)=ii
             END DO
             DO ii=i5,i6
                  cnt=cnt+1
                  block_radius(cnt)=ii
             END DO

             cone_r=0.1d0*cone_h

          100 CALL find_upper_layer_atom(block_radius,&
                 r_low,norm_low,cone_h,cone_r,item,r_item,b_item,&
                 angle_item,dist_item,noatom)
    
            IF(noatom.and.cone_r<5.0d0) THEN
!                 WRITE(6,*) 'cone_r',cone_r
                 cone_r=cone_r+0.15d0*cone_h
                 GOTO 100
            END IF
    
            IF (noatom) THEN
              not_found(1)=noatom
              !WRITE(6,*) 'Upper Hexagon',hexa_atom,'No atom found'
            ELSE
              !WRITE(6,*)'Atom_found', hexa_atom,item,r_item,b_item,&
              !          angle_item,dist_item,noatom
    
            item_to_find(1)=item
            r_to_find(1,1)=r_item(1)
            r_to_find(1,2)=r_item(2)
            r_to_find(1,3)=r_item(3)
            b_to_find(1)=b_item
            
            angle_to_find(1)=angle_item
            dist_to_find(1)=dist_item
            
            not_found(1)=noatom
           
        !   WRITE(6,*) item_to_find(1),item,r_to_find(1,:),r_item(:), &
        !              b_to_find(1),b_item,angle_to_find(1),angle_item,&
        !              dist_to_find(1),dist_item,not_found(1),noatom
   
            END IF

             DEALLOCATE(block_radius)
           
!Block 2
 
             ALLOCATE(block_radius(36))
             
             i1=4*J_dim*(i-3)+4*(j1-1)+1         
             i2=4*J_dim*(i-3)+4*(j1+2)
             i3=4*J_dim*(i-2)+4*(j1-1)+1
             i4=4*J_dim*(i-2)+4*(j1+2)
             i5=4*J_dim*(i-1)+4*(j1-1)+1
             i6=4*J_dim*(i-1)+4*(j1+2)

             !WRITE(6,*) 'Block2',i,j,i1,i2,i3,i4,i5,i6
             
             cnt=0
             DO ii=i1,i2
                  cnt=cnt+1
                  block_radius(cnt)=ii
             END DO
             DO ii=i3,i4
                  cnt=cnt+1
                  block_radius(cnt)=ii
             END DO
             DO ii=i5,i6
                  cnt=cnt+1
                  block_radius(cnt)=ii
             END DO

             cone_r=0.1d0*cone_h

          101 CALL find_upper_layer_atom(block_radius,&
                 r_low,norm_low,cone_h,cone_r,item,r_item,b_item,&
                 angle_item,dist_item,noatom)
    
           IF(noatom.and.cone_r<5.0d0) THEN
                cone_r=cone_r+0.15d0*cone_h
                GOTO 101
           END IF
    
           IF (noatom) THEN
             not_found(2)=noatom
             !WRITE(6,*) 'Upper Hexagon',hexa_atom,'No atom found'
           ELSE
            ! WRITE(6,*) 'Atom_found',hexa_atom,item,r_item,b_item,&
            !            angle_item,dist_item,noatom
    
           item_to_find(2)=item
           r_to_find(2,1)=r_item(1)
           r_to_find(2,2)=r_item(2)
           r_to_find(2,3)=r_item(3)
           b_to_find(2)=b_item
           
           angle_to_find(2)=angle_item
           dist_to_find(2)=dist_item
           
           not_found(2)=noatom
           
        !   WRITE(6,*) item_to_find(2),item,r_to_find(2,:),r_item(:), &
        !              b_to_find(2),b_item,angle_to_find(2),angle_item,&
        !              dist_to_find(2),dist_item,not_found(2),noatom
    
           END IF

            DEALLOCATE(block_radius)

!Block 3

             ALLOCATE(block_radius(36))
            
             i1=4*J_dim*(i-1)+4*(j1-3)+1         
             i2=4*J_dim*(i-1)+4*j1
             i3=4*J_dim*i+4*(j1-3)+1
             i4=4*J_dim*i+4*j1
             i5=4*J_dim*(i+1)+4*(j1-3)+1
             i6=4*J_dim*(i+1)+4*j1

             !WRITE(6,*) 'Block3',i,j,i1,i2,i3,i4,i5,i6

             cnt=0
             DO ii=i1,i2
                  cnt=cnt+1
                  block_radius(cnt)=ii
             END DO
             DO ii=i3,i4
                  cnt=cnt+1
                  block_radius(cnt)=ii
             END DO
             DO ii=i5,i6
                  cnt=cnt+1
                  block_radius(cnt)=ii
             END DO

             cone_r=0.1d0*cone_h

          102 CALL find_upper_layer_atom(block_radius,&
                 r_low,norm_low,cone_h,cone_r,item,r_item,b_item,&
                 angle_item,dist_item,noatom)
    
            IF(noatom.and.cone_r<5.0d0) THEN
                 cone_r=cone_r+0.15d0*cone_h
                 GOTO 102
            END IF
    
            IF (noatom) THEN
              not_found(3)=noatom
             ! WRITE(6,*) 'Upper Hexagon',hexa_atom,'No atom found'
            ELSE
             ! WRITE(6,*) 'Atom_found', hexa_atom,item,r_item,b_item,&
             !            angle_item,dist_item,noatom
    
            item_to_find(3)=item
            r_to_find(3,1)=r_item(1)
            r_to_find(3,2)=r_item(2)
            r_to_find(3,3)=r_item(3)
            b_to_find(3)=b_item
            
            angle_to_find(3)=angle_item
            dist_to_find(3)=dist_item
            
            not_found(3)=noatom
            
         !   WRITE(6,*) item_to_find(3),item,r_to_find(3,:),r_item(:), &
         !              b_to_find(3),b_item,angle_to_find(3),angle_item,&
         !              dist_to_find(3),dist_item,not_found(3),noatom
    
            END IF

            DEALLOCATE(block_radius)
           
!Block 4

            ALLOCATE(block_radius(36))

            i1=4*J_dim*(i-1)+4*(j1-1)+1         
            i2=4*J_dim*(i-1)+4*(j1+2)
            i3=4*J_dim*i+4*(j1-1)+1         
            i4=4*J_dim*i+4*(j1+2)
            i5=4*J_dim*(i+1)+4*(j1-1)+1
            i6=4*J_dim*(i+1)+4*(j1+2)

            !WRITE(6,*) 'Block4',i,j,i1,i2,i3,i4,i5,i6
            
            cnt=0
            DO ii=i1,i2
                 cnt=cnt+1
                 block_radius(cnt)=ii
            END DO
            DO ii=i3,i4
                 cnt=cnt+1
                 block_radius(cnt)=ii
            END DO
            DO ii=i5,i6
                 cnt=cnt+1
                 block_radius(cnt)=ii
            END DO

             cone_r=0.1d0*cone_h

          103 CALL find_upper_layer_atom(block_radius,&
                 r_low,norm_low,cone_h,cone_r,item,r_item,b_item,&
                 angle_item,dist_item,noatom)
    
           IF(noatom.and.cone_r<5.0d0) THEN
                cone_r=cone_r+0.15d0*cone_h
                GOTO 103
           END IF
    
           IF (noatom) THEN
             not_found(3)=noatom
             !WRITE(6,*) 'Upper Hexagon',hexa_atom,'No atom found'
           ELSE
             !WRITE(6,*) 'Atom_found',hexa_atom,item,r_item,b_item,&
             !            angle_item,dist_item,noatom
    
           item_to_find(4)=item
           r_to_find(4,1)=r_item(1)
           r_to_find(4,2)=r_item(2)
           r_to_find(4,3)=r_item(3)
           b_to_find(4)=b_item
           
           angle_to_find(4)=angle_item
           dist_to_find(4)=dist_item
           
           not_found(4)=noatom
           
           !WRITE(6,*) item_to_find(4),item,r_to_find(4,:),r_item(:), &
           !           b_to_find(4),b_item,angle_to_find(4),angle_item,&
           !           dist_to_find(4),dist_item,not_found(4),noatom
    
           END IF

           DEALLOCATE(block_radius)
           
!Block 5

            ALLOCATE(block_radius(36))

            i1=4*J_dim*(i-2)+4*(j1-2)+1         
            i2=4*J_dim*(i-2)+4*(j1+1)
            i3=4*J_dim*(i-1)+4*(j1-2)+1         
            i4=4*J_dim*(i-1)+4*(j1+1)
            i5=4*J_dim*i+4*(j1-2)+1
            i6=4*J_dim*i+4*(j1+1)

            !WRITE(6,*) 'Block5',i,j,i1,i2,i3,i4,i5,i6
            
            cnt=0
            DO ii=i1,i2
                 cnt=cnt+1
                 block_radius(cnt)=ii
            END DO
            DO ii=i3,i4
                 cnt=cnt+1
                 block_radius(cnt)=ii
            END DO
            DO ii=i5,i6
                 cnt=cnt+1
                 block_radius(cnt)=ii
            END DO

             cone_r=0.1d0*cone_h

          104 CALL find_upper_layer_atom(block_radius,&
                 r_low,norm_low,cone_h,cone_r,item,r_item,b_item,&
                 angle_item,dist_item,noatom)
    
            IF(noatom.and.cone_r<5.0d0) THEN
                 cone_r=cone_r+0.15d0*cone_h
                 GOTO 104
            END IF
    
            IF (noatom) THEN
              not_found(5)=noatom
             ! WRITE(6,*) 'Upper Hexagon',hexa_atom,'No atom found'
            ELSE
             ! WRITE(6,*) 'Atom_found',hexa_atom,item,r_item,b_item,&
             !            angle_item,dist_item,noatom
    
            item_to_find(5)=item
            r_to_find(5,1)=r_item(1)
            r_to_find(5,2)=r_item(2)
            r_to_find(5,3)=r_item(3)
            b_to_find(5)=b_item
            
            angle_to_find(5)=angle_item
            dist_to_find(5)=dist_item
            
            not_found(5)=noatom
            
         !   WRITE(6,*) item_to_find(5),item,r_to_find(5,:),r_item(:), &
         !              b_to_find(5),b_item,angle_to_find(5),angle_item,&
         !              dist_to_find(5),dist_item,not_found(5),noatom
    
            END IF
             min_dist_up=MINVAL(dist_to_find)
             pos_min_dist=MINLOC(dist_to_find,DIM=1)
     
            ! WRITE(6,*) min_dist_up, pos_min_dist
     
             item_up=item_to_find(pos_min_dist)
     
             r_up(1)=r_to_find(pos_min_dist,1)
             r_up(2)=r_to_find(pos_min_dist,2)
             r_up(3)=r_to_find(pos_min_dist,3)
     
             b_up=b_to_find(pos_min_dist)
             
             angle_up=angle_to_find(pos_min_dist)
             dist_up=dist_to_find(pos_min_dist)
     
             latom=not_found(pos_min_dist)
     
           !  WRITE(6,*) item_up,r_up,b_up,angle_up,dist_up,latom

        END SUBROUTINE search_atom_up_pbc
 
                !********************** 
