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
                    unit_output_1,unit_output_2,unit_output_3
                    
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
        REAL*8,DIMENSION(:,:,:),ALLOCATABLE :: hexa_norm_l1
        REAL*8,DIMENSION(:,:,:,:),ALLOCATABLE:: layer1_block_init_r,&
                                              layer1_block_fin_r, &
                                              layer2_block_init_r,&
                                              layer2_block_fin_r

        LOGICAL,DIMENSION(:,:),ALLOCATABLE:: mesh_block_l1_init,&
                                             mesh_block_l1_fin,&
                                             mesh_block_l2_init,&
                                             mesh_block_l2_fin    

        CHARACTER(len=100)::order_param_file,hexagon_data_dist_file, &
                             hexagon_data_xyz_file

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

!        REAL*8 FUNCTION average(array,tot_no)
!
!        INTEGER*4,INTENT(IN)::tot_no
!        REAL*8,INTENT(IN)::array
!        INTEGER*4::i
!        REAL*4::sum_no
!        
!        sum_no=0.0d0
!
!        DO i=1,tot_no
!                
!                sum_no=sum_no+array(i)
!
!        END DO
!
!        average=(sum_no/FLOAT(tot_no))
!
!        END FUNCTION average
!

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

        PROGRAM main
        USE global_variable
        USE my_functions        


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

!! Sort and identify the atoms in the 1st layer and commensurate 
!!the final file in accordance with the initial file

        CALL sort_index_l1

!! Find the nearest neighbors in the sorted layer 1

        CALL nearest_neighbor_l1

!! Find the hexagons in the layer 1

        CALL find_hexagon_lower_layer

! Separate out the second layers

        CALL Find_layer_l2

! Sort and identify the atoms in the 2nd layer and commensurate 
!the final file in accordance with the initial file

        CALL sort_index_l2

! Find the nearest neighbors in the sorted layer 2

        CALL nearest_neighbor_l2

! Find the order parameter

        CALL order_parameter

        DEALLOCATE(rx_init, ry_init, rz_init, rx_fin, ry_fin, rz_fin, &
                rx_sorted_l1,ry_sorted_l1, rz_sorted_l1,rx_sorted_l2, &
                ry_sorted_l2,rz_sorted_l2, atom_index_init, &
                atom_index_fin, layer1_init,layer1_fin,layer2_init,&
                layer2_fin,sorted_index_l1,sorted_index_l2,&
                nn4_sorted_l1,nn4_sorted_l2, hexa_l1,&
                cnt_atom_mb_l1_fin,cnt_atom_mb_l2_fin,&
                sum_atom_mb_l1_fin,sum_atom_mb_l2_fin,&
                layer1_block_fin,layer2_block_fin,b_fin,b_sorted_l1,&
                b_sorted_l2,layer1_block_fin_b,&
                layer2_block_fin_b,hexa_rcom_l1,hexa_rcom_l2,&
                hexa_norm_l1,layer1_block_fin_r,layer2_block_fin_r, &
                mesh_block_l1_fin,mesh_block_l2_fin) 
        
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
                   cnt5,cnt6,cnt7,cnt8,cnt9,cnt10,cnt11,l2_atom_num

        REAL*8::xmin,xmax,ymin,ymax,a0,xmin_adj,xmax_adj,ymin_adj,&
              ymax_adj,x_adjust,y_adjust,r_hexa_l2(3), &
              sum_z_l1,sum_z_l2,ave_z_l1,ave_z_l2,cone_height,&
              r_hexa_l1(3),norm_hexa_l1(3),rcom_hexa_l1(3),&
              rcom_hexa_l2(3), angle_two_lay,dist_l2,ord_pm_AA,&
              ord_pm_AA_1,ord_pm_AB,ord_pm_AB_1,ord_pm_AB_op2,&
              ord_pm_AB_inter,ord_pm_AB_no_up_atom,l1_sq,l2_sq,&
              cone_radius,temp_x,temp_y,temp_z,d_rcom
        
        LOGICAL::latoms,lsameatom,langlezero,ldisthalf,l1(6),l2(6),&
                 l3(6),l4(6),lhexagon
 
        CHARACTER::b_hexa_l1,b_hexa_l2,l2_atom_b

        INTEGER*4,DIMENSION(:,:),ALLOCATABLE::item_upper_match

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
               sum_z_l1=sum_z_l1+rz_sorted_l1(i)
         END DO
         write(6,*) sum_z_l1, ave_z_l1
         ave_z_l1=sum_z_l1/FLOAT(natoms_fin/2)
        
         sum_z_l2=0.d0
         DO i=1,tot_sort_l2
               sum_z_l2=sum_z_l2+rz_sorted_l2(i)
         END DO
         ave_z_l2=sum_z_l2/FLOAT(tot_sort_l2)

         cone_height=ABS(ave_z_l2-ave_z_l1)

         write(6,*) 'c_h',cone_height

!To find the order parameter divide the upper layer in M by M cells and
!find the atoms(atomtypes) in each cell.
! Now search the atoms on upper layer for each hexagon at the lower
! layer
         
         I_dim=x_cell
         J_dim=y_cell
         
         x_adjust=(xhi_fin-xlo_fin)/(xhi_init-xlo_init)

         y_adjust=(yhi_fin-ylo_fin)/(yhi_init-ylo_init)
         

         a0=2.504d0   ! B-N hex lattice constant

         ALLOCATE(b_upper_match(tot_hex,6),&
                 mesh_atoms(I_dim,J_dim),item_upper_match(tot_hex,6),&
                 r_upper_match(tot_hex,6,3),theta_up_down(tot_hex,6),&
                 dist_upper_match(tot_hex,6),no_atom_found(tot_hex,6),&
                 hexa_atom_match(tot_hex,6,3),hexa_op(tot_hex),&
                 hexa_rcom_l2(tot_hex,3))

        l1=(/.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE.,.TRUE./)
        l2=(/.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE./)
        l3=(/.FALSE.,.TRUE.,.FALSE.,.TRUE.,.FALSE.,.TRUE./)
        l4=(/.TRUE.,.FALSE.,.TRUE.,.FALSE.,.TRUE.,.FALSE./)
        
        OPEN(UNIT=unit_output_1,FILE=hexagon_data_dist_file)
        OPEN(UNIT=unit_output_2,FILE=order_param_file)
        OPEN(UNIT=unit_output_3,FILE=hexagon_data_xyz_file)

        cnt5=0
        cnt6=0
        cnt7=0
        cnt8=0
        cnt9=0        
        cnt10=0
        cnt11=0


        DO ii=1,tot_hex
              rcom_hexa_l1(:)=hexa_rcom_l1(ii,:) 
              !WRITE(6,*) rcom_hexa_l1
       inner: DO jj=1,6
! Lower layer hexagon points

              item=hexa_l1(ii,jj)
              r_hexa_l1(1)=rx_sorted_l1(item)
              r_hexa_l1(2)=ry_sorted_l1(item)
              r_hexa_l1(3)=rz_sorted_l1(item)
              b_hexa_l1=b_sorted_l1(item)  
       
              norm_hexa_l1(:)=hexa_norm_l1(ii,jj,:)     

             ! WRITE(6,*) ii, jj, r_hexa_l1, norm_hexa_l1, rcom_hexa_l1,&
             !            b_hexa_l1,hexa_l1(ii,jj)  
             ! WRITE(6,*) ii, jj, item
                      
              WRITE(6,*) 'Lower Hexagon', ii, jj, item, r_hexa_l1, &
                            b_hexa_l1             

! i=1 & j=1  
!              i=1;j=1
!              IF (item.le.sum_atom_mb_l1_fin(i,j))THEN
!                WRITE(6,*) ii,jj,item,i,j,sum_atom_mb_l1_fin(i,j)
!               
!               CALL search_atom_up_pbc(item,i,j,cone_height,&
!                        r_hexa_l1,norm_hexa_l1,l2_atom_num,r_hexa_l2,&
!                        l2_atom_b,angle_two_lay,dist_l2,lhexagon)
!
!                IF(lhexagon)THEN
!                WRITE(6,*) ii,jj, item,&
!                           'No atom found at the upper-layer',lhexagon
!                            no_atom_found(ii,jj)=lhexagon
!                WRITE(6,*) 'Upper Hexagon',ii,jj,item,'No atom found'
!                
!                ELSE
!
!                item_upper_match(ii,jj)=l2_atom_num
!                r_upper_match(ii,jj,1)=r_hexa_l2(1)
!                r_upper_match(ii,jj,2)=r_hexa_l2(2)
!                r_upper_match(ii,jj,3)=r_hexa_l2(3)
!                b_upper_match(ii,jj)=l2_atom_b        
!                
!                theta_up_down(ii,jj)=angle_two_lay
!                dist_upper_match(ii,jj)=dist_l2
!                
!                no_atom_found(ii,jj)=lhexagon
!                
!                WRITE(6,*) 'Upper Hexagon',ii,jj,item,&
!                item_upper_match(ii,jj),r_upper_match(ii,jj,:),&
!                b_upper_match(ii,jj),theta_up_down(ii,jj),&
!                dist_upper_match(ii,jj)
!                
!               END IF 
!              END IF
!                
!! i=1 & 2<=j<=J_dim-1
!              i=1
!              DO j=2,J_dim-1
!                IF(item.gt.sum_atom_mb_l1_fin(i,j-1) &
!                  .and.item.le.sum_atom_mb_l1_fin(i,j))THEN 
!                 WRITE(6,*)ii,jj,item,i,j,j-1,&
!                           sum_atom_mb_l1_fin(i,j-1),&
!                           sum_atom_mb_l1_fin(i,j)
!
!               CALL search_atom_up_pbc(item,i,j,cone_height,&
!                        r_hexa_l1,norm_hexa_l1,l2_atom_num,r_hexa_l2,&
!                        l2_atom_b,angle_two_lay,dist_l2,lhexagon)
!
!                IF(lhexagon)THEN
!                WRITE(6,*) ii,jj,item,&
!                           'No atom found at the upper-layer',lhexagon
!                            no_atom_found(ii,jj)=lhexagon
!                WRITE(6,*) 'Upper Hexagon',ii,jj,item,'No atom found'
!                
!                ELSE
!
!                item_upper_match(ii,jj)=l2_atom_num
!                r_upper_match(ii,jj,1)=r_hexa_l2(1)
!                r_upper_match(ii,jj,2)=r_hexa_l2(2)
!                r_upper_match(ii,jj,3)=r_hexa_l2(3)
!                b_upper_match(ii,jj)=l2_atom_b        
!                
!                theta_up_down(ii,jj)=angle_two_lay
!                dist_upper_match(ii,jj)=dist_l2
!                
!                no_atom_found(ii,jj)=lhexagon
!                
!                WRITE(6,*) 'Upper Hexagon',ii,jj,item,&
!                item_upper_match(ii,jj),r_upper_match(ii,jj,:),&
!                b_upper_match(ii,jj),theta_up_down(ii,jj),&
!                dist_upper_match(ii,jj)
!                
!               END IF 
!                 CYCLE inner
!              END IF
!             END DO
!                
!! 2<=i<=I_dim-1 & j=1
!              j=1
!              DO i=2,I_dim-1
!                IF(item.gt.sum_atom_mb_l1_fin(i-1,J_dim) &
!                  .and.item.le.sum_atom_mb_l1_fin(i,j))THEN 
!                 WRITE(6,*)ii,jj,item,i,i-1,j,&
!                           sum_atom_mb_l1_fin(i-1,J_dim),&
!                           sum_atom_mb_l1_fin(i,j)
!
!               CALL search_atom_up_pbc(item,i,j,cone_height,&
!                        r_hexa_l1,norm_hexa_l1,l2_atom_num,r_hexa_l2,&
!                        l2_atom_b,angle_two_lay,dist_l2,lhexagon)
!
!                IF(lhexagon)THEN
!                WRITE(6,*) ii,jj, item,&
!                           'No atom found at the upper-layer',lhexagon
!                            no_atom_found(ii,jj)=lhexagon
!                WRITE(6,*) 'Upper Hexagon',ii,jj,item,'No atom found'
!                
!                ELSE
!
!                item_upper_match(ii,jj)=l2_atom_num
!                r_upper_match(ii,jj,1)=r_hexa_l2(1)
!                r_upper_match(ii,jj,2)=r_hexa_l2(2)
!                r_upper_match(ii,jj,3)=r_hexa_l2(3)
!                b_upper_match(ii,jj)=l2_atom_b        
!                
!                theta_up_down(ii,jj)=angle_two_lay
!                dist_upper_match(ii,jj)=dist_l2
!                
!                no_atom_found(ii,jj)=lhexagon
!                
!                WRITE(6,*) 'Upper Hexagon',ii,jj,item,&
!                item_upper_match(ii,jj),r_upper_match(ii,jj,:),&
!                b_upper_match(ii,jj),theta_up_down(ii,jj),&
!                dist_upper_match(ii,jj)
!                
!               END IF 
!                 CYCLE inner
!              END IF
!             END DO
!                
! 2<=i<=I_dim-1 & 2<=j<=J_dim-1
              
             DO i=2,I_dim-1
              DO j=2,J_dim-1
              IF(item.gt.sum_atom_mb_l1_fin(i,j-1) &
                  .and.item.le.sum_atom_mb_l1_fin(i,j))THEN 
                  WRITE(6,*)ii,jj,item,i,j-1,j,&
                           sum_atom_mb_l1_fin(i,j-1),&
                           sum_atom_mb_l1_fin(i,j)

               CALL search_atom_up_pbc(item,i,j,cone_height,&
                        r_hexa_l1,norm_hexa_l1,l2_atom_num,r_hexa_l2,&
                        l2_atom_b,angle_two_lay,dist_l2,lhexagon)

                IF(lhexagon)THEN
                WRITE(6,*) ii,jj, item,&
                           'No atom found at the upper-layer',lhexagon
                            no_atom_found(ii,jj)=lhexagon
                WRITE(6,*) 'Upper Hexagon',ii,jj,item,'No atom found'
                
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
                
               END IF 
                 CYCLE inner
              END IF
             END DO
            END DO
                
! i=I_dim & 2<=j<=J_dim-1
!              
!              i=I_dim
!              DO j=2,J_dim-1
!                IF(item.gt.sum_atom_mb_l1_fin(i,j-1) &
!                  .and.item.le.sum_atom_mb_l1_fin(i,j))THEN 
!                 WRITE(6,*)ii,jj,item,i,j,j-1,&
!                           sum_atom_mb_l1_fin(i,j-1),&
!                           sum_atom_mb_l1_fin(i,j)
!
!                 CALL search_atom_up_pbc(item,i,j,cone_height,&
!                        r_hexa_l1,norm_hexa_l1,l2_atom_num,r_hexa_l2,&
!                        l2_atom_b,angle_two_lay,dist_l2,lhexagon)
!
!                IF(lhexagon)THEN
!                WRITE(6,*) ii,jj, item,&
!                           'No atom found at the upper-layer',lhexagon
!                            no_atom_found(ii,jj)=lhexagon
!                WRITE(6,*) 'Upper Hexagon',ii,jj,item,'No atom found'
!                
!                ELSE
!              
!                item_upper_match(ii,jj)=l2_atom_num
!                r_upper_match(ii,jj,1)=r_hexa_l2(1)
!                r_upper_match(ii,jj,2)=r_hexa_l2(2)
!                r_upper_match(ii,jj,3)=r_hexa_l2(3)
!                b_upper_match(ii,jj)=l2_atom_b        
!                
!                theta_up_down(ii,jj)=angle_two_lay
!                dist_upper_match(ii,jj)=dist_l2
!                
!                no_atom_found(ii,jj)=lhexagon
!                
!                WRITE(6,*) 'Upper Hexagon',ii,jj,item,&
!                item_upper_match(ii,jj),r_upper_match(ii,jj,:),&
!                b_upper_match(ii,jj),theta_up_down(ii,jj),&
!                dist_upper_match(ii,jj)
!                
!               END IF 
!                 CYCLE inner
!              END IF
!             END DO
!                
!! 2<=i<=I_dim-1 & j=J_dim
!              j=J_dim
!              DO i=2,I_dim-1
!                IF(item.gt.sum_atom_mb_l1_fin(i,j-1) &
!                  .and.item.le.sum_atom_mb_l1_fin(i,j))THEN 
!                 WRITE(6,*)ii,jj,item,i,j-1,j,&
!                           sum_atom_mb_l1_fin(i,j-1),&
!                           sum_atom_mb_l1_fin(i,j)
!
!                CALL search_atom_up_pbc(item,i,j,cone_height,&
!                        r_hexa_l1,norm_hexa_l1,l2_atom_num,r_hexa_l2,&
!                        l2_atom_b,angle_two_lay,dist_l2,lhexagon)
!
!                IF(lhexagon)THEN
!                WRITE(6,*) ii,jj, item,&
!                           'No atom found at the upper-layer',lhexagon
!                            no_atom_found(ii,jj)=lhexagon
!                WRITE(6,*) 'Upper Hexagon',ii,jj,item,'No atom found'
!                
!                ELSE
!              
!                item_upper_match(ii,jj)=l2_atom_num
!                r_upper_match(ii,jj,1)=r_hexa_l2(1)
!                r_upper_match(ii,jj,2)=r_hexa_l2(2)
!                r_upper_match(ii,jj,3)=r_hexa_l2(3)
!                b_upper_match(ii,jj)=l2_atom_b        
!                
!                theta_up_down(ii,jj)=angle_two_lay
!                dist_upper_match(ii,jj)=dist_l2
!                
!                no_atom_found(ii,jj)=lhexagon
!                
!                WRITE(6,*) 'Upper Hexagon',ii,jj,item,&
!                item_upper_match(ii,jj),r_upper_match(ii,jj,:),&
!                b_upper_match(ii,jj),theta_up_down(ii,jj),&
!                dist_upper_match(ii,jj)
!                
!               END IF 
!                 CYCLE inner
!              END IF
!             END DO
!                
!! i=I_dim & j=J_dim
!                
!                i=I_dim;j=J_dim 
!                IF(item.gt.sum_atom_mb_l1_fin(i,j-1) &
!                  .and.item.le.sum_atom_mb_l1_fin(i,j))THEN 
!                WRITE(6,*)ii,jj,item,i,j-1,j,&
!                           sum_atom_mb_l1_fin(i,j-1),&
!                           sum_atom_mb_l1_fin(i,j)
!                
!                CALL search_atom_up_pbc(item,i,j,cone_height,&
!                        r_hexa_l1,norm_hexa_l1,l2_atom_num,r_hexa_l2,&
!                        l2_atom_b,angle_two_lay,dist_l2,lhexagon)
!
!                IF(lhexagon)THEN
!                WRITE(6,*) ii,jj, item,&
!                           'No atom found at the upper-layer',lhexagon
!                            no_atom_found(ii,jj)=lhexagon
!                WRITE(6,*) 'Upper Hexagon',ii,jj,item,'No atom found'
!                
!                ELSE
!              
!                item_upper_match(ii,jj)=l2_atom_num
!                r_upper_match(ii,jj,1)=r_hexa_l2(1)
!                r_upper_match(ii,jj,2)=r_hexa_l2(2)
!                r_upper_match(ii,jj,3)=r_hexa_l2(3)
!                b_upper_match(ii,jj)=l2_atom_b        
!                
!                theta_up_down(ii,jj)=angle_two_lay
!                dist_upper_match(ii,jj)=dist_l2
!                
!                no_atom_found(ii,jj)=lhexagon
!                
!                WRITE(6,*) 'Upper Hexagon',ii,jj,item,&
!                item_upper_match(ii,jj),r_upper_match(ii,jj,:),&
!                b_upper_match(ii,jj),theta_up_down(ii,jj),&
!                dist_upper_match(ii,jj)
!                
!                END IF 
!              END IF
!                
       END DO inner

       hexa_rcom_l2(ii,:)=rcom(r_upper_match(ii,:,1),&
                          r_upper_match(ii,:,2),r_upper_match(ii,:,3),&
                          b_upper_match(ii,:))

       rcom_hexa_l2(:)=hexa_rcom_l2(ii,:) 

       DO i=1,3
              temp_x=rcom_hexa_l2(1)-rcom_hexa_l1(1)
              temp_y=rcom_hexa_l2(2)-rcom_hexa_l1(2)
              temp_z=rcom_hexa_l2(3)-rcom_hexa_l1(3)
              d_rcom=SQRT(temp_x**2+temp_y**2+temp_z**2)
       END DO

       WRITE(unit_output_1,*)'Hexagon#',ii 
       WRITE(unit_output_1,*)'Hexagon Info'

!       DO jj=1,6
!
!        item=hexa_l1(ii,jj)
!        b_hexa_l1=b_sorted_l1(item)
!        r_hexa_l1(1)=rx_sorted_l1(item)
!        r_hexa_l1(2)=ry_sorted_l1(item)
!        r_hexa_l1(3)=rz_sorted_l1(item)
!         
!        WRITE(unit_output_1,*) 'Lower Hexagon',ii,jj,item,r_hexa_l1,&
!                                b_hexa_l1
!
!        WRITE(unit_output_3,*) b_hexa_l1, r_hexa_l1, "#LH",ii,jj
!
!        IF(.NOT.no_atom_found(ii,jj))THEN
!
!        WRITE(unit_output_1,*) 'Upper Hexagon',ii,jj,item,&
!        item_upper_match(ii,jj),r_upper_match(ii,jj,:),&
!        b_upper_match(ii,jj),theta_up_down(ii,jj),&
!        dist_upper_match(ii,jj)
!
!        WRITE(unit_output_3,*) b_upper_match(ii,jj), &
!                              r_upper_match(ii,jj,:) , "#UH",ii,jj
!
!        IF(b_hexa_l1.eq.b_upper_match(ii,jj))THEN
!                lsameatom=.TRUE.
!                hexa_atom_match(ii,jj,1)=lsameatom
!        ELSE
!                lsameatom=.FALSE.
!                hexa_atom_match(ii,jj,1)=lsameatom
!        END IF
!
!        IF(0.d0.le.theta_up_down(ii,jj).and.&
!                   theta_up_down(ii,jj).le.10.d0) THEN
!                langlezero=.TRUE.
!                hexa_atom_match(ii,jj,2)=langlezero
!        ELSE
!                 langlezero=.FALSE.
!                 hexa_atom_match(ii,jj,2)=langlezero
!        ENDIF
!          
!        IF(0.d0.le.dist_upper_match(ii,jj).and.&
!                   dist_upper_match(ii,jj).le.0.5d0) THEN
!                ldisthalf=.TRUE.
!                hexa_atom_match(ii,jj,3)=ldisthalf
!        ELSE
!                 ldisthalf=.FALSE.
!                 hexa_atom_match(ii,jj,3)=ldisthalf
!        ENDIF
!          
!        WRITE(6,*) ii,jj,hexa_atom_match(ii,jj,1),'ATOM'
!        WRITE(6,*) ii,jj,hexa_atom_match(ii,jj,2),'ANGLE'
!        WRITE(6,*) ii,jj,hexa_atom_match(ii,jj,3),'DIST'
!
!        ELSE 
!
!       WRITE(unit_output_1,*) 'Upper Hexagon',ii,jj,item,'No atom found'
!          
!        END IF
!
!       END DO
! 
!       IF(ALL([(ALL(no_atom_found(ii,jj).eqv.l2),jj=1,6)]))THEN
!
!       WRITE(unit_output_1,*) ii,hexa_atom_match(ii,:,1),'ATOM'
!       WRITE(unit_output_1,*) ii,hexa_atom_match(ii,:,2),'ANGLE'
!       WRITE(unit_output_1,*) ii,hexa_atom_match(ii,:,3),'DIST'
!   
!       ELSE
!        
!       WRITE(unit_output_1,*) 'Atleast one of the atom at the lower',& 
!             '/','layer has no corresponding atom at the upper layer!!' 
!
!       END IF
!
!       !All hexagon has corresponding upper atom
!       IF(ALL([(ALL(no_atom_found(ii,jj).eqv.l2),jj=1,6)]))THEN
!
!          !All atom matches
!        IF(ALL([(ALL(hexa_atom_match(ii,jj,1).eqv.l1),jj=1,6)])) THEN 
!
!          !All atom at near zero degree angle
!         IF(ALL([(ALL(hexa_atom_match(ii,jj,2).eqv.l1),jj=1,6)])) THEN
!
!          !All atom at near zero plane distance
!          IF(ALL([(ALL(hexa_atom_match(ii,jj,3).eqv.l1),jj=1,6)])) THEN
!          
!                cnt5=cnt5+1
!                hexa_op(ii)='AA'
!
!          ELSE
!        
!                cnt5=cnt5+1
!                hexa_op(ii)='AA'
!
!          END IF
!
!         ELSEIF(ALL([(ALL(hexa_atom_match(ii,jj,2).eqv.l2),jj=1,6)])) &
!                                                                  THEN 
!          IF(ALL([(ALL(hexa_atom_match(ii,jj,3).eqv.l2),jj=1,6)])) THEN
!
!                cnt8=cnt8+1
!                hexa_op(ii)='AB_1'
!
!          END IF
!
!         ELSEIF(ANY([(ANY(hexa_atom_match(ii,jj,2).eqv.l2),jj=1,6)])) &
!                                                                  THEN 
!          IF(ANY([(ANY(hexa_atom_match(ii,jj,3).eqv.l2),jj=1,6)])) THEN
!
!                cnt8=cnt8+1
!                hexa_op(ii)='AB_1'
!
!          END IF
!
!         END IF
!
!       !ALL atoms are opposite
!         ELSEIF(ALL([(ALL(hexa_atom_match(ii,jj,1).eqv.l2),jj=1,6)]))THEN
!
!          !All atom at near zero degree angle
!         IF(ALL([(ALL(hexa_atom_match(ii,jj,2).eqv.l1),jj=1,6)])) THEN
!
!          !All atom at near zero distance
!          IF(ALL([(ALL(hexa_atom_match(ii,jj,3).eqv.l1),jj=1,6)])) THEN
!          
!                cnt6=cnt6+1
!                hexa_op(ii)='AA_1'
!
!          ELSE
!        
!                cnt6=cnt6+1
!                hexa_op(ii)='AA_1'
!
!          END IF
!
!         ELSEIF(ALL([(ALL(hexa_atom_match(ii,jj,2).eqv.l2),jj=1,6)])) &
!                                                                  THEN 
!
!          IF(ALL([(ALL(hexa_atom_match(ii,jj,3).eqv.l2),jj=1,6)])) THEN
!
!                cnt7=cnt7+1
!                hexa_op(ii)='AB'
!
!          END IF
!
!         ELSEIF(ANY([(ANY(hexa_atom_match(ii,jj,2).eqv.l2),jj=1,6)])) &
!                                                                  THEN
!       
!          IF(ANY([(ANY(hexa_atom_match(ii,jj,3).eqv.l2),jj=1,6)])) THEN
!
!                cnt7=cnt7+1
!                hexa_op(ii)='AB'
!
!          END IF
!
!         END IF
!
!        ELSE
!        
!        IF(ALL([(ALL(hexa_atom_match(ii,jj,1).eqv.l3),jj=1,6)]))THEN
!
!         IF(ANY([(ANY(hexa_atom_match(ii,jj,2).eqv.l2),jj=1,6)])) THEN
!          IF(ANY([(ANY(hexa_atom_match(ii,jj,3).eqv.l2),jj=1,6)])) THEN
!
!              cnt9=cnt9+1
!              hexa_op(ii)='AB_op2'
!           
!           END IF        
!          END IF
!
!        ELSEIF(ALL([(ALL(hexa_atom_match(ii,jj,1).eqv.l4),jj=1,6)]))&
!                                                                  THEN
!         IF(ANY([(ANY(hexa_atom_match(ii,jj,2).eqv.l2),jj=1,6)])) THEN
!          IF(ANY([(ANY(hexa_atom_match(ii,jj,3).eqv.l2),jj=1,6)])) THEN
!
!               cnt9=cnt9+1
!               hexa_op(ii)='AB_op2'
!
!           END IF        
!          END IF
!
!        ELSE
!               cnt10=cnt10+1
!               hexa_op(ii)='AB_inter'
!
!        END IF
!      
!       END IF          
!      
!       ELSEIF(ALL([(ALL(no_atom_found(ii,jj).eqv.l1),jj=1,6)]))THEN
!                 
!        WRITE(unit_output_1,*) "Atoms not found at upper layer for", &
!                   "every atom",'/',"in the hexagon at the lower layer"
!
!        cnt11=cnt11+1
!        hexa_op(ii)='No up atom'
!
!
!       ELSEIF(ANY([(ANY(no_atom_found(ii,jj).eqv.l1),jj=1,6)]))THEN
!             
!         WRITE(unit_output_1,*) "Atoms not found at upper layer for", &
!              "atleast one atom",'/',"in the hexagon at the lower layer"
!
!         cnt10=cnt10+1
!         hexa_op(ii)='AB_inter'
!        
!      END IF
!
        WRITE(unit_output_1,*)"Upper_hexagon_order:",'',hexa_op(ii)    
      END DO

      CLOSE(unit_output_1)

        WRITE(6,*) cnt5,cnt6,cnt7,cnt8,cnt9,cnt10,cnt11
  
 
        ord_pm_AA=cnt5*100.d0/FLOAT(tot_hex)

        ord_pm_AA_1=cnt6*100.d0/FLOAT(tot_hex)

        ord_pm_AB=cnt7*100.d0/FLOAT(tot_hex)

        ord_pm_AB_1=cnt8*100.d0/FLOAT(tot_hex)
        
        ord_pm_AB_op2=cnt9*100.d0/FLOAT(tot_hex)

        ord_pm_AB_inter=cnt10*100.d0/FLOAT(tot_hex)

        ord_pm_AB_no_up_atom=cnt11*100.d0/FLOAT(tot_hex)



        WRITE(6,*) 'The order parameter is represented in percentage'
        WRITE(6,*) '# of hexagon at 1st order(AA):', cnt5
        WRITE(6,*) 'The order parameter at 1st order(AA):',&
                                                           ord_pm_AA,'%'
        WRITE(6,*) '# of hexagon at 1st order(AA_1):', cnt6
        WRITE(6,*) 'The order parameter at 1st order(AA_1):',&
                                                         ord_pm_AA_1,'%'
        WRITE(6,*) '# of hexagon at 2nd order(AB):', cnt7
        WRITE(6,*) 'The order parameter at 2nd order(AB):',&
                                                           ord_pm_AB,'%'
        WRITE(6,*) '# of hexagon at 2nd order(AB_1):', cnt8
        WRITE(6,*) 'The order parameter at 2nd order(AB_1):',&
                                                        ord_pm_AB_1,'%'
        WRITE(6,*) '# of hexagon at for the order(AB_op2):', cnt9
        WRITE(6,*) 'The order parameter for AB_op2:',&
                                                       ord_pm_AB_op2,'%'
        WRITE(6,*) '# of hexagon at for the order(AB_inter):', cnt10
        WRITE(6,*) 'The order parameter for AB_inter:',&
                                                     ord_pm_AB_inter,'%'
        WRITE(6,*) '# of hexagon having no up atom:', cnt11
        WRITE(6,*) 'Atleast one of the atoms in lower hexagon having',&
                   'no up atom :',&
                                                ord_pm_AB_no_up_atom,'%'


        WRITE(unit_output_2,*) 'The order parameter is represented in',&
                                ' percentage'
        WRITE(unit_output_2,*) '# of hexagon at 1st order(AA):', cnt5
        WRITE(unit_output_2,*) 'The order parameter at 1st order(AA):',&
                                                           ord_pm_AA,'%'
        WRITE(unit_output_2,*) '# of hexagon at 1st order(AA_1):', cnt6
      WRITE(unit_output_2,*) 'The order parameter at 1st order(AA_1):',&
                                                         ord_pm_AA_1,'%'
        WRITE(unit_output_2,*) '# of hexagon at 2nd order(AB):', cnt7
        WRITE(unit_output_2,*) 'The order parameter at 2nd order(AB):',&
                                                           ord_pm_AB,'%'
        WRITE(unit_output_2,*) '# of hexagon at 2nd order(AB_1):', cnt8
      WRITE(unit_output_2,*) 'The order parameter at 2nd order(AB_1):',&
                                                        ord_pm_AB_1,'%'
     WRITE(unit_output_2,*) '# of hexagon with order_parameter AB_1:',&
                                                                   cnt9
        WRITE(unit_output_2,*) 'The order parameter for AB_op2:',&
                                                       ord_pm_AB_op2,'%'
     WRITE(unit_output_2,*)'# of hexagon with order_parameter AB_op2:',&
                                                                  cnt10
        WRITE(unit_output_2,*) 'The order parameter for AB_inter:',&
                                                     ord_pm_AB_inter,'%'
        WRITE(unit_output_2,*) '# of hexagon having no up atom:', cnt11
        WRITE(unit_output_2,*) 'Atleast one of the atoms in lower', &
                               'hexagon having no up atom :',&
                                                ord_pm_AB_no_up_atom,'%'
        CLOSE(unit_output_2)
        
         DEALLOCATE(b_upper_match,mesh_atoms,item_upper_match,&
                 r_upper_match,theta_up_down,hexa_atom_match,hexa_op)

        
        CALL cpu_time(finish)
    
        WRITE(6,*) '("Time taken for order_parameter = ",f6.3,&
                                    " seconds.")',finish-start
    
         END SUBROUTINE order_parameter
        
        SUBROUTINE nearest_neighbor_l2

        USE global_variable
        USE my_functions
        
        IMPLICIT NONE

        INTEGER*4::i,j,k,ii,jj,kk,atom_num(12),item1,item2,item3,cnt,&
                 I_dim,J_dim

        REAL*8:: rBN

        REAL*8 :: start, finish

        CALL cpu_time(start)

        ALLOCATE(nn4_sorted_l2(natoms_fin/2,12))
       
        rBN=rBN_l2+0.5d0

        I_dim=x_cell
        J_dim=y_cell

        nn4_sorted_l2=0
        
!For j.eq 1 and i.eq.1

         j=1
         i=1
         DO k=1,cnt_atom_mb_l2_fin(i,j)

            !WRITE(6,*) 'inner_n4','i',i,'j',j,'k',k
            CALL find_nn4_l2(i,j,k,i,i+1,j,j+1,rBN)

         END DO

!For i.eq 1 and 2<=j<=J_dim-1

         i=1
         DO j=2,J_dim-1
         DO k=1,cnt_atom_mb_l2_fin(i,j)

            !WRITE(6,*) 'inner_n4','i',i,'j',j,'k',k
            CALL find_nn4_l2(i,j,k,i,i+1,j-1,j+1,rBN)

         END DO
         END DO

!For j.eq 1 and 2<=i<=I_dim-1

         j=1
         DO i=2,I_dim-1
         DO k=1,cnt_atom_mb_l2_fin(i,j)
            
            !WRITE(6,*) 'inner_n4','i',i,'j',j,'k',k
            CALL find_nn4_l2(i,j,k,i-1,i+1,j,j+1,rBN)
                
         END DO
         END DO


! For 2 <= i <= I_dim-1 and 2 <= j <= J_dim-1

        DO i=2,I_dim-1
         DO j=2,J_dim-1
          DO k=1,cnt_atom_mb_l2_fin(i,j)

            !WRITE(6,*) 'inner_n4','i',i,'j',j,'k',k
            CALL find_nn4_l2(i,j,k,i-1,i+1,j-1,j+1,rBN)
                
           END DO
         END DO
        END DO
                                
!For i=I_dim and 2<=j<=J_dim-1

         i=I_dim
         DO j=2,J_dim-1
         DO k=1,cnt_atom_mb_l2_fin(i,j)

            !WRITE(6,*) 'inner_n4','i',i,'j',j,'k',k
            CALL find_nn4_l2(i,j,k,i-1,i,j-1,j+1,rBN)

         END DO
         END DO

!For 2<=i<=I_dim-1 and j=J_dim

         j=J_dim
         DO i=2,I_dim-1
         DO k=1,cnt_atom_mb_l2_fin(i,j)

            !WRITE(6,*) 'inner_n4','i',i,'j',j,'k',k
            CALL find_nn4_l2(i,j,k,i-1,i+1,j-1,j,rBN)
                
         END DO
         END DO


!For i=I_dim and j=J_dim
           
          i=I_dim
          j=J_dim
          DO k=1,cnt_atom_mb_l2_fin(i,j)

            !WRITE(6,*) 'inner_n4','i',i,'j',j,'k',k
            CALL find_nn4_l2(i,j,k,i-1,i,j-1,j,rBN)
                
          END DO

           DO i=1,tot_sort_l2
                                
               WRITE(6,*) i, nn4_sorted_l2(i,:)
             
           END DO 

!! Comment: If the program crashes here you must have something wrong ! with sorting the right atoms between snapshots
        
!        DO i=1,natoms_fin/2
!           item1=sorted_index_l2(i)
!           nn4_sorted_l2(i,1)=item1
!
!           cnt=0
!           DO j=1,natoms_fin/2
!              IF(i.ne.j) THEN
!                item2=sorted_index_l2(j)
!                IF(dist(rx_sorted_l2(item2),rx_sorted_l2(item1),&
!                   ry_sorted_l2(item2),ry_sorted_l2(item1),&
!                   rz_sorted_l2(item2),rz_sorted_l2(item1),&
!                   xlo_fin,xhi_fin,ylo_fin,yhi_fin,zlo_fin,zhi_fin) &
!                   .LE.rBN) THEN
!                           cnt=cnt+1
!                           nn4_sorted_l2(i,cnt+1)=item2
!       !                    WRITE(6,*) cnt,item1,item2
!                END IF
!              END IF
!           END DO
!           !WRITE(6,*) 'nn4_l2',i,nn4_sorted_l2(i,:)
!        END DO
        
        
        CALL cpu_time(finish)
    
        WRITE(6,*) '("Time taken for nearest_neighbor_l2 = ",f6.3,&
                                    " seconds.")',finish-start
    
        END SUBROUTINE nearest_neighbor_l2

        SUBROUTINE sort_index_l2

        USE global_variable
        USE my_functions

        IMPLICIT NONE

        INTEGER*4::i,j,k,cnt1,cnt2,cnt3,I_dim,J_dim,item,&
        norm_atom,ll,last,temp,cnt4, i1,i2,i3,i4
        

        REAL*8::norm,a0,rx_temp,rz_temp,ry_temp,sum_r_l2,dr

        REAL*8::xmin,xmax,ymin,ymax
        
        CHARACTER::b_temp

        LOGICAL::ldup,l_layer2_init=.False.,l_layer2_fin=.False.
        
        INTEGER count0,count1
        INTEGER cnt_rate,cnt_max,cnt_diff
        
        REAL*8 :: start, finish,wall_time
        
        CALL system_clock(COUNT_RATE=cnt_rate,COUNT_MAX=cnt_max)

        CALL system_clock(COUNT=count0)

        CALL cpu_time(start)


       ALLOCATE(layer2_block_init(x_cell,y_cell,12),&        
                layer2_block_fin(x_cell,y_cell,12),&
                cnt_atom_mb_l2_init(x_cell,y_cell),&
                cnt_atom_mb_l2_fin(x_cell,y_cell),&
                sum_atom_mb_l2_init(x_cell,y_cell),&
                sum_atom_mb_l2_fin(x_cell,y_cell),&
                layer2_block_init_r(x_cell,y_cell,12,3),&
                layer2_block_fin_r(x_cell,y_cell,12,3),&
                layer2_block_init_b(x_cell,y_cell,12),&
                layer2_block_fin_b(x_cell,y_cell,12),&
                mesh_block_l2_init(x_cell,y_cell),&
                mesh_block_l2_fin(x_cell,y_cell))
                         


! Layer2 sorted: STEP1: Identify the atoms in each block for initial and
! final frame
        
        I_dim=x_cell
 
        J_dim=y_cell

! For Static system
        a0=2.504d0   ! B-N hex lattice constant
! For MD we need to use the lattice vectors
! However, 1st thing to be sure that ordering of atom remains same
! throughout. Then define a1 and a2 locally, i.e 
! a1=(1/Sqrt(3))((x_2(4*(i-1)+2)-x_1(4*(i-1)+1)),&
! (y_2(4*(i-1)+2)-y_2(4*(i-1)+1)),(z_2(4*(i-1)+2)-z_1(4*(i-1)+1)))
! and
! a2=((x_2(4*(i-1)+5)-x_1(4*(i-1)+1)),(y_2(4*(i-1)+5)-y_1(4*(i-1)+1)),&
!    (z_2(4*(i-1)+5)-z_1(4*(i-1)+1)))
        cnt4=0
        cnt1=0

         DO i=1,I_dim
       
         DO j=1,J_dim

           IF (i==I_dim.AND.j==J_dim) THEN
             i1=4*J_dim*(i-1)+4*(j-1)+1
             i2=4*J_dim*(i-1)+4*j
             i3=4*J_dim*(i-1)+4*(j-1)+1
             i4=4*J_dim*(i-1)+4*j

             
             xmin=rx_init(layer2_init(i1))-0.5d0
             xmax=rx_init(layer2_init(i2))+0.1d0
             ymin=ry_init(layer2_init(i3))-0.5d0
             ymax=ry_init(layer2_init(i4))+0.5d0

           ELSEIF (i==I_dim.AND.j<J_dim) THEN
            
             i1=4*J_dim*(i-1)+4*(j-1)+1
             i2=4*J_dim*(i-1)+4*j
             i3=4*J_dim*(i-1)+4*(j-1)+1
             i4=4*J_dim*(i-1)+4*j+1

 
             xmin=rx_init(layer2_init(i1))-0.5d0
             xmax=rx_init(layer2_init(i2))+0.1d0
             ymin=ry_init(layer2_init(i3))-0.5d0
             ymax=ry_init(layer2_init(i4))

           ELSEIF (i<I_dim.AND.j==J_dim)THEN
             i1=4*J_dim*(i-1)+4*(j-1)+1
             i2=4*J_dim*i+4*(j-1)+1
             i3=4*J_dim*(i-1)+4*(j-1)+1
             i4=4*J_dim*(i-1)+4*j
        
             xmin=rx_init(layer2_init(i1))-0.5d0
             xmax=rx_init(layer2_init(i2))
             ymin=ry_init(layer2_init(i3))-0.5d0
             ymax=ry_init(layer2_init(i4))+0.1d0

           ELSEIF(i<I_dim.AND.j<J_dim) THEN 
            i1=4*J_dim*(i-1)+4*(j-1)+1
            i2=4*J_dim*i+4*(j-1)+1
            i3=4*J_dim*(i-1)+4*(j-1)+1
            i4=4*J_dim*(i-1)+4*j+1

            xmin=rx_init(layer2_init(i1))-0.5d0
            xmax=rx_init(layer2_init(i2))
            ymin=ry_init(layer2_init(i3))-0.5d0
            ymax=ry_init(layer2_init(i4))

           END IF

            cnt3=0
            WRITE(6,*) i,j,i1,i2,i3,i4,xmin,xmax,ymin,ymax

           DO k=1,natoms_init/2

            IF(rx_init(layer2_init(k)).gt.xmin.and.&
               rx_init(layer2_init(k)).le.xmax) THEN
            IF(ry_init(layer2_init(k)).gt.ymin.and.&
               ry_init(layer2_init(k)).le.ymax) THEN
                            
            
                item=layer2_init(k)  

            ! Removing the spurious atoms coming from MD simulations                   
               IF (i==I_dim .and. j==J_dim) THEN
                IF(item.ge.(i1+natoms_init/2).and.&
                   item.le.(i4+natoms_init/2)) THEN
               
               cnt3=cnt3+1                

               layer2_block_init(i,j,cnt3)=item
               layer2_block_init_b(i,j,cnt3)=b_init(item)
               layer2_block_init_r(i,j,cnt3,1)=rx_init(layer2_init(k))
               layer2_block_init_r(i,j,cnt3,2)=ry_init(layer2_init(k))
               layer2_block_init_r(i,j,cnt3,3)=rz_init(layer2_init(k))
               
                WRITE(6,*) i,j,k,cnt3,layer2_block_init(i,j,cnt3),&
                       layer2_block_init_b(i,j,cnt3),&
                       layer2_block_init_r(i,j,cnt3,:)

                END IF
 
              ELSEIF(i==I_dim .and. j<J_dim) THEN
                IF(item.ge.(i1+natoms_init/2).and.&
                   item.lt.(i4+natoms_init/2)) THEN

                cnt3=cnt3+1                

                layer2_block_init(i,j,cnt3)=item
                layer2_block_init_b(i,j,cnt3)=b_init(item)
                layer2_block_init_r(i,j,cnt3,1)=rx_init(layer2_init(k))
                layer2_block_init_r(i,j,cnt3,2)=ry_init(layer2_init(k))
                layer2_block_init_r(i,j,cnt3,3)=rz_init(layer2_init(k))
                
                 WRITE(6,*) i,j,k,cnt3,layer2_block_init(i,j,cnt3),&
                        layer2_block_init_b(i,j,cnt3),&
                        layer2_block_init_r(i,j,cnt3,:)

                END IF

              ELSEIF(i<I_dim .and. j==J_dim) THEN
                IF(item.ge.(i1+natoms_init/2).and.&
                   item.le.(i4+natoms_init/2)) THEN

                cnt3=cnt3+1                

                layer2_block_init(i,j,cnt3)=item
                layer2_block_init_b(i,j,cnt3)=b_init(item)
                layer2_block_init_r(i,j,cnt3,1)=rx_init(layer2_init(k))
                layer2_block_init_r(i,j,cnt3,2)=ry_init(layer2_init(k))
                layer2_block_init_r(i,j,cnt3,3)=rz_init(layer2_init(k))
                
                 WRITE(6,*) i,j,k,cnt3,layer2_block_init(i,j,cnt3),&
                        layer2_block_init_b(i,j,cnt3),&
                        layer2_block_init_r(i,j,cnt3,:)

                END IF

               ELSEIF(i<I_dim .and. j<J_dim) THEN
                IF(item.ge.(i1+natoms_init/2).and.&
                   item.lt.(i4+natoms_init/2)) THEN

                cnt3=cnt3+1                

                layer2_block_init(i,j,cnt3)=item
                layer2_block_init_b(i,j,cnt3)=b_init(item)
                layer2_block_init_r(i,j,cnt3,1)=rx_init(layer2_init(k))
                layer2_block_init_r(i,j,cnt3,2)=ry_init(layer2_init(k))
                layer2_block_init_r(i,j,cnt3,3)=rz_init(layer2_init(k))
                
                 WRITE(6,*) i,j,k,cnt3,layer2_block_init(i,j,cnt3),&
                        layer2_block_init_b(i,j,cnt3),&
                        layer2_block_init_r(i,j,cnt3,:)
             
                END IF 
               END IF
              END IF
            END IF
        END DO 

        IF (cnt3.eq.4) THEN
               l_layer2_init=.TRUE.
        ELSE 
               l_layer2_init=.FALSE.
        END IF
        mesh_block_l2_init(i,j)=l_layer2_init
        cnt_atom_mb_l2_init(i,j)=cnt3       
         
 
        cnt1=cnt1+cnt3
        
        sum_atom_mb_l2_init(i,j)=cnt1
       ! WRITE(6,*)'init', i,j,mesh_block_l2_init(i,j),&
       !               cnt_atom_mb_l2_init(i,j),&
       !            cnt1,sum_atom_mb_l2_init(i,j)
        WRITE(6,*) 'init',i,j,mesh_block_l2_init(i,j),&
                   cnt_atom_mb_l2_init(i,j),&
                   sum_atom_mb_l2_init(i,j)
       
           IF (i==I_dim.AND.j==J_dim) THEN
             i1=4*J_dim*(i-1)+4*(j-1)+1
             i2=4*J_dim*(i-1)+4*j
             i3=4*J_dim*(i-1)+4*(j-1)+1
             i4=4*J_dim*(i-1)+4*j

             
             xmin=rx_fin(layer2_fin(i1))-0.5d0
             xmax=rx_fin(layer2_fin(i2))+0.1d0
             ymin=ry_fin(layer2_fin(i3))-0.5d0
             ymax=ry_fin(layer2_fin(i4))+0.1d0

           ELSEIF (i==I_dim.AND.j<J_dim) THEN
            
             i1=4*J_dim*(i-1)+4*(j-1)+1
             i2=4*J_dim*(i-1)+4*j
             i3=4*J_dim*(i-1)+4*(j-1)+1
             i4=4*J_dim*(i-1)+4*j+1

 
             xmin=rx_fin(layer2_fin(i1))-0.5d0
             xmax=rx_fin(layer2_fin(i2))+0.1d0
             ymin=ry_fin(layer2_fin(i3))-0.5d0
             ymax=ry_fin(layer2_fin(i4))

           ELSEIF (i<I_dim.AND.j==J_dim)THEN
             i1=4*J_dim*(i-1)+4*(j-1)+1
             i2=4*J_dim*i+4*(j-1)+1
             i3=4*J_dim*(i-1)+4*(j-1)+1
             i4=4*J_dim*(i-1)+4*j
        
             xmin=rx_fin(layer2_fin(i1))-0.5d0
             xmax=rx_fin(layer2_fin(i2))
             ymin=ry_fin(layer2_fin(i3))-0.5d0
             ymax=ry_fin(layer2_fin(i4))+0.1d0

           ELSEIF(i<I_dim.AND.j<J_dim) THEN 
            i1=4*J_dim*(i-1)+4*(j-1)+1
            i2=4*J_dim*i+4*(j-1)+1
            i3=4*J_dim*(i-1)+4*(j-1)+1
            i4=4*J_dim*(i-1)+4*j+1

            xmin=rx_fin(layer2_fin(i1))-0.5d0
            xmax=rx_fin(layer2_fin(i2))
            ymin=ry_fin(layer2_fin(i3))-0.5d0
            ymax=ry_fin(layer2_fin(i4))

           END IF

            cnt3=0
            WRITE(6,*) i,j,i1,i2,i3,i4,xmin,xmax,ymin,ymax

           DO k=1,natoms_fin/2

            IF(rx_fin(layer2_fin(k)).gt.xmin.and.&
               rx_fin(layer2_fin(k)).le.xmax) THEN
            IF(ry_fin(layer2_fin(k)).gt.ymin.and.&
               ry_fin(layer2_fin(k)).le.ymax) THEN
                            
            
                item=layer2_fin(k)  

            ! Removing the spurious atoms coming from MD simulations                   
               IF (i==I_dim .and. j==J_dim) THEN
                IF(item.ge.(i1+natoms_fin/2).and.&
                   item.le.(i4+natoms_fin/2)) THEN
               
               cnt3=cnt3+1                

               layer2_block_fin(i,j,cnt3)=item
               layer2_block_fin_b(i,j,cnt3)=b_fin(item)
               layer2_block_fin_r(i,j,cnt3,1)=rx_fin(layer2_fin(k))
               layer2_block_fin_r(i,j,cnt3,2)=ry_fin(layer2_fin(k))
               layer2_block_fin_r(i,j,cnt3,3)=rz_fin(layer2_fin(k))
               
                WRITE(6,*) i,j,k,cnt3,layer2_block_fin(i,j,cnt3),&
                       layer2_block_fin_b(i,j,cnt3),&
                       layer2_block_fin_r(i,j,cnt3,:)

                END IF
 
              ELSEIF(i==I_dim .and. j<J_dim) THEN
                IF(item.ge.(i1+natoms_fin/2).and.&
                   item.lt.(i4+natoms_fin/2)) THEN

                cnt3=cnt3+1                

               layer2_block_fin(i,j,cnt3)=item
               layer2_block_fin_b(i,j,cnt3)=b_fin(item)
               layer2_block_fin_r(i,j,cnt3,1)=rx_fin(layer2_fin(k))
               layer2_block_fin_r(i,j,cnt3,2)=ry_fin(layer2_fin(k))
               layer2_block_fin_r(i,j,cnt3,3)=rz_fin(layer2_fin(k))
               
                WRITE(6,*) i,j,k,cnt3,layer2_block_fin(i,j,cnt3),&
                       layer2_block_fin_b(i,j,cnt3),&
                       layer2_block_fin_r(i,j,cnt3,:)

                END IF

              ELSEIF(i<I_dim .and. j==J_dim) THEN
                IF(item.ge.(i1+natoms_fin/2).and.&
                   item.le.(i4+natoms_fin/2)) THEN

                cnt3=cnt3+1                

               layer2_block_fin(i,j,cnt3)=item
               layer2_block_fin_b(i,j,cnt3)=b_fin(item)
               layer2_block_fin_r(i,j,cnt3,1)=rx_fin(layer2_fin(k))
               layer2_block_fin_r(i,j,cnt3,2)=ry_fin(layer2_fin(k))
               layer2_block_fin_r(i,j,cnt3,3)=rz_fin(layer2_fin(k))
               
                WRITE(6,*) i,j,k,cnt3,layer2_block_fin(i,j,cnt3),&
                       layer2_block_fin_b(i,j,cnt3),&
                       layer2_block_fin_r(i,j,cnt3,:)

                END IF

               ELSEIF(i<I_dim .and. j<J_dim) THEN
                IF(item.ge.(i1+natoms_fin/2).and.&
                   item.lt.(i4+natoms_fin/2)) THEN

                cnt3=cnt3+1                

               layer2_block_fin(i,j,cnt3)=item
               layer2_block_fin_b(i,j,cnt3)=b_fin(item)
               layer2_block_fin_r(i,j,cnt3,1)=rx_fin(layer2_fin(k))
               layer2_block_fin_r(i,j,cnt3,2)=ry_fin(layer2_fin(k))
               layer2_block_fin_r(i,j,cnt3,3)=rz_fin(layer2_fin(k))
               
                WRITE(6,*) i,j,k,cnt3,layer2_block_fin(i,j,cnt3),&
                       layer2_block_fin_b(i,j,cnt3),&
                       layer2_block_fin_r(i,j,cnt3,:)
             
             
                END IF 
               END IF
              END IF
            END IF
        END DO 

        IF (cnt3.eq.4) THEN
               l_layer2_fin=.TRUE.
        ELSE 
               l_layer2_fin=.FALSE.
        END IF
        mesh_block_l2_fin(i,j)=l_layer2_fin
        cnt_atom_mb_l2_fin(i,j)=cnt3

       cnt4=cnt4+cnt3
      
       sum_atom_mb_l2_fin(i,j)=cnt4

       ! WRITE(6,*) 'fin',i,j,mesh_block_l2_fin(i,j),&
       !                  cnt_atom_mb_l2_fin(i,j)
        WRITE(6,*) 'fin',i,j,mesh_block_l2_fin(i,j),&
                   cnt_atom_mb_l2_fin(i,j),&
                   sum_atom_mb_l2_fin(i,j)

        END DO 
       END DO

!! COMMENT: THIS PART OF THE CODE IS REDUNDENT


        ALLOCATE(b_sorted_l2(cnt4),rx_sorted_l2(cnt4),&
                 ry_sorted_l2(cnt4),rz_sorted_l2(cnt4),&
                 sorted_index_l2(cnt4))


!! Layer2 sorted: STEP2: Find the minimum of the norm between each atoms
!! in each block for consecutive 4 blocks for the two frames and record
!! the atoms at the final frame

       cnt2=0

!Case I
! i=1, j=1

        i=1
        j=1
        DO k=1,cnt_atom_mb_l2_fin(i,j)
       
           CALL sort_frame_by_frame_l2(i,j,k,i,i+1,j,j+1,norm_atom)

            cnt2=cnt2+1
       
            ll=norm_atom
            sorted_index_l2(cnt2)=ll
            rx_sorted_l2(cnt2)=rx_fin(layer2_fin(ll))
            ry_sorted_l2(cnt2)=ry_fin(layer2_fin(ll))
            rz_sorted_l2(cnt2)=rz_fin(layer2_fin(ll))
            b_sorted_l2(cnt2)=b_fin(layer2_fin(ll))

           
            WRITE(6,*) i,j,k,cnt2,norm_atom,ll,&
                     sorted_index_l2(cnt2),&
                     rx_sorted_l2(cnt2),ry_sorted_l2(cnt2),&
                     rz_sorted_l2(cnt2),b_sorted_l2(cnt2)
         END DO

!Case II
! i=1, j=J_dim

        i=1
        j=J_dim
        DO k=1,cnt_atom_mb_l2_fin(i,j)
       
           CALL sort_frame_by_frame_l2(i,j,k,i,i+1,j-1,j,norm_atom)

            cnt2=cnt2+1
       
            ll=norm_atom
            sorted_index_l2(cnt2)=ll
            rx_sorted_l2(cnt2)=rx_fin(layer2_fin(ll))
            ry_sorted_l2(cnt2)=ry_fin(layer2_fin(ll))
            rz_sorted_l2(cnt2)=rz_fin(layer2_fin(ll))
            b_sorted_l2(cnt2)=b_fin(layer2_fin(ll))

           
            WRITE(6,*) i,j,k,cnt2,norm_atom,ll,&
                     sorted_index_l2(cnt2),&
                     rx_sorted_l2(cnt2),ry_sorted_l2(cnt2),&
                     rz_sorted_l2(cnt2),b_sorted_l2(cnt2)
        END DO 

!Case III
! j=1, i=I_dim 

        i=I_dim
        j=1
        DO k=1,cnt_atom_mb_l2_fin(i,j)
       
           CALL sort_frame_by_frame_l2(i,j,k,i-1,i,j,j+1,norm_atom)

            cnt2=cnt2+1
       
            ll=norm_atom
            sorted_index_l2(cnt2)=ll
            rx_sorted_l2(cnt2)=rx_fin(layer2_fin(ll))
            ry_sorted_l2(cnt2)=ry_fin(layer2_fin(ll))
            rz_sorted_l2(cnt2)=rz_fin(layer2_fin(ll))
            b_sorted_l2(cnt2)=b_fin(layer2_fin(ll))

           
            WRITE(6,*) i,j,k,cnt2,norm_atom,ll,&
                     sorted_index_l2(cnt2),&
                     rx_sorted_l2(cnt2),ry_sorted_l2(cnt2),&
                     rz_sorted_l2(cnt2),b_sorted_l2(cnt2)
       
        END DO 

!Case IV
!i=I_dim,j=J_dim

        i=I_dim
        j=J_dim
        DO k=1,cnt_atom_mb_l2_fin(i,j)
       
           CALL sort_frame_by_frame_l2(i,j,k,i-1,i,j-1,j,norm_atom)

            cnt2=cnt2+1
       
            ll=norm_atom
            sorted_index_l2(cnt2)=ll
            rx_sorted_l2(cnt2)=rx_fin(layer2_fin(ll))
            ry_sorted_l2(cnt2)=ry_fin(layer2_fin(ll))
            rz_sorted_l2(cnt2)=rz_fin(layer2_fin(ll))
            b_sorted_l2(cnt2)=b_fin(layer2_fin(ll))

           
            WRITE(6,*) i,j,k,cnt2,norm_atom,ll,&
                     sorted_index_l2(cnt2),&
                     rx_sorted_l2(cnt2),ry_sorted_l2(cnt2),&
                     rz_sorted_l2(cnt2),b_sorted_l2(cnt2)
       
        END DO 

!Case V
!i=1,2<=j<=J_dim-1

        i=1
        DO j=2,J_dim-1
        DO k=1,cnt_atom_mb_l2_fin(i,j)
       
           CALL sort_frame_by_frame_l2(i,j,k,i,i+1,j-1,j+1,norm_atom)

            cnt2=cnt2+1
       
            ll=norm_atom
            sorted_index_l2(cnt2)=ll
            rx_sorted_l2(cnt2)=rx_fin(layer2_fin(ll))
            ry_sorted_l2(cnt2)=ry_fin(layer2_fin(ll))
            rz_sorted_l2(cnt2)=rz_fin(layer2_fin(ll))
            b_sorted_l2(cnt2)=b_fin(layer2_fin(ll))

           
            WRITE(6,*) i,j,k,cnt2,norm_atom,ll,&
                     sorted_index_l2(cnt2),&
                     rx_sorted_l2(cnt2),ry_sorted_l2(cnt2),&
                     rz_sorted_l2(cnt2),b_sorted_l2(cnt2)
       
        END DO 
        END DO

!Case VI
!i=I_dim,2<=j<=J_dim-1

        i=I_dim
        DO j=2,J_dim-1
        DO k=1,cnt_atom_mb_l2_fin(i,j)
       
           CALL sort_frame_by_frame_l2(i,j,k,i-1,i,j-1,j+1,norm_atom)

            cnt2=cnt2+1
       
            ll=norm_atom
            sorted_index_l2(cnt2)=ll
            rx_sorted_l2(cnt2)=rx_fin(layer2_fin(ll))
            ry_sorted_l2(cnt2)=ry_fin(layer2_fin(ll))
            rz_sorted_l2(cnt2)=rz_fin(layer2_fin(ll))
            b_sorted_l2(cnt2)=b_fin(layer2_fin(ll))

           
            WRITE(6,*) i,j,k,cnt2,norm_atom,ll,&
                     sorted_index_l2(cnt2),&
                     rx_sorted_l2(cnt2),ry_sorted_l2(cnt2),&
                     rz_sorted_l2(cnt2),b_sorted_l2(cnt2)
       
        END DO 
        END DO

!Case VII
!j=1,2<=i<=I_dim-1

        j=1
        DO i=2,I_dim-1
        DO k=1,cnt_atom_mb_l2_fin(i,j)
       
           CALL sort_frame_by_frame_l2(i,j,k,i-1,i+1,j,j+1,norm_atom)

            cnt2=cnt2+1
       
            ll=norm_atom
            sorted_index_l2(cnt2)=ll
            rx_sorted_l2(cnt2)=rx_fin(layer2_fin(ll))
            ry_sorted_l2(cnt2)=ry_fin(layer2_fin(ll))
            rz_sorted_l2(cnt2)=rz_fin(layer2_fin(ll))
            b_sorted_l2(cnt2)=b_fin(layer2_fin(ll))

           
            WRITE(6,*) i,j,k,cnt2,norm_atom,ll,&
                     sorted_index_l2(cnt2),&
                     rx_sorted_l2(cnt2),ry_sorted_l2(cnt2),&
                     rz_sorted_l2(cnt2),b_sorted_l2(cnt2)
       
        END DO 
        END DO

!Case VIII
!j=J_dim,2<=i<=I_dim-1

        j=J_dim
        DO i=2,I_dim-1
        DO k=1,cnt_atom_mb_l2_fin(i,j)
       
           CALL sort_frame_by_frame_l2(i,j,k,i-1,i+1,j-1,j,norm_atom)

            cnt2=cnt2+1
       
            ll=norm_atom
            sorted_index_l2(cnt2)=ll
            rx_sorted_l2(cnt2)=rx_fin(layer2_fin(ll))
            ry_sorted_l2(cnt2)=ry_fin(layer2_fin(ll))
            rz_sorted_l2(cnt2)=rz_fin(layer2_fin(ll))
            b_sorted_l2(cnt2)=b_fin(layer2_fin(ll))

           
            WRITE(6,*) i,j,k,cnt2,norm_atom,ll,&
                     sorted_index_l2(cnt2),&
                     rx_sorted_l2(cnt2),ry_sorted_l2(cnt2),&
                     rz_sorted_l2(cnt2),b_sorted_l2(cnt2)
       
        END DO 
        END DO

!Case IX
!2<=j<=J_dim,2<=i<=I_dim-1

        DO i=2,I_dim-1
        DO j=2,J_dim-1
        DO k=1,cnt_atom_mb_l2_fin(i,j)
       
           CALL sort_frame_by_frame_l2(i,j,k,i-1,i+1,j-1,j+1,norm_atom)

            cnt2=cnt2+1
       
            ll=norm_atom
            sorted_index_l2(cnt2)=ll
            rx_sorted_l2(cnt2)=rx_fin(layer2_fin(ll))
            ry_sorted_l2(cnt2)=ry_fin(layer2_fin(ll))
            rz_sorted_l2(cnt2)=rz_fin(layer2_fin(ll))
            b_sorted_l2(cnt2)=b_fin(layer2_fin(ll))

           
            WRITE(6,*) i,j,k,cnt2,norm_atom,ll,&
                     sorted_index_l2(cnt2),&
                     rx_sorted_l2(cnt2),ry_sorted_l2(cnt2),&
                     rz_sorted_l2(cnt2),b_sorted_l2(cnt2)
       
        END DO
        END DO 
        END DO

         tot_sort_l2=cnt2

!!Check Duplicates
   
!          WRITE(6,*) 'tot_atom_l2',SIZE(sorted_index_l2)

          DO i=1,SIZE(sorted_index_l2)
                WRITE(6,*) i,sorted_index_l2(i)
          END DO

         sum_r_l2=0.d0
         DO i=1,SIZE(sorted_index_l2)-1
               dr=dist(rx_sorted_l2(i+1),rx_sorted_l2(i),&
                       ry_sorted_l2(i+1),ry_sorted_l2(i),&
                       rz_sorted_l2(i+1),rz_sorted_l2(i),&
                       xlo_fin,xhi_fin,ylo_fin,yhi_fin,&
                       zlo_fin,zhi_fin)
               sum_r_l2=sum_r_l2+SQRT(dr)
         END DO
         rBN_l2=sum_r_l2/FLOAT((SIZE(sorted_index_l2)-1))
         write(6,*) sum_r_l2, rBN_l2

         
!          ldup=.FALSE.
!          DO i=1,tot_sort_l2
!           DO j=1,tot_sort_l2
!              IF(i.ne.j) THEN
!                 WRITE(6,*) i,sorted_index_l2(i),j,sorted_index_l2(j)
!                IF(sorted_index_l2(i).eq.sorted_index_l2(j)) THEN
!                        ldup=.TRUE.
!                        WRITE(6,*) 'Duplicates', i,j
!                END IF 
!             END IF
!           END DO  
!          END DO      
         
!          WRITE(6,*) 'ldup=',ldup

!!Layer2_sorted: STEP 3: Now sort the sorted_index. This is just a bubble sort (O(N^2)). The
!! time could be reduced to O(NlogN) using quicksort/heapsort. 
!
          WRITE(6,*) 'Total atoms in the 2nd Layer', cnt2, cnt4
            
          tot_sort_l2=cnt2       
          
          last=tot_sort_l2
!          CALL heap_sort_nr(last,sorted_index_l2) 

          DO i=last-1,1,-1
            DO j=1,i
                IF (sorted_index_l2(j+1).LT.sorted_index_l2(j)) THEN
                    temp=sorted_index_l2(j+1)
                    rx_temp=rx_sorted_l2(j+1)
                    ry_temp=ry_sorted_l2(j+1)
                    rz_temp=rz_sorted_l2(j+1)
                    b_temp=b_sorted_l2(j+1)
                    sorted_index_l2(j+1)=sorted_index_l2(j)
                    rx_sorted_l2(j+1)=rx_sorted_l2(j)
                    ry_sorted_l2(j+1)=ry_sorted_l2(j)
                    rz_sorted_l2(j+1)=rz_sorted_l2(j)
                    b_sorted_l2(j+1)=b_sorted_l2(j)
                    sorted_index_l2(j)=temp
                    rx_sorted_l2(j)=rx_temp
                    ry_sorted_l2(j)=ry_temp
                    rz_sorted_l2(j)=rz_temp
                    b_sorted_l2(j)=b_temp
                END IF
             END DO
           END DO
       
         OPEN(UNIT=1,FILE='layer2.dat')

          DO i=1,last
!           WRITE(6,*)'Layer2_sorted',i,sorted_index_l2(i), &
!                         rx_sorted_l2(sorted_index_l2(i)), &
!                         ry_sorted_l2(sorted_index_l2(i)), &
!                         rz_sorted_l2(sorted_index_l2(i)), &
!                         b_sorted_l2(sorted_index_l2(i))
          WRITE(6,*)'Layer2_sorted',i,sorted_index_l2(i),&
                     rx_sorted_l2(i),ry_sorted_l2(i),rz_sorted_l2(i),&
                     b_sorted_l2(i)
          WRITE(1,*) i,rx_sorted_l2(i),ry_sorted_l2(i),rz_sorted_l2(i)
          END DO

          CLOSE(1)
                 
        DEALLOCATE(layer2_block_init,cnt_atom_mb_l2_init,&
                   sum_atom_mb_l2_init,layer2_block_init_r,&
                   layer2_block_init_b, mesh_block_l2_init)

        CALL cpu_time(finish)
    
        WRITE(6,*) '("Time taken for sort_index_l2 = ",f6.3,&
                                    " seconds.")',finish-start
        
        CALL system_clock(COUNT=count1)
        cnt_diff=count1-count0
        IF(count1<count0) cnt_diff=cnt_diff+cnt_max

        wall_time=REAL(cnt_diff)/cnt_rate

        WRITE(6,*) '("Walltime = ",f6.3,"seconds")',wall_time
        
        END SUBROUTINE sort_index_l2
        
        SUBROUTINE Find_layer_l2

        USE global_variable
        USE my_functions

        IMPLICIT NONE
      
        INTEGER*4::i,j,k
        REAL*8 :: dx,dy,dz,dr 
        REAL*8 :: start, finish

        INTEGER*4::nlines,io

        INTEGER*4,DIMENSION(:),ALLOCATABLE::indices
        CHARACTER(len=9),DIMENSION(:),ALLOCATABLE::box_pos

        CALL cpu_time(start)

        nlines = 0
        OPEN (1, file = 'lay2_info_1p0_AA.xyz')
        DO
        READ(1,*,iostat=io)
        IF (io/=0) EXIT
        nlines = nlines + 1
        END DO

        WRITE(6,*) nlines

        ALLOCATE(layer2_init(natoms_init/2),layer2_fin(natoms_fin/2), &
                indices(nlines),box_pos(nlines)) 

        REWIND(1)

        DO i=1,nlines
         READ(1,*) indices(i)
         WRITE(6,*) indices(i)
        END DO  

        REWIND(1)

        DO i=1,nlines
         READ(1,110) box_pos(i)
     110 FORMAT(13X, A9)  
         WRITE(6,*) box_pos(i)
        END DO  

        CLOSE(1)


        DO i=1,natoms_init/2
              
             k=i+natoms_init/2             
  
             layer2_init(i)=atom_index_init(k)
        
                DO j=1,nlines

                  IF(layer2_init(i)==indices(j)) THEN
                     IF(box_pos(j)=='neg_box_x') THEN
                       rx_init(layer2_init(i))=rx_init(layer2_init(i))-&
                                               box_init(1)
                       write(6,*) 'i-term, -x',j 
                     ELSEIF(box_pos(j)=='pos_box_x') THEN
                       rx_init(layer2_init(i))=rx_init(layer2_init(i))+&
                                               box_init(1) 
                       write(6,*) 'i-term, +x',j 
                     ELSEIF(box_pos(j)=='neg_box_y') THEN
                       ry_init(layer2_init(i))=ry_init(layer2_init(i))-&
                                               box_init(2) 
                       write(6,*) 'i-term, -y',j 
                     ELSEIF(box_pos(j)=='pos_box_y') THEN
                       ry_init(layer2_init(i))=ry_init(layer2_init(i))+&
                                               box_init(2) 
                       write(6,*) 'i-term, +y',j 
                     END IF    
                  END IF            

                END DO

         
             IF(i.ge.1.and.i.le.x_cell*4/2) THEN
                IF(rx_init(layer2_init(i))>(xhi_init-1.d0))THEN
                   rx_init(layer2_init(i))=ABS(rx_init(layer2_init(i)) &
                                      - box_init(1))
                   WRITE(6,*)'i-term',i
                END IF
             END IF
         
             IF(i.ge.1.and.i.le.y_cell*4/2) THEN
                IF(ry_init(layer2_init(i))>(yhi_init-1.d0))THEN
                   ry_init(layer2_init(i))=ABS(ry_init(layer2_init(i)) &
                                      - box_init(2))
                   WRITE(6,*)'i-term',i
                END IF
             END IF
              
             WRITE(6,*) 'Layer2_initial',layer2_init(i),&
             rx_init(layer2_init(i)),ry_init(layer2_init(i)),&
             rz_init(layer2_init(i)),xhi_init,yhi_init,zhi_init

        END DO
           
        DO i=1,natoms_fin/2
              
             k=i+natoms_fin/2             
 
             layer2_fin(i)=atom_index_fin(k)
        
                DO j=1,nlines

                  IF(layer2_fin(i)==indices(j)) THEN
                     IF(box_pos(j)=='neg_box_x') THEN
                       rx_fin(layer2_fin(i))=rx_fin(layer2_fin(i))-&
                                               box_fin(1) 
                       write(6,*) 'i-term, -x',j 
                     ELSEIF(box_pos(j)=='pos_box_x') THEN
                       rx_fin(layer2_fin(i))=rx_fin(layer2_fin(i))+&
                                               box_fin(1) 
                       write(6,*) 'i-term, +x',j 
                     ELSEIF(box_pos(j)=='neg_box_y') THEN
                       ry_fin(layer2_fin(i))=ry_fin(layer2_fin(i))-&
                                               box_fin(2) 
                       write(6,*) 'i-term, -y',j 
                     ELSEIF(box_pos(j)=='pos_box_y') THEN
                       ry_fin(layer2_fin(i))=ry_fin(layer2_fin(i))+&
                                               box_fin(2) 
                       write(6,*) 'i-term, +y',j 
                     END IF    
                  END IF            

                END DO

             IF(i.ge.1.and.i.le.x_cell*4/2) THEN
                IF(rx_fin(layer2_fin(i))>(xhi_fin-1.d0))THEN
                   rx_fin(layer2_fin(i))=ABS(rx_fin(layer2_fin(i)) &
                                      - box_fin(1))
                   WRITE(6,*)'i-term',i
                END IF
             END IF
         
             IF(i.ge.1.and.i.le.y_cell*4/2) THEN
                IF(ry_fin(layer2_fin(i))>(yhi_fin-1.d0))THEN
                   ry_fin(layer2_fin(i))=ABS(ry_fin(layer2_fin(i)) &
                                      - box_fin(2))
                   WRITE(6,*)'i-term',i
                END IF
             END IF
                
             WRITE(6,*) 'Layer2_final',layer2_fin(i),&
             rx_fin(layer2_fin(i)),ry_fin(layer2_fin(i)),&
             rz_fin(layer2_fin(i)),xhi_fin,yhi_fin,zhi_fin
          

        END DO

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
             END IF
             
                 
         END DO 

        DO i=1,natoms_fin/2
         
             dx=rx_fin(layer2_fin(i))-rx_init(layer1_fin(i))
             dy=ry_fin(layer2_fin(i))-ry_init(layer1_fin(i))
             dz=rz_fin(layer2_fin(i))-rz_init(layer1_fin(i))
             dr=SQRT(dx*dx+dy*dy+dz*dz)

             WRITE(6,*) 'Lay1_Lay2_dist_comp',layer2_fin(i), &
                                              layer1_fin(i), &
                                              dx,dy,dz,dr
        
         END DO 

        
        CALL cpu_time(finish)
    
        WRITE(6,*) '("Time taken for find_layer_l2 = ",f6.3,&
                                    " seconds.")',finish-start
    
        END SUBROUTINE Find_layer_l2

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

        CHARACTER,DIMENSION(:),ALLOCATABLE::b_hexa
 
        REAL*8 :: start, finish

        CALL cpu_time(start)
        
        !M=99  ! number of unit cell in x direction
      
        !N=99  ! number of unit cell in y direction     

        J_dim=y_cell

        I_dim=x_cell

        J_hex=INT(FLOOR((J_dim)/3.d0))+INT(FLOOR((J_dim-1)/3.d0))
 
        J_hex_prim=INT(FLOOR(J_dim/3.d0))

        tot_hex=J_hex*(I_dim-2)+J_hex_prim
                
        
        WRITE(6,*) J_hex,J_hex_prim,tot_hex

        cnt=0

        ALLOCATE(hexa_l1(tot_hex,6),hexa_rcom_l1(tot_hex,3))

!In J direction we look for the hexagons. For I=1 and I=M-1, consider
!every isolated hexagon at lower layer after 3 unit cells each 

         DO i=1,I_dim-1

             IF(i.le.(I_dim-2))THEN

                WRITE(6,*) i                
              
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

             ELSEIF(i.eq.(I_dim-1))THEN

                WRITE(6,*) i                

                DO j=0,J_hex_prim-1

                  cnt=cnt+1
         
                  hexa_l1(cnt,1)=4*J_dim*(i-1)+3+12*j
                  hexa_l1(cnt,2)=4*J_dim*(i-1)+4+12*j
                  hexa_l1(cnt,3)=4*J_dim*i+5+12*j
                  hexa_l1(cnt,4)=4*J_dim*(i-1)+8+12*j
                  hexa_l1(cnt,5)=4*J_dim*(i-1)+7+12*j
                  hexa_l1(cnt,6)=4*J_dim*(i-1)+6+12*j

                  WRITE(6,*)cnt,i,j+1,hexa_l1(cnt,1:6)

                END DO
                
             END IF

         END DO  
      

! Find the COM for each hexagon

         DO i=1,tot_hex
           
           ALLOCATE(x_hexa(6),y_hexa(6),z_hexa(6),b_hexa(6))    

           DO j=1,6
        
                item=hexa_l1(i,j)
                
                x_hexa(j)=rx_sorted_l1(item)
                y_hexa(j)=ry_sorted_l1(item)
                z_hexa(j)=rz_sorted_l1(item)
                b_hexa(j)=b_sorted_l1(item)

               WRITE(6,*) i,j,hexa_l1(i,j),item,x_hexa(j),&
                          y_hexa(j),z_hexa(j),b_hexa(j)
                 
           END DO

           hexa_rcom_l1(i,:)=rcom(x_hexa(:),y_hexa(:),z_hexa(:),&
                                                         b_hexa(:))

           WRITE(6,*) 'hexagon_com',i, hexa_l1(i,:), hexa_rcom_l1(i,:)


           DEALLOCATE(x_hexa,y_hexa,z_hexa,b_hexa)                

         END DO
        
!Find normal for each of the points of the hexagon

          ALLOCATE(hexa_norm_l1(tot_hex,6,3))

          DO i=1,tot_hex

             DO j=1,6

                item=hexa_l1(i,j)
                item1=nn4_sorted_l1(item,1)
                item2=nn4_sorted_l1(item,2)
                item3=nn4_sorted_l1(item,3)
                item4=nn4_sorted_l1(item,4)

                v1(1)=rx_sorted_l1(item1)
                v1(2)=ry_sorted_l1(item1)
                v1(3)=rz_sorted_l1(item1)
                
                v2(1)=rx_sorted_l1(item2)
                v2(2)=ry_sorted_l1(item2)
                v2(3)=rz_sorted_l1(item2)

                v3(1)=rx_sorted_l1(item3)
                v3(2)=ry_sorted_l1(item3)
                v3(3)=rz_sorted_l1(item3)
                
                v4(1)=rx_sorted_l1(item4)
                v4(2)=ry_sorted_l1(item4)
                v4(3)=rz_sorted_l1(item4)

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


        SUBROUTINE nearest_neighbor_l1

        USE global_variable
        USE my_functions
        
        IMPLICIT NONE

        INTEGER*4::i,j,k,ii,jj,kk,atom_num(12),item1,item2,item3,cnt,&
                 I_dim,J_dim

        REAL*8:: rBN

        REAL*8 :: start, finish

        CALL cpu_time(start)

        ALLOCATE(nn4_sorted_l1(natoms_fin/2,12))
       
        rBN=rBN_l1+0.5d0

        I_dim=x_cell
        J_dim=y_cell

        nn4_sorted_l1=0

!For j.eq 1 and i.eq.1

         j=1
         i=1
         DO k=1,cnt_atom_mb_l1_fin(i,j)

            CALL find_nn4_l1(i,j,k,i,i+1,j,j+1,rBN)

         END DO

!For i.eq 1 and 2<=j<=J_dim-1

         i=1
         DO j=2,J_dim-1
         DO k=1,cnt_atom_mb_l1_fin(i,j)

            CALL find_nn4_l1(i,j,k,i,i+1,j-1,j+1,rBN)

         END DO
         END DO

!For j.eq 1 and 2<=i<=I_dim-1

         j=1
         DO i=2,I_dim-1
         DO k=1,cnt_atom_mb_l1_fin(i,j)

            CALL find_nn4_l1(i,j,k,i-1,i+1,j,j+1,rBN)
                
         END DO
         END DO


! For 2 <= i <= I_dim-1 and 2 <= j <= J_dim-1

        DO i=2,I_dim-1
         DO j=2,J_dim-1
          DO k=1,cnt_atom_mb_l1_fin(i,j)

            CALL find_nn4_l1(i,j,k,i-1,i+1,j-1,j+1,rBN)
                
           END DO
         END DO
        END DO
                                
!For i=I_dim and 2<=j<=J_dim-1

         i=I_dim
         DO j=2,J_dim-1
         DO k=1,cnt_atom_mb_l1_fin(i,j)

            CALL find_nn4_l1(i,j,k,i-1,i,j-1,j+1,rBN)

         END DO
         END DO

!For 2<=i<=I_dim-1 and j=J_dim

         j=J_dim
         DO i=2,I_dim-1
         DO k=1,cnt_atom_mb_l1_fin(i,j)

            CALL find_nn4_l1(i,j,k,i-1,i+1,j-1,j,rBN)
                
         END DO
         END DO


!For i=I_dim and j=J_dim
           
          i=I_dim
          j=J_dim
          DO k=1,cnt_atom_mb_l1_fin(i,j)

            CALL find_nn4_l1(i,j,k,i-1,i,j-1,j,rBN)
                
          END DO

           DO i=1,natoms_fin/2  

               WRITE(6,*)'nearest_neighbor', i, nn4_sorted_l1(i,:)
             
           END DO 

       ! DO i=1,natoms_fin/2
       ! WRITE(6,*) i,sorted_index_l1(i),rx_sorted_l1(i),&
       !            ry_sorted_l1(i),rz_sorted_l1(i),&
       !            b_sorted_l1(i)
       ! END DO

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


        SUBROUTINE sort_index_l1

        USE global_variable
        USE my_functions

        IMPLICIT NONE

        INTEGER*4::i,j,k,cnt1,cnt2,cnt3,I_dim,J_dim,item,&
        norm_atom,ll,last,temp,cnt4,i1,i2,i3,i4
        

        REAL*8::norm,a0,rx_temp,rz_temp,ry_temp,sum_r_l1,dr

        REAL*8::xmin,xmax,ymin,ymax
        
        CHARACTER::b_temp

        LOGICAL::ldup,l_layer1_init=.False.,l_layer1_fin=.False.

        INTEGER count0,count1

        INTEGER cnt_rate,cnt_max,cnt_diff

        REAL*8 :: start,finish,wall_time
        
        CALL system_clock(COUNT_RATE=cnt_rate,COUNT_MAX=cnt_max)

        CALL system_clock(COUNT=count0)

        CALL cpu_time(start)
        

       ALLOCATE(b_sorted_l1(natoms_fin/2),&
                rx_sorted_l1(natoms_fin/2),ry_sorted_l1(natoms_fin/2),&
                rz_sorted_l1(natoms_fin/2),&
                sorted_index_l1(natoms_fin/2),&
                layer1_block_init(x_cell,y_cell,12),&        
                layer1_block_fin(x_cell,y_cell,12),&
                cnt_atom_mb_l1_init(x_cell,y_cell),&
                cnt_atom_mb_l1_fin(x_cell,y_cell),&
                sum_atom_mb_l1_init(x_cell,y_cell),&
                sum_atom_mb_l1_fin(x_cell,y_cell),&
                layer1_block_init_r(x_cell,y_cell,12,3),&
                layer1_block_fin_r(x_cell,y_cell,12,3),&
                layer1_block_init_b(x_cell,y_cell,12),&
                layer1_block_fin_b(x_cell,y_cell,12),&
                mesh_block_l1_init(x_cell,y_cell),&
                mesh_block_l1_fin(x_cell,y_cell))
                         


! Layer1 sorted: STEP1: Identify the atoms in each block for initial and
! final frame

        I_dim=x_cell
 
        J_dim=y_cell

! For Static system
        a0=2.504d0   ! B-N hex lattice constant
! For MD we need to use the lattice vectors
! However, 1st thing to be sure that ordering of atom remains same
! throughout. Then define a1 and a2 locally, i.e 
! a1=(1/Sqrt(3))((x_2(4*(i-1)+2)-x_1(4*(i-1)+1)),&
! (y_2(4*(i-1)+2)-y_2(4*(i-1)+1)),(z_2(4*(i-1)+2)-z_1(4*(i-1)+1)))
! and
! a2=((x_2(4*(i-1)+5)-x_1(4*(i-1)+1)),(y_2(4*(i-1)+5)-y_1(4*(i-1)+1)),&
!    (z_2(4*(i-1)+5)-z_1(4*(i-1)+1)))
        cnt4=0
        cnt1=0

         DO i=1,I_dim
       
         DO j=1,J_dim

           IF (i==I_dim.AND.j==J_dim) THEN
             i1=4*J_dim*(i-1)+4*(j-1)+1
             i2=4*J_dim*(i-1)+4*j
             i3=4*J_dim*(i-1)+4*(j-1)+1
             i4=4*J_dim*(i-1)+4*j

             
             xmin=rx_init(layer1_init(i1))-0.5d0
             xmax=rx_init(layer1_init(i2))+0.1d0
             ymin=ry_init(layer1_init(i3))-0.5d0
             ymax=ry_init(layer1_init(i4))+0.5d0

           ELSEIF (i==I_dim.AND.j<J_dim) THEN
            
             i1=4*J_dim*(i-1)+4*(j-1)+1
             i2=4*J_dim*(i-1)+4*j
             i3=4*J_dim*(i-1)+4*(j-1)+1
             i4=4*J_dim*(i-1)+4*j+1

 
             xmin=rx_init(layer1_init(i1))-0.5d0
             xmax=rx_init(layer1_init(i2))+0.1d0
             ymin=ry_init(layer1_init(i3))-0.5d0
             ymax=ry_init(layer1_init(i4))

           ELSEIF (i<I_dim.AND.j==J_dim)THEN
             i1=4*J_dim*(i-1)+4*(j-1)+1
             i2=4*J_dim*i+4*(j-1)+1
             i3=4*J_dim*(i-1)+4*(j-1)+1
             i4=4*J_dim*(i-1)+4*j
        
             xmin=rx_init(layer1_init(i1))-0.5d0
             xmax=rx_init(layer1_init(i2))
             ymin=ry_init(layer1_init(i3))-0.5d0
             ymax=ry_init(layer1_init(i4))+0.1d0

           ELSEIF(i<I_dim.AND.j<J_dim) THEN 
            i1=4*J_dim*(i-1)+4*(j-1)+1
            i2=4*J_dim*i+4*(j-1)+1
            i3=4*J_dim*(i-1)+4*(j-1)+1
            i4=4*J_dim*(i-1)+4*j+1

            xmin=rx_init(layer1_init(i1))-0.5d0
            xmax=rx_init(layer1_init(i2))
            ymin=ry_init(layer1_init(i3))-0.5d0
            ymax=ry_init(layer1_init(i4))

           END IF

            cnt3=0
            WRITE(6,*) i,j,i1,i2,i3,i4,xmin,xmax,ymin,ymax

           DO k=1,natoms_init/2

            IF(rx_init(layer1_init(k)).gt.xmin.and.&
               rx_init(layer1_init(k)).le.xmax) THEN
            IF(ry_init(layer1_init(k)).gt.ymin.and.&
               ry_init(layer1_init(k)).le.ymax) THEN
                            
            
                item=layer1_init(k)  

            ! Removing the spurious atoms coming from MD simulations                   
               IF (i==I_dim .and. j==J_dim) THEN
                IF(item.ge.i1.and.item.le.i4) THEN
               
               cnt3=cnt3+1                

               layer1_block_init(i,j,cnt3)=item
               layer1_block_init_b(i,j,cnt3)=b_init(item)
               layer1_block_init_r(i,j,cnt3,1)=rx_init(layer1_init(k))
               layer1_block_init_r(i,j,cnt3,2)=ry_init(layer1_init(k))
               layer1_block_init_r(i,j,cnt3,3)=rz_init(layer1_init(k))
               
                WRITE(6,*) i,j,k,cnt3,layer1_block_init(i,j,cnt3),&
                       layer1_block_init_b(i,j,cnt3),&
                       layer1_block_init_r(i,j,cnt3,:)

                END IF
 
              ELSEIF(i==I_dim .and. j<J_dim) THEN
                IF(item.ge.i1.and.item.lt.i4) THEN

                cnt3=cnt3+1                

                layer1_block_init(i,j,cnt3)=item
                layer1_block_init_b(i,j,cnt3)=b_init(item)
                layer1_block_init_r(i,j,cnt3,1)=rx_init(layer1_init(k))
                layer1_block_init_r(i,j,cnt3,2)=ry_init(layer1_init(k))
                layer1_block_init_r(i,j,cnt3,3)=rz_init(layer1_init(k))
                
                 WRITE(6,*) i,j,k,cnt3,layer1_block_init(i,j,cnt3),&
                        layer1_block_init_b(i,j,cnt3),&
                        layer1_block_init_r(i,j,cnt3,:)

                END IF

              ELSEIF(i<I_dim .and. j==J_dim) THEN
                IF(item.ge.i1.and.item.le.i4) THEN

                cnt3=cnt3+1                

                layer1_block_init(i,j,cnt3)=item
                layer1_block_init_b(i,j,cnt3)=b_init(item)
                layer1_block_init_r(i,j,cnt3,1)=rx_init(layer1_init(k))
                layer1_block_init_r(i,j,cnt3,2)=ry_init(layer1_init(k))
                layer1_block_init_r(i,j,cnt3,3)=rz_init(layer1_init(k))
                
                 WRITE(6,*) i,j,k,cnt3,layer1_block_init(i,j,cnt3),&
                        layer1_block_init_b(i,j,cnt3),&
                        layer1_block_init_r(i,j,cnt3,:)

                END IF

               ELSEIF(i<I_dim .and. j<J_dim) THEN
                IF(item.ge.i1.and.item.lt.i4) THEN

                cnt3=cnt3+1                

                layer1_block_init(i,j,cnt3)=item
                layer1_block_init_b(i,j,cnt3)=b_init(item)
                layer1_block_init_r(i,j,cnt3,1)=rx_init(layer1_init(k))
                layer1_block_init_r(i,j,cnt3,2)=ry_init(layer1_init(k))
                layer1_block_init_r(i,j,cnt3,3)=rz_init(layer1_init(k))
                
                 WRITE(6,*) i,j,k,cnt3,layer1_block_init(i,j,cnt3),&
                        layer1_block_init_b(i,j,cnt3),&
                        layer1_block_init_r(i,j,cnt3,:)
             
                END IF 
               END IF
              END IF
            END IF
        END DO 

        IF (cnt3.eq.4) THEN
               l_layer1_init=.TRUE.
        ELSE 
               l_layer1_init=.FALSE.
        END IF
        mesh_block_l1_init(i,j)=l_layer1_init
        cnt_atom_mb_l1_init(i,j)=cnt3       
         
 
        cnt1=cnt1+cnt3
        
        sum_atom_mb_l1_init(i,j)=cnt1
       ! WRITE(6,*)'init', i,j,mesh_block_l1_init(i,j),&
       !               cnt_atom_mb_l1_init(i,j),&
       !            cnt1,sum_atom_mb_l1_init(i,j)
        WRITE(6,*) 'init',i,j,mesh_block_l1_init(i,j),&
                   cnt_atom_mb_l1_init(i,j),&
                   sum_atom_mb_l1_init(i,j)
       
           IF (i==I_dim.AND.j==J_dim) THEN
             i1=4*J_dim*(i-1)+4*(j-1)+1
             i2=4*J_dim*(i-1)+4*j
             i3=4*J_dim*(i-1)+4*(j-1)+1
             i4=4*J_dim*(i-1)+4*j

             
             xmin=rx_fin(layer1_fin(i1))-0.5d0
             xmax=rx_fin(layer1_fin(i2))+0.1d0
             ymin=ry_fin(layer1_fin(i3))-0.5d0
             ymax=ry_fin(layer1_fin(i4))+0.1d0

           ELSEIF (i==I_dim.AND.j<J_dim) THEN
            
             i1=4*J_dim*(i-1)+4*(j-1)+1
             i2=4*J_dim*(i-1)+4*j
             i3=4*J_dim*(i-1)+4*(j-1)+1
             i4=4*J_dim*(i-1)+4*j+1

 
             xmin=rx_fin(layer1_fin(i1))-0.5d0
             xmax=rx_fin(layer1_fin(i2))+0.1d0
             ymin=ry_fin(layer1_fin(i3))-0.5d0
             ymax=ry_fin(layer1_fin(i4))

           ELSEIF (i<I_dim.AND.j==J_dim)THEN
             i1=4*J_dim*(i-1)+4*(j-1)+1
             i2=4*J_dim*i+4*(j-1)+1
             i3=4*J_dim*(i-1)+4*(j-1)+1
             i4=4*J_dim*(i-1)+4*j
        
             xmin=rx_fin(layer1_fin(i1))-0.5d0
             xmax=rx_fin(layer1_fin(i2))
             ymin=ry_fin(layer1_fin(i3))-0.5d0
             ymax=ry_fin(layer1_fin(i4))+0.1d0

           ELSEIF(i<I_dim.AND.j<J_dim) THEN 
            i1=4*J_dim*(i-1)+4*(j-1)+1
            i2=4*J_dim*i+4*(j-1)+1
            i3=4*J_dim*(i-1)+4*(j-1)+1
            i4=4*J_dim*(i-1)+4*j+1

            xmin=rx_fin(layer1_fin(i1))-0.5d0
            xmax=rx_fin(layer1_fin(i2))
            ymin=ry_fin(layer1_fin(i3))-0.5d0
            ymax=ry_fin(layer1_fin(i4))

           END IF

            cnt3=0
            WRITE(6,*) i,j,i1,i2,i3,i4,xmin,xmax,ymin,ymax

           DO k=1,natoms_fin/2

            IF(rx_fin(layer1_fin(k)).gt.xmin.and.&
               rx_fin(layer1_fin(k)).le.xmax) THEN
            IF(ry_fin(layer1_fin(k)).gt.ymin.and.&
               ry_fin(layer1_fin(k)).le.ymax) THEN
                            
            
                item=layer1_fin(k)  

            ! Removing the spurious atoms coming from MD simulations                   
               IF (i==I_dim .and. j==J_dim) THEN
                IF(item.ge.i1.and.item.le.i4) THEN
               
               cnt3=cnt3+1                

               layer1_block_fin(i,j,cnt3)=item
               layer1_block_fin_b(i,j,cnt3)=b_fin(item)
               layer1_block_fin_r(i,j,cnt3,1)=rx_fin(layer1_fin(k))
               layer1_block_fin_r(i,j,cnt3,2)=ry_fin(layer1_fin(k))
               layer1_block_fin_r(i,j,cnt3,3)=rz_fin(layer1_fin(k))
               
                WRITE(6,*) i,j,k,cnt3,layer1_block_fin(i,j,cnt3),&
                       layer1_block_fin_b(i,j,cnt3),&
                       layer1_block_fin_r(i,j,cnt3,:)

                END IF
 
              ELSEIF(i==I_dim .and. j<J_dim) THEN
                IF(item.ge.i1.and.item.lt.i4) THEN

                cnt3=cnt3+1                

               layer1_block_fin(i,j,cnt3)=item
               layer1_block_fin_b(i,j,cnt3)=b_fin(item)
               layer1_block_fin_r(i,j,cnt3,1)=rx_fin(layer1_fin(k))
               layer1_block_fin_r(i,j,cnt3,2)=ry_fin(layer1_fin(k))
               layer1_block_fin_r(i,j,cnt3,3)=rz_fin(layer1_fin(k))
               
                WRITE(6,*) i,j,k,cnt3,layer1_block_fin(i,j,cnt3),&
                       layer1_block_fin_b(i,j,cnt3),&
                       layer1_block_fin_r(i,j,cnt3,:)

                END IF

              ELSEIF(i<I_dim .and. j==J_dim) THEN
                IF(item.ge.i1.and.item.le.i4) THEN

                cnt3=cnt3+1                

               layer1_block_fin(i,j,cnt3)=item
               layer1_block_fin_b(i,j,cnt3)=b_fin(item)
               layer1_block_fin_r(i,j,cnt3,1)=rx_fin(layer1_fin(k))
               layer1_block_fin_r(i,j,cnt3,2)=ry_fin(layer1_fin(k))
               layer1_block_fin_r(i,j,cnt3,3)=rz_fin(layer1_fin(k))
               
                WRITE(6,*) i,j,k,cnt3,layer1_block_fin(i,j,cnt3),&
                       layer1_block_fin_b(i,j,cnt3),&
                       layer1_block_fin_r(i,j,cnt3,:)

                END IF

               ELSEIF(i<I_dim .and. j<J_dim) THEN
                IF(item.ge.i1.and.item.lt.i4) THEN

                cnt3=cnt3+1                

               layer1_block_fin(i,j,cnt3)=item
               layer1_block_fin_b(i,j,cnt3)=b_fin(item)
               layer1_block_fin_r(i,j,cnt3,1)=rx_fin(layer1_fin(k))
               layer1_block_fin_r(i,j,cnt3,2)=ry_fin(layer1_fin(k))
               layer1_block_fin_r(i,j,cnt3,3)=rz_fin(layer1_fin(k))
               
                WRITE(6,*) i,j,k,cnt3,layer1_block_fin(i,j,cnt3),&
                       layer1_block_fin_b(i,j,cnt3),&
                       layer1_block_fin_r(i,j,cnt3,:)
             
             
                END IF 
               END IF
              END IF
            END IF
        END DO 

        IF (cnt3.eq.4) THEN
               l_layer1_fin=.TRUE.
        ELSE 
               l_layer1_fin=.FALSE.
        END IF
        mesh_block_l1_fin(i,j)=l_layer1_fin
        cnt_atom_mb_l1_fin(i,j)=cnt3

       cnt4=cnt4+cnt3
      
       sum_atom_mb_l1_fin(i,j)=cnt4

       ! WRITE(6,*) 'fin',i,j,mesh_block_l1_fin(i,j),&
       !                  cnt_atom_mb_l1_fin(i,j)
        WRITE(6,*) 'fin',i,j,mesh_block_l1_fin(i,j),&
                   cnt_atom_mb_l1_fin(i,j),&
                   sum_atom_mb_l1_fin(i,j)

        END DO 
       END DO


! Layer1 sorted: STEP2: Find the minimum of the norm between each atoms
! in each block for consecutive 4 blocks for the two frames and record
! the atoms at the final frame

       cnt2=0

!Case I
! i=1, j=1

        i=1
        j=1
        DO k=1,cnt_atom_mb_l1_fin(i,j)
       
           CALL sort_frame_by_frame_l1(i,j,k,i,i+1,j,j+1,norm_atom)

            cnt2=cnt2+1
       
            ll=norm_atom
            sorted_index_l1(cnt2)=ll
            rx_sorted_l1(cnt2)=rx_fin(layer1_fin(ll))
            ry_sorted_l1(cnt2)=ry_fin(layer1_fin(ll))
            rz_sorted_l1(cnt2)=rz_fin(layer1_fin(ll))
            b_sorted_l1(cnt2)=b_fin(layer1_fin(ll))

           
            WRITE(6,*) i,j,k,cnt2,norm_atom,ll,&
                     sorted_index_l1(cnt2),&
                     rx_sorted_l1(cnt2),ry_sorted_l1(cnt2),&
                     rz_sorted_l1(cnt2),b_sorted_l1(cnt2)
         END DO

!Case II
! i=1, j=J_dim

        i=1
        j=J_dim
        DO k=1,cnt_atom_mb_l1_fin(i,j)
       
           CALL sort_frame_by_frame_l1(i,j,k,i,i+1,j-1,j,norm_atom)

            cnt2=cnt2+1
       
            ll=norm_atom
            sorted_index_l1(cnt2)=ll
            rx_sorted_l1(cnt2)=rx_fin(layer1_fin(ll))
            ry_sorted_l1(cnt2)=ry_fin(layer1_fin(ll))
            rz_sorted_l1(cnt2)=rz_fin(layer1_fin(ll))
            b_sorted_l1(cnt2)=b_fin(layer1_fin(ll))

           
            WRITE(6,*) i,j,k,cnt2,norm_atom,ll,&
                     sorted_index_l1(cnt2),&
                     rx_sorted_l1(cnt2),ry_sorted_l1(cnt2),&
                     rz_sorted_l1(cnt2),b_sorted_l1(cnt2)
       
        END DO 

!Case III
! j=1, i=I_dim 

        i=I_dim
        j=1
        DO k=1,cnt_atom_mb_l1_fin(i,j)
       
           CALL sort_frame_by_frame_l1(i,j,k,i-1,i,j,j+1,norm_atom)

            cnt2=cnt2+1
       
            ll=norm_atom
            sorted_index_l1(cnt2)=ll
            rx_sorted_l1(cnt2)=rx_fin(layer1_fin(ll))
            ry_sorted_l1(cnt2)=ry_fin(layer1_fin(ll))
            rz_sorted_l1(cnt2)=rz_fin(layer1_fin(ll))
            b_sorted_l1(cnt2)=b_fin(layer1_fin(ll))

           
            WRITE(6,*) i,j,k,cnt2,norm_atom,ll,&
                     sorted_index_l1(cnt2),&
                     rx_sorted_l1(cnt2),ry_sorted_l1(cnt2),&
                     rz_sorted_l1(cnt2),b_sorted_l1(cnt2)
       
        END DO 

!Case IV
!i=I_dim,j=J_dim

        i=I_dim
        j=J_dim
        DO k=1,cnt_atom_mb_l1_fin(i,j)
       
           CALL sort_frame_by_frame_l1(i,j,k,i-1,i,j-1,j,norm_atom)

            cnt2=cnt2+1
       
            ll=norm_atom
            sorted_index_l1(cnt2)=ll
            rx_sorted_l1(cnt2)=rx_fin(layer1_fin(ll))
            ry_sorted_l1(cnt2)=ry_fin(layer1_fin(ll))
            rz_sorted_l1(cnt2)=rz_fin(layer1_fin(ll))
            b_sorted_l1(cnt2)=b_fin(layer1_fin(ll))

           
            WRITE(6,*) i,j,k,cnt2,norm_atom,ll,&
                     sorted_index_l1(cnt2),&
                     rx_sorted_l1(cnt2),ry_sorted_l1(cnt2),&
                     rz_sorted_l1(cnt2),b_sorted_l1(cnt2)
       
        END DO 

!Case V
!i=1,2<=j<=J_dim-1

        i=1
        DO j=2,J_dim-1
        DO k=1,cnt_atom_mb_l1_fin(i,j)
       
           CALL sort_frame_by_frame_l1(i,j,k,i,i+1,j-1,j+1,norm_atom)

            cnt2=cnt2+1
       
            ll=norm_atom
            sorted_index_l1(cnt2)=ll
            rx_sorted_l1(cnt2)=rx_fin(layer1_fin(ll))
            ry_sorted_l1(cnt2)=ry_fin(layer1_fin(ll))
            rz_sorted_l1(cnt2)=rz_fin(layer1_fin(ll))
            b_sorted_l1(cnt2)=b_fin(layer1_fin(ll))

           
            WRITE(6,*) i,j,k,cnt2,norm_atom,ll,&
                     sorted_index_l1(cnt2),&
                     rx_sorted_l1(cnt2),ry_sorted_l1(cnt2),&
                     rz_sorted_l1(cnt2),b_sorted_l1(cnt2)
       
        END DO 
        END DO

!Case VI
!i=I_dim,2<=j<=J_dim-1

        i=I_dim
        DO j=2,J_dim-1
        DO k=1,cnt_atom_mb_l1_fin(i,j)
       
           CALL sort_frame_by_frame_l1(i,j,k,i-1,i,j-1,j+1,norm_atom)

            cnt2=cnt2+1
       
            ll=norm_atom
            sorted_index_l1(cnt2)=ll
            rx_sorted_l1(cnt2)=rx_fin(layer1_fin(ll))
            ry_sorted_l1(cnt2)=ry_fin(layer1_fin(ll))
            rz_sorted_l1(cnt2)=rz_fin(layer1_fin(ll))
            b_sorted_l1(cnt2)=b_fin(layer1_fin(ll))

           
            WRITE(6,*) i,j,k,cnt2,norm_atom,ll,&
                     sorted_index_l1(cnt2),&
                     rx_sorted_l1(cnt2),ry_sorted_l1(cnt2),&
                     rz_sorted_l1(cnt2),b_sorted_l1(cnt2)
       
        END DO 
        END DO

!Case VII
!j=1,2<=i<=I_dim-1

        j=1
        DO i=2,I_dim-1
        DO k=1,cnt_atom_mb_l1_fin(i,j)
       
           CALL sort_frame_by_frame_l1(i,j,k,i-1,i+1,j,j+1,norm_atom)

            cnt2=cnt2+1
       
            ll=norm_atom
            sorted_index_l1(cnt2)=ll
            rx_sorted_l1(cnt2)=rx_fin(layer1_fin(ll))
            ry_sorted_l1(cnt2)=ry_fin(layer1_fin(ll))
            rz_sorted_l1(cnt2)=rz_fin(layer1_fin(ll))
            b_sorted_l1(cnt2)=b_fin(layer1_fin(ll))

           
            WRITE(6,*) i,j,k,cnt2,norm_atom,ll,&
                     sorted_index_l1(cnt2),&
                     rx_sorted_l1(cnt2),ry_sorted_l1(cnt2),&
                     rz_sorted_l1(cnt2),b_sorted_l1(cnt2)
       
        END DO 
        END DO

!Case VIII
!j=J_dim,2<=i<=I_dim-1

        j=J_dim
        DO i=2,I_dim-1
        DO k=1,cnt_atom_mb_l1_fin(i,j)
       
           CALL sort_frame_by_frame_l1(i,j,k,i-1,i+1,j-1,j,norm_atom)

            cnt2=cnt2+1
       
            ll=norm_atom
            sorted_index_l1(cnt2)=ll
            rx_sorted_l1(cnt2)=rx_fin(layer1_fin(ll))
            ry_sorted_l1(cnt2)=ry_fin(layer1_fin(ll))
            rz_sorted_l1(cnt2)=rz_fin(layer1_fin(ll))
            b_sorted_l1(cnt2)=b_fin(layer1_fin(ll))

           
            WRITE(6,*) i,j,k,cnt2,norm_atom,ll,&
                     sorted_index_l1(cnt2),&
                     rx_sorted_l1(cnt2),ry_sorted_l1(cnt2),&
                     rz_sorted_l1(cnt2),b_sorted_l1(cnt2)
       
        END DO 
        END DO

!Case IX
!2<=j<=J_dim,2<=i<=I_dim-1

        DO i=2,I_dim-1
        DO j=2,J_dim-1
        DO k=1,cnt_atom_mb_l1_fin(i,j)
       
           CALL sort_frame_by_frame_l1(i,j,k,i-1,i+1,j-1,j+1,norm_atom)

            cnt2=cnt2+1
       
            ll=norm_atom
            sorted_index_l1(cnt2)=ll
            rx_sorted_l1(cnt2)=rx_fin(layer1_fin(ll))
            ry_sorted_l1(cnt2)=ry_fin(layer1_fin(ll))
            rz_sorted_l1(cnt2)=rz_fin(layer1_fin(ll))
            b_sorted_l1(cnt2)=b_fin(layer1_fin(ll))

           
            WRITE(6,*) i,j,k,cnt2,norm_atom,ll,&
                     sorted_index_l1(cnt2),&
                     rx_sorted_l1(cnt2),ry_sorted_l1(cnt2),&
                     rz_sorted_l1(cnt2),b_sorted_l1(cnt2)
       
        END DO
        END DO 
        END DO

           WRITE(6,*) 'cnt2', cnt2
           DO i =1,cnt2
               WRITE(6,*) 'rx_sorted_l1', rx_sorted_l1(i), &
                          'ry_sorted_l1', ry_sorted_l1(i), &
                          'rz_sorted_l1', rz_sorted_l1(i)

           END DO

         sum_r_l1=0.d0
         DO i=1,SIZE(sorted_index_l1)-1
               dr=dist(rx_sorted_l1(i+1),rx_sorted_l1(i),&
                       ry_sorted_l1(i+1),ry_sorted_l1(i),&
                       rz_sorted_l1(i+1),rz_sorted_l1(i),&
                       xlo_fin,xhi_fin,ylo_fin,yhi_fin,&
                       zlo_fin,zhi_fin)
               sum_r_l1=sum_r_l1+SQRT(dr)
         END DO
         rBN_l1=sum_r_l1/FLOAT((SIZE(sorted_index_l1)-1))
         write(6,*) sum_r_l1, rBN_l1


!Check Duplicates
   
          ldup=.FALSE.
          DO i=1,SIZE(sorted_index_l1)
           DO j=1,SIZE(sorted_index_l1)
              IF(i.ne.j) THEN
                IF(sorted_index_l1(i).eq.sorted_index_l1(j)) THEN
                        ldup=.TRUE.
                        WRITE(6,*) 'Duplicates', i,j
                END IF 
             END IF
           END DO  
          END DO      
         
          WRITE(6,*) 'ldup=',ldup

!Layer1_sorted: STEP 3: Now sort the sorted_index. This is just a bubble sort (O(N^2)). The
! time could be reduced to O(NlogN) using quicksort/heapsort. 
   
          
          last=size(sorted_index_l1)
!          CALL heap_sort_nr(last,sorted_index_l1) 

          DO i=last-1,1,-1
            DO j=1,i
                IF (sorted_index_l1(j+1).LT.sorted_index_l1(j)) THEN
                    temp=sorted_index_l1(j+1)
                    rx_temp=rx_sorted_l1(j+1)
                    ry_temp=ry_sorted_l1(j+1)
                    rz_temp=rz_sorted_l1(j+1)
                    b_temp=b_sorted_l1(j+1)
                    sorted_index_l1(j+1)=sorted_index_l1(j)
                    rx_sorted_l1(j+1)=rx_sorted_l1(j)
                    ry_sorted_l1(j+1)=ry_sorted_l1(j)
                    rz_sorted_l1(j+1)=rz_sorted_l1(j)
                    b_sorted_l1(j+1)=b_sorted_l1(j)
                    sorted_index_l1(j)=temp
                    rx_sorted_l1(j)=rx_temp
                    ry_sorted_l1(j)=ry_temp
                    rz_sorted_l1(j)=rz_temp
                    b_sorted_l1(j)=b_temp
                END IF
             END DO
           END DO
       
           DO i=1,last
!           WRITE(6,*)'Layer1_sorted',i,sorted_index_l1(i), &
!                         rx_sorted_l1(sorted_index_l1(i)), &
!                         ry_sorted_l1(sorted_index_l1(i)), &
!                         rz_sorted_l1(sorted_index_l1(i)), &
!                         b_sorted_l1(sorted_index_l1(i))
           WRITE(6,*)'Layer1_sorted',i,sorted_index_l1(i),&
                     rx_sorted_l1(i),ry_sorted_l1(i),rz_sorted_l1(i),&
                     b_sorted_l1(i)
           END DO

                 
        DEALLOCATE(layer1_block_init,cnt_atom_mb_l1_init,&
                   sum_atom_mb_l1_init,layer1_block_init_r,&
                   layer1_block_init_b, mesh_block_l1_init)
        
        CALL cpu_time(finish)
    
        WRITE(6,*) '("Time taken for sort_index_l1 = ",f6.3,&
                                    " seconds.")',finish-start
    

        CALL system_clock(COUNT=count1)
        cnt_diff=count1-count0
        IF(count1<count0) cnt_diff=cnt_diff+cnt_max

        wall_time=REAL(cnt_diff)/cnt_rate

        WRITE(6,*) '("Walltime = ",f6.3,"seconds")',wall_time
        
        END SUBROUTINE sort_index_l1


        SUBROUTINE Find_layer_l1

        USE global_variable
        USE my_functions

      
        INTEGER*4::i,k
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
              
        END DO

! Similarity check

        DO i=1,natoms_fin/2
         
             dx=rx_fin(layer1_fin(i))-rx_init(layer1_init(i))
             dy=ry_fin(layer1_fin(i))-ry_init(layer1_init(i))
             dz=rz_fin(layer1_fin(i))-rz_init(layer1_init(i))
             dr=SQRT(dx*dx+dy*dy+dz*dz)

             WRITE(6,*) 'Lay1_dist_comp',layer1_fin(i),layer1_init(i), &
                                      dx,dy,dz,dr
        
             IF (dr > 1.0d0) THEN

               WRITE(6,*) 'The atom',i,'at the 1st layer for initial',&
                          'frame is not associated with final frame '
             END IF
             
                 
         END DO 

         WRITE(6,*) 'rz_test',rz_fin
! Find the boundaries of the lower layer
        
!        lay1_fin_inn1(1)=xlo_fin
!        lay1_fin_inn1(2)=ylo_fin
!        lay1_fin_inn1(3)=average(rz_fin,natoms_fin/2)
!        
!        lay1_fin_inn2(1)=xlo_fin
!        lay1_fin_inn2(2)=yhi_fin
!        lay1_fin_inn2(3)=average(rz_fin,natoms_fin/2)
! 
!        lay2_fin_out1(1)=xhi_fin
!        lay2_fin_out1(2)=ylo_fin
!        lay2_fin_out1(3)=average(rz_fin,natoms_fin/2)
!
!        lay2_fin_out2(1)=xhi_fin
!        lay2_fin_out2(2)=yhi_fin
!        lay2_fin_out2(3)=average(rz_fin,natoms_fin/2)
!       
!        WRITE(6,*)'Boundaries:',lay1_fin_inn1(:),lay1_fin_inn2(:),&
!                   lay1_fin_out1(:),lay1_fin_out2(:)
!         
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
        READ(5,*) hexagon_data_xyz_file
 
        WRITE(6,*)"Look for output file 1","\t", order_param_file,"\n",&
                       'output file 2',"\t",hexagon_data_dist_file
        READ(5,*)unit_input_1,unit_input_2,unit_output_1,unit_output_2,&
                 unit_output_3
        WRITE(6,*)unit_input_1,unit_input_2,unit_output_1,&
                  unit_output_2,unit_output_3
       
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
 
        cone_slope=r*r/h

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

        SUBROUTINE sort_frame_by_frame_l1(i,j,k,imin,imax,jmin,jmax,&
                                                          norm_atom)
        
        USE global_variable
        USE my_functions

        IMPLICIT NONE

        INTEGER*4,INTENT(IN)::i,j,k,imin,imax,jmin,jmax
        
        INTEGER*4,INTENT(OUT)::norm_atom

        INTEGER*4:: ii,jj,kk,cnt,cnt1,atom_num(12),min_pos,&
                  min_norm_pos(12),min_atom_num(12),norm_act
                  

        REAL*8::norm1(12),min_val,min_norm_val(12) 

        cnt=0
         DO ii=imin,imax
           DO jj=jmin,jmax
              !WRITE(6,*) i,j,k,jj,ii
              cnt1=0
            DO kk=1,cnt_atom_mb_l1_init(ii,jj)
              cnt1=cnt1+1
              norm1(kk)=dist(layer1_block_fin_r(i,j,k,1),&
                             layer1_block_init_r(ii,jj,kk,1),&
                             layer1_block_fin_r(i,j,k,2),&
                             layer1_block_init_r(ii,jj,kk,2),&
                             layer1_block_fin_r(i,j,k,3),&
                             layer1_block_init_r(ii,jj,kk,3),&
                             xlo_fin,xhi_fin,ylo_fin,yhi_fin,&
                             zlo_fin,zhi_fin)
               IF(ii.eq.1.and.jj.eq.1) THEN
                        atom_num(kk)=kk
               ELSE
                        atom_num(kk)=kk+sum_atom_mb_l1_init(ii,jj) &
                                        -cnt_atom_mb_l1_init(ii,jj)
               END IF
               
             !  WRITE(6,*) 'inner2',ii,jj,kk,sum_atom_mb_l1_init(ii,jj),&
             !              cnt_atom_mb_l1_init(ii,jj),norm1(kk),&
             !              atom_num(kk)
            END DO     
            cnt=cnt+1
        
            CALL find_min(norm1(:),cnt1,min_val,min_pos)
            min_norm_val(cnt)=min_val
            min_norm_pos(cnt)=min_pos
            min_atom_num(cnt)=atom_num(min_norm_pos(cnt))
           ! WRITE(6,*) 'inner1',ii,jj,cnt,min_norm_val(cnt),&
           !           MINLOC(norm1(:)),min_norm_pos(cnt),&
           !           min_atom_num(cnt)
           END DO
         END DO
             
         CALL find_min(min_norm_val(:),cnt,min_val,min_pos)
         norm_act=min_pos
         norm_atom=min_atom_num(norm_act)
        
        END SUBROUTINE sort_frame_by_frame_l1

        SUBROUTINE sort_frame_by_frame_l2(i,j,k,imin,imax,jmin,jmax,&
                                                          norm_atom)
        
        USE global_variable
        USE my_functions

        IMPLICIT NONE

        INTEGER*4,INTENT(IN)::i,j,k,imin,imax,jmin,jmax
        
        INTEGER*4,INTENT(OUT)::norm_atom

        INTEGER*4:: ii,jj,kk,cnt,cnt1,atom_num(12),min_pos,&
                  min_norm_pos(12),min_atom_num(12),norm_act
                  

        REAL*8::norm1(12),min_val,min_norm_val(12) 

        cnt=0
         DO ii=imin,imax
           DO jj=jmin,jmax
         !      WRITE(6,*) jj,ii
              cnt1=0
            DO kk=1,cnt_atom_mb_l2_init(ii,jj)
              cnt1=cnt1+1
              norm1(kk)=dist(layer2_block_fin_r(i,j,k,1),&
                             layer2_block_init_r(ii,jj,kk,1),&
                             layer2_block_fin_r(i,j,k,2),&
                             layer2_block_init_r(ii,jj,kk,2),&
                             layer2_block_fin_r(i,j,k,3),&
                             layer2_block_init_r(ii,jj,kk,3),&
                             xlo_fin,xhi_fin,ylo_fin,yhi_fin,&
                             zlo_fin,zhi_fin)
               IF(ii.eq.1.and.jj.eq.1) THEN
                        atom_num(kk)=kk
               ELSE
                        atom_num(kk)=kk+sum_atom_mb_l2_init(ii,jj) &
                                        -cnt_atom_mb_l2_init(ii,jj)
               END IF
                
              ! WRITE(6,*) i,j,k,ii,jj,'inner2',kk,norm1(kk),atom_num(kk)
            END DO     
            cnt=cnt+1
        
            CALL find_min(norm1(:),cnt1,min_val,min_pos)
            min_norm_val(cnt)=min_val
            min_norm_pos(cnt)=min_pos
            min_atom_num(cnt)=atom_num(min_norm_pos(cnt))
        !    WRITE(6,*) 'inner1',ii,jj,cnt,min_norm_val(cnt),&
        !              MINLOC(norm1(:)),min_norm_pos(cnt),&
        !              min_atom_num(cnt)
           END DO
         END DO
             
         CALL find_min(min_norm_val(:),cnt,min_val,min_pos)
         norm_act=min_pos
         norm_atom=min_atom_num(norm_act)
        
        END SUBROUTINE sort_frame_by_frame_l2

      SUBROUTINE find_nn4_l1(i,j,k,imin,imax,jmin,jmax,rBN)

        USE global_variable
        USE my_functions

        IMPLICIT NONE
        
        INTEGER*4,INTENT(IN)::i,j,k,imin,imax,jmin,jmax

        INTEGER*4::ii,jj,kk,item1,item2,item3,atom_num(12),cnt

        REAL*8,INTENT(IN):: rBN
        
              
        item1=layer1_block_fin(i,j,k)
        IF (i.eq.1.and.j.eq.1) THEN 
           atom_num(k)=k
           item3=atom_num(k)
        ELSE
           atom_num(k)=k+sum_atom_mb_l1_fin(i,j) &
                        -cnt_atom_mb_l1_fin(i,j)
           item3=atom_num(k)
        END IF
        nn4_sorted_l1(atom_num(k),1)=item1
!             ! WRITE(6,*) 'i',i,'j',j,'k',k,atom_num(k),item1,item3,&
!             !             nn4_sorted_l1(atom_num(k),1)
        cnt=1
        DO ii=imin,imax
          DO jj=jmin,jmax
           DO kk=1,cnt_atom_mb_l1_fin(ii,jj)
             item2=layer1_block_fin(ii,jj,kk)
             IF(item1.ne.item2) THEN
              IF(dist(layer1_block_fin_r(i,j,k,1),&
                      layer1_block_fin_r(ii,jj,kk,1),&
                      layer1_block_fin_r(i,j,k,2),&
                      layer1_block_fin_r(ii,jj,kk,2),&
                      layer1_block_fin_r(i,j,k,3),&
                      layer1_block_fin_r(ii,jj,kk,3),&
                      xlo_fin,xhi_fin,ylo_fin,yhi_fin,&
                      zlo_fin,zhi_fin).LT.rBN) THEN
                      cnt=cnt+1 
                      nn4_sorted_l1(atom_num(k),cnt)=item2
      !           WRITE(6,*)'inner_inner_n4','ii',ii,'jj',jj,'kk',&
      !                      kk,cnt,item1,item2,item3
               END IF
              ! WRITE(6,*) 'inner_n4','i',ii,'j',jj,'k',kk,&
              !                        item1,item2
             END IF
           END DO
         END DO
       END DO
 
       END SUBROUTINE find_nn4_l1

      SUBROUTINE find_nn4_l2(i,j,k,imin,imax,jmin,jmax,rBN)

        USE global_variable
        USE my_functions

        IMPLICIT NONE
        
        INTEGER*4,INTENT(IN)::i,j,k,imin,imax,jmin,jmax

        INTEGER*4::ii,jj,kk,item1,item2,item3,atom_num(12),cnt,i1,i2

        REAL*8,INTENT(IN):: rBN
       
     
        item1=layer2_block_fin(i,j,k)
        IF (i.eq.1.and.j.eq.1) THEN 
           atom_num(k)=k
           item3=atom_num(k)
        ELSE
           atom_num(k)=k+sum_atom_mb_l2_fin(i,j) &
                        -cnt_atom_mb_l2_fin(i,j)
           item3=atom_num(k)
        END IF
        !WRITE(6,*)'atom_num',i,j,k,atom_num(k),&
        !           layer2_block_fin(i,j,k)-natoms_fin/2
        i1 = (item1 - natoms_fin/2)
        nn4_sorted_l2(atom_num(k),1)=i1
!              WRITE(6,*) atom_num(k),item1,item3,&
!                         nn4_sorted_l1(atom_num(k),1)
        cnt=1
        DO ii=imin,imax
          DO jj=jmin,jmax
           DO kk=1,cnt_atom_mb_l2_fin(ii,jj)
             item2=layer2_block_fin(ii,jj,kk)
             IF(item1.ne.item2) THEN
              IF(dist(layer2_block_fin_r(i,j,k,1),&
                      layer2_block_fin_r(ii,jj,kk,1),&
                      layer2_block_fin_r(i,j,k,2),&
                      layer2_block_fin_r(ii,jj,kk,2),&
                      layer2_block_fin_r(i,j,k,3),&
                      layer2_block_fin_r(ii,jj,kk,3),&
                      xlo_fin,xhi_fin,ylo_fin,yhi_fin,&
                      zlo_fin,zhi_fin).LT.rBN) THEN
                      cnt=cnt+1 
                      i2 = (item2 - natoms_fin/2)
                      nn4_sorted_l2(atom_num(k),cnt)=i2
                ! WRITE(6,*)'inner_inner_n4','ii',ii,'jj',jj,'kk',&
                !            kk,cnt,i1,i2,item3
               END IF
             END IF
              ! WRITE(6,*) 'inner_n4','i',ii,'j',jj,'k',kk,&
              !            nn4_sorted_l2(atom_num(k),:)
           END DO
         END DO
       END DO
 
       END SUBROUTINE find_nn4_l2

      SUBROUTINE find_upper_layer_atom(imin,imax,jmin,jmax,r_hexa_l1,&
                         norm_hexa_l1,cone_height,cone_radius,atom_num,&
                         atom_pos,atom_type,angle,dist_plane,lnohexagon)

      USE global_variable
      USE my_functions

      IMPLICIT NONE
 
      INTEGER*4,INTENT(IN)::imin,imax,jmin,jmax 
      REAL*8,INTENT(IN)::r_hexa_l1(3),norm_hexa_l1(3),cone_height,&
                         cone_radius
      INTEGER*4,INTENT(OUT)::atom_num
      REAL*8,INTENT(OUT)::angle,dist_plane,atom_pos(3)
      CHARACTER,INTENT(OUT)::atom_type
      LOGICAL,INTENT(OUT)::lnohexagon
 
      INTEGER*4:: i,j,k,cnt,cnt1,diff,imodmin,imodmax,jmodmin,jmodmax
      REAL*8::dist_pt(3),r_hexa_l2(3),norm_hexa_tr_l1(3),angle_cone,&
              r_hexa_tr_l1(3),r_hexa_tr_l2(3),atom_tr_pos(3),nx,ny,nz
      LOGICAL::ptinsidecone
 
      LOGICAL,DIMENSION(:),ALLOCATABLE::lupperatom(:),l1(:),l2(:),l3(:)

 
      !WRITE(6,*) imin,imax,jmin,jmax,r_hexa_l1,norm_hexa_l1,cone_height
        
        cnt1=0      
       
        WRITE(6,*) 'tg2',imin,imax,jmin,jmax
        IF (imax<imin.and.jmax<jmin) THEN
                imodmin=imax-1
                jmodmin=jmax-1
                imodmax=imax
                jmodmax=jmax
        ELSE IF(jmax<jmin.and.imax>imin) THEN
                imodmax=imax
                imodmin=imin
                jmodmax=jmax
                jmodmin=jmax-1
        ELSE IF(imax<imin.and.jmax>jmin) THEN
                imodmax=imax
                imodmin=imax-1
                jmodmax=jmax
                jmodmin=jmin
        ELSE
                imodmax=imax
                imodmin=imin
                jmodmax=jmax
                jmodmin=jmin
        END IF
                
 
        diff=(imodmax-imodmin+1)*(jmodmax-jmodmin+1)
        
       ! WRITE(6,*) diff

        ALLOCATE(l2(diff),l3(diff))

        WRITE(6,*) imodmin,imodmax,jmodmin,jmodmax

        DO i=imodmin,imodmax
        inner: DO j=jmodmin,jmodmax
           cnt=0
           write(6,*) imodmin,imodmax,jmodmin,jmodmax
           IF(cnt_atom_mb_l2_fin(i,j).eq.0) THEN
                cnt1=cnt1+1
                l3(cnt1)=.FALSE.
                l2(cnt1)=.FALSE.
                !WRITE(6,*) 'No atom at the upper layer', i,j
                CYCLE inner
           ELSE
           cnt1=cnt1+1
           l3(cnt1)=.FALSE.
           ALLOCATE(lupperatom(cnt_atom_mb_l2_fin(i,j)),&
                    l1(cnt_atom_mb_l2_fin(i,j)))
           DO k=1,cnt_atom_mb_l2_fin(i,j)
                
                 l1(k)=.FALSE.
                 
                 cnt=cnt+1

                IF((norm_hexa_l1(1).eq.0.d0).AND.&
                   (norm_hexa_l1(2).eq.0.d0).AND.&
                   (norm_hexa_l1(3).eq.1.d0))THEN
                   

                !PBC kicks in
                IF(imodmin.ne.imax.or.jmodmin.ne.jmax) THEN
                        r_hexa_l2(1)=layer2_block_fin_r(i,j,cnt,1) &
                                                               -xhi_fin 
                        r_hexa_l2(2)=layer2_block_fin_r(i,j,cnt,2)&
                                                               -xlo_fin
                        r_hexa_l2(3)=layer2_block_fin_r(i,j,cnt,3)
                ELSE
                        r_hexa_l2(1)=layer2_block_fin_r(i,j,cnt,1)
                        r_hexa_l2(2)=layer2_block_fin_r(i,j,cnt,2)
                        r_hexa_l2(3)=layer2_block_fin_r(i,j,cnt,3)
                END IF
                

                CALL inside_cone(r_hexa_l1,norm_hexa_l1,cone_height,&
                                cone_radius,r_hexa_l2,angle_cone,&
                                ptinsidecone)

                lupperatom(k)=ptinsidecone
 
                IF(ptinsidecone) THEN
                        atom_num=layer2_block_fin(i,j,cnt)-natoms_fin/2
                        atom_type=layer2_block_fin_b(i,j,cnt)
                        atom_pos(1)=r_hexa_l2(1) 
                        atom_pos(2)=r_hexa_l2(2) 
                        atom_pos(3)=r_hexa_l2(3) 
                
             
                        dist_pt(1)=r_hexa_l2(1)-r_hexa_l1(1)
                        dist_pt(2)=r_hexa_l2(2)-r_hexa_l1(2) 
                        dist_pt(3)=r_hexa_l2(3)-r_hexa_l1(3) &
                                  -cone_height*norm_hexa_l1(3)
             
                        angle=angle_cone
                        
                        dist_plane=SQRT(DOT_PRODUCT(dist_pt,dist_pt))
                                             
             !           WRITE(6,*) 'Inside the cone', ptinsidecone,&
             !           i,j,k,cnt,atom_type,atom_num,atom_pos,angle,&
             !           dist_plane
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

               !PBC kicks in
                IF(imodmin.ne.imax.or.jmodmin.ne.jmax) THEN
                        r_hexa_l2(1)=layer2_block_fin_r(i,j,cnt,1) &
                                                               -xhi_fin 
                        r_hexa_l2(2)=layer2_block_fin_r(i,j,cnt,2)&
                                                               -xlo_fin
                        r_hexa_l2(3)=layer2_block_fin_r(i,j,cnt,3)
                ELSE
                        r_hexa_l2(1)=layer2_block_fin_r(i,j,cnt,1)
                        r_hexa_l2(2)=layer2_block_fin_r(i,j,cnt,2)
                        r_hexa_l2(3)=layer2_block_fin_r(i,j,cnt,3)
                END IF

               r_hexa_tr_l2=MATMUL(ROT_TRANSFORM(nx,ny,nz),& 
                                                  r_hexa_l2)
               
                !WRITE(6,*)r_hexa_l2(:),r_hexa_tr_l2(:)
                 

             CALL inside_cone(r_hexa_tr_l1,norm_hexa_tr_l1,cone_height,&
                              cone_radius,r_hexa_tr_l2,angle_cone,&
                              ptinsidecone)

                 lupperatom(k)=ptinsidecone
 
                 IF(ptinsidecone) THEN
                         atom_num=layer2_block_fin(i,j,cnt)-natoms_fin/2
                         atom_type=layer2_block_fin_b(i,j,cnt)
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
                                              
              !           WRITE(6,*) 'Inside the cone', ptinsidecone,&
              !           i,j,k,cnt,atom_type,atom_num,atom_pos,angle,&
              !           dist_plane
                 END IF
                END IF
                
              END DO

             !WRITE(6,*) 'lupperatom',lupperatom,'l1',l1

            ! WRITE(6,*) 'cnt',cnt,'cnt1',cnt1  
             IF(ALL([(ALL(lupperatom(i).eqv.l1),i=1,cnt)])) THEN
                       l2(cnt1)=.FALSE.
            ! WRITE(6,*) l2(cnt1)  
             ELSE
                       l2(cnt1)=.TRUE.
            ! WRITE(6,*) l2(cnt1)  
             END IF

             END IF             
             
             DEALLOCATE(lupperatom,l1) 

            END DO inner
           END DO
           !WRITE(6,*) 'l2',l2,'l3',l3
           IF(ALL([(ALL(l3(i).eqv.l2),i=1,cnt1)])) THEN
                  lnohexagon=.TRUE.
            ! WRITE(6,*) lnohexagon 
           ELSE 
                  lnohexagon=.FALSE.
            ! WRITE(6,*) lnohexagon 
           END IF
                  
           !WRITE(6,*) ' I am at the end of the subroutine',&
           !           'find_upper_layer_atom'
           DEALLOCATE(l2,l3) 
        END SUBROUTINE find_upper_layer_atom
        
        SUBROUTINE search_atom_up_pbc(hexa_atom,ibot,jbot,cone_h,&
                    r_low,norm_low,item_up,r_up,b_up,angle_up,&
                    dist_up,latom)

        USE global_variable
        USE my_functions

        IMPLICIT NONE

        INTEGER*4,INTENT(IN)::hexa_atom,ibot,jbot
        REAL*8,INTENT(IN)::cone_h,r_low(3),norm_low(3)

        INTEGER*4,INTENT(OUT)::item_up
        REAL*8,INTENT(OUT)::r_up(3),angle_up,dist_up
        CHARACTER,INTENT(OUT)::b_up
        LOGICAL,INTENT(OUT)::latom

        INTEGER*4::item,imod,imodpl,imodmin,jmod,jmodpl,jmodmin,&
                   I_dim,J_dim,pos_min_dist
        REAL*8::r_item(3),angle_item,dist_item,cone_r,min_dist_up
        CHARACTER::b_item
        LOGICAL::noatom

        INTEGER*4,DIMENSION(:),ALLOCATABLE::item_to_find
        REAL*8,DIMENSION(:),ALLOCATABLE::angle_to_find,dist_to_find
        REAL*8,DIMENSION(:,:),ALLOCATABLE::r_to_find
        CHARACTER,DIMENSION(:),ALLOCATABLE::b_to_find
        LOGICAL,DIMENSION(:),ALLOCATABLE::not_found
       
        I_dim=x_cell
        J_dim=y_cell

        IF(ibot.eq.1.and.jbot.eq.1)THEN
           imodmin=I_dim;imod=ibot;imodpl=ibot+1
           jmodmin=J_dim;jmod=jbot;jmodpl=jbot+1
           WRITE(6,*) imodmin,imod,imodpl,jmodmin,jmod,jmodpl
        ELSEIF(ibot.eq.1.and.jbot.ne.J_dim.and.jbot.ne.1)THEN
            imodmin=I_dim;imod=ibot;imodpl=ibot+1
            jmodmin=jbot-1;jmod=jbot;jmodpl=jbot+1
           WRITE(6,*) imodmin,imod,imodpl,jmodmin,jmod,jmodpl
        ELSEIF(ibot.ne.1.and.jbot.eq.1.and.ibot.ne.I_dim)THEN
            imodmin=ibot-1;imod=ibot;imodpl=ibot+1
            jmodmin=J_dim;jmod=jbot;jmodpl=jbot+1
           WRITE(6,*) imodmin,imod,imodpl,jmodmin,jmod,jmodpl
        ELSEIF(ibot.eq.I_dim.and.jbot.eq.J_dim)THEN
            imodmin=I_dim-1;imod=I_dim;imodpl=1
            jmodmin=J_dim-1;jmod=J_dim;jmodpl=1
           WRITE(6,*) imodmin,imod,imodpl,jmodmin,jmod,jmodpl
        ELSEIF(ibot.eq.I_dim.and.jbot.ne.1.and.jbot.ne.J_dim)THEN
            imodmin=I_dim-1;imod=I_dim;imodpl=1
            jmodmin=jbot-1;jmod=jbot;jmodpl=jbot+1
           WRITE(6,*) imodmin,imod,imodpl,jmodmin,jmod,jmodpl
        ELSEIF(ibot.ne.I_dim.and.jbot.eq.J_dim.and.ibot.ne.1)THEN
            imodmin=ibot-1;imod=ibot;imodpl=ibot+1
            jmodmin=J_dim-1;jmod=J_dim;jmodpl=1
           WRITE(6,*) imodmin,imod,imodpl,jmodmin,jmod,jmodpl
        ELSEIF(ibot.eq.1.and.jbot.eq.J_dim)THEN
            imodmin=I_dim;imod=ibot;imodpl=ibot+1
            jmodmin=J_dim-1;jmod=J_dim;jmodpl=1
           WRITE(6,*) imodmin,imod,imodpl,jmodmin,jmod,jmodpl
        ELSEIF(ibot.eq.I_dim.and.jbot.eq.J_dim)THEN
            imodmin=I_dim-1;imod=I_dim;imodpl=1
            jmodmin=J_dim-1;jmod=J_dim;jmodpl=1
           WRITE(6,*) imodmin,imod,imodpl,jmodmin,jmod,jmodpl
        ELSE
            imodmin=ibot-1;imod=ibot;imodpl=ibot+1
            jmodmin=jbot-1;jmod=jbot;jmodpl=jbot+1
           WRITE(6,*) imodmin,imod,imodpl,jmodmin,jmod,jmodpl
        END IF

        cone_r=0.75d0*cone_h             

      WRITE(6,*)'tg',ibot,jbot,imodmin,imod,imodpl,jmodmin,jmod,jmodpl 

      100 CALL find_upper_layer_atom(imod,imodpl,jmod,jmodpl,&
             r_low,norm_low,cone_h,cone_r,item,r_item,b_item,&
             angle_item,dist_item,noatom)

        IF(noatom.and.cone_r<10.0d0) THEN
             cone_r=cone_r+0.25d0*cone_h
             GOTO 101
        END IF

        IF (noatom) THEN
        WRITE(6,*) 'No atom found at the upper-layer',latom
                    not_found(1)=noatom
        WRITE(6,*) 'Upper Hexagon',hexa_atom,'No atom found'
        ELSE
        WRITE(6,*) hexa_atom,item,r_item,b_item,angle_item,&
                   dist_item,noatom

        item_to_find(1)=item
        r_to_find(1,1)=r_item(1)
        r_to_find(1,2)=r_item(2)
        r_to_find(1,3)=r_item(3)
        b_to_find(1)=b_item
        
        angle_to_find(1)=angle_item
        dist_to_find(1)=dist_item
        
        not_found(1)=noatom
        
        WRITE(6,*) item_to_find(1),item,r_to_find(1,:),r_item(:), &
                   b_to_find(1),b_item,angle_to_find(1),angle_item,&
                   dist_to_find(1),dist_item,not_found(1),noatom

        END IF

        cone_r=0.75d0*cone_h             

      101 CALL find_upper_layer_atom(imod,imodpl,jmodmin,jmod,&
             r_low,norm_low,cone_h,cone_r,item,r_item,b_item,&
             angle_item,dist_item,noatom)

        IF(noatom.and.cone_r<10.0d0) THEN
             cone_r=cone_r+0.25d0*cone_h
             GOTO 101
        END IF

        IF (noatom) THEN
        WRITE(6,*) 'No atom found at the upper-layer',latom
                    not_found(2)=noatom
        WRITE(6,*) 'Upper Hexagon',hexa_atom,'No atom found'
        ELSE
        WRITE(6,*) hexa_atom,item,r_item,b_item,angle_item,&
                   dist_item,noatom

        item_to_find(2)=item
        r_to_find(2,1)=r_item(1)
        r_to_find(2,2)=r_item(2)
        r_to_find(2,3)=r_item(3)
        b_to_find(2)=b_item
        
        angle_to_find(2)=angle_item
        dist_to_find(2)=dist_item
        
        not_found(2)=noatom

        WRITE(6,*) item_to_find(2),item,r_to_find(2,:),r_item(:), &
                   b_to_find(2),b_item,angle_to_find(2),angle_item,&
                   dist_to_find(2),dist_item,not_found(2),noatom
        END IF

        cone_r=0.75d0*cone_h             

      102 CALL find_upper_layer_atom(imodmin,imod,jmod,jmodpl,r_low,&
             norm_low,cone_h,cone_r,item,r_item,b_item,angle_item,&
             dist_item,noatom)

        IF(noatom.and.cone_r<10.0d0) THEN
             cone_r=cone_r+0.25d0*cone_h
             GOTO 102
        END IF

        IF (noatom) THEN
        WRITE(6,*) 'No atom found at the upper-layer',latom
                    not_found(1)=noatom
        WRITE(6,*) 'Upper Hexagon',hexa_atom,'No atom found'
        ELSE
        WRITE(6,*) hexa_atom,item,r_item,b_item,angle_item,&
                   dist_item,noatom

        item_to_find(3)=item
        r_to_find(3,1)=r_item(1)
        r_to_find(3,2)=r_item(2)
        r_to_find(3,3)=r_item(3)
        b_to_find(3)=b_item
        
        angle_to_find(3)=angle_item
        dist_to_find(3)=dist_item
        
        not_found(3)=noatom

        WRITE(6,*) item_to_find(3),item,r_to_find(3,:),r_item(:), &
                   b_to_find(3),b_item,angle_to_find(3),angle_item,&
                   dist_to_find(3),dist_item,not_found(3),noatom
        END IF

        cone_r=0.75d0*cone_h             

      103 CALL find_upper_layer_atom(imodmin,imod,jmodmin,jmod,r_low,&
             norm_low,cone_h,cone_r,item,r_item,b_item,angle_item,&
             dist_item,noatom)

        IF(noatom.and.cone_r<10.0d0) THEN
             cone_r=cone_r+0.25d0*cone_h
             GOTO 103
        END IF

        IF (noatom) THEN
        WRITE(6,*) 'No atom found at the upper-layer',latom
                    not_found(4)=noatom
        WRITE(6,*) 'Upper Hexagon',hexa_atom,'No atom found'
        ELSE
        WRITE(6,*) hexa_atom,item,r_item,b_item,angle_item,&
                   dist_item,noatom

        item_to_find(4)=item
        r_to_find(4,1)=r_item(1)
        r_to_find(4,2)=r_item(2)
        r_to_find(4,3)=r_item(3)
        b_to_find(4)=b_item
        
        angle_to_find(4)=angle_item
        dist_to_find(4)=dist_item
        
        not_found(4)=noatom
        
        WRITE(6,*) item_to_find(4),item,r_to_find(4,:),r_item(:), &
                   b_to_find(4),b_item,angle_to_find(4),angle_item,&
                   dist_to_find(4),dist_item,not_found(4),noatom

        END IF

        min_dist_up=MINVAL(dist_to_find)
        pos_min_dist=MINLOC(dist_to_find,DIM=1)

        WRITE(6,*) min_dist_up, pos_min_dist

        item_up=item_to_find(pos_min_dist)

        r_up(1)=r_to_find(pos_min_dist,1)
        r_up(2)=r_to_find(pos_min_dist,2)
        r_up(3)=r_to_find(pos_min_dist,3)

        b_up=b_to_find(pos_min_dist)
        
        angle_up=angle_to_find(pos_min_dist)
        dist_up=dist_to_find(pos_min_dist)

        latom=not_found(pos_min_dist)

        WRITE(6,*) item_up,r_up(:),b_up,angle_up,dist_up,latom

        END SUBROUTINE search_atom_up_pbc
 
         !********************** 
