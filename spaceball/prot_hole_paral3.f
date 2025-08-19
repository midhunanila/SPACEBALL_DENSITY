
      program prot_hole
      implicit none
      
      integer*4 i,j,k,counter,ii,jj,kk,counter_t,counter_z,l,xw,yw,zw
      integer*4 nx,ny,nz,nres,nat,tot_at,tot,nxr,nyr,nzr,dles,iii,jjj
      integer*1, dimension(:,:,:), allocatable :: xl
      integer*4, dimension(:,:,:), allocatable :: pa
      integer*4, dimension(:,:,:), allocatable :: pr,xl2
      integer*4 omp_get_thread_num,ncell,tempt,omp_get_num_threads
      integer*4 surf,vol,edg,edg1,edg2,surfq,vt,scz,scy,scx,scxt
      integer*4 scyt,sczt,sum,imax,cl_number,rep,aver,nncell,kkk
      real*4 xmin,xmax,ymin,ymax,zmin,zmax,aa,often1
      real*4 gridx,gridy,gridz,rextend,volv,surfv,totv
      real*4 r_wall,r_water,r,r_in,vx,vy,vz,sx,sy,sz
      real*4 std_v,std_s,std_t,median_t,median_v,median_s 
      real*4, dimension(:,:), allocatable :: xp,yp,zp
      real*4, dimension(:), allocatable :: xs,ys,zs,vol_av,surf_av
      real*4, dimension(:), allocatable :: tot_av
      real*4, dimension(:,:), allocatable :: vrad
      character*3, dimension(:), allocatable :: seq
      character*2, dimension(:,:), allocatable :: at
      integer*4, dimension(:), allocatable :: na
      character*80 pdb_name,folder,d1,d2,arg,desc,text1,text2,text3
      character*4 v_name
      integer*1 machine
      integer*2 mark,mark_temp,mark1,maxval,a,b,mark_temp2
      integer*4, dimension(:,:), allocatable :: cl
      integer i_arg,rr,often
!      logical test

      do i_arg=1,iargc()
         call getarg(i_arg,arg)
      enddo
                     
      rr=1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      r_wall=1.42 !rs 1.4
      r_water=1.42 !rw 1.4
      r_in=1.0 !bylo 15.0
      gridx=0.2 !0.6
      gridy=0.2 ! 0.6
      gridz=0.2 ! 0.6
      rextend=0.0 !wall ball radius extension. Don't change this parameter!
      dles=1 !grid division
      cl_number=5 !number of clusters written into file <=85
      machine=1 !1 for PC, 2 for CLUSTER
      aver=25 !average over N values
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      

      open(rr,file=trim(arg),status='old')
11    read(rr,'(a)',end=12, err=12) desc
       if (desc(1:14) .eq. 'pdb_structure:' ) then
        read(desc(15: ), *) pdb_name
       endif
       if (desc(1:22) .eq. 'wall_probe_radius_[A]:' ) then
        read(desc(23: ), *) r_wall
       endif
       if (desc(1:22) .eq. 'wall_probe_radius_[A]:' ) then
        read(desc(23: ), *) r_wall
       endif
       if (desc(1:23) .eq. 'water_probe_radius_[A]:' ) then
        read(desc(24: ), *) r_water
       endif
       if (desc(1:11) .eq. 'grid_X_[A]:' ) then
        read(desc(12: ), *) gridx
       endif
       if (desc(1:11) .eq. 'grid_Y_[A]:' ) then
        read(desc(12: ), *) gridy
       endif
       if (desc(1:11) .eq. 'grid_Z_[A]:' ) then
        read(desc(12: ), *) gridz
       endif
       if (desc(1:41) .eq. 'number_of_clusters_written_to_the_output:' ) 
     &  then
        read(desc(42: ), *) cl_number
       endif
       if (desc(1:25) .eq. 'machine_[1_PC/2_CLUSTER]:' ) then
        read(desc(26: ), *) machine
       endif
       if (desc(1:20) .eq. 'number_of_rotations:' ) then
        read(desc(21: ), *) aver
       endif
       
      goto 11
12    close(rr)      

      tempt=0 ! don't bother of this parameter
      volv=0.0
      surfv=0.0
      totv=0.0
      std_v=0.0
      std_s=0.0
      std_t=0.0
      
      nxr=int(rextend/gridx) !additional gridpoindst forbidden for the atoms 
      nyr=int(rextend/gridy) 
      nzr=int(rextend/gridz)
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
!      do i_arg=1,iargc()
!         call getarg(i_arg,pdb_name)
!      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                           
      
      v_name=trim(pdb_name)
      !folder=v_name//'_clusters'
      folder='clusters'
      open(86, file=v_name//'.out',form='formatted')
      open(87, file=v_name//'_water.pdb',form='formatted')
      open(88, file=v_name//'_surface.pdb',form='formatted')
      !open(48, file=v_name//'_surfacerrr.pdb',form='formatted') 
      !open(49, file=v_name//'_surfaceDUPA.pdb',form='formatted')      
!      open(90, file=v_name//'_shell.pdb',form='formatted')
      open(91, file=v_name//'_x_cross-section.pdb',form='formatted')
      open(92, file=v_name//'_y_cross-section.pdb',form='formatted')              
      open(93, file=v_name//'_z_cross-section.pdb',form='formatted')
      open(94, file=v_name//'_xt_cross-section.pdb',form='formatted')
      open(95, file=v_name//'_yt_cross-section.pdb',form='formatted')              
      open(96, file=v_name//'_zt_cross-section.pdb',form='formatted')
      open(97, file=v_name//'_structure_test.pdb',form='formatted')
      if (machine.eq.2) then
       open(98, file=v_name//'_screen.txt',form='formatted')
      endif
      
      call system ('rm -r '//folder)
      
      call system ('mkdir '//folder)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
     
      

      
      allocate(vol_av(aver))
      allocate(surf_av(aver))
      allocate(tot_av(aver))
      
      do rep=1,aver !(B)
      
      !write(*,*)'Zaczalem',rep
      
      call get_nres(trim(pdb_name),nres,nat)
       
!      write(*,*)'nat',nat,nres
      
      !nres=nres-1   
          
      allocate(xp(nres,nat))
      allocate(yp(nres,nat))
      allocate(zp(nres,nat))     
      allocate(seq(nres))
      allocate(at(nres,nat))
      allocate(na(nres))
      allocate(vrad(nres,nat))
  
      call read_all_atom(trim(pdb_name),nres,nat,seq,na,xp,yp,zp,at,
     & tot_at)
      
       !do i=1,nres
       !Bwrite(*,*)i,na(i)
       !enddo
       !stop
       
       
      !do i=1,nres
      !do j=1,na(i)
      !write(*,*)'P',i,j,xp(i,j),yp(i,j),zp(i,j)
      !enddo
      !enddo
      
!      write(*,*)'nat2',nat,nres
      
      
      if (rep.gt.1) then 
      call rotation(xp,yp,zp,nres,nat,na,rep)
      endif
      
      
      
      !write(*,*)'xp: ',xp(1,1),xp(nres,1)
       
!       write(*,*)'na',nres,xp(nres,1),xp(nres-1,1),na(nres-1)
       !stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! testing correct reading
      if (rep.eq.1) then    
            counter=0
      do i=1,nres
         do j=1,na(i)
            counter=counter+1
      write(97,'(a,i7,2x,a,2x,a,i6,f12.3,2f8.3)')
     & 'ATOM',counter,at(i,j),'XXX',i,xp(i,j),yp(i,j),zp(i,j)   
         enddo
      enddo 
      close(97)
      !stop
      endif
      
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

      call assign_vdw_radii(trim(pdb_name),nres,nat,seq,na,at,vrad) 
      
      !box borders evaluation
      
      xmin=1d10
      ymin=1d10
      zmin=1d10
      xmax=-1d10
      ymax=-1d10
      zmax=-1d10
      
      do i=1,nres
         do j=1,na(i)
            if (xp(i,j).gt.xmax) xmax=xp(i,j)
            if (xp(i,j).lt.xmin) xmin=xp(i,j)
            if (yp(i,j).gt.ymax) ymax=yp(i,j)
            if (yp(i,j).lt.ymin) ymin=yp(i,j)
            if (zp(i,j).gt.zmax) zmax=zp(i,j)
            if (zp(i,j).lt.zmin) zmin=zp(i,j)
            !vrad(i,j)=vrad(i,j)*1.244455099105835d0
            !write(*,*)i,j,at(i,j),vrad(i,j)
          enddo
       enddo

         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
   
      
      r=max(r_wall,r_water)
      if (r.lt.2.20d0) r=2.20d0  !!!!!!!!!!!!!!! 2.20
                
      nx=int((xmax-xmin+8*r)/gridx)+1
      ny=int((ymax-ymin+8*r)/gridy)+1
      nz=int((zmax-zmin+8*r)/gridz)+1
      
      nx=int(nx/dles)+dles-1
      ny=int(ny/dles)+dles-1
      nz=int(nz/dles)+dles-1
      
      nxr=int(nxr/dles)+dles-1
      nyr=int(nyr/dles)+dles-1
      nzr=int(nzr/dles)+dles-1
     
      gridx=gridx*dles
      gridy=gridy*dles
      gridz=gridz*dles
     
      if (machine.eq.1) write(*,*)'Need to calculate ',nx,ny,nz,'steeps'
      if (machine.eq.2) write(98,*)'Need to calculate ',nx,ny,nz,
     & 'steeps'
      
      
      allocate(xs(nx))
      allocate(ys(ny))
      allocate(zs(nz))    
      allocate(xl(nx+nxr,ny+nyr,nz+nzr))
      allocate(pr(nx+nxr,ny+nyr,nz+nzr))
      allocate(pa(nx+nxr,ny+nyr,nz+nzr))
      
     
      do i=1,nx
         xs(i)=xmin-4*r+(i-1)*gridx
      enddo         
      do i=1,ny     
         ys(i)=ymin-4*r+(i-1)*gridy
      enddo
      do i=1,nz         
         zs(i)=zmin-4*r+(i-1)*gridz
      enddo   
      
      do i=1,nx
         do j=1,ny
            do k=1,nz
            xl(i,j,k)=0
            enddo   
         enddo
      enddo
      
      
       
       do ii=1,nres
          do jj=1,na(ii)
          !write(*,*)'Tu tez jestem',xp(ii,jj)
          xw=int((xp(ii,jj)-xmin+4*r)/gridx)+1
          yw=int((yp(ii,jj)-ymin+4*r)/gridy)+1
          zw=int((zp(ii,jj)-zmin+4*r)/gridz)+1 
          
          if (xl(xw,yw,zw).eq.0) then
              xl(xw,yw,zw)=10
              pr(xw,yw,zw)=ii
              pa(xw,yw,zw)=jj
          else
         if (machine.eq.1) write(*,*)'Two atoms at the same gridpoint!'
         if (machine.eq.2) write(98,*)'Two atoms at the same gridpoint!'
          !stop
          endif          
          enddo
       enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
! Number of cells tested around the moving direction
          
          ncell=int((r_wall+2.20)/gridx)+1 
          ncell=ncell+int(0.2*ncell)
          
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
! Number of cells tested around the moving direction - border correction
         
          nncell=int(r_wall/gridx)+1 
          !write(*,*)nncell
          nncell=nncell+int(0.2*nncell)
          !write(*,*)nncell
          !stop
          !nncell=1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!          
         
!$omp parallel
       if (machine.eq.1) write(*,*)'Processors: ',omp_get_num_threads()
       if (machine.eq.2) write(98,*)'Processors: ',omp_get_num_threads()
!$omp end parallel
        

!$omp parallel default(shared), private(i,j,k,l,ii,jj,kk,tempt)

!$omp do schedule(dynamic)         
      do j=1,ny
      if (machine.eq.1) write(*,*)'I ',j,'/',ny,omp_get_thread_num(),rep
      if (machine.eq.2) write(98,*)'I ',j,'/',ny,omp_get_thread_num(),
     & rep 
         do k=1,nz
            do i=1,nx      
                
                do ii=-ncell,ncell
             if (((j+ii).le.0).or.((j+ii).gt.ny)) goto 120
 
                   do jj=-ncell,ncell
              if (((k+jj).le.0).or.((k+jj).gt.nz)) goto 130
                   
                      do kk=0,ncell
               if (((i+kk).le.0).or.((i+kk).gt.nx)) goto 140
                      
                      if (xl(i+kk,j+ii,k+jj).gt.9) then  
            
               if (sqrt((xs(i)-xp(pr(i+kk,j+ii,k+jj),
     &          pa(i+kk,j+ii,k+jj)))**2+
     &          (ys(j)-yp(pr(i+kk,j+ii,k+jj),pa(i+kk,j+ii,k+jj)))**2+
     &          (zs(k)-zp(pr(i+kk,j+ii,k+jj),pa(i+kk,j+ii,k+jj)))**2)
     &       .le.(vrad(pr(i+kk,j+ii,k+jj),pa(i+kk,j+ii,k+jj))+r_wall)) 
     &          then
               if ((xl(i,j,k).ne.10)
     &             .and.(xl(i,j,k).ne.13)) then
                       xl(i,j,k)=1
                   else
                       xl(i,j,k)=11
                   endif
                       
                   do l=1,nxr
                   if ((xl(i+l,j,k).ne.10).and.
     &                  (xl(i+l,j,k).ne.11)) then
                       xl(i+l,j,k)=3
                   else
                       xl(i+l,j,k)=13
                   endif
                   enddo
                   
                   tempt=1
                   goto 110
               
               endif
                     endif
140                    continue
                    enddo                                     
130                    continue
                    enddo
120                     continue      
                  enddo
!               if (tempt.eq.0) xl(i,j,k)=2!xl(i,j,k)+2
                if ((tempt.eq.0).and.(xl(i,j,k).ne.1).and.
     &           (xl(i,j,k).ne.11)) xl(i,j,k)=2
               tempt=0          
            enddo
110         continue
           tempt=0   
         enddo         
      enddo
!$omp enddo nowait

!$omp do schedule(dynamic)         
      do j=1,ny
      if (machine.eq.1) write(*,*)'II ',j,'/',ny,omp_get_thread_num()
     & ,rep 
      if (machine.eq.2) write(98,*)'II ',j,'/',ny,omp_get_thread_num()
     & ,rep 
         do k=1,nz
            do i=nx,1,-1      
                
                do ii=-ncell,ncell
             if (((j+ii).le.0).or.((j+ii).gt.ny)) goto 220
 
                   do jj=-ncell,ncell
              if (((k+jj).le.0).or.((k+jj).gt.nz)) goto 230
                      
                      do kk=-ncell,0
               if (((i+kk).le.0).or.((i+kk).gt.nx)) goto 240 
                      
                      if (xl(i+kk,j+ii,k+jj).gt.9) then  
            
               if (sqrt((xs(i)-xp(pr(i+kk,j+ii,k+jj),
     &          pa(i+kk,j+ii,k+jj)))**2+
     &          (ys(j)-yp(pr(i+kk,j+ii,k+jj),pa(i+kk,j+ii,k+jj)))**2+
     &          (zs(k)-zp(pr(i+kk,j+ii,k+jj),pa(i+kk,j+ii,k+jj)))**2)
     &       .le.(vrad(pr(i+kk,j+ii,k+jj),pa(i+kk,j+ii,k+jj))+r_wall)) 
     &           then
               if ((xl(i,j,k).ne.10)
     &             .and.(xl(i,j,k).ne.13)) then
                       xl(i,j,k)=1
                   else
                       xl(i,j,k)=11
                   endif
                       
                 do l=1,nxr
                   if ((xl(i+l,j,k).ne.10).and.
     &                  (xl(i+l,j,k).ne.11)) then
                       xl(i+l,j,k)=3
                   else
                       xl(i+l,j,k)=13
                   endif
                   enddo
                   
                   
                   tempt=1
                   goto 210
               
               endif
                     endif
240                    continue
                    enddo                                     
230                    continue
                    enddo
220                     continue      
                  enddo
!               if (tempt.eq.0) xl(i,j,k)=2!xl(i,j,k)+2
                if ((tempt.eq.0).and.(xl(i,j,k).ne.1).and.
     &           (xl(i,j,k).ne.11)) xl(i,j,k)=2
               tempt=0          
            enddo
210         continue
           tempt=0   
         enddo         
      enddo
!$omp enddo nowait

!$omp do schedule(dynamic)         
      do i=1,nx
      if (machine.eq.1) write(*,*)'III ',i,'/',nx,omp_get_thread_num()
     & ,rep 
      if (machine.eq.2) write(98,*)'III ',i,'/',nx,omp_get_thread_num()
     & ,rep
      
         do k=1,nz
            do j=1,ny      
                
                do ii=-ncell,ncell
             if (((i+ii).le.0).or.((i+ii).gt.nx)) goto 320
 
                   do jj=-ncell,ncell
              if (((k+jj).le.0).or.((k+jj).gt.nz)) goto 330
               
                      do kk=0,ncell
               if (((j+kk).le.0).or.((j+kk).gt.ny)) goto 340
                      
                      if (xl(i+ii,j+kk,k+jj).gt.9) then  
            
               if (sqrt((xs(i)-xp(pr(i+ii,j+kk,k+jj),
     &          pa(i+ii,j+kk,k+jj)))**2+
     &          (ys(j)-yp(pr(i+ii,j+kk,k+jj),pa(i+ii,j+kk,k+jj)))**2+
     &          (zs(k)-zp(pr(i+ii,j+kk,k+jj),pa(i+ii,j+kk,k+jj)))**2)
     &       .le.(vrad(pr(i+ii,j+kk,k+jj),pa(i+ii,j+kk,k+jj))+r_wall)) 
     &           then
               if ((xl(i,j,k).ne.10)
     &             .and.(xl(i,j,k).ne.13)) then
                       xl(i,j,k)=1
                   else
                       xl(i,j,k)=11
                   endif
                       
                 do l=1,nxr
                   if ((xl(i,j+l,k).ne.10).and.
     &                  (xl(i,j+l,k).ne.11)) then
                       xl(i+l,j+l,k)=3
                   else
                       xl(i,j+l,k)=13
                   endif
                   enddo
                   
                   
                   tempt=1
                   goto 310
               
               endif
                     endif
340                    continue
                    enddo                                     
330                    continue
                    enddo
320                     continue      
                  enddo
!               if (tempt.eq.0) xl(i,j,k)=2!xl(i,j,k)+2
                if ((tempt.eq.0).and.(xl(i,j,k).ne.1).and.
     &           (xl(i,j,k).ne.11)) xl(i,j,k)=2
               tempt=0          
            enddo
310         continue
           tempt=0   
         enddo         
      enddo
!$omp enddo nowait

!$omp do schedule(dynamic)         
      do i=1,nx
      if (machine.eq.1) write(*,*)'IV ',i,'/',nx,omp_get_thread_num()
     & ,rep
      if (machine.eq.2) write(98,*)'IV ',i,'/',nx,omp_get_thread_num()
     & ,rep
      
         do k=1,nz
            do j=ny,1,-1      
                
                do ii=-ncell,ncell
             if (((i+ii).le.0).or.((i+ii).gt.nx)) goto 420
 
                   do jj=-ncell,ncell
              if (((k+jj).le.0).or.((k+jj).gt.nz)) goto 430
                    
                      do kk=-ncell,0
              if (((j+kk).le.0).or.((j+kk).gt.ny)) goto 440
                      
                      if (xl(i+ii,j+kk,k+jj).gt.9) then  
            
               if (sqrt((xs(i)-xp(pr(i+ii,j+kk,k+jj),
     &          pa(i+ii,j+kk,k+jj)))**2+
     &          (ys(j)-yp(pr(i+ii,j+kk,k+jj),pa(i+ii,j+kk,k+jj)))**2+
     &          (zs(k)-zp(pr(i+ii,j+kk,k+jj),pa(i+ii,j+kk,k+jj)))**2)
     &       .le.(vrad(pr(i+ii,j+kk,k+jj),pa(i+ii,j+kk,k+jj))+r_wall)) 
     &           then
                if ((xl(i,j,k).ne.10)
     &             .and.(xl(i,j,k).ne.13)) then
                       xl(i,j,k)=1
                   else
                       xl(i,j,k)=11
                   endif
                       
                 do l=1,nxr
                   if ((xl(i,j+l,k).ne.10).and.
     &                  (xl(i,j+l,k).ne.11)) then
                       xl(i,j+l,k)=3
                   else
                       xl(i,j+l,k)=13
                   endif
                   enddo
                   
                   
                   tempt=1
                   goto 410
               
               endif
                     endif
440                    continue
                    enddo                                     
430                    continue
                    enddo
420                     continue      
                  enddo
!               if (tempt.eq.0) xl(i,j,k)=2!xl(i,j,k)+2
                if ((tempt.eq.0).and.(xl(i,j,k).ne.1).and.
     &           (xl(i,j,k).ne.11)) xl(i,j,k)=2
               tempt=0          
            enddo
410         continue
           tempt=0   
         enddo         
      enddo
!$omp enddo nowait

!$omp do schedule(dynamic)         
      do i=1,nx
      if (machine.eq.1) write(*,*)'V ',i,'/',nx,omp_get_thread_num()
     & ,rep
      if (machine.eq.2) write(98,*)'V ',i,'/',nx,omp_get_thread_num()
     & ,rep
         
         do j=1,ny
            do k=1,nz      
                
                do ii=-ncell,ncell
             if (((i+ii).le.0).or.((i+ii).gt.nx)) goto 520
 
                   do jj=-ncell,ncell
              if (((j+jj).le.0).or.((j+jj).gt.ny)) goto 530
                    
                      do kk=0,ncell
              if (((k+kk).le.0).or.((k+kk).gt.nz)) goto 540
                      
                      if (xl(i+ii,j+jj,k+kk).gt.9) then  
            
               if (sqrt((xs(i)-xp(pr(i+ii,j+jj,k+kk),
     &          pa(i+ii,j+jj,k+kk)))**2+
     &          (ys(j)-yp(pr(i+ii,j+jj,k+kk),pa(i+ii,j+jj,k+kk)))**2+
     &          (zs(k)-zp(pr(i+ii,j+jj,k+kk),pa(i+ii,j+jj,k+kk)))**2)
     &       .le.(vrad(pr(i+ii,j+jj,k+kk),pa(i+ii,j+jj,k+kk))+r_wall)) 
     &           then
                if ((xl(i,j,k).ne.10)
     &             .and.(xl(i,j,k).ne.13)) then
                       xl(i,j,k)=1
                   else
                       xl(i,j,k)=11
                   endif
                       
                 do l=1,nxr
                   if ((xl(i,j,k+l).ne.10).and.
     &                  (xl(i,j,k+l).ne.11)) then
                       xl(i,j,k+l)=3
                   else
                       xl(i,j,k+l)=13
                   endif
                   enddo
                   
                   tempt=1
                   goto 510
               
               endif
                     endif
540                    continue
                    enddo                                     
530                    continue
                    enddo
520                     continue      
                  enddo
!               if (tempt.eq.0) xl(i,j,k)=2!xl(i,j,k)+2
                if ((tempt.eq.0).and.(xl(i,j,k).ne.1).and.
     &           (xl(i,j,k).ne.11)) xl(i,j,k)=2
               tempt=0          
            enddo
510         continue
           tempt=0   
         enddo         
      enddo
!$omp enddo nowait

!$omp do schedule(dynamic)         
      do i=1,nx
      if (machine.eq.1) write(*,*)'VI ',i,'/',nx,omp_get_thread_num()
     & ,rep
      if (machine.eq.2) write(98,*)'VI ',i,'/',nx,omp_get_thread_num()
     & ,rep
      
         do j=1,ny
            do k=nz,1,-1      
                
                do ii=-ncell,ncell
             if (((i+ii).le.0).or.((i+ii).gt.nx)) goto 620
 
                   do jj=-ncell,ncell
              if (((j+jj).le.0).or.((j+jj).gt.ny)) goto 630
                   
                      do kk=-ncell,0
              if (((k+kk).le.0).or.((k+kk).gt.nz)) goto 640
                      
                      if (xl(i+ii,j+jj,k+kk).gt.9) then  
            
               if (sqrt((xs(i)-xp(pr(i+ii,j+jj,k+kk),
     &          pa(i+ii,j+jj,k+kk)))**2+
     &          (ys(j)-yp(pr(i+ii,j+jj,k+kk),pa(i+ii,j+jj,k+kk)))**2+
     &          (zs(k)-zp(pr(i+ii,j+jj,k+kk),pa(i+ii,j+jj,k+kk)))**2)
     &       .le.(vrad(pr(i+ii,j+jj,k+kk),pa(i+ii,j+jj,k+kk))+r_wall)) 
     &           then
                if ((xl(i,j,k).ne.10)
     &             .and.(xl(i,j,k).ne.13)) then
                       xl(i,j,k)=1
                   else
                       xl(i,j,k)=11
                   endif
                       
                 do l=1,nxr
                   if ((xl(i,j,k+l).ne.10).and.
     &                  (xl(i,j,k+l).ne.11)) then
                       xl(i,j,k+l)=3
                   else
                       xl(i,j,k+l)=13
                   endif
                   enddo
                   
                   tempt=1
                   goto 610
               
               endif
                     endif
640                    continue
                    enddo                                     
630                    continue
                    enddo
620                     continue      
                  enddo
!               if (tempt.eq.0) xl(i,j,k)=2!xl(i,j,k)+2
                if ((tempt.eq.0).and.(xl(i,j,k).ne.1).and.
     &           (xl(i,j,k).ne.11)) xl(i,j,k)=2
               tempt=0          
            enddo
610         continue
           tempt=0   
         enddo         
      enddo
!$omp enddo nowait

!$omp do schedule(dynamic)         
      do i=1,nx
      if (machine.eq.1) write(*,*)'TOTAL ',i,'/',nx,omp_get_thread_num()
     & ,rep
      if (machine.eq.2) write(98,*)'TOTAL ',i,'/',
     & nx,omp_get_thread_num(),rep
      
         do j=1,ny
            do k=1,nz     
                
                if (xl(i,j,k).eq.0) then
                
                do ii=-ncell,ncell
             if (((i+ii).le.0).or.((i+ii).gt.nx)) goto 720
 
                   do jj=-ncell,ncell
              if (((j+jj).le.0).or.((j+jj).gt.ny)) goto 730
                      
                      do kk=-ncell,ncell
              if (((k+kk).le.0).or.((k+kk).gt.nz)) goto 740
                         
                      if (xl(i+ii,j+jj,k+kk).gt.9) then  
            
               if (sqrt((xs(i)-xp(pr(i+ii,j+jj,k+kk),
     &          pa(i+ii,j+jj,k+kk)))**2+
     &          (ys(j)-yp(pr(i+ii,j+jj,k+kk),pa(i+ii,j+jj,k+kk)))**2+
     &          (zs(k)-zp(pr(i+ii,j+jj,k+kk),pa(i+ii,j+jj,k+kk)))**2)
     &       .le.(vrad(pr(i+ii,j+jj,k+kk),pa(i+ii,j+jj,k+kk))+r_wall)) 
     &          then
!                  write(*,*)'DUPA',sqrt((xs(i)-xp(pr(i+ii,j+jj,k+kk),
!     &          pa(i+ii,j+jj,k+kk)))**2+
!     &          (ys(j)-yp(pr(i+ii,j+jj,k+kk),pa(i+ii,j+jj,k+kk)))**2+
!     &          (zs(k)-zp(pr(i+ii,j+jj,k+kk),pa(i+ii,j+jj,k+kk)))**2),
!     &          vrad(pr(i+ii,j+jj,k+kk),pa(i+ii,j+jj,k+kk))
     
                  if ((xl(i,j,k).ne.10)
     &             .and.(xl(i,j,k).ne.11)
     &             .and.(xl(i,j,k).ne.13)) then
                       xl(i,j,k)=7
                   else
                       xl(i,j,k)=17
                   endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
                   !xl(i,j,k)=7
                   !do l=1,nxr
                   !   xl(i,j,k+l)=xl(i,j,k+l)+3
                   !enddo
                   !tempt=1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                   
                   goto 710
               
               endif
                     endif
740                    continue
                    enddo                                     
730                    continue
                    enddo
720                     continue      
                  enddo
               !if (tempt.eq.0) xl(i,j,k)=2!xl(i,j,k)+2
               !tempt=0          
710         continue 
             endif           
            enddo

           !tempt=0   
         enddo         
      enddo
!$omp enddo nowait




!          ncell=int((r_in*r_wall+1.88)/gridx)+1 
!          ncell=ncell+int(0.2*ncell)



! !$omp do schedule(dynamic)         
!      do i=1,nx
!       write(*,*)'TOTAL ',i,'/',nx,omp_get_thread_num()
!         do j=1,ny
!           do k=1,nz     
                
!                if (xl(i,j,k).eq.0) then
                
!                do ii=-ncell,ncell
!             if (((i+ii).le.0).or.((i+ii).gt.nx)) goto 820
 
!                   do jj=-ncell,ncell
!              if (((j+jj).le.0).or.((j+jj).gt.ny)) goto 830
                   
!                   do kk=-ncell,ncell
!              if (((k+kk).le.0).or.((k+kk).gt.nz)) goto 840   
                     
!                      if (xl(i+ii,j+jj,k+kk).gt.9) then  
            
!               if (sqrt((xs(i)-xp(pr(i+ii,j+jj,k+kk),
!     &          pa(i+ii,j+jj,k+kk)))**2+
!     &          (ys(j)-yp(pr(i+ii,j+jj,k+kk),pa(i+ii,j+jj,k+kk)))**2+
!     &          (zs(k)-zp(pr(i+ii,j+jj,k+kk),pa(i+ii,j+jj,k+kk)))**2)
!     &       .le.(vrad(pr(i+ii,j+jj,k+kk),
!     &          pa(i+ii,j+jj,k+kk))+r_in*r_wall)) then
!                   xl(i,j,k)=8
                   !do l=1,nxr
                   !   xl(i,j,k+l)=xl(i,j,k+l)+3
                   !enddo
                   !tempt=1
!                   goto 810
               
!               endif
!                     endif
!840                    continue
!                    enddo                                     
!830                    continue
!                    enddo
!820                     continue      
!                  enddo
               !if (tempt.eq.0) xl(i,j,k)=2!xl(i,j,k)+2
               !tempt=0          
!810         continue 
!             endif           
!            enddo

           !tempt=0   
!         enddo         
!      enddo
! !$omp enddo nowait 




!$omp end parallel     
         
         surf=0
         surfq=0
         vol=0
         edg=0
         edg1=0
         edg2=0
         vt=0
         
         do i=1,nx
         if (machine.eq.1) write(*,*)'Calculations...',i,'/',nx,rep
         if (machine.eq.2) write(98,*)'Calculations...',i,'/',nx,rep
         
            do j=1,ny
               do k=1,nz
          if ((xl(i,j,k).eq.1).or.(xl(i,j,k).eq.11)) surf=surf+1 ! OK surface
          if (xl(i,j,k).eq.0) vol=vol+1 !Water volume
          !if (xl(i,j,k).eq.8) edg=edg+1
          !if ((xl(i,j,k).eq.7).or.(xl(i,j,k).eq.17)) edg1=edg1+1
          if (xl(i,j,k).eq.2) edg2=edg2+1 !OK zewnetrze
          !if ((xl(i,j,k).eq.3).or.(xl(i,j,k).eq.13)) surfq=surfq+1
          if (xl(i,j,k).ne.2) vt=vt+1 !OK Total volume      
               enddo
            enddo
         enddo      


        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      
      if (rep.eq.1) then
       scx=0
         scy=0
         scz=0
         scxt=0
         scyt=0
         sczt=0
      do j=1,ny
       if (machine.eq.1) write(*,*)'cross-section ',j,'/',ny,rep
       if (machine.eq.2) write(98,*)'cross-section ',j,'/',ny,rep
       
         do k=1,nz
            do i=1,nx       
               
               if (xl(i,j,k).eq.0) then
                 scx=scx+1
                 write(91,'(a,i7,2x,a,3x,a,i6,f12.3,2f8.3)')
     & 'ATOM',counter,'F','AAA',i,xs(1),ys(j),zs(k)
               goto 191
               endif
            enddo
191            continue      
          enddo         
      enddo
      
      do i=1,nx
       if (machine.eq.1) write(*,*)'cross-section ',i,'/',nx,rep
       if (machine.eq.2) write(98,*)'cross-section ',i,'/',nx,rep
       
         do k=1,nz
            do j=1,ny       
               
               if (xl(i,j,k).eq.0) then
                 scy=scy+1
                 write(92,'(a,i7,2x,a,3x,a,i6,f12.3,2f8.3)')
     & 'ATOM',counter,'F','AAA',i,xs(i),ys(1),zs(k)
               goto 192
               endif
            enddo
192       continue                  
          enddo         
      enddo      


      do j=1,ny
       if (machine.eq.1) write(*,*)'cross-section ',j,'/',ny,rep
       if (machine.eq.2) write(98,*)'cross-section ',j,'/',ny,rep
       
         do i=1,nx
            do k=1,nz       
               
               if (xl(i,j,k).eq.0) then
                 scz=scz+1
                 write(93,'(a,i7,2x,a,3x,a,i6,f12.3,2f8.3)')
     & 'ATOM',counter,'F','AAA',i,xs(i),ys(j),zs(1)
               goto 193
               endif
            enddo      
193        continue          
          enddo         
      enddo
      
            do j=1,ny,10
       if (machine.eq.1) write(*,*)'cross-section ',j,'/',ny,rep
       if (machine.eq.2) write(98,*)'cross-section ',j,'/',ny,rep
         
         do k=1,nz,10
            do i=1,nx,10       
               
               if (xl(i,j,k).ne.2) then
                 scxt=scxt+1
                 write(94,'(a,i7,2x,a,3x,a,i6,f12.3,2f8.3)')
     & 'ATOM',counter,'F','AAA',i,xs(1),ys(j),zs(k)
               goto 194
               endif
            enddo
194            continue      
          enddo         
      enddo
      
      do i=1,nx,10
       if (machine.eq.1) write(*,*)'cross-section ',i,'/',nx,rep
       if (machine.eq.2) write(98,*)'cross-section ',i,'/',nx,rep
       
         do k=1,nz,10
            do j=1,ny,10       
               
               if (xl(i,j,k).ne.2) then
                 scyt=scyt+1
                 write(95,'(a,i7,2x,a,3x,a,i6,f12.3,2f8.3)')
     & 'ATOM',counter,'F','AAA',i,xs(i),ys(1),zs(k)
               goto 195
               endif
            enddo
195       continue                  
          enddo         
      enddo      


      do j=1,ny,10
       if (machine.eq.1) write(*,*)'cross-section ',j,'/',ny,rep
       if (machine.eq.2) write(98,*)'cross-section ',j,'/',ny,rep
       
         do i=1,nx,10
            do k=1,nz,10       
               
               if (xl(i,j,k).ne.2) then
                 sczt=sczt+1
                 write(96,'(a,i7,2x,a,3x,a,i6,f12.3,2f8.3)')
     & 'ATOM',counter,'F','AAA',i,xs(i),ys(j),zs(1)
               goto 196
               endif
            enddo      
196        continue          
          enddo         
      enddo

    





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        vx=4/3*3.14159265*((xmax-xmin)/2)**3
        vy=4/3*3.14159265*((ymax-ymin)/2)**3
        vz=4/3*3.14159265*((zmax-zmin)/2)**3
        sx=4*3.14159265*((xmax-xmin)/2)**2
        sy=4*3.14159265*((ymax-ymin)/2)**2
        sz=4*3.14159265*((zmax-zmin)/2)**2
        
        
     
      counter=0                  
      do i=1,nx
        if (machine.eq.1) write(*,*)i,'/',nx,rep
        if (machine.eq.2) write(98,*)i,'/',nx,rep
        
        do j=1,ny
           do k=1,nz
              if ((xl(i,j,k).eq.0)) then
              !if ((xl(i,j,k).eq.2)) then
!                 if (xl(i,j,k).eq.3) then
                 counter=counter+1
      write(87,'(a,i7,2x,a,3x,a,i6,f12.3,2f8.3)')
     & 'ATOM',counter,'F','AAA',i,xs(i),ys(j),zs(k)
              endif
           enddo
        enddo
      enddo
      


      endif

!      counter=0                  
!      do i=1,nx  
!        write(*,*)i,'/',nx
!        do j=1,ny
!           do k=1,nz
!              if ((xl(i,j,k).eq.0)) then
!                 counter=counter+1
!      write(89,'(a,i7,2x,a,3x,a,i6,f12.3,2f8.3)')
!     & 'ATOM',counter,'F','AAA',i,xs(i),ys(j),zs(k)
!              endif
!           enddo
!        enddo
!      enddo
      

!SHELL WRITTING HAS TO BE CORRECTED
!      counter=0                  
!      do i=1,nx,10
!         if (machine.eq.1) write(*,*)i,'/',nx
!         if (machine.eq.2) write(98,*)i,'/',nx
       
!        do j=1,ny,10
!           do k=1,nz,10
!              if ((xl(i,j,k).eq.8).or.(xl(i,j,k).eq.7)) then
              !if ((xl(i,j,k).eq.2)) then
                 !if (xl(i,j,k).eq.3) then
!                 counter=counter+1
!      write(90,'(a,i7,2x,a,3x,a,i6,f12.3,2f8.3)')
!     & 'ATOM',counter,'F','AAA',i,xs(i),ys(j),zs(k)
!              endif
!           enddo
!        enddo
!      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! CAPSULES !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        SURFACE CORRECTION      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do i=1,nx
       do j=1,ny
        do k=1,nz
         if ((xl(i,j,k).eq.1).or.(xl(i,j,k).eq.11)) then
           do ii=-ncell,ncell
            if (((i+ii).le.0).or.((i+ii).gt.nx)) goto 922
             do jj=-ncell,ncell
              if (((j+jj).le.0).or.((j+jj).gt.ny)) goto 932            
               do kk=-ncell,ncell
                if (((k+kk).le.0).or.((k+kk).gt.nz)) goto 942
            
           if ((xl(i+ii,j+jj,k+kk).ne.2).and.(xl(i+ii,j+jj,k+kk).lt.9)
!     $   .and.(xl(i+ii,j+jj,k+kk).ne.11).and.(xl(i+ii,j+jj,k+kk).ne.1)
     $     .and.(xl(i+ii,j+jj,k+kk).ne.6)) 
     $     then
            if (sqrt((xs(i)-xs(i+ii))**2+(ys(j)-ys(j+jj))**2+
     $       (zs(k)-zs(k+kk))**2).le.r_wall) then

            xl(i+ii,j+jj,k+kk)=6
     
           endif             
            endif
                 
            
942              continue
             enddo
932              continue             
            enddo
922              continue            
           enddo
         endif    
        enddo
       enddo
      enddo
      
!            if (rep.eq.1) then

!       counter=0                  
!      do i=1,nx
!        if (machine.eq.1) write(*,*)i,'/',nx
!        if (machine.eq.2) write(98,*)i,'/',nx
!        do j=1,ny
!           do k=1,nz
!              if ((xl(i,j,k).eq.6).or.(xl(i,j,k).eq.6)) then
!                 counter=counter+1
!              endif
!           enddo
!        enddo
!      enddo
      
!      often1=counter/1000000
!      often1=often1**(1./3.)
!      often=often1
!      often=often+1
      
!      counter=0                  
!      do i=1,nx,often
!        if (machine.eq.1) write(*,*)i,'/',nx
!        if (machine.eq.2) write(98,*)i,'/',nx
!        do j=1,ny,often
!           do k=1,nz,often
!              if ((xl(i,j,k).eq.6).or.(xl(i,j,k).eq.6)) then
              !if ((xl(i,j,k).eq.2)) then
!                 if (xl(i,j,k).eq.3) then
!                 counter=counter+1
!      write(49,'(a,i7,2x,a,3x,a,i6,f12.3,2f8.3)')
!     & 'ATOM',counter,'S','AAA',i,xs(i),ys(j),zs(k)
!              endif
!           enddo
!        enddo
!      enddo

!      endif      
      

      
 
      ! ncell=2*ncell             
            
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!$omp parallel default(shared), private(i,j,k,l,ii,jj,kk,tempt)

!$omp do schedule(dynamic)         
      do j=1,ny
      if (machine.eq.1) write(*,*)'7 ',j,'/',ny,omp_get_thread_num(),rep
      if (machine.eq.2) write(98,*)'7 ',j,'/',ny,omp_get_thread_num(),
     & rep 
         do k=1,nz
            do i=1,nx      
               if (xl(i,j,k).eq.6) then 
                do ii=-ncell,ncell
             if (((j+ii).le.0).or.((j+ii).gt.ny)) goto 122
 
                   do jj=-ncell,ncell
              if (((k+jj).le.0).or.((k+jj).gt.nz)) goto 132
                   
                      do kk=-ncell,ncell
               if (((i+kk).le.0).or.((i+kk).gt.nx)) goto 142
                      
                      if (xl(i+kk,j+ii,k+jj).gt.9) then  
            
               if (sqrt((xs(i)-xp(pr(i+kk,j+ii,k+jj),
     &          pa(i+kk,j+ii,k+jj)))**2+
     &          (ys(j)-yp(pr(i+kk,j+ii,k+jj),pa(i+kk,j+ii,k+jj)))**2+
     &          (zs(k)-zp(pr(i+kk,j+ii,k+jj),pa(i+kk,j+ii,k+jj)))**2)
     &       .le.(vrad(pr(i+kk,j+ii,k+jj),pa(i+kk,j+ii,k+jj)))) 
     &          then
                       
                     xl(i,j,k)=7
                      !if (xl(i-1,j,k).gt.9) xl(i-1,j,k)=11  
                   
                   
                   !tempt=1
                   goto 112
               
               endif
                     endif
142                    continue
                    enddo                                     
132                    continue
                    enddo
122                     continue      
                  enddo
!               if (tempt.eq.0) xl(i,j,k)=2!xl(i,j,k)+2
!                if ((tempt.eq.0).and.(xl(i,j,k).ne.1).and.
!     &           (xl(i,j,k).ne.11)) xl(i,j,k)=2
               !tempt=0
               endif          
            enddo
112         continue
           !tempt=0   
         enddo         
      enddo
!$omp enddo nowait
!$omp end parallel


!            if (rep.eq.1) then

!B       counter=0                  
!      do i=1,nx
!        if (machine.eq.1) write(*,*)i,'/',nx
!        if (machine.eq.2) write(98,*)i,'/',nx
!        do j=1,ny
!           do k=1,nz
!              if ((xl(i,j,k).eq.7).or.(xl(i,j,k).eq.7)) then
!                 counter=counter+1
!              endif
!           enddo
!        enddo
!      enddo
      
!      often1=counter/1000000
!      often1=often1**(1./3.)
!      often=often1
!      often=often+1
      
!      counter=0                  
!      do i=1,nx,often
!        if (machine.eq.1) write(*,*)i,'/',nx
!        if (machine.eq.2) write(98,*)i,'/',nx
!        do j=1,ny,often
!           do k=1,nz,often
!              if ((xl(i,j,k).eq.7).or.(xl(i,j,k).eq.7)) then
              !if ((xl(i,j,k).eq.2)) then
!                 if (xl(i,j,k).eq.3) then
!                 counter=counter+1
!      write(48,'(a,i7,2x,a,3x,a,i6,f12.3,2f8.3)')
!     & 'ATOM',counter,'R','AAA',i,xs(i),ys(j),zs(k)
!              endif
!           enddo
!        enddo
!      enddo

!      endif


      do i=1,nx
       do j=1,ny
        do k=1,nz
             
           if (xl(i,j,k).eq.6) then
      
      if ((xl(i+1,j,k).eq.0).or.(xl(i+1,j,k).eq.7).or.
     $     (xl(i-1,j,k).eq.0).or.(xl(i-1,j,k).eq.7).or.
     $     (xl(i,j+1,k).eq.0).or.(xl(i,j+1,k).eq.7).or.
     $     (xl(i,j-1,k).eq.0).or.(xl(i,j-1,k).eq.7).or.
     $     (xl(i,j,k+1).eq.0).or.(xl(i,j,k+1).eq.7).or.
     $     (xl(i,j,k-1).eq.0).or.(xl(i,j,k-1).eq.7)) then
      
           xl(i,j,k)=1
           else
           xl(i,j,k)=6
       endif              
       endif           
           
!           do ii=-1,1
!            if (((i+ii).le.0).or.((i+ii).gt.nx)) goto 9221
!             do jj=-1,1
!              if (((j+jj).le.0).or.((j+jj).gt.ny)) goto 9321            
!               do kk=-1,1
!                if (((k+kk).le.0).or.((k+kk).gt.nz)) goto 9421
            
!           if ((xl(i+ii,j+jj,k+kk).ne.6).and.
!     $           (xl(i+ii,j+jj,k+kk).ne.2).and.
!     $           (xl(i+ii,j+jj,k+kk).ne.1)) then
!             xl(i,j,k)=1
!             goto 9222
!             else
!             xl(i,j,k)=6
!            endif
                 
            
!9421              continue
!             enddo
!9321              continue             
!            enddo
!9221              continue            
!           enddo           
!         endif
!9222              continue           
        enddo        
       enddo
      enddo



            if (rep.eq.1) then

       counter=0                  
      do i=1,nx
        if (machine.eq.1) write(*,*)i,'/',nx
        if (machine.eq.2) write(98,*)i,'/',nx
        do j=1,ny
           do k=1,nz
              if ((xl(i,j,k).eq.1).or.(xl(i,j,k).eq.11)) then
                 counter=counter+1
              endif
           enddo
        enddo
      enddo
      
      often1=counter/1000000
      often1=often1**(1./3.)
      often=often1
      often=often+1
      
      counter=0                  
      do i=1,nx,often
        if (machine.eq.1) write(*,*)i,'/',nx
        if (machine.eq.2) write(98,*)i,'/',nx
        do j=1,ny,often
           do k=1,nz,often
              if ((xl(i,j,k).eq.11).or.(xl(i,j,k).eq.1)) then
              !if ((xl(i,j,k).eq.2)) then
!                 if (xl(i,j,k).eq.3) then
                 counter=counter+1
      write(88,'(a,i7,2x,a,3x,a,i6,f12.3,2f8.3)')
     & 'ATOM',counter,'F','AAA',i,xs(i),ys(j),zs(k)
              endif
           enddo
        enddo
      enddo

      endif      
      







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do i=1,nx
       do j=1,ny
        do k=1,nz
         if (xl(i,j,k).eq.0) then
           do ii=-nncell,nncell
            if (((i+ii).le.0).or.((i+ii).gt.nx)) goto 920
             do jj=-nncell,nncell
              if (((j+jj).le.0).or.((j+jj).gt.ny)) goto 930            
               do kk=-nncell,nncell
                if (((k+kk).le.0).or.((k+kk).gt.nz)) goto 940
            
           if ((xl(i+ii,j+jj,k+kk).ne.0).and.(xl(i+ii,j+jj,k+kk).lt.9)
     $     .and.(xl(i+ii,j+jj,k+kk).ne.5)) 
     $     then
            if (sqrt((xs(i)-xs(i+ii))**2+(ys(j)-ys(j+jj))**2+
     $       (zs(k)-zs(k+kk))**2).le.r_wall) then

            xl(i+ii,j+jj,k+kk)=5
     
           endif             
            endif
                 
            
940              continue
             enddo
930              continue             
            enddo
920              continue            
           enddo
         endif    
        enddo
       enddo
      enddo
      



!$omp parallel default(shared), private(i,j,k,l,ii,jj,kk,tempt)
!$omp do schedule(dynamic)         
      do i=1,nx
      if (machine.eq.1) write(*,*)'TOT2 ',i,'/',nx,omp_get_thread_num()
     & ,rep
      if (machine.eq.2) write(98,*)'TOT2 ',i,'/',
     & nx,omp_get_thread_num(),rep
      
         do j=1,ny
            do k=1,nz     
                
                if (xl(i,j,k).eq.5) then
                
                do ii=-ncell,ncell
             if (((i+ii).le.0).or.((i+ii).gt.nx)) goto 820
 
                   do jj=-ncell,ncell
              if (((j+jj).le.0).or.((j+jj).gt.ny)) goto 830
                      
                      do kk=-ncell,ncell
              if (((k+kk).le.0).or.((k+kk).gt.nz)) goto 840
                         
                      if (xl(i+ii,j+jj,k+kk).gt.9) then  
            
               if (sqrt((xs(i)-xp(pr(i+ii,j+jj,k+kk),
     &          pa(i+ii,j+jj,k+kk)))**2+
     &          (ys(j)-yp(pr(i+ii,j+jj,k+kk),pa(i+ii,j+jj,k+kk)))**2+
     &          (zs(k)-zp(pr(i+ii,j+jj,k+kk),pa(i+ii,j+jj,k+kk)))**2)
     &       .gt.(vrad(pr(i+ii,j+jj,k+kk),pa(i+ii,j+jj,k+kk)))) 
     &          then
!                  write(*,*)'DUPA',sqrt((xs(i)-xp(pr(i+ii,j+jj,k+kk),
!     &          pa(i+ii,j+jj,k+kk)))**2+
!     &          (ys(j)-yp(pr(i+ii,j+jj,k+kk),pa(i+ii,j+jj,k+kk)))**2+
!     &          (zs(k)-zp(pr(i+ii,j+jj,k+kk),pa(i+ii,j+jj,k+kk)))**2),
!     &          vrad(pr(i+ii,j+jj,k+kk),pa(i+ii,j+jj,k+kk))
     
!                  if ((xl(i,j,k).ne.10)
!     &             .and.(xl(i,j,k).ne.11)
!     &             .and.(xl(i,j,k).ne.13)) then
!                       xl(i,j,k)=7
!                   else
                       xl(i,j,k)=0
!                   endif


                   goto 810
               
               endif
                     endif
840                    continue
                    enddo                                     
830                    continue
                    enddo
820                     continue      
                  enddo
               !if (tempt.eq.0) xl(i,j,k)=2!xl(i,j,k)+2
               !tempt=0          
810         continue 
             endif           
            enddo

           !tempt=0   
         enddo         
      enddo
!$omp enddo nowait
!$omp end parallel     
!!!!!!!!!!! THE LAST PART !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
      mark=1  
       
      do i=1,nx
       do j=1,ny
        do k=1,nz
         if (xl(i,j,k).eq.0) xl(i,j,k)=9
         if ((xl(i,j,k).ne.0).and.(xl(i,j,k).ne.9)) xl(i,j,k)=4    
        enddo
       enddo
      enddo 
      
      do i=1,nx
       do j=1,ny
        do k=1,nz
         if (xl(i,j,k).eq.9) xl(i,j,k)=1
         if (xl(i,j,k).eq.4) xl(i,j,k)=0      
        enddo
       enddo
      enddo
      
      open(99, file='temporary_file',form='formatted')

      
      do i=1,nx
       if (machine.eq.1) write(*,*)'Writing down',i,'/',nx,rep
       if (machine.eq.2) write(98,*)'Writing down',i,'/',nx,rep
       if (machine.eq.2) flush(98)
       do j=1,ny
        do k=1,nz
         write(99,*) xl(i,j,k)      
        enddo
       enddo
      enddo
      
      close(99)
      deallocate(xl)
      deallocate(pa)
      deallocate(pr)
!      deallocate(xs)
!      deallocate(ys)
!      deallocate(zs)
      deallocate(xp)
      deallocate(yp)
      deallocate(zp)
      deallocate(vrad)
      deallocate(na)
      allocate(xl2(nx,ny,nz))
      open(99, file='temporary_file',form='formatted')
       
      do i=1,nx
       if (machine.eq.1) write(*,*)'Reading',i,'/',nx,rep
       if (machine.eq.2) write(98,*)'Reading',i,'/',nx,rep
       if (machine.eq.2) flush(98) 
       do j=1,ny
        do k=1,nz
         read(99,*) xl2(i,j,k)      
        enddo
       enddo
      enddo
      
      
      call system ('rm '//'temporary_file')
      
      
      do i=1,nx
       if(machine.eq.1) write(*,*)'Changing the table: ',i,'/',nx,rep
       if(machine.eq.2) write(98,*)'Changing the table: ',i,'/',nx,rep
       if (machine.eq.2) flush(98)
       do j=1,ny
        do k=1,nz
         if (xl2(i,j,k).eq.1) then
          if (xl2(i+1,j,k).gt.1) mark1=xl2(i+1,j,k)
          if (xl2(i-1,j,k).gt.1) mark1=xl2(i-1,j,k)    
          if (xl2(i,j-1,k).gt.1) mark1=xl2(i,j-1,k)
          if (xl2(i,j+1,k).gt.1) mark1=xl2(i,j+1,k)
          if (xl2(i,j,k-1).gt.1) mark1=xl2(i,j,k-1)
          if (xl2(i,j,k+1).gt.1) mark1=xl2(i,j,k+1)   
           if ((xl2(i+1,j,k).le.1).and.(xl2(i-1,j,k).le.1).and.    
     &        (xl2(i,j-1,k).le.1).and.(xl2(i,j+1,k).le.1).and.
     &        (xl2(i,j,k-1).le.1).and.(xl2(i,j,k+1).le.1)) then  
              mark=mark+1
              mark1=mark
           endif
             
             xl2(i,j,k)=mark1
          
          if (xl2(i+1,j,k).eq.1) xl2(i+1,j,k)=mark1
          
          if (xl2(i-1,j,k).eq.1) xl2(i-1,j,k)=mark1
          
          if (xl2(i,j+1,k).eq.1) xl2(i,j+1,k)=mark1
          
          if (xl2(i,j-1,k).eq.1) xl2(i,j-1,k)=mark1
          
          if (xl2(i,j,k-1).eq.1) xl2(i,j,k-1)=mark1
          
          if (xl2(i,j,k+1).eq.1) xl2(i,j,k+1)=mark1
         
         endif      
        enddo
       enddo
      enddo 
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                  THE LAST CONNECTING




      do i=1,nx
      if (machine.eq.1) write(*,*)'The last corrections: ',i,'/',nx,rep
      if (machine.eq.2) write(98,*)'The last corrections: ',i,'/',nx,
     & rep 
      if (machine.eq.2) flush(98) 
       do j=1,ny
        do k=1,nz
            
             if (xl2(i,j,k).gt.1) then
             
             mark_temp=xl2(i,j,k)
          
          if ((xl2(i+1,j,k).gt.1).and.(xl2(i+1,j,k).ne.mark_temp)) then
              mark_temp2=xl2(i+1,j,k)   

!$omp parallel default(shared),private(ii,jj,kk)
!$omp do schedule(dynamic)              
              do ii=1,nx
               do jj=1,ny
               do kk=1,nz
               if (xl2(ii,jj,kk).eq.mark_temp2) xl2(ii,jj,kk)=mark_temp
                enddo
               enddo
              enddo   
!$omp enddo nowait
!$omp end parallel
          endif
          
          if ((xl2(i-1,j,k).gt.1).and.(xl2(i-1,j,k).ne.mark_temp)) then
              mark_temp2=xl2(i-1,j,k)   
!$omp parallel default(shared),private(ii,jj,kk)
!$omp do schedule(dynamic)
              do ii=1,nx
               do jj=1,ny
                do kk=1,nz
               if (xl2(ii,jj,kk).eq.mark_temp2) xl2(ii,jj,kk)=mark_temp
                enddo
               enddo
              enddo
!$omp enddo nowait
!$omp end parallel
          endif
          
          if ((xl2(i,j+1,k).gt.1).and.(xl2(i,j+1,k).ne.mark_temp)) then
              mark_temp2=xl2(i,j+1,k)   
!$omp parallel default(shared),private(ii,jj,kk)
!$omp do schedule(dynamic) 
              do ii=1,nx
               do jj=1,ny
                do kk=1,nz
               if (xl2(ii,jj,kk).eq.mark_temp2) xl2(ii,jj,kk)=mark_temp
                enddo
               enddo
              enddo
!$omp enddo nowait
!$omp end parallel
          endif
          
          if ((xl2(i,j-1,k).gt.1).and.(xl2(i,j-1,k).ne.mark_temp)) then
              mark_temp2=xl2(i,j-1,k)   
!$omp parallel default(shared),private(ii,jj,kk)
!$omp do schedule(dynamic) 
              do ii=1,nx
               do jj=1,ny
                 do kk=1,nz
               if (xl2(ii,jj,kk).eq.mark_temp2) xl2(ii,jj,kk)=mark_temp
                enddo
               enddo
              enddo
!$omp enddo nowait
!$omp end parallel
          endif
          
          if ((xl2(i,j,k-1).gt.1).and.(xl2(i,j,k-1).ne.mark_temp)) then
              mark_temp2=xl2(i,j,k-1)   
!$omp parallel default(shared),private(ii,jj,kk)
!$omp do schedule(dynamic) 
              do ii=1,nx
               do jj=1,ny
                do kk=1,nz
               if (xl2(ii,jj,kk).eq.mark_temp2) xl2(ii,jj,kk)=mark_temp
                enddo
               enddo
              enddo
!$omp enddo nowait
!$omp end parallel
          endif
          
          if ((xl2(i,j,k+1).gt.1).and.(xl2(i,j,k+1).ne.mark_temp)) then
             mark_temp2=xl2(i,j,k+1)   
!$omp parallel default(shared),private(ii,jj,kk)
!$omp do schedule(dynamic) 
              do ii=1,nx
               do jj=1,ny
                do kk=1,nz
               if (xl2(ii,jj,kk).eq.mark_temp2) xl2(ii,jj,kk)=mark_temp
                enddo
               enddo
              enddo
!$omp enddo nowait                
!$omp end parallel
          endif
       
         endif
        enddo
       enddo   
      enddo

      
 
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
     
      allocate(cl(mark,2))
      
      sum=0
      maxval=-100
      
      do i=1,mark
      cl(i,1)=0
      cl(i,2)=0
      enddo
     
!$omp parallel default(shared),private(i,sum,ii,jj,kk)
!$omp do schedule(dynamic)
      do i=1,mark
      if (machine.eq.1) write(*,*)i,'/',mark,rep
      if (machine.eq.2) write(98,*)i,'/',mark,rep
      if (machine.eq.2) flush(98)
      do ii=1,nx
       do jj=1,ny
        do kk=1,nz
         if (xl2(ii,jj,kk).eq.i) sum=sum+1
        enddo
       enddo
      enddo      
      cl(i,1)=i
      cl(i,2)=sum
      sum=0
      enddo
!$omp enddo nowait      
!$omp end parallel       
      
      do ii=mark,2,-1
      if (machine.eq.1) write(*,*)ii,'/',mark,rep
      if (machine.eq.2) write(98,*)ii,'/',mark,rep
      if (machine.eq.2) flush(98)
       do jj=1,ii-1
        if (cl(jj,2).lt.cl(jj+1,2)) then
         a=cl(jj,2)
         b=cl(jj,1)
         cl(jj,2)=cl(jj+1,2)
         cl(jj,1)=cl(jj+1,1)
         cl(jj+1,1)=b
         cl(jj+1,2)=a
        endif
       enddo
      enddo
      
      
      
      if (rep.eq.1) then !(A) 
     
      if (cl_number.le.mark) then
      
      do ii=1,cl_number
       write(d1,'(i8)')ii
       write(d2,'(i8)')cl(ii,2)
       d1=adjustl(d1)
       d2=adjustl(d2)
       open(ii, file='./'//trim(folder)//'/'//trim(d1)//'_'//trim(d2)//
     &  '.pdb',form='formatted')
!       open(ii, file=trim(d1)//'cav'//
!     &  '.pdb',form='formatted')
       
       counter=0                  
      do i=1,nx
        if(machine.eq.1) write(*,*)i,'/',nx,ii,'/',cl_number,rep
        if(machine.eq.2) write(98,*)i,'/',nx,ii,'/',cl_number,rep
        if (machine.eq.2) flush(98)
        do j=1,ny
           do k=1,nz
              if (xl2(i,j,k).eq.cl(ii,1)) then
             counter=counter+1
!      write(ii,'(a,i7,2x,a,3x,a,i6,f12.3,2f8.3)')
!     & 'ATOM',counter,'F','AAA',i,xs(i),ys(j),zs(k)
              endif
           enddo
        enddo
      enddo
      
      often1=counter/1000000
      often1=often1**(1./3.)
      often=often1
      often=often+1
      
       counter=0                  
      do i=1,nx,often
        if(machine.eq.1) write(*,*)i,'/',nx,ii,'/',cl_number,rep
        if(machine.eq.2) write(98,*)i,'/',nx,ii,'/',cl_number,rep
        if (machine.eq.2) flush(98)
        do j=1,ny,often
           do k=1,nz,often
              if (xl2(i,j,k).eq.cl(ii,1)) then
             counter=counter+1
      write(ii,'(a,i7,2x,a,3x,a,i6,f12.3,2f8.3)')
     & 'ATOM',counter,'F','AAA',i,xs(i),ys(j),zs(k)
              endif
           enddo
        enddo
      enddo       
       close(ii)
      enddo
      
      else
      
      do ii=1,mark
       write(d1,'(i8)')ii
       write(d2,'(i8)')cl(ii,2)
       d1=adjustl(d1)
       d2=adjustl(d2)
       open(ii, file='./'//trim(folder)//'/'//trim(d1)//'_'//trim(d2)//
     &  '.pdb',form='formatted')
       
       counter=0                  
      do i=1,nx
        if (machine.eq.1) write(*,*)i,'/',nx,ii,'/',cl_number,rep
        if (machine.eq.2) write(98,*)i,'/',nx,ii,'/',cl_number,rep
        if (machine.eq.2) flush(98)
        do j=1,ny
           do k=1,nz
              if (xl2(i,j,k).eq.cl(ii,1)) then
             counter=counter+1
      write(ii,'(a,i7,2x,a,3x,a,i6,f12.3,2f8.3)')
     & 'ATOM',counter,'F','AAA',i,xs(i),ys(j),zs(k)
              endif
           enddo
        enddo
      enddo 
       close(ii)
      enddo
      
      endif      
      
         

!        write(86,*)'Cavity volume: ',vol,'[pt]',vol*gridx*gridy*gridz,
!     &  '[A^3]'
        write(86,*)'PDB_file_name: ',trim(pdb_name)
        write(86,'(a24,1x,f5.2,1x,a1)')'The radius of the probe:', 
     $   r_wall,'A'
        !write(86,*)'Radius of a probe of water', r_water,'A'
        write(86,'(a33,1x,f5.2,1x,a1)')
     $   'The grid size in the x direction: ',gridx,'A'
        write(86,'(a33,1x,f5.2,1x,a1)')
     $   'The grid size in the y direction: ',gridy,'A'
        write(86,'(a33,1x,f5.2,1x,a1)')
     $   'The grid size in the z direction: ',gridz,'A'
        
!        write(86,*)'Volume avaliable for water: ',vol,'[pt]',
!     &   (vol)*gridx*gridy*gridz,'[A^3]'
!        write(86,*)'Shell volume: ',nx*ny*nz-edg2-cl(1,2),'[pt]',
!     &   (nx*ny*nz-edg2-cl(1,2))*gridx*gridy*gridz,'[A^3]'
!        write(86,*)'To powinno byc zero, jesli nie powiekszam: ',surfq
        write(86,'(a37,1x,f8.2,1x,a1)')
     $   'Maximum dimension in the x direction: '
     $   ,xmax-xmin,'A' 
        write(86,'(a37,1x,f8.2,1x,a1)')
     $   'Maximum dimension in the y direction: '
     $   ,ymax-ymin,'A' 
        write(86,'(a37,1x,f8.2,1x,a1)')
     $   'Maximum dimension in the z direction: '
     $   ,zmax-zmin,'A'
!        write(86,*)'Volume evaluated with rx: ',vx,'A^3'
!        write(86,*)'Volume evaluated with ry: ',vy,'A^3'
!        write(86,*)'Volume evaluated with rz: ',vz,'A^3'   
!        write(86,*)'Real volume: ',vt*gridx*gridy*gridz,'A^3'
!        write(86,*)'Surface evaluated with rx: ',sx,'A^2'
!        write(86,*)'Surface evaluated with ry: ',sy,'A^2'
!        write(86,*)'Surface evaluated with rz: ',sz,'A^2' 
!        write(86,*)'Real surface: ',surf,'[pt]',surf*gridx*gridx,'[A^2]'
!        write(86,*)'x cross-section surface: ',scx,'[pt]',
!     &   scx*gridy*gridz,'[A^2]'
!        write(86,*)'y cross-section surface: ',scy,'[pt]',
!     &   scy*gridy*gridx,'[A^2]'
!        write(86,*)'z cross-section surface: ',scz,'[pt]',
!     &   scz*gridy*gridx,'[A^2]'
!        write(86,*)'x total cross-section surface: ',scxt,'[pt]',
!     &   scxt*gridy*gridz,'[A^2]',3.14159265*((xmax-xmin)/2)**2
!        write(86,*)'y total cross-section surface: ',scyt,'[pt]',
!     &   scyt*gridy*gridx,'[A^2]',3.14159265*((ymax-ymin)/2)**2
!        write(86,*)'z total cross-section surface: ',sczt,'[pt]',
!     &   sczt*gridy*gridx,'[A^2]',3.14159265*((zmax-zmin)/2)**2
        write(86,*)'##### of cavities have been detected'
        write(86,'(a29,1x,i3,1x,a13)')
     $   'The results are averaged over' ,aver, ' of rotations.'
        text1='The other enclosed file, contains the coordinates'
        text2=' of points (in the PDB format) belonging to the' 
        text3=' largest chamber.'
        write(86,*)trim(text1),trim(text2),trim(text3)
        write(86,*)'=========================================='
        write(86,*)'|          Clusters volumes              |'
        write(86,*)'=========================================='
        write(86,*)'1',cl(1,2),'[pt]',cl(1,2)*gridx*gridy*gridz,'[A^3]'
        write(86,*)'------------------------------------------'
        do ii=2,mark
         if (cl(ii,2).gt.0) write(86,*)ii,cl(ii,2),'[pt]',
     &    cl(ii,2)*gridx*gridy*gridz,'[A^3]'
        enddo
       
       endif !(A)
       
       !write(*,*)'TO przeszedlem'
       
       vol_av(rep)=cl(1,2)*gridx*gridy*gridz
       !write(*,*)cl(1,2)*gridx*gridy*gridz,vol_av(rep)
       volv=volv+cl(1,2)*gridx*gridy*gridz
       surf_av(rep)=surf*gridx*gridx
       surfv=surfv+surf*gridx*gridx
       tot_av(rep)=vt*gridx*gridy*gridz
       totv=totv+vt*gridx*gridy*gridz
       deallocate(seq)
       deallocate(at)
       deallocate(xl2)
       deallocate(cl)
       deallocate(xs)
       deallocate(ys)
       deallocate(zs)
       enddo !(B)
       volv=volv/real(aver)
       surfv=surfv/real(aver)
       totv=totv/real(aver)  
       
       do ii=aver,2,-1
       do jj=1,ii-1
        if (vol_av(jj).gt.vol_av(jj+1)) then
         aa=vol_av(jj)
         vol_av(jj)=vol_av(jj+1)
         vol_av(jj+1)=aa
        endif
       enddo
      enddo
       do ii=aver,2,-1
       do jj=1,ii-1
        if (surf_av(jj).gt.surf_av(jj+1)) then
         aa=surf_av(jj)
         surf_av(jj)=surf_av(jj+1)
         surf_av(jj+1)=aa
        endif
       enddo
      enddo

       do ii=aver,2,-1
       do jj=1,ii-1
        if (tot_av(jj).gt.tot_av(jj+1)) then
         aa=tot_av(jj)
         tot_av(jj)=tot_av(jj+1)
         tot_av(jj+1)=aa
        endif
       enddo
      enddo
      
        do ii=1,aver
        std_v=std_v+(vol_av(ii)-volv)*(vol_av(ii)-volv)
        enddo
        do ii=1,aver
        std_s=std_s+(surf_av(ii)-surfv)*(surf_av(ii)-surfv)
        enddo
        do ii=1,aver
        std_t=std_t+(tot_av(ii)-totv)*(tot_av(ii)-totv)
        enddo
        
        std_v=sqrt(std_v/aver/(aver-1))
        std_s=sqrt(std_s/aver/(aver-1))
        std_t=sqrt(std_t/aver/(aver-1))
        
        if (mod(aver,2).eq.0) then
         median_v=(vol_av(aver/2)+vol_av(aver/2+1))/2
         median_s=(surf_av(aver/2)+surf_av(aver/2+1))/2
         median_t=(tot_av(aver/2)+tot_av(aver/2+1))/2
        else 
          median_v=vol_av((aver-1)/2+1)
          median_s=surf_av((aver-1)/2+1)
          median_t=tot_av((aver-1)/2+1)
        endif
 
!        write(86,*)'=========================================='
!        write(86,*)'|          Averaged values               |'
!        write(86,*)'=========================================='
        if (aver.gt.1) then
        write(86,'(a31,1x,f12.2,1x,a3,1x,f12.2,1x,a5)')
     $   'Outer surface of the structure: '
     $   ,surfv,'+/-',std_s,'[A^2]'
        else
        write(86,'(a31,1x,f12.2,1x,a5)')
     $   'Outer surface of the structure: '
     $   ,surfv,'[A^2]'
        endif
        write(86,'(a4,1x,f12.2,1x,a5)')'Min: ',surf_av(1),'[A^2]'
        write(86,'(a4,1x,f12.2,1x,a5)')'Max: ',surf_av(aver),'[A^2]'
        write(86,'(a7,1x,f12.2,1x,a5)')'Median: ',median_s,'[A^3]'        
        write(86,*)'------------------------------------------'
        if (aver.gt.1) then
        write(86,'(a50,1x,f12.2,1x,a3,1x,f12.2,1x,a5)')
     $   'Volume of the whole structure (with the cavities): '
     $   ,totv,'+/-',std_t,'[A^3]'
        else
        write(86,'(a50,1x,f12.2,1x,a5)')
     $   'Volume of the whole structure (with the cavities): '
     $   ,totv,'[A^3]'
        endif
        write(86,'(a4,1x,f12.2,1x,a5)')'Max: ',tot_av(1),'[A^3]'
        write(86,'(a4,1x,f12.2,1x,a5)')'Min: ',tot_av(aver),'[A^3]'
        write(86,'(a7,1x,f12.2,1x,a5)')'Median: ',median_t,'[A^3]'
        write(86,*)'------------------------------------------'
        if (aver.gt.1) then
        write(86,'(a29,1x,f12.2,1x,a3,1x,f12.2,1x,a5)')
     $   'Volume of the largest cavity: ',volv,'+/-',std_v,
     $   '[A^3]'
        else
        write(86,'(a29,1x,f12.2,1x,a5)')
     $   'Volume of the largest cavity: ',volv,
     $   '[A^3]'
        endif
        write(86,'(a4,1x,f12.2,1x,a5)')'Min: ',vol_av(1),'[A^3]'
        write(86,'(a4,1x,f12.2,1x,a5)')'Max: ',vol_av(aver),'[A^3]'
        write(86,'(a7,1x,f12.2,1x,a5)')'Median: ',median_v,'[A^3]'
        write(86,*)'------------------------------------------'
        write(86,*)'Results obtained for the individual orientations'
        write(86,*)'=========================================='
        write(86,*)'|         Volume of the cavity           |'
        write(86,*)'=========================================='        
        do ii=1,aver
        write(86,'(9x,i3,1x,f12.2,1x,a5)')ii,vol_av(ii),'[A^3]'
        enddo
        write(86,*)'=========================================='
        write(86,*)'|      Surface of the structure          |'
        write(86,*)'=========================================='
        do ii=1,aver
        write(86,'(9x,i3,1x,f12.2,1x,a5)')ii,surf_av(ii),'[A^2]'
        enddo
        write(86,*)'=========================================='
        write(86,*)'|       Total volume of the structure    |'
        write(86,*)'=========================================='
        do ii=1,aver
        write(86,'(9x,i3,1x,f12.2,1x,a5)')ii,tot_av(ii),'[A^3]'
        enddo
      
       
      end     
