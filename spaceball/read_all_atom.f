      subroutine read_all_atom (pdbfile, ires, iat, seq, na, 
     $  x, y, z, at,tot_at)
C THIS SURROUTINE READS IN THE STRUCTURE FROM THE PDB FILE
C IT IS USED IN THE CONTACT MAP CALCULATION
      implicit none
      character (len=*) pdbfile ! PDB FILE
      integer*4 nres		! NUMBER OF RESIDUES
      integer*4 nat,tot_at		! MAX NUMBER OF ATOMS IN A RESIDUE
      character*3 seq(ires)	! AMINO ACID SEQUENCE
      integer*4 na(ires)		! NUMBER OF ATOMS IN RESIDUES
      character*2 at(ires,iat)	! HEAVY ATOMS N, S, O, C
      real*4 x(ires, iat)	! X-COORDINATES OF ATOMS
      real*4 y(ires, iat)	! Y-COORDINATES OF ATOMS
      real*4 z(ires, iat)	! Z-COORDINATES OF ATOMS
      integer*4 iresn(ires)	! RESIDUR NUMBERS IN PDB FILE
      logical struc_error	! PROBLEMS WITH PDB FILE
      integer*4 ires, iat
      integer*4 i,j
      character*80 line, outputfile
      character*1 a,b,c,d
      integer test1
      
      
      tot_at=0
      iat=0
      test1=0  !to zmienilem
      
C OPEN PDB FILE
      open(20, file=trim(pdbfile),form='formatted',status='old')
       
      ires=1
       
 1    read(20, '(a)', end = 100, err = 100) line
      
      if ((line(1:4).eq.'ATOM')) then
        !if (line (14:16) .eq. 'N  ' ) then
            
         !   ires = ires + 1
            test1=1 !to zmienilem
            read( line(18:20), '(a3)' ) seq(ires)
            read( line(24:26), '(i3)' ) iresn(ires)
         !   iat = 0
         !endif
         read( line(14:14), '(a1)' ) a
         read( line(15:15), '(a1)' ) b
         read( line(13:13), '(a1)' ) d
     
       if ((a.eq.'H').or.
     $      (a.eq.'N').or.  
     $      (a.eq.'P').or.
     $     ((a.eq.'A').and.(b.eq.'S')).or.
     $     ((a.eq.'S').and.(b.eq.'B')).or.
     $      (a.eq.'O').or.
     $      (a.eq.'S').or.
     $     ((a.eq.'S').and.(b.eq.'E')).or.
     $     ((a.eq.'T').and.(b.eq.'E')).or.
     $      (a.eq.'F').or.
     $     ((a.eq.'C').and.(b.eq.'L')).or.
     $     ((a.eq.'B').and.(b.eq.'R')).or.
     $     (a.eq.'I').or.
     $      (a.eq.'C').or.
     $      (d.eq.'H')) then
                                                                              
        iat=iat+1
        tot_at=tot_at+1
        na(ires)=iat
        c=' '
        
        if ((a.eq.'A').and.(b.eq.'S')) c='S'
        if ((a.eq.'S').and.(b.eq.'B')) c='B'
        if ((a.eq.'S').and.(b.eq.'E')) c='E'
        if ((a.eq.'T').and.(b.eq.'E')) c='E'
        if ((a.eq.'C').and.(b.eq.'L')) c='L'
        if ((a.eq.'B').and.(b.eq.'R')) c='R'
        if (d.eq.'H') then 
         c=' '
         a='H'
        endif
        
!        if ((((a.eq.'A').and.(b.ne.'S')).or.
!     $       ((a.eq.'S').and.(b.ne.'B')).or.
!     $       ((a.eq.'S').and.(b.ne.'E')).or.
!     $       ((a.eq.'T').and.(b.ne.'E')).or.
!     $       ((a.eq.'C').and.(b.ne.'L')).or.
!     $      ((a.eq.'B').and.(b.ne.'R'))).or.
!     $       ((b.ne.'S').and.(b.ne.'B').and.
!     $       (b.ne.'E').and.(b.ne.'L').and.
!     $       (b.ne.'R'))) b=' ' 
        at(ires, iat) =a//c
         
         read( line(31:54), '(3f8.3)' ) x(ires, iat), 
     $   y(ires, iat), z(ires, iat)
        if (mod(iat,10).eq.0) then
         iat=0
         ires=ires+1
         test1=0 !to zmienilem  
        endif
       endif
      endif      
      
      goto 1
 100  close (20)
       
       if (test1.eq.0) ires=ires-1 ! to zmieniłem
       
       iat=10                    
      return
      end 
