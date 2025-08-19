      subroutine get_nres (pdbfile, nres, nat)
C THIS SUBROUTINE READS THE PDB FILE AND FINDS 
C THE NUMBER OF RESIDUES (nres) AS WELL AS THE MAXIMAL 
C NUMBER OF HEAVY ATOMS IN A RESIDUE (nat) WHICH IS USUALLY 14
      implicit none
      character (len=*) pdbfile 	! PDB FILE
      integer*4 nres		! NUMBER OR RESIDUES
      integer*4 nat
      integer*4 iat
      integer*4 ss
      character*1 a,b,d
      character*80 line
      integer test1

      
      ss=13 ! SPECIFICATION OF A LOGICAL UNIT FOR INPUT FILE
      test1=0  !to zmienilem
      
      
      open (ss, file=trim(pdbfile),form='formatted',status='old')
      
      nres=1
      nat=0
      iat=0
 1    read (ss, '(a)', end = 100, err = 100 ) line
      
      if ((line(1:4).eq.'ATOM')) then
       !if ( line (14:16) .eq. 'N  ' ) then
        ! nres = nres + 1
        ! if ( iat .gt. Nat ) then
        !  nat = iat
        ! endif
        ! iat = 0
       !endif
       test1=1 !to zmienilem
       a=line(14:14)
       b=line(15:15)
       d=line(13:13)
       !b=line(33:39)
       if ((a.eq.'H').or.
     $      (a.eq.'N').or.
     $      (a.eq.'P').or.
     $      ((a.eq.'A').and.(b.eq.'S')).or.
     $      ((a.eq.'S').and.(b.eq.'B')).or.
     $      (a.eq.'O').or.
     $      (a.eq.'S').or.
     $      ((a.eq.'S').and.(b.eq.'E')).or.
     $      ((a.eq.'T').and.(b.eq.'E')).or.
     $      (a.eq.'F').or.
     $      ((a.eq.'C').and.(b.eq.'L')).or.
     $      ((a.eq.'B').and.(b.eq.'R')).or.
     $     (a.eq.'I').or.
     $      (a.eq.'C').or.
     $      (d.eq.'H')) then
        iat = iat + 1
           if (mod(iat,10).eq.0) then
             nres=nres+1
             iat=0
             test1=0 !to zmienilem
        endif
       endif
      endif
       
      
!      if (line(1:4).eq.'ATOM') then
!       if ( line (14:16) .eq. 'P  ' ) then
!         nres = nres + 1
!         if ( iat .gt. Nat ) then
!          nat = iat
!         endif
!         iat = 0
!       endif
!       a=line(14:14)
!       if ((a.eq.'N').or.(a.eq.'S').or.(a.eq.'O').or.(a.eq.'C')
!     $  .or.(a.eq.'P').or.(a.eq.'F').or.(a.eq.'Cl')) then
!        iat = iat + 1
!       endif
!      endif
            
      goto 1

 100  close (ss)
 
      nat=10
      if (test1.eq.0) nres=nres-1 ! to zmieniłem
      
      !write(*,*)nat,nres
!      if ( nat .lt. 14 ) then
!        nat = 14
!      endif
      return
      end
