      subroutine assign_vdw_radii (pdb_name, nres, nat, 
     $                             seq, na, at, vrad )
C THIS SUBROUTINE ASSIGNS VAN DER WAALS RADII
      character (len=*) pdb_name
      integer*4 nres		! MAX NUMBER OF RESIDUES
      integer*4 nat		! MAX NUMBER OF ATOMS IN A RESIDUE
      character*3 seq(nres)	! AMINO ACID SEQUENCE
      integer*4 na(nres)		! NUMBER OF ATOMS IN RESIDUES
      character*2 at(nres, nat)	! HEAVY ATOMS N, S, O, C
      real*4 vrad(nres, nat)	! VAN DER WAALS RADII
      integer*4 ib, iat, j
      integer*4 nb(nres)
      character*3 ares
      character*1 ana
      character*80 outputfile
      real*8 rad
      do ib = 1, nres
       do j=1,na(ib)
        if (at(ib,j).eq.'H ') vrad(ib,j)=1.00d0
        if (at(ib,j).eq.'N ') vrad(ib,j)=1.50d0
        if (at(ib,j).eq.'P ') vrad(ib,j)=1.90d0
        if (at(ib,j).eq.'AS') vrad(ib,j)=2.00d0
        if (at(ib,j).eq.'SB') vrad(ib,j)=2.20d0
        if (at(ib,j).eq.'O ') vrad(ib,j)=1.40d0
        if (at(ib,j).eq.'S ') vrad(ib,j)=1.85d0
        if (at(ib,j).eq.'SE') vrad(ib,j)=2.00d0
        if (at(ib,j).eq.'TE') vrad(ib,j)=2.20d0
        if (at(ib,j).eq.'F ') vrad(ib,j)=1.35d0
        if (at(ib,j).eq.'CL') vrad(ib,j)=1.80d0
        if (at(ib,j).eq.'BR') vrad(ib,j)=1.95d0
        if (at(ib,j).eq.'I ') vrad(ib,j)=2.15d0
        if (at(ib,j).eq.'C ') vrad(ib,j)=1.70d0
       enddo
      enddo
      !stop 
      return
      end

