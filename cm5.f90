!   The MIT License (MIT)
!  
!   Copyright (c) 2016 Holger Kruse
!  
!   Permission is hereby granted, free of charge, to any person obtaining a copy
!   of this software and associated documentation files (the "Software"), to deal
!   in the Software without restriction, including without limitation the rights
!   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
!   copies of the Software, and to permit persons to whom the Software is
!   furnished to do so, subject to the following conditions:
!  
!   The above copyright notice and this permission notice shall be included in all
!   copies or substantial portions of the Software.
!  
!   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
!   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
!   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
!   SOFTWARE.
!  
program cm5charge
implicit none
integer nat,i,j,k,l,ii,jj,kdelta
integer, allocatable :: iat(:)
real(8), allocatable :: xyz(:,:),qcm5(:),qh(:)
real(8) dx,dy,dz,rij,s,alpha,tkk,bkk,fscale,dipole(3)
character(120) infile,chrgfile
character(2) esym
real(kind=8) rvdw(94),rconv(94),bohr,quad(3,3)
real(8), parameter :: AMBER_ELECTROSTATIC = 18.2223d0
real(8), parameter :: INV_AMBER_ELECTROSTATIC = 1.0d0/AMBER_ELECTROSTATIC


! atomic radii from Mantina, Valero, Cramer, Truhlar "Atomic radii of elements"
! Copyed by hand,  maybe contain typos...
!            H       He
data rvdw /1.10d0,1.40d0, &
    ! Li     Be     B     C       N      O     F      Ne
    1.82d0,1.53d0,1.92d0,1.70d0,1.54d0,1.52d0,1.47d0,1.54d0, &
    ! Na    Mg     Al     Si     P      S       Cl     Ar
    2.27d0,1.73d0,1.84d0,2.10d0,1.80d0,1.80d0,1.75d0,1.88d0, &
    ! K      Ca     Sc     Ti     V      Cr      Mn     Fe     Co    Ni     Cu     Zn     Ga     Ge     As     Se     Br    Kr
    2.75d0,2.31d0,2.15d0,2.11d0,2.07d0,2.06d0,2.05d0,2.04d0,2.00d0,1.97d0,1.96d0,2.01d0,1.87d0,2.11d0,1.85d0,1.90d0,1.85d0,2.02d0, &
    !  Rb    Sr     Y      Zr      Nb     Mo    Tc     Ru     Rh     Pd     Ag     Cd     In    Sn      Sb      Te     I     Xe
    3.03d0,2.49d0,2.26d0,2.23d0,2.18d0,2.17d0,2.16d0,2.13d0,2.10d0,2.10d0,2.11d0,2.18d0,1.93d0,2.17d0,2.06d0,2.06d0,1.98d0,2.16d0, &
    ! Cs Ba
    3.32d0,2.68d0, &
    ! La-Lu
    2.43d0,2.42d0,2.40d0,2.46d0,2.38d0,2.36d0,2.35d0,2.34d0,2.33d0,2.31d0,2.30d0,2.29d0,2.27d0,2.26d0,2.24d0, &
    ! Hf     Ta     W       Re     Os    Ir     Pt     Au     Hg     Ti     Pb     Bi     Po     At     Rn
    2.23d0,2.22d0,2.18d0,2.16d0,2.16d0,2.13d0,2.13d0,2.23d0,2.23d0,2.11d0,2.02d0,2.07d0,1.97d0,2.02d0,2.20d0, &
    ! Fr-Pu
    3.48d0,2.83d0,2.47d0,2.45d0,2.43d0,2.41d0,2.39d0,2.43d0/


data rconv /0.32d0,0.37d0, &
    ! Li     Be     B     C       N      O     F      Ne
    1.30d0,0.99d0,0.84d0,0.75d0,0.71d0,0.64d0,0.60d0,0.62d0,  &
    ! Na    Mg     Al     Si     P      S       Cl     Ar
    1.60d0,1.40d0,1.24d0,1.14d0,1.09d0,1.04d0,1.00d0,1.01d0,  &
    ! K      Ca     Sc     Ti     V      Cr      Mn     Fe     Co    Ni     Cu     Zn     Ga     Ge     As     Se     Br    Kr
    2.00d0,1.74d0,1.59d0,1.48d0,1.44d0,1.30d0,1.29d0,1.24d0,1.18d0,1.17d0,1.22d0,1.20d0,1.23d0,1.20d0,1.20d0,1.18d0,1.17d0,1.24d0, &
    !  Rb    Sr     Y      Zr      Nb     Mo    Tc     Ru     Rh     Pd     Ag     Cd     In    Sn      Sb      Te     I     Xe
    2.15d0,1.90d0,1.78d0,1.64d0,1.56d0,1.46d0,1.38d0,1.36d0,1.34d0,1.30d0,1.36d0,1.40d0,1.42d0,1.40d0,1.40d0,1.37d0,1.32d0,1.36d0, &
    ! Cs Ba
    2.38d0,2.06d0,  &
    ! La-Lu
     1.94d0,1.84d0,1.90d0,1.73d0,1.86d0,1.85d0,1.83d0,1.82d0,1.81d0,1.80d0,1.79d0,1.77d0,1.77d0,1.78d0,1.74d0,  &
    ! Hf     Ta     W       Re     Os    Ir     Pt     Au     Hg     Ti     Pb     Bi     Po     At     Rn
    1.64d0,1.58d0,1.50d0,1.41d0,1.36d0,1.32d0,1.30d0,1.64d0,1.88d0,1.48d0,1.45d0,1.50d0,1.42d0,1.47d0,1.46d0,  &
    ! Fr-Pu
    2.42d0,2.11d0,2.01d0,1.90d0,1.84d0,1.83d0,1.80d0,1.80d0/



print*,'|-----------------------------------------|'
print*,'| CM5 charges calculator from             |'
print*,'| ORCA Hirsheld charges                   |'
print*,'| set: "%output print[P_Hirshfeld] 1 end" |'
print*,'| v1.0 Holger Kruse                       |'
print*,'|-----------------------------------------|'
print*,''

bohr=0.52917726d0

alpha=2.474  ! angstrom

call getarg(1,infile)
call getarg(2,chrgfile)

if(chrgfile=='') then
  fscale=1.20
else
 read(chrgfile,*) fscale
endif

if(infile=='') then
print*,'USAGE: cm5 <orca output> [<scale factor]'
print*,''
stop 
endif

write(*,*) ''
write(*,'(a,F4.2)'),'* scaling factor for CM5 *: ', fscale
write(*,*) ''

nat=1
call tmolrd(nat,xyz,iat,infile,.false.,.true.)
allocate(xyz(3,nat),iat(nat))
allocate(qh(nat),qcm5(nat))
call tmolrd(nat,xyz,iat,infile,.true.,.false.)

xyz=xyz*bohr ! angstrom

write(*,'(a)'),' data in electron/angstrom '
call read_hirsh(qh,nat,trim(infile))
!print*,' Hirshfeld charges '
!do i=1,nat
!write(*,'(I3,a,F5.2)') i,esym(iat(i)),qh(i)
!enddo

qcm5=0d0
do i=1,nat
ii=iat(i)
 s=0
 do j=1,nat
  if(i==j) cycle
   jj=iat(j)
   dx=xyz(1,i)-xyz(1,j)
   dy=xyz(2,i)-xyz(2,j)
   dz=xyz(3,i)-xyz(3,j)
   rij=sqrt(dx*dx+dy*dy+dz*dz)
   bkk=exp(-alpha*(rij-rconv(ii)-rconv(jj))) ! eq.2
   s=s+tkk(ii,jj)*bkk
!  print*,s,tkk(ii,jj)*bkk,tkk(ii,jj),bkk,rij
 enddo
 qcm5(i)=qh(i)+s
enddo

print*,'CHARGES:'
print*,'           CM5      CM5scale     HPA         X            Y            Z'
do i=1,nat
write(*,'(I3,2x,a2,2x,3(F9.5,x),3(F12.7,x))') i,esym(iat(i)),qcm5(i),qcm5(i)*fscale,qh(i),xyz(1:3,i)
enddo

write(*,*) ''
write(*,*) 'electrostatic moments'
write(*,*) ' ref. point: 0.0  0.0  0.0'

!CM5 dipole moment
dipole=0
do i=1,nat
 do j=1,3
dipole(j)=dipole(j)+qcm5(i)*xyz(j,i)
enddo
enddo

dipole=dipole*4.802889778d0
dx= dipole(1)
dy= dipole(2)
dz= dipole(3)
s=sqrt(dipole(1)**2+dipole(2)**2+dipole(3)**2)

write(*,'(a,F10.6)'     ) '  Total CM5 charge       : ',sum(qcm5)
write(*,*)                '                             X          Y          Z      TOTAL'
write(*,'(a,4(F10.6,x))') '  CM5 dipole [Debye]  : ', dx,dy,dz,s


! scaled 
! could have just scaled the dipole directly, oh well...
dipole=0
do i=1,nat
 do j=1,3
dipole(j)=dipole(j)+qcm5(i)*xyz(j,i)*fscale
enddo
enddo



dipole=dipole*4.802889778d0
dx= dipole(1)
dy= dipole(2)
dz= dipole(3)
s=sqrt(dipole(1)**2+dipole(2)**2+dipole(3)**2)

write(*,'(a,4(F10.6,x))') '  CM5 scaled [Debye]  : ', dx,dy,dz,s


! Quadrupole moment
quad=0
  do k=1,nat
   dx=xyz(1,k)/bohr
   dy=xyz(2,k)/bohr
   dz=xyz(3,k)/bohr
   quad(1,1)=quad(1,1)+ dx*dx*qcm5(k)
   quad(2,2)=quad(2,2)+ dy*dy*qcm5(k)
   quad(3,3)=quad(3,3)+ dz*dz*qcm5(k)
   quad(1,2)=quad(1,2)+ dx*dy*qcm5(k)
   quad(1,3)=quad(1,3)+ dx*dz*qcm5(k)
   quad(2,3)=quad(2,3)+ dy*dz*qcm5(k)
  enddo


!quad=quad*4.802889778d0
print*, ' quadrupole moment (au):'
write(*,'(2x,3(a,F12.6,x))') 'XX= ',quad(1,1),'YY= ',quad(2,2),'ZZ= ',quad(3,3)
write(*,'(2x,3(a,F12.6,x))') 'XY= ',quad(1,2),'XZ= ',quad(1,3),'YZ= ',quad(2,3)
write(*,'(2x,a,F12.6)') '  1/3 trace= ',(quad(1,1)+quad(2,2)+quad(3,3)/3d0)

! write data
print*,'writing: '//trim(infile)//'_cm5.dat with normal and  scaled CM5 charges'
open(99,file=trim(infile)//'_cm5.dat')
do i=1,nat
write(99,'(2(F10.6,x))') qcm5(i),qcm5(i)*fscale
enddo
close(99)


! write topology file section
qcm5=qcm5*fscale/INV_AMBER_ELECTROSTATIC
print*,'writing: '//trim(infile)//'_cm5.top (amber topology) with scaled CM5 charges'
open(111,file=trim(infile)//'_cm5.top')
  write(111,'(a)')'%FLAG CHARGE'
  write(111,'(a)')'%FORMAT(5E16.8)'
  write(111,'(5E16.8)') (qcm5(i),i = 1,nat)
close(111)



end program

! form Tkk' eq.3,4
real(kind=8) function tkk(i,j)
implicit none
integer i,j
real(8) dz(54)
real(8) dzz(6)

data dz /0.0056,-0.1543, &
 0.0, 0.0333, -0.1030, -0.0446, -0.1072, -0.0802, -0.0629, -0.1088, & ! Li-Ne
0.0184, 0.0 ,-0.0726, -0.0790, -0.0756, -0.0565, -0.0444, -0.07676,  & !Na-Ar
0.0130, 0.0,10*0.0,-0.0557,-0.0533,-0.0399, -0.0313,0.0,0.0, & !K-Kr
16*0.0, -0.0220,0.0/ !Rb - Xe

!H-C; H-N H-O C-N C-O N-O
data dzz /0.0502d0,0.1747d0,0.1671d0,0.0556d0, 0.0234d0,-0.0346d0/

tkk=dz(i)-dz(j)

! catch special cases
 if(i==1.and.j==6) tkk=dzz(1)
 if(i==1.and.j==7) tkk=dzz(2)
 if(i==1.and.j==8) tkk=dzz(3)
 if(i==6.and.j==7) tkk=dzz(4)
 if(i==6.and.j==8) tkk=dzz(5)
 if(i==7.and.j==8) tkk=dzz(6)

 if(i==6.and.j==1) tkk=-dzz(1)
 if(i==7.and.j==1) tkk=-dzz(2)
 if(i==8.and.j==1) tkk=-dzz(3)
 if(i==7.and.j==6) tkk=-dzz(4)
 if(i==8.and.j==6) tkk=-dzz(5)
 if(i==8.and.j==7) tkk=-dzz(6)


end function

subroutine read_hirsh(qh,nat,fname)
implicit none
integer nat,i
real(8) qh(nat)
character(*) fname
character(80) aa
real(8) xx,yy

write(*,*) 'reading ORCA file :  ', trim(fname)
open(99,file=fname)
do
read(99,'(a)') aa
if(index(aa,'HIRSHFELD ANALYSIS').ne.0) then
! 6 dummy line, then read
read(99,'(a)') aa
read(99,'(a)') aa
read(99,'(a)') aa
read(99,'(a)') aa
read(99,'(a)') aa
read(99,'(a)') aa
do i=1,nat
read(99,*) yy,aa,qh(i),xx
enddo
exit
endif
enddo

! read charges file for testing
!open(99,file='charges')
!do i=1,nat
!read(99,*) qh(i)
!enddo




close(99)



end



!************************************************************
!* reads a turbomole (bohr) or xmol (angst)rom file.        *
!* Tests if xmol starts with "number of atoms + blank" or   *
!* directly with the coordinates.                           *
!************************************************************
subroutine tmolrd(nat,xyz,iat,infile,echo,c_nat)
!use parm
implicit none

integer i,j,k,l,nat
real(8) xyz(3,*),energy,s2r
integer ifrez(nat),iat(nat)

character*2 cc,ff
character(200)  atmp
character*(*) infile
real(kind=8) txyz(3,nat),xx(10)
real(kind=8) bohr
integer tiat(nat),nn,tifrez(nat),istat,iff
logical da,c_nat,echo,doORCA
bohr=0.52917726d0
i=0
tifrez=0
iff=0
doORCA=.false.

inquire(file=infile,exist=da)
select case (da)
case (.true.)
      if(echo) write(*,'('' reading...'',$)')

open(unit=33,file=infile)
! test for tmol or txyz file

! test for ORCA
do i=1,10
 read(33,'(a)') atmp ! $coord
if(index(atmp,' O   R   C   A ').ne.0) doORCA=.true.
enddo
rewind(33)

 read(33,'(a)') atmp ! $coord
rewind(33)
if(index(atmp,'$coord').ne.0) then

 ! count number of atoms
 do while (da)
  read(33,'(a)',end=100) atmp ! $coord
   if(index(atmp,'$coord').ne.0) cycle
   if(index(atmp,'$').ne.0) exit
   i=i+1
  enddo
 nat=i
 100 continue
 if(c_nat) then
  close(33)
  return  ! just return number of atoms
 endif
 rewind(unit=33)

 ! read TMOL file
 read(33,*) atmp ! $coord
 do j=1,nat
    read(33,'(a)') atmp ! $coord
    backspace(33)
    if(index(atmp,' f ').ne.0) then
    read(33,*) txyz(1,j),txyz(2,j),txyz(3,j),cc,ff
     tifrez(j)=1
    iff=iff+1
    else ! default
     read(33,*) txyz(1,j),txyz(2,j),txyz(3,j),cc
   endif
   call elem(cc,tiat(j))
   txyz(1:3,j)=txyz(1:3,j)
  enddo
 if(echo) write(*,*) ' Turbomole file [bohr] :  ', trim(infile)

 close(33)

!**********
! PDB FILE
!**********
elseif(index(atmp,'REMARK').ne.0.or.index(atmp,'ATOM').ne.0) then
rewind(33)

 if(c_nat) then ! just count atoms
  do
   read (33,'(a)',end=777) atmp
   if(index(atmp,'ATOM').ne.0) i=i+1
  enddo
  777 continue
  nat=i
  close(33)
  return
  endif

   i=0
   do
     read (33,'(a)',end=778) atmp
     if(index(atmp,'ATOM').ne.0) then
      i=i+1
      txyz(1,i)=s2r(atmp(31:38))
      txyz(2,i)=s2r(atmp(39:46))
      txyz(3,i)=s2r(atmp(47:54))
     endif
   enddo
  txyz=0
 778 continue
   txyz=txyz*1d0/bohr
   close(33)
   if(echo) write(*,'(5x,'' PDB file [angst]: '',a)')  trim(infile)
!******************
!* ORCA output file
!******************
elseif(doORCA) then
i=0
if(c_nat) then ! just count atoms
do
 read(33,'(a)',end=801) atmp
 if(index(atmp,'CARTESIAN COORDINATES (A.U.)').ne.0) then
  read(33,'(a)') atmp
 do
  read(33,'(a)') atmp
  i=i+1
  if(index(atmp,'------').ne.0) goto 800
 enddo 
 endif
enddo
800 continue
nat=i-3
close(33)
801 continue
if(nat<=0) stop 'no atoms found!'
return
endif

i=0
do
 read(33,'(a)',end=779) atmp
 if(index(atmp,'CARTESIAN COORDINATES (A.U.)').ne.0) then
 read(33,'(a)') atmp ! dummy line
 read(33,'(a)') atmp ! dummy line
 do i=1,nat
   read(33,'(a)') atmp ! dummy line
   call readl(atmp,xx,nn)
   call elem(atmp,tiat(i))
   txyz(1:3,i)=xx(5:7)
 enddo
 endif
enddo
779 continue
close(33)
 if(echo) write(*,'(5x,'' ORCA output file [A.U.]: '',a)')  trim(infile)

!******************
!* TMOL output file
!******************
else ! txyz file
       read(33,'(a)',end=101) atmp
! check for first two lines
       call readl(atmp,xx,nn)
        if(nn.gt.1) then   ! more than one argument found, assuming they are coords
           do
            nat=nat+1
            read(33,'(a)',end=123) atmp
           enddo
            if(c_nat) then
              close(33)
              return  ! just return number of atoms
            endif
          else
            nat=idint(xx(1))
            if(c_nat) then
             close(33)
             return  ! just return number of atoms
            endif
           read(33,'(a)',end=101) atmp  !titel line
        endif
 123   if(nn.gt.1) rewind(3)
       do i=1,nat
            read(33,'(a)') atmp
            call readl(atmp,xx,nn)
            call elem(atmp,tiat(i))
            txyz(1:3,i)=xx(1:3)*1d0/bohr
       enddo
 101  close(33)
      if(echo) write(*,'(5x,'' XYZ file [angst]: '',a)')  trim(infile)
      endif

if(maxval(tifrez,nat).eq.1) then
 if(echo) then
  write(*,'(a,x,I4,x,a)') '  found ',iff,' frozen cart. coordinates'
  if(iff.lt.50) then ! dont spam the output to much ...
   write(*,'(a,$)') '  atom nr: '
   do i=1,nat
     if(tifrez(i)==1) write(*,'(x,I2,$)') i
   enddo
   print*,''
  endif
 endif
endif

case (.false.)
  write(*,*) ' no input file <',trim(infile) ,'> found !! '
end select

do i=1,nat
xyz(1:3,i)=txyz(1:3,i)
iat(i)=tiat(i)
ifrez(i)=tifrez(i)
enddo
return

end subroutine


 SUBROUTINE ELEM(KEY1, NAT)
 IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 CHARACTER*(*) KEY1
 CHARACTER*2 ELEMNT(94),E

 DATA ELEMNT/'h ','he',                                      &
 'li','be','b ','c ','n ','o ','f ','ne',                    &
 'na','mg','al','si','p ','s ','cl','ar',                    &
 'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu',     &
 'zn','ga','ge','as','se','br','kr',                         &
 'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag',     &
 'cd','in','sn','sb','te','i ','xe',                         &
 'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy',&
 'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt',&
 'au','hg','tl','pb','bi','po','at','rn',                    &
 'fr','ra','ac','th','pa','u ','np','pu'/

 nat=0
 e='  '
 k=1
 DO J=1,len(key1)
    if (k.gt.2)exit
    N=ICHAR(key1(J:J))
    if(n.ge.ichar('A') .and. n.le.ichar('Z') )then
       e(k:k)=char(n+ICHAR('a')-ICHAR('A'))
       k=k+1
    endif
    if(n.ge.ichar('a') .and. n.le.ichar('z') )then
       e(k:k)=key1(j:j)
       k=k+1
    endif
 enddo
 DO I=1,107
         if(e.eq.elemnt(i))then
            NAT=I
            RETURN
         ENDIF
      ENDDO

      end

character*2 FUNCTION ESYM(I)
CHARACTER*2 ELEMNT(94)
DATA ELEMNT/'h ','he',                                           &
  'li','be','b ','c ','n ','o ','f ','ne',                       &
  'na','mg','al','si','p ','s ','cl','ar',                       &
  'k ','ca','sc','ti','v ','cr','mn','fe','co','ni','cu',        &
  'zn','ga','ge','as','se','br','kr',                            &
  'rb','sr','y ','zr','nb','mo','tc','ru','rh','pd','ag',        &
  'cd','in','sn','sb','te','i ','xe',                            &
  'cs','ba','la','ce','pr','nd','pm','sm','eu','gd','tb','dy',   &
  'ho','er','tm','yb','lu','hf','ta','w ','re','os','ir','pt',   &
  'au','hg','tl','pb','bi','po','at','rn',                       &
  'fr','ra','ac','th','pa','u ','np','pu'/
  ESYM=ELEMNT(I)
  RETURN
END

      SUBROUTINE READL(A1,X,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*(*) A1
      DIMENSION X(*)
      I=0
      IS=1
  10  I=I+1
      X(I)=READAA(A1,IS,IB,IE)
      IF(IB.GT.0 .AND. IE.GT.0) THEN
                                IS=IE
                                GOTO 10
      ENDIF
      N=I-1
      RETURN
      END


      FUNCTION READAA(A,ISTART,IEND,IEND2)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 READAA
      CHARACTER*(*) A
      NINE=ICHAR('9')
      IZERO=ICHAR('0')
      MINUS=ICHAR('-')
      IDOT=ICHAR('.')
      ND=ICHAR('D')
      NE=ICHAR('E')
      IBL=ICHAR(' ')
      IEND=0
      IEND2=0
      IDIG=0
      C1=0
      C2=0
      ONE=1.D0
      X = 1.D0
      NL=LEN(A)
      DO 10 J=ISTART,NL-1
         N=ICHAR(A(J:J))
         M=ICHAR(A(J+1:J+1))
         IF(N.LE.NINE.AND.N.GE.IZERO .OR.N.EQ.IDOT)GOTO 20
         IF(N.EQ.MINUS.AND.(M.LE.NINE.AND.M.GE.IZERO.OR. M.EQ.IDOT)) GOTO 20

   10 CONTINUE
      READAA=0.D0
      RETURN
   20 CONTINUE
      IEND=J
      DO 30 I=J,NL
         N=ICHAR(A(I:I))
         IF(N.LE.NINE.AND.N.GE.IZERO) THEN
            IDIG=IDIG+1
            IF (IDIG.GT.10) GOTO 60
            C1=C1*10+N-IZERO
         ELSEIF(N.EQ.MINUS.AND.I.EQ.J) THEN
            ONE=-1.D0
         ELSEIF(N.EQ.IDOT) THEN
            GOTO 40
         ELSE
            GOTO 60
         ENDIF
   30 CONTINUE
   40 CONTINUE
      IDIG=0
      DO 50 II=I+1,NL
         N=ICHAR(A(II:II))
         IF(N.LE.NINE.AND.N.GE.IZERO) THEN
            IDIG=IDIG+1
            IF (IDIG.GT.10) GOTO 60
            C2=C2*10+N-IZERO
            X = X /10
         ELSEIF(N.EQ.MINUS.AND.II.EQ.I) THEN
            X=-X
         ELSE
            GOTO 60
         ENDIF
   50 CONTINUE
!C
!C PUT THE PIECES TOGETHER
!C
   60 CONTINUE
      READAA= ONE * ( C1 + C2 * X)
      DO 55 J=IEND,NL
         N=ICHAR(A(J:J))
         IEND2=J
         IF(N.EQ.IBL)RETURN
   55 IF(N.EQ.ND .OR. N.EQ.NE)GOTO 57
      RETURN

   57 C1=0.0D0
      ONE=1.0D0
      DO 31 I=J+1,NL
         N=ICHAR(A(I:I))
         IEND2=I
         IF(N.EQ.IBL)GOTO 70
         IF(N.LE.NINE.AND.N.GE.IZERO) C1=C1*10.0D0+N-IZERO
         IF(N.EQ.MINUS)ONE=-1.0D0
   31 CONTINUE
   61 CONTINUE
   70 READAA=READAA*10**(ONE*C1)
      RETURN
      END
!********************************
!* convert a word to lower case *
!********************************
      subroutine lower_case(word)
      character (len=*) , intent(in out) :: word
      integer :: i,ic,nlen
      nlen = len(word)
      do i=1,nlen
      ic = ichar(word(i:i))
      if (ic >= 65 .and. ic < 90) word(i:i) = char(ic+32)
      end do
      end subroutine lower_case


!********************************************
!* split a string s into n separate words w *
!********************************************
subroutine charsplit(s,n,w)
implicit none
integer i,n,k
character*80, intent(in) :: s
character*80, intent(out) :: w(n)
character*80   a,aa

aa=adjustl(s)
do i=1,n
  a=aa
  k=index(a,' ')
  w(i)=trim(a(:(k-1)))
  aa=adjustl(trim(a((k+1):)))
enddo
return
end subroutine



!***************************************************
!* split a string s iand the return the x'ths word *
!***************************************************
subroutine charXsplit(s,wx,x)
implicit none
integer i,n,k,x
character*80, intent(in) :: s
character*80, intent(out) ::wx
character*80   w(20)
character*80   a,aa

aa=adjustl(s)
i=0
do while (index(aa,' ').ne.1)
  i=i+1
  a=aa
  k=index(a,' ')
  w(i)=trim(a(:(k-1)))
  aa=adjustl(trim(a((k+1):)))
!  print*,'AA',trim(aa),index(aa,' ')
  if(i.gt.50) stop '!* string split error: subroutine charXsplit *!'
enddo
wx=w(x)
!print*, trim(wx)
return
end subroutine


!***********************************
!*  count words n in string s      *
!***********************************
subroutine cstring(s,n)
implicit none
integer i,n,k,x
character(*), intent(in) :: s
character*80   a,aa
aa=adjustl(s)
i=0
do while (index(aa,' ').ne.1)
  i=i+1
  a=aa
  k=index(a,' ')
!  w(i)=trim(a(:(k-1)))
  aa=adjustl(trim(a((k+1):)))
!  print*,'AA',trim(aa),index(aa,' ')
  if(i.gt.50) stop '!* string split error: subroutine charXsplit *!'
enddo
n=i
return
end subroutine

integer function nwords(s)
implicit none
integer i,k,x
character(*), intent(in) :: s
character*80   a,aa
aa=adjustl(s)
i=0
do while (index(aa,' ').ne.1)
  i=i+1
  a=aa
  k=index(a,' ')
!  w(i)=trim(a(:(k-1)))
  aa=adjustl(trim(a((k+1):)))
!  print*,'AA',trim(aa),index(aa,' ')
  if(i.gt.50) stop '!* string split error: subroutine charXsplit *!'
enddo
nwords=i
return
end function



! String to integer
integer pure function s2i(a)
implicit none
character(*), intent(in):: a
read(a,*) s2i
return
end function


! String to real(8)
real(8) pure function s2r(a)
implicit none
character(*), intent(in):: a
read(a,*) s2r
return
end function




!******************
!* write xyz      *
!******************
subroutine wrxyz(iat,nat,xyz,infile)
implicit none
integer i,j,k,l,nat,iat(nat)
real(8) xyz(3,*)
integer ierr
character(*) infile
character(2) esym
real(8) f

!f=0.5291770d0
f=1d0
open(unit=55,file=infile,iostat=ierr,status='replace')
if(ierr.ne.0) stop 'cannot write xopt.xyz'
write(55,'(I5)') nat
write(55,*)
!write(55,'(2F16.8)')energy,gnorm
do i=1,nat
 write(55,'(a2,5x,3(F18.14,3x))') esym(iat(i)), xyz(1,i)*f,xyz(2,i)*f,xyz(3,i)*f
enddo
close(55)
end subroutine wrxyz



integer function kdelta(i,j)
implicit none
integer i,j

kdelta=0
if(i==j) kdelta=1

end function


