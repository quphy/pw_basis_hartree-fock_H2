program hf_pw_h2
  implicit none
  integer, allocatable :: n_G(:,:)
  integer nmax,ncount,nx,ny,nz,i,j,k,l,o
  complex(kind=8), allocatable :: T(:,:),Vn(:,:),ge(:,:,:,:),Hc(:,:),F(:,:),work(:),C(:,:),P(:,:)
  real, allocatable :: ep(:),eps(:),rwork(:)
  real , dimension(3,2) :: R
  real , dimension(2) :: Z
  real pi, E, celllength,minep
  complex zi
  integer :: lwork, info,minint
 !according to your requirment to change the nmax and cellength
 nmax=2
 celllength=5

 ncount=0
 do nx=-nmax,nmax
   do ny=-nmax,nmax
    do nz=-nmax,nmax
      if(nx**2+ny**2+nz**2<=nmax**2) then 
      ncount=ncount+1
      end if
    end do
   end do
  end do
allocate(n_G(ncount,3),T(ncount,ncount),Vn(ncount,ncount),ge(ncount,ncount,ncount,ncount))
allocate(Hc(ncount,ncount),F(ncount,ncount),C(ncount,ncount),P(ncount,ncount),ep(ncount),eps(ncount))
i=1
do nx=-nmax,nmax
   do ny=-nmax,nmax
    do nz=-nmax,nmax
      if(nx**2+ny**2+nz**2<=nmax**2) then 
      n_G(i,1)=nx
      n_G(i,2)=ny
      n_G(i,3)=nz
      i=i+1
      end if
    end do
   end do
  end do

 R(1,1)=0
 R(2,1)=0
 R(3,1)=0.7
 R(1,2)=0
 R(2,2)=0
 R(3,2)=0.7

 T=0
 Vn=0
 Z=1
 pi=dacos(-1.D0)
 zi=(0.0,1.0)
 ep=0
 lwork=-1
 allocate(work(1))
 allocate(rwork(3*ncount-2))

 !calculation of the kinetic energy matrix T and nuclear attraction matrix Vn
 do i =1,ncount
  T(i,i)= 0.5 * (2.0 * pi / celllength)**2 * real(n_G(i,1)**2 + n_G(i,2)**2 + n_G(i,3)**2)
 end do

 !The diagonal elements of the nuclear attraction matrix Vndiverge, but we are concerned with relative energy,
 ! so we can set the diagonal elements to 0.
 do i=1,ncount
  do j=1,ncount
   if (i==j)then
   Vn(i,j)=0
   else 
   do k=1,2
   Vn(i,j)=Vn(i,j)-celllength**2/(pi*((n_G(i,1)-n_G(j,1))**2+(n_G(i,2)-n_G(j,2))**2+(n_G(i,3)-n_G(j,3))**2))*&
   Z(k)*exp(zi*(2*pi/celllength)*((n_G(j,1)-n_G(i,1))*R(1,k)+(n_G(j,2)-n_G(i,2))*R(2,k)+(n_G(j,3)-n_G(i,3))*R(3,k)))
   end do
   end if
  end do
 end do

! the same situation to the 2-e integral，but we can set it to 0 
 do i=1,ncount
  do j=1,ncount
    do k=1,ncount
      do l=1,ncount
        if ((n_G(i,1)+n_G(j,1) == n_G(k,1)+n_G(l,1)) .and. &
            (n_G(i,2)+n_G(j,2) == n_G(k,2)+n_G(l,2)) .and. &
            (n_G(i,3)+n_G(j,3) == n_G(k,3)+n_G(l,3))) then
          if ((n_G(i,1)-n_G(k,1) == 0) .and. (n_G(i,2)-n_G(k,2) == 0) .and. (n_G(i,3)-n_G(k,3) == 0)) then
            ge(i,j,k,l) = 0
          else
            ge(i,j,k,l) = (1) / ((celllength*pi)*real((n_G(i,1)-n_G(k,1))**2 + (n_G(i,2)-n_G(k,2))**2 + (n_G(i,3)-n_G(k,3))**2))
          end if
        else
          ge(i,j,k,l) = 0
        end if
      end do
    end do
  end do
end do

Hc=T+Vn


call zheev('V','U', ncount, Hc, ncount, eps, work, lwork, rwork, info)
if (info /= 0) then
  print*, 'Error in first zheev call, info=', info
  stop
end if
lwork = int(real(work(1)))
deallocate(work)
allocate(work(lwork))


call zheev('V','U', ncount, Hc, ncount, eps, work, lwork, rwork, info)
 C=Hc
 P=0

minep=100.0
do o=1,ncount 
    if (eps(o)<minep)then
    minep=eps(o)
    minint=o
    end if
end do

 do i=1,ncount
   do j=1,ncount
  P(i,j)=P(i,j)+2*C(i,minint)*conjg(C(j,minint))
   end do
  end do

DO WHILE( abs(eps(minint)-ep(minint))>0.000001  )
ep=eps
F=T+Vn
 do i=1,ncount
  do j=1,ncount
   do k=1,ncount
    do l=1,ncount
    F(i,j)=F(i,j)+P(k,l)*(ge(i,k,j,l)-0.5*ge(i,l,j,k))
    end do
   end do
  end do
 end do
call zheev('V','U', ncount,F, ncount, eps, work, lwork, rwork, info)
 C=F
 P=0

 minep=100
do o=1,ncount 
    if (eps(o)<minep)then
    minep=eps(o)
    minint=i
     end if
end do

  do i=1,ncount
   do j=1,ncount
  P(i,j)=P(i,j)+2*C(i,minint)*conjg(C(j,minint))
  end do
 end do
 F=0
end DO
