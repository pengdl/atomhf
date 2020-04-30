Program AtomHF ! Program for atomic closed-shell Hatree-Fock calculation
  implicit real(8) (a-h,o-z)
  character(250) :: tempc
  integer,parameter :: ml=4,maxln=999,mn=29,maxscf=39
  integer :: st(0:ml),nbas(0:ml),nocc(0:ml)
  real(8) :: ex(maxln),norm(maxln),occ(9,0:ml),fs(0:mn),fd(mn)
  real(8),allocatable :: cof(:,:,:),mas(:,:),mat(:,:),mavc(:,:),mavn(:,:),mavx(:,:)
  real(8),allocatable :: maf(:,:),mad(:,:),dold(:,:),ff(:,:),ss(:,:),w(:),work(:)
  sqrpi = dsqrt(atan(1.d0))*2
  fs(0) = 1.d0
  fd(1) = 1.d0
  do i=1,mn
    fs(i) = fs(i-1) * i
    if(i.ge.3.and.mod(i,2).eq.1) fd(i) = fd(i-2) * i
  end do
  read(*,*) tempc,nuc
  read(*,*) tempc,lmax
  j = 0
  do l=0,lmax
    st(l) = j
    read(*,*) tempc,nbas(l),m
    do i=1,nbas(l)
      j = j + 1
      read(*,*) ex(j)
      norm(j) = dsqrt(2**(l+2) * dsqrt((2*ex(j))**(2*l+3)) / fd(2*l+1) / sqrpi)
    end do
    nocc(l) = (m+4*l+1) /(4*l+2)
    occ(1:nocc(l)-1,l) = 2*l+1
    occ(nocc(l),l) = m*0.5d0 - (2*l+1)*(nocc(l)-1)
  end do
  nf = j 
  allocate(cof(0:lmax,0:lmax,0:lmax),mas(nf,nf),mat(nf,nf),mavn(nf,nf),mavc(nf,nf),mavx(nf,nf),maf(nf,nf),mad(nf,nf),dold(nf,nf))
  do l1=0,lmax
    do l2=0,lmax
      do k=0,min(l1,l2)
        l3 = abs(l1-l2)+k+k
        cof(k,l1,l2) = fs((l1+l2+l3)/2)**2 * fs(l1+l2-l3) * fs(l2+l3-l1) * fs(l3+l1-l2) / fs(l1+l2+l3+1) &
                     / fs((l1+l2-l3)/2)**2 / fs((l2+l3-l1)/2)**2 / fs((l3+l1-l2)/2)**2 
      end do
    end do
  end do
  do i=0,lmax
    do k1=st(i)+1,st(i)+nbas(i)
      do k2=st(i)+1,st(i)+nbas(i)
        mas(k2,k1) =  fd(2*i+1)*sqrpi/dsqrt(2**(2*i+4)*(ex(k1)+ex(k2))**(2*i+3))*norm(k1)*norm(k2)
        mat(k2,k1) = (fd(2*i+1)*sqrpi/dsqrt(2**(2*i+4)*(ex(k1)+ex(k2))**(2*i+3))*ex(k2)*(2*i+3) & 
         - fd(2*i+3)*sqrpi/dsqrt(2**(2*i+6)*(ex(k1)+ex(k2))**(2*i+5))*2*ex(k2)*ex(k2))*norm(k1)*norm(k2)
        mavn(k2,k1) = -nuc*fs(i)/(ex(k1)+ex(k2))**(i+1)/2*norm(k1)*norm(k2)
      end do
    end do
  end do
  maf = mat + mavn
  eold = 0.d0
  do iter=1,maxscf
    mad = 0.d0
    do l=0,lmax
      n = nbas(l)
      allocate(ff(n,n),ss(n,n),w(n),work(n*8))
      ff(1:n,1:n) = maf(st(l)+1:st(l)+n,st(l)+1:st(l)+n)
      ss(1:n,1:n) = mas(st(l)+1:st(l)+n,st(l)+1:st(l)+n)
      call dsygv(1,'V','L',n,ff,n,ss,n,w,work,n*8,info)
      do i1=1,n
        do i2=1,n
          do k=1,nocc(l)
            mad(st(l)+i2,st(l)+i1) = mad(st(l)+i2,st(l)+i1) + ff(i2,k)*occ(k,l)*ff(i1,k)
          end do
        end do
      end do
      deallocate(ff,ss,w,work)
    end do
    if(iter.ge.2) mad = 0.85d0*mad + 0.15d0*dold
    dold = mad
    mavc = 0.d0
    mavx = 0.d0
    etot = 0.d0
    do l1=0,lmax
      do j1=st(l1)+1,st(l1)+nbas(l1)
        do j2=st(l1)+1,st(l1)+nbas(l1)
          do l2=0,lmax
            do j3=st(l2)+1,st(l2)+nbas(l2)
              do j4=st(l2)+1,st(l2)+nbas(l2)
                zz = norm(j1) * norm(j2) * norm(j3) * norm(j4) * sqrpi * mad(j3,j4)
                mavc(j1,j2) = mavc(j1,j2) + zz*2*g2e(l2*2,l1*2,0,ex(j3)+ex(j4),ex(j1)+ex(j2),fs,fd,mn)
                do i=0,min(l1,l2)
                  mavx(j1,j2)=mavx(j1,j2)-zz*cof(i,l1,l2)*g2e(l1+l2,l1+l2,abs(l1-l2)+i+i,ex(j3)+ex(j2),ex(j1)+ex(j4),fs,fd,mn)
                end do
              end do
            end do
          end do
          etot = etot + (mat(j1,j2)*2+mavn(j1,j2)*2+mavc(j1,j2)+mavx(j1,j2))*mad(j1,j2)
        end do
      end do
    end do
    maf = mat + mavn + mavc + mavx
    write(6,"(i4,a,f17.9,d13.3)") iter," : ",etot,dabs(etot-eold)
    if(dabs(etot-eold).lt.1.d-9) exit
    eold = etot
  end do
  deallocate(cof,mas,mat,mavn,mavc,mavx,maf,mad,dold)
End Program AtomHF
!----------------------------------------------------------------------|
real(8) function g2e(n1,n2,l,a1,a2,fs,fd,m)
  implicit real(8) (a-h,o-z)
  real(8) :: fs(0:m),fd(m)
  g2e = ri((n1+l)/2+1,(n2-l)/2,a1,a2,fs,fd,m) + ri((n2+l)/2+1,(n1-l)/2,a2,a1,fs,fd,m)
end function g2e
!----------------------------------------------------------------------|
real(8) function ri(n1,n2,a1,a2,fs,fd,m)
  implicit real(8) (a-h,o-z)
  real(8) :: c(0:n2),fs(0:m),fd(m)
  q = a2/(a1+a2)/2
  c(0) = 1.d0
  do k=1,n2
    c(k) = c(k-1) * q
  end do
  do k=0,n2
    c(k) = c(k) * fd(n1+n1+k+k-1) / fs(k)
  end do
  ri = sum(c)*fs(n2)/(2**(n1+2)*a2**(n2+1)*(a1+a2)**(n1)*dsqrt(a1+a2))
end function ri
