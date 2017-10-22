! test program for FastScape Course

program FastScape0

implicit none

! declaring arrays

real, dimension(:), allocatable :: h,a,length
real, dimension(:,:), allocatable :: x,y,z
integer, dimension(:), allocatable :: rec,ndon,stack
integer, dimension(:,:), allocatable :: donor

integer nx,ny,nn,nstep,nfreq,nstack
integer i,j,ij,ii,jj,iii,jjj,ijk,ijr,istep
real xl,yl,dx,dy,dt,k,n,m,u,l,slope,smax
real diff,fact,h0,hp,tol

! defining size of the problem

nx=501;ny=501
nx=501;ny=501
nn=nx*ny

! allocating memory

allocate (h(nn),a(nn),length(nn),rec(nn),ndon(nn),stack(nn),donor(8,nn))
allocate (x(nx,ny),y(nx,ny),z(nx,ny))

! defining geometrical and temporal constants

xl=100.e3;yl=100.e3
dx=xl/(nx-1);dy=yl/(ny-1)
dt=10000.
dt=1000.
nstep=1
!nstep=1000
nfreq=10

! generating initial topography

call random_number (h)
  do j=1,ny
    do i=1,nx
    ij=i+(j-1)*nx
    if (i.eq.1.or.i.eq.nx.or.j.eq.1.or.j.eq.ny) h(ij)=0.d0
    enddo
  enddo

! initializing erosional parameters

k=2.e-6;n=2;m=0.8
u=2.e-3
tol=1.e-3

! begining of time stepping

  do istep=1,nstep

! initializing rec and length

    do ij=1,nn
    rec(ij)=ij
    length(ij)=0.
    enddo

! computing receiver array

    do j=2,ny-1
      do i=2,nx-1
!      do i=1,nx
      ij=i+(j-1)*nx
      smax=tiny(smax)
        do jj=-1,1
          do ii=-1,1
          iii=i+ii
!          iii=modulo(iii-1,nx)+1
          jjj=j+jj
          ijk=iii+(jjj-1)*nx
            if (ijk.ne.ij) then
            l=sqrt((dx*ii)**2+(dy*jj)**2)
            slope=(h(ij)-h(ijk))/l
              if (slope.gt.smax) then
              smax=slope
              rec(ij)=ijk
              length(ij)=l
              endif
            endif
          enddo
        enddo
      enddo
    enddo

! initialising number of donors per node to 0

  ndon=0

! computing donor arrays

    do ij=1,nn
      if (rec(ij).ne.ij) then
      ijk=rec(ij)
      ndon(ijk)=ndon(ijk)+1
      donor(ndon(ijk),ijk)=ij
      endif
    enddo

! computing stack

  nstack=0
    do ij=1,nn
      if (rec(ij).eq.ij) then
      nstack=nstack+1
      stack(nstack)=ij
      call find_stack (ij,donor,ndon,nn,stack,nstack)
      endif
    enddo

! computing drainage area

  a=dx*dy
    do ij=nn,1,-1
    ijk=stack(ij)
      if (rec(ijk).ne.ijk) then
      a(rec(ijk))=a(rec(ijk))+a(ijk)
      endif
    enddo

! adding uplift to landscape

    do j=2,ny-1
      do i=2,nx-1
!      do i=1,nx
      ij=i+(j-1)*nx
      h(ij)=h(ij)+u*dt
      enddo
    enddo

! computing erosion

    do ij=1,nn
    ijk=stack(ij)
    print*,ijk
    ijr=rec(ijk)
      if (ijr.ne.ijk) then
      fact=k*dt*a(ijk)**m/length(ijk)**n
      h0=h(ijk)
      hp=h0
      diff=tol*2.
        do while (abs(diff).gt.tol)
        h(ijk)=h(ijk)-(h(ijk)-h0+ &
        fact*(h(ijk)-h(ijr))**n)/(1.+fact*n*(h(ijk)-h(ijr))**(n-1))
        diff=h(ijk)-hp
        hp=h(ijk)
        enddo
      endif
    enddo

! outputing summary of results every nfreq steps

    if ((istep/nfreq)*nfreq.eq.istep) then
    print*,minval(h),sum(h)/nn,maxval(h)
    endif

  enddo

! output final topography to ascii file

open (7,file='FinalTopo.txt',status='unknown')
  do j=1,ny
  write (7,*) (h(i+(j-1)*nx),i=1,nx)
  enddo
close (7)

! output largest river profile to ascii file

open (7,file='Profile.txt',status='unknown')
i=maxloc(a,1)
  do while (ndon(i).ne.0)
  i=donor(maxloc(a(donor(1:ndon(i),i)),1),i)
  write (7,*) h(i),a(i),(h(i)-h(rec(i)))/length(i)
  enddo
close (7)

end

!----------

! recursive routine to compute the stack

recursive subroutine find_stack	(ij,donor,ndon,nn,stack,nstack)

implicit none

integer ij,k,ijk,nn,nstack
integer donor(8,nn),ndon(nn),stack(nstack)

  do k=1,ndon(ij)
  ijk=donor(k,ij)
  nstack=nstack+1
  stack(nstack)=ijk
  call find_stack (ijk,donor,ndon,nn,stack,nstack)
  enddo

return
end
