
      subroutine getrow( a,lda ,n,m ,i ,row ,l1,l2)
      real a(lda,n) ,row(n)

c     --- gets the ith row of a symmetric banded matrix ---

c      if ( numarg().eq.8 ) then
         k1=l1
         k2=l2
c      else
c         k1=1
c         k2=n
c      endif

         do 7 j=k1,k2
7        row(j)=0.

         do 2 j=0,min0(m,i-1)  
2        row(i-j)=a(m+1-j,i)


         do 3 j=1,min0(m,n-i)  
3        row(i+j)=a(m+1-j,i+j)

      return
      end

      subroutine putval( a,lda ,n,m ,i,j ,val )
      real a(lda,n)

c     --- stores an element into a symmetric band store matrix ---

**?   if ( iabs(i-j).gt.m ) return
**?   if ( i.lt.1 .or. j.lt.1 .or. i.gt.n .or. j.gt.n ) return
      ii=max0(i,j)
c     jj=iabs(i-j)
      a(m+1-iabs(i-j),ii)=val
      return
      end

      subroutine symmsq( a,lda,n,m ,as,ldas,ms )
      real a(lda,n) ,as(ldas,n)
     *  ,r(10000),c(10000)
      common /qqq/ aa(200,200),bb(200,200)

c   --- squares the symmetric matrix a ---

      ms=2*m
         write(6,*)'debug: symmsq: do 1'
         do 1 j=1,ms+1
         do 1 i=1,n
1        as(j,i)=0.

         write(6,*)'debug: symmsq: do 10'
         do 10 i=1,n
         call getrow( a,lda,n,m ,i,r )
      k1=max(1,i-m-1)
      k2=min(n,i+m+1)

         do 10 j=max0(1,i-ms),i
         call getrow( a,lda,n,m ,j,c ,k1,k2)

         s=0.
         write(6,*)'debug: symmsq: do 5'
            do 5 k=k1,min(k2,j+m+1)
cccc        do 5 k=k1,k2
5           s=s+r(k)*c(k)

         call putval( as,ldas,n,ms ,i,j ,s )
10       continue

      return
      end
c   file bandc forfor

      subroutine bandc ( a,lda ,n,m ,info )
      real a(lda,n)
*dir$ vfunction sxmups ]fortran only version.
      real rowi(1000)

c     --- matrix decomposition for banded positive definite matrices
c     --- ( this routine, and it's partner bndchls can be used to
c     --- accomplish the same thing as linpack's spbfa,spbsl -note
c     --- however that bandc&bndchls store d**(-1/2) on the diagnol
c     --- whereas linpack stores d**(1/2)
c     --- note: data is most frequently accessed stride (lda-1)
c     --- avoid lda = l*2**q+1 : q >2

         do 50 i=1,m

c        --- get rowi ---

         j1=max0(m+2-i,1)
            do 10 j=j1,m
10          rowi(j)=-a(j,i)

         nrows=min0(m+1,n+1-i)
            do 30 jj=j1,m
cdir$ ivdep
            do 30 ii=1,min0(jj,nrows)
30          a(m+2-ii,i+ii-1)=a(m+2-ii,i+ii-1)
     *      +rowi(jj)*a(jj+1-ii,i-1+ii)

      info=i
      if(a(m+1,i).le.0. ) go to 60
         a(m+1,i)=1./sqrt(a(m+1,i))
            do 40 ii=1,nrows-1
40          a(m+1-ii,i+ii)=a(m+1-ii,i+ii)*a(m+1,i)

50       continue

      j1=1
      nrows=     m+1

      if ( m.le.64 ) then
         do 100 i=m+1,n-m
        a(m+1,i)=sxmupf( a(1,i),a(1,i),a(m+1,i),lda-1,m )
c         a(m+1,i)=sxmups( loc(a(1,i)),loc(a(1,i)),loc(a(m+1,i))
c     *                    ,lda-1,m )
         info=i
         if ( a(m+1,i).le.0. ) goto 60
         a(m+1,i)=1./a(m+1,i)
cdir$ ivdep
            do 100 ii=1,m
100         a(m+1-ii,i+ii)=a(m+1-ii,i+ii)*a(m+1,i)
      else
         do 150 i=m+1,n-m
c        --- get rowi ---

            do 110 j=j1,m
110         rowi(j)=-a(j,i)

c            do 130 jj=j1,m
c            do 130 ii=1,nrows
c130         if ( jj+1-ii .gt. 0 )
c     *      a(m+2-ii,i+ii-1)=a(m+2-ii,i+ii-1)
c     *      +rowi(jj)*a(jj+1-ii,i-1+ii)
         call mxmup( a(1,i),lda-1 ,rowi,1 ,a(m+1,i),lda-1 ,m )

      info=i
      if(a(m+1,i).le.0. ) go to 60
         a(m+1,i)=1./sqrt(a(m+1,i))
            do 140 ii=1,nrows-1
140         a(m+1-ii,i+ii)=a(m+1-ii,i+ii)*a(m+1,i)

150      continue
      endif

         do 250 i=n+1-m,n

c        --- get rowi ---

         j1=max0(m+2-i,1)
            do 210 j=j1,m
210         rowi(j)=-a(j,i)

         nrows=min0(m+1,n+1-i)
            do 230 jj=j1,m
cdir$ ivdep
            do 230 ii=1,min0(jj,nrows)   
230         if ( jj+1-ii .gt. 0 )
     *      a(m+2-ii,i+ii-1)=a(m+2-ii,i+ii-1)
     *      +rowi(jj)*a(jj+1-ii,i-1+ii)

      info=i
      if(a(m+1,i).le.0. ) go to 60
         a(m+1,i)=1./sqrt(a(m+1,i))
            do 240 ii=1,nrows-1
240         a(m+1-ii,i+ii)=a(m+1-ii,i+ii)*a(m+1,i)

250      continue
      info=0
60    return
      end
      subroutine bndchls( a,lda ,n,m ,x,b )
      real a(lda,n) ,x(n),b(n)

      if ( loc(x).ne.loc(b) ) then
         do 1 j=1,n
1        x(j)=b(j)
      endif

         do 20 k=1,n 
         x(k)=x(k)*a(m+1,k)
cdir$ ivdep
            do 10 i=k+1,min0(k+m,n)
10          x(i)=x(i)-x(k)*a(m+k+1-i,i)
20       continue

         do 40 k=n,1,-1
         x(k)=x(k)*a(m+1,k)
cdir$ ivdep
            do 30 i=1,min0(k-1,m)
30          x(k-i)=x(k-i)-x(k)*a(m+1-i,k)
40       continue

      return
      end

c /// fortran versions of assembly language kernals
c      function sxmups ( pq,pr,py,ldy ,m )
c      integer pq,pr,py
c      real p(*),r(*),q(*)
c      pointer ( pq,q ),( pr,r ),( py,y )
c      sxmups = sxmupf ( q,r,y,ldy ,m )
c      return
c      end
      function sxmupf( q,r,y,ldy ,m )
      real q(ldy,m),r(m),y(ldy,m)
c /// specialized m in range (1,64) only (incorrect results otherwise
         do 10 j=m,1,-1
         do 10 i=1,j
10       y(1,i)=y(1,i) -r(j)*q(j,i)
      sxmupf=y(1,1)
      if ( sxmupf.gt.0 ) sxmupf=sqrt(sxmupf)
      return
      end
      subroutine mxmup ( a,lda ,r,ldr ,y,ldy ,m )
      real a(lda,m),r(ldr,m),y(ldy,m)

         do 10 i=1,m
         do 10 j=1,m
10       if ( j.ge.i ) y(1,i)=y(1,i)+ r(1,j)*a(j,i)

      return
      end
      function ranf()
      real ranf
c
      data init /1325/
            init = mod(3125*init,65536)
            ranf = (init - 32768.0)/16384.0
      return
      end
