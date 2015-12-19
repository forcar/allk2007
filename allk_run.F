      
      subroutine allk_run(q2,w,eps,costh,phi,k,iupd,iwndx)
      
      include 'mp.inc'
      real m_p,mpi2,m_pi(2),q2,qcm
      real*8 win,qin
      logical ok_wsave,ok_qsave,ok_csave,skip,iupd
      data m_p/0.93828/
      data m_pi/0.1349764,0.13957/
      data pi/3.14159/
      save wsave,qsave,csave
      
      if (iupd) then
        wsave = 0.0
        qsave = 0.0
        csave = -2
        iupd = .false.
      endif
      
      ok_wsave	= w.eq.wsave
      ok_qsave	= q2.eq.qsave
      ok_csave 	= costh.eq.csave
      skip      = ok_wsave.and.ok_qsave
      
      wsave 	= w
      qsave	= q2
      csave	= costh
      
      win        = w*1000
      qin	 = q2
      
      if (.not.skip) then
        if (iwndx.eq.0) then	! Recalculate multipoles
          call allkqw(win,qin,1)
          n=1
        else			! Use precalculated table
          n=iwndx
        endif
        call mul_amp(n,k)
      endif
      
      skip = skip.and.ok_csave
      if (skip) go to 10
      
      mpi2 = m_pi(k)**2
      ecm  = (w*w + mpi2 - m_p**2)/(2.*w)
      pcm  = ecm**2 - mpi2
      pcm  = sqrt(pcm)
      fkt  = 2.0*w*pcm/(w**2-m_p**2)
      qcm  = ((w**2 + q2 + m_p**2)/2.0/w)**2 - m_p**2
      ekin = sqrt(q2/qcm)
      
      vl   = eps
      vtt  = eps
      vlt  = sqrt(2*eps*(1+eps))
      vltp = sqrt(2*eps*(1-eps))
            
      call hel_amp(costh)
      
      rsp(1) = fkt*sig_t
      rsp(2) = fkt*sig_l*vl*ekin**2
      rsp(3) = fkt*sig_tt
      rsp(4) = fkt*sig_lt
      rsp(5) = fkt*sig_ltp 
      
10    phir = phi*pi/180.0      
      cs2f = cos(2*phir)
      csf  = cos(phir)
      snf  = sin(phir)
      
      sigma_u  =     rsp(1)  
     1         +     rsp(2)  
     1         + vtt*rsp(3)*cs2f
     1         + vlt*rsp(4)*ekin*csf 
     
      sigma_lt  =  vlt*rsp(4)*ekin*csf 
      sigma_ltp = vltp*rsp(5)*ekin*snf 
      
      robs(1,k) = sigma_u
      robs(2,k) = sigma_ltp/sigma_u
      robs(3,k) = sigma_lt/sigma_u
      robs(4,k) = -sig_tt/sig_t
      robs(5,k) = rsp(1)
      robs(6,k) = rsp(1)+rsp(2)
      robs(7,k) = rsp(3)
      robs(8,k) = rsp(4)*ekin
      robs(9,k) = rsp(5)*ekin
      robs(10,k) = (rsp(1)+rsp(2))/fkt	! Bates R0
      robs(11,k) = fkt*(sig_t32-sig_t12)
      robs(12,k) =  sigt_t_y/sig_t	! T
      robs(13,k) = -sigt_tt_y/sig_t	! P
      robs(14,k) = -sigt_tt_z/sig_t	! G
      robs(15,k) = -sigt_tt_x/sig_t	! H
      return
      end
      
      subroutine mul_amp(n,k)
      
      implicit real*8(a-b,d-h,o-z)
      implicit complex*16(c)
      parameter(ndim=20,maxpot=8,nedim=50,npdim=20,nbmas=10,
     &          ladim=10)
      common /mult/ cdmul(nedim,4,0:ladim,8),cbmul(nedim,4,0:ladim,8)
     & ,cbrmul(nedim,4,0:ladim,8),cbdmul(nedim,4,0:ladim,8)
     & ,ctmul (nedim,4,0:ladim,8)
     & ,cxb(4,0:ladim,8),cxbr(4,0:ladim,8),cxbd(4,0:ladim,8)
     & ,cxd(4,0:ladim,8),cxt(4,0:ladim,8)
      common / c2001 / wcmmev(1000),qsgev(1000),fact(1000),
     &ebeam,qmin,wmin,delq,delw,meshq,meshw,imaxdt
      include 'mp.inc'
      
      factor = 0.141383e0     ! Convert from 10**-3/m_pi+ to sqrt(microbarns)
      
      do i=0,5
        sp(i)=ctmul(n,k,i,7)*fact(n)*factor
        sm(i)=ctmul(n,k,i,8)*fact(n)*factor
        ep(i)=ctmul(n,k,i,1)*fact(n)*factor
        em(i)=ctmul(n,k,i,2)*fact(n)*factor
        mp(i)=ctmul(n,k,i,3)*fact(n)*factor
        mm(i)=ctmul(n,k,i,4)*fact(n)*factor
      enddo
      
c      c12=(ctmul(n,1,1,4)-sqrt(2.)*ctmul(n,2,1,4))/3.
c      c32=(ctmul(n,1,1,4)+ctmul(n,2,1,4)/sqrt(2.))
c      c32=c32*0.8
c      mm(1)=(c12+2*c32/3.)*fact(n)*factor

      end
      
      subroutine hel_amp(x)

      include 'mp.inc'

      real x
      data root2/1.41421356/

      pol(0,1) = 0.
      pol(1,1) = 1.
      pol(2,1) = 3.*x 
      pol(3,1) = (15.*x*x-3.)/2.
      pol(4,1) = (35.*x*x*x-15.*x)/2.
      pol(5,1) = (315.*x*x*x*x-210.*x*x+15.)/8.
      pol(6,1) = (693.*x*x*x*x*x-630.*x*x*x +105.*x)/8.
      pol(7,1) = (3003.*x*x*x*x*x*x-3465.*x*x*x*x+945*x*x-35.)/16.
      pol(0,2) = 0.
      pol(1,2) = 0.
      pol(2,2) = 3.
      pol(3,2) = 15.*x 
      pol(4,2) = (105.*x*x-15.)/2.
      pol(5,2) = (315.*x*x*x-105.*x)/2.
      pol(6,2) = (3465.*x*x*x*x-1890.*x*x+105.)/8.  
      pol(7,2) = (9009.*x*x*x*x*x-6930.*x*x*x +945.*x)/8.
      
      ff1 = 0.0
      ff2 = 0.0
      ff3 = 0.0
      ff4 = 0.0
      ff5 = 0.0
      ff6 = 0.0
      
      do l = 0, mwave_L
        if(l.lt.2) then
          ff1 = ff1 + (    l*mp(l) +  ep(l))*pol(l+1,1)
        else
          ff1 = ff1 + (    l*mp(l) +  ep(l))*pol(l+1,1)
     1              + ((l+1)*mm(l) +  em(l))*pol(l-1,1)
        endif
      enddo

      do l = 1, mwave_L
        ff2 = ff2 + ((l+1)*mp(l) + l*mm(l))*pol(l,1)
        if(l.lt.2) then
           ff3 = ff3 + (   ep(l) -   mp(l))*pol(l+1,2)
        else
           ff3 = ff3 + (   ep(l) -   mp(l))*pol(l+1,2) +
     1                 (   em(l) +   mm(l))*pol(l-1,2)
        endif
      enddo

      do l = 2, mwave_L
        ff4 = ff4 + (mp(l) - ep(l) -mm(l) - em(l))* pol(l,2)
      enddo

      do l = 0, mwave_L
        if(l.lt.2) then
          ff5 = ff5 +  (l+1)*sp(l)*pol(l+1,1)
        else
          ff5 = ff5 + ((l+1)*sp(l)*pol(l+1,1) - l*sm(l)*pol(l-1,1))
        endif
      enddo

      do l = 1, mwave_L
        ff6 = ff6 + (l*sm(l) - (l+1)*sp(l))*pol(l,1)
      enddo

      theta_cm	= acos(x)
      s		= sin(theta_cm)
      c		= cos(theta_cm)
      s2	= sin(theta_cm/2.0)
      c2	= cos(theta_cm/2.0)
      
      hh1 = -s*c2*(ff3 + ff4)/root2
      hh2 =  c2*((ff2 - ff1) + 0.5*(1 - c)*(ff3 - ff4))*root2
      hh3 =  s*s2*(ff3 - ff4)/root2
      hh4 =  s2*((ff2 + ff1) + 0.5*(1 + c)*(ff3 + ff4))*root2
      hh5 =  c2*(ff6 + ff5)
      hh6 =  s2*(ff6 - ff5)
      
      sig_t  = (cabs(hh1)**2+cabs(hh2)**2
     1         +cabs(hh3)**2+cabs(hh4)**2)/2.
      sig_l  =  cabs(hh5)**2+cabs(hh6)**2
      sig_tt =  real(-conjg(hh1)*hh4+conjg(hh2)*hh3)
      
      sig_lt =    real(conjg(hh5)*(hh1-hh4)  + 
     1                 conjg(hh6)*(hh2+hh3))/root2      
      sig_ltp = -aimag(conjg(hh5)*(hh1-hh4)  + 
     1                 conjg(hh6)*(hh2+hh3))/root2
     
      sig_t12 = (cabs(hh2)**2+cabs(hh4)**2)
      sig_t32 = (cabs(hh1)**2+cabs(hh3)**2)
      
      sigt_t_y   = -aimag(conjg(hh1)*hh2+conjg(hh3)*hh4)	! T
      sigt_tt_y  = -aimag(conjg(hh1)*hh3+conjg(hh2)*hh4)	! P
      sigt_tt_z  = -aimag(conjg(hh1)*hh4+conjg(hh2)*hh3)	! G
      sigt_tt_x  = -aimag(conjg(hh1)*hh3-conjg(hh2)*hh4)	! H
           
      end
      
      subroutine sigtot(sigt,sig32)
      
      include 'mp.inc'
      
      sum1 = 0.
      sum2 = 0.
      do l = 0, mwave_L
        a=(l+2)*(cabs(ep(l))**2+cabs(mm(l+1))**2)
        b=    l*(cabs(mp(l))**2+cabs(em(l+1))**2)
        c=real(conjg(ep(l))*mp(l))-real(conjg(em(l+1))*mm(l+1))
        sigt=(a+b)*(l+1)**2
        sig32=(b-a-2*l*(l+2)*c)*(l+1)
        sum1=sum1+sigt
        sum2=sum2+sig32
      enddo
      sigt = sum1*fkt*2*3.14159
      sig32= sum2*fkt*4*3.14159
      
      end

      function ran(idum)
c     Random number generator used because I can't find one in the
c     library.

      implicit none
      integer*4 idum
      integer*4 mbig,mseed, mz
      real ran,fac
      parameter (mbig=1000000000, mseed=161803398,mz=0,fac=1./mbig)
      integer*4 i,ii,inext,inextp,k
      integer*4 mj,mk,ma(55)
      save inext,inextp,ma

c     Initialization section:
      if (idum .lt. 0.)then
         mj=mseed-idum
         mj=mod(mj,mbig)
         ma(55)=mj
         mk=1
         do  i=1,54
            ii=mod(21*i,55)
            ma(ii)=mk
            mk=mj-mk
            if(mk .lt. mz)mk=mk+mbig
            mj=ma(ii)
         enddo
         do k=1,4
            do i=1,55
               ma(i)=ma(i)-ma(1+mod(i+30,55))
               if(ma(i) .lt. mz)ma(i)=ma(i)+mbig
            enddo
         enddo
         inext=0
         inextp=31
         idum=1
      endif
 25   inext=inext+1
      if(inext .eq. 56)inext=1
      inextp=inextp+1
      if(inextp .eq. 56)inextp=1
      mj=ma(inext)-ma(inextp)
      if(mj .lt. mz)mj=mj+mbig
      ma(inext)=mj
      idum=mj
      ran=mj*fac
      if (ran .eq. 0. .or. ran .eq. 1.)go to 25
      if (ran .lt. 0. .or. ran .gt. 1.)then
         write(6,*)' random error, ran =',ran
         go to 25
      endif

      return
      end
