      subroutine allk_fit_run(qinp,qout,eout)
      
      double precision arg1,arg2,arg3,arg4,ar_list(35)
      real fputil,qinp,qout,eout
      
      include 'allk_fit.inc'
      
      integer i,erf
      external allk_fmin
      data lmini/5/lmino/66/lminf/11/linput/21/fputil/0./f_out/49/
      data chimax/200000./
      
      print *, 'ALLK_FIT: Reading slfit_inp.dat'
      open(unit=linput,file='slfit_inp.dat',status='OLD')
      
      read (linput,*) inffile
      read (linput,*) maxcalls
      read (linput,*) tolorance 
      read (linput,*) errorup
      read (linput,*) parms
      read (linput,*) obsv2fit     
      read (linput,*) ebeam
      read (linput,*) wmin,wmax,wstep
      read (linput,*) q2min,q2max,q2step
      read (linput,*) cosmin,cosmax,cosstep
      read (linput,*) phimin,phimax,phistep
      read (linput,*) parofile      
      read (linput,*) datform
      read (linput,*) rndmize
      read (linput,*) fname(1)
      read (linput,*) fname(2)
      read (linput,*) fname(3)
      read (linput,*) fname(4)
     
      close(linput)
      
      chimax = 200000.
      qout = 0.5*(q2min+q2max)
      eout = ebeam
      if (qinp.eq.-4) return
      
      print *, 'ALLK_FIT: Initializing Minuit'     
c     Initialize Minuit, Define I/O
      CALL mninit(lmini,lmino,lminf) 
                
c     Define the title of the Job
      CALL mnseti('Sato-Lee Analysis of Single Pion Production')
          
      print *, 'ALLK_FIT: Defining parameters'
c     Define parameters, Set the initial values
      fix_param=0.
      do i=1,parms
        arg1 = pval(i)
        arg2 = pstp(i)
        arg3 = pmin(i)
        arg4 = pmax(i)
        call mnparm(plist(i),pname(i),arg1,arg2,arg3,arg4,erf,fputil)
        if(erf.ne.0)then
        write(6,*)'Unable to define parameter  ',pname(i),
     +  i,' Error',erf
        endif
      enddo
      
      if (qout.eq.0) pstat(3)=2	! Override fitting status of gdc2 if Q2=0
      
      ar_list(1)=1.
      CALL mnexcm(allk_fmin,'CALL UFC',ar_list,1,erf,fputil)
      
      print *, 'ALLK_FIT: SCAN or FIX parameters'
      do i=1,parms
        ar_list(1)=i
c       ar_list(2)=100
        ar_list(2)=10
        if(abs(pstat(i)).eq.1)then
          print *,'ALLK_FIT: Scaning parameter ',i
          call mnexcm(allk_fmin,'SCAN',ar_list,2,erf,fputil)
        endif
        if(pstat(i).gt.1)then
          print *,'ALLK_FIT: Fixing parameter ',i
          fix_param=fix_param+1.
          ar_list(1)=plist(i)
          CALL mnexcm(allk_fmin,'FIX',ar_list,1,erf,fputil)
          if(erf.ne.0)then
            write(6,*)'Unable to Fix the parameter  ',pname(i),
     +      i,' Error',erf
          endif
        endif
      enddo
   
      print *,'ALLK_FIT: Number of parameters = ',parms
      print *,'ALLK_FIT: Number of fixed parameters = ',fix_param

      fputil=0.
      ar_list(1)=1.
      CALL mnexcm(allk_fmin,'SET PRINT',ar_list,1,erf,fputil)
     
      ar_list(1)=errorup
      CALL mnexcm(allk_fmin,'SET ERRORDEF',ar_list,1,erf,fputil)
     
      ar_list(1)=2.
      CALL mnexcm(allk_fmin,'CALL UFC',ar_list,1,erf,fputil)
     
      ar_list(1)=maxcalls
      ar_list(2)=tolorance
      write(6,*)'ALLK_FIT: Starting fit', ar_list(1), ar_list(2)
      CALL mnexcm(allk_fmin,'MIG',ar_list,2,erf,fputil)
      
      write(6,*),'ALLK_FIT: MIGRAD',erf
      write(lmino,*)'MIGRAD',erf
     
      if(erf.ne.0)then
        ar_list(1)=tolorance
        ar_list(2)=1
        ar_list(3)=2
        ar_list(4)=3
        CALL mnexcm(allk_fmin,'MINOS',ar_list,1,erf,fputil)
        write(6,*),'MINOS',erf
        write(lmino,*),'MINOS',erf
      endif
  
      print *,'ALLK_FIT: Sending STOP to MINUI '
      CALL mnexcm(allk_fmin,'STOP',0,0,erf,fputil)
      
c   77 call hjanr_plot 
c      call hjanr_pick
c      go to 77
      
      end
      
      subroutine allk_fmin(npar,grad,fval,xval,iflag,futil)
      save
      parameter (maxpoints=2000)
      double precision grad,xval,fval 
      dimension grad(37),xval(37)
      integer obsv,fitop,point_stat(maxpoints),debug,f_out
      integer*4 idum
      real xx(maxpoints,10),nb,q2,grn(2)
      logical bad_xsec,bad_data,bad_acc,bad_stat,iupd,tst_obs
      data lmini/5/lmino/66/lminf/11/linput/21/fputil/0./f_out/49/
      data chimax/200000./
      
      include 'mp.inc'
      include 'allk_fit.inc'
      
      print *, 'ALLK_FMIN: iflag',iflag

c Read in data file  

      IF ( iflag .eq. 1 ) THEN
        open(unit=f_out,file=inffile,status='UNKNOWN')
        call timex(timdd2)
        call getunixtime(idum)
        idum       =-idum
        time_add   = 0.
        fitop      = 1
        itr_number = 0
        itr_count  = 0
        Npoints_max= 0
        N_zero     = 0
        N_bad      = 0
        N_good     = 0
        nhic       = 0
        
        do if=1,4
        if (fname(if).ne.'JUNK'.and.fname(if).ne.'junk') then
          print *, 'ALLK_FMIN: Opening file:',fname(if)
          
          open(unit=44,file=fname(if),status='OLD')
          
          do i=1,maxpoints
            
            if(datform.eq.1)then
            
              READ(44,*,err=999,end=998)W,Q2,epsilon,costh,phi,sec,err1
              bad_data    = sec.le.0..and.err1.le.0
                            
            elseif(datform.eq.2)then	! e1e data file
            
              READ(44,*,err=999,end=998)W,Q2,epsilon,costh,phi,sec,err1,
     +        nb,ac,rad,dum,xx1,xx2
              if (costh.gt.1.0) costh=cos(costh*3.1415926/180.)
              bad_xsec  = sec.le.0.or.err1.le.0
c              bad_acc  = ac.le.0.05.or.(ac.lt.0.1.and.nb.lt.10)
              bad_acc   = ac.le.0.05
              tst_obs   = pi_type.eq.1.and.obsv.eq.1
              bad_data	= tst_obs.and.(bad_xsec.or.bad_acc)
              pi_type   = xx1
              obsv      = xx2
              if (rndmize.eq.1) then
                print *, 'ALLK_FMIN: Randomizing data points'
                rr1 = ran(idum)
                rr2 = ran(idum)
                z = sin(2*3.14159*rr1)*sqrt(-2.*alog(rr2))
                sec = sec + err1*z
              endif
            elseif(datform.eq.3)then	! e1b data file

              READ(44,*,err=999,end=998)W,Q2,epsilon,costh,phi,sec,err1,
     +        nb,ac,rad,dum,pi_type,obsv
     
              bad_xsec     = sec.le.0..and.err1.le.0
              bad_acc      = ac.lt.0.05.or.(ac.lt.0.1.and.nb.lt.10)
              bad_data	   = pi_type.eq.1.and.(bad_xsec.or.bad_acc)
                            
            endif
            
            Npoints_max  = Npoints_max+1
            
            if(bad_data)then
              point_stat(Npoints_max) = 1
              N_zero = N_zero+1
            else
              point_stat(Npoints_max) = 0
            endif
              
            if(w.lt.wmin.or.w.gt.wmax)             point_stat(Npoints_max)=1
            if(q2.lt.q2min.or.q2.gt.q2max)         point_stat(Npoints_max)=1
            if(costh.lt.cosmin.or.costh.gt.cosmax) point_stat(Npoints_max)=1
            if(phi.lt.phimin.or.phi.gt.phimax)     point_stat(Npoints_max)=1
            
            if (.not.btest(obsv2fit,pi_type-1))    point_stat(Npoints_max)=1
            if (.not.btest(obsv2fit,obsv+1))       point_stat(Npoints_max)=1
            
            xx(Npoints_max,1)     = w
            xx(Npoints_max,2)     = q2
            xx(Npoints_max,3)     = epsilon
            xx(Npoints_max,4)     = costh
            xx(Npoints_max,5)     = phi
            xx(Npoints_max,6)     = sec
            xx(Npoints_max,7)     = err1
            xx(Npoints_max,8)     = pi_type
            xx(Npoints_max,9)     = obsv
            
            if (point_stat(Npoints_max).eq.0) N_good = N_good+1
            
          enddo   
        endif
         
 998    close(44)
 
        enddo
          
        OPEN (UNIT=53,FILE='start_tab',status='UNKNOWN')
          
        write(lmino,*)'Number of data points =  ',Npoints_max
        write(lmino,*)'Number of good data points =  ',Npoints_max-N_zero
        write(lmino,*)'Number of fitted data points =  ',N_good
          
        write(f_out,*)'Kinematics:   W  =',wmin,' - ',wmax
        write(f_out,*)'              Q2 =',q2min,' - ',q2max
        write(f_out,*)'Number of data points =  ',Npoints_max
        write(f_out,*)'Number of good data points =  ',Npoints_max-N_zero
        write(f_out,*)'Number of fitted data points =  ',N_good
          
        write(lmino,*)'Number of data points for first iteration =  ',Npoints_max-N_zero 
          
        write(6,111)(pname(i),xval(i),i=1,parms)
111     format('Starting parameter values  ',/,
     +  37(a9,'=',f10.4/))
     
      ENDIF	! IFLAG.eq.1
      
      IF ( iflag .eq. 2 ) THEN
        print *,'Changing option'
        fitop=0
      ELSEIF ( iflag .EQ. 3 ) THEN
        write(f_out,*)'Number of data points to fit =  ',ndf
        write(f_out,*)'Number of iterations = ',itr_number,'  Chi2/point',fval_last/float(ndf)
        rq2=0.5*(q2min+q2max)
        call prterr(rq2)
        itr_number=0
      ENDIF

      IF(iflag.eq.1.or.iflag.ge.3)THEN
        timdd1   = timdd2
        call timex(timdd2)
        timdd    = timdd2-timdd1
        time_add = time_add+timdd
        write(6,*)
        write (*,*)'Iteration time = ',timdd,'  Total time = ',time_add
        fitted_points=0
        fval=0.
        
        call parput(xval)
        
        iupd = .true.
        
        DO i_point=1,Npoints_max
        
          if(point_stat(i_point).eq.0)then
           
            W       = xx(i_point,1)
            Q2      = xx(i_point,2)
            epsilon = xx(i_point,3)
            costh   = xx(i_point,4)
            phi     = xx(i_point,5)
            pi_type = xx(i_point,8)
            obsv    = xx(i_point,9)
            
            call allk_run(q2,w,epsilon,costh,phi,pi_type,iupd,0)
            
            IF(ierr.eq.0) THEN
            
              var_calc = robs(obsv,pi_type)
              c2weight = 1.
                            
              diff=c2weight*((var_calc-xx(i_point,6))/xx(i_point,7))**2
              if (debug.eq.1) write(6,120) obsv,w,costh,var_calc,xx(i_point,6),diff
120      format(i3,5f10.3)
              
              if(diff.gt.chimax.and.itr_number.eq.0.and.fitop.eq.0)then
                nhic = nhic+1
                point_stat(i_point) = 1
                write(55,134)nhic,i_point,W,Q2,costh,phi,xx(i_point,6),var_calc,diff
                itr_bad = nhic
              else
                if(itr_number.eq.0.and.(iflag.eq.1.or.iflag.eq.3))then
                write(53,133)i_point,W,Q2,epsilon,costh,phi,xx(i_point,6),
     +xx(i_point,7),var_calc,diff
                endif
                
                fitted_points = fitted_points + 1.
                fval          = fval + diff
                fval_last     = fval
              endif
            Else
              write(6,*)'Bad calculation',w,q2,costh,phi,i_point
            ENDIF
          endif
        ENDDO

c *** lcs begin hjanr plugin for interactive plots        
c        do i=1,30
c          xnew(i)=sngl(xval(i))
c        enddo
c        call hjanr_plot 
c *** lcs end hjanr plugin for interactive plots        
        
        if(itr_number.eq.0)then
          if(iflag.eq.1) then
            close(53)
            ndf=fitted_points
          else
            ndf=fitted_points-(parms-fix_param)
          endif
          write(lmino,*)'Number of degrees of freedom for fit =  ',ndf
          print *, 'ALLK_FMIN: Number of degrees of freedom for fit =  ',ndf
        endif
         
        if(iflag.ne.1)itr_number=itr_number+1
        if(iflag.ne.1)itr_count=itr_count+1
        
        write(6,*)'Iteration = ',itr_number,' Value to minimize',fval,' Reduced Chi2',fval/float(ndf)
        write(6,*)'No. fitted param.',int(parms-fix_param),' No. fitted points',int(fitted_points)
        write(26,*)'Iteration = ',itr_number,' Value to minimize',fval,' Reduced Chi2',fval/float(ndf)
        
        debug = 0
        if (itr_count.eq.50) then
          itr_count=0
          debug=1
        endif
        
        do ip=1,parms
          if((pstat(ip).le.0.and.fitop.eq.0).or.(abs(pstat(ip)).eq.1.and.fitop.eq.1))then
            write(6,*)pname(ip),'  ',sngl(xval(ip))
            write(26,*)pname(ip),'  ',xval(ip)
          endif
        enddo
      
      ENDIF
 133  format(i7,8f9.4,f10.5)
 134  format(2i7,7f8.3)
        
      RETURN
      
 999  WRITE(6,*) 'Error while reading the data file',W,Q2,epsilon,costh,phi,sec,err1,
     +        nb,ac,rad,dum,pi_type,obsv      
      end
      
      subroutine prterr(rq2)

      double precision fmin,fedm,errdef,eplus,eminus,eparab,globcc,parv,parerr,bn1,bn2
      integer npari,nparx,istat,i,iv,f_out
      character*10 pnamei(35)
      
      include 'allk_fit.inc'

      OPEN (UNIT=27,FILE=parofile,STATUS='UNKNOWN')
      
      call mnstat(fmin,fedm,errdef,npari,nparx,istat)
      write(6,*)
      write(6,*) 'Minimized function=',fmin,'Reduced Chi2=',fmin/ndf
      write(6,*) 'Estimated distance to minimum',fedm
      write(6,*)
      write(6,*) 'Parameters and errors'

      do i=1,parms
         call mnpout(i,pnamei(i),parv,parerr,bn1,bn2,iv)
         if(pstat(i).le.0)then
            call mnerrs(-i,eplus,eminus,eparab,globcc)
            write(6,*)i,'  ',pnamei(i),parv,parerr,bn1,bn2
            write(6,*)i,eplus,eminus,eparab,globcc
            write(f_out,*)i,'  ',pnamei(i),parv,parerr
            write(27,53)plist(i),pname(i),parv,parerr,
     +         pstat(i),pmin(i),pmax(i)
         else
            write(27,53)plist(i),pname(i),pval(i),pstp(i),
     +         pstat(i),pmin(i),pmax(i)
         endif
      enddo
      
      close(27)
      
53    format(i2,1x,A10,2f10.4,i10,2f10.4)

      return
      end

