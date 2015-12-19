      integer parms,datform,pi_type,obsv2fit,plist,pstat,fix_param
      integer rndmize
      real q2min,q2max,q2step,pval,pstp,pmin,pmax
      character pname*20,parofile*20,inffile*20,fname(4)*60
      common /parc_out/wmin,wmax,wstep,
     &q2min,q2max,q2step,cosmin,cosmax,cosstep,phimin,phimax,phistep,
     &parms,datform,rndmize,obsv2fit,fix_param,parofile,inffile,fname
      common / fit /pval(37),pstp(37),pmin(37),pmax(37),
     &              plist(37),pstat(37),pname(37)