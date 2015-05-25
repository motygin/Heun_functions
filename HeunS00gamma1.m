% the second local Heun function at z=0
% the case when gamma = 1
%
% evaluation by power series
%
% Usage:
% [val,dval,err,numb,wrnmsg] = HeunS00gamma1(a,q,alpha,beta,delta,z)
%
% |z| should not exceed the convergency radius min{1,|a|}
%
% The function uses a parameter Heun_klimit which can be changed by HeunOpts: 
% Heun_klimit is the maximum number of series' terms
%
% Returned parameters:
% val is the value of the Heun function
% dval is the value of z-derivative of the Heun function
% err is the estimated error
% numb is the number of power series terms needed for the evaluation
% wrnmsg is a warning message:
%   it is empty if computations are ok
%   otherwise it is a diagnostic message
% if wrnmsg is not empty, then the function returns val, dval = NaN
%
% Oleg V. Motygin, copyright 2015, license: GNU GPL v3
%
% 28 April 2015
%
function [val,dval,err,numb,wrnmsg] = HeunS00gamma1(a,q,alpha,beta,delta,z)
  
  global Heun_klimit;
  
  if isempty(Heun_klimit)
    HeunOpts();
  end

  R = min(1,abs(a));
  
  wrnmsg = '';

  if (abs(z)>=R)

    wrnmsg = strcat('HeunS00gamma1: z is out of the convergence radius = ',num2str(R),'; ');
    val = NaN; dval = NaN; err = NaN; numb = NaN;

  else

    epsilon = alpha+beta-delta;
  
    recur0 = @(k,ckm1,ckm2) ...
      (ckm1*z*(q/k+(1-1/k)*((a+1)+epsilon+a*delta)+(k+2/k-3)*(a+1)) - ...
       ckm2*z^2*((1-2/k)*(1+epsilon+delta)+1/k*alpha*beta+k+6/k-5))/(a*k);

    recur1 = @(k,ckm1,ckm2,dkm0,dkm1,dkm2) recur0(k,ckm1,ckm2) + ...
      (-dkm0*a*2 + dkm1*z*((epsilon+a*delta)/k+(a+1)*2*(1-1/k)) + ...
       dkm2*z^2*((4-alpha-beta)/k-2))/(a*k);

    L1 = 0; dL1 = 0; ddL1 = 0; dm1 = 0; dm2 = NaN; ckm0 = NaN; ckm1 = 0; ckm2 = 0; 

    L2 = 1; dL2 = 0; ddL2 = 0; skm2 = 0; skm1 = 1;
    
    dsm1 = 0; dsm2 = NaN; skm0 = NaN;
    
    k = 1;

    while ( (k<=Heun_klimit) && ( (dsm2~=dsm1) || (abs(skm0)>eps) || ...
      (dm2~=dm1) || (abs(ckm0)>eps) ) )

      skm0 = recur0(k,skm1,skm2);
      ckm0 = recur1(k,ckm1,ckm2,skm0,skm1,skm2);

      L1 = L1+ckm0; dL1 = dm1+k*ckm0/z; ddL1 = ddL1+k*(k-1)*ckm0/z^2;
      ckm2 = ckm1; ckm1 = ckm0;
      dm2 = dm1; dm1 = dL1;

      L2 = L2+skm0; dL2 = dsm1+k*skm0/z; ddL2 = ddL2+k*(k-1)*skm0/z^2;
      skm2 = skm1; skm1 = skm0;
      dsm2 = dsm1; dsm1 = dL2;

      k = k+1;

    end
    
    numb = k-1;

    val = L1 + log(z) * L2;
    dval = dL1 + log(z) * dL2 + L2/z;
    ddval = ddL1 - L2/z^2 + 2*dL2/z + log(z) * ddL2;

    if ( isinf(val) || isinf(dval) || isnan(val) || isnan(dval) )

      wrnmsg = 'HeunS00gamma1: failed convergence of recurrence and summation; '; 
      val = NaN; dval = NaN; err = NaN;
      
    else
    
      val2 = ( z*(z-1)*(z-a)*ddval+((z-1)*(z-a)+delta*z*(z-a)+epsilon*z*(z-1))*dval ) / ...
              (q-alpha*beta*z);
      err1 = abs(val-val2);

      if abs(q-alpha*beta*z)<0.01
      
        err2 = abs(ckm0)*sqrt(numb) + abs(L1)*eps*numb + ...
               abs(log(z)) * ( abs(skm0)*sqrt(numb) + abs(L2)*eps*numb );           
        err =  min(err1,err2);
      
      else
      
        err = err1;
      
      end
       
    end

  end

end

