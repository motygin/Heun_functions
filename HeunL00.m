% local Heun function, by power series at z=0 for Hl(0)=1, Hl'(0)=q/(a*gamma)
% the case when gamma is not equal to 0, -1, -2, ...
%
% Usage:
% [val,dval,err,numb,wrnmsg] = HeunL00(a,q,alpha,beta,gamma,delta,z)
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
function [val,dval,err,numb,wrnmsg] = HeunL00(a,q,alpha,beta,gamma,delta,z)
  
  global Heun_klimit;
  
  if isempty(Heun_klimit)
    HeunOpts();
  end

  R = min(1,abs(a));
  
  epsilon = alpha+beta+1-gamma-delta;
  
  wrnmsg = '';

  isnonpositivegamma = abs(ceil(gamma-5*eps)+abs(gamma))<5*eps;

  if isnonpositivegamma

    wrnmsg = 'HeunL00: gamma is a non-positive integer, use HeunL00log; ';
    val = NaN; dval = NaN; err = NaN; numb = NaN;

  elseif (abs(z)>=R)

    wrnmsg = strcat('HeunL00: z is out of the convergence radius = ',num2str(R),'; ');
    val = NaN; dval = NaN; err = NaN; numb = NaN;

  elseif (z==0)
  
    val = 1; dval = q/(a*gamma);
    err = 0; numb = 1;
  
  else

    recur = @(k,ckm1,ckm2) ...
      (ckm1*z*(q/k+(1-1/k)*((a+1)*gamma+epsilon+a*delta)+(k+2/k-3)*(a+1))-...
       ckm2*z^2*((1-2/k)*(gamma+epsilon+delta)+1/k*alpha*beta+k+6/k-5))/(a*gamma+(k-1)*a);

    ckm2 = 1; ckm1 = z*q/(a*gamma);

    val = ckm2+ckm1; vm1 = val; vm2 = NaN;
    dval = q/(a*gamma); dm1 = dval; dm2 = NaN;
    ddval = 0;
    
    k = 2; ckm0=1;

    while ( (k<=Heun_klimit) && ( (vm2~=vm1) || (dm2~=dm1) || (abs(ckm0)>eps) ) )
      
      ckm0 = recur(k,ckm1,ckm2);
      val = val+ckm0; dval = dm1+k*ckm0/z;
      ddval = ddval+k*(k-1)*ckm0/z^2;
      ckm2 = ckm1; ckm1 = ckm0;
      vm2 = vm1; vm1 = val;
      dm2 = dm1; dm1 = dval;
      k = k+1;
    
    end

    numb = k-1;

    if ( isinf(val) || isinf(dval) || isnan(val) || isnan(dval) )
    
      wrnmsg = 'HeunL00: failed convergence of recurrence and summation; '; 
      val = NaN; dval = NaN; err = NaN;
      
    else
    
      val2 = ( z*(z-1)*(z-a)*ddval+(gamma*(z-1)*(z-a)+delta*z*(z-a)+epsilon*z*(z-1))*dval ) / (q-alpha*beta*z);
      err1 = abs(val-val2);
      
      if abs(q-alpha*beta*z)<0.01
        err2 = abs(ckm0) * sqrt(numb) + eps * numb * abs(val);
        err =  min(err1,err2);
      else
        err = err1;
      end
      
    end

  end

end
