% local Heun function, by power series at z=0, for Hl(0)=1
% the case when gamma = 0, -1, -2, ...
%
% Usage:
% [val,dval,err,numb,wrnmsg] = HeunL00log(a,q,alpha,beta,gamma,delta,z)
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
function [val,dval,err,numb,wrnmsg] = HeunL00log(a,q,alpha,beta,gamma,delta,z)
  
  global Heun_klimit;
  
  if isempty(Heun_klimit)
    HeunOpts();
  end

  R = min(1,abs(a));
  
  wrnmsg = '';

  isnonpositivegamma = abs(ceil(gamma-5*eps)+abs(gamma))<5*eps;

  if ~isnonpositivegamma

    wrnmsg = 'HeunL00log: gamma is not a non-positive integer; ';
    val = NaN; dval = NaN; err = NaN; numb = NaN;

  elseif (abs(z)>=R)

    wrnmsg = strcat('HeunL00log: z is out of the convergence radius = ',num2str(R),'; ');
    val = NaN; dval = NaN; err = NaN; numb = NaN;

  else

    epsilon = alpha+beta+1-gamma-delta;
  
    N = 1-gamma;
  
    recur0 = @(k,ckm1,ckm2) ...
      (ckm1*z*(q/k+(1-1/k)*((a+1)*gamma+epsilon+a*delta)+(k+2/k-3)*(a+1)) - ...
       ckm2*z^2*((1-2/k)*(gamma+epsilon+delta)+1/k*alpha*beta+k+6/k-5))/(a*gamma+(k-1)*a);

    recur1 = @(k,ckm1,ckm2,dkm0,dkm1,dkm2) recur0(k,ckm1,ckm2) + ...
      (dkm0*a*(1-gamma-2*k) + dkm1*z*(epsilon+a*delta+(a+1)*(gamma+2*k-3)) + ...
       dkm2*z^2*(4-2*k-alpha-beta))/(a*k*(gamma+k-1));

    L1 = 1; dL1 = 0; ddL1 = 0; ckm0 = 1; ckm1 = 1; ckm2 = 0;
    
    for k=1:N-1   
      ckm0 = recur0(k,ckm1,ckm2);
      L1 = L1+ckm0; dL1 = dL1+k*ckm0/z; ddL1 = ddL1+k*(k-1)*ckm0/z^2;
      ckm2 = ckm1; ckm1 = ckm0; 
    end
       
    sN = (ckm1*z*(q-gamma*(epsilon+a*delta-a-1)) - ...
         ckm2*z^2*(alpha*beta-(gamma+1)*(epsilon+delta-2)))/(a*(1-gamma));
    
    L2 = 0; dL2 = 0; ddL2 = 0; dm1 = dL2; dm2 = NaN;
    ckm1 = 0; ckm2 = ckm0; 

    L3 = sN; skm2 = 0; skm1 = sN;
    dL3 = N*sN/z; ddL3 = N*(N-1)*sN/z^2; 
    dsm1 = dL3; dsm2 = NaN; skm0 = NaN;
    
    k = N+1;

    while ( (k<=Heun_klimit) && ( (dsm2~=dsm1) || (abs(skm0)>eps) || ...
      (dm2~=dm1) || (abs(ckm0)>eps) ) )

      skm0 = recur0(k,skm1,skm2);
      ckm0 = recur1(k,ckm1,ckm2,skm0,skm1,skm2);

      L2 = L2+ckm0; dL2 = dm1+k*ckm0/z; ddL2 = ddL2+k*(k-1)*ckm0/z^2;
      ckm2 = ckm1; ckm1 = ckm0;
      dm2 = dm1; dm1 = dL2;

      L3 = L3+skm0; dL3 = dsm1+k*skm0/z; ddL3 = ddL3+k*(k-1)*skm0/z^2;
      skm2 = skm1; skm1 = skm0;
      dsm2 = dsm1; dsm1 = dL3;

      k = k+1;

    end
    
    numb = k-1;

    val = L1 + L2 + log(z) * L3;
    dval = dL1 + dL2 + log(z) * dL3 + L3/z;
    ddval = ddL1 + ddL2 - L3/z^2 + 2*dL3/z + log(z) * ddL3;

    if ( isinf(val) || isinf(dval) || isnan(val) || isnan(dval) )

      wrnmsg = 'HeunL00log: failed convergence of recurrence and summation; '; 
      val = NaN; dval = NaN; err = NaN;
      
    else
    
      val2 = ( z*(z-1)*(z-a)*ddval+(gamma*(z-1)*(z-a)+delta*z*(z-a)+epsilon*z*(z-1))*dval ) / ...
              (q-alpha*beta*z);
      err1 = abs(val-val2);

      if abs(q-alpha*beta*z)<0.01
        err2 = abs(L1)*eps*N + abs(ckm0)*sqrt(numb-N+1) + abs(L2)*eps*(numb-N+1) + ...
               abs(log(z)) * ( abs(skm0)*sqrt(numb-N+1) + abs(L3)*eps*(numb-N+1) );
        err =  min(err1,err2);
      else
        err = err1;
      end
       
    end

  end

end

