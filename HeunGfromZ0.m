% Heun function, by power series at Z0 for the given values H(Z0)=H0, H'(Z0)=dH0 
%
% Usage:
% [val,dval,ddval,err,numb,wrnmsg] = HeunGfromZ0(a,q,alpha,beta,gamma,delta,z,Z0,H0,dH0)
%
% it is assumed that z, Z0 are not equal to 0, 1, a
%   and |z-Z0| < min{|Z0|,|Z0-1|,|Z0-a|}
%
% The function uses a parameter Heun_klimit which can be changed by HeunOpts: 
% Heun_klimit is the maximum number of series' terms
%
% Returned parameters:
% val is the value of the Heun function at point z
% dval is the value of z-derivative of the Heun function at point z
% err is the estimated error
% numb is the number of power series terms needed for the evaluation
% wrnmsg is empty if computations are ok
%   otherwise it is a diagnostic message
%
% if wrnmsg is not empty, then the function returns val, dval = NaN
%
% Oleg V. Motygin, copyright 2015, license: GNU GPL v3
%
% 04 June 2015
%
function [val,dval,err,numb,wrnmsg] = HeunGfromZ0(a,q,alpha,beta,gamma,delta,z,Z0,H0,dH0)
  
  global Heun_klimit;
  
  if isempty(Heun_klimit)
    HeunOpts();
  end

  R = min([abs(Z0),abs(Z0-1),abs(Z0-a)]);
  
  wrnmsg = '';

  epsilon = alpha+beta+1-gamma-delta;
  
  if (abs(z-Z0)>=R)

    wrnmsg = strcat('HeunGfromZ0: z is out of the convergence radius = ',num2str(R)); val = NaN; dval = NaN; err = NaN; numb = NaN;

  elseif ((min(abs(z-[1 a]))<eps) || (min(abs(Z0-[1 a]))<eps))

    wrnmsg = 'HeunGfromZ0: z or Z0 is too close to the singular points'; 
    val = NaN; dval = NaN; err = NaN; numb = NaN;

  elseif (z==Z0)

    val = H0; dval = dH0; 
    err = 0; numb = 0;

  else
  
    zeta = z-Z0;
  
    recur = @(k,ckm1,ckm2,ckm3) ...
     - (ckm1*zeta*(1-1/k)*((gamma+delta+epsilon+3*(k-2))*Z0^2 - ...
         ((a+1)*(gamma+2*k-4)+epsilon+a*delta)*Z0 + a*(gamma+k-2)) + ...
        ckm2*zeta^2*((2*(1-2/k)*(gamma+delta+epsilon+3/2*(k-3))+alpha*beta/k)*Z0 - ...
          q/k - (1-2/k)*((a+1)*(gamma+k-3)+epsilon+a*delta)) + ...
        ckm3*zeta^3*((1-3/k)*(gamma+epsilon+delta+k-4)+alpha*beta/k) ) / ...
        ((k-1)*Z0*(Z0-1)*(Z0-a));

    ckm3 = H0; ckm2 = dH0*zeta; ckm1 = recur(2,ckm2,ckm3,0);

    val = ckm3+ckm2+ckm1; vm1 = val; vm2 = NaN;
    dm2 = dH0; dm1 = dH0+2*ckm1/zeta; dval = dm1;
    ddval = 2*ckm1/zeta^2; 
  
    k = 3; ckm0 = 1;
    
    while (k<=Heun_klimit) && ( ( vm2~=vm1 ) || ( dm2~=dm1 ) || (abs(ckm0)>eps) )
      ckm0 = recur(k,ckm1,ckm2,ckm3);
      val = val+ckm0; dval = dm1+k*ckm0/zeta;
      ddval = ddval+k*(k-1)*ckm0/zeta^2;
      ckm3 = ckm2; ckm2 = ckm1; ckm1 = ckm0;
      vm2 = vm1; vm1 = val;
      dm2 = dm1; dm1 = dval;
      k=k+1;
    end

    numb = k-1;

    if ( isinf(val) || isinf(dval) || isnan(val) || isnan(dval) )
    
      wrnmsg = 'failed convergence of recurrence and summation'; 
      val = NaN; dval = NaN; err = NaN; 
    
    else
    
      val2 = ( z*(z-1)*(z-a)*ddval+(gamma*(z-1)*(z-a)+delta*z*(z-a)+epsilon*z*(z-1))*dval ) / (q-alpha*beta*z);
      err1 = abs(val-val2);
      
      if abs(q-alpha*beta*z)<0.01
        err2 = abs(ckm0) * sqrt(numb) + abs(val) * eps * numb;
        err =  min(err1,err2);
      else
        err = err1;
      end

    end

  end  
end
