% the second local Heun function at z=0
%
% evaluation by analytic continuation from z=0
% to z using a consequence of power expansions
%
% Usage:
% [val,dval,err,numb,wrnmsg] = HeunS0(a,q,alpha,beta,gamma,delta,z)
%
% it is assumed that z is not equal to 1 or a;
%
% The function uses a parameter Heun_cont_coef which can be changed by HeunOpts: 
% for each power expansion of the analytic continuation procedure Heun_cont_coef
% is the relative distance from the centre to the calculated point 
%
% Returned parameters:
% val is the value of the Heun function
% dval is the value of z-derivative of the Heun function
% err is the estimated error
% numb is the total number of power series terms needed for the evaluation
% wrnmsg is empty if computations are ok
%   otherwise it is a diagnostic message
%
% if wrnmsg is not empty, then the function returns val, dval = NaN
%
% Oleg V. Motygin, copyright 2015, license: GNU GPL v3
%
% 01 May 2015
%
function [val,dval,err,numb,wrnmsg] = HeunS0(a,q,alpha,beta,gamma,delta,z)

  global Heun_cont_coef;
  
  if isempty(Heun_cont_coef)
    HeunOpts();
  end

  if (real(z/a)>=1)&&(imag(z/a)==0)||(real(z)>=1)&&(imag(z)==0)
    
    wrnmsg = 'HeunS0: z belongs to a possible branch cut; '; 
    val = NaN; dval = NaN; err = NaN; numb = NaN;
    
  elseif (min(abs(z-[1 a]))<eps)

    wrnmsg = 'HeunS0: z is too close to one of the singular points; '; 
    val = NaN; dval = NaN; err = NaN; numb = NaN;

  else

    if (gamma==1)
    
      [val,dval,err,numb,wrnmsg] = HeunS0gamma1(a,q,alpha,beta,delta,z);
    
    else
    
      epsilon = alpha+beta+1-gamma-delta;
      [H0w,dH0w,errw,numb,wrnmsg] = HeunL0(a,q-(gamma-1)*(epsilon+a*delta),beta-gamma+1,alpha-gamma+1,2-gamma,delta,z);
      val = z^(1-gamma)*H0w;
      dval = (1-gamma)*z^(-gamma)*H0w + z^(1-gamma)*dH0w;
      err = abs(z^(1-gamma))*errw;
      
    end
  
  end
  
end

