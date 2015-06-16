% (local) Heun function, based on HeunL0 (analytic continuation from z=0)
% but with improvements near singular points z = 1, a, infinity 
%
% Usage:
% [val,dval,err,numb,wrnmsg] = HeunL(a,q,alpha,beta,gamma,delta,z,memlimit)
%
% the optional parameter memlimit (500 as default) is the maximum number of already
%   computed matching data which are kept in memory
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
% 15 June 2015
%
function [val,dval,err,numb,wrnmsg] = HeunL(a,q,alpha,beta,gamma,delta,z,varargin)

  isnonpositivegamma = abs(ceil(gamma-5*eps)+abs(gamma))<5*eps;

  if isnonpositivegamma
    gamma = ceil(gamma-5*eps);
  end

  if (real(z/a)>=1)&&(angle(z/a)==0)||(real(z)>=1)&&(angle(z)==0)||...
    (isnonpositivegamma&&(angle(z)==0)&&(sign(real(z))<0))
    
    wrnmsg = 'HeunL: z belongs to a possible branch cut; '; 
    val = NaN; dval = NaN; err = NaN; numb = NaN;
    
  elseif abs(angle(z))==pi
    
    [val,dval,err,numb,wrnmsg] = HeunL0(a,q,alpha,beta,gamma,delta,z);
    
  else

    if length(varargin)>0
      [val,dval,err,numb,wrnmsg] = HeunLS(1,a,q,alpha,beta,gamma,delta,z,varargin);
    else
      [val,dval,err,numb,wrnmsg] = HeunLS(1,a,q,alpha,beta,gamma,delta,z);
    end

  end

end
