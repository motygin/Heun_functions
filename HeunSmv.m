% the second local Heun function at z=0
%
% evaluation by analytic continuation from z=0
% to z using a consequence of power expansions
%
% Usage:
% [val,dval,err,numb,wrnmsg] = HeunSmv(a,q,alpha,beta,gamma,delta,path2z)
%
% path2z is a list of points from 0 to z
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
% 04 June 2015
%
function [val,dval,err,numb,wrnmsg] = HeunSmv(a,q,alpha,beta,gamma,delta,path2z)

  global Heun_cont_coef;
  
  if isempty(Heun_cont_coef)
    HeunOpts();
  end

  if (path2z(1)~=0)

    wrnmsg = 'HeunSmv: path2z should start from zero; '; 
    val = NaN; dval = NaN; err = NaN; numb = NaN;    

  elseif (dist2sing(path2z,[1,a])<10*eps)||(dist2sing(path2z(2:end),[0])<10*eps)

    wrnmsg = 'HeunSmv: path2z is too close to one of the singular points; '; 
    val = NaN; dval = NaN; err = NaN; numb = NaN;

  else

    if (gamma==1)
    
      [val,dval,err,numb,wrnmsg] = HeunSmvgamma1(a,q,alpha,beta,delta,path2z);
    
    else
    
      z = path2z(end);

      epsilon = alpha+beta+1-gamma-delta;
      [H0w,dH0w,errw,numb,wrnmsg] = HeunLmv(a,q-(gamma-1)*(epsilon+a*delta),beta-gamma+1,alpha-gamma+1,2-gamma,delta,path2z);

      zp0 = path2z(2)^(-gamma);
% increment of argument on the line segment (z0,z1)
      addarg = @(z0,z1) asin((real(z0)*imag(z1)-real(z1)*imag(z0))/(abs(z0)*abs(z1)));

      for k=3:length(path2z)

        zp1 = zp0 * exp(-gamma*(log(abs(path2z(k)))-log(abs(path2z(k-1)))+1i*addarg(path2z(k-1),path2z(k))));
        zp0 = zp1;

      end

      val = z * zp0 * H0w;
      dval = (1-gamma)*zp0*H0w + z*zp0*dH0w;
      err = abs(z*zp0)*errw;
      
    end
  
  end
  
end

% minimal distance from path to singular points
function dmin = dist2sing(path,singpts)

  dmin = NaN;
  for k=1:length(path)-1
    for n=1:length(singpts)
      d = dist(singpts(n), path(k), path(k+1));
      if ((k==1)&&(n==1))||(d<dmin)
        dmin = d;
      end
    end
  end 
end

% distance from point P to AB
function [d,isinside] = dist(P, A, B)
    
    isinside = true;
    a = abs(P-A);
    b = abs(P-B);
    c = abs(A-B);
    
    if(a^2>=b^2+c^2)
      isinside = false; d = b; return;
    end
    if(b^2>=a^2+c^2)
      isinside = false; d = b; return;
    end
    
    p = (a+b+c)/2;
    s = sqrt((p-a)*(p-b)*(p-c)*p);
    
    d = s*2/c;

end
