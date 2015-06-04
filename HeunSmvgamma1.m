% the second local Heun function at z=0
% the case when gamma = 1
%
% evaluation by analytic continuation from z=0
% to z using consequence of power expansions
%
% Usage:
% [val,dval,err,numb,wrnmsg] = HeunS0gamma1(a,q,alpha,beta,delta,path2z)
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
function [val,dval,err,numb,wrnmsg] = HeunSmvgamma1(a,q,alpha,beta,delta,path2z)

  global Heun_cont_coef;
  
  if isempty(Heun_cont_coef)
    HeunOpts();
  end

  if (path2z(1)~=0)

    wrnmsg = 'HeunSmvgamma1: path2z should start from zero; '; 
    val = NaN; dval = NaN; err = NaN; numb = NaN;    

  elseif (dist2sing(path2z,[1,a])<10*eps)||(dist2sing(path2z(2:end),[0])<10*eps)

    wrnmsg = 'HeunSmvgamma1: path2z is too close to one of the singular points; '; 
    val = NaN; dval = NaN; err = NaN; numb = NaN;

  else

    R0 = min([1,abs(a)]);
    z = path2z(end);
  
    if abs(z) <= R0*(Heun_cont_coef+0.01)
  
      [val,dval,err,numb,wrnmsg] = HeunS00gamma1(a,q,alpha,beta,delta,z);
    
    else

      errsum = 0; numbsum = 0;

      failure = false;

      for k=1:(length(path2z)-1)
        z0 = path2z(k);
        theta = angle(path2z(k+1)-z0);
        insearch = true;
    
        while (insearch && ~failure)

          if (z0==0)
            R = R0;
          else
            R = min(abs([z0,z0-1,z0-a]));
          end

          if abs(path2z(k+1)-z0) <= R*Heun_cont_coef
            z1 = path2z(k+1); insearch = false;
          else
            z1 = z0 + Heun_cont_coef * R * exp(1i*theta);
          end

          if (z0==0)
            
              [H0,dH0,err,numb,wrnmsg] = HeunS00gamma1(a,q,alpha,beta,delta,z1);
            
            if ~isempty(wrnmsg)
              wrnmsg = strcat('HeunSmvgamma1: problem invoking HeunS00gamma1(',num2str(a),',',...
                num2str(q),',',num2str(alpha),',',num2str(beta),',',...
                num2str(delta),',',num2str(z1),...
                '); warning: ',wrnmsg,'; ');
              failure = true;
            end

          else
          
            [H0,dH0,err,numb,wrnmsg] = HeunGfromZ0(a,q,alpha,beta,1,delta,z1,z0,H0,dH0);
            if ~isempty(wrnmsg)
              wrnmsg = strcat('HeunSmvgamma1: problem invoking HeunGfromZ0(',num2str(a),',',...
                num2str(q),',',num2str(alpha),',',num2str(beta),',',...
                num2str(gamma),',',num2str(delta),',',num2str(z1),',',...
                num2str(z0),',',num2str(H0),',',num2str(dH0),...
                '); warning: ',wrnmsg,'; ');
              failure = true;
            end
           
          end
          
          errsum = errsum + err;
          numbsum = numbsum + numb;

          z0 = z1;

        end
        
        if failure
          break;
        end
        
      end
      
      val = H0; dval = dH0; numb = numbsum; err = errsum;
      
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
