% the second local Heun function at z=0
% the case when gamma = 1
%
% evaluation by analytic continuation from z=0
% to z using consequence of power expansions
%
% Usage:
% [val,dval,err,numb,wrnmsg] = HeunS0gamma1(a,q,alpha,beta,delta,z)
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
% 28 April 2015
%
function [val,dval,err,numb,wrnmsg] = HeunS0gamma1(a,q,alpha,beta,delta,z)

  global Heun_cont_coef;
  
  if isempty(Heun_cont_coef)
    HeunOpts();
  end

  if (real(z/a)>=1)&&(imag(z/a)==0)||(real(z)>=1)&&(imag(z)==0)||...
    (imag(z)==0)&&(sign(z)<0)
    
    wrnmsg = 'HeunS0gamma1: z belongs to a possible branch cut; '; 
    val = NaN; dval = NaN; err = NaN; numb = NaN;
    
  elseif (min(abs(z-[1 a]))<eps)

    wrnmsg = 'HeunS0gamma1: z is too close to one of the singular points; '; 
    val = NaN; dval = NaN; err = NaN; numb = NaN;

  else

    R0 = min([1,abs(a)]);
  
    if abs(z) <= R0*(Heun_cont_coef+0.01)
  
      [val,dval,err,numb,wrnmsg] = HeunS00gamma1(a,q,alpha,beta,delta,z);
    
    else

      if (abs(a)<1)
        P1 = a; P2 = 1;
      else
        P1 = 1; P2 = a;
      end

      R1 = min([abs(P1),abs(a-1)]);
      R2 = min([abs(P2),abs(a-1)]);

      addzz = [];
      
      [d,isinside] = dist(P1, 0, z);
      if isinside && (d < R1/2)
        addzz(1) = P1+exp(1i*(pi/2+angle(z)))*min(R1/2,abs(z-P1))*sign(imag(z/P1));
      end
  
      [d,isinside] = dist(P2, 0, z);
      if isinside && (d < R2/2)
        addzz(end+1) = P2+exp(1i*(pi/2+angle(z)))*min(R2/2,abs(z-P2))*sign(imag(z/P2));
      end

      zz = [0,addzz,z];

      errsum = 0; numbsum = 0;

      failure = false;

      for k=1:(length(zz)-1)
        z0 = zz(k);
        theta = angle(zz(k+1)-z0);
        insearch = true;
    
        while (insearch && ~failure)

          if (z0==0)
            R = R0;
          else
            R = min(abs([z0,z0-1,z0-a]));
          end

          if abs(zz(k+1)-z0) <= R*Heun_cont_coef
            z1 = zz(k+1); insearch = false;
          else
            z1 = z0 + Heun_cont_coef * R * exp(1i*theta);
          end

          if (z0==0)
            
              [H0,dH0,err,numb,wrnmsg] = HeunS00gamma1(a,q,alpha,beta,delta,z1);
            
            if ~isempty(wrnmsg)
              wrnmsg = strcat('HeunS0gamma1: problem invoking HeunS00gamma1(',num2str(a),',',...
                num2str(q),',',num2str(alpha),',',num2str(beta),',',...
                num2str(delta),',',num2str(z1),...
                '); warning: ',wrnmsg,'; ');
              failure = true;
            end

          else
          
            [H0,dH0,err,numb,wrnmsg] = HeunGfromZ0(a,q,alpha,beta,1,delta,z1,z0,H0,dH0);
            if ~isempty(wrnmsg)
              wrnmsg = strcat('HeunS0gamma1: problem invoking HeunGfromZ0(',num2str(a),',',...
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
