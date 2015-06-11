% local Heun function, based on HeunL0 or HeunS0
% by analytic continuation from z=0
% with improvements near singular points z = 1, a, infinity 
%
% Usage:
% [val,dval,err,numb,wrnmsg] = HeunLS(numfunc,a,q,alpha,beta,gamma,delta,z,memlimit)
%
% if numfunc=1 HeunLS computes the value of the first local solution HeunL
% otherwise HeunLS computes the value of the second local solution HeunS
%
% the optional parameter memlimit (500 as default) is the maximum number of already
%   computed matching data which are kept in memory
%
% The function uses parameters Heun_proxco and Heun_proxcoinf which can be changed
% by HeunOpts: they specify relative proximity to singular point where we use the
% special representations
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
% 10 June 2015
%
function [val,dval,err,numb,wrnmsg] = HeunLS(numfunc,a,q,alpha,beta,gamma,delta,z,varargin)

  if length(varargin)>0
    memlimit = varargin{1};
  else
    memlimit = 500;
  end

  global Heun_proxco1st Heun_proxcoinf1st Heun_proxco Heun_proxcoinf;
  
  if isempty(Heun_proxco1st)||isempty(Heun_proxco)||isempty(Heun_proxcoinf1st)||isempty(Heun_proxcoinf)
    HeunOpts();
  end
  
  Heun0 = @HeunL0;

  if numfunc~=1
    numfunc = 2;
    Heun0 = @HeunS0;
  end

  Ra = min([abs(a),abs(a-1)]);
  R1 = min([1,abs(a-1)]);
  Rinf = max([1,abs(a)]);

  wrnmsg = '';
  failed = true;
  
  if abs(z-a)<Heun_proxco*Ra

    singpt = 'a'; sdir = a*sign(imag(a))/abs(a);

    if (imag(a)==0)&&((a<0)||(a>1))
      if imag(z)>0
        singpt = 'aU'; sdir = 1;
      else
        singpt = 'aL'; sdir = -1;
      end
    end 

    midpoint = a*0.5 + 0.70710678i*sdir;
    impev = abs(z-a)<Heun_proxco1st*Ra;

    [val,dval,err,numb,wrnmsg,failed] = HeunLSnearsing(numfunc,Heun0,@HeunG_near_a,singpt,midpoint,memlimit,impev,a,q,alpha,beta,gamma,delta,z);

  elseif abs(z-1)<Heun_proxco*R1

    singpt = '1'; sdir = -sign(imag(a));

    if (imag(a)==0)&&(a>0)&&(a<1)
      if imag(z)>0
        singpt = '1U'; sdir = 1;
      else
        singpt = '1L'; sdir = -1;
      end
    end 

    midpoint = 0.5 + 0.70710678i*sdir;
    impev = abs(z-1)<Heun_proxco1st*R1;

    [val,dval,err,numb,wrnmsg,failed] = HeunLSnearsing(numfunc,Heun0,@HeunG_near_1,singpt,midpoint,memlimit,impev,a,q,alpha,beta,gamma,delta,z);

  elseif abs(z)>Heun_proxcoinf*Rinf

    ls = sort([-pi,0,pi,angle(a)]);
    idx = sum(angle(z)>ls);
    singpt = strcat('I',num2str(idx));
    midarg = (ls(idx)+ls(idx+1))/2;

    midpoint = Heun_proxcoinf * Rinf * exp(1i*midarg);
    impev = abs(z)>Heun_proxcoinf1st*Rinf;

    [val,dval,err,numb,wrnmsg,failed] = HeunLSnearsing(numfunc,Heun0,@HeunG_near_infty,singpt,midpoint,impev,memlimit,a,q,alpha,beta,gamma,delta,z);

  end

  if failed
    [val,dval,err,numb,wrnmsg] = Heun0(a,q,alpha,beta,gamma,delta,z);
  end

end

function [val,dval,err,numb,wrnmsg,failed] = HeunLSnearsing(numfunc,Heun0,HeunG_nearsing,singpt,midpoint,memlimit,impev,a,q,alpha,beta,gamma,delta,z)

  val = NaN; dval = NaN; err = NaN; numb = NaN; wrnmsg = '';
  
  [A1,A2,errco,consts_known] = extrdatfromsav(numfunc,a,q,alpha,beta,gamma,delta,singpt);
  
  if ~consts_known && impev

    [val1,dval1,err1,numb1,wrnmsg1] = HeunG_nearsing(a,q,alpha,beta,gamma,delta,midpoint,1,0);
    [val2,dval2,err2,numb2,wrnmsg2] = HeunG_nearsing(a,q,alpha,beta,gamma,delta,midpoint,0,1);
    [val0,dval0,err0,numb0,wrnmsg0] = Heun0(a,q,alpha,beta,gamma,delta,midpoint);
    M = [val1 val2; dval1 dval2];
    dcs = M \ [val0; dval0];
    A1 = dcs(1); A2 = dcs(2);

    errco = err0 + err1 + err2;
    numb0 = numb0 + numb1 + numb2;

    if min(svd(M))/max(svd(M))<10^(-6)     
      A1 = NaN; A2 = NaN;
    end
      
    keepdattosav(numfunc,a,q,alpha,beta,gamma,delta,errco,A1,A2,singpt,memlimit);
    wrnmsg = strcat(wrnmsg0,wrnmsg1,wrnmsg2);

  end
  
  failed = isnan(A1)||isnan(A2)||(~consts_known&&~impev);
  
  if ~failed
    
    [val1,dval1,err1,numb1,wrnmsg1] = HeunG_nearsing(a,q,alpha,beta,gamma,delta,z,1,0);
    if ~isempty(wrnmsg1)
       if abs(A1)<errco
          wrnmsg1 = ''; val1 = 0; dval1 = 0; err1 = 0; 
       end
    end
    [val2,dval2,err2,numb2,wrnmsg2] = HeunG_nearsing(a,q,alpha,beta,gamma,delta,z,0,1);
    if ~isempty(wrnmsg2)
       if abs(A2)<errco
          wrnmsg2 = ''; val2 = 0; dval2 = 0; err2 = 0; 
       end
    end
    wrnmsg = strcat(wrnmsg,wrnmsg1,wrnmsg2);

    val = A1*val1 + A2*val2; dval = A1*dval1 + A2*dval2;
    err = abs(A1)*err1 + abs(A2)*err2 + errco;
    numb = numb1 + numb2;    

    if ~consts_known
      numb = numb + numb0;
    end
  end
  
  failed = failed||isnan(val);
  
end
  
function [A1,A2,errco,consts_known] = extrdatfromsav(numfunc,a,q,alpha,beta,gamma,delta,singpt)

  global savdata;

  A1 = NaN; A2 = NaN; errco=NaN; consts_known = false;
  if length(savdata)~=0
    for k=1:length(savdata)
      if (savdata(k).numfunc==numfunc)&&strcmp(savdata(k).singpt,singpt)&&(savdata(k).a==a)&&...
         (savdata(k).q==q)&&(savdata(k).alpha==alpha)&&(savdata(k).beta==beta)&&...
         (savdata(k).gamma==gamma)&&(savdata(k).delta==delta)
        A1=savdata(k).A1;
        A2=savdata(k).A2;
        errco=savdata(k).errco;
        consts_known = true;
        break;
      end
    end
  end
end

function keepdattosav(numfunc,a,q,alpha,beta,gamma,delta,errco,A1,A2,singpt,memlimit)
 
  global savdata;

  if length(savdata)==0
    savdata=struct('numfunc',numfunc,'singpt',singpt,'a',a,'q',q,'alpha',alpha,...
      'beta',beta,'gamma',gamma,'delta',delta,'errco',errco,'A1',A1,'A2',A2);
  else
    if length(savdata)<=memlimit
      savdata(end+1)=struct('numfunc',numfunc,'singpt',singpt,'a',a,'q',q,'alpha',alpha,...
        'beta',beta,'gamma',gamma,'delta',delta,'errco',errco,'A1',A1,'A2',A2);
    else
      savdata(1).numfunc=numfunc; savdata(1).singpt=singpt; savdata(1).a=a;
      savdata(1).q=q; savdata(1).alpha=alpha; savdata(1).beta=beta;
      savdata(1).gamma=gamma; savdata(1).delta=delta;
      savdata(1).errco=errco; savdata(1).A1=A1; savdata(1).A2=A2;
      savdata = shift(savdata,-1);
    end  
  end
end

