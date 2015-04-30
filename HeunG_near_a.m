% Heun function for z close to a (|z-a|<max{|a|,|a-1|})
%
% Usage:
% [val,dval,ddval,err,numb,wrnmsg] = HeunG_near_a(a,q,alpha,beta,gamma,delta,z,C1,C2)
%
% Returned parameters:
% val is the combination with weights C1 and C2 of two linearly independent solutions
%   of the Heun equation
% dval is the value of z-derivative of the expression
% err is the estimated error
% numb is the total number of series terms needed for the evaluation
% wrnmsg is empty if computations are ok
%   otherwise it is a diagnostic message from HeunL0 or HeunS0
%     and the function returns val, dval = NaN
% 
% Oleg V. Motygin, copyright 2015, license: GNU GPL v3
%
% 28 April 2015
%
function [val,dval,err,numb,wrnmsg] = HeunG_near_a(a,q,alpha,beta,gamma,delta,z,C1,C2)

  epsilon = alpha+beta+1-gamma-delta;

  val = 0; dval = 0; err = 0; numb = 0; wrnmsg = '';
  
  if C1~=0
%[a_+0_+1_+][\infty_+] in Table 2, Maier, 2007, The 192 solutions of the Heun equation
    [H0w,dH0w,errw,numbw,wrnmsg] = HeunL0((a-1)/a,alpha*beta-q/a,alpha,beta,epsilon,...
                              gamma,(a-z)/a);
    val = C1*H0w;
    dval = -C1*dH0w/a;
    err = abs(C1)*errw;
    numb = numbw;
  end
  
  if C2~=0
    [H0w,dH0w,errw,numbw,wrnmsg2] = HeunS0((a-1)/a,alpha*beta-q/a,alpha,beta,epsilon,...
                              gamma,(a-z)/a);
    val = val+C2*H0w;
    dval = dval-C2*dH0w/a;
    err = err+abs(C2)*errw;
    numb = numb+numbw;
    wrnmsg = strcat(wrnmsg,wrnmsg2);
  end
  
end
