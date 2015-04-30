% Change internal parameters of HeunL, HeunL00, HeunL00log, HeunGfromZ0
%
% Usage: HeunOpts('cont_coef',Heun_cont_coef,'klimit',Heun_klimit,
%                 'proxco',Heun_proxco,'proxcoinf',Heun_proxcoinf)
%
% parameters are optional
%
% if a pair of parameters is omitted, but the corresponding global variable
% is empty then it is set to the default value
%
% Use [] as the value to reset coefficient to its default
% e.g. HeunOpts('cont_coef',[])
%
% Heun_cont_coef is used in HeunL0 and HeunS0gamma1; for each power expansion
% of the analytic continuation procedure Heun_cont_coef is the relative
% (to the radius of convergence) distance from the centre to the calculated point.
% By default Heun_cont_coef = 0.38 (golden ratio)
%
% Heun_klimit is used in HeunL00, HeunL00log, HeunGfromZ0, HeunS0gamma1, 
% HeunS00gamma1; it is the maximum number of power series' terms.
% Default value is 1000 
%
% Heun_proxco and Heun_proxcoinf are used in HeunLS; they specify relative proximity
% to singular point where special representation is used
% By default Heun_proxco = 0.05; Heun_proxcoinf = 5 
%
% Oleg V. Motygin, copyright 2015, license: GNU GPL v3
%
% 14 April 2015
%
function HeunOpts(varargin)

  global Heun_cont_coef Heun_klimit Heun_proxco Heun_proxcoinf;
  
  [reg, props] = parseparams(varargin);
  opts = cell2struct(props(2:2:end),props(1:2:end),2);
  
  if isfield(opts,'cont_coef')
    if isempty(opts.cont_coef)
      Heun_cont_coef = 0.38;
    else
      Heun_cont_coef = opts.cont_coef;
    end
  else
    if isempty(Heun_cont_coef)
      Heun_cont_coef = 0.38;
    end
  end

  if isfield(opts,'klimit')
    if isempty(opts.klimit)
      Heun_klimit = 1000;
    else
      Heun_klimit = opts.klimit;
    end
  else
    if isempty(Heun_klimit)
      Heun_klimit = 1000;
    end
  end
  
  if isfield(opts,'proxco')
    if isempty(opts.proxco)
      Heun_proxco = 0.05;
    else
      Heun_proxco = opts.proxco;
    end
  else
    if isempty(Heun_proxco)
      Heun_proxco = 0.05;
    end
  end
  
  if isfield(opts,'proxcoinf')
    if isempty(opts.proxcoinf)
      Heun_proxcoinf = 5;
    else
      Heun_proxcoinf = opts.proxcoinf;
    end
  else
    if isempty(Heun_proxcoinf)
      Heun_proxcoinf = 5;
    end
  end
  
end
