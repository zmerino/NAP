function [a99, a90, a75, c50, b75, b90, b99] = getTargets(Ns)
% getTargets calculates position-dependent threshold values about the mean,
% according to a beta distribution with parameters k and (n + 1 - k), where
% k is the position and n is the total number of positions.  These beta
% distributions represent probability per position for sort order
% statistics for a uniform distribution.
%
% [a99, a90, a75, c50, b75, b90, b99] = getTargets(N) Computes thresholds
% for a sample size of N.  Output format is as follows:
%   a{nn} = a contour line representing nn% of points above the mean
%   b{nn} = a contour line representing nn% of points below the mean

%   a99 = contour line ABOVE mean & 99% of scaled residuals reside below it
%   a90 = contour line ABOVE mean & 90% of scaled residuals reside below it
%   a75 = contour line ABOVE mean & 75% of scaled residuals reside below it
%   c50 = CENTER contour line = median. 50% scaled residuals above & below
%   b75 = contour line BELOW mean & 75% of scaled residuals reside above it
%   b90 = contour line BELOW mean & 90% of scaled residuals reside above it
%   b99 = contour line BELOW mean & 99% of scaled residuals reside above it

% --------------------------------- initialize for subsequent calculations
a99 = zeros(1,Ns);          % 99% of all residuals reside below this line
a90 = zeros(1,Ns);          % 90% of all residuals reside below this line
a75 = zeros(1,Ns);          % 75% of all residuals reside below this line
c50 = zeros(1,Ns);    % 50% of all residuals reside below/above this line
b75 = zeros(1,Ns);          % 75% of all residuals reside above this line
b90 = zeros(1,Ns);          % 90% of all residuals reside above this line
b99 = zeros(1,Ns);          % 99% of all residuals reside above this line
%                   NOTE: 99% gives too much numerical error to be useful
   if( Ns > 100000000 )
   error('maximum Ns = 100000000 = 100M. Must extend the limit');
   end
Ns1 = Ns + 1;
% ---------------------------------------------------- set up du - spacing
np = ceil( log10(Ns) );
du1000 = 0.001;
duMIN = du1000*(0.2)^np;                                        % heristic
nipW = 1000*(np+1);       % # of integration points within WINDOW: heristic
%       du1000 = duMIN*(gRate^nipW) => nipW*log(gRate) = log(du1000/duMIN)
gRate = exp( log(du1000/duMIN)/nipW );               % gRate = growth rate
duRelative0 = zeros(1,nipW);
du = duMIN;
   for i=1:nipW
   duRelative0(i) = du;
   du = du*gRate;
   end
nipMAX = 1000 + nipW + nipW + 1000;                     % => four sections
kStart = 1000 + nipW;
uList = zeros(1,nipMAX);
duRelative = zeros(1,nipMAX);
duRelative(1:1000) = du1000;
duRelative(1001:kStart) = duRelative0(nipW:-1:1);
duRelative(kStart+1:nipMAX-1000) = duRelative0;
duRelative(kStart+1 + nipW:nipMAX) = du1000;
%disp( duRelative );                                   % quick debug/check
bd0 = zeros(1,nipMAX);                         % => beta distribution: CDF
bd1 = zeros(1,nipMAX);                  % => beta distribution: CDF' = PDF
bd2 = zeros(1,nipMAX);                % => beta distribution: CDF'' = PDF'
% ------------------------- select the sort order statistic (SOS) k-values
nHalf = ceil(Ns/2);                            % employ refection symmetry
nksos = nHalf;              % # of k-values for SOS: True when nHalf <= 50
   if( nHalf <= 50 )     % note: 100*(0.025) = 2.5 => floor(100*0.025) = 2
   ksosArray = 1:nHalf;     % k-th SOS index, SOS => sort order statistics
   else                       % this part of the code requires nHalf >= 50 
   nksos = 0;
   k075 = floor(0.075*Ns);
   k = 0;
   dk = 1;
      while(k < k075)
         for i16times=1:16
         k = k + dk;
            if( k > k075 )
            break;
            else
            nksos = nksos + 1;
            end
         end
      dk = 2*dk;
      end
   nksos = nksos + 9;
   ksosArray = zeros(1,nksos);
   % --------------------------------------------------- populate k-values
   jk = 0;
   k = 0;
   dk = 1;
      while(k < k075)
         for i16times=1:16
         k = k + dk;
            if( k > k075 )
            break;
            else
            jk = jk + 1;
            ksosArray(jk) = k;
            end
         end
      dk = 2*dk;
      end
   jk1 = jk + 1;
   jkLast = nksos - 1;
   ff = 0.10;                                          % initialize at 10%
      for jk=jk1:jkLast
      ksosArray(jk) = floor(ff*Ns); 
      ff = ff + 0.05;                                    % increment by 5%
      end
   ksosArray(nksos) = nHalf;   % finish at 50% since nHalf = ceil(0.50*Ns)
   end
%disp( ksosArray );   
% error('stop here for now');
tag = -ones(1,nHalf);
   for jjj=1:nksos
   ksos = ksosArray(jjj);
   tag(ksos) = 1;
   alpha = ksos;
   beta = Ns - ksos + 1;
   alpha1 = alpha - 1;
   beta1 = beta - 1;
   %alpha2 = alpha1 - 1;                              % <= implicitly used
   %beta2 = beta1 - 1;                                % <= implicitly used
%  PDF(u) = u^(alpha-1)*(1 - u)^(beta - 1);                  % unormalized
%  PDF(u) = u^alpha1*(1 - u)^beta1;                          % unormalized
   uMode = (alpha - 1)/(alpha + beta - 2);
   uMode = max(1.0e-10,uMode);
% -------------------------------------------- generate integration points
   k = kStart;
   uList(k) = uMode;
      while( uList(k) < 1 )
      k = k + 1;
      uList(k) = uList(k-1) + duRelative(k); 
      end
   kMAX = k;
   uList(kMAX) = 1;
   % --------------
   k = kStart + 1;
      while( uList(k) > 0 )
      k = k - 1;
      uList(k-1) = uList(k) - duRelative(k); 
      end
   kMIN = k;
   uList(kMIN) = 1.0e-25;
   u = uList(kMIN:kMAX);
   sMode = kStart - kMIN + 1;
   %disp( [u(sMode),uMode] );
   %disp( sMode );
   %pause
   nip = kMAX - kMIN;
% ========================================= determine range of integration
   sL = sMode;
      if( Ns < 100 )
      dsL = 1;
      else
      dsL = 20;
      end
      while( sL > 1 )
      sL = sL - dsL;
         if( sL < 1 )
         sL = 1;
         break;
         end
      rA = u(sL)/uMode;
      rB = ( 1 - u(sL) )/( 1 - uMode );
      log10TestRatio = alpha1*log10(rA) + beta1*log10(rB);
         if( log10TestRatio < -20 )
         break;
         end
         if( log10TestRatio > -15 )
         dsL = 1;
         end
      end
%    disp(['-------------------------------------------------------', ... 
%          '----------------------- ksos = ',num2str(ksos)]);
   sR = sMode;
      if( Ns < 100 )
      dsR = 1;
      else
      dsR = 20;
      end
      while( sR < nip )
      sR = sR + dsR;
         if( sR > nip )
         sR = nip;
         break;
         end
      rA = u(sR)/uMode;
      rB = ( 1 - u(sR) )/( 1 - uMode );
      log10TestRatio = alpha1*log10(rA) + beta1*log10(rB);
      %disp( num2str( [sR,rA,rB,log10TestRatio,nip] ) );
         if( log10TestRatio < -20 )
         break;
         end
         if( log10TestRatio > -10 )
         dsR = 1;
         end
      end
% ------------------------------------------------------ right propagation
   ln_rA = log( u(sL)/uMode );
   ln_rB = log( ( 1 - u(sL) )/(1 - uMode) );
   ln_bd1 = alpha1*ln_rA + beta1*ln_rB;               % using natural logs
   bd1(sL) = exp(ln_bd1);
   ln1 = ln_bd1 - ln_rA;         % => ln1 = alpha2*log(rA) + beta1*log(rB)
   ln2 = ln_bd1 - ln_rB;         % => ln2 = alpha1*log(rA) + beta2*log(rB)
   bd2(sL) = (alpha1/uMode)*exp(ln1) - ( beta1/(1 - uMode) )*exp(ln2);
   bd0(sL) = 0;
%    disp( [bd0(sL),bd1(sL),bd2(sL),uMode] );
%    pause
% -------------------------------------------------------------- integrate
   sLess1 = sL - 1;
   sL1 = sL + 1;
      for s=sL1:sMode 
      sLess1 = sLess1 + 1;
      du = u(s) - u(sLess1);
      duHalf = du/2;
      ln_rA = log( u(s)/uMode );
      ln_rB = log( ( 1 - u(s) )/(1 - uMode) );
      ln_bd1 = alpha1*ln_rA + beta1*ln_rB;            % using natural logs
      bd1(s) = exp(ln_bd1);
      bd2(s) = -bd2(sLess1)  + ( bd1(s) - bd1(sLess1) )/duHalf;
      bd0(s) =  bd0(sLess1)  + ( bd1(s) + bd1(sLess1) )*duHalf ...
                             + ( bd2(sLess1) - bd2(s) )*du*du/8;
      end 
   %bd0(sL:sR)
   topValue = bd0(sMode);
   %pause;
% ------------------------------------------------------ left propagation
   ln_rA = log( u(sR)/uMode );
   ln_rB = log( ( 1 - u(sR) )/(1 - uMode) );
   ln_bd1 = alpha1*ln_rA + beta1*ln_rB;               % using natural logs
   bd1(sR) = exp(ln_bd1);
   ln1 = ln_bd1 - ln_rA;         % => ln1 = alpha2*log(rA) + beta1*log(rB)
   ln2 = ln_bd1 - ln_rB;         % => ln2 = alpha1*log(rA) + beta2*log(rB)
   bd2(sR) = (alpha1/uMode)*exp(ln1) - ( beta1/(1 - uMode) )*exp(ln2);
   bd0(sR) = 0;
%    disp( [bd0(sR),bd1(sR),bd2(sR),uMode] );
%    pause
% -------------------------------------------------------------- integrate
   sLess1 = sR;
   sMode1 = sMode + 1;
      for s=sR:-1:sMode1
      sLess1 = sLess1 - 1;
      du = u(s) - u(sLess1);
      duHalf = du/2;
      ln_rA = log( u(sLess1)/uMode );
      ln_rB = log( ( 1 - u(sLess1) )/(1 - uMode) );
      ln_bd1 = alpha1*ln_rA + beta1*ln_rB;            % using natural logs
      bd1(sLess1) = exp(ln_bd1);
      bd2(sLess1) = -bd2(s)  + ( bd1(s) - bd1(sLess1) )/duHalf;
      bd0(sLess1) =  bd0(s)  - ( bd1(s) + bd1(sLess1) )*duHalf ...
                             + ( bd2(s) - bd2(sLess1) )*du*du/8;
      end
% ----------------------------------------------------- end of integration

% -------------------------------- match left and right integrals at sMode
   shift_bd0 = topValue - bd0(sMode);
   bd0(sMode:sR) = bd0(sMode:sR) + shift_bd0;
   %disp(shift_bd0)
   %disp( bd0(sL:sR) )
% ---------------------------------------------------------- normalize CDF
   scaleFactor = 1/bd0(sR);
   bd0 = bd0*scaleFactor;   


% ------------------------------------------------ calculate contour lines
   F99L = 0.01;
   F90L = 0.10;
   F75L = 0.25;
   F50B = 0.50;   % B => both lower (L) and upper (U)
   F75U = 0.75;
   F90U = 0.90;
   F99U = 0.99;
% -----------------
   max_dF99L = 100;
   max_dF90L = 100;
   max_dF75L = 100;
   max_dF50B = 100;   % B => both lower (L) and upper (U)
   max_dF75U = 100;
   max_dF90U = 100;
   max_dF99U = 100;
   s99L = sL;
   s90L = sL;
   s75L = sL;
   s50B = sL;
   s75U = sR;
   s90U = sR;
   s99U = sR;
% ------------------------------------------------------------------------
      for s=sL:sR
      F = bd0(s);
      dF99L = abs(F - F99L);
         if( dF99L < max_dF99L )
         max_dF99L = dF99L;
         s99L = s;
         end
      dF90L = abs(F - F90L);
         if( dF90L < max_dF90L )
         max_dF90L = dF90L;
         s90L = s;
         end
      dF75L = abs(F - F75L);
         if( dF75L < max_dF75L )
         max_dF75L = dF75L;
         s75L = s;
         end
      dF50B = abs(F - F50B);
         if( dF50B < max_dF50B )
         max_dF50B = dF50B;
         s50B = s;
         end
      dF75U = abs(F - F75U);
         if( dF75U < max_dF75U )
         max_dF75U = dF75U;
         s75U = s;
         end
      dF90U = abs(F - F90U);
         if( dF90U < max_dF90U )
         max_dF90U = dF90U;
         s90U = s;
         end
      dF99U = abs(F - F99U);
         if( dF99U < max_dF99U )
         max_dF99U = dF99U;
         s99U = s;
         end
      end
% --------------------------- define residuals to be w.r.t. mean quantiles
   mu = (1:Ns)/Ns1;                                % => mean SOS positions
   z = ( u - mu(ksos) )*sqrt(Ns+2);                   % => scaled residual
% ------------------------------------------------ quadratic interpolation
   x = 0.01;
   s = s99L;
   x0 = bd0(s-1);
   x1 = bd0(s);
   x2 = bd0(s+1);
   y0 = z(s-1);
   y1 = z(s);
   y2 = z(s+1);
   c0 = (x - x1)*(x - x2)/( (x0 - x1)*(x0 - x2) );
   c1 = (x - x0)*(x - x2)/( (x1 - x0)*(x1 - x2) );
   c2 = (x - x0)*(x - x1)/( (x2 - x0)*(x2 - x1) );
   b99(ksos) = y0*c0 + y1*c1 + y2*c2;
% ------------------------------------------------------------------------
   x = 0.10;
   s = s90L;
   x0 = bd0(s-1);
   x1 = bd0(s);
   x2 = bd0(s+1);
   y0 = z(s-1);
   y1 = z(s);
   y2 = z(s+1);
   c0 = (x - x1)*(x - x2)/( (x0 - x1)*(x0 - x2) );
   c1 = (x - x0)*(x - x2)/( (x1 - x0)*(x1 - x2) );
   c2 = (x - x0)*(x - x1)/( (x2 - x0)*(x2 - x1) );
   b90(ksos) = y0*c0 + y1*c1 + y2*c2;
% ------------------------------------------------------------------------
   x = 0.25;
   s = s75L;
   x0 = bd0(s-1);
   x1 = bd0(s);
   x2 = bd0(s+1);
   y0 = z(s-1);
   y1 = z(s);
   y2 = z(s+1);
   c0 = (x - x1)*(x - x2)/( (x0 - x1)*(x0 - x2) );
   c1 = (x - x0)*(x - x2)/( (x1 - x0)*(x1 - x2) );
   c2 = (x - x0)*(x - x1)/( (x2 - x0)*(x2 - x1) );
   b75(ksos) = y0*c0 + y1*c1 + y2*c2; 
% ------------------------------------------------------------------------
   x = 0.50;
   s = s50B;
   x0 = bd0(s-1);
   x1 = bd0(s);
   x2 = bd0(s+1);
   y0 = z(s-1);
   y1 = z(s);
   y2 = z(s+1);
   c0 = (x - x1)*(x - x2)/( (x0 - x1)*(x0 - x2) );
   c1 = (x - x0)*(x - x2)/( (x1 - x0)*(x1 - x2) );
   c2 = (x - x0)*(x - x1)/( (x2 - x0)*(x2 - x1) );
   c50(ksos) = y0*c0 + y1*c1 + y2*c2;
% ------------------------------------------------------------------------
   x = 0.75;
   s = s75U;
   x0 = bd0(s-1);
   x1 = bd0(s);
   x2 = bd0(s+1);
   y0 = z(s-1);
   y1 = z(s);
   y2 = z(s+1);
   c0 = (x - x1)*(x - x2)/( (x0 - x1)*(x0 - x2) );
   c1 = (x - x0)*(x - x2)/( (x1 - x0)*(x1 - x2) );
   c2 = (x - x0)*(x - x1)/( (x2 - x0)*(x2 - x1) );
   a75(ksos) = y0*c0 + y1*c1 + y2*c2;
% ------------------------------------------------------------------------
   x = 0.90;
   s = s90U;
   x0 = bd0(s-1);
   x1 = bd0(s);
   x2 = bd0(s+1);
   y0 = z(s-1);
   y1 = z(s);
   y2 = z(s+1);
   c0 = (x - x1)*(x - x2)/( (x0 - x1)*(x0 - x2) );
   c1 = (x - x0)*(x - x2)/( (x1 - x0)*(x1 - x2) );
   c2 = (x - x0)*(x - x1)/( (x2 - x0)*(x2 - x1) );
   a90(ksos) = y0*c0 + y1*c1 + y2*c2;
% ------------------------------------------------------------------------
   x = 0.99;
   s = s99U;
   x0 = bd0(s-1);
   x1 = bd0(s);
   x2 = bd0(s+1);
   y0 = z(s-1);
   y1 = z(s);
   y2 = z(s+1);
   c0 = (x - x1)*(x - x2)/( (x0 - x1)*(x0 - x2) );
   c1 = (x - x0)*(x - x2)/( (x1 - x0)*(x1 - x2) );
   c2 = (x - x0)*(x - x1)/( (x2 - x0)*(x2 - x1) );
   a99(ksos) = y0*c0 + y1*c1 + y2*c2;
   end 
% ------------------------------------------ interpolate lemon drop points
j2 = 3;
   for k=1:nHalf
      if( tag(k) < 0 )
      x = mu(k);
         while( ksosArray(j2) < k )
         j2 = j2 + 1;
         j1 = j2 - 1;
         j0 = j1 - 1;
         k0 = ksosArray(j0);
         k1 = ksosArray(j1);
         k2 = ksosArray(j2);
         dk21 = k2 - k1;
            while( (k1 - k0) < floor(0.92*dk21) )
            j0 = j0 - 1;
            k0 = ksosArray(j0);
            end
% ------------------------------------------------- generate anchor points
         x0 = mu(k0);
         x1 = mu(k1);
         x2 = mu(k2);
         %-----------------
         y0a99 = a99(k0);
         y0a90 = a90(k0);
         y0a75 = a75(k0);
         y0c50 = c50(k0);
         y0b75 = b75(k0);
         y0b90 = b90(k0);
         y0b99 = b99(k0);
         %-----------------
         y1a99 = a99(k1);
         y1a90 = a90(k1);
         y1a75 = a75(k1);
         y1c50 = c50(k1);
         y1b75 = b75(k1);
         y1b90 = b90(k1);
         y1b99 = b99(k1);
         %-----------------
         y2a99 = a99(k2);
         y2a90 = a90(k2);
         y2a75 = a75(k2);
         y2c50 = c50(k2);
         y2b75 = b75(k2);
         y2b90 = b90(k2);
         y2b99 = b99(k2);
         end
      c0 = (x - x1)*(x - x2)/( (x0 - x1)*(x0 - x2) );
      c1 = (x - x0)*(x - x2)/( (x1 - x0)*(x1 - x2) );
      c2 = (x - x0)*(x - x1)/( (x2 - x0)*(x2 - x1) );
      % ---------------------------------------------
      b99(k) = c0*y0b99 + c1*y1b99 + c2*y2b99;
      b90(k) = c0*y0b90 + c1*y1b90 + c2*y2b90;
      b75(k) = c0*y0b75 + c1*y1b75 + c2*y2b75;
      c50(k) = c0*y0c50 + c1*y1c50 + c2*y2c50;
      a75(k) = c0*y0a75 + c1*y1a75 + c2*y2a75;
      a90(k) = c0*y0a90 + c1*y1a90 + c2*y2a90;
      a99(k) = c0*y0a99 + c1*y1a99 + c2*y2a99;
      end
   end
% --------------------------------------------------- copy reflected parts
nHalf1 = nHalf + 1;
   for k=nHalf1:Ns
   %disp( num2str([aLine(1,Ns1-k), k, Ns1 - k]) );
   a99(k) = -b99(Ns1-k);
   a90(k) = -b90(Ns1-k);
   a75(k) = -b75(Ns1-k);
   c50(k) = -c50(Ns1-k);
   b75(k) = -a75(Ns1-k);
   b90(k) = -a90(Ns1-k);
   b99(k) = -a99(Ns1-k);
   end
end
