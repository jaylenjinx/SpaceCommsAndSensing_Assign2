function s = LFMCW(t,t0,phi0,B,T)
% Function LFMCW() generates linear frequency modulated continuous wave signal with bandwidth B and period T.
% Inputs:
% t: vector of time in second;
% t0: initial time (delay);
% phi0: initial phase;
% B: bandwidth in Hz;
% T: period.
% Output:
% s: vector of LFMCW signal
    
    s = zeros(1,length(t));
    
    for i = 1:length(t)
        ti = t(i) - t0;
        while ti > T
            ti = ti - T;
        end
        
        s(i) = exp(1j*pi*B/T*(ti-T/2).^2);
    end
    
    s = s*exp(1j*phi0);
end