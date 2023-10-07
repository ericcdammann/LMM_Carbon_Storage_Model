function wave_power(amp,wind,fetch,hs)

    depth = hs
    fac = min(1,depth/(2*amp))
    D = (depth+(depth-fac*2*amp))/2
    Hs = wave_height(fetch,wind,D)
    Tp = wave_period(fetch,wind,D)
    kk = wave_number(1/Tp,D)
    cg = 2*pi/kk/Tp*0.5*(1+2*kk*D/(sinh(2*kk*D)))
    W = cg*9800/16*abs(Hs)^2

    return W
    
end