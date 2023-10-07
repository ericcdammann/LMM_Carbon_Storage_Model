function marsh_bed_shear_stress(fetch,wind,Hm,Bfrac)

    if Bfrac==0  
    
        ko = 0.001;
        Hs = wave_height(fetch,wind,Hm)
        Tp = wave_period(fetch,wind,Hm)
        k = wave_number(1/Tp,Hm)
        Um = (pi*Hs/Tp/sinh(k*Hm))
        aw = Tp*Um/(2*pi)
        fw = 0.4*(aw/ko)^-0.75
        tau_2 = 1/2*1020*fw*Um^2
    
    else
    
        tau_2 = 0;
    
    end

    return tau_2

end