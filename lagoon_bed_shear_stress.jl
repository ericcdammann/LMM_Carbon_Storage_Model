function lagoon_bed_shear_stress(fetch,wind,Hl)

    ko = 0.001;
    Hs = wave_height(fetch,wind,Hl)
    Tp = wave_period(fetch,wind,Hl)
    k = wave_number(1/Tp,Hl)
    Um = (pi*Hs/Tp/sinh(k*Hl))
    aw = Tp*Um/(2*pi)
    fw = 0.4*(aw/ko)^-0.75
    tau_1 = 1/2*1020*fw*Um^2

    return tau_1

end
