function wave_period(fetch,wind,h)

    g = 9.8;
    delta = h*g/wind^2
    chi = fetch*g./wind^2
    ni = 0.133*(tanh(0.331*delta^1.01)*tanh(5.215*10^-4*chi^0.73/tanh(0.331*delta^1.01)))^-0.37
    Tp = wind/ni/g

    return Tp

end