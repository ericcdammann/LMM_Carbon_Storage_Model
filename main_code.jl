using Plots
using Statistics

include("wave_number.jl")
include("wave_height.jl")
include("wave_period.jl")
include("wave_power.jl")
include("lagoon_bed_shear_stress.jl")
include("marsh_bed_shear_stress.jl")

# Computational Parameters

years = 150
n_iter = 1000
t_init = 0
t_fin = years*365*24*60*60
h = (t_fin-t_init)/n_iter
t_sec = range(t_init,t_fin,n_iter)
t = t_sec / (60*60*24*365)

# Inital Conditions

wm_init = 1000
wl_init = 9000
hm_init = 0.2602
hl_init = 3.2786
wb_init = wl_init+wm_init

ctot_init = 0
clab_init = 1.222
cref_init = 1.222
cdec_init = 0

# Model Paramters 

rhos = 1000
P = 12.5*3600*1
ws = 0.5*10^-3
tcr = 0.1
Co = 0.05
wind = 6
ka = 2
ke = 0.1/(365*24*3600)
amp = 1.4/2
RSLR = [1 3 5 7]*(10^-3)/(3600*24*365)
rhom = 1000
lamda = 0.0001
beta = 10^10
kk = 0.012/(24*60*60)
wtidal = 15*wm_init

# Preallocating Arrays

WL = zeros(n_iter,length(RSLR))
HM = zeros(n_iter,length(RSLR))
HL = zeros(n_iter,length(RSLR))
WB = zeros(n_iter,length(RSLR))

CTOT = zeros(n_iter,length(RSLR))
CLAB = zeros(n_iter,length(RSLR))
CREF = zeros(n_iter,length(RSLR))

MFLUX = zeros(n_iter,length(RSLR))
MTOT = zeros(n_iter,length(RSLR))

for i in eachindex(RSLR)

    WL[1,i] = wl_init
    HM[1,i] = hm_init
    HL[1,i] = hl_init
    WB[1,i] = wb_init

    CTOT[1,i] = ctot_init
    CLAB[1,i] = clab_init
    CREF[1,i] = cref_init

end

for j in eachindex(RSLR)

    for k in 1:n_iter-1

        fetch = WL[k,j]
        hl = HL[k,j]
        hm = HM[k,j]
        wb = WB[k,j]
        wm = wb-fetch

        ctot = CTOT[k,j]
        clab = CLAB[k,j]
        cref = CREF[k,j]
        mtot = MTOT[k,j]

        # Organic Sediment Production

        BMax = 2.500
        Dmin = 0
        Dmax = 0.7167*2*amp-0.483
        AA = 0.25*(-Dmin-Dmax)*(Dmax-3*Dmin)
        Bpeak = BMax*(hm-Dmax)*(hm-Dmin)/AA
    
        if Bpeak <= 1e-3
        
            Bpeak = 0
    
        end
    
        Bfrac = (Bpeak/BMax)
        nuGp = 0.0138
        AMC = (365/2)*Bpeak*(nuGp)/(365*24*60*60)
        por = 1000/2650
        chiref = 0.15
        Rref = AMC*chiref
        po = 1000 
        O = (1/por)*(Rref/po)
    
        # Average Depths
    
        fac = min(1,hl/(2*amp))
        fac2 = min(1,hm/(2*amp))
        Hf = (hl+(hl-fac*2*amp))/2
        Hm = (hm+(hm-fac2*2*amp))/2
    
        # Sediment Concentrations
    
            # Sediment Concentration in Lagoon
    
        tau_1 =  lagoon_bed_shear_stress(fetch,wind,Hf)
        tau_lagoon = max((tau_1-tcr)/tcr,0)*lamda
        Cr = rhos*tau_lagoon/(1+tau_lagoon)
    
            # Sediment Concentration in Marsh
    
        if Hm > 1e-4
    
            tau_2 =  marsh_bed_shear_stress(fetch,wind,Hm,Bfrac)
    
        else 
    
            tau_2 = 0
    
        end
    
        tau_marsh = max((tau_2-tcr)/tcr,0)*lamda
        Cm = rhos*tau_marsh/(1+tau_marsh)
    
        # Sediment Fluxes
    
            # Flux from Lagoon to Marsh
    
        Fm = (Cr-Cm)*min(2*amp,hm)/P/rhom
    
            # Flux from Open Ocean to Lagoon
    
        Fc = (Cr-Co)*(fac*2*amp)/P/rhom
    
        # Lagoon/Marsh Edge
    
            # Marsh Erosion Rate
    
        dist = 10
        hs = hm+(hl-hm)*(1-exp(-dist*0.1/hl));
        W = wave_power(amp,wind,fetch,hs)
        Be = ke*W/(hs-hm)
    
            # Marsh Progradation Rate
    
        Ba = ka*Cr*ws/rhom
    
        # Upland Migration
    
        Umig = (Fm+O)/beta

        # Carbon Decomposition

        chilab = 1-chiref
        Rlab  = AMC*chilab
        phi = wm/(2*wtidal)
        
        # Carbon Dynamics Outputs

        dCTOTdt = AMC-(kk*clab)
        dCLABdt = Rlab-(kk*clab)
        dCREFdt = Rref
        mflux = phi*kk*clab

        # Morphodynamic Outputs
    
        dWLdt = Be-Ba
        dHMdt = -Fm-O+RSLR[j]
        dHLdt = -(Be-Ba)*(hl-hm)/fetch+Fm*wm/fetch+Fc+RSLR[j]
        dWBdt = Umig

        # Foward Euler Scheme

        WL[k+1,j] = fetch+h*dWLdt
        HM[k+1,j] = hm+h*dHMdt
        HL[k+1,j] = hl+h*dHLdt
        WB[k+1,j] = wb+h*dWBdt

        CTOT[k+1,j] = ctot+h*dCTOTdt
        CLAB[k+1,j] = clab+h*dCLABdt
        CREF[k+1,j] = cref+h*dCREFdt
        MTOT[k+1,j] = mtot+h*mflux*wm

        MFLUX[k,j] = mflux 
        MFLUX[k+1,j] = mflux

    end

end

mass_balance = CTOT-(CLAB+CREF).+(clab_init+cref_init)

# Morphodynamics

WM = WB - WL 
p1_md = plot(t,WL)
p2_md = plot(t,WM)
p3_md = plot(t,HL)
p4_md = plot(t,HM)

fig1=plot(p1_md,p2_md,p3_md,p4_md, layout = 4, 
title=["Width of Lagoon (wl)" "Width of Marsh (wl)" "Depth of Lagoon (hl)" "Depth of Marsh below MHW (hm)" ] ,
label=["RSLR = 1 mm/year" "RSLR = 3 mm/year" "RSLR = 5 mm/year" "RSLR = 7 mm/year"],
xlabel="Years",
ylabel="Meters",
size=(1000,1000)
)
display(fig1)

p1_cd = plot(t,CREF)
p2_cd = plot(t,CLAB)

fig2 = plot(p1_cd,p2_cd, layout=2,
title=["Refractory Carbon Pool (Cref)" "Labile Carbon Pool (Clab)"],
label=["RSLR = 1 mm/year" "RSLR = 3 mm/year" "RSLR = 5 mm/year" "RSLR = 7 mm/year"],
xlabel="Years",
ylabel="kgs of Carbon per m^2",
legend=:outerbottom,
size=(700,500)
)
display(fig2)

MFLUX = MFLUX*(1.0*10^6*3600)
leg_ent = ["RSLR = $(RSLR[l]/((10^-3)/(3600*24*365))) mm/yr" for l in eachindex(RSLR)]
fig3=plot(title="CH4 Flux over Time", legend=:outerright)
for m in eachindex(RSLR)
    plot!(fig3,t,MFLUX[:,m],label=leg_ent[m])
end
xlabel!("Years")
ylabel!("mgs of CH4 per m^2 per hour")
up_bound = t*0 .+27.4
low_bound = t*0 .+8
plot!(fig3,t,up_bound,linecolor=:black,linestyle=:dash,label="Comer-Warner et al. 2022")
plot!(fig3,t,low_bound,linecolor=:black,linestyle=:dash,label=false)
display(fig3)

fig4=plot(t,MTOT,
title="Total Methane Emissions (Mtot)",
label=["RSLR = 1 mm/year" "RSLR = 3 mm/year" "RSLR = 5 mm/year" "RSLR = 7 mm/year"]
)
xlabel!("Years")
ylabel!("kgs of CH4")
display(fig4)

average_mflux = zeros(1,length(RSLR))
for n in eachindex(RSLR)
    average_mflux[n] = mean(MFLUX[:,n])
end
print(average_mflux)