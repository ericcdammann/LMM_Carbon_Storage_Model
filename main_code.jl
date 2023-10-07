using Plots

include("wave_number.jl")
include("wave_height.jl")
include("wave_period.jl")
include("wave_power.jl")
include("lagoon_bed_shear_stress.jl")
include("marsh_bed_shear_stress.jl")

# Computational Parameters

years = 150
n_iter = 100000
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

Ctot_init = 0
ctot_init = 0
Clab_init = 0
clab_init = 0
Cref_init = 0
cref_init = 0
torg_init = 0
Mtot_init = 0
mtot_init = 0
closs_init = 0

# Model Paramters 

rhos = 1000
P = 12.5*3600*1
ws = 0.5*10^-3
tcr = 0.1
Co = 0.1
wind = 6
ka = 2
ke = 0.1/(365*24*3600)
amp = 1.4/2
RSLR = 5*(10^-3)/(3600*24*365)
rhom = 1000
lamda = 0.0001
beta = 10^-3
k = 1.3*10^-6
wtidal = 25000

# Preallocating Arrays

WL = zeros(n_iter)
WL[1] = wl_init
HM = zeros(n_iter)
HM[1] = hm_init
HL = zeros(n_iter)
HL[1] = hl_init
WB = zeros(n_iter)
WB[1] = wb_init

CTOT = zeros(n_iter)
CTOT[1] = Ctot_init
cTOT = zeros(n_iter)
cTOT[1] = ctot_init
CLAB = zeros(n_iter)
CLAB[1] = Clab_init
cLAB = zeros(n_iter)
cLAB[1] = clab_init
CREF = zeros(n_iter)
CREF[1] = Cref_init
cREF = zeros(n_iter)
cREF[1] = cref_init
TORG = zeros(n_iter)
TORG[1] = torg_init
MTOT = zeros(n_iter)
MTOT[1] = Mtot_init
mTOT = zeros(n_iter)
mTOT[1] = mtot_init
cLOSS = zeros(n_iter)
cLOSS[1] = closs_init

for i in 1:n_iter-1

    fetch = WL[i]
    hl = HL[i]
    hm = HM[i]
    wb = WB[i]
    wm = wb-fetch

    Ctot = CTOT[i]
    ctot = cTOT[i]
    Clab = CLAB[i]
    clab = cLAB[i]
    Cref = CREF[i]
    cref = cREF[i]
    torg = TORG[i]
    Mtot = MTOT[i]
    mtot = mTOT[i]
    closs = cLOSS[i]

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
    AMC = (180)*Bpeak*(nuGp)/(365*24*3600)
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

    # Carbon Budget

    Rae = ((Ba-Be)+Umig)*po*por*torg
    chilab = 1-chiref
    Rlab  = AMC*chilab
    phi = wm/(2*wtidal)
    
    # Defining Equations

        # Morphodynamic Outputs
    
    dWLdt = Be-Ba
    dHMdt = -Fm-O+RSLR
    dHLdt = -(Be-Ba)*(hl-hm)/fetch+Fm*wm/fetch+Fc+RSLR
    dWBdt = Umig

        # Carbon Dynamics Outputs

    dCTOTdt = AMC*wm
    dcTOTdt = AMC
    dCLABdt = Rlab*wm-k*Clab - Rae
    dcLABdt = Rlab-k*clab
    dCREFdt = Rref*wm-Rae
    dcREFdt = Rref
    dTORGdt = O
    dMTOTdt = k*Clab*phi
    dmTOTdt = k*clab*phi
    dcLOSSdt = k*clab

    # Foward Euler Scheme

    WL[i+1] = fetch+h*dWLdt
    HM[i+1] = hm+h*dHMdt
    HL[i+1] = hl+h*dHLdt
    WB[i+1] = wb+h*dWBdt

    CTOT[i+1] = Ctot+h*dCTOTdt
    cTOT[i+1] = ctot+h*dcTOTdt
    CLAB[i+1] = Clab+h*dCLABdt
    cLAB[i+1] = clab+h*dcLABdt
    CREF[i+1] = Cref+h*dCREFdt
    cREF[i+1] = cref+h*dcREFdt
    TORG[i+1] = torg+h*dTORGdt
    MTOT[i+1] = Mtot+h*dMTOTdt
    mTOT[i+1] = mtot+h*dmTOTdt
    cLOSS[i+1] = closs+h*dcLOSSdt

end

# Width of the Marsh

WM = WB-WL

# Mass Balance

mass_balance = cTOT-(cLAB+cREF+cLOSS)

# Plotting Results

# Morphodynamic Outputs [L]

p1_md = plot(t,WL)
p2_md = plot(t,WM)
p3_md = plot(t,HL)
p4_md = plot(t,HM)

display(plot(p1_md,p2_md,p3_md,p4_md, layout = 4, 
plot_title="Morphodynamic Outputs", 
label=["wl" "wm" "hl" "hm"],
xlabel="Years",
ylabel="Meters"
))

# Carbon Dynamics Outputs [M * L^-1]

p1_cd1 = plot(t,CTOT)
p2_cd1 = plot(t,CLAB)
p3_cd1 = plot(t,CREF)
p4_cd1 = plot(t,MTOT)

display(plot(p1_cd1,p2_cd1,p3_cd1,p4_cd1, layout = 4,
plot_title="Carbon Dynamics Outputs [M*L^-1]", 
label=["Ctot" "Clab" "Cref" "Mtot"],
xlabel="Years",
ylabel="kg per meter"
))

# Carbon Dynamics Outputs [M * L^-2]

p1_cd2 = plot(t,cTOT)
p2_cd2 = plot(t,cLAB)
p3_cd2 = plot(t,cREF)
p4_cd2 = plot(t,mTOT)

display(plot(p1_cd2,p2_cd2,p3_cd2,p4_cd2, layout = 4,
plot_title="Carbon Dynamics Outputs [M*L^-2]", 
label=["ctot" "clab" "cref" "mtot"],
xlabel="Years",
ylabel="kg per meter^2"
))

# Mass Balance

display(plot(t,mass_balance,
plot_title="Mass Balance", 
label=false,
xlabel="Years",
ylabel="Error"
))