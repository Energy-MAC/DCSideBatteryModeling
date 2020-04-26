using PowerSystems
using NLsolve

omib_file_dir= joinpath(pwd(), "data/OMIB.raw")
# Data in raw file only contains partial network information. Checks are disabled since we know that the data can't pass all checks.
omib_sys = System(PowerModelsData(omib_file_dir), runchecks=false)
slack_bus = [c for c in get_components(Bus, omib_sys) if c.bustype == BusTypes.REF][1]

battery = GenericBattery(
    name = "Battery",
    primemover = PrimeMovers.BA,
    available = true,
    bus = slack_bus,
    energy = 5.0,
    capacity = (min = 5.0, max = 100.0),
    rating = 70,
    activepower = 10.0,
    inputactivepowerlimits = (min = 0.0, max = 50.0),
    outputactivepowerlimits = (min = 0.0, max = 50.0),
    reactivepower = 0.0,
    reactivepowerlimits = (min = -50.0, max = 50.0),
    efficiency = (in = 0.80, out = 0.90),
)
add_component!(omib_sys, battery)

res = solve_powerflow!(omib_sys, nlsolve)

###### Converter Data ######
converter() = AverageConverter(
    690.0, #Rated Voltage
    2.75,  #Rated MVA
)
###### DC Source Data ######
dc_source() = FixedDCSource(600.0) #Needs to be improved for this paper.

###### Filter Data ######
filter() = LCLFilter(
    0.08, #Series inductance lf in pu
    0.003, #Series resitance rf in pu
    0.074, #Shunt capacitance cf in pu
    0.2, #Series ractance rg to grid connection
    0.01, #Series resistance lg to grid connection
)

###### Outer Control ######

# Needs to be implemented

######## Inner Control ######
inner_control() = CurrentControl(
    0.59, #kpv:: Voltage controller proportional gain
    736.0, #kiv:: Voltage controller integral gain
    0.0, #kffv:: Binary variable enabling the voltage feed-forward in output of current controllers
    0.0, #rv:: Virtual resistance in pu
    0.2, #lv: Virtual inductance in pu
    1.27, #kpc:: Current controller proportional gain
    14.3, #kiv:: Current controller integral gain
    0.0, #kffi:: Binary variable enabling the current feed-forward in output of current controllers
    50.0, #ωad:: Active damping low pass filter cut-off frequency
    0.2, #kad:: Active damping gain
)

inverter = PSY.DynamicInverter(
        1, #Number
        "Storage", #name
        slack_bus, #bus
        1.0, # ω_ref,
        get_voltage(bus), #V_ref
        get_activepower(battery), #P_ref
        get_reactivepower(battery), #Q_ref
        2.75, #MVABase
        converter(), #converter
        outer_control(), #outer control
        inner_control(), #inner control voltage source
        dc_source(), #dc source
        no_pll(), #pll
        filter(),
    ) #filter
