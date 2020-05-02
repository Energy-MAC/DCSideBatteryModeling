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
    v_rated = 690.0,
    s_rated = 2.75,
)
###### DC Source Data ######
dc_source() = FixedDCSource(voltage = 600.0) #Not in the original data, guessed.

###### Filter Data ######
filter() = LCLFilter(
    lf = 0.08,
    rf = 0.003,
    cf = 0.074,
    lg = 0.2,
    rg = 0.01,
)

###### Outer Control ######
function outer_control()
    #Need to implement proper outer control Data in PSY
    return OuterControl(virtual_inertia(), reactive_droop())
end

######## Inner Control ######
inner_control() = CurrentControl(
    kpv = 0.59,     #Voltage controller proportional gain
    kiv = 736.0,    #Voltage controller integral gain
    kffv = 0.0,     #Binary variable enabling the voltage feed-forward in output of current controllers
    rv = 0.0,       #Virtual resistance in pu
    lv = 0.2,       #Virtual inductance in pu
    kpc = 1.27,     #Current controller proportional gain
    kic = 14.3,     #Current controller integral gain
    kffi = 0.0,     #Binary variable enabling the current feed-forward in output of current controllers
    ωad = 50.0,     #Active damping low pass filter cut-off frequency
    kad = 0.2,
)
inverter = DynamicInverter(
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

to_json(omib_sys, joinpath(pwd(), "data/OMIB_DCBattery.json"))
