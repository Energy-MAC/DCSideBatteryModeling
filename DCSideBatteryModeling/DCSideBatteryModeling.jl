# Cant use import, see https://github.com/SciML/ModelingToolkit.jl/issues/319
using ModelingToolkit
import PowerSystems
import OrdinaryDiffEq

const MTK = ModelingToolkit
const PSY = PowerSystems

include("model.jl")
#include("parameters.jl")
#include("utils")
