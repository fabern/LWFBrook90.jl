# test script to set up Pfynwald model
include("LWFBrook90.jl")

input_prefix="pfynwald";
input_path="./examples/PFY2025-capacitance/";

model = LWFBrook90.loadSPAC(input_path, input_prefix; simulate_isotopes = true);

simulation = LWFBrook90.setup(model; requested_tspan = (0,100));

LWFBrook90.simulate!(simulation)