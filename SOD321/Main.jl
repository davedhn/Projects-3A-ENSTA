include("distance.jl")
include("Poly_form.jl")
include("Exp_form.jl")

n,d,f,Amin,Nr,R,regions,coords,D = readInstance("instances/aerodrome_70_1.txt")
Prob_w_poly_const(n,d,f,Amin,Nr,R,regions,D)
Prob_w_exp_const(n,d,f,Amin,Nr,R,regions,D)
