## 3-species competition in a self-cycling fermenter
#=
The code contained in this Julia file constructs Figure 2 from 
Smith?, Meadows, Wolkowicz (2023) Competition in the Nutrient Driven Self Cycling Fermentation Process
=#
using DifferentialEquations, Parameters
using Plots, LaTeXStrings; Plots.backend(:gr)
using ColorBrewer
theme(:wong,
     size = (450,200),
     fontfamily = "computer modern")
pars = (m = [2.142653, 1.0, 7.0], K = [6.33, 1.0, 32.5], Sin = 20, sbar = 0.1, r = 0.5)

f(s,p,j) = p.m[j]*s/(p.K[j]+s)
sbar_plus = pars.Sin*pars.r+(1-pars.r)*pars.sbar
## Plot response functions
plot(0:0.01:20.0,[x->f(x,pars,1),x-> f(x,pars,2), x-> f(x,pars,3)],
                lw=2,
                label = [L"f_1(s)" L"f_2(s)" L"f_3(s)"],
                xlabel = L"s")
plot!([pars.sbar,pars.sbar],[0.0,2.6],
    lw = 2,
    linestyle = :dash,
    c = 4,
    label = nothing)
plot!([sbar_plus,sbar_plus],[0.0,2.6],lw =2 ,linestyle = :dash, c = 4,label = nothing)
annotate!([0.5],[2.5],L"\bar{s}")
annotate!([11.3],[2.5],L"\bar{s}^+")
savefig("Figures/response_functions.pdf")
#Continuous Phase
function continuous!(du,u,p,t)
    s = u[1]
    n = length(u)-1
    du[1] = ds = -sum([f(s,pars,j)*u[j+1] for j in 1:n])
    for j in 1:n
        du[j+1] = f(s,p,j)*u[j+1]
    end
end

#Decanting condition
condition(u,t,integrator) =  u[1] -integrator.p.sbar

#Decanting phase
function decant!(integrator)
    @unpack Sin,r,sbar = integrator.p
    integrator.u[1] = Sin*r +(1.0-r)*sbar
    integrator.u[2:end] .= (1.0-r).*integrator.u[2:end]
end

cb = ContinuousCallback(condition,decant!)

##Timespan and initial conditions

tspan = (0.0,40.0)
ICs = [sbar_plus,2.5,1.0,5.5]
prob = ODEProblem(continuous!,ICs,tspan,pars)

sol = solve(prob,Tsit5(), saveat=0.01, callback = cb)

## Plot time series
p1 = plot(sol, idxs = [2,3,4], 
    labels = [L"x_1(t)" L"x_2(t)" L"x_3(t)"], 
    lw = 1,
    legend = :topright,
    xlabel = L"t")
savefig("Figures/Three_Species_all.pdf")

## Calculate Floquet Multipliers

function Λ(pars,j,k)
    @unpack sbar, r,Sin,m,K = pars
    sbarp = Sin*r+(1.0-r)*sbar
    return ((1.0/(1.0-r))^(m[k]*(K[j]+Sin)/(m[j]*(K[k]+Sin))-1))*((K[k]+sbarp)/(K[k]+sbar))^(m[k]*(K[j]-K[k])/(m[j]*(K[k]+Sin)))
end

Λ(pars,1,2)
Λ(pars,2,1)
Λ(pars,2,3)
Λ(pars,3,2)
Λ(pars,1,3)
Λ(pars,3,1)

## Pairwise Plots

# x1, x2
ICs12 = [sbar_plus,6.0,1.0,0.0]
prob12 = remake(prob, u0 = ICs12,tspan = (0.0,100.0))
sol = solve(prob12,Tsit5(),saveat=0.01,callback = cb)
p2 = plot(sol, idxs = [2 3],
    lw = 1, 
    c = [1 2], 
    label = [L"x_1(t)" L"x_2(t)"],
    legend = false,
    xlabel = "",
    )
savefig("Figures/Three_species_x1x2.pdf")
# x1, x3
ICs13 = [sbar_plus,7.0,0.0,4.0]
prob13 =remake(prob, u0 = ICs13)
sol = solve(prob13,Tsit5(), saveat=0.01,callback = cb)
p3 = plot(sol, idxs = [2 4],
    lw = 1,
    c= [1 3],
    label = [L"x_1(t)" L"x_3(t)"],
    legend = false,
    xlabel = "")
savefig("Figures/Three_species_x1x3.pdf")
# x2, x3
ICs23 = [sbar_plus,0.0,1.5,7.0]
prob23 = remake(prob,u0 = ICs23)
sol = solve(prob23,Tsit5(),saveat=0.01,callback = cb)
p4 = plot(sol, idxs = [3 4],
    lw = 1, 
    c = [2 3], 
    label = [L"x_2(t)" L"x_3(t)"],
    xlabel = "",
    legend = false)
savefig("Figures/Three_Species_x2x3.pdf")

l = @layout [a{0.6w} [grid(3,1)]]
plot(p1,p2,p3,p4,layout = l,
    size = (450,350),
    title = ["A" "B" "C" "D"],
    titlefont = font(10,"computer modern"),
    titleloc = :left,
    top_margin = (0.5,:mm))
savefig("Figures/comp4.pdf")