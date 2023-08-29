## 2-species competition in a self-cycling fermenter
#=
The code contained in this Julia file constructs Figure 1 from 
Smith?, Meadows, Wolkowicz (2023) Competition in the Nutrient Driven Self Cycling Fermentation Process
=#
using DifferentialEquations, Parameters
using Plots, LaTeXStrings
using ColorBrewer
theme(:wong,
     size = (450,200),
     fontfamily = "computer modern")
pars = (m = [1.75, 1.0], K = [5.0, 1.0], Sin = 20, sbar = 0.1, r = 0.5)

f(s,p,j) = p.m[j]*s/(p.K[j]+s)


#Continuous Phase
function continuous!(du,u,p,t)
    s,x1,x2,x3 = u
    du[1] = ds = -f(s,p,1)*x1-f(s,p,2)*x2
    du[2] = dx1 = f(s,p,1)*x1
    du[3] = dx2 = f(s,p,2)*x2
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
sbar_plus = pars.Sin*pars.r+(1-pars.r)*pars.sbar
tspan = (0.0,30.0)
ICs = [sbar_plus,1.0,1.0,1.0]

prob = ODEProblem(continuous!,ICs,tspan,pars)

sol = solve(prob,Tsit5(), saveat=20, callback = cb)

plot(sol, idxs = [2,3], labels = [L"x_1(t)" L"x_2(t)"], lw = 2)
savefig("Figures/coexistence_2Species.pdf")

## Calculate Floquet Multipliers
using Setfield


function Λ(pars,j,k)
    @unpack sbar, r,Sin,m,K = pars
    sbarp = Sin*r+(1.0-r)*sbar
    return ((1.0/(1.0-r))^(m[k]*(K[j]+Sin)/(m[j]*(K[k]+Sin))-1))*((K[k]+sbarp)/(K[k]+sbar))^(m[k]*(K[j]-K[k])/(m[j]*(K[k]+Sin)))
end

Λ(pars,1,2)
Λ(pars,2,1)
## Plot Coexistence Region
using DataFrames

df = DataFrame(m1=Float64[],
               K1=Float64[],
               color = Int64[])

mrange = 0.01:0.001:2.0
krange = 0.01:0.005:6.0
for m1 in mrange
    @set! pars.m[1] = m1
    for K1 = krange
        @set! pars.K[1] = K1
        if f(sbar_plus,pars,1)-f(sbar_plus,pars,2) >0
            if f(pars.sbar,pars,1)-f(pars.sbar,pars,2) ≥0
                push!(df,[m1,K1,0])
            elseif (Λ(pars,1,2)>1)&(Λ(pars,2,1)>1)
                    push!(df,[m1,K1,1])
            elseif (Λ(pars,1,2)>1)&(Λ(pars,2,1)<1)
                    push!(df,[m1,K1,2])
            elseif (Λ(pars,1,2)<1)&(Λ(pars,2,1)>1)
                    push!(df,[m1,K1,3])
            else
                    push!(df,[m1,K1,5])
            end
        else 
            if f(pars.sbar,pars,1)-f(pars.sbar,pars,2) ≤ 0
                push!(df,[m1,K1,4])
            elseif (Λ(pars,1,2)>1)&(Λ(pars,2,1)>1)
                    push!(df,[m1,K1,1])
            elseif (Λ(pars,1,2)>1)&(Λ(pars,2,1)<1)
                    push!(df,[m1,K1,2])
            elseif (Λ(pars,1,2)<1)&(Λ(pars,2,1)>1)
                    push!(df,[m1,K1,3])
            else
                    push!(df,[m1,K1,5])
            end
        end
    end
end

heatmap(mrange,krange, Matrix(unstack(df,:K1,:m1,:color)[:,2:end]),
        size = (650,400),
        c = cgrad(:seaborn_colorblind6,6,categorical=true,rev = false),
        legend = false,
        )
xlabel!(L"m_1")
ylabel!(L"K_1")
annotate!(0.5,4.0,L"A")
annotate!([0.6,1.45],[0.25,4.5],[L"B",L"B"])
annotate!([0.9,1.7], [0.25,2.8],L"D")
annotate!(1.5,0.3,L"E")
annotate!([1.9,0.75],[5.8,-0.2],L"C")
savefig("Figures/coexistence_regions.pdf")