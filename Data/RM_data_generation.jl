using DifferentialEquations
using CSV, Random

r = 4
K = 1
a2 = 3
c2 = 1
d2 = 0.1
a3 = 4
d3 = 0.05
c3 = 0.1
sigma = 0.5

p = [r,K,a2,c2,d2,a3,d3,c3,sigma]

function RM_three_chain!(dx,x,p,t)
    r = p[1]
    K = p[2]
    a2 = p[3]
    c2 = p[4]
    d2 = p[5]
    a3 = p[6]
    d3 = p[7]
    c3 = p[8]

    dx[1] = r*x[1]*(1-x[1]/K) - a2*x[1]*x[2]/(1+x[1])
    dx[2] = c2*a2*x[1]*x[2]/(1+x[1]) - d2*x[2] - a3*x[2]*x[3]/(1+x[2])
    dx[3] = c3*a3*x[2]*x[3]/(1+x[2]) - d3*x[3]
end

g(du,u,p,t) = (du.=sigma*u)

Random.seed!(556)

x0 = rand(3)
tspan = (0.0,100.0)
W = WienerProcess(0.0,0.0,0.0)
problem = SDEProblem(RM_three_chain!,g,x0,tspan,p,noise=W)
data = solve(problem)

CSV.write("RM_data.csv",data)