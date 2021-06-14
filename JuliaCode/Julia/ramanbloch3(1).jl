### A Pluto.jl notebook ###
# v0.12.21

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 2be91432-6181-11eb-1715-19c0b8f1a75a
using Plots, LinearAlgebra

# ╔═╡ 5540104e-6181-11eb-14cd-0bbbe92838a8
# define Bloch sphere
function bloch()
	N=100
	a=0.075
	b=100
	lw=1
	theta=LinRange(0.0,2*pi,N)
	phi=LinRange(0.0,2*pi,N)
	z0=0.0
	plot3d(sin.(theta).*cos.(0.0),sin.(theta).*sin.(0.0),cos.(theta),legend=false,line = (:black, lw, abs.(a.-sin.(phi)/b)),size=(1.0*600,1.0*520),axis=([], false),box=:off,camera = (40, 30))
	for n in 2:17
		long=n*pi/18
		plot3d!(sin.(theta).*cos.(long),sin.(theta).*sin.(long),cos.(theta),line = (:black, lw, abs.(a.-sin.(phi)/b)))
	end
#	plot3d!(sin.(theta),z0*ones(N),cos.(theta),line = (:black, lw, abs.(a.-z0/b)))
#	plot3d!(z0*ones(N),sin.(theta),cos.(theta),line = (:black, lw, abs.(a.-sin.(theta)/b)))
	plot3d!(sin.(theta),cos.(theta),zeros(N),line = (:black, lw, abs.(a.-cos.(theta)/b)))
	z1=0.66
	rho1=sqrt(1-z1*z1)
	x1=rho1*sin.(theta)
	y1=rho1*cos.(theta)
	plot3d!(x1,y1,z1*ones(N),line = (:black, lw, abs.(a.-cos.(theta)/b)))
	z2=-0.66
	rho2=sqrt(1-z2*z2)
	x2=rho2*sin.(theta)
	y2=rho2*cos.(theta)
	plot3d!(x2,y2,z2*ones(N),line = (:black, lw, abs.(a.-cos.(theta)/b)))
end

# ╔═╡ 870e18a2-6181-11eb-2515-cfc29e9e4c6b
@bind Delta html"<input type=range min=0.01 max=20.0 step=0.05>"

# ╔═╡ 883bebbc-6181-11eb-3ce4-3f1f5c5e177d
@bind Omega1 html"<input type=range min=0.01 max=5.0 step=0.05>"

# ╔═╡ 8d743c92-6181-11eb-029a-37f9b6fb9c8f
@bind Omega2 html"<input type=range min=0.01 max=5.0 step=0.05>"

# ╔═╡ 91aa876c-6181-11eb-1dc7-c7936a8ea7ab
#@bind phiL1 html"<input type=range min=0.00 max=20.0 step=0.05>"

# ╔═╡ 96a6251e-6181-11eb-0ad9-6b4b56c1ece3
#@bind phiL2 html"<input type=range min=0.00 max=20.0 step=0.05>"

# ╔═╡ 9c1dec66-6181-11eb-2a5f-131548a29058
@bind tt html"<input type=range min=0.00 max=100.0 step=0.005>"

# ╔═╡ 64f75592-6181-11eb-14e2-63fbc7265f80
begin

	phiL1 = 0.0
	phiL2 = 0.0
	delta = (abs(Omega2)^2-abs(Omega1)^2)/(4*Delta)
	#define Hamiltonian
	H = 0.5*[-delta 0 Omega1*exp(1im*phiL1);0 delta Omega2*exp(1im*phiL2);conj(Omega1)*exp(-1im*phiL1) conj(Omega2)*exp(-1im*phiL2) 2*Delta]
	evals, evecs = eigen(H)

	#define time parameter
	num = 1000
	tmax = tt
	t = LinRange(0,tmax,num)
	dt = t[2]-t[1]

	#define initial state
	rho0 = [1.0;0.0;0.0]
	state1 = [1.0;0.0;0.0]
	state2 = [0.0;1.0;0.0]
	state3 = [0.0;0.0;1.0]

	#evolution of quantum state
	s1 = zeros(Complex,num)
	s2 = zeros(Complex,num)
	s3 = zeros(Complex,num)
	theta1 = zeros(num)
	phi1 = zeros(num)
	theta2 = zeros(num)
	phi2 = zeros(num)
	rho_all = rho0
	for i in 1:num
		global rho_all = evecs*Diagonal(exp.(-1im*evals*dt))*conj(evecs')*rho_all
		s1[i] = abs(state1'*rho_all)^2
		s2[i] = abs(state2'*rho_all)^2
		s3[i] = abs(state3'*rho_all)^2
		theta1[i] = 2*atan(sqrt(abs(rho_all[2]^2)+abs(rho_all[3]^2))/abs(rho_all[1]))
		phi1[i] = atan(imag(rho_all[2])/real(rho_all[2])) + atan(imag(rho_all[1])/real(rho_all[1]))
		theta2[i] = atan(abs(rho_all[3]/rho_all[2]))
		phi2[i] = atan(imag(rho_all[3]/rho_all[2])/real(rho_all[3]/rho_all[2]))
	end

	u1 = sin.(theta1).*cos.(phi1)
	v1 = sin.(theta1).*sin.(phi1)
	w1 = cos.(theta1)
	u2 = sin.(theta2).*cos.(phi2)
	v2 = sin.(theta2).*sin.(phi2)
	w2 = cos.(theta2)

	l = @layout [a;b]
	p1 = plot(t,real.(s1),label="Prob1")
	p1 = plot!(t,real.(s2),label="Prob2")
	p1 = plot!(t,real.(s3),label="Prob3")
	p2 = bloch()
	p2 = scatter3d!(u1,v1,w1,ms=4,zcolor=t/tmax,m =(:heat, 0.25),markerstrokewidth=0)
	p2 = scatter3d!(u2,v2,w2,ms=4,zcolor=t/tmax,m =(:blues, 0.25),markerstrokewidth=0)
	p3 = plot(p1,p2,layout=l,size=(450,700))
	#savefig("ramantran")


end



# ╔═╡ 9e8def0c-618b-11eb-35c8-6171de7267c1
[1 2 3].*[2 3 4]

# ╔═╡ Cell order:
# ╠═2be91432-6181-11eb-1715-19c0b8f1a75a
# ╠═5540104e-6181-11eb-14cd-0bbbe92838a8
# ╠═870e18a2-6181-11eb-2515-cfc29e9e4c6b
# ╠═883bebbc-6181-11eb-3ce4-3f1f5c5e177d
# ╠═8d743c92-6181-11eb-029a-37f9b6fb9c8f
# ╠═91aa876c-6181-11eb-1dc7-c7936a8ea7ab
# ╠═96a6251e-6181-11eb-0ad9-6b4b56c1ece3
# ╠═9c1dec66-6181-11eb-2a5f-131548a29058
# ╠═64f75592-6181-11eb-14e2-63fbc7265f80
# ╠═9e8def0c-618b-11eb-35c8-6171de7267c1
