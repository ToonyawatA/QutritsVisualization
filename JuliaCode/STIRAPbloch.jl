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

# ╔═╡ 3c2c4258-756a-11eb-0bf3-0ba77da71717
using Plots, LinearAlgebra

# ╔═╡ 3ee11d1a-75af-11eb-345a-434ac61ba9ec
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

# ╔═╡ 45170b54-75af-11eb-2e7c-f159673105e7
#define dissipator operator
function D(sigma)
	d = size(sigma)[1]
	I0 = Matrix(1.0I,d,d)
	kron(conj(sigma),sigma)-0.5*kron(I0,transpose(conj(sigma))*sigma)-0.5*kron(transpose(transpose(conj(sigma))*sigma),I0)
end

# ╔═╡ 913c74e6-75b0-11eb-3f34-7f2a3d384a79
@bind E_0 html"<input type=range min=0.01 max=50.0 step=0.05>"

# ╔═╡ 93007b4c-75b0-11eb-3e28-15c6c7881a6b
@bind tau1 html"<input type=range min=0.01 max=5.0 step=0.05>"

# ╔═╡ 14bd2b56-75b3-11eb-2ccb-f50700357871
@bind t1 html"<input type=range min=0.00 max=50.0 step=0.05>"

# ╔═╡ 1d8cc94e-75b3-11eb-1eee-890a481c1784
@bind t2 html"<input type=range min=0.00 max=50.0 step=0.05>"

# ╔═╡ 9d264e58-75b0-11eb-12fc-dd4ee11ae3f0
@bind tt html"<input type=range min=0.00 max=100.0 step=0.05>"

# ╔═╡ fd4a2d8a-75b1-11eb-3f7e-f7cf4aa20bef
@bind Delta html"<input type=range min=0.01 max=20.0 step=0.05>"

# ╔═╡ 05da0350-75b2-11eb-3ce0-3b2f1e20c493


# ╔═╡ f782f312-75b3-11eb-3d7f-ff524bf48ece
b = [1 2 3]


# ╔═╡ fddffb74-75b3-11eb-17e5-7f641e5114f9
a = [0 0 0;1 1 1]

# ╔═╡ 6fac2970-75b0-11eb-1c89-d1493aa0c889
begin 
	#define sigma_plus and sig_minus
	g1 = [1.0;0.0;0.0]
	g2 = [0.0;1.0;0.0]
	e = [0.0;0.0;1.0]
	sigma1 = g1*e'
	sigma2 = g2*e'
	sigmax = Complex[0 1 0;1 0 0;0 0 0]
	sigmay = Complex[0 -1im 0;1im 0 0;0 0 0]
	sigmaz = Complex[1 0 0;0 -1 0;0 0 0]
	sigma0 = Complex[0 0 0;0 0 0;0 0 1]
	
	#define Rabi pulse (Gaussian bean profile)
	#E_0 = sqrt(0.5*pi)*(1/tau1)
    #tau =2
	#omega = 1
	#alpha = 1
	num = 1001
	tmax = tt
	t = LinRange(0,tmax,num)
	dt = t[2]-t[1]
	Omega1 = abs.(real(E_0*exp.((-1/(2*tau1^2))*(t.-t1).^2)))
	Omega2 = abs.(real(E_0*exp.((-1/(2*tau1^2))*(t.-t2).^2)))
	
	
	#define Detuning
	delta = 0.0
	
	#define Hamiltonian
	H = zeros(Complex,3,3,num)
	H[1,1,:] .= Delta
	H[2,2,:] .= Delta-delta
	H[1,3,:] .= 0.5*Omega1
	H[2,3,:] .= 0.5*Omega2
	H[3,1,:] .= 0.5*conj(Omega1)
	H[3,2,:] .= 0.5*conj(Omega2)
	
	
	#define initial state
	psi_0 = [1.0;0.0;0.0]
	
	#quantum state evolution 
	u1 = zeros(num)
	v1 = zeros(num)
	w1 = zeros(num)
	u2 = zeros(num)
	v2 = zeros(num)
	w2 = zeros(num)
	rho_ee = zeros(num)
	rho_gg1 = zeros(num)
	rho_gg2 = zeros(num)
	psi_all = psi_0
	for i in 1:num
		evals,evecs = eigen(H[:,:,i])
		global psi_all = evecs*Diagonal(exp.(-1im*dt*evals))*conj(evecs')*psi_all
		rho_all = psi_all*psi_all'
		phi = atan(-1*imag(rho_all[3])/real(rho_all[3]))
		u1[i] = real(tr(rho_all*sigmax))
		v1[i] = real(tr(rho_all*sigmay))
		w1[i] = real(tr(rho_all*sigmaz))
		u2[i] = real(tr(sigma0*rho_all)*cos(phi))
		v2[i] = real(tr(sigma0*rho_all)*sin(phi))
		rho_ee[i] = real(rho_all[9])
		rho_gg2[i] = real(rho_all[5])
		rho_gg1[i] = real(rho_all[1])
	end 
	
	
	
	l = @layout [a;b;c]
	p1 = plot(t,rho_ee,label="Prob in e state (Numerical)",legend=:bottomright)
	p1 = plot!(t,rho_gg2,label="Prob in g2 state (Numerical)")
	p1 = plot!(t,rho_gg1,label="Prob in g1 state (Numerical)")
	p2 = plot(t,Omega1,label="Rabi Pulse1")
	p2 = plot!(t,Omega2,label="Rabi Pulse2")
	p3 = bloch()
	p3 = scatter3d!(u1,v1,w1,ms=3,zcolor=t/tmax,m =(:heat, 0.25),markerstrokewidth=0)
	p3 = scatter3d!(u2,v2,w2,ms=3,zcolor=t/tmax,m =(:algae, 0.25),markerstrokewidth=0)
	p4 = plot(p2,p1,p3,layout=l,size=(500,1200))
	
end
	
	
	



# ╔═╡ ae1c0712-75b4-11eb-2e58-570e887323ca
Omega1

# ╔═╡ 53e86464-75b4-11eb-2050-d5560cf22ec3
H

# ╔═╡ 92ce58de-75b4-11eb-368a-8fc4b0c4be44
conj(a)

# ╔═╡ 146569e2-75b4-11eb-340f-5d12067a8e76
a[1,:] .= b

# ╔═╡ 1a6c44e6-75b4-11eb-3a94-6575eea73c03
a

# ╔═╡ Cell order:
# ╠═3c2c4258-756a-11eb-0bf3-0ba77da71717
# ╠═3ee11d1a-75af-11eb-345a-434ac61ba9ec
# ╠═45170b54-75af-11eb-2e7c-f159673105e7
# ╠═913c74e6-75b0-11eb-3f34-7f2a3d384a79
# ╠═93007b4c-75b0-11eb-3e28-15c6c7881a6b
# ╠═14bd2b56-75b3-11eb-2ccb-f50700357871
# ╠═1d8cc94e-75b3-11eb-1eee-890a481c1784
# ╠═9d264e58-75b0-11eb-12fc-dd4ee11ae3f0
# ╠═fd4a2d8a-75b1-11eb-3f7e-f7cf4aa20bef
# ╠═05da0350-75b2-11eb-3ce0-3b2f1e20c493
# ╠═6fac2970-75b0-11eb-1c89-d1493aa0c889
# ╠═ae1c0712-75b4-11eb-2e58-570e887323ca
# ╠═53e86464-75b4-11eb-2050-d5560cf22ec3
# ╠═f782f312-75b3-11eb-3d7f-ff524bf48ece
# ╠═fddffb74-75b3-11eb-17e5-7f641e5114f9
# ╠═92ce58de-75b4-11eb-368a-8fc4b0c4be44
# ╠═146569e2-75b4-11eb-340f-5d12067a8e76
# ╠═1a6c44e6-75b4-11eb-3a94-6575eea73c03
