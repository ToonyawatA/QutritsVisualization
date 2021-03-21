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

# ╔═╡ d8440102-72af-11eb-35fc-fd767329d18f
using Plots, LinearAlgebra

# ╔═╡ 53612950-72b0-11eb-1580-99530d5661d5
# define Bloch sphere
function bloch()
	N=100
	a=0.150
	b=100
	lw=1
	lw1 = 5
	theta=LinRange(0.0,2*pi,N)
	phi=LinRange(0.0,2*pi,N)
	z0=0.0
	r = 1.0
	plot3d(r*sin.(theta).*cos.(0.0),r*sin.(theta).*sin.(0.0),r*cos.(theta),legend=false,line = (:black, lw, abs.(a.-sin.(phi)/b)),size=(1.0*600,1.0*520),axis=([], false),box=:on,camera = (40,40),xlim=(-1.0,1.0), ylim=(-1.0,1.0),zlim=(-1.0,1.0))
	for n in 2:17
		long=n*pi/18
		plot3d!(r*sin.(theta).*cos.(long),r*sin.(theta).*sin.(long),r*cos.(theta),line = (:black, lw, abs.(a.-sin.(phi)/b)))
	end
#	plot3d!(sin.(theta),z0*ones(N),cos.(theta),line = (:black, lw, abs.(a.-z0/b)))
#	plot3d!(z0*ones(N),sin.(theta),cos.(theta),line = (:black, lw, abs.(a.-sin.(theta)/b)))
	plot3d!(r*sin.(theta),r*cos.(theta),zeros(N),line = (:black, lw, abs.(a.-cos.(theta)/b)))
	z1= 2*r/3
	rho1= sqrt(r^2-z1*z1)
	x1=rho1*sin.(theta)
	y1=rho1*cos.(theta)
	plot3d!(x1,y1,z1*ones(N),line = (:black, lw, abs.(a.-cos.(theta)/b)))
	z2= -2*r/3
	rho2= sqrt(r^2-z2*z2)
	x2=rho2*sin.(theta)
	y2=rho2*cos.(theta)
	plot3d!(x2,y2,z2*ones(N),line = (:black, lw, abs.(a.-cos.(theta)/b)))
#   plot3d! unit circle behind
	alpha = 0
	x3 = r*cos.(theta)
	z3 = zeros(N)
	y3 = r*sin.(theta)
	#tune the angle of the circle by rotating with the rotation matrix
	x3 = cos(alpha)*x3 - sin(alpha)*y3
	y3 = sin(alpha)*x3 + cos(alpha)*y3

	plot3d!(x3,y3,z3,line = (:blue, lw1, abs.(a.-cos.(theta)/b)))
end

# ╔═╡ 531a52e4-72b0-11eb-2c01-79a3cb6155c5
bloch()

# ╔═╡ 72c17c30-72b0-11eb-3565-9580aab1e3ea
#define dissipator operator
function D(sigma)
	d = size(sigma)[1]
	I0 = Matrix(1.0I,d,d)
	kron(conj(sigma),sigma)-0.5*kron(I0,transpose(conj(sigma))*sigma)-0.5*kron(transpose(transpose(conj(sigma))*sigma),I0)
end

# ╔═╡ 62d76ff4-72b1-11eb-13e6-bf5c50fed494
#@bind Delta1 html"<input type=range min=0.01 max=20.0 step=0.05>"

# ╔═╡ 64bfbc04-72b1-11eb-3258-7dae266c6f17
@bind Omega1 html"<input type=range min=0.01 max=1.0 step=0.05>"

# ╔═╡ 6a2b4bec-72b1-11eb-0fed-e7b9ea504a48
@bind Omega2 html"<input type=range min=0.01 max=1.0 step=0.05>"

# ╔═╡ 70a96022-72b1-11eb-116d-f59e801b68d8
@bind Gamma1 html"<input type=range min=0.01 max=3.0 step=0.05>"

# ╔═╡ 7c7ebe28-72b1-11eb-0bc5-179f8c6c9cac
@bind Gamma2 html"<input type=range min=0.01 max=3.0 step=0.05>"

# ╔═╡ a9dae890-72b1-11eb-2665-896ab46586d0
@bind gammaG html"<input type=range min=0.01 max=3.0 step=0.05>"

# ╔═╡ 7f909044-72b1-11eb-2579-2f95f9065d99
@bind tt html"<input type=range min=0.01 max=20.0 step=0.05>"

# ╔═╡ 84173ff6-72b0-11eb-2da6-0555b1306cb3
begin
	#define sigma_plus and sig_minus
	g1 = [1.0;0.0;0.0]
	g2 = [0.0;1.0;0.0]
	e = [0.0;0.0;1.0]
	sigma1 = g1*e'
	sigma2 = g2*e'
	sigmaG = g2*g2' - g1*g1'
	sigmax = Complex[0 1 0;1 0 0;0 0 0]
	sigmay = Complex[0 -1im 0;1im 0 0;0 0 0]
	sigmaz = Complex[1 0 0;0 -1 0;0 0 0]
	sigma0 = Complex[0 0 0;0 0 0;0 0 1]


	#define detuning array and susceptibility
	numd = 1000
	Dmax = 5.0
	Delta1 = 0.0
	Delta2 = 0.0
	Delta = LinRange(-Dmax,Dmax,numd)
    Chi_r = zeros(numd) # real part of susceptibility
	Chi_i = zeros(numd) # imaginary part of susceptibility

	#define time parameter
	numt = 1000
	tmax = tt
	t = LinRange(0,tmax,numt)
	dt = t[2]-t[1]

	#define initial state
	rho_0 = zeros(9)
	rho_0[1] = 1.0


	#define Hamiltonian
	H_a = Delta1*g1*g1' + Delta2*g2*g2'
	H_af = 0.5*Omega1*(sigma1 + conj(sigma1')) + 0.5*Omega2*(sigma2 + conj(sigma2'))
	H = H_a + H_af

	#define superoperator
	I3 = Matrix(1.0I,3,3)
	H_eff = -1im*(kron(I3,H)-kron(conj(H'),I3))
	L_eff = Gamma1*D(sigma1)+Gamma2*D(sigma2)+gammaG*D(sigmaG)
	S = H_eff+L_eff #superoperator
	evals, evecs = eigen(S)


	# (A) Plot Delta against linear susceptibility

	for i in 1:numd
		Ha = Delta1*g1*g1' + Delta[i]*g2*g2'
	    Haf = 0.5*Omega1*(sigma1 + conj(sigma1')) + 0.5*Omega2*(sigma2 + conj(sigma2'))
	    HH = Ha + Haf
		Heff = -1im*(kron(I3,HH)-kron(conj(HH'),I3))
	    SS = Heff+L_eff #superoperator
		evals1, evecs1 = eigen(SS)
		rho_all = evecs1*Diagonal(exp.(evals1*tt))*inv(evecs1)*rho_0
		Chi_r[i] = real(rho_all[6])
		Chi_i[i] = -imag(rho_all[6])
	end



	# (B) quantum state evolution

	u1 = zeros(numt) #
	v1 = zeros(numt) ### Bloch vector for c1 and c2
	w1 = zeros(numt) #

	u2 = zeros(numt) #
	v2 = zeros(numt) ## circular vector for c3
	w2 = zeros(numt) #


	rho_ee = zeros(numt)
	rho_gg1 = zeros(numt)
	rho_gg2 = zeros(numt)

	for i in 1:numt
		rho_all = evecs*Diagonal(exp.(evals*t[i]))*inv(evecs)*rho_0
		phi = atan(-1*imag(rho_all[3])/real(rho_all[3]))
		rho_ee[i] = real(rho_all[9])
		rho_gg2[i] = real(rho_all[5])
		rho_gg1[i] = real(rho_all[1])
		DM = reshape(rho_all,3,3)
		u1[i] = real(tr(sigmax*DM))
		v1[i] = real(tr(sigmay*DM))
		w1[i] = real(tr(sigmaz*DM))
		u2[i] = real(tr(sigma0*DM)*cos(phi))
		v2[i] = real(tr(sigma0*DM)*sin(phi))
	end

	l = @layout [a;b;c]
	p1 = plot(t,rho_ee,xlabel="time",ylabel="Probability",label="Prob in e state (OBE)")
	p2 = plot!(t,rho_gg2,label="Prob in g2 state (OBE)")
	p2 = plot!(t,rho_gg1,label="Prob in g1 state (OBE)")
	p1 = plot(Delta,Chi_r,xlabel="Delta",ylabel="Chi",label="Dispersion")
	p1 = plot!(Delta,Chi_i,label="Absorption")
	p3 = bloch()
	p3 = scatter3d!(u1,v1,w1,ms=3,zcolor=t/tmax,m =(:heat, 0.25),markerstrokewidth=0)
	p3 = scatter3d!(u2,v2,w2,ms=3,zcolor=t/tmax,m =(:algae, 0.25),markerstrokewidth=0)
	p4 = plot(p1,p2,p3,layout=l,size=(600,1400))



end


# ╔═╡ Cell order:
# ╠═d8440102-72af-11eb-35fc-fd767329d18f
# ╠═53612950-72b0-11eb-1580-99530d5661d5
# ╠═531a52e4-72b0-11eb-2c01-79a3cb6155c5
# ╠═72c17c30-72b0-11eb-3565-9580aab1e3ea
# ╠═62d76ff4-72b1-11eb-13e6-bf5c50fed494
# ╠═64bfbc04-72b1-11eb-3258-7dae266c6f17
# ╠═6a2b4bec-72b1-11eb-0fed-e7b9ea504a48
# ╠═70a96022-72b1-11eb-116d-f59e801b68d8
# ╠═7c7ebe28-72b1-11eb-0bc5-179f8c6c9cac
# ╠═a9dae890-72b1-11eb-2665-896ab46586d0
# ╠═7f909044-72b1-11eb-2579-2f95f9065d99
# ╠═84173ff6-72b0-11eb-2da6-0555b1306cb3
