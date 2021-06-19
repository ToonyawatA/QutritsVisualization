### A Pluto.jl notebook ###
# v0.14.8

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

# ╔═╡ 9814729a-6080-11eb-0de8-bb08c7111606
using Plots, LinearAlgebra

# ╔═╡ ee16d9ee-6080-11eb-1564-07708d618050
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

# ╔═╡ f14477f2-6080-11eb-171f-f90d82b45c5d
#define dissipator operator
function D(sigma)
	d = size(sigma)[1]
	I0 = Matrix(1.0I,d,d)
	kron(conj(sigma),sigma)-0.5*kron(I0,transpose(conj(sigma))*sigma)-0.5*kron(transpose(transpose(conj(sigma))*sigma),I0)
end

# ╔═╡ fcb70890-6081-11eb-04d2-8dcd42374330
@bind Delta1 html"<input type=range min=0.01 max=10.0 step=0.05>"

# ╔═╡ 0b27b3b0-6082-11eb-0001-1fbf739823b5
@bind Delta2 html"<input type=range min=0.01 max=10.0 step=0.05>"

# ╔═╡ 0e81ae2c-6082-11eb-1bd3-ffbda94217d5
@bind Omega1 html"<input type=range min=0.01 max=10.0 step=0.05>"

# ╔═╡ 16160df6-6082-11eb-09aa-c584d4e325ea
@bind Omega2 html"<input type=range min=0.01 max=10.0 step=0.05>"

# ╔═╡ 60c36bdc-6082-11eb-19c8-050b6b1a6af0
@bind Gamma1 html"<input type=range min=0.01 max=3.0 step=0.05>"

# ╔═╡ 687dc818-6082-11eb-3666-d37237ed6373
@bind Gamma2 html"<input type=range min=0.01 max=3.0 step=0.05>"

# ╔═╡ 3956808e-6084-11eb-32c4-8fe5099d4234
@bind tt html"<input type=range min=0.01 max=30.0 step=0.05>"

# ╔═╡ 3f5e82e8-69fa-11eb-2739-79b806a14f89
Delta1

# ╔═╡ 434f404a-69fa-11eb-2bb0-012f828753ba
Delta2

# ╔═╡ 466a7074-69fa-11eb-1eef-33159384848c
Omega1

# ╔═╡ 4a00664e-69fa-11eb-3946-bd05562017d8
Omega2

# ╔═╡ 4c536b26-69fa-11eb-2a50-67221de3bf76
Gamma1

# ╔═╡ 508bdcf0-69fa-11eb-1755-45c0159b9275
Gamma2

# ╔═╡ 538fdf50-69fa-11eb-37b1-aba39acb965b
tt

# ╔═╡ 46e00ac4-608a-11eb-356d-15580d5653f8
a = Array{Float64}(undef,2,2,3)

# ╔═╡ 58525dba-6081-11eb-2f4a-2947c4998d6d
begin

	#define sigma_plus and sig_minus
	g1 = [1.0;0.0;0.0]
	g2 = [0.0;1.0;0.0]
	e = [0.0;0.0;1.0]
	sigma1 = g1*e'
	sigma2 = g2*e'

	#define Hamiltonian
	H_a = Delta1*g1*g1' + Delta2*g2*g2'
	H_af = 0.5*Omega1*(sigma1 + conj(sigma1')) + 0.5*Omega2*(sigma2 + conj(sigma2'))
	H = H_a + H_af

	#define superoperator
	I3 = Matrix(1.0I,3,3)
	H_eff = -1im*(kron(I3,H)-kron(H',I3))
	L_eff = Gamma1*D(sigma1)+Gamma2*D(sigma2)
	S = H_eff+L_eff #superoperator
	evals, evecs = eigen(S)

	#define time parameter
	num = 1000
	tmax = tt
	t = LinRange(0,tmax,num)
	dt = t[2]-t[1]

	#define initial state
	rho_0 = zeros(9)
	rho_0[9] = 1.0

	#quantum state evolution
	rho_ee = zeros(num)
	rho_gg1 = zeros(num)
	rho_gg2 = zeros(num)
	for i in 1:num
		rho_all = evecs*Diagonal(exp.(evals*t[i]))*inv(evecs)*rho_0
		rho_ee[i] = real(rho_all[9])
		rho_gg2[i] = real(rho_all[5])
		rho_gg1[i] = real(rho_all[1])
	end


	###################### quantum jump method ###############################

	#define jump operator
	c1 = sqrt(Gamma1)*sigma1
	C1 = c1'*c1

	c2 = sqrt(Gamma2)*sigma2
	C2 = c2'*c2

	#define effective non-Hermitain Hamiltonian
	nH_eff = H - (1im/2)*(C1 + C2)
	nojumpE = (I3 - (1im*nH_eff*dt))


	#define initial state
	Ns =300
	psi0 = e

	psi_en = zeros(Complex,3,1,Ns) #ensemble of states
	psi_en[:,:,:] .= psi0


	rhoj_t = zeros(Complex,3,3,num) #storing density matrix for each time steps

	# quantum jump method
	for i in 1:num

		rhoj = zeros(3,3) #using for update density matrix for each time step
		epsilon = rand(1,Ns)

		for j in 1:Ns
			dp1 = real(dt*psi_en[:,:,j]'*C1*psi_en[:,:,j])
			dp2 = real(dt*psi_en[:,:,j]'*C2*psi_en[:,:,j])
			dp = dp1+dp2

			#No jump evolution
			if epsilon[j] > dp[1]
				norm = real(psi_en[:,:,j]'*nojumpE'*nojumpE*psi_en[:,:,j])
				#psi_en[:,:,j] = (1/sqrt(1-dp[1]))*nojumpE*psi_en[:,:,j]
				psi_en[:,:,j] = (1/sqrt(norm[1]))*nojumpE*psi_en[:,:,j]
				rhoj += psi_en[:,:,j]*conj(psi_en[:,:,j]')

			#jump evolution
			else
				epsilon1 = rand()
				if epsilon1 < (dp1[1]/dp[1])
					psi_en[:,:,j] = sqrt(dt/dp1[1])*c1*psi_en[:,:,j]
				else
					psi_en[:,:,j] = sqrt(dt/dp2[1])*c2*psi_en[:,:,j]
				end
				rhoj += psi_en[:,:,j]*conj(psi_en[:,:,j]')
			end
		end
		rhoj_t[:,:,i] = (1/Ns)*rhoj
	end

	p_ee = abs.(rhoj_t[3,3,:])
	p_g2 = abs.(rhoj_t[2,2,:])
	p_g1 = abs.(rhoj_t[1,1,:])
	#for i in 1:num
	#	p_ee[i] = real(e'*rhoj_t[:,:,i]*e)
	#	p_g2[i] = real(g2'*rhoj_t[:,:,i]*g2)
	#	p_g1[i] = real(g1'*rhoj_t[:,:,i]*g1)
	#end

	l = @layout [a;b;c]
	p1 = plot(t,rho_ee,label="Prob in e state (OBE)")
	p1 = plot!(t,p_ee,label="Prob in e state (Jump)")
	p2 = plot(t,rho_gg2,label="Prob in g2 state (OBE)")
	p2 = plot!(t,p_g2,label="Prob in g2 state (Jump)")
	p3 = plot(t,rho_gg1,label="Prob in g1 state (OBE)")
	p3 = plot!(t,p_g1,label="Prob in g1 state (Jump)")

	p4 = plot(p1,p2,p3,layout=l,size=(450,900))
	#savefig("QuantunJump")

end


# ╔═╡ 69f62638-65b4-11eb-18ac-9f70b39ea388
a[:,:,1] = [1 1;1 1]

# ╔═╡ 7369bc1e-65b4-11eb-30f5-0d4b12bbf065
a[:,:,2] = [2 2;2 2]

# ╔═╡ 77178fa8-65b4-11eb-23a0-d3da6649ed32
a[:,:,3] = [3 3;3 3]

# ╔═╡ 81ca4650-65b4-11eb-3164-1f19d6e7e2b8
a

# ╔═╡ 6aeae542-657b-11eb-1620-15184a26b490
dp1

# ╔═╡ f9eba016-659d-11eb-056e-5192939c3ab2
(psi_en[:,:,1]'*C1*psi_en[:,:,1])[1]

# ╔═╡ 3da00b62-65b2-11eb-1944-6592aa4b11e9
zeros(1,2,3)'

# ╔═╡ 0fcfe3dc-657c-11eb-383a-5b5167f0889f
1.0 > 1.0 + 1im

# ╔═╡ ee2325b6-65b0-11eb-0c02-f730d25e8be3
yy

# ╔═╡ 0ca35ed6-6a64-11eb-134f-bf7082584ffb
xx = [1 1;1 -1]

# ╔═╡ 2115815c-6a64-11eb-2b95-37c0b7c358ba
eigva, eigve = eigen(xx)

# ╔═╡ 4b7846b2-6a64-11eb-1294-b7168ff2f973
dd = eigve[:,2]*eigve[:,2]'

# ╔═╡ 7eaeb9d0-6a64-11eb-284e-1d4940a88d01
conj(eigve')*dd*eigve

# ╔═╡ 6005472c-40f7-41e3-9eb3-91844365050b
aaa = [1 2 3 4]

# ╔═╡ 5f0a8311-a483-4aad-9e1e-eede4cb8df4b
aaa'*aaa

# ╔═╡ Cell order:
# ╠═9814729a-6080-11eb-0de8-bb08c7111606
# ╠═ee16d9ee-6080-11eb-1564-07708d618050
# ╠═f14477f2-6080-11eb-171f-f90d82b45c5d
# ╠═fcb70890-6081-11eb-04d2-8dcd42374330
# ╠═0b27b3b0-6082-11eb-0001-1fbf739823b5
# ╠═0e81ae2c-6082-11eb-1bd3-ffbda94217d5
# ╠═16160df6-6082-11eb-09aa-c584d4e325ea
# ╠═60c36bdc-6082-11eb-19c8-050b6b1a6af0
# ╠═687dc818-6082-11eb-3666-d37237ed6373
# ╠═3956808e-6084-11eb-32c4-8fe5099d4234
# ╠═3f5e82e8-69fa-11eb-2739-79b806a14f89
# ╠═434f404a-69fa-11eb-2bb0-012f828753ba
# ╠═466a7074-69fa-11eb-1eef-33159384848c
# ╠═4a00664e-69fa-11eb-3946-bd05562017d8
# ╠═4c536b26-69fa-11eb-2a50-67221de3bf76
# ╠═508bdcf0-69fa-11eb-1755-45c0159b9275
# ╠═538fdf50-69fa-11eb-37b1-aba39acb965b
# ╠═58525dba-6081-11eb-2f4a-2947c4998d6d
# ╠═46e00ac4-608a-11eb-356d-15580d5653f8
# ╠═69f62638-65b4-11eb-18ac-9f70b39ea388
# ╠═7369bc1e-65b4-11eb-30f5-0d4b12bbf065
# ╠═77178fa8-65b4-11eb-23a0-d3da6649ed32
# ╠═81ca4650-65b4-11eb-3164-1f19d6e7e2b8
# ╠═6aeae542-657b-11eb-1620-15184a26b490
# ╠═f9eba016-659d-11eb-056e-5192939c3ab2
# ╠═3da00b62-65b2-11eb-1944-6592aa4b11e9
# ╠═0fcfe3dc-657c-11eb-383a-5b5167f0889f
# ╠═ee2325b6-65b0-11eb-0c02-f730d25e8be3
# ╠═0ca35ed6-6a64-11eb-134f-bf7082584ffb
# ╠═2115815c-6a64-11eb-2b95-37c0b7c358ba
# ╠═4b7846b2-6a64-11eb-1294-b7168ff2f973
# ╠═7eaeb9d0-6a64-11eb-284e-1d4940a88d01
# ╠═6005472c-40f7-41e3-9eb3-91844365050b
# ╠═5f0a8311-a483-4aad-9e1e-eede4cb8df4b
