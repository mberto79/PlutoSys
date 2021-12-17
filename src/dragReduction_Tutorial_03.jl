### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# ╔═╡ f5851eb7-93ba-474f-ab32-5c1a79a4caf5
begin
	using PlutoUI
	using Downloads: download
	using Plots
	using LaTeXStrings
	using Latexify; set_default(fmt = FancyNumberFormatter(4), env=:equation)
	aerofluids = "http://aerofluids.org/notebooks/"
	download(aerofluids*"src/myPlutoUtils.jl") |> include
	styling = get_plot_style()
	toc = TableOfContents()
	md"*Cell: loading packages*"
end

# ╔═╡ 544fc8f0-5429-11ec-2ae3-a33713d9c95e
md"""
# DRAG REDUCTION $(html"<br>") *Tutorial sheet 3 - Lifting line (estimation)*

#### MMME3080 - Advanced Aerodynamics
"""

# ╔═╡ b09c0be0-42f3-4f81-9696-3f2be96d4e96
md"""
Q1. The distribution of circulation, ``\Gamma``, along a finite wing is required when using Prandtl's lifting line theory in order to determine the lift and induced drag characteristics. The distribution of ``\Gamma`` can be determined by assuming that it can represented with a Fourier series of the form

$$\Gamma(\theta) = (2bV_\infty) \sum_{n=1}^\infty A_n \sin n\theta$$

When the expression above is incoporated into the general wing equation, it can be shown that the induced drag coefficient is given by

```math
C_{D_i} =
\pi AR (A_1^2) \left[ 1 + \sum_{n=2}^\infty n \left(\frac{A_n}{A_1} \right)^2\right]
```

Use the equations above (and some mathematicl reasoning) to show that the induced drag coefficient is a minimum for an elliptical distribution of ``\Gamma`` and it can be represented as

```math
\frac{\Gamma^2}{\Gamma_0^2} + \frac{y^2}{(b/2)^2} = 1
```

Where ``\Gamma_0`` is the value of the circulation at the centre of the wing and ``b`` is the wing span.


"""

# ╔═╡ f18defae-abe6-4ee3-bd55-a22cdd66c7a7
md"""
**SOLUTION**

This can be found in the literature. The following book (included in the Bibliogrpahy for the lecture notes) uses a nice approach (however, you are encouraged to find alternative sources since this solution is presented differently depending on the author's preferences. Perhaps an alternative approach may suit your preferences better)

$(bibItem("Flandro, Gary A. McMahon, Howard M. Roach, Robert L.", 2012, "Basic Aerodynamics - Incompressible Flow", "Cambridge University Press"))
"""
# md"""
# (answer: proof required)
# """

# ╔═╡ a362cf2d-7879-4727-a4cb-b9ead98de139
md"""
Q2. For the special case of an elliptical wing, the downwash is constant and equal to

$$w = -\frac{\Gamma_0}{2b}$$

Where ``\Gamma_0`` is the value of the circulation at the centre of the wing and ``b`` is the wing span

* Show that the induced angle can be related to the lift coefficent, ``C_L``, and the apspect ration, ``AR``, by the following expression

$$\alpha_i = \frac{C_L}{\pi AR}$$

* Show that the induced lift coefficient for the wing is given by

$$C_{D_i} = \frac{C_L^2}{\pi AR}$$

* What are the implications of this result for aircraft designers?
"""

# ╔═╡ d52964b0-a58d-47ed-911d-c7da8fb18288
md"""
**SOLUTION**

This can be found in the literature. The following book (included in the Bibliogrpahy for the lecture notes) uses a nice approach (however, you are encouraged to find alternative sources since this solution is presented differently depending on the author's preferences. Perhaps an alternative approach may suit your preferences better)

$(bibItem("Flandro, Gary A. McMahon, Howard M. Roach, Robert L.", 2012, "Basic Aerodynamics - Incompressible Flow", "Cambridge University Press"))
"""
# md"""
# (answer: proof required)
# """

# ╔═╡ db1828a2-79ee-4ef3-819c-202253deea5b
begin	
	# Input information
	# V 		= 10
	alpha0_degrees = 1.5
	alpha_degrees = 15
	alpha0 	= deg2rad(alpha0_degrees)
	alpha 	= deg2rad(alpha_degrees)
	chord 	= 1.2
	b 		= 12
	m0 		= round(2*π, digits= 4)
	nPoints = 4
	n = [2i-1 for i ∈ 1:nPoints]

	# Define solution points
	bh = b/2
	DY = b/2/nPoints
	y = [0+DY/2:DY:b/2;]
	# bh = b/2
	# DY = b/nPoints
	# y = [-bh+DY/2:DY:bh;]
md"""
Q3. Using Prandtl's lifting line theory, determing the lift and drag coefficient for a finite wing for which the following properties are known

|Property   |  Value |
|:---|---|
| Zero angle of attack  | $(alpha0_degrees) degrees |
| Angle of attack  | $(alpha_degrees) degrees |
| Chord  | $(chord) metres |
| Wing span  | $(b) metres |
| Lit-curve slope  | $(m0) metres |

Perform the calculations at the following spanwise locations, ``y_0 \in`` $(" $(y)") along the positive semi-span.
"""
end

# ╔═╡ da72e4b2-89ee-4bd7-b974-efba6ea3ab60
begin
	# Input information
	# V 		= 10
	alpha0Tip_degrees = 1.0 	#1.5
	alphaTip_degrees = 10 		#15
	alpha0Tip 	= deg2rad(alpha0Tip_degrees)
	alphaTip 	= deg2rad(alphaTip_degrees)
	chordTip 	= 1.0 		#1.2
md"""
Q4. Using Prandtl's lifting line theory, determing the lift and drag coefficient for a finite wing for which the following properties are known

|Property   |  Value |
|:---|---|
| Zero angle of attack  | $(alpha0_degrees) degrees |
| Zero angle of attack (tip)  | $(alpha0Tip_degrees) degrees |
| Angle of attack  | $(alpha_degrees) degrees |
| Angle of attack (tip) | $(alphaTip_degrees) degrees |
| Chord  | $(chord) metres |
| Chord (tip) | $(chordTip) metres |
| Wing span  | $(b) metres |
| Lit-curve slope  | $(m0) metres |

Perform the calculations at the following spanwise locations, ``y_0 \in`` $(" $(y)") along the positive semi-span.

"""
end

# ╔═╡ 9ce1500a-601b-41b4-92fc-210aa993db5b
begin
	# Function definitions
	Theta(y, b) = acos(2*(y)/b) 

	function Chord(y, b, chordR; chordT=chordR)
	    m = (chordT - chordR)/(b/2)
	    cy = y*m + chordR
	    return cy
	end
	
	function Alpha(y, b, alphaR; alphaT=alphaR)
	    m = (alphaT - alphaR)/(b/2)
	    ay = y*m + alphaR
	    return ay
	end
	
	function Alpha0(y, b, alpha0R; alpha0T=alpha0R)
	    m = (alpha0T - alpha0R)/(b/2)
	    a0y = y*m + alpha0R
	    return a0y
	end
	
	function expansion_coeff(y, N, theta, c, m)
	    t1 = (4*b)/(m*c)
	    t2 = N/sin(theta)
	    t3 = t1 + t2
	    return t3*sin(N*theta)
	end

	function build_system(y, a, a0, theta, c, m)
		# Pre-allocate memory
		nPoints = length(y)
		C = zeros(nPoints, nPoints)
		D = zeros(nPoints)
		
		# Assign coefficient matrix
		for j ∈ 1:nPoints
		    for i ∈ 1:nPoints
		        C[i,j] = expansion_coeff(y[i], 2j-1, theta(y[i]), c(y[i]), m(y[i]))
		    end
		    D[j] = (a(y[j]) - a0(y[j]))
		end
		return C, D
	end
	
	function Γ(A, y, theta)
		nPoints = length(y)
	    G = zeros(nPoints)
	    for i ∈ 1:nPoints
	        for n ∈ 1:nPoints
	            G[i] += A[n]*sin((2*n-1)*theta(y[i]))
	        end
	    end
	    return G
	end

	# function δ(A)
	# 	N = length(A)
	# 	sum = 0.0
	# 	for n ∈ 2:N
	# 		sum += n*(A[n]/A[1])^2
	# 	end
	# 	return sum
	# end

	function δ(A)
		N = length(A)
		sum = 0.0
		for n ∈ 2:N
			sum += (2n-1)*(A[n]/A[1])^2
		end
		return sum
	end

	e(δ) = 1/(1+δ)
	AR(b, c; c_tip=c) = (b^2)/( b*(c + c_tip)/2 ) # Need to extend (to handle taper)

	function plot_gamma_distribution(y, G, CL, CD; styling=styling)
		gammaDistribution = plot(y, G, label="Circulation", 
			title="CL = $(latexify(CL; env=:raw)), CD = $(latexify(CD; env=:raw))",
			xlims=(0.0, b/2); styling...)
		xlabel!("Wing semi-span [m]")
		ylabel!("Circulation [m²/s]")
		scatter!(y, 0.0*y, color=:red, label="Grid points")
		
		# savefig(gammaDistribution, "fig_solution_example_20.svg")
		gammaDistribution
	end
md"*Cell: function definitions*"
end

# ╔═╡ b2acdc1b-0048-486c-8075-00ff0c6cc5d1
let
	# Function redefinitions
	θ(y) = Theta(y, b)
	c(y) = chord 
	a(y) = alpha
	m(y) = m0
	a0(y) = alpha0
	
	C, D = build_system(y, a, a0, θ, c, m)

	# Solve system
	A = similar(D)
	try
		A .= C\D
	catch
		print("NOT DOING THIS!")
	end
	
	# Calculate circulation and C_L from coefficients
	G = Γ(A, y, θ)
	CL = π*A[1]*AR(b,chord)
	CD = CL^2/(π*e(δ(A))*AR(b,chord))
	
	# Plot result
	p3 = plot_gamma_distribution(y, G, CL, CD)

md"""
**SOLUTION**

The aim is to find the distribution of ``\Gamma`` for the specific wing geometry given. This is achived using a Fourier series approximation for ``\Gamma`` and representing the fundamental wing equation in matrix form

```math
CA=D \qquad \text{(equation 1)}
```

Where ``C`` is the coefficient matrix, ``A`` is the matrix representing the unkonwn Fourier series expansion coefficients and ``D`` is a vector matrix of the known angles (geometric and aerodynamic).

The wing is to be evaluated at the given spanwise locations

``y_0 = `` $(latexify(y))

The following coordinate transformation is used

$$y = \frac{b}{2} \cos \theta$$

or in terms of the transformation angle (``\theta``)

$$\theta(y) = \cos^{-1} \left( \frac{2}{b} y \right)$$

Thus, evaluating for each element in ``y_0``, the corresponding ``\theta`` locations are

``\theta(y_0) = `` $(latexify(θ.(y)))

Since 4 stations are being considered (and only odd coefficients are required for symmetrical wing loading), 

``n = `` $(latexify(n))

Corresponding to the following unknown expansion coefficients (i.e. matrix ``A`` in equation 1)

``A = \begin{bmatrix} A_1 \\ A_3 \\ A_5 \\ A_7\end{bmatrix}`` 

The matrix of the known angles, ``D`` in equation 1, can be calculated as (requires conversion to radians)

```math
D(\theta(y)) = \alpha(\theta(y)) - \alpha_0(\theta(y))
```

Since ``\alpha =`` $(latexify(alpha)) and ``\alpha_0 =`` $(latexify(alpha0)), and both are constant along the span, the vector ``D`` becomes

``D = `` $(latexify(D))

The coefficient matrix, ``C`` can be calcuated using the equation

```math
C(\theta(y), n) = 
\left( 
\frac{4b}{m_0(\theta(y))c(\theta(y))}   
+
\frac{n}{\sin \theta(y)}
\right) \sin n\theta(y)
```

Since the gradient of the lift curve (``m_0``) and the chord (c), on this ocassion, are constant along the span. The equation above can be written as

```math
C(\theta(y_0), n) = 
\left( 
\frac{4b}{m_0c}   
+
\frac{n}{\sin \theta(y_0)}
\right) \sin n\theta(y_0)
```

Now, each location ``y_0``, corresponds to a value of ``\theta`` (as calculated above following the coordinate transformation). Also, for each location ``y_0`` there are ``n`` coefficients, where ``n =`` $(" $n "). It is now a matter of substituting these values into the equationa above. As an illustration, let us evaluate the first row of the coefficient matrix ``C``

``C(\theta(y_{0_1}), `` $(latexify(n[1])) ``) = C(`` $(latexify(θ.(y)[1])), $(latexify(n[1])) ``) = `` $(latexify(C[1,1]))

``C(\theta(y_{0_1}), `` $(latexify(n[2])) ``) = C(`` $(latexify(θ.(y)[2])), $(latexify(n[2])) ``) = `` $(latexify(C[1,2]))

``C(\theta(y_{0_1}), `` $(latexify(n[3])) ``) = C(`` $(latexify(θ.(y)[3])), $(latexify(n[3])) ``) = `` $(latexify(C[1,3]))

``C(\theta(y_{0_1}), `` $(latexify(n[4])) ``) = C(`` $(latexify(θ.(y)[4])), $(latexify(n[4])) ``) = `` $(latexify(C[1,4]))

Applying the same process to the remaining rows the full coefficient matrix ``C`` can be found. Therefore, the fundamental wing equation in matrix form becomes (fully expanded)

 $(latexify(C)) 
``\begin{bmatrix} A_1 \\ A_3 \\ A_5 \\ A_7\end{bmatrix} =`` 
$(latexify(D))

Which represents a system of linear equations and can be easily solved (e.g. using Gaussian elimination). Therefore, the coefficients we seek are

``A = `` $(latexify(A))

The *lift coefficient* can be found using the equation

``C_L = \pi A_1 AR``

Where ``AR = b^2/s`` is the aspect ratio, that is ``AR = ``$(latexify(AR(b,chord))). It is a matter of substituting for the first coefficient in matrix ``A`` in the equation above, we get

``C_L = \pi *`` $(latexify(A[1])) ``*`` $(latexify(AR(b,chord))) ``=``  $(latexify(CL) )

The *drag coefficient* is given by 

$$C_{D_i} = \frac{C_L^2}{\pi e AR}$$

Where the *span efficiency factor* is defined as

$$e = \frac{1}{(1 + \delta)}$$

and ``\delta`` is given by

``\delta =  \sum_{n=3}^N n \left(\frac{A_n}{A_1} \right)^2 = `` $(latexify(δ(A)))

Therefore, the *span efficiency factor* is

``e = `` $(latexify("1/(1 + $(δ(A))) = $(e(δ(A)))"))

Finally, the induced drag coefficient is

``C_{D_i} = \frac{C_L^2}{\pi e AR} =`` $(latexify(" $(CL)^2/(π* $(e(δ(A))) * $(AR(b,chord))) = $(CD)"))

Below a plot of the underlying distribution of ``\Gamma`` (calculated by subsituting back into the original Fourier series definition for ``\Gamma``)

$p3

"""
# md"""
# (answer: ``C_L = `` $(latexify(CL) ), ``C_D = `` $(latexify(CD) )
# """
end

# ╔═╡ 773517c7-c7b8-4f4a-971a-542462464ba1
let
# Function redefinitions
	θ(y) = Theta(y, b)
	c(y) = Chord(y, b, chord; chordT=chordTip)
	a(y) = Alpha(y, b, alpha; alphaT=alphaTip)
	m(y) = m0
	a0(y) = Alpha0(y, b, alpha0; alpha0T=alpha0Tip)
	
	C, D = build_system(y, a, a0, θ, c, m)

	# Solve system
	A = similar(D)
	try
		A .= C\D
	catch
		print("NOT DOING THIS!")
	end
	
	# Calculate circulation and C_L from coefficients
	G = Γ(A, y, θ)
	ARv = AR(b,chord, c_tip=chordTip)
	CL = π*A[1]*ARv
	δv = δ(A)
	ev = e(δv)
	CD = CL^2/(π*ev*ARv)
	
	# Plot result
	p4 = plot_gamma_distribution(y, G, CL, CD)
	
md"""
**SOLUTION**

The aim is to find the distribution of ``\Gamma`` for the specific wing geometry given. This is achived using a Fourier series approximation for ``\Gamma`` and representing the fundamental wing equation in matrix form

```math
CA=D \qquad \text{(equation 1)}
```

Where ``C`` is the coefficient matrix, ``A`` is the matrix representing the unkonwn Fourier series expansion coefficients and ``D`` is a vector matrix of the known angles (geometric and aerodynamic).

*NOTE* the solution proceeds in a similar fashion as shown previously. The key difference in this question is that some properties are changing linearly with the spanwise coordinate. Therefore, those wing properties that change along the span will not be treated as a vector. For example, let us consider the chord, a suitable expression must be found to represent the change

$$c(y) = c_{root} + \frac{(c_{tip} - c_{root})}{b/2}$$

The wing is to be evaluated at the given spanwise locations

``y_0 = `` $(latexify(y))

At each location we can evaluate the chord, giving the vector

``c(y) =`` $(latexify(c.(y)))

Similarly, the geometry and zero-lift angles change with the spanwise direction, after working out the linear equation and converting to radians the corresponding vectors are

``\alpha = `` $(latexify(a.(y))) and ``\alpha_0 = `` $(latexify(a0.(y)))

With the information above, we can proceed the analysis. The following coordinate transformation is used

$$y = \frac{b}{2} \cos \theta$$

or in terms of the transformation angle (``\theta``)

$$\theta(y) = \cos^{-1} \left( \frac{2}{b} y \right)$$

Thus, evaluating for each element in ``y_0``, the corresponding ``\theta`` locations are

``\theta(y_0) = `` $(latexify(θ.(y)))

Since 4 stations are being considered (and only odd coefficients are required for symmetrical wing loading), 

``n = `` $(latexify(n))

Corresponding to the following unknown expansion coefficients (i.e. matrix ``A`` in equation 1)

``A = \begin{bmatrix} A_1 \\ A_3 \\ A_5 \\ A_7\end{bmatrix}`` 

The matrix of the known angles, ``D`` in equation 1, can be calculated as (requires conversion to radians)

```math
D(\theta(y)) = \alpha(\theta(y)) - \alpha_0(\theta(y))
```

Since ``\alpha`` and ``\alpha_0 =`` are no longer constant along the span, the vector ``D`` becomes

``D = `` $(latexify(D))

The coefficient matrix, ``C`` can be calcuated using the equation

```math
C(\theta(y), n) = 
\left( 
\frac{4b}{m_0(\theta(y))c(\theta(y))}   
+
\frac{n}{\sin \theta(y)}
\right) \sin n\theta(y)
```

Since the gradient of the lift curve (``m_0``) constant along the span, but the chord, on this ocassion, is changeingalong the span. The equation above can be re-written as

```math
C(\theta(y_0), n) = 
\left( 
\frac{4b}{m_0c(\theta(y_0))}   
+
\frac{n}{\sin \theta(y_0)}
\right) \sin n\theta(y_0)
```

Now, each location ``y_0``, corresponds to a value of ``\theta`` (as calculated above following the coordinate transformation). Also, for each location ``y_0`` there are ``n`` coefficients, where ``n =`` $(" $n "). It is now a matter of substituting these values into the equationa above. As an illustration, let us evaluate the first row of the coefficient matrix ``C``

``C(\theta(y_{0_1}), `` $(latexify(n[1])) ``) = C(`` $(latexify(θ.(y)[1])), $(latexify(n[1])) ``) = `` $(latexify(C[1,1]))

``C(\theta(y_{0_1}), `` $(latexify(n[2])) ``) = C(`` $(latexify(θ.(y)[2])), $(latexify(n[2])) ``) = `` $(latexify(C[1,2]))

``C(\theta(y_{0_1}), `` $(latexify(n[3])) ``) = C(`` $(latexify(θ.(y)[3])), $(latexify(n[3])) ``) = `` $(latexify(C[1,3]))

``C(\theta(y_{0_1}), `` $(latexify(n[4])) ``) = C(`` $(latexify(θ.(y)[4])), $(latexify(n[4])) ``) = `` $(latexify(C[1,4]))

Applying the same process to the remaining rows the full coefficient matrix ``C`` can be found. Therefore, the fundamental wing equation in matrix form becomes (fully expanded)

 $(latexify(C)) 
``\begin{bmatrix} A_1 \\ A_3 \\ A_5 \\ A_7\end{bmatrix} =`` 
$(latexify(D))

Which represents a system of linear equations and can be easily solved (e.g. using Gaussian elimination). Therefore, the coefficients we seek are

``A = `` $(latexify(A))

The *lift coefficient* can be found using the equation

``C_L = \pi A_1 AR``

Where ``AR = b^2/s`` is the aspect ratio, that is ``AR = ``$(latexify(ARv)). It is a matter of substituting for the first coefficient in matrix ``A`` in the equation above, we get

``C_L = \pi *`` $(latexify(A[1])) ``*`` $(latexify(ARv)) ``=``  $(latexify(CL) )

The *drag coefficient* is given by 

$$C_{D_i} = \frac{C_L^2}{\pi e AR}$$

Where the *span efficiency factor* is defined as

$$e = \frac{1}{(1 + \delta)}$$

and ``\delta`` is given by

``\delta =  \sum_{n=3}^N n \left(\frac{A_n}{A_1} \right)^2 = `` $(latexify(δ(A)))

Therefore, the *span efficiency factor* is

``e = `` $(latexify("1/(1 + $(δv)) = $(ev)"))

Finally, the induced drag coefficient is

``C_{D_i} = \frac{C_L^2}{\pi e AR} =`` $(latexify(" $(CL)^2/(π* $(ev) * $(ARv)) = $(CD)"))

Below a plot of the underlying distribution of ``\Gamma`` (calculated by subsituting back into the original Fourier series definition for ``\Gamma``)

$p4

"""
# md"""
# (answer: ``C_L = `` $(latexify(CL) ), ``C_D = `` $(latexify(CD) )
# """
end

# ╔═╡ 4fba3872-9dba-4a3e-9441-d90b234884e6
# begin
# 	Markdown.html(io::IO, ls::LaTeXString) = 
# 		Markdown.html(io, Markdown.LaTeX(
# 			repr(MIME"text/latex"(), ls)
# 		))
	
# 	Markdown.htmlinline(io::IO, ls::LaTeXString) = 
# 		Markdown.htmlinline(io, Markdown.LaTeX(
# 			repr(MIME"text/latex"(), ls)
# 		))
# md"*Cell: Markdown sorcery*"
# end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Downloads = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
Latexify = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
LaTeXStrings = "~1.3.0"
Latexify = "~0.15.9"
Plots = "~1.25.2"
PlutoUI = "~0.7.23"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "abb72771fd8895a7ebd83d5632dc4b989b022b5b"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.2"

[[Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "84918055d15b3114ede17ac6a7182f68870c16f7"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.1"

[[ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "f2202b55d816427cd385a9a4f3ffb226bee80f99"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+0"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "4c26b4e9e91ca528ea212927326ece5918a04b47"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.11.2"

[[ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "bf98fa45a0a4cee295de98d4c1462be26345b9a1"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.2"

[[ColorSchemes]]
deps = ["ColorTypes", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "a851fec56cb73cfdf43762999ec72eff5b86882a"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.15.0"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "024fe24d83e4a5bf5fc80501a314ce0d1aa35597"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.0"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "44c37b4636bc54afac5c574d2d02b625349d6582"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.41.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[DataAPI]]
git-tree-sha1 = "cc70b17275652eb47bc9e5f81635981f13cea5c8"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.9.0"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "3daef5523dd2e769dad2365274f760ff5f282c7d"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.11"

[[DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b3bfd02e98aedfa5cf885665493c5598c350cd2f"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.2.10+0"

[[FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d8a578692e3077ac998b50c0217dfd67f21d1e5f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.0+0"

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "0c603255764a1fa0b61752d2bec14cfbd18f7fe8"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.5+1"

[[GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "30f2b340c2fff8410d89bfcdc9c0a6dd661ac5f7"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.62.1"

[[GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "fd75fa3a2080109a2c0ec9864a6e14c60cca3866"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.62.0+0"

[[GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "58bcdf5ebc057b085e58d95c138725628dd7453c"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.1"

[[Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "74ef6288d071f58033d54fd6708d4bc23a8b8972"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+1"

[[Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "0fa77022fe4b511826b39c894c90daf5fce3334a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.17"

[[HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[HypertextLiteral]]
git-tree-sha1 = "2b078b5a615c6c0396c77810d92ee8c6f470d238"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.3"

[[IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[IniFile]]
deps = ["Test"]
git-tree-sha1 = "098e4d2c533924c921f9f9847274f2ad89e018b8"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.0"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "a7254c0acd8e62f1ac75ad24d5db43f5f19f3c65"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.2"

[[IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "642a199af8b68253517b80bd3bfd17eb4e84df6e"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.3.0"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "8076680b162ada2a031f707ac7b4953e30667a37"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.2"

[[JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d735490ac75c5cb9f1b00d8b5509c11984dc6943"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.0+0"

[[LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "a8f4f279b6fa3c3c4f1adadd78a621b13a506bce"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.9"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "340e257aada13f95f98ee352d316c3bed37c8ab9"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.3.0+0"

[[Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "be9eef9f9d78cecb6f262f3c10da151a6c5ab827"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.5"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[NaNMath]]
git-tree-sha1 = "bfe47e760d60b82b66b61d2d44128b62e3a369fb"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.5"

[[NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7937eda4681660b4d6aeeecc2f7e1c81c8ee4e2f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+0"

[[OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "15003dcb7d8db3c6c857fda14891a539a8f2705a"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.10+0"

[[Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "ae4bbcadb2906ccc085cf52ac286dc1377dceccc"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.1.2"

[[Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PlotThemes]]
deps = ["PlotUtils", "Requires", "Statistics"]
git-tree-sha1 = "a3a964ce9dc7898193536002a6dd892b1b5a6f1d"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "2.0.1"

[[PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "e4fe0b50af3130ddd25e793b471cb43d5279e3e6"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.1.1"

[[Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun"]
git-tree-sha1 = "65ebc27d8c00c84276f14aaf4ff63cbe12016c70"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.25.2"

[[PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "5152abbdab6488d5eec6a01029ca6697dff4ec8f"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.23"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00cfd92944ca9c760982747e9a1d0d5d86ab1e5a"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.2.2"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "ad368663a5e20dbb8d6dc2fddeefe4dae0781ae8"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+0"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "7ad0dfa8d03b7bcf8c597f59f5292801730c55b8"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.4.1"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "8f82019e525f4d5c669692772a6f4b0a58b06a6a"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.2.0"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "3c76dde64d03699e074ac02eb2e8ba8254d428da"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.2.13"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StatsAPI]]
git-tree-sha1 = "0f2aa8e32d511f758a2ce49208181f7733a0936a"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.1.0"

[[StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "2bb0cb32026a66037360606510fca5984ccc6b75"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.13"

[[StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "2ce41e0d042c60ecd131e9fb7154a3bfadbf50d3"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.3"

[[TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "TableTraits", "Test"]
git-tree-sha1 = "bb1064c9a84c52e277f1096cf41434b675cd368b"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.6.1"

[[Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "66d72dc6fcc86352f01676e8f0f698562e60510f"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.23.0+0"

[[XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "1acf5bdf07aa0907e0a37d3718bb88d4b687b74a"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.12+0"

[[XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "cc4bf3fdde8b7e3e9fa0351bdeedba1cf3b7f6e6"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.0+0"

[[libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "c45f4e40e7aafe9d086379e5578947ec8b95a8fb"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+0"

[[nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╟─544fc8f0-5429-11ec-2ae3-a33713d9c95e
# ╟─b09c0be0-42f3-4f81-9696-3f2be96d4e96
# ╟─f18defae-abe6-4ee3-bd55-a22cdd66c7a7
# ╟─a362cf2d-7879-4727-a4cb-b9ead98de139
# ╟─d52964b0-a58d-47ed-911d-c7da8fb18288
# ╟─db1828a2-79ee-4ef3-819c-202253deea5b
# ╟─b2acdc1b-0048-486c-8075-00ff0c6cc5d1
# ╟─da72e4b2-89ee-4bd7-b974-efba6ea3ab60
# ╟─773517c7-c7b8-4f4a-971a-542462464ba1
# ╟─9ce1500a-601b-41b4-92fc-210aa993db5b
# ╟─f5851eb7-93ba-474f-ab32-5c1a79a4caf5
# ╟─4fba3872-9dba-4a3e-9441-d90b234884e6
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
