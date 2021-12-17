### A Pluto.jl notebook ###
# v0.17.3

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 2500b290-30fb-11ec-2224-ed2674139124
begin
	using PlutoUI
	using Downloads: download
	using Plots
	using LaTeXStrings
	aerofluids = "http://aerofluids.org/notebooks/"
	download(aerofluids*"src/myPlutoUtils.jl") |> include
	styling = get_plot_style()
	toc = TableOfContents()
	md"*Cell: loading packages*"
end

# ╔═╡ e3b13a64-4c36-46b3-8532-575a71ad3dee
md"""
# DRAG REDUCTION $(html"<br>") *Induced drag (estimation)*

#### MMME3080 - Advanced Aerodynamics
Dr Humberto Medina (*revision version 1.0*)
"""

# ╔═╡ 81548881-1de0-4c19-898d-4b85cb71cff1
toc

# ╔═╡ 0267d4be-7a71-4b21-a7a3-48e0bb86bb64
md"""
## Introduction
"""

# ╔═╡ 620122a7-7623-43ba-9ee1-a765786d1e93
md"""
Previously, the building blocks and general formulation of the lifting line theory have been introduced. It has been shown that, using the Biot-Savart law, it is possible to estimate the downwash generated at any point in a wing when it is modelled as a horseshoe vortex (and vortex sheet). However, the fundamental model that was formulated can only be of practical use if the distribution of the circulation along the bound vortex (lifting line) can be determined. The purpose of this lecture is to present a general method for estimating the distribution of the circulation for an arbiratry wing geometry. Thus, key performance infomation about the wing planform can be estimated. 

To this end, in this lecture the fundamental wing equation will be deduced. The next step will be formulate a general approach to estimate the distribution of circulation for an arbitrary wing, along with a general solution method (in the form of a fourier series approximation). This will allow to estimate the lift, induced drag, induced angle and downwash velocity distribution of any wing geometry (provided it adheres to the limits imposed by the theory).
"""

# ╔═╡ b25aff1f-d0f4-4003-95fc-f512f4113f2b
md"""
## Lecture video
"""

# ╔═╡ fe2a2cd8-99ba-4b5f-ae0c-c0f096002108
html"""
<iframe height="400" width="640" allowfullscreen frameborder=0 src="https://echo360.org.uk/media/02930a83-8f2b-45f8-8f22-74124acbf8a7/public?autoplay=false&automute=false"></iframe>
"""

# ╔═╡ c0bd19e4-a932-41c3-8a7e-45a326593e27
md"""
## Prandtl Lifting Line Theory (application)
"""

# ╔═╡ 98576fe8-cd61-4d8b-9628-1ff8e83ccffb
md"""
### Fundamental wing equation
"""

# ╔═╡ 8d6f9c8e-5b66-4d42-8836-a4e484dd7cbf
md"""
The objective in this section is to use fundamental aerodynamic concepts in order to modify the basic horseshoe and vortex sheet model to enable its application to finite wings of arbitrary geometry. To do this, the fundamental equation will have to be manipulated to allow the calculation of circulation over the bound vortex. The modifications start by writing the equation for the lift coefficient in terms of the circulation. The lift coefficient at a given spanwise location ``y_0`` can be written as
"""

# ╔═╡ f80f03cb-f5cb-40eb-9483-48d7ac6b33bc
md"""
```math
C_l(y_0) = \frac{L^\prime(y_0)}{1/2\rho V_\infty^2 c(y_0)} = \frac{2\Gamma(y_0)}{V_\infty c(y_0)} \qquad \text{(equation 1)}
```
"""

# ╔═╡ d13bdb89-4c57-4a08-b764-698b8bece8b8
md"""
In the equation above, the lift per unit span, ``L^\prime``, is replaced by the expression given by the Kutta-Joukousky theorem, ``L^\prime = 2V_\infty\Gamma``. Notice, that the equation above is evaluated at a fixed wing span location, ``y_0``. 
"""

# ╔═╡ 9a637330-680c-43db-881f-cf2198913640
md"""
A general wing planform may exhibit *geometric* and/or *aerodynamic twist* along its span. Thus, they could both be functions of the wing span, ``y``. *Geometric twist* refers to changes in the physical angle at which the local aerofoil section is positioned on the wing. It will be denoted as ``\alpha(y_0)`` to highlight the fact that it can change at different locations ``y_0``. The *aerodynamic twist* relates to changes to the aerodynamic response of the aerofoil with angle of attack (e.g. using a different aerofoil profile). Changes in aerodynamic twist can be incorporated using the zero lift angle for the given aerofoil profile/section, and it will be denoted as ``\alpha_0(y_0)`` (again, as a function of ``y_0`` to stress the fact that this angle can also change along the span of the wing). With these definitions in place, and assuming that the *lift slope* (lift coefficient per radians) of the aerfoil, ``m_0(y_0)``, is also known, it is possible to write an additional expression for the lift coefficient in terms of the lift slope and the characteristic angles (geometric, aerodynamic and induced). Recall that the effective angle of attack (from geometric inspection) is given as ``\alpha_{eff}(y_0) = \alpha(y_0) - \alpha_i(y_0)``, where ``\alpha_i(y_0)`` is the induced angle. It is now possible to express the lift coefficient as
"""

# ╔═╡ 0215c2f4-892f-40c4-acdd-f0e90d850c41
md"""
$$C_l(y_0) = m_0(y_0)[\alpha(y_0) - \alpha_i(y_0) - \alpha_0(y_0)]$$
"""

# ╔═╡ 1d221e4c-a532-4bcd-9b39-9e708daf0aa7
md"""
The equation above can be rewritten as
"""

# ╔═╡ d946d4ed-847e-4f1f-ae27-357435aa68b6
md"""
$$\alpha(y_0) = \frac{C_l(y_0)}{m_0(y_0)} + \alpha_i(y_0) + \alpha_0(y_0)$$
"""

# ╔═╡ 56f1ea76-f726-4327-87b2-b43192874401
md"""
Substituting in for the expression of the lift coefficient (equation 1) and for the expression for the induced angle derived in the previous lecture. We obtain the fundamental wing equation
"""

# ╔═╡ 8389b67f-86ee-4f86-90da-90ccf4dc450a
md"""
```math
\alpha(y_0) = 
\frac{2\Gamma(y_0)}{m_0(y_0) V_\infty c(y_0)} 
-
\frac{1}{4\pi V_\infty} \int_{-b/2}^{b/2}\frac{(d\Gamma/dy)|_{y_0}}{ (y - y_0)}dy 
+ 
\alpha_0(y_0) \qquad \text{(equation 2)}
```
"""

# ╔═╡ 67b75aa2-be8b-49f3-be61-bb5e75477b42
md"""
This is a very elegant reformulation of the theory. By some careful manipulation of existing results, it has been possible to arrive at an equation where the only unknown is not the distribution of the circulation, ``\Gamma``. Now, this differential-integral equation has no simple solution. Nonetheless, it will now be shown that it can be readily "solved" using Fourier series. 
"""

# ╔═╡ 55e66d16-8759-49e5-8c18-de0d41319f48
md"""
### Spanwise lift distribution for an arbitrary wing
"""

# ╔═╡ ef93499e-d302-45bd-b3d1-43f4a768ae41
md"""
The key to solving the fundamental wing equation shown earlier is to assume that the distribution of the circulation, ``\Gamma``, along the wing can be expressed as a Fourier series. It will be give the following form
"""

# ╔═╡ 1f5cda8a-d26f-4d5c-b1a5-29f3c9cf2ff7
md"""
$$\Gamma(\theta) = (2bV_\infty) \sum_{n=1}^\infty A_n \sin n\theta \qquad \text{(equation 3)}$$
"""

# ╔═╡ 7ce50ac3-de76-4044-a3bc-d9046d60ac53
md"""
In the equation above the contant term ``(2bV_\infty)`` is included to ensure units are consistent. Effectively, the problem has been reduced from trying to find a continuous solution for ``\Gamma`` to a discrete problem in terms of a new spatial coordinate ``\theta``, where the expansion coeffients ``A_n`` are needed. First, the transformation used is illustrated in the figure below
"""

# ╔═╡ e1d220ba-4e19-42fd-9797-67b3fecde01d
Resource("http://aerofluids.org/notebooks/dragReduction/03_theta_transformation.svg")

# ╔═╡ f86f3db6-a93c-45f6-8250-6b77c2cd3157
md"""
Examining the figure above, it can be shown that the relationship between ``y`` and ``\theta`` is as follows
"""

# ╔═╡ e710ee7d-1110-4a91-b4a2-073afa13e218
md"""
$$y = \frac{b}{2} \cos \theta$$
"""

# ╔═╡ acba8da7-e635-416f-9105-892344953e35
md"""
Now, we turn our attention to the second term in equation 2 (the integral-differential term). It includes the derivative of the circulation with respect to the spanwise coordinate so it can be expressed as
"""

# ╔═╡ 5d62a58f-2c10-4997-a0ef-87f8b99a026a
md"""
$$\frac{d\Gamma}{dy} = \frac{d\Gamma}{d\theta}\frac{d\theta}{dy}$$
"""

# ╔═╡ eeed306f-c59b-4586-9e7f-70caf49e39f9
md"""
After incorporating the coordinate transformation, the second term in equation 2 can be written as
"""

# ╔═╡ 0d085718-6bfb-4e0d-a08b-acab43870a4b
md"""
```math
\frac{1}{4\pi V_\infty} \int_{-b/2}^{b/2} \frac{d\Gamma/dy}{(y_0 - y)} dy = 
\frac{1}{\pi}
\int_\pi^0 
\frac{d\Gamma/d\theta}{\frac{b}{2}(\cos \theta_0 - \cos\theta)} 
\frac{d\theta}{dy}dy
```
"""

# ╔═╡ f47d77a9-652d-4d54-a12c-a7252fb6cfe5
md"""
Evaluating the derivative ``d\Gamma/d\theta`` from our assumed Fourier series (equation 3), changing the order of the integration limits, integrating (requires using the principal value approach due to the ill-behaved denominator for ``\theta=\theta_0``), and rearranging, we finally arrive at an equivalent expression for the second term on the right-hand side of equation 2. Given as
"""

# ╔═╡ 6c77ad6b-606a-4050-b136-312f667b84a1
md"""
$$\sum_{n=1}^\infty nA_n \frac{\sin n\theta}{\sin \theta} \qquad \text{(equation 4)}$$
"""

# ╔═╡ c0930aa3-b436-4ec4-bdba-62f044298a2f
md"""
Substituting equation 3 for ``\Gamma(y)`` in the first term on the right-hand side of equation 2, and replacing the second term of equation 2 for equation 4, we arrive at a new form of the fundamental wing equation in terms of the transformation variable theta, as shown below
"""

# ╔═╡ 06aafb04-1347-4653-abe6-6c6f0c2fbf93
md"""
```math
\alpha(\theta) = 
\frac{4b}{m_0(\theta)c(\theta)} \sum_{n=1}^\infty A_n \sin n\theta
+
\sum_{n=1}^\infty nA_n \frac{\sin n\theta}{\sin \theta}
+
\alpha_0(\theta)
\qquad \text{(equation 5)}
```
"""

# ╔═╡ ba451210-b849-424d-ac77-37a0ef98ee1a
md"""
With this equation it is now possible to approximate the required distribution of circulation to satisfy the wings geometry and aerodynamic constraints. A solution requires finding the unknown coefficients ``A_n`` where ``n`` is an integer values in the range ``n \in [1,\infty]``. However, for *symmetrical* wing loading, it can be shown that only odd values of ``n`` are non-zero. 
"""

# ╔═╡ 0212b76d-2304-45b1-b23f-0228e8ab39b8
md"""
**Lift coefficient**
"""

# ╔═╡ 4164f7fa-96df-4cf5-bd79-4fcff7a5f660
md"""
From the Kutta-Joukosky theorem, an expression for the lift per unit span is give. Thus, it is possible to obtain the total lift force byintegrating over the entire wing
"""

# ╔═╡ 531b2b1a-903b-4cab-b4dd-82827a6721be
md"""
$$L = \int_{-b/2}^{b/2} L^\prime dy = \rho V_\infty \int_{-b/2}^{b/2} \Gamma(y) dy$$
"""

# ╔═╡ e53dbc90-4840-44e1-8e4b-7c9835f10fb2
md"""
Therefore, the lift coefficient for the wing (notice the use of capital L for the subscript to denote a wing)
"""

# ╔═╡ 0405524a-7915-41e4-a264-d545560f44a3
md"""
$$C_L = \frac{L}{\frac{1}{2}\rho V_\infty^2 S} = \frac{2}{V_\infty S} \int_{-b/2}^{b/2} \Gamma(y) dy$$
"""

# ╔═╡ e806a0a9-dcdf-40f3-8ba6-870866c979ed
md"""
Applying the coordinate transformation and integrating
"""

# ╔═╡ 751648da-3723-44ae-8a19-683437969922
md"""
$$C_L = \frac{2}{V_\infty S} \int_{\pi}^{0} \left[ (2bV_\infty) \sum_{n=1}^\infty A_n \sin n\theta\right] \left[ (\frac{b}{2})(-\sin \theta d\theta)\right]$$
"""

# ╔═╡ 58026876-f278-4211-ae68-042fab586bdc
md"""
On inspection, the product ``(\sin n\theta)(\sin \theta)`` is zero when integrated between ``\pi`` and ``0``. The only remaining term constains the ``A_1`` coefficient. Therefore,
"""

# ╔═╡ e64bb91c-4105-479f-bb82-e2ddff6eda23
md"""
$$C_L = \pi A_1 AR$$
"""

# ╔═╡ 436845d8-72fb-49f7-abaa-45e141a8b3b5
md"""
The observation that ``(\sin n\theta)(\sin theta)`` is zero when integrated ``\in [\pi, 0]`` can be explored visually below.
"""

# ╔═╡ 31d527d1-03ad-4ffc-83a0-59b3ee9f8d84
md"""
*Select integer* (n) $(html"<br>") 
$(@bind n Slider(1:1:20, show_value=true))
"""

# ╔═╡ 81127462-87c7-4c58-8ef8-097101f2fc5d
let
	theta = [0:pi/150:pi;]
	a = sin.(theta)
	b = sin.(n.*theta)
	plot(theta, a, ylims=(-1,1), label=L"\sin \theta", legend=:bottomleft;styling...)
	plot!(theta, b, label=L"\sin n\theta";styling...)
	plot!(theta, a.*b, label=L"\sin \theta \sin n \theta";styling...)
end

# ╔═╡ d63c8e12-def8-4d56-85e9-8457a4a79c09
md"""
**Drag coefficient (induced)**
"""

# ╔═╡ f5cdc961-feeb-478a-8b2a-9b45c649fbe2
md"""
From geometry inspection, it has already been established that the induced drag per unit span is ``D^\prime_i = L^\prime\alpha_i``. Therefore, the total induced drag can now be calculated as follows
"""

# ╔═╡ c89ce8fb-0ac3-4360-b7c2-0fb6a05fcb93
md"""
$$D_i = \int_{-b/2}^{b/2} L^\prime\alpha_i dy$$
"""

# ╔═╡ 95be0392-e4ae-42c2-910b-e7950afb3f51
md"""
After substitution using the Kutta-Joukousky theorem gives
"""

# ╔═╡ 7cb6536a-66e5-4156-bfbc-3fe3d61f9966
md"""
$$D_i = \int_{-b/2}^{b/2} \rho V_\infty \Gamma\alpha_i dy$$
"""

# ╔═╡ 4cea8ab8-05d9-4509-acec-e22442aa21ae
md"""
Recal that both the circulation and the induced angle are both functions of ``y``, ``\Gamma(y)`` and ``\alpha_i(y)``, respectively. Substituting for their corresponding Fourier series expressions (equations 3 and 4) and applying the previous coordinate transformation for ``dy``, we obtain the expression for the induced drag
"""

# ╔═╡ ee4ce3ca-5870-4292-864d-24ddfcbaa249
md"""
```math
\int_\pi^0 
\left[\rho V_\infty (2bV_\infty) \sum_{n=1}^\infty A_n \sin n\theta \right]
\left[\sum_{n=1}^\infty nA_n \frac{\sin n\theta}{\sin \theta} \right]
\left[-\frac{b}{2} \sin\theta d\theta \right]
```
"""

# ╔═╡ af702e5a-2692-47e4-aee7-a334f8c163fe
md"""
And the induced drag coefficient becomes
"""

# ╔═╡ 7af5aa5a-6412-4dac-835b-c21475fc04b3
md"""
```math
C_{D_i} = \frac{D_i}{\frac{1}{2}\rho V_\infty^2 S} = 
2\frac{b^2}{S} \int_0^\pi 
\left[ \sum_{n=1}^\infty A_m \sin m\theta \right]
\left[\sum_{n=1}^\infty nA_n \sin n\theta \right]
d\theta
```
"""

# ╔═╡ 81de1481-3c19-4112-aa8e-ea26231dbdd3
md"""
The integrand is the product of two summations with the following properties
"""

# ╔═╡ 3cd73b58-7045-42a4-bbe2-e9213d6b0432
md"""
```math
\int_0^\pi \sin n\theta\sin m\theta d\theta = 
\begin{cases}
	0 \quad\text{ for } n\neq m\\
	\frac{\pi}{2} \quad\text{for } n=m
\end{cases}
```
"""

# ╔═╡ 6bdff171-e866-4059-a0f5-b6863e7b1729
md"""
This allows the expection with only one index, say ``n``, so it is possible to expand and evaluate the expression above for the two Fourier expansions which results in
"""

# ╔═╡ f92b9592-51e3-45b1-a497-d08b108fcb4e
md"""
```math
C_{D_i} =
2AR\left( \frac{\pi}{2}\right) \left[ A_1^2 + \sum_{n=2}^\infty nA_n^2\right]
```
"""

# ╔═╡ c2662e31-63af-4147-8c7a-096595925d81
md"""
Or alternatively
"""

# ╔═╡ 9c3c4d8e-01f8-4419-a380-24d95bfd4b0b
md"""
```math
C_{D_i} =
\pi AR (A_1^2) \left[ 1 + \sum_{n=2}^\infty n \left(\frac{A_n}{A_1} \right)^2\right]
```
"""

# ╔═╡ 724fdab4-bc9a-4118-b848-acc6eab6a2e7
md"""
Which can be written as
"""

# ╔═╡ 42c90cf9-b04d-4166-8757-965cbee75447
md"""
```math
C_{D_i} =
\frac{C_L^2}{\pi AR} \left( 1 + \delta\right)
```
"""

# ╔═╡ 45724bb4-08fe-4c75-8b7e-40bfbb1ba704
md"""
This results is interesting. It highlights the link between lift generation and the emergence of a drag as a result. In the equation above ``\delta`` is given as
"""

# ╔═╡ 8ec13686-818d-4030-a222-ee9d22e023e8
md"""
$$\delta =  \sum_{n=2}^\infty n \left(\frac{A_n}{A_1} \right)^2$$
"""

# ╔═╡ 851e21a9-ba3e-423a-ac28-a6b434e8649c
md"""
Notice that ``delta`` will always be positive because the ratios of Fourier coefficients is squared. Generally, the induced drag equation is written in terms of the wing *span efficiency factor*, ``e``, such that
"""

# ╔═╡ f2078958-27c3-4b11-8265-2443ccdec2d4
md"""
$$C_{D_i} = \frac{C_L^2}{\pi e AR}$$
"""

# ╔═╡ abe2ac68-a282-4a38-9733-213eb52b3a6c
md"""
Where the *span efficiency factor* is defined as
"""

# ╔═╡ d650a4ba-58f6-4b48-943f-a7dc46ecc553
md"""
$$e = \frac{1}{(1 + \delta)}$$
"""

# ╔═╡ 31763a87-f25f-463c-bb9a-5bf8992298c8
md"""
Exploring the implications of the induced drag equation shown in this section. It can be deduced that, in order to reduce the induced drag, it is desirable to increase the wing's aspect ratio, ``AR``. The reader is encouraged to revisit the definition of the aspect ratio to further explore the implications. What design variables would influence induce drag generation?
"""

# ╔═╡ c523c738-1e5b-42ab-9c7e-de51cfdddda0
md"""
### The special case of elliptical loading
"""

# ╔═╡ 6c8151ed-fe79-4faf-848b-708dbad1b171
md"""
Whilst the details will not be included here, it is a worthwhile exercise to explore the results shown above and to demonstrate that, in the special case of an elliptical distribution of circulation along the wing span, the induced drag is a minimum. It can also be shown that the induced velocity (or downwash) is constant. For convenience, the relationships for elliptical wing loading are given below.
"""

# ╔═╡ 34d51d30-f395-4843-a8cd-f006f1f62db5
md"""
```math
\Gamma(y) = \Gamma_0 \left[ 1 - \left(\frac{2y}{b}\right)^2 \right]^\frac{1}{2}
```
"""

# ╔═╡ 746057fb-28ff-4161-83eb-e040caa0ccd1
md"""
```math
\frac{d\Gamma(y)}{dy} = -\frac{4\Gamma_0}{b}\frac{(\frac{2y}{b})}{ \left[ 1 - \left(\frac{2y}{b}\right)^2 \right]^\frac{1}{2}}
```
""";

# ╔═╡ 96fee0b0-40da-4b7f-9bc2-a52393db2e8c
let
	b = 40
	b2 = b/2
	gamma0 = 100
	y = [-b2:b2/200:b2;]
	gamma(y) = gamma0*(1-(2y/b)^2)^(0.5)
	plot(y, gamma.(y), title="Elliptical distribution of circulation",legend=false; styling...)
	xlabel!("Wing span [m]")
	ylabel!("Circulation, Γ [m²/s]")
end

# ╔═╡ e4c009b8-92cc-46f2-b7f0-9638a4f42991
md"""
### Numerical solution procedure
"""

# ╔═╡ 91648402-73d3-4c5f-852b-4376809cdcb7
md"""
Now that the fundamental wing equation has been modified in terms of a Fourier series, the question that remains is how to use this result to find an estimation for the distribution of circulation over an arbitrary wing geometry. 
"""

# ╔═╡ cafe001c-e33b-4f3a-8bb8-c8f24e63db73
md"""
The solution requires access to all the geometric characteristics of the wing as a function of the spanwise coordinate ``y``. In particular,

* The geometric twist, ``\alpha(y)``
* The aerodynamic twist, ``\alpha_0(y)``
* The lift slope, ``m_0(y)`` (although, commonly it is assumed to be constant)
* Taper ratio, ``\lambda``
* Wing geometry (root and tip chord lengths, and wing span)
"""

# ╔═╡ 39e89cc0-8b1a-4b99-8c97-1187451cbe4d
md"""
The solution involves the evaluation of equation 5 (the transformed fundamental wing equation) at discrete points along the span. This will result in a linear systems of equations with unkown expansion coefficients ``A_n``, where ``n`` is the number of stations evaluated and only odd integers are required. That is, ``n_i = 2i-1`` where ``i = 1, 2, 3 \dots m``. 
"""

# ╔═╡ cfe864ae-4372-4c13-927c-bc909d33cd4d
md"""
It is convenient to rewrite the fundamental wing equation in matrix form, such that
"""

# ╔═╡ 666623e2-fa52-4871-bf00-a40b75d7b08b
md"""
```math
C(\theta(y), n) = 
\left( 
\frac{4b}{m_0(\theta(y))c(\theta(y))}   
+
\frac{n}{\sin \theta(y)}
\right) \sin n\theta(y)
```
"""

# ╔═╡ 86128bbd-7ee6-445a-a52a-8dd0db6384be
md"""
$$A(n) = A_n$$
"""

# ╔═╡ ffdc3be4-e2e9-4fae-995b-cfb764abb2c8
md"""
$$D(\theta(y)) = \alpha(\theta(y)) - \alpha_0(\theta(y))$$
"""

# ╔═╡ a0595f1d-b9da-40b0-9248-2faded70cea0
md"""
Resulting in the system
"""

# ╔═╡ f7779cd6-6a40-43d6-806b-69059da906c5
md"""
$$\sum_{n=1}^N C(\theta(y),n)A(n) = D(\theta(y))$$
"""

# ╔═╡ 936c6242-3313-44f0-be44-4da19783e443
md"""
The above is a system of linear equations that can then be solved using any familiar method for linear systems, such as the Guassian Elimination Method.
"""

# ╔═╡ 0f1c1126-5b51-4a82-8450-b8082530aaaf
md"""
### Computational solution (example)
"""

# ╔═╡ a4dbac35-2973-454e-9005-7b7068fb16d0
begin
	# Input information
	velocity 	= @bind V NumberField(1:100, default=20.0)
	f_m0 		= @bind m0 NumberField(0.01:0.01:3π, default=round(2π, digits=2))
	max_a 		= round(20, digits=2)
	max_a4 		= round(20/4, digits=2)
	max_a2 		= round(20/2, digits=2)
	s_alpha0 	= @bind alpha0 NumberField(0.0:max_a2/100:max_a2, default=max_a4)
	s_alpha 	= @bind alpha NumberField(0.0:max_a/100:max_a, default=max_a2)
	f_chord 	= @bind chord NumberField(0.5:0.5:10, default=5)
	f_b 	 	= @bind b NumberField(1:100, default=40)
	f_nPoints 	= @bind nPoints NumberField(1:50, default=4)
md"""
 $(velocity) *Velocity* [m/s] $(html"<br>")
 $(f_m0) *Lift slope* [/rad] $(html"<br>")
 $(s_alpha0) *Zero-lift angle* [deg] $(html"<br>")
 $(s_alpha) *Geometric angle* [deg] $(html"<br>")
 $(f_chord) *Chord* [m] $(html"<br>")
 $(f_b) *Span* [m] $(html"<br>")
 $(f_nPoints) *Points* $(html"<br>")
"""
end

# ╔═╡ 52e09fd6-beb4-48d4-98c0-7cafd9b6e230
let
	# Function definitions
	theta(y, b) = acos(2*abs(y)/b) 

	function c(y, chordR, chordT, b)
	    m = (chordT - chordR)/(b/2)
	    cy = y*m + chordR
	    return cy
	end
	
	function a(y, alphaR, alphaT, b)
	    m = (alphaT - alphaR)/(b/2)
	    ay = y*m + alphaR
	    return ay
	end
	
	function a0(y, alpha0R, alpha0T, b)
	    m = (alpha0T - alpha0R)/(b/2)
	    a0y = y*m + alpha0R
	    return a0y
	end
	
	function expansion_coeff(y, N, theta, c, m)
	    t1 = (4*b)/(m(y)*c(y))
	    t2 = N/sin(theta(y))
	    t3 = t1 + t2
	    return t3*sin(N*theta(y))
	end
	
	function Gamma0(A, nPoints::Int, y; theta=theta)
	    G = zeros(nPoints)
	    for i ∈ 1:nPoints
	        for n ∈ 1:nPoints
	            G[i] += A[n]*sin((2*n-1)*theta(y[i]))
	        end
	    end
	    return G
	end
	
	# Function redefinitions
	theta(y) = theta(y, b)
	c(y) = chord #c(y, chordR, chordT, b)
	a(y) = alpha #a(y, alphaR, alphaT, b)
	m(y) = m0
	
	
	
	# Main solution logic
	bh = b/2
	DY = b/2/nPoints
	y = [0+DY/2:DY:b/2;]
	
	C = zeros(nPoints, nPoints)
	D = zeros(nPoints)
	
	for j ∈ 1:nPoints
	    for i ∈ 1:nPoints
	        C[i,j] = expansion_coeff(y[i], 2j-1, theta, c, m)
	    end
	    D[j] = deg2rad(a(y[j]) - alpha0)
	end

	# Solve system
	A = similar(D)
	try
		A = C\D
	catch
		print("NOT DOING THIS!")
	end
	

	# Calculate circulation and C_L from coefficients
	G = Gamma0(A, nPoints, y)
	CL = round(π*A[1]*b/chord, digits=6)
	
	# Plot result
	gammaDistribution = plot(y, G, label="Circulation",
		title="Lift coefficient = $(CL)", 
		xlims=(0.0, b/2); styling...)
	xlabel!("Wing semi-span [m]")
	ylabel!("Circulation [m²/s]")
	scatter!(y, 0.0*y, color=:red, label="Grid points")
	if alpha0 >= alpha
		yLoc = (ylims(gammaDistribution)[1] + ylims(gammaDistribution)[2])/2
		annotate!(gammaDistribution, b/4, yLoc, text("Check angles!", :red))
	end
	# savefig(gammaDistribution, "fig_solution_example_20.svg")
	gammaDistribution
end

# ╔═╡ 2e464a75-4241-497c-8ac8-77a86eb389cc
md"""
!!! warning
	Some combinations of variables can cause the solution to become unbounded. In most cases this can be handled by Julia. However, in some cases it will result in a memory overflow and the kernel will become unresponsive. This is not a reason for concern, but if it happens just relaunch the notebook. I will include some numerical checks in the future to prevent this from happening. As an educational opportunity... from the theory we have explored, you should know why this is the case, can you think of the reason?
"""

# ╔═╡ 22d47c0e-566c-4766-bc00-2f63a7d85a3b
md"""
### Final observations
"""

# ╔═╡ ba9eeaff-3c45-4838-8c9d-8cc98dc84250
md"""
A method to calculate the aerodynamic performance of a finite wing has been explored, the lifting line theory. This method requires that basic information about a wing plantform is available, from which the 3D aerodynamic behaviour can be estimated. It surprising that despite the assumptions made, an the ingenious manipulations required to arrived at a "solvable" form of the general wing equation, the lifting line theory can return resonably accurate results. A major limitation of the lifting line theory is that in only consider inviscid flow. A more complete model of the drag of an aircraft wing can be formulated by including viscous drag (the topic of the next lectures in this lecture series). For example, the drag polar model is given as
"""

# ╔═╡ 15e79076-22f0-4996-84ab-9fda367db05d
md"""
$$C_D = C_{D_0} + C_{D_i}$$
"""

# ╔═╡ 598b953f-6656-4c73-892b-1432f2102581
md"""
Where, ``C_D`` is the total drag.``C_{D_0}`` and ``C_{D_i}`` are the *viscous* and *induced* drag, respectively. Subtituting the expression for induced drag (shown earlier)
"""

# ╔═╡ 9ffc04eb-ff4a-4b3f-9dad-fd969a77950f
md"""
$$C_D = C_{D_0} + \frac{C_L^2}{\pi e AR}$$
"""

# ╔═╡ 60f4e04b-b53d-4b05-b0ef-21bcde6f7ebb
md"""
This is a useful result that aerodynamicists know well. Lift generation comes at a cost (``\propto C_L^2``), but by increasing the aspect ratio of a wing and choosing an optimal platfrom the induced drag component can be reduced. The viscous drag ``C_{D_0}`` is, in fact, rather difficult to estimate theoretically and it is often necessary to carry out experiments or used correlations i.e. a certain level of empirism is generally required. This topic will be explored in the next lectures.
"""

# ╔═╡ aecce740-7ef7-4311-b7c1-467c6a21cec1
md"""
## Bibliography
"""

# ╔═╡ c45da147-72b5-4347-b2e0-aa6f7ef42ca3
md"""
$(bibItem("E.L. Houghton, P.W. Carpenter, Steven H. Collicott, Daniel T. Valentine",2015,"Aerodynamics for Engineering Students, Seventh Edition","Elsevier"))
$(bibItem("Flandro, Gary A. McMahon, Howard M. Roach, Robert L.", 2012, "Basic Aerodynamics - Incompressible Flow", "Cambridge University Press"))
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Downloads = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
LaTeXStrings = "~1.3.0"
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
# ╟─e3b13a64-4c36-46b3-8532-575a71ad3dee
# ╟─81548881-1de0-4c19-898d-4b85cb71cff1
# ╟─0267d4be-7a71-4b21-a7a3-48e0bb86bb64
# ╟─620122a7-7623-43ba-9ee1-a765786d1e93
# ╟─b25aff1f-d0f4-4003-95fc-f512f4113f2b
# ╟─fe2a2cd8-99ba-4b5f-ae0c-c0f096002108
# ╟─c0bd19e4-a932-41c3-8a7e-45a326593e27
# ╟─98576fe8-cd61-4d8b-9628-1ff8e83ccffb
# ╟─8d6f9c8e-5b66-4d42-8836-a4e484dd7cbf
# ╟─f80f03cb-f5cb-40eb-9483-48d7ac6b33bc
# ╟─d13bdb89-4c57-4a08-b764-698b8bece8b8
# ╟─9a637330-680c-43db-881f-cf2198913640
# ╟─0215c2f4-892f-40c4-acdd-f0e90d850c41
# ╟─1d221e4c-a532-4bcd-9b39-9e708daf0aa7
# ╟─d946d4ed-847e-4f1f-ae27-357435aa68b6
# ╟─56f1ea76-f726-4327-87b2-b43192874401
# ╟─8389b67f-86ee-4f86-90da-90ccf4dc450a
# ╟─67b75aa2-be8b-49f3-be61-bb5e75477b42
# ╟─55e66d16-8759-49e5-8c18-de0d41319f48
# ╟─ef93499e-d302-45bd-b3d1-43f4a768ae41
# ╟─1f5cda8a-d26f-4d5c-b1a5-29f3c9cf2ff7
# ╟─7ce50ac3-de76-4044-a3bc-d9046d60ac53
# ╟─e1d220ba-4e19-42fd-9797-67b3fecde01d
# ╟─f86f3db6-a93c-45f6-8250-6b77c2cd3157
# ╟─e710ee7d-1110-4a91-b4a2-073afa13e218
# ╟─acba8da7-e635-416f-9105-892344953e35
# ╟─5d62a58f-2c10-4997-a0ef-87f8b99a026a
# ╟─eeed306f-c59b-4586-9e7f-70caf49e39f9
# ╟─0d085718-6bfb-4e0d-a08b-acab43870a4b
# ╟─f47d77a9-652d-4d54-a12c-a7252fb6cfe5
# ╟─6c77ad6b-606a-4050-b136-312f667b84a1
# ╟─c0930aa3-b436-4ec4-bdba-62f044298a2f
# ╟─06aafb04-1347-4653-abe6-6c6f0c2fbf93
# ╟─ba451210-b849-424d-ac77-37a0ef98ee1a
# ╟─0212b76d-2304-45b1-b23f-0228e8ab39b8
# ╟─4164f7fa-96df-4cf5-bd79-4fcff7a5f660
# ╟─531b2b1a-903b-4cab-b4dd-82827a6721be
# ╟─e53dbc90-4840-44e1-8e4b-7c9835f10fb2
# ╟─0405524a-7915-41e4-a264-d545560f44a3
# ╟─e806a0a9-dcdf-40f3-8ba6-870866c979ed
# ╟─751648da-3723-44ae-8a19-683437969922
# ╟─58026876-f278-4211-ae68-042fab586bdc
# ╟─e64bb91c-4105-479f-bb82-e2ddff6eda23
# ╟─436845d8-72fb-49f7-abaa-45e141a8b3b5
# ╟─31d527d1-03ad-4ffc-83a0-59b3ee9f8d84
# ╟─81127462-87c7-4c58-8ef8-097101f2fc5d
# ╟─d63c8e12-def8-4d56-85e9-8457a4a79c09
# ╟─f5cdc961-feeb-478a-8b2a-9b45c649fbe2
# ╟─c89ce8fb-0ac3-4360-b7c2-0fb6a05fcb93
# ╟─95be0392-e4ae-42c2-910b-e7950afb3f51
# ╟─7cb6536a-66e5-4156-bfbc-3fe3d61f9966
# ╟─4cea8ab8-05d9-4509-acec-e22442aa21ae
# ╟─ee4ce3ca-5870-4292-864d-24ddfcbaa249
# ╟─af702e5a-2692-47e4-aee7-a334f8c163fe
# ╟─7af5aa5a-6412-4dac-835b-c21475fc04b3
# ╟─81de1481-3c19-4112-aa8e-ea26231dbdd3
# ╟─3cd73b58-7045-42a4-bbe2-e9213d6b0432
# ╟─6bdff171-e866-4059-a0f5-b6863e7b1729
# ╟─f92b9592-51e3-45b1-a497-d08b108fcb4e
# ╟─c2662e31-63af-4147-8c7a-096595925d81
# ╟─9c3c4d8e-01f8-4419-a380-24d95bfd4b0b
# ╟─724fdab4-bc9a-4118-b848-acc6eab6a2e7
# ╟─42c90cf9-b04d-4166-8757-965cbee75447
# ╟─45724bb4-08fe-4c75-8b7e-40bfbb1ba704
# ╟─8ec13686-818d-4030-a222-ee9d22e023e8
# ╟─851e21a9-ba3e-423a-ac28-a6b434e8649c
# ╟─f2078958-27c3-4b11-8265-2443ccdec2d4
# ╟─abe2ac68-a282-4a38-9733-213eb52b3a6c
# ╟─d650a4ba-58f6-4b48-943f-a7dc46ecc553
# ╟─31763a87-f25f-463c-bb9a-5bf8992298c8
# ╟─c523c738-1e5b-42ab-9c7e-de51cfdddda0
# ╟─6c8151ed-fe79-4faf-848b-708dbad1b171
# ╟─34d51d30-f395-4843-a8cd-f006f1f62db5
# ╟─746057fb-28ff-4161-83eb-e040caa0ccd1
# ╟─96fee0b0-40da-4b7f-9bc2-a52393db2e8c
# ╟─e4c009b8-92cc-46f2-b7f0-9638a4f42991
# ╟─91648402-73d3-4c5f-852b-4376809cdcb7
# ╟─cafe001c-e33b-4f3a-8bb8-c8f24e63db73
# ╟─39e89cc0-8b1a-4b99-8c97-1187451cbe4d
# ╟─cfe864ae-4372-4c13-927c-bc909d33cd4d
# ╟─666623e2-fa52-4871-bf00-a40b75d7b08b
# ╟─86128bbd-7ee6-445a-a52a-8dd0db6384be
# ╟─ffdc3be4-e2e9-4fae-995b-cfb764abb2c8
# ╟─a0595f1d-b9da-40b0-9248-2faded70cea0
# ╟─f7779cd6-6a40-43d6-806b-69059da906c5
# ╟─936c6242-3313-44f0-be44-4da19783e443
# ╟─0f1c1126-5b51-4a82-8450-b8082530aaaf
# ╟─a4dbac35-2973-454e-9005-7b7068fb16d0
# ╟─52e09fd6-beb4-48d4-98c0-7cafd9b6e230
# ╟─2e464a75-4241-497c-8ac8-77a86eb389cc
# ╟─22d47c0e-566c-4766-bc00-2f63a7d85a3b
# ╟─ba9eeaff-3c45-4838-8c9d-8cc98dc84250
# ╟─15e79076-22f0-4996-84ab-9fda367db05d
# ╟─598b953f-6656-4c73-892b-1432f2102581
# ╟─9ffc04eb-ff4a-4b3f-9dad-fd969a77950f
# ╟─60f4e04b-b53d-4b05-b0ef-21bcde6f7ebb
# ╟─aecce740-7ef7-4311-b7c1-467c6a21cec1
# ╟─c45da147-72b5-4347-b2e0-aa6f7ef42ca3
# ╟─2500b290-30fb-11ec-2224-ed2674139124
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
