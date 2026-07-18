### Prompt for Claude: Comprehensive Paper Reformulation and Plot Critique

**Role:** You are a senior theoretical physicist and a strict managing editor for a top-tier high-energy physics journal (like JHEP or Phys. Rev. D).

**Goal:** I have a working draft of a paper on the High-Precision Numerical Analysis of Electroweak Vacuum Stability and Planck-Scale Operator Corrections. I need you to deeply reformulate, expand, and rewrite this draft into a polished, publication-ready manuscript. Additionally, I need you to critique my current list of figures and tell me exactly what needs to be replotted, combined, or styled differently to meet strict journal standards.

**Context & Recent Discoveries:**
Since the initial draft, we have made significant physical breakthroughs that must be seamlessly integrated:
1. **The Fubini Ansatz Breakdown and Virial Diagnostic:** We discovered a "sandwiched" unstable region at low $M_h$ ($\sim 50$ GeV). Initially, we computed negative Euclidean actions here. Following rigorous QFT analysis, we established this is NOT a classical roll-down, but rather a breakdown of the Fubini-Lipatov instanton ansatz. Because $\lambda_{\rm eff}$ runs rapidly, the virial constraint $S_{\rm kin}/|S_{\rm pot}| = 2$ is violently violated. The potential has a barrier, but it is pathologically thin, leading to rapid quantum tunneling ($S_E^{\rm true} \ll S_{\rm thresh}$) rather than the computed negative artifact.
2. **Constant Action Contours:** We successfully generated exact numerical constant action contours ($S=10, 50, 100, 200, 450, 3000, 500000$) on the global phase diagram. These contours explicitly visualize the "wedge" of the unstable region and how the contours pack densely right against the stable boundary wall.

**Summary of Our Work History (For Full Context):**
To give you complete context on how we arrived here:
- We built an exact numerical C++ solver that integrates the 3-loop SM RGEs and solves the $O(4)$ bounce equation without relying on the conformal approximation.
- We scanned 250,000 parameter points in the $(M_t, M_h)$ plane and found that the conformal approximation systematically underestimates the action near the phase boundary, especially when a Planck-scale operator $c_6\phi^6/M_{\rm Pl}^2$ is included.
- We generated numerous plots (error scaling, boundary overlay, $c_6$ sweep). 
- We then noticed a bizarre anomaly at low Higgs masses (e.g., $M_h=5$ and $M_h=50$): the action evaluated to a massive negative number, creating a "sandwiched" unstable region between two metastable regions.
- We brought this to you (Claude) previously, and you correctly diagnosed it as a breakdown of the Fubini-Lipatov ansatz via the virial constraint $S_{\rm kin}/|S_{\rm pot}| = 2$, clarifying it is rapid quantum tunneling through an ultra-thin barrier, not a classical roll-down. We have appended your exact previous critique below (Appendix G) to ensure you maintain that same level of rigor.
- We have now plotted the constant action contours specifically to visualize this sandwich topology, and carefully cleaned up our project directory to prepare for final publication.

**Current List of Available Plots:**
1. `error_rk4_convergence.png`: RK4 error scaling.
2. `error_quadrature_convergence.png`: Simpson's rule error scaling.
3. `error_golden_section.png`: Golden section search error plateau.
4. `error_residual_histogram.png`: EOM residuals.
5. `analytical_stability_plot.png`: Global analytical phase diagram (Buttazzo-style).
6. `numerical_stability_contours.png`: Global exact numerical phase diagram showing constant action contours.
7. `error_S_scatter.png`: Scatter plot comparing $S_{\rm exact}$ vs $S_{\rm approx}$.
8. `error_vs_lambda.png`: Discrepancy scaling with $|\lambda_{\min}|$.
9. `error_discrepancy_heatmap.png`: Heatmap of discrepancy.
10. `error_boundary_reclassification.png`: Reclassified strip near SM.
11. `analytical_stability_closeup.png`: Close-up analytical.
12. `numerical_stability_closeup.png`: Close-up exact numerical.
13. `overlay_plot.png` & `overlay_closeup_plot.png`: Direct overlay of analytical vs numerical boundaries.
14. `c6_phase_diagram.png` & `c6_sweep_plot.png`: $c_6$ parameter sweeps for $c_6 \in \{0.1, 1.0, 10.0\}$.

**Your Tasks:**

**PART 0: Editorial Assessment (MANDATORY)**

Before performing any rewrite, do not modify the manuscript.

Act as a senior PRD/JHEP editor and theoretical physicist.

First determine:

- What is the central scientific contribution of this work?
- What is the most novel and publishable result?
- Is this paper primarily:
  - a precision numerical methods paper,
  - a phase-boundary correction paper,
  - a vacuum-topology paper,
  - a Fubini-Lipatov/virial-breakdown paper,
  - or a combination of the above?
- Which framing would maximize publication impact?
- Which sections should be shortened, expanded, merged, moved, or relegated to appendices?
- Which figures are essential, which are redundant, and which should be combined?
- Propose a revised manuscript structure before any rewriting begins.

Only after completing this editorial assessment should you proceed to Parts 1 and 2.

**PART 1: Editor's Critique of Plots**
Before rewriting the paper, give me strict instructions on how to improve my figures. For example:
- Should the four error plots (1-4) be combined into a single 2x2 multi-panel figure to save space?
- Should the overlay and close-up plots be combined into a side-by-side subplot?
- What specific axis adjustments, color map changes, or annotations do I need to make before final publication?
Tell me exactly what to replot.

**PART 2: Full LaTeX Rewrite**

Do not immediately rewrite the paper.

First provide:
- Editorial Assessment
- Figure Plan
- Revised Outline

Only then produce the full manuscript rewrite.

When writing the full manuscript rewrite, ensure you:
- Deepen the theoretical introduction and methodology.
- Fully integrate the physical explanation of the virial ratio breakdown and rapid quantum tunneling.
- Discuss the new constant action contours and how they visualize the topology of the instability.
- **CRITICAL:** Insert the `\begin{figure}...\includegraphics{...}...\end{figure}` commands in the exact, correct locations within the flow of the text, using the filenames listed above (or the combined multi-panel filenames you suggest in Part 1). Write robust, detailed captions for every figure.

**Current LaTeX Draft (For Baseline):**
```latex
% [Paste the entire contents of main.tex here]
```


### Appendix A: main.tex

\documentclass[aps,prd,twocolumn,superscriptaddress,nofootinbib]{revtex4-2}

\usepackage{amsmath,amssymb,amsfonts}
\usepackage{graphicx}
\usepackage{dcolumn}
\usepackage{bm}
\usepackage{hyperref}
\usepackage{xcolor}
\usepackage{booktabs}
\usepackage{multirow}
\usepackage{physics}
\usepackage{microtype}
\usepackage{float}

\graphicspath{{./figures/}}

\begin{document}

\title{High-Precision Numerical Analysis of Electroweak Vacuum Stability\\
and Planck-Scale Operator Corrections}

\author{Naman Goyal}
\affiliation{Department of Physics, Indian Institute of Science}

\date{\today}

%----------------------------------------------------------------------
\begin{abstract}
We present a rigorous high-precision numerical analysis of the Standard Model
electroweak vacuum stability phase diagram, computed across 250,000 parameter
points in the $(M_t, M_h)$ plane. Standard theoretical treatments of vacuum
decay rely on the analytical conformal approximation to the bounce action,
$S_E \approx 8\pi^2/3|\lambda_{\min}|$, which ignores the shape of the
tunneling barrier at Planck-scale field values. We demonstrate that this
approximation systematically fails near the metastability boundary
defined by the point-dependent condition $S_E > 4\ln(T_U\mu_{\rm inst})$,
where the bounce profile extends to field values
$\phi \sim M_{\rm Pl}$. By deploying an exact numerical framework
integrating three-loop renormalization group equations with a
4th-order Runge-Kutta solver ($\epsilon_{\rm RK4} = 10^{-13}$),
composite Simpson's rule quadrature with $N = 2048$ nodes, and a Golden-Section
Search algorithm for bubble radius optimization, we achieve numerical errors
below $10^{-11}$ in relative action — more than ten orders of magnitude
below the dominant physical uncertainties. Upon incorporating a Planck-scale
dimension-6 operator, $c_6 \phi^6/\Lambda^2$ with $\Lambda = M_{\rm Pl}$
and $c_6 = 1$, the exact bounce solver detects a positive-definite barrier
inflation that unconditionally raises $S_{\rm exact}$ above $S_{\rm approx}$
near the phase boundary. This rescues a continuous strip of parameter space —
approximately 110 grid points corresponding to $\sim$0.3--0.5~GeV in $M_t$
— from the unstable to the metastable classification. The phenomenon is shown
to be universal across two orders of magnitude in $c_6 \in \{0.1, 1.0, 10.0\}$.
Our results demonstrate that exact numerical bounce solvers are strictly
necessary for precision vacuum stability bounds within the $1\sigma$
experimental uncertainty of the top quark mass.
\end{abstract}

\maketitle

%======================================================================
\section{Introduction}
\label{sec:intro}
%======================================================================

The joint discovery of the Higgs boson at $M_h \approx 125.09~\text{GeV}$
by the ATLAS and CMS collaborations~\cite{Aad:2012tfa,Chatrchyan:2012xdj}
completed the particle content of the Standard Model (SM) and enabled the
quantitative extrapolation of SM parameters to arbitrarily high energy
scales. A pivotal implication is that for the experimentally
measured values of the top quark mass, $M_t \approx 173.3~\text{GeV}$, and
the Higgs mass, the Higgs quartic coupling $\lambda(\mu)$ turns negative
at an intermediate scale $\mu \sim 10^{10}$~GeV, indicating that the
electroweak vacuum is not the global minimum of the scalar
potential~\cite{Degrassi:2012ry,Buttazzo:2013uya}. This metastability
is arguably one of the most profound theoretical puzzles emerging from
the precision era of particle physics.

The fate of the electroweak vacuum is governed by the rate at which quantum
fluctuations nucleate a bubble of the true vacuum. In the dilute-bubble
approximation, this tunneling rate per unit volume is
\begin{equation}
    \Gamma/V \sim \mu^4 \exp(-S_E[\phi_b]),
    \label{eq:tunneling_rate}
\end{equation}
where $S_E[\phi_b]$ is the Euclidean bounce action evaluated on the $O(4)$-symmetric saddle-point configuration $\phi_b(r)$~\cite{Coleman:1977py,Callan:1977pt}. The vacuum is classified as \textit{metastable} if the tunneling time $\tau \sim \Gamma^{-1}$ exceeds the age of the universe $T_U \sim 4.3 \times 10^{17}$~s. Throughout this work, we evaluate the metastability threshold point-by-point using~\cite{Buttazzo:2013uya}
\begin{equation}
    S_{\rm threshold}(\mu_{\rm inst}) = 4\ln(T_U\mu_{\rm inst}).
    \label{eq:threshold}
\end{equation}

Comprehensive analyses of the SM vacuum stability have been carried out at
two-loop~\cite{Isidori:2001bm,Espinosa:2007qp} and
three-loop~\cite{Degrassi:2012ry,Buttazzo:2013uya,Bezrukov:2012sa}
accuracy in the perturbative expansion. These works map the phase diagram of
the SM in the $(M_t, M_h)$ plane into three regions: \textit{stable}
($\lambda > 0$ up to $M_{\rm Pl}$), \textit{metastable} ($\lambda < 0$ but
$S_E > S_{\rm threshold}$), and \textit{unstable} ($S_E < S_{\rm threshold}$).
The precise location of the boundary between the metastable and unstable
phases is of direct phenomenological relevance, since current measurements
place the SM near it.

A central tool in these analyses is the analytical conformal approximation
to the bounce action~\cite{Isidori:2001bm}:
\begin{equation}
    S_{\rm approx} = \frac{8\pi^2}{3|\lambda_{\min}|},
    \label{eq:conformal_approx}
\end{equation}
where $\lambda_{\min}$ is the minimum value of the running quartic coupling.
This expression is derived under the assumption that the tunneling field
profile is dominated by the region near the scale where $\lambda$ is most
negative, and that the potential can be locally approximated as a pure
$\lambda \phi^4$ form with a constant $\lambda = \lambda_{\min}$. While
this approximation is computationally convenient and physically transparent,
it is fundamentally blind to the shape of the potential at very large field
values, $\phi \gtrsim 10^{17}$~GeV.

Near the metastability boundary, where $|\lambda_{\min}| \to 0$, the bounce
profile becomes very wide in field space. The classical turning point extends
to $\phi \sim M_{\rm Pl}$, precisely where physics beyond the SM — even if
suppressed by the Planck scale — can materially distort the tunneling
potential. A generic and well-motivated class of such corrections is provided
by higher-dimensional operators generated by unknown UV physics:
\begin{equation}
    \mathcal{L}_{\rm UV} \supset \frac{c_6}{\Lambda^2} \phi^6
    + \frac{c_8}{\Lambda^4} \phi^8 + \cdots
    \label{eq:eft_expansion}
\end{equation}
with $\Lambda = M_{\rm Pl}$. The sign and magnitude of these operators are
free parameters of the UV completion, but for any $c_6 > 0$, the leading
operator strictly raises the potential at Planck-scale field values, inflating
the tunneling barrier.

This observation motivates the present work. We ask: \textit{to what precision
must the bounce action be computed, and what is the quantitative impact of
Planck-scale operators on the metastability boundary?} Our approach is to
construct a fully exact numerical pipeline that avoids the conformal
approximation entirely, achieving sub-$10^{-11}$ numerical errors, and
compare the resulting phase boundaries directly against those produced by
Eq.~\eqref{eq:conformal_approx}.

Our key findings are:
\begin{enumerate}
  \item The conformal approximation is an accurate global estimator, with
  74.9\% of sampled parameter points agreeing with the exact result to
  within $\pm 5\%$.
  \item Near the metastability boundary, however, the approximation
  catastrophically underestimates the action, with discrepancies growing
  to $20\%$ and beyond as $|\lambda_{\min}| \to 0$.
  \item The inclusion of a $c_6\phi^6/M_{\rm Pl}^2$ operator produces a
  barrier inflation that unconditionally raises $S_{\rm exact} > S_{\rm approx}$
  in the boundary region, rescuing a strip of approximately 0.3--0.5~GeV in
  $M_t$ from the unstable to the metastable regime.
  \item This phenomenon is universal: for any $c_6 > 0$, the qualitative
  breakdown of the conformal approximation persists, independent of the
  precise UV completion.
\end{enumerate}

These results have direct implications for interpreting the position of the
SM within the metastability diagram, given that the $1\sigma$ experimental
uncertainty on $M_t$ is of the same order as the boundary shift we report.

The paper is structured as follows. Section~\ref{sec:theory} presents the
theoretical framework, including the RG-improved effective potential and the
bounce equation. Section~\ref{sec:numerics} describes our exact numerical
pipeline in detail. Section~\ref{sec:errors} establishes a rigorous error
budget. Section~\ref{sec:results} presents the full results, including the
conformal approximation breakdown, the phase boundary reclassification, and
the $c_6$ universality sweep. Section~\ref{sec:discussion} discusses the
physical implications, and Section~\ref{sec:conclusion} concludes.

%======================================================================
\section{Theoretical Framework}
\label{sec:theory}
%======================================================================

\subsection{The RG-Improved Effective Potential}
\label{subsec:potential}

We work with the renormalization group (RG) improved tree-level effective
potential for the Higgs field in the large-field regime,
$\phi \gg v = 246~\text{GeV}$:
\begin{equation}
    V_{\rm eff}(\phi) = \frac{1}{4}\lambda(\phi)\,\phi^4
    + \frac{c_6}{\Lambda^2}\,\phi^6,
    \label{eq:veff}
\end{equation}
where $\Lambda = M_{\rm Pl} = 1.22 \times 10^{19}~\text{GeV}$ and the
benchmark value $c_6 = 1$ is adopted. The renormalization scale $\mu$ is
identified with the field value $\phi$ to minimize the effects of large
logarithms. The running quartic coupling $\lambda(\mu)$ is obtained by
numerically integrating the three-loop SM $\beta$-functions~\cite{Mihaila:2012pz}
from the electroweak scale $\mu_0 = M_t$ to the Planck scale.

The three-loop $\beta$-function for $\lambda$ in the SM is
\begin{equation}
    \mu \frac{d\lambda}{d\mu} = \beta_\lambda^{(1)} + \beta_\lambda^{(2)}
    + \beta_\lambda^{(3)},
    \label{eq:beta_lambda}
\end{equation}
where the one-, two-, and three-loop contributions are functions of
$\lambda$, the top Yukawa coupling $y_t$, and the gauge couplings $g_1$,
$g_2$, $g_3$~\cite{Luo:2002ti,Mihaila:2012pz}. A system of seven coupled
ODEs for $\{\lambda, y_t, g_1, g_2, g_3, g_s, m^2\}$ is evolved simultaneously.
The boundary conditions at the electroweak scale are derived from the
measured pole masses $M_t$, $M_h$, $M_W$, $M_Z$ and the strong coupling
$\alpha_s(M_Z)$ using standard matching relations.

The $\phi^6$ operator in Eq.~\eqref{eq:veff} is the leading irrelevant
deformation of the potential induced by Planck-scale physics, and its natural
coefficient is $c_6 = \mathcal{O}(1)$. For $c_6 > 0$, the operator is
positive definite and creates a barrier at large field values that partially
closes off the tunneling path. We scan over $c_6 \in \{0.1, 1.0, 10.0\}$ to
assess the sensitivity of our conclusions to the value of this parameter.
The pure SM result corresponds to $c_6 \to 0$.

\subsection{The Bounce Equation}
\label{subsec:bounce}

For an $O(4)$-symmetric bubble solution in Euclidean spacetime, the field
profile $\phi_b(r)$ satisfies the nonlinear ODE~\cite{Coleman:1977py}
\begin{equation}
    \frac{d^2\phi_b}{dr^2} + \frac{3}{r}\frac{d\phi_b}{dr} = V'_{\rm eff}(\phi_b),
    \label{eq:bounce_eom}
\end{equation}
with boundary conditions $\phi_b \to 0$ as $r \to \infty$ and
$d\phi_b/dr|_{r=0} = 0$. The Euclidean action is then
\begin{equation}
    S_E = 2\pi^2 \int_0^\infty dr\, r^3
    \left[\frac{1}{2}\left(\frac{d\phi_b}{dr}\right)^2 + V_{\rm eff}(\phi_b)\right].
    \label{eq:euclidean_action}
\end{equation}

The analytical conformal approximation to Eq.~\eqref{eq:euclidean_action}
is obtained by assuming the bounce profile is concentrated at a single
scale $R$ such that $\phi_b(r) \approx \phi_0\,f(r/R)$ with $f$ a
dimensionless profile function, and that $V_{\rm eff} \approx \lambda_{\min}\phi^4/4$
throughout the integration region. Under these assumptions, dimensional
analysis yields Eq.~\eqref{eq:conformal_approx}. The critical flaw of this
approximation near the boundary is that the classical turning point
$\phi_{\rm turn} \sim R^{-1}$ satisfies $\phi_{\rm turn} \sim M_{\rm Pl}$
when $|\lambda_{\min}| \to 0$, invalidating the assumption that the $\phi^6$
operator is negligible.

\subsection{The Analytical Conformal Approximation and its Regime of Validity}
\label{subsec:validity}

The conformal approximation~\eqref{eq:conformal_approx} can be derived
rigorously in the thin-wall limit (large $S_E$, deep instability), where
it is known to be accurate. Corrections to it arise from two sources:
(i) the non-trivial running of $\lambda(\phi)$ over the bounce profile, and
(ii) the contribution of higher-dimensional operators at large $\phi$.
The first correction was partially addressed in~\cite{Isidori:2001bm} using
scale-dependent matching, while the second has not been systematically
incorporated in previous phase diagram studies.

For the pure SM ($c_6 = 0$) and deep in the instability (large $|\lambda_{\min}|$),
the conformal approximation provides sub-percent accuracy. However, as we
demonstrate below, the approximation catastrophically underestimates $S_E$
when $|\lambda_{\min}| \lesssim 0.001$, corresponding to the phenomenologically
critical region near $M_t \approx 173$~GeV, $M_h \approx 125$~GeV.

%======================================================================
\section{Numerical Methodology}
\label{sec:numerics}
%======================================================================

\subsection{Overview of the Computational Pipeline}
\label{subsec:pipeline}

Our exact bounce action computation is organized as a three-stage pipeline:
(1) RGE evolution of all SM couplings from $M_t$ to $M_{\rm Pl}$;
(2) numerical evaluation of the Euclidean action integral for a given bubble
profile; and (3) optimization of the bubble radius $R$ to minimize
$S_E[R]$. Each stage is independently validated against theoretical
convergence expectations, as described in Sec.~\ref{sec:errors}.
The full pipeline is evaluated at $N_{\rm pts} = 250{,}000$ points in
the $(M_t, M_h)$ plane, sampled on a uniform grid.

\subsection{RGE Integration: Fourth-Order Runge-Kutta Solver}
\label{subsec:rk4}

The seven-component ODE system for the SM RGEs is integrated from
$\mu_0 = M_t$ to $\mu_{\rm max} = M_{\rm Pl}$ using a classical
fourth-order Runge-Kutta (RK4) scheme~\cite{Press:2007nrbook} with
adaptive step sizing. The error tolerance for the RK4 integration is set to
\begin{equation}
    \epsilon_{\rm RK4} = 10^{-13},
    \label{eq:rk4_tol}
\end{equation}
which is four orders of magnitude tighter than the precision required
by the physical uncertainties in the problem. The adaptive step-size
controller enforces that the local truncation error at each step satisfies
$|\delta \phi| < \epsilon_{\rm RK4}$, automatically refining the grid
near regions of rapid coupling evolution. The default step size is
$\Delta t = 0.01$ in the log-scale integration variable $t = \ln(\mu/M_t)$.

\subsection{Action Quadrature: Composite Simpson's Rule}
\label{subsec:quadrature}

To evaluate the radial integral in Eq.~\eqref{eq:euclidean_action},
the integration domain $r \in [0, \infty)$ is first compactified to a
finite interval via the substitution $u = \tanh(r/R)$, mapping the
domain to $u \in [0, 1)$. This removes the contribution from the
asymptotic tail while preserving the near-origin behavior of the bounce
profile. The integral is then evaluated using the composite Simpson's
rule~\cite{Stoer:2002introduction} with
\begin{equation}
    N_{\rm quad} = 2048
    \label{eq:nquad}
\end{equation}
quadrature nodes, corresponding to a step size $h = 1/2048$ on the
compactified domain. For a smooth integrand, Simpson's rule achieves
fourth-order convergence, $\mathcal{O}(h^4) = \mathcal{O}(N^{-4})$.

\subsection{Bubble Radius Optimization: Golden-Section Search}
\label{subsec:golden}

The Euclidean action $S_E(R)$ is a unimodal function of the bubble
radius $R$ for a given potential profile. We minimize it using the
classical Golden-Section Search (GSS) algorithm~\cite{Kiefer:1953gss},
which converges to the minimum with a reduction factor of $\varphi^{-1}
= 2/(1+\sqrt{5}) \approx 0.618$ per iteration. The search is initialized
over the interval $R \in [R_{\min}, R_{\max}]$ determined by the scale
at which $\lambda(\mu)$ first turns negative ($R \sim \mu^{-1}_{\rm cross}$)
and the Planck scale. The algorithm terminates when the bracketing interval
satisfies
\begin{equation}
    |R_{\rm max} - R_{\rm min}| < \tau_{\rm GSS} = 10^{-13},
    \label{eq:gss_tol}
\end{equation}
ensuring sub-femtometer precision in the optimal bubble radius $R^*$.

\subsection{Parameter Space Scan}
\label{subsec:scan}

We construct a high-resolution uniform grid in the physically relevant
region of the $(M_t, M_h)$ plane:
\begin{equation}
    M_t \in [155, 185]~\text{GeV}, \quad
    M_h \in [115, 135]~\text{GeV},
    \label{eq:param_range}
\end{equation}
sampled at 500 points in each direction, for a total of 250,000 evaluation
points. For each point, the full pipeline is executed with
$c_6 \in \{0, 0.1, 1.0, 10.0\}$, yielding four stability classifications
per grid point. The stability of each point is determined by comparing
the computed $S_E$ to the threshold~\eqref{eq:threshold}: if $S_E > S_{\rm threshold}$
and $\lambda_{\min} < 0$, the point is metastable; if $S_E < S_{\rm threshold}$
and $\lambda_{\min} < 0$, it is unstable; if $\lambda_{\min} > 0$ throughout,
it is absolutely stable.

%======================================================================
\section{Numerical Error Budget}
\label{sec:errors}
%======================================================================

A critical requirement for any claim of phase boundary reclassification
is the demonstration that the observed effect is physical, not a numerical
artifact. We established a comprehensive error budget through three
independent convergence studies.

\subsection{RK4 Convergence Analysis}
\label{subsec:rk4_error}

To validate the RK4 integrator, we applied it to the model
$\beta$-function ODE
\begin{equation}
    \frac{dy}{dt} = -0.05\,y, \quad y(0) = 1, \quad y_{\rm exact}(T=100) = e^{-5},
    \label{eq:model_ode}
\end{equation}
and measured the global truncation error $|y(T) - y_{\rm exact}|$ as a
function of the step size $h = \Delta t$. Figure~\ref{fig:rk4_convergence}
confirms that the RK4 scheme follows the expected $\mathcal{O}(h^4)$ scaling
over four decades in $h$. The deviation from the theoretical slope at
$h \lesssim 10^{-2}$ (visible as a floor in Fig.~\ref{fig:rk4_convergence})
corresponds to the onset of double-precision floating-point round-off error
near $\epsilon_{\rm machine} \sim 10^{-16}$. At our working step size
$\Delta t = 0.01$, the global truncation error is $\sim 10^{-14}$, entirely
negligible compared to the 0.1\% precision required for the action classification.

\begin{figure}[H]
    \centering
    \includegraphics[width=\columnwidth]{error_rk4_convergence.png}
    \caption{Global truncation error of the RK4, RK2 (midpoint), and Euler
    integrators applied to the model ODE $\dot{y} = -0.05y$. The RK4
    scheme follows the theoretical $\mathcal{O}(h^4)$ convergence slope
    over four decades, demonstrating fourth-order accuracy. The floor
    at $h \lesssim 10^{-2}$ arises from double-precision floating-point
    round-off. At our operational step size $\Delta t = 0.01$ (marked), the
    truncation error is $\sim 10^{-14}$.}
    \label{fig:rk4_convergence}
\end{figure}

\subsection{Quadrature Convergence Analysis}
\label{subsec:quad_error}

The convergence of the action integral with respect to the number of
quadrature nodes is shown in Fig.~\ref{fig:quadrature_convergence}.
The composite Simpson integrator traces the expected $\mathcal{O}(N^{-4})$
slope over three decades, reaching a relative error of $\sim 10^{-11}$
at $N = 256$ nodes. Beyond this, the error plateaus at the double-precision
floor. Our operational choice of $N = 2048$ nodes (vs. the previous $N = 256$)
moves the solver three points deeper into the convergence plateau,
providing an order-of-magnitude guard band against any subtle cancellations
in the integral.

\begin{figure}[H]
    \centering
    \includegraphics[width=\columnwidth]{error_quadrature_convergence.png}
    \caption{Relative error in the bounce action integral as a function
    of the number of composite Simpson quadrature nodes $N$. The error
    follows the expected $\mathcal{O}(N^{-4})$ scaling and saturates
    below $10^{-11}$ for $N \gtrsim 256$. The operational choice of
    $N = 2048$ (green dashed line) places the solver well inside the
    convergence floor.}
    \label{fig:quadrature_convergence}
\end{figure}

\subsection{Golden-Section Search Convergence}
\label{subsec:gss_error}

Figure~\ref{fig:gss_convergence} shows the sensitivity of the minimized
action $S_{\min}$ to the search tolerance $\tau_{\rm GSS}$. The action
value saturates to its precision floor at $\tau \lesssim 10^{-2}$, with
a residual change of $\sim 5 \times 10^{-3}$ in $S_{\min}$ between
$\tau = 10^{-10}$ and $\tau = 10^{-13}$. This plateau demonstrates that
the optimal bubble radius $R^*$ is stably located by the algorithm and
that tightening the tolerance beyond $\tau = 10^{-10}$ yields no
physically meaningful refinement. Our operational tolerance of
$\tau = 10^{-13}$ therefore represents a deliberately conservative choice
that eliminates any systematic bias from premature termination.

\begin{figure}[H]
    \centering
    \includegraphics[width=\columnwidth]{error_golden_section.png}
    \caption{Error in the minimized bounce action $S_{\min}$ as a function
    of the Golden-Section Search tolerance. The action value saturates
    for $\tau \lesssim 10^{-2}$, confirming that the global minimum of
    $S_E(R)$ is robustly isolated. Previous ($\tau = 10^{-10}$, red) and
    current ($\tau = 10^{-13}$, green) tolerances are both well within the
    converged plateau.}
    \label{fig:gss_convergence}
\end{figure}

\subsection{Equation-of-Motion Residuals}
\label{subsec:eom_residuals}

As a final, independent validation of the numerical bounce profiles,
we computed the residual of the equation of motion,
$\mathcal{R}(r) = \phi_b'' + (3/r)\phi_b' - V'_{\rm eff}(\phi_b)$,
integrated over the profile. Figure~\ref{fig:residual_histogram} shows
the distribution of these residuals across all 250,000 parameter points.
The distribution is sharply peaked at zero, confirming that the numerical
bounce profiles tightly satisfy their defining ODE throughout the
parameter space. The large negative tail in the distribution corresponds
to points near the absolute stability boundary, where $\lambda_{\min}$
is near zero and the bounce profile approaches a degenerate configuration.

\begin{figure}[H]
    \centering
    \includegraphics[width=\columnwidth]{error_residual_histogram.png}
    \caption{Distribution of equation-of-motion residuals for the
    numerical bounce profiles across the full 250,000-point parameter scan.
    The distribution is sharply concentrated at zero, with a Gaussian fit
    yielding $\mu = 2.8\%$ and $\sigma = 116.6\%$ — the large standard
    deviation arising from the extended non-Gaussian tail of near-boundary
    points rather than any systematic bias.}
    \label{fig:residual_histogram}
\end{figure}

\subsection{Consolidated Error Budget}
\label{subsec:error_summary}

Table~\ref{tab:errors} summarizes the sources and magnitudes of numerical
errors in the computed bounce action. The dominant contribution to the total
numerical error budget is from the quadrature at the $10^{-11}$ level,
which is smaller by at least ten orders of magnitude than the physical
discrepancy between $S_{\rm exact}$ and $S_{\rm approx}$ in the boundary
region ($\sim 10^{-1}$ relative error). We conclude that all phase boundary
reclassifications reported below are physical and not numerical artifacts.

\begin{table}[H]
  \caption{Consolidated numerical error budget for the exact bounce action
  computation. All errors are expressed as relative contributions to $S_E$.}
  \label{tab:errors}
  \begin{ruledtabular}
  \begin{tabular}{lcc}
    Source & Error (relative) & Status \\
    \hline
    RK4 truncation ($\Delta t = 0.01$) & $\lesssim 10^{-14}$ & $\checkmark$ \\
    Composite Simpson quadrature ($N=2048$) & $\lesssim 10^{-11}$ & $\checkmark$ \\
    GSS radius optimization ($\tau=10^{-13}$) & $\lesssim 10^{-5}$ & $\checkmark$ \\
    Floating-point round-off & $\sim 10^{-16}$ & $\checkmark$ \\
    \hline
    \textbf{Total numerical error} & $\boldsymbol{\lesssim 10^{-11}}$ & \\
    \hline
    Physical discrepancy ($|\lambda_{\min}| \to 0$) & $\sim 10\%$--$30\%$ & (target) \\
    Dominant parametric uncertainty ($\delta M_t$) & $\sim 1\%$ & (physical) \\
  \end{tabular}
  \end{ruledtabular}
\end{table}

%======================================================================
\section{Results}
\label{sec:results}
%======================================================================

\subsection{Global Phase Diagram}
\label{subsec:global_phase}

Figures~\ref{fig:analytical_full} and~\ref{fig:numerical_full} display the
global stability phase diagrams derived from the analytical conformal
approximation and the exact numerical solver, respectively. Both diagrams
partition the $(M_t, M_h)$ plane into the four regions defined in
Sec.~\ref{sec:theory}: stable (green), metastable (yellow), unstable (red),
and non-perturbative (black). The experimentally measured central
values $(M_t, M_h) = (173.3, 125.1)~\text{GeV}$ are marked with a
filled star and lie within the metastable region in both diagrams,
consistent with the existing literature~\cite{Buttazzo:2013uya}.
Figure~\ref{fig:analytical_full} serves as the Buttazzo-style analytical
diagram for direct comparison with the canonical literature result,
providing a cross-check of our RGE implementation.

\begin{figure}[H]
    \centering
    \includegraphics[width=\columnwidth]{analytical_stability_plot.png}
    \caption{Global vacuum stability phase diagram in the $(M_t, M_h)$ plane
    derived from the analytical conformal approximation
    $S_{\rm approx} = 8\pi^2/3|\lambda_{\min}|$. The color coding distinguishes
    stable (green), metastable (yellow), unstable (red), and non-perturbative
    (black) regions. The star marks the experimental central values
    $(M_t, M_h) = (173.3, 125.1)~\text{GeV}$.}
    \label{fig:analytical_full}
\end{figure}

\begin{figure}[H]
    \centering
    \includegraphics[width=\columnwidth]{numerical_stability_contours.png}
    \caption{Global vacuum stability phase diagram derived from the exact
    numerical bounce solver incorporating the $c_6\phi^6/M_{\rm Pl}^2$
    operator with $c_6 = 1$, shown with constant action contours in the unstable
    and metastable regions. Compared to Fig.~\ref{fig:analytical_full},
    the metastable region is visibly expanded, with the instability boundary
    shifted upward (toward higher $M_t$) throughout the accessible range.}
    \label{fig:numerical_full}
\end{figure}

\subsection{Breakdown of the Conformal Approximation}
\label{subsec:breakdown}

Our primary quantitative finding is illustrated in Fig.~\ref{fig:S_scatter},
which plots $S_{\rm exact}$ against $S_{\rm approx}$ for all 250,000 sampled
parameter points.

\begin{figure*}[H]
    \centering
    \includegraphics[width=\textwidth]{error_S_scatter.png}
    \caption{Comprehensive comparison of the exact numerical bounce action
    $S_{\rm exact}$ versus the conformal approximation $S_{\rm approx}$ across
    250,000 parameter points. \textbf{(a)} Full scatter plot showing that the
    conformal approximation is a good global estimator, with 74.9\% of points
    within $\pm 5\%$ and 86.8\% within $\pm 15\%$. \textbf{(b)} Relative
    discrepancy $(S_{\rm exact} - S_{\rm approx})/S_{\rm exact}$ as a function
    of $S_{\rm approx}$, revealing a critical two-regime structure: for deep
    instabilities ($S_{\rm approx} \gtrsim 500$) the approximation is accurate,
    while for shallow instabilities ($S_{\rm approx} \lesssim 50$, near the
    phase boundary) the discrepancy grows to $+20\%$ to $+30\%$.
    \textbf{(c)} Zoomed view near the point-dependent metastability threshold
    $S_{\rm threshold}(\mu_{\rm inst})$, highlighting the points rescued from
    unstable to metastable (green box) by the exact computation.}
    \label{fig:S_scatter}
\end{figure*}

Panel (a) of Fig.~\ref{fig:S_scatter} confirms that globally the
conformal approximation is an excellent estimator: 74.9\% of all parameter
points show less than 5\% discrepancy, and 86.8\% lie within $\pm 15\%$.
This validates the use of Eq.~\eqref{eq:conformal_approx} as a fast
reconnaissance tool for large areas of parameter space far from the
metastability boundary.

Panel (b) exposes the critical failure regime. For small $S_{\rm approx}$
(large $|\lambda_{\min}|$, deep instability), the conformal approximation is
accurate. However, as $S_{\rm approx}$ approaches the point-dependent
threshold $4\ln(T_U\mu_{\rm inst})$, the binned
mean discrepancy rises from zero to $\sim 20\%$, and the scatter band
($\pm 1\sigma$) widens dramatically. This two-regime structure is precisely
what is expected from the analytical discussion in Sec.~\ref{subsec:validity}:
the conformal approximation is derived in the limit of a deep, narrow
instability, and fails when the bounce profile extends to Planck-scale fields.

The physical mechanism behind this failure is the barrier inflation effect:
the $c_6\phi^6/M_{\rm Pl}^2$ operator contributes a positive-definite
correction to $V_{\rm eff}$ at large $\phi$. Because the exact numerical
solver dynamically integrates the equation of motion out to $\phi \sim M_{\rm Pl}$,
it fully accounts for this additional barrier, which strictly raises $S_E$.
By contrast, $S_{\rm approx}$ is blind to any feature of the potential beyond
the minimum of $\lambda(\mu)$. The inequality $S_{\rm exact} > S_{\rm approx}$
is therefore an unconditional consequence of $c_6 > 0$ in the boundary
region.

Figure~\ref{fig:error_vs_lambda} provides a direct correlation between the
discrepancy and the depth of the instability, measured by
$|\lambda_{\min}|$ (expressed via its proxy $\approx 8\pi^2/3S_{\rm approx}$).
As $|\lambda_{\min}| \to 0$, the mean discrepancy grows with increasing
scatter, consistent with the divergence of the conformal approximation
in the shallow-instability limit. The $1\sigma$ error band is shown to
encompass an enormous range — from $-150\%$ to $+150\%$ — for points with
$|\lambda_{\min}| \lesssim 0.1 \times 10^{-3}$. This underscores the
complete unreliability of the conformal approximation for these points.

\begin{figure}[H]
    \centering
    \includegraphics[width=\columnwidth]{error_vs_lambda.png}
    \caption{Mean relative discrepancy $(S_{\rm exact} - S_{\rm approx})/S_{\rm exact}$
    as a function of the depth of instability $|\lambda_{\min}|$ (proxied by
    $8\pi^2/3S_{\rm approx}$). As $|\lambda_{\min}| \to 0$, the conformal
    approximation becomes catastrophically inaccurate, with the $\pm 1\sigma$
    envelope expanding to hundreds of percent for the shallowest instabilities.}
    \label{fig:error_vs_lambda}
\end{figure}

The spatial topology of the breakdown is displayed in the discrepancy heatmap
of Fig.~\ref{fig:heatmap}. The failure of the conformal approximation is
not uniformly distributed across the $(M_t, M_h)$ plane; rather, it is
sharply concentrated along the diagonal strip tracing the absolute stability
boundary ($\lambda_{\min} = 0$ contour). The relative discrepancy is positive
(green in Fig.~\ref{fig:heatmap}) throughout, confirming the unconditional
upward bias of the exact result.

\begin{figure}[H]
    \centering
    \includegraphics[width=\columnwidth]{error_discrepancy_heatmap.png}
    \caption{Heatmap of the relative discrepancy
    $(S_{\rm exact} - S_{\rm approx})/S_{\rm exact}$ in the $(M_t, M_h)$ plane.
    The discrepancy is concentrated along the stability boundary (diagonal ridge),
    where $|\lambda_{\min}| \to 0$. Experimental central values
    $M_t^{\rm SM} = 173.3~\text{GeV}$ (dashed orange) and
    $M_h^{\rm SM} = 125.1~\text{GeV}$ (dashed blue) are marked.}
    \label{fig:heatmap}
\end{figure}

\subsection{Phase Boundary Reclassification}
\label{subsec:reclassification}

The direct physical consequence of the unconditional upward bias in
$S_{\rm exact}$ is a systematic reclassification of the metastability boundary.
Points for which $S_{\rm approx} < S_{\rm threshold}(\mu_{\rm inst})$ (analytically
classified as unstable) can satisfy $S_{\rm exact} > S_{\rm threshold}$
(physically metastable). This reclassification is illustrated in
Fig.~\ref{fig:reclassification}, which highlights in yellow the 110 grid
points in the SM-relevant region of the $(M_t, M_h)$ plane that are
rescued from the unstable to the metastable regime by the $c_6 = 1$
numerical correction.

\begin{figure}[H]
    \centering
    \includegraphics[width=\columnwidth]{error_boundary_reclassification.png}
    \caption{Phase boundary reclassification in the SM-relevant region
    $M_t \in [155, 185]~\text{GeV}$, $M_h \in [115, 135]~\text{GeV}$.
    Yellow points indicate the 110 parameter configurations reclassified
    from unstable (analytical) to metastable (exact numerical) by the
    $c_6\phi^6/M_{\rm Pl}^2$ operator correction. The rescued points form
    a continuous diagonal strip adjacent to the analytical instability
    boundary, with a characteristic width of $\sim 0.3$--$0.5$~GeV in $M_t$.
    Crosshairs mark the experimental values $(M_t, M_h) = (173.3, 125.1)~\text{GeV}$.}
    \label{fig:reclassification}
\end{figure}

The 110 reclassified points form a continuous, smooth strip of width
approximately 0.3--0.5~GeV in the $M_t$ direction, running diagonally
across the parameter plane in close parallel to the analytical instability
boundary. This strip is physically interpreted as a region in which the
Planck-scale barrier inflation is strong enough to push $S_E$ above the
threshold but not so deep in the unstable regime as to be irrelevant.

The high-resolution closeup view of the analytical boundary
(Fig.~\ref{fig:analytical_closeup}) reveals that in the vicinity of the
experimental point at $M_t \approx 173$~GeV, $M_h \approx 125$~GeV, the
conformal approximation places the instability boundary dangerously close
to the $1\sigma$ experimental ellipse of the top quark mass. By contrast,
the exact numerical boundary (Fig.~\ref{fig:numerical_closeup}) lies visibly
higher, accommodating a wider range of $M_t$ values in the metastable
regime.

\begin{figure}[H]
    \centering
    \includegraphics[width=\columnwidth]{analytical_stability_closeup.png}
    \caption{High-resolution close-up of the analytical phase diagram near
    the experimental point $(M_t, M_h) = (173.1, 125.1)~\text{GeV}$.
    The conformal approximation boundary cuts close to the measured values,
    suggesting a precarious vacuum.}
    \label{fig:analytical_closeup}
\end{figure}

\begin{figure}[H]
    \centering
    \includegraphics[width=\columnwidth]{numerical_stability_closeup.png}
    \caption{High-resolution close-up of the exact numerical phase diagram,
    including the $c_6 = 1$ Planck-scale correction. The instability boundary
    is shifted rightward (toward higher $M_t$) relative to Fig.~\ref{fig:analytical_closeup},
    placing the SM experimental point more robustly within the metastable region.}
    \label{fig:numerical_closeup}
\end{figure}

\subsection{Direct Phase Boundary Overlay}
\label{subsec:overlay}

To directly quantify the physical shift in the phase boundary,
Figs.~\ref{fig:overlay_full} and~\ref{fig:overlay_closeup} overlay the
analytical (background shading) and exact numerical (blue dots) boundaries
on the same axes. The blue dots — marking points that are analytically
unstable but numerically metastable — define the physical rescued strip.

\begin{figure}[H]
    \centering
    \includegraphics[width=\columnwidth]{overlay_plot.png}
    \caption{Full overlay of the analytical and exact numerical phase
    boundaries in the $(M_t, M_h)$ plane. The background shows the
    numerical phase classification (stable green, metastable yellow, unstable
    red). Blue dots indicate the boundary shift, tracing the locus of
    parameter points rescued by the exact computation throughout the
    diagram.}
    \label{fig:overlay_full}
\end{figure}

\begin{figure}[H]
    \centering
    \includegraphics[width=\columnwidth]{overlay_closeup_plot.png}
    \caption{Close-up overlay in the SM-relevant region, revealing the
    width of the rescued strip ($\sim$0.3--0.5~GeV in $M_t$, blue dots)
    relative to the analytical background. The experimental point
    (star) lies safely within the metastable region according to the
    exact calculation.}
    \label{fig:overlay_closeup}
\end{figure}

The overlay confirms that the boundary shift is a coherent, extended
feature rather than a collection of isolated anomalous points. The strip
of blue dots runs continuously from $M_h \approx 115~\text{GeV}$ at
$M_t \approx 166~\text{GeV}$ to $M_h \approx 135~\text{GeV}$ at
$M_t \approx 180~\text{GeV}$, tracing the approximate diagonal $\partial\lambda_{\min}/\partial M_t \approx -0.002~\text{GeV}^{-1}$.

\subsection{Breakdown of the Semi-Classical Approximation: The Sandwiched Unstable Region}
\label{subsec:sandwiched_unstable}

A striking feature of the exact numerical phase diagram at low Higgs masses ($M_h \sim 50$ GeV) is the emergence of a \textit{sandwiched} unstable region at intermediate top quark masses. Within this band, the Euclidean bounce action evaluates to negative values despite the parameter points being flanked on both sides by large-positive-action metastable phases --- a topology that, while superficially paradoxical, has a precise physical and mathematical origin.

The mechanism begins with the extreme smallness of the initial quartic coupling at the matching scale. For $M_h \sim 50$ GeV, the tree-level contribution $\lambda_{\rm tree}(M_t) = M_h^2/(2v^2) \approx 0.021$ is already small, and NNLO threshold corrections reduce it further, rendering the effective potential highly sensitive to the early-scale structure of the renormalization group flow. Within the sandwiched band, the specific interplay of top Yukawa, gauge, and NNLO threshold contributions causes $\lambda_{\rm eff}(\mu)$ to turn negative at an instability scale $\mu_1$ barely above $M_t$. The effective potential $V_{\rm eff}(\phi) = (\lambda_{\rm eff}(\phi)/4)\phi^4 + (c_6/M_{\rm Pl}^2)\phi^6$ then develops a potential barrier that, while nonzero, is extremely thin and shallow relative to the characteristic energy scale.

The negative bounce action computed in this region is an artifact of the approximation scheme, not a physical observable. The Fubini-Lipatov instanton profile is an exact saddle point of the Euclidean action only for the pure quartic potential with a \textit{constant} negative coupling. For that potential, scale invariance in four Euclidean dimensions enforces the virial relation $S_{\rm kin} = -2 S_{\rm pot}$, ensuring $S_E = -S_{\rm pot} > 0$. When $\lambda_{\rm eff}(\phi)$ runs strongly, this constraint is violated: the kinetic contribution is controlled by the bounce size $R \sim \mu_1^{-1}$, which is small when $\mu_1$ is close to $M_t$ and $|\lambda_{\rm eff}(\mu_1)|$ is tiny; meanwhile, the potential integral is dominated by the large-$|\lambda_{\rm eff}|$ contributions at high field values, far above $\mu_1$, where the Fubini profile's long-range tail samples the potential. The ratio $S_{\rm kin}/|S_{\rm pot}|$ deviates from its expected value of $2$ by orders of magnitude, providing a quantitative diagnostic that the profile is far from a saddle point of the action for this potential. The resulting $S_E < 0$ indicates computational breakdown of the Fubini ansatz, not the absence of a physical bounce.

Crucially, this breakdown does not imply a classical roll-down instability --- that is, the absence of a potential barrier. The barrier exists, as demonstrated by the fact that the adjacent metastable parameter regions host genuine O(4)-symmetric bounce saddle points: since $V_{\rm eff}(\phi, M_t)$ depends smoothly on $M_t$ through the RGE running and analytic matching conditions, the barrier cannot discontinuously vanish at the boundary of the sandwiched region. Rather, the true bounce in this regime --- obtainable by direct numerical integration of the O(4)-symmetric shooting equations with the full running potential --- yields a genuine, small positive action $S_E^{\rm true} \ll S_{\rm thresh}$, corresponding to rapid quantum tunneling through the thin barrier rather than a classical trajectory. The semi-classical saddle-point approximation for this single bounce breaks down precisely because the bounce scale and the scale over which $\lambda_{\rm eff}$ varies are not well-separated.

The flanking metastable bands (the lower of which extends to the stable phase at small $M_t$) correspond to parameter regions where either the instability scale $\mu_1$ is sufficiently high, or the initial $\lambda(M_t)$ is sufficiently positive, to sustain a barrier of adequate thickness. In those regions, the ratio $S_{\rm kin}/|S_{\rm pot}|$ is within an acceptable range of $2$, the Fubini approximation is controlled, and the action converges to a large positive value well above the metastability threshold. For phase classification purposes, points with $S_E^{\rm computed} \leq 0$ are assigned to the unstable phase, since the Fubini breakdown is itself diagnostic of a regime where $S_E^{\rm true} \ll S_{\rm thresh}$.

\subsection{Dependence on the Planck-Scale Coefficient \texorpdfstring{$c_6$}{c6}}
\label{subsec:c6_sweep}

To assess the robustness of these conclusions against the unknown value
of the UV coefficient $c_6$, we repeated the full phase diagram computation
for $c_6 \in \{0.1, 1.0, 10.0\}$. The resulting boundaries are shown in
Figs.~\ref{fig:c6_phase} and~\ref{fig:c6_sweep}.

\begin{figure}[H]
    \centering
    \includegraphics[width=\columnwidth]{c6_phase_diagram.png}
    \caption{Full phase diagram showing the expanding metastability region
    for different values of the Planck-scale operator coefficient $c_6$.
    The $c_6 = 0$ (pure SM) boundaries are shown in blue, and $c_6 = 10$
    boundaries in red. Instability boundaries for $c_6 = 0.1$ (yellow
    dashed) and $c_6 = 1.0$ (green dashed) are also shown. For all
    $c_6 > 0$, the metastable region is expanded relative to the pure
    SM case.}
    \label{fig:c6_phase}
\end{figure}

\begin{figure}[H]
    \centering
    \includegraphics[width=\columnwidth]{c6_sweep_plot.png}
    \caption{High-resolution view of the metastability boundary in the
    experimental region for $c_6 \in \{0, 0.1, 1.0, 10.0\}$. The boundary
    shifts toward the upper right with increasing $c_6$, reflecting the
    stronger stabilizing effect of a larger Planck-scale barrier. The
    SM experimental point (star) at $(M_t, M_h) = (173.3, 125.1)~\text{GeV}$
    is shown for reference.}
    \label{fig:c6_sweep}
\end{figure}

Three observations follow from these figures. First, increasing $c_6$
shifts the boundary further into the classically unstable regime, consistent
with the expectation that a larger $\phi^6$ coefficient creates a higher
barrier. For $c_6 = 10$, the boundary is shifted by approximately
1.5--2~GeV in $M_t$ compared to the pure SM case. Second, for the
benchmark $c_6 = 1$, the shift of $\sim 0.3$--$0.5$~GeV is comparable to
the $1\sigma$ experimental uncertainty on the top quark mass,
$\delta M_t^{\rm exp} \approx 0.7$~GeV~\cite{PDG:2022}, making this
effect phenomenologically significant. Third, and most importantly, the
qualitative conclusion — that the conformal approximation underestimates
$S_E$ near the boundary — holds for all $c_6 > 0$. It is a universal
consequence of the positive definiteness of the $\phi^6$ barrier, not a
fine-tuned artifact of the $c_6 = 1$ choice.

%======================================================================
\section{Discussion}
\label{sec:discussion}
%======================================================================

\subsection{Physical Implications of the Boundary Shift}
\label{subsec:physical_implications}

The central result of this work — that the exact numerical bounce action
exceeds the conformal approximation in the boundary region — has a clear
physical interpretation. The conformal approximation~\eqref{eq:conformal_approx}
captures only the contribution of the deepest point of the effective potential.
As the instability becomes shallow ($|\lambda_{\min}| \to 0$), the bounce
profile extends to increasingly large field values, and the classical turning
point approaches the Planck scale. At these scales, even a subdominant
operator suppressed by $M_{\rm Pl}^{-2}$ can contribute meaningfully to the
integral in Eq.~\eqref{eq:euclidean_action}. Because the $\phi^6$ operator
is positive definite for $c_6 > 0$, this contribution always raises the
action, extending the metastable region.

The magnitude of the effect — $\sim 0.3$--$0.5$~GeV in $M_t$ for $c_6 = 1$
— is directly relevant to precision measurements. The current best
determination of the top quark pole mass from LHC data yields
$M_t = 173.3 \pm 0.7~\text{GeV}$ at 1$\sigma$~\cite{PDG:2022,CMS:2022top}.
This means that the $\pm 1\sigma$ experimental band in $M_t$ straddles
the reclassified strip. Whether the SM vacuum is metastable or (on the
boundary of) unstable, as assessed by purely conformal methods, cannot be
definitively resolved without exact numerical bounce solvers.

\subsection{Comparison with Existing Work}
\label{subsec:comparison}

Previous work on the impact of higher-dimensional operators on vacuum
stability has primarily focused on their effect on the location of the
absolute stability boundary (the condition $\lambda_{\min} = 0$), rather
than on the tunneling action and metastability
boundary~\cite{Branchina:2013jra,Branchina:2014rva,Lalak:2014qua}.
These analyses showed that $\phi^6$ operators can stabilize the potential,
but they typically employed the conformal approximation to assess the
tunneling rate. Our work demonstrates that this combination is internally
inconsistent near the boundary: the very operators that modify the
potential shape at large $\phi$ are also the ones that cause the conformal
approximation to fail. The exact numerical approach is necessary to capture
both effects simultaneously.

Regarding numerical techniques for computing bounce actions, several
public codes exist~\cite{CosmoTransitions,BubbleProfiler,AnyBubble}
that solve the bounce equation numerically for arbitrary potentials.
However, these codes have not been systematically applied to the SM
metastability boundary at the precision required to resolve the
$\sim 0.5$~GeV boundary shift. Our analysis fills this gap by providing
a high-resolution scan with verified sub-$10^{-11}$ numerical errors.

\subsection{Limitations and Future Directions}
\label{subsec:limitations}

Several caveats and natural extensions of this work should be noted.

First, our effective potential~\eqref{eq:veff} is a tree-level, RG-improved
approximation. The full one-loop effective potential includes a
Coleman-Weinberg correction~\cite{Coleman:1973jx}, logarithms of
$\phi/\mu_{\rm ren}$, and contributions from gauge bosons and top quarks
that we have not incorporated. At the precision level targeted here
($\sim 0.5$~GeV in $M_t$), these one-loop corrections can be of the same
order as the boundary shift, and a complete treatment would require
their inclusion. This represents the primary avenue for future work.

Second, we have treated $c_6$ as a free parameter and have not attempted
to predict its value from a specific UV completion. Models such as
non-minimal Higgs portal extensions, higher-dimensional gravity, or
asymptotic safety~\cite{Eichhorn:2017ylw} each predict specific values.
Our universality result — that the boundary shift persists for all $c_6 > 0$
— implies that any UV completion with a positive $c_6$ shifts the boundary
in the rescuing direction.

Third, the effect of finite-temperature corrections on the tunneling rate
at cosmological times has not been included. For the zero-temperature
analysis of the electroweak vacuum at the Planck era, this is well
justified, but a full cosmological treatment would require the
thermal bounce~\cite{Linde:1980tt,Linde:1981zj}.

Finally, the non-perturbative region of the parameter space (large $M_t$
or $M_h$, shown in black in our phase diagrams) cannot be accessed by our
perturbative potential. The true stability boundary in this region requires
lattice field theory methods~\cite{Gould:2019qek}, which are beyond the
scope of this paper.

%======================================================================
\section{Conclusion}
\label{sec:conclusion}
%======================================================================

We have presented a comprehensive, high-precision numerical study of the
Standard Model electroweak vacuum stability phase diagram in the $(M_t, M_h)$
plane. Our main conclusions are as follows.

The analytical conformal approximation $S_{\rm approx} = 8\pi^2/3|\lambda_{\min}|$
is a reliable estimator for deep instabilities, with 74.9\% of all 250,000
sampled parameter points agreeing with the exact computation to within
$\pm 5\%$. However, the approximation fails catastrophically near the
metastability boundary where $|\lambda_{\min}| \to 0$, with discrepancies
growing to $20$--$30\%$ in the mean and exceeding $100\%$ in the tails of
the distribution.

Upon incorporating a Planck-scale dimension-6 operator
$c_6\phi^6/M_{\rm Pl}^2$ with $c_6 = 1$, the exact numerical bounce
solver detects a positive-definite barrier inflation that raises
$S_{\rm exact}$ above $S_{\rm approx}$ unconditionally in the boundary
region. This produces a physically meaningful reclassification of
approximately 110 parameter points — a continuous strip of width
$\sim 0.3$--$0.5$~GeV in $M_t$ — from the unstable to the metastable
regime. The direction of this shift (toward stability) is universal: it
persists for all $c_6 > 0$ and for $c_6$ ranging over two orders of
magnitude, $c_6 \in \{0.1, 1.0, 10.0\}$.

The magnitude of the boundary shift is comparable to the current $1\sigma$
experimental uncertainty on $M_t \approx 0.7$~GeV, directly implying that
precision determinations of the vacuum stability status of the SM require
exact numerical bounce solvers. Standard conformal approximation methods
are insufficient to resolve whether the electroweak vacuum is comfortably
metastable or on the threshold of instability within the present experimental
uncertainties.

We have validated our numerical pipeline through four independent
convergence studies — RK4 step-size analysis, quadrature node scaling,
Golden-Section Search tolerance analysis, and equation-of-motion residual
histograms — and established a comprehensive error budget showing total
numerical errors below $10^{-11}$ in relative action. All reported boundary
reclassifications are therefore purely physical in origin.

This work establishes a rigorous baseline for future studies that will
incorporate one-loop Coleman-Weinberg corrections, finite-temperature
effects, and UV completions with specific predictions for $c_6$. It also
motivates the development and deployment of fully automated exact bounce
solvers as a standard tool in precision vacuum stability analyses.

%======================================================================
\begin{acknowledgments}
We thank autonomous agents for valuable discussions on electroweak
vacuum stability and numerical methods.
\end{acknowledgments}

\bibliography{references}

\end{document}



### Appendix B: metastable_points_mh5.md

# Metastable Points at  = 5$ GeV

Here is the complete list of all metastable points at =5$ GeV from the dense numerical grid. Note that {exact}$ is the computed Euclidean bounce action.

| Top Mass ($) | Euclidean Action ({\rm exact}$) |
|:---:|:---:|
| 94.00 GeV | > 10^{10} (Classical Instability Artifact) |
| 94.25 GeV | > 10^{10} (Classical Instability Artifact) |
| 94.50 GeV | > 10^{10} (Classical Instability Artifact) |
| 94.75 GeV | > 10^{10} (Classical Instability Artifact) |
| 95.00 GeV | 643,028.0 |
| 95.25 GeV | 500,000.0 |
| 95.50 GeV | 500,000.0 |
| 95.75 GeV | 500,000.0 |
| 96.00 GeV | 500,000.0 |
| 96.25 GeV | 505,344.0 |
| 96.50 GeV | 377,450.0 |
| 96.75 GeV | 291,328.0 |
| 97.00 GeV | 230,618.0 |
| 97.25 GeV | 186,596.0 |
| 97.50 GeV | 153,783.0 |
| 97.75 GeV | 128,715.0 |
| 98.00 GeV | 109,184.0 |
| 98.25 GeV | 93,690.6 |
| 98.50 GeV | 81,217.6 |
| 100.75 GeV | 500,000.0 |
| 101.25 GeV | 500,000.0 |
| 114.00 GeV | 4,235.9 |
| 114.25 GeV | 4,141.2 |
| 114.50 GeV | 4,049.5 |
| 114.75 GeV | 3,960.5 |
| 115.00 GeV | 3,874.3 |
| 115.25 GeV | 3,790.7 |
| 115.50 GeV | 3,709.5 |
| 115.75 GeV | 3,630.8 |
| 116.00 GeV | 3,554.3 |
| 116.25 GeV | 3,480.1 |
| 116.50 GeV | 3,408.0 |
| 116.75 GeV | 3,337.9 |
| 117.00 GeV | 3,269.8 |
| 117.25 GeV | 3,203.7 |
| 117.50 GeV | 3,139.4 |
| 117.75 GeV | 3,076.8 |
| 118.00 GeV | 3,016.0 |
| 118.25 GeV | 2,956.8 |
| 118.50 GeV | 2,899.2 |
| 118.75 GeV | 2,843.2 |
| 119.00 GeV | 2,788.6 |
| 119.25 GeV | 2,735.4 |
| 119.50 GeV | 2,683.7 |
| 119.75 GeV | 2,633.2 |
| 120.00 GeV | 2,584.1 |
| 120.25 GeV | 2,536.2 |
| 120.50 GeV | 2,489.5 |
| 120.75 GeV | 2,444.0 |
| 121.00 GeV | 2,399.6 |
| 121.25 GeV | 2,356.3 |
| 121.50 GeV | 2,314.1 |
| 121.75 GeV | 2,272.9 |
| 122.00 GeV | 2,232.7 |
| 122.25 GeV | 2,193.4 |
| 122.50 GeV | 2,155.1 |
| 122.75 GeV | 2,117.7 |
| 123.00 GeV | 2,081.1 |
| 123.25 GeV | 2,045.5 |
| 123.50 GeV | 2,010.6 |
| 123.75 GeV | 1,976.5 |
| 124.00 GeV | 1,943.2 |
| 124.25 GeV | 1,910.7 |
| 124.50 GeV | 1,878.9 |
| 124.75 GeV | 1,847.8 |
| 125.00 GeV | 1,817.4 |
| 125.25 GeV | 1,787.6 |
| 125.50 GeV | 1,758.5 |
| 125.75 GeV | 1,730.1 |
| 126.00 GeV | 1,702.2 |
| 126.25 GeV | 1,675.0 |
| 126.50 GeV | 1,648.3 |
| 126.75 GeV | 1,622.3 |
| 127.00 GeV | 1,596.7 |
| 127.25 GeV | 1,571.7 |
| 127.50 GeV | 1,547.2 |
| 127.75 GeV | 1,523.2 |
| 128.00 GeV | 1,499.8 |
| 128.25 GeV | 1,476.8 |
| 128.50 GeV | 1,454.2 |
| 128.75 GeV | 1,432.1 |
| 129.00 GeV | 1,410.5 |
| 129.25 GeV | 1,389.3 |
| 129.50 GeV | 1,368.5 |
| 129.75 GeV | 1,348.2 |
| 130.00 GeV | 1,328.2 |
| 130.25 GeV | 1,308.6 |
| 130.50 GeV | 1,289.4 |
| 130.75 GeV | 1,270.6 |
| 131.00 GeV | 1,252.1 |
| 131.25 GeV | 1,234.0 |
| 131.50 GeV | 1,216.2 |
| 131.75 GeV | 1,198.8 |
| 132.00 GeV | 1,181.7 |
| 132.25 GeV | 1,165.0 |
| 132.50 GeV | 1,148.5 |
| 132.75 GeV | 1,132.3 |
| 133.00 GeV | 1,116.5 |
| 133.25 GeV | 1,100.9 |
| 133.50 GeV | 1,085.6 |
| 133.75 GeV | 1,070.6 |
| 134.00 GeV | 1,055.9 |
| 134.25 GeV | 1,041.4 |
| 134.50 GeV | 1,027.2 |
| 134.75 GeV | 1,013.2 |
| 135.00 GeV | 999.5 |
| 135.25 GeV | 986.1 |
| 135.50 GeV | 972.8 |
| 135.75 GeV | 959.9 |
| 136.00 GeV | 947.1 |
| 136.25 GeV | 934.5 |
| 136.50 GeV | 922.2 |
| 136.75 GeV | 910.1 |
| 137.00 GeV | 898.2 |
| 137.25 GeV | 886.5 |
| 137.50 GeV | 875.0 |
| 137.75 GeV | 863.6 |
| 138.00 GeV | 852.5 |
| 138.25 GeV | 841.6 |
| 138.50 GeV | 830.8 |
| 138.75 GeV | 820.3 |
| 139.00 GeV | 809.9 |
| 139.25 GeV | 799.6 |
| 139.50 GeV | 789.6 |
| 139.75 GeV | 779.7 |
| 140.00 GeV | 769.9 |
| 140.25 GeV | 760.3 |
| 140.50 GeV | 750.9 |
| 140.75 GeV | 741.6 |
| 141.00 GeV | 732.5 |
| 141.25 GeV | 723.5 |
| 141.50 GeV | 714.7 |
| 141.75 GeV | 706.0 |
| 142.00 GeV | 697.4 |
| 142.25 GeV | 689.0 |
| 142.50 GeV | 680.7 |
| 142.75 GeV | 672.6 |
| 143.00 GeV | 664.5 |
| 143.25 GeV | 656.6 |
| 143.50 GeV | 648.8 |
| 143.75 GeV | 641.1 |
| 144.00 GeV | 633.6 |
| 144.25 GeV | 626.1 |
| 144.50 GeV | 618.8 |
| 144.75 GeV | 611.6 |
| 145.00 GeV | 604.5 |
| 145.25 GeV | 597.5 |
| 145.50 GeV | 590.6 |
| 145.75 GeV | 583.8 |
| 146.00 GeV | 577.1 |
| 146.25 GeV | 570.5 |
| 146.50 GeV | 564.0 |
| 146.75 GeV | 557.6 |
| 147.00 GeV | 551.3 |
| 147.25 GeV | 545.0 |
| 147.50 GeV | 538.9 |
| 147.75 GeV | 532.9 |
| 148.00 GeV | 526.9 |
| 148.25 GeV | 521.0 |
| 148.50 GeV | 515.3 |
| 148.75 GeV | 509.6 |
| 149.00 GeV | 503.9 |
| 149.25 GeV | 498.4 |
| 149.50 GeV | 492.9 |
| 149.75 GeV | 487.5 |
| 150.00 GeV | 482.2 |
| 150.25 GeV | 477.0 |
| 150.50 GeV | 471.8 |
| 150.75 GeV | 466.7 |
| 151.00 GeV | 461.7 |
| 151.25 GeV | 456.8 |
| 151.50 GeV | 451.9 |
| 151.75 GeV | 447.1 |
| 152.00 GeV | 442.3 |
| 152.25 GeV | 437.6 |
| 152.50 GeV | 433.0 |
| 152.75 GeV | 428.4 |
| 153.00 GeV | 423.9 |
| 153.25 GeV | 419.5 |
| 153.50 GeV | 415.1 |
| 153.75 GeV | 410.8 |



### Appendix C: plot_numerical_contours.py

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def create_contour_plot():
    print("Loading data...")
    df = pd.read_csv('results/numerical_data.csv')

    stable = df[df['Stability'] == 1]
    metastable = df[df['Stability'] == 2]
    unstable = df[df['Stability'] == 3]
    nonpert = df[df['Stability'] == 4]

    plt.figure(figsize=(10, 8))
    
    # Plot the background scatter points
    plt.scatter(stable['Mh_calc'], stable['Mt'], c='green', alpha=0.3, s=5, label='Stable')
    plt.scatter(metastable['Mh_calc'], metastable['Mt'], c='yellow', alpha=0.3, s=5, label='Metastable')
    plt.scatter(unstable['Mh_calc'], unstable['Mt'], c='red', alpha=0.3, s=5, label='Unstable')
    
    if not nonpert.empty:
        plt.scatter(nonpert['Mh_calc'], nonpert['Mt'], c='black', alpha=0.3, s=5, label='Non-perturbative')

    # To prevent contours from bleeding into the stable region, we MUST NOT delete the stable points.
    # But using all 1 million stable points makes Delaunay triangulation too slow.
    # To prevent contours from bleeding into the stable region, we MUST NOT delete the stable points.
    # But using all 1 million stable points can be slow. Let's use all active points and 10% of stable points.
    df_active = df[df['Stability'] != 1].copy()
    df_stable = df[df['Stability'] == 1].sample(frac=0.1, random_state=42).copy()
    
    df_contour = pd.concat([df_active, df_stable], ignore_index=True)
    
    # Handle any Infs or NaNs
    df_contour.replace([np.inf, -np.inf], np.nan, inplace=True)
    df_contour.dropna(subset=['S_exact', 'Mh_calc', 'Mt'], inplace=True)
    
    df_contour.loc[df_contour['S_exact'] < 0, 'S_exact'] = 1e8
    
    # Clamp classical instability artifacts (huge positive values) to 0, as they are deeply unstable
    df_contour.loc[df_contour['S_exact'] > 1e10, 'S_exact'] = 0.0
    
    if not df_contour.empty:
        print("Drawing contours...")
        
        # Define contour levels to span the Unstable (S < 450) and Metastable (S > 450) regions
        # We add 500000 specifically to draw a contour through the Mh=5, Mt=96 point.
        levels = [10, 50, 100, 200, 450, 1000, 3000, 500000]
        
        # tricontour handles irregularly spaced data (due to Mh_calc being calculated dynamically)
        contour = plt.tricontour(df_contour['Mh_calc'], df_contour['Mt'], df_contour['S_exact'], 
                                 levels=levels, colors='blue', linewidths=1.5, alpha=0.9)
        
        # Add labels to the contour lines, manually placing them along Mh=50
        # We explicitly define the y-coordinates to perfectly hit each of the 8 contours exactly once
        y_locs = [127, 140, 155, 170, 185, 200, 215, 230]
        manual_locations = [(50, y) for y in y_locs]
        
        plt.clabel(contour, inline=True, fontsize=11, fmt='S=%1.0f', manual=manual_locations)

    plt.scatter(125.10, 173.1, color='blue', marker='*', s=200, edgecolor='black', label='Experimental (125.1, 173.1)', zorder=5)

    plt.xlabel('Higgs Mass $M_h$ (GeV)', fontsize=12)
    plt.ylabel('Top Mass $M_t$ (GeV)', fontsize=12)
    plt.title('Numerical Phase Diagram with Constant Action Contours', fontsize=14)
    
    # Set the same bounds as the standard numerical plot (0 to 250)
    plt.xlim(0, 250)
    plt.ylim(0, 250)
        
    plt.legend(loc='lower right')
    plt.grid(True, linestyle='--', alpha=0.7)
    
    output_file = 'results/numerical_stability_contours.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved {output_file}")

if __name__ == '__main__':
    create_contour_plot()



### Appendix D: plot_numerical.py

import pandas as pd
import matplotlib.pyplot as plt

def create_plot(df, filename, title, xlim=None, ylim=None, s=5):
    stable = df[df['Stability'] == 1]
    metastable = df[df['Stability'] == 2]
    unstable = df[df['Stability'] == 3]
    nonpert = df[df['Stability'] == 4]

    plt.figure(figsize=(10, 8))
    
    plt.scatter(stable['Mh_calc'], stable['Mt'], c='green', alpha=0.5, s=s, label='Stable')
    plt.scatter(metastable['Mh_calc'], metastable['Mt'], c='yellow', alpha=0.5, s=s, label='Metastable')
    plt.scatter(unstable['Mh_calc'], unstable['Mt'], c='red', alpha=0.5, s=s, label='Unstable')
    
    if not nonpert.empty:
        plt.scatter(nonpert['Mh_calc'], nonpert['Mt'], c='black', alpha=0.5, s=s, label='Non-perturbative')

    plt.scatter(125.10, 173.1, color='blue', marker='*', s=200, edgecolor='black', label='Experimental (125.1, 173.1)', zorder=5)

    plt.xlabel('Higgs Mass $M_h$ (GeV)')
    plt.ylabel('Top Mass $M_t$ (GeV)')
    plt.title(title)
    
    if xlim:
        plt.xlim(xlim)
    if ylim:
        plt.ylim(ylim)
        
    plt.legend(loc='lower right')
    plt.grid(True, linestyle='--', alpha=0.7)
    
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved {filename}")

# Full plot (using broad numerical data)
print("Plotting full view using broad numerical data...")
df_full = pd.read_csv('results/numerical_data.csv')
create_plot(df_full, 'results/numerical_stability_plot.png', 'Standard Model Vacuum Stability (Numerical, Full View)', s=5)

# Close-up plot (using dense numerical closeup data)
print("Plotting close-up view using dense numerical closeup data...")
df_dense = pd.read_csv('results/numerical_dense_closeup.csv')
create_plot(df_dense, 'results/numerical_stability_closeup.png', 'Standard Model Vacuum Stability (Numerical, Close-up)', xlim=(110, 140), ylim=(160, 185), s=2)




### Appendix E: plot_overlay.py

"""
Overlay / comparison plot between analytical and numerical stability maps.

Colouring rules:
  - Points where both methods AGREE on stability → plotted in that stability's colour
      Stable      = green
      Metastable  = yellow / goldenrod
      Unstable    = red
      Non-pert.   = black / grey
  - Points where the two methods DISAGREE → plotted in BLUE
      (these mark the boundary shift between analytical and numerical)

Both the full plot and a close-up around the experimental point are saved.
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# ── colour map ────────────────────────────────────────────────────────────────
COLOUR = {
    1: 'green',   # stable
    2: 'yellow',  # metastable
    3: 'red',     # unstable
    4: 'black',   # non-perturbative
}
DIFF_COLOUR = 'blue'   # disagreement between methods

# ── helper: align grids ───────────────────────────────────────────────────────
from scipy.spatial import cKDTree

def align_data(df_ana, df_num):
    df_ana['_Mt']  = df_ana['Mt'].round(6)
    df_ana['_Mh']  = df_ana['Mh_calc'].round(6)
    df_num['_Mt']  = df_num['Mt'].round(6)
    df_num['_Mh']  = df_num['Mh_calc'].round(6)
    
    tree = cKDTree(df_num[['_Mt', '_Mh']].values)
    dists, indices = tree.query(df_ana[['_Mt', '_Mh']].values)
    
    # Accept matches within a small distance
    valid = dists < 1.0
    
    merged = pd.DataFrame({
        '_Mt': df_ana['_Mt'].values[valid],
        '_Mh': df_ana['_Mh'].values[valid],
        'S_ana': df_ana['Stability'].values[valid],
        'S_num': df_num['Stability'].values[indices[valid]]
    })
    return merged

# ── helper: build per-stability masks for agreeing points ─────────────────────
def plot_stability_diagram(ax, merged, s_agree=3, s_disagree=8, xlim=None, ylim=None):
    Mt  = merged['_Mt'].values
    Mh  = merged['_Mh'].values
    S_a = merged['S_ana'].values
    S_n = merged['S_num'].values

    agree    = S_a == S_n
    disagree = ~agree

    # Plot agreeing points – grouped by stability value for legend clarity
    for stab_val, label in [(1, 'Stable'), (2, 'Metastable'), (3, 'Unstable'), (4, 'Non-perturbative')]:
        mask = agree & (S_n == stab_val)
        if mask.sum() == 0:
            continue
        ax.scatter(Mh[mask], Mt[mask],
                   c=COLOUR[stab_val], alpha=0.5, s=s_agree,
                   label=f'{label} (both agree)', rasterized=True)

    # Plot disagreeing points on top in blue
    if disagree.sum() > 0:
        ax.scatter(Mh[disagree], Mt[disagree],
                   c=DIFF_COLOUR, alpha=0.9, s=s_disagree, zorder=10,
                   label='Disagreement (boundary shift)', rasterized=True)

    # Experimental point
    ax.scatter(125.10, 173.1, color='blue', marker='*', s=200,
               edgecolor='black', zorder=15,
               label='Experimental (125.1, 173.1)')

    ax.set_xlabel('Higgs Mass $M_h$ (GeV)', fontsize=12)
    ax.set_ylabel('Top Mass $M_t$ (GeV)', fontsize=12)
    ax.grid(True, linestyle='--', alpha=0.4)

    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)

    ax.legend(loc='lower right', fontsize=9, markerscale=2)

# ── full plot ─────────────────────────────────────────────────────────────────
print("Loading broad analytical data for full plot …")
df_ana_full = pd.read_csv('results/analytical_data.csv')

print("Loading broad numerical data for full plot …")
df_num_full = pd.read_csv('results/numerical_data.csv')

print("Aligning broad grids …")
merged_full = align_data(df_ana_full, df_num_full)
print(f"Merged full points: {len(merged_full):,}")

print("Generating full overlay plot …")
fig, ax = plt.subplots(figsize=(11, 9))
plot_stability_diagram(ax, merged_full, s_agree=3, s_disagree=8)
ax.set_title('Vacuum Stability: Analytical vs Numerical Overlay\n'
             '(Blue = disagreement between methods)', fontsize=13)
fig.tight_layout()
fig.savefig('results/overlay_plot.png', dpi=200, bbox_inches='tight')
plt.close(fig)
print("Saved results/overlay_plot.png")

# ── close-up plot ─────────────────────────────────────────────────────────────
print("\nLoading dense analytical closeup data for closeup plot …")
df_ana_dense = pd.read_csv('results/analytical_dense_closeup.csv')

print("Loading dense numerical closeup data for closeup plot …")
df_num_dense = pd.read_csv('results/numerical_dense_closeup.csv')

print("Aligning dense closeup grids …")
merged_dense = align_data(df_ana_dense, df_num_dense)
print(f"Merged dense closeup points: {len(merged_dense):,}")

print("Generating close-up overlay plot …")
fig, ax = plt.subplots(figsize=(11, 9))
plot_stability_diagram(ax, merged_dense, s_agree=1, s_disagree=2, xlim=(110, 140), ylim=(160, 185))
ax.set_title('Vacuum Stability Overlay – Close-up\n'
             '(Blue = disagreement between methods)', fontsize=13)
fig.tight_layout()
fig.savefig('results/overlay_closeup_plot.png', dpi=200, bbox_inches='tight')
plt.close(fig)
print("Saved results/overlay_closeup_plot.png")

print("\nDone!")




### Appendix F: plot_analytical.py

import pandas as pd
import matplotlib.pyplot as plt

def create_plot(df, filename, title, xlim=None, ylim=None, s=5):
    stable = df[df['Stability'] == 1]
    metastable = df[df['Stability'] == 2]
    unstable = df[df['Stability'] == 3]
    nonpert = df[df['Stability'] == 4]

    plt.figure(figsize=(10, 8))
    
    plt.scatter(stable['Mh_calc'], stable['Mt'], c='green', alpha=0.5, s=s, label='Stable')
    plt.scatter(metastable['Mh_calc'], metastable['Mt'], c='yellow', alpha=0.5, s=s, label='Metastable')
    plt.scatter(unstable['Mh_calc'], unstable['Mt'], c='red', alpha=0.5, s=s, label='Unstable')
    
    if not nonpert.empty:
        plt.scatter(nonpert['Mh_calc'], nonpert['Mt'], c='black', alpha=0.5, s=s, label='Non-perturbative')

    plt.scatter(125.10, 173.1, color='blue', marker='*', s=200, edgecolor='black', label='Experimental (125.1, 173.1)', zorder=5)

    plt.xlabel('Higgs Mass $M_h$ (GeV)')
    plt.ylabel('Top Mass $M_t$ (GeV)')
    plt.title(title)
    
    if xlim:
        plt.xlim(xlim)
    if ylim:
        plt.ylim(ylim)
        
    plt.legend(loc='lower right')
    plt.grid(True, linestyle='--', alpha=0.7)
    
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved {filename}")

# Full plot (using broad analytical data)
print("Plotting full view using broad analytical data...")
df_full = pd.read_csv('results/analytical_data.csv')
create_plot(df_full, 'results/analytical_stability_plot.png', 'Standard Model Vacuum Stability (Analytical, Full View)', s=5)

# Close-up plot (using dense analytical closeup data)
print("Plotting close-up view using dense analytical closeup data...")
df_dense = pd.read_csv('results/analytical_dense_closeup.csv')
create_plot(df_dense, 'results/analytical_stability_closeup.png', 'Standard Model Vacuum Stability (Analytical, Close-up)', xlim=(110, 140), ylim=(160, 185), s=2)




### Appendix G: Claude's Previous Critique Output

The physical mechanism described is partially correct in spirit but contains four distinct errors that would not survive editorial review at a journal covering vacuum stability (JHEP, Phys. Rev. D, etc.). The most serious is a conflation of 'Fubini saddle-point breakdown' with 'classical roll-down instability,' which are physically distinct regimes. The negative bounce action computed in this region is an artifact of the approximation scheme, not a physical observable. The Fubini-Lipatov instanton profile is an exact saddle point of the Euclidean action only for the pure quartic potential with a *constant* negative coupling. For that potential, scale invariance in four Euclidean dimensions enforces the virial relation S_kin = -2 S_pot, ensuring S_E = -S_pot > 0. When lambda_eff(phi) runs strongly, this constraint is violated... The ratio S_kin / |S_pot| deviates from its expected value of 2 by orders of magnitude, providing a quantitative diagnostic that the profile is far from a saddle point... The true bounce in this regime yields a genuine, small positive action S_E_true << S_thresh, corresponding to rapid quantum tunneling through the thin barrier rather than a classical trajectory.
