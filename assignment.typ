#import "@preview/ctheorems:1.1.3": *
#import "@preview/plotst:0.2.0": *
#import "@preview/codly:1.2.0": *
#import "@preview/codly-languages:0.1.1": *
#codly(languages: codly-languages)

#show: codly-init.with()
#show: thmrules.with(qed-symbol: $square$)

#set heading(numbering: "1.1.")
#set page(numbering: "1")
#set heading(numbering: "1.")
#set math.equation(
  numbering: "1",
  supplement: none,
)
#show ref: it => {
  // provide custom reference for equations
  if it.element != none and it.element.func() == math.equation {
    // optional: wrap inside link, so whole label is linked
    link(it.target)[eq.~(#it)]
  } else {
    it
  }
}

#let theorem = thmbox("theorem", "Theorem", fill: rgb("#ffeeee")) //theorem color
#let corollary = thmplain(
  "corollary",
  "Corollary",
  base: "theorem",
  titlefmt: strong
)
#let definition = thmbox("definition", "Definition", inset: (x: 1.2em, top: 1em))

#set math.mat(delim: "[")
#set math.vec(delim: "[")
#let example = thmplain("example", "Example").with(numbering: "1.")
#let proof = thmproof("proof", "Proof")

//shortcuts
#let span = "span"

#align(center, text(21pt)[
  Assignment 2 - Numerical Linear Algebra
])

#align(center, text(17pt)[
  Arthur Rabello Oliveira

  #datetime.today().display("[day]/[month]/[year]")
])

#align(center)[
  #set par(justify: true)
  *Abstract*\
  We derive linear and polynomial regression in subsets of $RR$ and discuss the condition number of the associated matrices, numerical algorithms for the SVD and QR factorization are built and used on an efficiency analysis of the 3 methods to do linear or polynomial regression, stability of these algorirths is discussed.
]

#outline()

#pagebreak()

= Introduction

Given $D subset RR^2$, a dataset, approximating this set through a _continuous_ $f: RR -> RR$ is a very important problem in statistics, we will derive the 2 most important and most used methods to do this: linear and polynomial regression. Both are based on the least squares minimization problem. We will also discuss the conditioning number of the problems shown. A computational approach to regression is shown as well. We discuss how the condition number changes when the matrix is QR or SVD decomposed, and the algorithms for such decompositions are built.

= Norms and Problems

A #text(weight: "bold")[norm] is a way to quantify _size_ on a vector space. A norm is a class of functions that satisfy some properties, which we define below:

#definition[(Norm)
  Given $E$ a vector space over a field $KK$, a #text(weight: "bold")[norm] is a function $norm(dot):E -> RR$ that satisfies:

  - $norm(x) > 0_E, forall x in E^*$, and $norm(x) = 0_E <=> x = 0_E$\
  - $norm(x + y) <= norm(x) + norm(y)$
  - $norm(phi x) = abs(phi) norm(x), forall phi in KK$
]

Throughout this document we will use the most famous class of norms, the p-norms defined below:

#definition[(p-norm)
  Given $p in RR$, the #text(weight: "bold")[p-norm] of $x in CC^m$ is:

  $
    norm(x)_p = (sum_(i = 1)^m abs(x_i)^p)^(1/p)
  $
]

Some famous cases are:

$
  norm(x)_1 = sum_(i = 1)^m abs(x_i)\
  norm(x)_2 = (sum_(i = 1)^m abs(x_i)^2)^(1/2)\
  norm(x)_infinity = max_(1 <= i <= m) abs(x_i)
$

Now we proceed with the main topic of this section, the condition number of a problem. It is useful to see a problem as a function $f: X -> Y$ from a _normed_ vector space $X$ of data and a space $Y$ of solutions, $f$ is not always a well-behaved continuous function, which is why we are interested in #text(weight: "bold")[well-conditioned] problems and not in  #text(weight: "bold")[ill-conditioned] problems.

#definition[(Well-Conditioned Problem)
  A problem $f:X -> Y$ is _well-conditioned_ at $x_0 in X$ $<=> forall epsilon > 0, exists delta > 0 |  norm(x - x_0) < delta => f(x) - norm(f(x_0)) < epsilon$.
] <definition_of_stable_problem>

This means that small perturbations in $x$ lead to small changes in $f(x)$, a problem is #text(weight: "bold")[ill-conditioned] if $f(x)$ can suffer huge changes with small changes in $x$.

We usually say $f$ is well-conditioned if it is well-conditioned $forall x in X$. If there is at least one $x_i$ in which the problem is ill-conditioned, then we can use that whole problem is ill-conditioned.


== The Condition number of a problem

Condition numbers are a tool to quantify how well/ill conditioned a problem is:

#definition[(Absolute Conditioning Number)
  Let $delta x$ be a small pertubation of $x$, so $delta f = f(x + delta x) - f(x)$. The #text(weight: "bold")[absolute] conditioning number of $f$ is:

  $
    hat(kappa) = lim_(delta -> 0) sup_(norm(delta x) <= delta) (norm(delta f)) / (norm(delta x))
  $
] <definition_absolute_conditioning_number>

The limit of the supremum can be seen as the supremum of all _infinitesimal_ perturbations, so this can be rewritten as:

$
  hat(kappa) = sup_(delta x) (norm(delta f)) / (norm(delta x))
$

If $f$ is differentiable, we can evaluate the abs.conditioning number using its derivative. If $J$ is the matrix whose $i times j$ entry is the derivative $(diff f_i) / (diff x_j)$ (jacobian of $f$), then we know that $delta f approx J(x) delta x$, with equality in the limit $norm(delta x) -> 0$. So the absolute conditioning number of $f$ becomes:

$
  hat(kappa) = norm(J(x))
$ <absolute_conditioning_number__jacobian>

== The Relative Condition Number
When, instead of analyzing the whole set $X$ of data, we are interested in _relative_ changes, we use the #text(weight: "bold")[relative condition number:]

#definition[(Relative Condition Number)
  Given $f:X -> Y$ a problem, the _relative condition number_ $kappa(x)$ at $x in X$ is:

  $
    kappa(x) = lim_(delta -> 0) sup_(norm(delta x) <= delta) ((norm(delta f)) / (norm(f(x)))) dot ((norm(delta x)) / (norm(x)))^(-1)
  $
] <definition_relative_condition_number>

Or, as we did in @definition_absolute_conditioning_number, assuming that $delta f$ and $delta x$ are infinitesimal:

$
  kappa(x) = sup_(delta x) ((norm(delta f)) / (norm(f(x)))) dot ((norm(delta x)) / (norm(x)))^(-1)
$

If $f$ is differentiable:

$
  kappa(x) = (norm(J(x))) dot ((norm(f(x))) / (norm(x)))^(-1)
$ <relative_condition_number_through_jacobian>

Relative condition numbers are more useful than absolute conditioning numbers because the #text(weight: "bold")[floating point arithmetic] used in many computers produces _relative_ errors, the latter is not a highlight of this discussion.

Here are some examples of the definitions above:

#example[
  Consider the problem of obtaining the scalar $x/2$ from $x in RR$. The function $f(x) = x/2$ is differentiablle, so by @relative_condition_number_through_jacobian:

  $
    kappa(x) = (norm(J)) dot ((norm(f(x))) / (norm(x)))^(-1) = (1/2) dot ((x/2) / x)^(-1) = 1.
  $
]

This problem is well-conditioned ($kappa$ is small).

#example[
  Consider the problem of computing the scalar $x_1 - x_2$ from $(x_1, x_2) in RR^2$ (Use the $infinity$-norm in $RR^2$ for simplicity). The function associated is differentiable and the jacobian is:

  $
    J = mat((diff f) / (diff x_1) , (diff f) / (diff x_2)) = mat(1,-1)
  $
  With $norm(J)_infinity = 2$, so the condition number is:

  $
    kappa = (norm(J)_infinity) dot ((norm(f(x))) / (norm(x)))^(-1) = 2 / (|x_1 - x_2| dot max{|x_1|, |x_2|})
  $

  This problem can be ill-conditioned if $|x_1 - x_2| approx 0$ ($kappa$ gets huge), and well-conditioned otherwise
].

== Condition Number of Matrices

We will deduce the conditioning number of a matrix from the conditioning number of _matrix-vector_ multiplication:

Consider the problem of obtaining $A x$ given $A in CC^(m times n)$. We will calculate the relative condition number with respect to perturbations on $x$. Directly from @definition_relative_condition_number, we have:

$
  kappa = sup_(delta x) (norm(A (x + delta x) - A x)) / (norm(A x)) dot ((norm(delta x)) / (norm(x)))^(-1) = sup_(delta x) (norm(A delta x)) / (norm(delta x)) dot ((norm(A x)) / (norm(x)))^(-1)
$ <beginning_proof>

Since $sup_(forall x) (norm(A delta x)) / (norm(delta x)) = norm(A)$, we have:

$
  kappa = norm(A) dot (norm(x)) / (norm(A x))
$ <condition_number_of_matrix_vector_multiplication>

This is a precise formula as a function of $(A, x)$.

The following theorem will be useful in a near future:

#theorem[
  $forall x in CC^n, A in CC^(n times n), det(A) != 0$, the following holds:

  $
      (norm(x)) / (norm(A x)) <= norm(A^(-1))
  $
] <theorem_inequality_norms>

#proof[
  Since $norm(A B ) <= norm(A) norm(B)$, we have:

  $
    norm(A A^(-1) x) <= norm(A x) norm(A^(-1)) <=> (norm(x)) / (norm(A x)) <= norm(A^(-1))
  $
]

So using this in @condition_number_of_matrix_vector_multiplication, we can write:

$
  kappa <= norm(A) dot norm(A^(-1))
$

Or:

$
  kappa = alpha norm(A) dot norm(A^(-1))
$

With

$
  alpha = (norm(x)) / (norm(A x)) dot (norm(A^(-1)))^(-1)
$ <end_proof>

From @theorem_inequality_norms, we can choose $x$ to make $alpha = 1$, and therefore $kappa = norm(A) dot norm(A^(-1))$.

Consider now the problem of calculating $A^(-1) b$ given $A in CC^(n times n)$. This is mathematically identical to the problem we just analyzed, so the following theorem has already been proven:

#theorem[
  Let $A in CC^(n times n), det(A) != 0$, and consider the problem of computing $b$, from $A x = b$, by perturbating $x$. Then the following holds:

  $
    kappa = norm(A) (norm(x)) / (norm(b)) <= norm(A) dot norm(A^(-1))
  $
  
  Where $kappa$ is the condition number of the problem.
]

#proof[
  Read from @beginning_proof to @end_proof.
]

Finally, $norm(A) dot norm(A^(-1))$ is so useful it has a name: #text(weight: "bold")[the condition number of A] (relative to the norm $norm(dot)$)

If $A$ is singular, $kappa(A) = oo$. Notice that if $norm(dot) = norm(dot)_2$, then $norm(A) = sigma_1$ and $norm(A^(-1)) = 1/sigma_m$, so:

$
  kappa(A) = sigma_1 / sigma_m
$

This is the condition number of $A$ with respect to the $2$-norm, which is the most used norm in practice. The condition number of a matrix is a measure of how sensitive the solution of a system of equations is to perturbations in the data. A large condition number indicates that the matrix is ill-conditioned, meaning that small changes in the input can lead to large changes in the output.

= Linear Regression (1a)
<section_simple_linear_regression>

Given the dataset:

$
  D := {t_i = i/m}, i = 0, 1, dots, m in RR  
$ <dataset_linear_regression>

of equally spaced points , linear regression consists of finding the best _line_ $f(t) = alpha + beta t$ that approximates the points $(t_i, b_i) in RR^2$, where $b_i$ are arbitrary.

Approximating $2$ points in $RR^2$ by a line is trivial, now approximating more points is a task that requires linear algebra. To see this, we will analyze the following example to build intuition for the general case:

#figure(
  image("plots/least_squares_idea.png", width: 80%),
  caption: [
    A good approximation for the 3 points shown
  ],
) <figure_simples_linear_regression>

@figure_simples_linear_regression is a glimpse onto what we are about to produce.


Given the points $(1, 1), (2, 2), (3, 2) in RR^2$, we have $(t_1, b_1) = (1, 1), (t_2, b_2) = (2, 2), (t_3, b_3) = (3, 2) $ we would like a _line_ $f(t) = y(t) = alpha + beta t$ that best approximates $(t_i, b_i)$. In other words, since we know that the line does not pass through all 3 points, we would like to find the _closest_ line to #text(weight: "bold")[each point] of the dataset $D$. So the system:

$
  f(1) = alpha + beta = 1\
  f(2) = alpha + 2 beta = 2\
  f(3) = alpha + 3 beta = 2\
$

Which is:

$
  underbrace(mat(
    1, 1;
    1, 2;
    1, 3
  ), A) dot underbrace(mat(alpha;beta), x) = underbrace(mat(1;2;2), b)
$

Clearly has no solution. But it has a _closest solution_ which we can find through #text(weight: "bold")[minimizing] the errors produced by this approximation.

Let $x^* != x$ be a solution to the system. And let the error produced by approximating the points through a line be $e = A x - b$. Minimizing the error requires a _norm_, which is defined

$
  e_1^2 + e_2^2 + e_3^2
$

Is what we want to minimize, where $e_i$ is the error (distance) from the ith point to the line:

#figure(
  image("plots/least_squares_with_errors.png", width: 80%),
  caption: [
    The errors (distances)
  ],
)<figure_errors_simple_linear_regression>

So we will project $b$ into $C(A)$, giving us the closest solution, and the least squares solutions is when $hat(x)$ minimizes $norm(A x - b)^2$, this occurs when the residual $e = A x - b$ is orthogonal to $C(A)$. Since $N(A^*) perp C(A)$ and the dimensions sum up the left dimension of the matrix. so by the well-known projection formula, we have:

$
  A^* A hat(x) = A ^* b\

  = mat(
    1, 1, 1;
    1, 2, 3
  ) dot mat(
    1, 1;
    1, 2;
    1, 3
  ) dot mat(alpha;beta) = mat(
    3, 6;
    6, 14
  ) dot mat(alpha;beta)\

  = mat(
    3, 6;
    6, 14
  ) dot mat(1;2;3) = mat(5;11)
$

So the system to find $hat(x) = [hat(alpha), hat(beta)]$ becomes:

$ 
  3 alpha + 5 beta = 5\
  6 alpha + 14 beta = 11
$ <system_1>

Notice that with the _errors_ $e_i^2$ as:

$
  e_1^2 = (f(t_1) - b_1)^2 = (f(1) - 1)^2 = (alpha + beta - 1)^2\
  e_2^2 = (f(t_2) - b_2)^2 = (f(2) - 2)^2 = (alpha + 2 beta - 2)^2\
  e_3^2 = (f(t_3) - b_2)^2 = (f(3) - 2)^2 = (alpha + 3 beta - 2)^2
$

This is consistent with what we see in @figure_errors_simple_linear_regression. Notice that @system_1 is _precisely_ what is obtained after using partial derivatives to minimize the erros sum as a function of $(alpha, beta)$:

$
  f(alpha, beta) = (alpha + beta - 1)^2 + (alpha + 2 beta - 2)^2 + (alpha + 3 beta - 2)^2\
  = 3 alpha^2 + 14 beta^2 + 12 alpha beta - 10 alpha- 22 beta + 9,\

  (diff f) / (diff alpha) = (diff f) / (diff beta) = 0 <=> 6 alpha + 12 d - 10 = 28 beta + 12 alpha - 22 = 0 <=> cases(
    3c + 6d = 11, 6c + 14d = 11
  )
$ <system_2>

This new system has a solution in $hat(alpha) = 2/3 , hat(beta) = 1/2$, so the equation of the optimal line, obtained through _linear regression_ (or least squares) is:

$
  y(t) = 2/3 + 1/2 t.
$

If we have $n > 3$ points to approximate with linear regression, the reasoning is analogous:

Going back to $D$ in @dataset_linear_regression , we want to find the extended system as we did in @system_1, so let:
$
  f(t) = alpha + beta t
$

Be the linear regression line, for $D = {(0, b_0), (1/m, b_1), dots , (1, b_m)}$. The system is:

$
  f(0/m =   0) = b_0 = alpha,\
  f(1/m) = b_1 = alpha + beta/m,\
  f(2/m) = b_2 = alpha + 2/m beta\
  dots\
  f(m/m = 1) = b_m = alpha + beta
$

Or:

$
  underbrace(mat(
    1, 0;
    1, 1/m;
    dots.v, dots.v;
    1, 1
  ), A) dot underbrace(mat(alpha;beta), x) = underbrace(mat(b_0; dots.v; b_m), b)
$

Projecting into $C(A)$, we have:

$
  A^* A x = A^* b\
  = mat(1, 1, dots, 1; 0, 1/m, dots, 1) dot mat(
    1, 0;
    1, 1/m;
    dots.v, dots.v;
    1, 1
  ) = mat(
    m + 1, (m+1)/2;
    (m+1)/2 , ((m+1)(2m+2))/(6 m)
  ) dot mat(hat(alpha);hat(beta))\

  = mat(
    1, 1, dots, 1;
    0, 1/m, dots, 1
  ) dot mat(b_0; dots.v; b_m) = mat(
    b_0 + b_2 + dots + b_m;
    1/m [b_1 + 2 b_2 + dots + (m-1) b_(m-1) + b_m ]
  )
$

So the system to find the optimal vector $vec(hat(alpha), hat(beta))$ is:

$
  mat(
    m + 1, (m+1)/2;
    (m+1)/2 , ((m+1)(2m+2))/(6 m)
  ) dot mat(hat(alpha);hat(beta)) = mat(
    b_0 + b_1 + dots + b_m;
    1/m [b_1 + 2 b_2 + dots + (m-1) b_(m-1) + b_m ]
  )
$ <equation_with_AtransposeA>

Or, as a function of $t_i, b_i$ and $m$:

$
  underbrace(mat(
    m + 1, sum_(i = 1)^m t_i;
    sum_(i = 1)^m t_i, sum_(i = 1)^m t_i^2;
  ), hat(A)) dot mat(hat(alpha); hat(beta)) = mat(
    sum_(i = 1)^m b_i;
    sum_(i = 1)^m i/m dot b_i
  )
$ <matrix_system_least_squares>

This provides the optimal vector $hat(x)$ that minimizes the least squares error, which is the solution to the linear regression problem.

= How the condition number of A changes (1b)
<condition_number_linear_regression_system>

We are interested in the condition number of $hat(A) = A^* A$ shown in @equation_with_AtransposeA. We will analyze how the condition number of $hat(A)$ changes with respect to perturbations on $m$, the number of points in the dataset. A computational approach is appropriate.

Here is a python code that numerically calculates many values of $kappa(hat(A)_m)$ as a function of $m$:

```python
import numpy as np
import matplotlib.pyplot as plt

def cond_number(m):

    """
    This function computes the condition number of the matrix A(m) in the 2-norm. The matrix A is defined.

    Args:
        m (float): parameter for the matrix A(m)
    Returns:
        float: condition number of A(m)
    Raises:
        ZeroDivisionError: if m = 0
        np.linalg.LinAlgError: if A(m) is not invertible
    """

    A = np.array([
        [m + 1,          (m + 1) / 2],
        [(m + 1) / 2,  (m + 1)**2 / (3 * m)]
    ])
    A_inv = np.linalg.inv(A)
    return np.linalg.norm(A, 2) * np.linalg.norm(A_inv, 2)

M = float(input("Enter maximum m (M > 0): "))
N = int(input("Enter number of sample points: ")) #however the user wants to plot 

m_vals = np.linspace(0, M, N)
conds  = []

for m in m_vals:
    try:
        conds.append(cond_number(m))
    except (ZeroDivisionError, np.linalg.LinAlgError):
        conds.append(np.inf) #if it is not invertible

plt.figure()
plt.plot(m_vals, conds)
plt.xlabel('m')
plt.ylabel('Condition number $κ(A^* A)$')
plt.title('Condition number of $A^* A(m)$')
plt.grid(True)
plt.tight_layout()
plt.show()
```

Good visualizations of this are:

#figure(
  image("plots/condition_number.png", width: 80%),
  caption: [
    Condition number of $hat(A)$ over [0, 100]
  ],
) <condition_number_plot>

#figure(
  image("plots/condition_number_2.png", width: 80%),
  caption: [
    Condition number of $hat(A)$ over [0, 10000]
  ],
) <condition_number_plot_2>

@condition_number_plot and @condition_number_plot_2 show us that apparently $kappa(hat(A)_m)$ converges to a real number. We will evaluate this hypothesis below:

Using $norm(dot)_2$, the condition number of $hat(A)$ is:

$
  kappa(hat(A)) = norm(hat(A))_2 dot norm(hat(A)^(-1))_2 = sigma_1 / sigma_m 
$

Singular Values are better explored in @section_SVD_decompositon. Now we calculate the singular values of $hat(A)$, which are the square roots of the eigenvalues of $hat(A)$ (see @theorem_eigencalculation_of_svd). So we have:
$
  det(hat(A) - lambda I) = 0 <=> det(mat(
    m+1 - lambda , (m + 1) / 2;
    (m + 1) / 2 , ((m + 1) (2m + 1)) / (6m) - lambda
  )) = 0\ 

  <=> (m + 1 - lambda) [((m + 1) (2m + 1)) / (6m) - lambda] - ((m + 1)/2)^2 = 0\

  <=> lambda^2
    - ((m + 1)(8m + 1)) / (6m) lambda + ((m + 1)^2 (m + 2)) / (12m)
    = 0\

  <=> lambda = (m+1)/(12m) [(8m + 1) plus.minus sqrt(52m^2-  8m + 1)]
$

And the singular values are:

$
  sigma_1 = sqrt(lambda_1) = sqrt((m+1)/(12m) [(8m + 1) + sqrt(52m^2-  8m + 1)]),\

  sigma_2 = sqrt(lambda_2) = sqrt((m+1)/(12m) [(8m + 1) - sqrt(52m^2-  8m + 1)])
$

This gives:

$
  kappa(hat(A)) = sigma_1 / sigma_m = sqrt((m+1)/(12m) [(8m + 1) + sqrt(52m^2-  8m + 1)]) / sqrt((m+1)/(12m) [(8m + 1) - sqrt(52m^2-  8m + 1)])\

  = sqrt(((8m + 1) + sqrt(52m^2-  8m + 1)) / ((8m + 1) - sqrt(52m^2-  8m + 1)))
$ <condition_number_problem_b>

And the limit as $m -> oo$:

$
  lim_(m -> oo) kappa(hat(A)) = lim_(m -> oo) sqrt(((8m + 1) + sqrt(52m^2-  8m + 1)) / ((8m + 1) - sqrt(52m^2-  8m + 1)))\
$

Multiplying by the conjugate of the denominator and ignoring the square root (it is irelevant for the limit):

$
  = lim_(m -> oo) [ ((8m + 1) + sqrt(52m^2-  8m + 1)) / ((8m + 1) - sqrt(52m^2 - 8m + 1)) dot ((8m + 1) + sqrt(52m^2 - 8m + 1)) / ((8m + 1) + sqrt(52m^2 - 8m + 1))]\

  = lim_(m -> oo) (((8m + 1) + sqrt(52m^2 - 8m + 1))^2) / ((8m + 1)^2 - 52m^2 - 8m + 1)\

  = lim_(m -> oo) ((8m + 1)^2 + 2(8m + 1)sqrt(52m^2 - 8m + 1) + (52m^2 - 8m + 1)) / ((8m + 1)^2 - (52m^2 - 8m + 1))\

  = lim_(m -> oo) (64 m^2 + 16 m + 1 + (16m + 1) sqrt(52m^2 - 8m + 1) + 52 m^2 - 8 m + 1) / (64 m^2 + 16 m + 1 - 52m^2 + 8m - 1)\
$

Regretting having ignored the square root, and putting it back, we have:

$
  = lim_(m -> oo) sqrt((((8m + 1) + sqrt(52m^2 - 8m + 1))^2) / (12 m^2 + 24 m))\

  = lim_(m -> oo) ((8m + 1) + sqrt(52m^2 - 8m + 1)) / (sqrt(12 m^2 + 24 m))\

  = lim_(m -> oo) (8m + 1 + m sqrt(52 - 8/m + 1/m^2)) / (m sqrt(12 + 24/m))\

  = lim_(m -> oo) (m(8 + 1/m + sqrt(52 - 8/m + 1/m^2))) / (m sqrt(12 + 24/m))\

  = lim_(m -> oo) (8 + 1/m + sqrt(52 - 8/m + 1/m^2)) / (sqrt(12 + 24/m))\
$

And finally:

$
  lim_(m -> oo) kappa(hat(A)_m) = (8 + sqrt(52)) / sqrt(12) = (4 + sqrt(13)) / sqrt(3)  
$

A very good visualization of this is:

#figure(
  image("plots/desmos_condition_number.png", width: 60%),
  caption: [
    The purple line is the limit and the red is the function @condition_number_problem_b
  ],
) <condition_number_desmos_plot>


@condition_number_desmos_plot shows the function approaching the limit. One could say that this problem is well conditioned, for $kappa(hat(A))_m < (4 + sqrt(13)) / sqrt(3), forall m > 0$, and $(4 + sqrt(13)) / sqrt(3)$ is not a very big number. We will not go deep into the discussion of how well-condition this problem is, but we can say that the condition number of $A$ is not a problem for the linear regression algorithm.

= Polynomial Regression (1c)

In this section we will discuss what changes when we decide to use #text(weight: "bold")[polynomials]  instead of #text(weight: "bold")[lines] to approximate our dataset:

$
  f(t) = alpha + beta t -> p(t) = phi_0 + phi_1 t + dots + phi_n t^n
$

From a first perspective, it seems way more efficient to describe a dataset with many variables then to do so with a simple line $alpha + beta t$, so let's use the same dataset $S := {(t_i,b_i), t_i = i/m}, i = 0, 1, dots, m$. Where $b_i$ is arbitrary. As we did in @section_simple_linear_regression, finding the new system to be solved gives us:

$
  p(t_0 = 0) = b_0 = phi_0,\

  p(t_1 = 1/m) = b_1 = phi_0 + phi_1 1/m + dots + phi_n (1/m)^n\

  p(t_2 = 2/m) = b_2 = phi_0 + phi_1 2/m + phi_2 (2/m)^2 + dots + phi_n (2/m)^n\

  dots.v\
  
  p(t_m = 1) = b_m = phi_0 + dots + phi_n
$

Or:

$
 underbrace(mat(
  1, 0, 0, dots, 0;
  1, 1/m, (1/m)^2, dots, (1/m)^n;
  1, 2/m, (2/m)^2, dots, (2/m)^n;
  dots.v, dots.v, dots.v, dots.v, dots.v;
  1, 1, 1, dots, 1
 ),A_(m+1 times n+1)) dot underbrace(mat(
  phi_0; phi_1; phi_2; dots.v; phi_n
 ), Phi_(n+1 times 1)) = underbrace(mat(
  b_0; b_1; b_2; dots.v; b_m
 ), b_(m + 1 times 1))
$

Projecting into $C(A)$:

$
  A^* A hat(Phi) = mat(
    1, 1, 1, dots, 1;
    0, 1/m, 2/m, dots, 1;
    0, (1/m)^2, (2/m)^2, dots, 1;
    dots.v, dots.v, dots.v, dots, dots.v;
    0, (1/m)^n, (2/m)^n, dots, 1
  ) dot mat(
  1, 0, 0, dots, 0;
  1, 1/m, (1/m)^2, dots, (1/m)^n;
  1, 2/m, (2/m)^2, dots, (2/m)^n;
  dots.v, dots.v, dots.v, dots.v, dots.v;
  1, 1, 1, dots, 1
 ) dot mat(hat(phi_0); hat(phi_1); hat(phi_2); dots.v; hat(phi_n))\

 = mat(
  m + 1, sum_(i = 1)^m i/m, sum_(i = 1)^m (i/m)^2, dots, sum_(i = 1)^m (i/m)^n;
  sum_(i = 1)^m i/m, sum_(i = 1)^m (i/m)^2, sum_(i = 1)^m (i/m)^3, dots, sum_(i = 1)^m (i/m)^(n+1);
  sum_(i = 1)^m (i/m)^2, sum_(i = 1)^m (i/m)^3, sum_(i = 1)^m (i/m)^4, dots, sum_(i = 1)^m (i/m)^(n+2);
  dots.v, dots.v, dots.v, dots.v, dots.v;
  sum_(i = 1)^m (i/m)^n, sum_(i = 1)^m (i/m)^(n+1), sum_(i = 1)^m (i/m)^(n+2), dots, sum_(i = 1)^m (i/m)^(2n)
 ) dot mat(hat(phi_0); hat(phi_1); hat(phi_2); dots.v; hat(phi_n))\


 = mat(
    1, 1, 1, dots, 1;
    0, 1/m, 2/m, dots, 1;
    0, (1/m)^2, (2/m)^2, dots, 1;
    dots.v, dots.v, dots.v, dots, dots.v;
    0, (1/m)^n, (2/m)^n, dots, 1
  ) dot mat(
  b_0; b_1; b_2; dots.v; b_m 
  ) = mat(
    sum_(i = 0)^m b_i; sum_(i = 0)^m (i b_i) / m; sum_(i = 0)^m (i/m)^2 m; dots.v; sum_(i = 0)^m (i/m)^n b_i
  )
$

So the system to be solved is:

$
  underbrace(mat(
    m+ 1, sum_(i = 1)^m i/m, dots, sum_(i = 1)^m (i/m)^n;
    sum_(i = 1)^m i/m, sum_(i = 1)^m (i/m)^2, dots, sum_(i = 1)^m (i/m)^(n+1);
    sum_(i = 1)^m (i/m)^2, sum_(i = 1)^m (i/m)^3, dots, sum_(i = 1)^m (i/m)^(n+2);
    dots.v, dots.v, dots.v, dots.v;
    sum_(i = 1)^m (i/m)^n, sum_(i = 1)^m (i/m)^(n+1), dots, sum_(i = 1)^m (i/m)^(2n)
 ), hat(A)) dot mat(hat(phi_0); hat(phi_1); hat(phi_2); dots.v; hat(phi_n)) = mat(
    sum_(i = 0)^m b_i; sum_(i = 0)^m (i b_i) / m; sum_(i = 0)^m (i/m)^2 m; dots.v; sum_(i = 0)^m (i/m)^n b_i
  )
$ <matrix_system_polynomial_least_squares>

This gives the optimal vector $hat(Phi)$ that solves the least squares problem. We will use computational methods to analyze this system in some of the next sections.

= Computing the polynomial regression matrix, given (m,n) (1d)
<poly_ls_section>

Here is a python function that calculates the polynomial regression matrix $hat(A) = A^* A$ from @matrix_system_polynomial_least_squares, given the dimensions $(m, n)$:

```python
import numpy as np

def poly_ls(m, n):

    """
    Builds the (n+1) x (n+1) matrix A^T A for least-squares polynomial fitting.
    
    Args:
        m (int): number of subintervals (m >= 0)
        n (int): polynomial degree (n >= 0)
    Returns:
        np.ndarray: shape (n+1, n+1) Gram matrix
    Raises:
        ValueError: if m or n is negative or not integer
    """

    if not isinstance(m, int) or not isinstance(n, int):
        raise ValueError("m and n must be integers")
    if m < 0 or n < 0:
        raise ValueError("m and n must be non-negative")

    x = np.linspace(0, 1, m+1) #sample space

    A = np.zeros((n+1, n+1), dtype=float) #intializes 0 matrix to be filled
    np.set_printoptions(precision=3, suppress=True)
    for j in range(n+1): #THIS IS NOT A, IT IS A^* A
        for k in range(n+1):
            A[j, k] = np.sum(x**(j + k)) #fills each entry

    return A

for m, n in [(1, 1), (2, 2), (2, 3)]: #trivial examples
    M = poly_ls(m, n)
    print(f"m = {m}, n = {n}:")
    print(M, end="\n\n")

```

Some simple cases are:

$ 
  hat(A)(1, 1) = mat(
    2, 1;
    1, 1
  )\

  hat(A)(2, 2) = mat(
    3, 1.5, 1.25;
    1.5,   1.25,  1.125;
    1.25,  1.125, 1.062
  )\

  hat(A)(2, 3) = mat(
    3, 1.5, 1.25, 1.125;
    1.5,   1.25,  1.125, 1.062;
    1.25,  1.125, 1.062, 1.031;
    1.125, 1.062, 1.031, 1.016
  )
$

= How Perturbations Affect The Condition Number (1e)
<section_perturbations_polynomial_regression>

Still on polynomial regression, in this section we analyze what happens to $kappa(hat(A))$, when $hat(A)$ is perturbated with $m = 100$ and $n = 1, dots, 20$.

We will run _poly_ls(m, n)_ built in @poly_ls_section for $m = 100$ and $n = 1, dots, 20$. Given the sad fact that polynomial rootfinding is not a well-conditioned problem, we will then numerically calculate the condition number of $hat(A)$ for each output of _poly_ls(m, n)_. See the code:


```python
import numpy as np
import matplotlib.pyplot as plt

def format_scientific(x, sig=3):

    """
    Formats a number in scientific notation with a specified number of significant digits.

    Args:
        x (float): number to format
        sig (int): number of significant digits (default: 3)
    Returns:
        str: formatted string in scientific notation
    """

    if x == 0:
        return "0"
    exp = int(np.floor(np.log10(abs(x))))
    mant = x / 10**exp
    return f"{mant:.{sig}f} * 10^{exp}"

def compute_condition_numbers(m, max_n):

    """
    Returns a list of the condition numbers of the polynomial least-squares matrix A(m) for degrees n = 1 to max_n.

    Args:
        m (int): number of subintervals (m >= 0)
        max_n (int): maximum polynomial degree (max_n >= 0)
    Returns:
        list: condition numbers of A(m) for degrees n = 1 to max_n
    """

    conds = []
    for n in range(1, max_n + 1):
        A = poly_ls(m, n) 
        sv = np.linalg.svd(A, compute_uv=False) #computes singular values
        conds.append(sv[0] / sv[-1]) #condition number is the ratio of the largest to smallest singular value.
    return conds

if __name__ == "__main__":
    m = 100
    max_n = 20

    cond_nums = compute_condition_numbers(m, max_n)
    n_values = np.arange(1, max_n + 1)

    print(f"Condition numbers of A (m={m}) for degree n:")
    for n, c in zip(n_values, cond_nums):
        print(f"  n = {n:2d} → κ₂(A) = {format_scientific(c)}")

    plt.figure()
    plt.semilogy(n_values, cond_nums, marker="o", linestyle="-")
    plt.xlabel("Polynomial degree $n$")
    plt.ylabel("Condition number $\\kappa(A)$")
    plt.title(f"Growth of Condition Number, $m={m}$")
    plt.grid(True, which="both", ls="--")
    plt.tight_layout()
    plt.show()
```

A good plot of the growth of the condition number is:

#figure(
  image("plots/(m,n)condition_perturbation.png", width: 80%),
  caption: [
    Growth of the condition number of A(m) for polynomial regression
  ],
) <condition_number_perturbation_plot>

@condition_number_perturbation_plot Shows that with higher degrees of polynomials, $kappa(hat(A)_m)$ grows #text(weight: "bold")[exponentially]. This can make polynomial regression a #text(weight: "bold")[bad choice] for approximating datasets with many points. Which contradicts our initial assumption! 


One could notice that the graph wiggles up and down at $12 <= n <= 20$. To understand this we must take into account that $epsilon_m approx 2.22 dot 10^(-16)$, and at $n = 12$, $kappa(hat(A))$ is already past $10^16$. Notice as well that we are calculating on $hat(A)$, not $A$, so:

$
  kappa(hat(A) = A^* A) = sigma_max(A^* A) / sigma_min(A^* A) = (sigma_max(A) / sigma_min(A))^2 = kappa(A)^2
$ <equation_weird_condition_plot>

Squaring $kappa(A)$ pushes numbers past the #link("https://en.wikipedia.org/wiki/Double-precision_floating-point_format")[IEEE 754 double precision] quickly. So once the trye $kappa(hat(A) = A^* A)$ exceeds $1/epsilon_m$, the smallest singular value underflows to $0$. So everything above that is #text(weight: "bold")[numerically indistinguishable]. So the line stops at $approx 10^16$ and wiggles for a while.

= Polynomial Regression with a Different Dataset 
== A Different Dataset

If we change $S := {(t_i, b_i) | t_i = i / m, i = 0, 1, dots, m}$ to $hat(S) = {(t_i, b_i) | t_i = i/m - 1/2}$, the polynomial regression becomes:

$
  p(t_0 = 0 - 1/2) = phi_0 + phi_1 (-1/2) + dots + phi_n (-1/2)^n = b_0\
  p(t_1 = 1/m - 1/2) = phi_0 + phi_1 (1/m - 1/2) + phi_2 (1/m - 1/2)^2 + dots + phi_n (1/m - 1/2)^n\
  dots.v\
  p(t_m = 1 - 1/2) = phi_0 + phi_1 (1 - 1/2) + dots + phi_n (1 - 1/2)^n
$

So:

$
  underbrace(mat(
    1, -1/2, (-1/2)^2, dots, (-1/2)^n;
    1, (1/m - 1/2), (1/m - 1/2)^2, dots, (1/m - 1/2)^n;
    1, (2/m - 1/2), (2/m - 1/2)^2, dots, (2/m - 1/2)^n;
    dots.v, dots.v, dots.v, dots.v, dots.v;
    1, (-1/2), (-1/2)^2, dots, (-1/2)^n
  ), A) dot underbrace(mat(
    phi_0; phi_1; dots.v; phi_n
  ), Phi) = underbrace(mat(
    b_0; b_1; dots.v; b_m
  )
  , b)
$

Projecting onto $C(A)$:

$
  underbrace(mat(
    1, 1, 1, dots, 1;
    -1/2, (1/m - 1/2), (2/m - 1/2), dots, -1/2;
    (-1/2)^2, (1/m, - 1/2)^2, (2/m - 1/2)^2, dots, (-1/2)^2;
    dots.v, dots.v, dots.v, dots.v, dots.v;
    (-1/2)^n, (1/m - 1/2)^n , (2/m - 1/2)^n, dots, (-1/2)^n
  ), A^*) dot underbrace(mat(
    1, -1/2, (-1/2)^2, dots, (-1/2)^n;
    1, (1/m - 1/2), (1/m - 1/2)^2, dots, (1/m - 1/2)^n;
    1, (2/m - 1/2), (2/m - 1/2)^2, dots, (2/m - 1/2)^n;
    dots.v, dots.v, dots.v, dots.v, dots.v;
    1, -1/2, (-1/2)^2, dots, (-1/2)^n
  ), A) dot underbrace(mat(
    hat(phi_0); hat(phi_1); dots.v; hat(phi_n)
  ), hat(Phi))

$

Notice that to calculate $A^* A$ we can do:

$
  (A^* A)_(i j) = angle.l l_i^(A^*), c_j^A angle.r = angle.l c_i^A, c_j^A angle.r = sum_(k = 0)^m (k/m - 1/2)^(i + j - 2)
$

So we have:

$
  mat(
    n + 1, sum_(i = 0)^n (i/m - 1/2), sum_(i = 0)^n (i/m - 1/2)^2, dots,  sum_(i = 0)^n (i/m - 1/2)^n;
    sum_(i = 0)^n i/m - 1/2,sum_(i = 0)^n (i/m - 1/2)^2, sum_(i = 0)^n (i/m - 1/2)^3, dots, sum_(i = 0)^n (i/m - 1/2)^(n+1);
    dots.v, dots.v, dots.v, dots.v, dots.v;
    sum_(i = 0)^n (i/m - 1/2)^n, sum_(i = 0)^n (i/m - 1/2)^(n+1), sum_(i = 0)^n (i/m - 1/2)^(n+2), dots, sum_(i = 0)^n (i/m - 1/2)^(2n)
  ) dot mat(
    hat(phi_0);
    hat(phi_1);
    dots.v;
    hat(phi_n)
  )
$

And doing $A^* b$ gives:

$

  = mat(
    1, 1, 1, dots, 1;
    -1/2, (1/m - 1/2), (2/m - 1/2), dots, -1/2;
    (-1/2)^2, (1/m, - 1/2)^2, (2/m - 1/2)^2, dots, (-1/2)^2;
    dots.v, dots.v, dots.v, dots.v, dots.v;
    (-1/2)^n, (1/m - 1/2)^n , (2/m - 1/2)^n, dots, (-1/2)^n
  ) dot mat(
    b_0; b_1; b_2; dots.v; b_m
  ) = mat(
    sum_(i = 0)^n b_i;
    sum_(i = 0)^n (i/m - 1/2) b_i;
    sum_(i = 0)^n (i/m - 1/2)^2 b_i;
    dots.v;
    sum_(i = 0)^n (i/m - 1/2)^n b_i
  )
$

So the system to be solved is:

$
  underbrace(mat(
    m + 1, sum_(i = 0)^n (i/m - 1/2), dots,  sum_(i = 0)^n (i/m - 1/2)^n;
    sum_(i = 0)^n i/m - 1/2, sum_(i = 0)^n (i/m - 1/2)^2,dots, sum_(i = 0)^n (i/m - 1/2)^(n+1);
    dots.v, dots.v, dots.v, dots.v;
    sum_(i = 0)^n (i/m - 1/2)^n, sum_(i = 0)^n (i/m - 1/2)^(n+1), dots, sum_(i = 0)^n (i/m - 1/2)^(2n)
  ), hat(A)) dot vec(hat(phi_0), hat(phi_1), hat(phi_2), dots.v, hat(phi_n)) = dot mat(
    sum_(i = 0)^n b_i;
    sum_(i = 0)^n (i/m - 1/2) b_i;
    sum_(i = 0)^n (i/m - 1/2)^2 b_i;
    dots.v;
    sum_(i = 0)^n (i/m - 1/2)^n b_i
  )
$ <new_dataset_polynomial_regression>

The following code calculates the new matrix $hat(A)$ in @new_dataset_polynomial_regression:

```python
import numpy as np
import matplotlib.pyplot as plt

def poly_ls_2(m, n):

    """
    Builds the (n+1) x (n+1) matrix for least-squares polynomial fitting.

    Args:
        m (int): number of subintervals (m >= 0)
        n (int): polynomial degree (n >= 0)
    Returns:
        np.ndarray: shape (n+1, n+1) Gram matrix
    Raises:
        ValueError: if m or n is negative or not integer
    """

    if not (isinstance(m, int) and isinstance(n, int)) or m < 0 or n < 0:
        raise ValueError("m and n must be non-negative integers")

    t = np.linspace(0, 1, m + 1) - 0.5 
    powers = t[:, None] ** np.arange(2 * n + 1) 
    col_sums = powers.sum(axis=0)
    M = np.empty((n + 1, n + 1))
    for i in range(n + 1):
        for j in range(n + 1):
            M[i, j] = col_sums[i + j] #fills each entry

    return M

#examples:
m_1 = poly_ls_2(2, 1)
m_2 = poly_ls_2(2, 2)
m_3 = poly_ls_2(2, 3)
print("m = 2, n = 1:")
print(m_1)
print("\nm = 2, n = 2:")
print(m_2)
print("\nm = 2, n = 3:")
print(m_3)
```

The examples are:

#example[
  $
    M(2,1) = mat(
      3, 0;
      0, 0.5
    )
  $
]

#example[
  $
    M(2,2) = mat(
      3, 0, 0.5;
      0, 0.5, 0;
      0.5, 0, 0.125
    )
  $
]

#example[
  $
    M(2,3) = mat(
      3, 0, 0.5, 0;
      0, 0.5, 0, 0.125;
      0.5, 0, 0.125, 0;
      0, 0.125, 0, 0.031
    )
  $
]


== How Conditioning changes (1f)

Here we will analyze how the condition number of $hat(A)$ shown in the previous section changes with perturbations on the degree $n$. We will use the same method used in @section_perturbations_polynomial_regression. $m = 100$ and $n = 1, dots, 20$. The following code is used:

```python
import numpy as np
import matplotlib.pyplot as plt

def compute_condition_numbers_centered(m: int, max_n: int):

    """
    Computes the condition numbers of the polynomial least-squares matrix M(m) for degrees n = 1 to max_n.

    Args:
        m (int): number of subintervals (m >= 0)
        max_n (int): maximum polynomial degree (max_n >= 0)
    Returns:
        list: condition numbers of M(m) for degrees n = 1 to max_n
    """

    conds = []
    for n in range(1, max_n + 1):
        M  = poly_ls_2(m, n)
        s  = np.linalg.svd(M, compute_uv=False) #computes singular values
        conds.append(s[0] / s[-1]) #κ = σ_max / σ_min
    return conds

m, max_n = 100, 20

cond_nums = compute_condition_numbers_centered(m, max_n)
n_values  = np.arange(1, max_n + 1)

print(f"Condition numbers at (m = {m})")
for n, κ in zip(n_values, cond_nums):
    print(f"  n = {n:2d} → κ(G) = {format_scientific(κ)}")

plt.figure()
plt.semilogy(n_values, cond_nums, marker="o")
plt.xlabel("Polynomial degree $n$")
plt.ylabel(r"Condition number $\kappa_2(G)$")
plt.title(fr"Growth of $\kappa$, $m={m}$")
plt.grid(True, which="both", ls="--")
plt.tight_layout()
plt.show()
```

The expected output is:

#figure(
  image("plots/growth_of_condition_number_polynomial_regression_new_dataset.png", width: 80%),
  caption: [
    Growth of the condition number of $hat(A)$ with a new dataset
  ]
) <growth_of_condition_number_polynomial_regression_new_dataset_plot>


The weird behavior noticed in @condition_number_perturbation_plot is not seen in @growth_of_condition_number_polynomial_regression_new_dataset_plot. This is due to the fact that the dataset has been scaled to be centered at $0$, reducing the dynamic range in each column and preventing the singular values from the underflow zone (at least for $20 >= n$)

= Comparing the Condition Number
<section_graphically_comparing_condition_numbers>

Here we graphically compare both condition numbers seen in @condition_number_perturbation_plot and @growth_of_condition_number_polynomial_regression_new_dataset_plot. The following code is used:

```python 
def compute_condition_numbers(m, max_n):

    """
    Compute k(A^* A).

    Each matrix A is constructed by
    poly_ls(m, n)

    Args:
        m (int): Number of subintervals in the sample grid t_i = i/m.
        max_n (int): Maximum polynomial degree to evaluate.

    Returns:
        list[float]: The list  
        [ k(A_i^* A_i)) ],  
        where k(A^* A) = k(A)^2 and κ is the 2-norm condition number.
    """

    conds = []
    for n in range(1, max_n + 1):
        A = poly_ls(m, n)
        sigma = np.linalg.svd(A, compute_uv=False)
        κ2 = sigma[0] / sigma[-1]  #k_2(A)
        conds.append(κ2 ** 2) #κ_2(A^* A) = κ_2(A)^2
    return conds


def compute_condition_numbers_centered(m, max_n):

    """
    Compute k(M) for a centred / scaled polynomial basis.
    Each matrix M is produced by poly_ls_2(m, n)

    Args:
        m (int): Number of subintervals in the sample grid t_i = i/m.
        max_n (int): Maximum polynomial degree to evaluate.

    Returns:
        list[float]: The list  
        [ k(M_i) ] where  
        k is the 2-norm condition number of the centred design matrix.
    """

    conds = []
    for n in range(1, max_n + 1):
        M = poly_ls_2(m, n)
        sigma = np.linalg.svd(M, compute_uv=False)
        conds.append(sigma[0] / sigma[-1])
    return conds


def plot_condition_numbers(m = 100, max_n = 20):

    """
    Plot the growth of condition numbers for raw and centred Vandermonde bases.

    Args:
        m (int, optional): Number of subintervals in the sample grid
            t_i = i/m. Defaults to 100.

        max_n (int, optional): Maximum polynomial degree to display.
            Defaults to 20.

    Returns:
        None.  The function displays a semilog plot comparing  
        k(A^* A) (raw basis) and k(M) (centred / scaled basis).
    """

    n_vals = np.arange(1, max_n + 1)
    κ_raw = compute_condition_numbers(m, max_n)
    κ_centered = compute_condition_numbers_centered(m, max_n)

    plt.figure(figsize=(7, 5))
    plt.semilogy(
        n_vals,
        κ_raw,
        marker="o",
        linestyle="-",
        linewidth=1.4,
        markersize=5,
        label=r"$\kappa(A^{\top}A)$ (raw basis)",
    )
    plt.semilogy(
        n_vals,
        κ_centered,
        marker="s",
        linestyle="--",
        linewidth=1.4,
        markersize=5,
        label=r"$\kappa(M)$ (centred / scaled basis)",
    )

    plt.xlabel("Polynomial degree $n$")
    plt.ylabel("Condition number")
    plt.title(fr"Condition-number growth, $m={m}$")
    plt.grid(True, which="both", ls=":", lw=0.7)
    plt.legend()
    plt.tight_layout()
    plt.show()

plot_condition_numbers(m=100, max_n=20)
```

The expected output is the plot below:

#figure(
  image("plots/graphical_comparison_condition_numbers.png", width: 80%),
  caption: [
    A comparison of both datasets
  ]
) <figure_comparison_datasets>

@figure_comparison_datasets confirms the hypothesis that the condition number of the non-centered dataset tends to produce bigger condition numbers.

= Least Squares with QR and SVD decompositions

We have shown the solutions to the least squares problem $A x = b$, but this problem could be solved with factorizations of $A$, such as the QR and SVD, in the following sections we will show these factorizations and use them to solve the least squares problem.


== QR
<section_QR_decomposition>


The QR factorization of a full-rank $A in CC^(m times n)$, $m >= n$ an consists of finding orthonormal vectors $q_1, dots, q_n$ such that $q_1, dots, q_i$ spans $a_1, dots, q_1$, where $a_i$ is the ith-column of $A$. So we want:

$
  span(a_1) = span(q_1)\
  span(a_1, a_2) = span(q_1, q_2)\
  dots.v\
  span(a_1, dots, a_n) = span(q_1, dots, q_n)
$

This is equivalent to:

$
  A = mat(
    ,,;
    q_1, dots, q_n;
    ,,
  ) dot mat(
    r_11, r_12, dots, r_(1 n);
    ,r_22, dots, r_(2 n);
    , , dots.v, dots.v;
    ,,,r_(n n)
  )
$ <equation_open_reduced_QR>

Where $r_(i i) != 0$, because $a_i$ will be expressed as a linear combination of $q_i$, and since the triangular matrix is invertible, $q_i$ can be expressed as a linear combination of $a_i$. Therefore @equation_open_reduced_QR is:

$
  a_1 = q_1 r_11,\
  a_2 = r_12 q_1 + r_22 q_2,\
  dots.v\
  a_n = r_(1 n) q_1 + r_(2 n) q_2 + dots + r_(n n) q_n.
$

Or:

$
  A = hat(Q) hat(R)
$ <equation_closed_reduced_QR>

Is the _reduced_ QR decomposition of $A$.

The _full_ QR decomposition of $A in CC^(m times n)$ not of full-rank is analogous to the reduced, but $|m - n|$ 0-columns are appended to $hat(Q)$ to make it a unitary $m times m$ matrix $Q$, and 0-rows are aded to $hat(R)$ to make it a $m times n$ still triangular matrix:

#figure(
  image("plots/full_qr.png", width: 60%),
  caption: [
    Full QR factorization
  ],
)
<figure_full_qr>

And the decomposition becomes:

$
  A = Q R
$

Here are some examples:

  #example[
    $
      A = mat(
      3, 0, 0;
      0, 4, 0;
      0, 0, 5
      )
    $

    This is a diagonal matrix, so its QR factorization is particularly simple:
    
    $
      Q = mat(
      1, 0, 0;
      0, 1, 0;
      0, 0, 1
      ),

      R = mat(
      3, 0, 0;
      0, 4, 0;
      0, 0, 5
      )
    $

    With diagonal matrices, $Q$ is the identity matrix and $R = A$.
  ]

  #example[
    $
      A = mat(
      1, 1;
      1, 0;
      0, 1
      )
    $

    For this $3 times 2$ matrix, we compute the reduced QR factorization:

    $
      hat(Q) = mat(
      1/sqrt(2), 0;
      1/sqrt(2), -1/sqrt(2);
      0, 1/sqrt(2)
      ),

      hat(R) = mat(
      sqrt(2), 1/sqrt(2);
      0, sqrt(2)/2
      )
    $

    This is a reduced QR factorization where $hat(Q)$ is $3 times 2$. The full QR factorization would require extending $hat(Q)$ to a $3 times 3$ orthogonal matrix and adding a row of zeros to $hat(R)$ as shown in @figure_full_qr.
  ]

== SVD
<section_SVD_decompositon>

The _singular value decomposition_  of a matrix is based on the fact that the image of the unit sphere under a $m times n$ matrix is a #text(weight: "bold")[hyperellipse:]

#figure(
  image("plots/hyperellipse.png", width: 80%),
  caption: [
    SVD of a $2 times 2$ matrix
  ],
)

So the independent directions $v_1, v_2$ have been mapped to another set of orthogonal directions $sigma_1 v_1, sigma_2 v_2$, so with $S:= {v in CC^n | norm(v) = 1}$ as the unit ball,  let's define:

#definition[(Singular Values)
  The $n$ _singular values_ $sigma_i$ of $A in CC(m times n)$ are the lengths of the $n$ new axes of $A S$, written in non-crescent order $sigma_1 >= dots >= sigma_n$.
]<definition_singular_values>

#definition[(Left Singular Vectors)
  The $n$ #text(weight: "bold")[left] singular vectors of $A$ are the unit vectors $u_i$ laying in $A S$, oriented to correspond and  number  the singular values $sigma_i$, respectively
]<definition_left_singular_vectors>

#definition[(Right Singular Vectors)
  The #text(weight: "bold")[right] singular vectors of $A$ are the $v_i$ in $S$ that are the preimages of $sigma_i u_i in A S$, such that $A v_i = sigma_i u_i$
]<definition_right_singular_vectors>

The equation $A v_i = sigma_i u_i$ is equivalent to:

$
 A dot mat(
    ,,,;
    v_1, v_2, dots, v_n;
    ,
  ) = mat(
    ,,,;
    sigma_1 u_1, sigma_2 u_2, dots, sigma_n u_n;
    ,
  )
$

Better:

$
  A dot mat(
    ,,,;
    v_1, v_2, dots, v_n;
    ,
  ) = mat(
    ,,,;
    u_1, u_2, dots, u_n;
    ,
  ) dot mat(
    sigma_1,0 ,dots,0;
    0, sigma_2, dots, 0;
    dots.v, dots.v, dots.v, dots.v;
    0, 0, dots, sigma_n
  )
$

Or simple $A V = U Sigma$, but since $V$ has orthonormal columns:

$
  A = U Sigma V^*
$<SVD_equation>


The SVD is a very particular factorization for matrices, as the following theorem states:

#theorem[(Existence of SVD)
  _Every_ matrix $A in CC^(m times n)$ has a singular value decomposition
]

#proof[
  We prove the existane by fixing the largest image of $A$ and using induction on the dimension of $A$:

  Let $sigma_1 = norm(A)_2$. There must exist unitary vectors $u_1, v_1 in CC^n$ such that $A v_1 = sigma_1 u_1$, with $norm(v_1)_2 = norm(u_1)_2 = 1.$ Let ${v_j}$ and ${u_j}$ be 2 orthonormal bases of $CC^n$. These column vectors form the unitary matrices $V_1$ and $U_1$. We will compute:

  $
    Phi = U_1^* A V_1
  $ <equation_Phi_Ut_A_V_matrix>

  Notice that the first column of $Phi$ is $U_1^* A v_1 = sigma_1 U_1^* v_1 = sigma_1 e_1$, since $u_1$ is the first column of $U_1$. So $Phi$ looks like:

  $
    Phi = mat(
      sigma_1, w^*;
      0, B
    )
  $

  Where $w^*$ is the rest of the first row, the action of $A$ onto the remaining columns $v_j$. $B$ acts on the subspace orthogonal to $v_1$.

  We want $w = 0$, we can force this by using the norm. We know that:

  $
    norm(mat(
      sigma_1, w^*;
      0, B
    ) dot vec(sigma_1, w))_2 = norm(mat(
      sigma_1^2 +  w^* w;
      B w; 
    ))_2 = sqrt(|sigma_1^2 + w^* w|^2 + norm(B w)_2^2)  
  $

  And:

  $
    sqrt(|sigma_1^2 + w^* w|^2 + norm(B w)_2^2) >= sigma_1^2 + w^* w
  $

  We also know:
  
  $
   norm(Phi)_2 = sup_(norm(y) = 1) norm(Phi y)_2
  $
  
  For the specific $x = [sigma_1, w]$ scaled to the unit ball, and knowing $norm(Phi)_2 = sigma_1$, we have:

  $
    norm(Phi)_2 >= (norm(Phi x)_2) / (norm(x)_2) >= (sigma_1^2 + w^* w) / sqrt(sigma_1^2 + w^* w) = sqrt(sigma_1^2 + w^* w) <=> sigma_1 >= sqrt(sigma_1^2 + w^* w)\

    <=> sigma_1^2 >= sigma_1^2 + w^* w <=> w^* w = 0 <=> w = 0.
  $

  If $m = 1$ or $n = 1$, we are done, If not, $B$ has an SVD decomposition $B = U_2 Sigma_2 V_2^*$ by the induction hypothesis, so from @equation_Phi_Ut_A_V_matrix we have that the following is a SVD decomposition of $A$, completing the proof:

  $
    A = U_1 mat(1, 0; 0, U_2) mat(sigma_1, 0; 0, Sigma_2) mat(1, 0; 0, V_2)^* V_1^*
  $
]

Is the SVD factorization of A. There are more about the SVD on computing $U, Sigma, V^*$, as we will show below:

#theorem[
  $forall A in CC^(m times n)$, the following holds:
  
  - The eigenvalues of $A^*A$ are the singular values _squared_ of $A$, and the column-eigenvectors of $A^*A$ form the matrix $V$.

  - The eigenvalues of $A A^*$ are the singular values _squared_ of $A$, and the column-eigenvectors of $A A^*$ form the matrix $U$.
]<theorem_eigencalculation_of_svd>

#proof[
  Let $U Sigma V^* = A$ be the SVD of $A$, then computing $A^* A$, knowing $U,V$ are unitary matrices, we have:

  $
    A^* A = (U Sigma V^*)^* (U Sigma V^*) = V Sigma^* U^* U Sigma V^* = V Sigma^* Sigma V^* = V Sigma^2 V^*
  $

  This is an _eigenvalue_ decomposition of $A^* A$, where the eigenvalues are the entries of $Sigma^2$, which are the singular values of $A$ #text(weight: "bold")[squared], and the eigenvectors are the columns of $V$.

  For $A A^*$, we have:

  $
    A A^* = (U Sigma V^*) (U Sigma V^*)^* = U Sigma V^* V Sigma^* U^* = U Sigma Sigma^* U^* = U Sigma^2 U^*
  $

  The reasoning here is analogous. So the proof is complete.
]

By @theorem_eigencalculation_of_svd, calculating the SVD of $A$ has been reduced to calculating the eigenvalues and eigenvectors of $A^* A$ and $A A^*$, here are some examples of singular value decompositions:

#example[
  Consider $A = mat(3, 2; 2, 3)$. Computing the SVD:
  
  First, find $A^* A = mat(13, 12; 12, 13)$ and calculate its eigenvalues:
  $lambda_1 = 25, lambda_2 = 1$
  
  The singular values are $sigma_1 = 5, sigma_2 = 1$.
  
  The right singular vectors (eigenvectors of $A^* A$):
  $V = mat(1/sqrt(2), 1/sqrt(2); 1/sqrt(2), -1/sqrt(2))$
  
  The left singular vectors (obtained from $A v_i = sigma_i u_i$):
  $U = mat(1/sqrt(2), 1/sqrt(2); 1/sqrt(2), -1/sqrt(2))$
  
  Therefore, the SVD is:
  $A = mat(1/sqrt(2), 1/sqrt(2); 1/sqrt(2), -1/sqrt(2)) dot mat(5, 0; 0, 1) dot mat(1/sqrt(2), 1/sqrt(2); 1/sqrt(2), -1/sqrt(2))^*$
]

#example[
  Consider a non-square matrix $A = mat(1, 1, 0; 0, 1, 1)$. For this $2 times 3$ matrix, for the SVD we do:
  
  
  $
    A^* A = mat(1, 1, 0; 1, 2, 1; 0, 1, 1)
  $
  
  The eigenvalues of $A^* A$ are $lambda_1 = 3, lambda_2 = 1, lambda_3 = 0$, so the singular values are $sigma_1 = sqrt(3), sigma_2 = 1, sigma_3 = 0$
  
  The right singular vectors (eigenvectors of $A^* A$) are:

  $
    V = mat(1/2, -1/sqrt(2), 1/2; 1/sqrt(2), 0, -1/sqrt(2); 1/2, 1/sqrt(2), 1/2)
  $
  
  And now for $A A^*$:

  $
    A A^* = mat(2, 1; 1, 2)
  $ 
  
  The eigenvalues are $lambda_1 = 3, lambda_2 = 1$, so the singular values are $sigma_1 = sqrt(3), sigma_2 = 1$. The eigenvectors are:

  $
    U = mat(1/sqrt(2), 1/sqrt(2); 1/sqrt(2), -1/sqrt(2))
  $
  
  Therefore, the full SVD is:

  $
    A = mat(1/sqrt(2), 1/sqrt(2); 1/sqrt(2), -1/sqrt(2)) dot mat(sqrt(3), 0, 0; 0, 1, 0) dot mat(1/2, -1/sqrt(2), 1/2; 1/sqrt(2), 0, -1/sqrt(2); 1/2, 1/sqrt(2), 1/2)^*
  $
]

== Least Squares with QR and SVD
<section_python_qr_svd>

Here we will write code that solves the least squares problem usig the 2 factorizations shown in @section_QR_decomposition and @section_SVD_decompositon, as well as the ordinary approach to least squares shown in @section_simple_linear_regression.

The following code has functions that solve the least squares problem using the QR and SVD decompositions, as well as through the normal approach shown in @section_simple_linear_regression. A function to build the matrix $A$ is also included. With a boolean parameter to alternate betweeen a centered ($t_i = i/m$) or a non-centered $(t_i = i/m - 1/2)$ dataset:

```python
from numpy.linalg import qr, svd, solve, pinv, cond #easier

def ls_qr(A, b):

    """
    Solves the least squares problem Ax = b using QR decomposition.

    Args:
        A : system matrix
        b : right-hand side vector
    Returns:
        x : solution vector
        y : fitted values (Ax)
    """

    Q, R = qr(A, mode='reduced')
    x = solve(R, Q.T @ b)
    y = A @ x
    
    return x, y

def ls_svd(A, b):

    """
    Solves the least squares problem Ax = b using SVD decomposition.

    Args:
        A : system matrix
        b : right-hand side vector
    Returns:
        x : solution vector
        y : fitted values (Ax)
    """

    U, S, Vt = svd(A, full_matrices=False)
    
    S_inv = np.zeros_like(S) #calculate the pseudo-inverse using SVD, x = V * S⁻¹ * U.T * b

    tol = np.finfo(float).eps * max(A.shape) * S[0] #tolarance for singular values
    
    for i in range(len(S)):
        if S[i] > tol:
            S_inv[i] = 1.0 / S[i] #inverts if above tolarance
    
    x = (Vt.T @ np.diag(S_inv) @ U.T) @ b
    y = A @ x
    
    return x, y

def ls_normal(A, b):

    """
    Solves the least squares problem Ax = b using normal equations.

    Args:
        A : system matrix
        b : right-hand side vector
    Returns:
        x : solution vector
        y : fitted values (Ax)
    """

    ATA = A.T @ A
    ATb = A.T @ b
    
    x = solve(ATA, ATb)
    y = A @ x
    
    return x, y

def build_A_matrix(m, n, centralized=False):

    """
    Creates the matrix A for polynomial regression of degree n. 

    Args:
        m (int): number of subintervals (m >= 0)
        n (int): polynomial degree (n >= 0)
        centralized (bool): if True, use centralized points
    Returns:
        A (np.ndarray): shape (m+1, n+1) matrix for polynomial regression
        t (np.ndarray): array of points used to create the matrix
    """

    if centralized:
        t = np.array([i/m - 1/2 for i in range(m+1)])
    else:
        t = np.array([i/m for i in range(m+1)])
    
    A = np.zeros((m+1, n+1))
    
    for i in range(n+1):
        A[:, i] = t**i  #clever
        
    return A, t
``` 

== Examples (2b)

We will use the functions defined on the previous sections to do regression on the following functions:

$
  f(t) = sin(t)\
  g(t) = e^t\
  h(t) = cos(3t)
$ <functions_to_be_numerically_analysed>

The results are shown:

#figure(
  image("plots/linear_regression_cos.png", width: 100%),
  caption: [
    Linear regression on $cos(t)$
  ]
) <linear_regression_cos>

#figure(
  image("plots/linear_regression_et.png", width: 100%),
  caption: [
    Linear regression on $e^t$
  ]
) <linear_regression_et>

#figure(
  image("plots/linear_regression_sin.png", width: 100%),
  caption: [
    Linear regression on $sin(t)$
  ]
) <linear_regression_sin>

Plots @linear_regression_cos, @linear_regression_et, @linear_regression_sin shows us that all 3 methods result in the same line. This is expected, for we are using $100$ equally sparse points on a small portion of the image of the functions.




== How good are the approximations? (2c)

In @section_perturbations_polynomial_regression, we have seen that polynomial regression is not a very well-conditioned problem, so here we will analyze the errors produced by such method. We use what has been done in the previous section.

The following code plots the errors:

```python
m = 100
t = np.linspace(0, 1, m+1)

function_data = [(res['name'], res['values']) for res in results] #data here

max_degree = 15
all_errors = {name: {'qr': [], 'svd': [], 'normal': []} for name, _ in function_data}
all_cond_numbers = []

for n in range(1, max_degree + 1):
    print(f"\nPolynomial degree: {n}")
    
    A, _ = build_A_matrix(m, n) #matrix for regression
    cond_num = cond(A)
    all_cond_numbers.append(cond_num)
    print(f"Condition number of A: {format_scientific(cond_num)}")
    
    for name, values in function_data: #forall funcs, solve using all methods
        try:
            x_qr, y_qr = ls_qr(A, values)
            error_qr = np.linalg.norm(y_qr - values)
            all_errors[name]['qr'].append(error_qr)
        except Exception as e:
            print(f"Error with QR for {name}, degree {n}: {e}")
            all_errors[name]['qr'].append(np.nan)
            
        try:
            x_svd, y_svd = ls_svd(A, values)
            error_svd = np.linalg.norm(y_svd - values)
            all_errors[name]['svd'].append(error_svd)
        except Exception as e:
            print(f"Error with SVD for {name}, degree {n}: {e}")
            all_errors[name]['svd'].append(np.nan)
            
        try:
            x_normal, y_normal = ls_normal(A, values)
            error_normal = np.linalg.norm(y_normal - values)
            all_errors[name]['normal'].append(error_normal)
        except Exception as e:
            print(f"Error with Normal Equations for {name}, degree {n}: {e}")
            all_errors[name]['normal'].append(np.nan)
        
        #show errors
        print(f"{name} - Errors: QR: {format_scientific(all_errors[name]['qr'][-1])}, "
              f"SVD: {format_scientific(all_errors[name]['svd'][-1])}, "
              f"Normal: {format_scientific(all_errors[name]['normal'][-1])}")

degrees = list(range(1, max_degree + 1))

for name, _ in function_data:
    plt.figure(figsize=(6, 5), dpi=150) 
    plt.semilogy(degrees, all_errors[name]['qr'],     'ro-', label='QR')
    plt.semilogy(degrees, all_errors[name]['svd'],    'gs-', label='SVD')
    plt.semilogy(degrees, all_errors[name]['normal'], 'bx-', label='Normal Eqs.')

    plt.title(f'Error vs. Polynomial Degree for {name}')
    plt.xlabel('Polynomial degree $n$')
    plt.ylabel('Residual error (log scale)')
    plt.grid(True, which='both', ls=':')
    plt.legend()
    plt.tight_layout()
    plt.show()

plt.figure(figsize=(6, 5), dpi=150)
plt.semilogy(degrees, all_cond_numbers, 'mo-')
plt.xlabel('Polynomial degree $n$')
plt.ylabel('Condition number (log scale)')
plt.title('Condition Number vs. Polynomial Degree')
plt.grid(True, which='both', ls=':')
plt.tight_layout()
plt.show()
```

The expected output are the plots:

#figure(
  image("plots/error_analysis_cost.png", width: 100%),
  caption: [
    Error on $cos(3t)$
  ]
) <figure_errors_cos3t>

#figure(
  image("plots/error_analysis_et.png", width: 100%),
  caption: [
    Error on $e^t$
  ]
) <figure_error_et>

#figure(
  image("plots/error_analysis_sint.png", width: 100%),
  caption: [
    Error on $sin(t)$
  ]
) <figure_error_sint>

We know that the functions $e^t$ and $sin(t)$ are crescent on $[0,1]$ and shifted up, so the linear model would require higher degree to better describe this curve.

- #text(weight: "bold")[Normal System:] The error falls considerably until $approx 10^(-8)$, which is expected. Since the condition numbr is squared, it climbs sooner than the 2 other approaches. 

- #text(weight: "bold")[SVD:] After $approx n = 8$ it starts climbing, due to the loss of significance from the ill-conditioned Vandermonde Matrix (huge condition number due to a big $n$). We see that the normal system of equations climbs sooner than the SVD, due to the difference in magnitude of the condition numbers

- #text(weight: "bold")[QR:] Climbs all the way down until machine precision, which is expected from the QR factorization.

Now the function $cos(3t)$ behaves differently. Expanding it on the Taylor Series:

$
  cos(3t) = 1 - (3t)^2 / 2! + (3t)^4 / 4! - (3t)^6 / 6! + dots
$

Only _even_ powers appear. So when the degree increases from an odd $phi$ to the next odd $phi + 1$, no even power is added. The expansion gains a term that has zero coefficient, as a consequence, the error does not improve. This explains the nearly horizontal parts between consecutive degrees.
