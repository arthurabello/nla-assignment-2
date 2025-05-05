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
  We derive linear and polynomial regression in subsets of $RR$ and discuss the condition number of the associated matrices, numerical algorithms for the SVD and QR factorization are built and used on an efficiency analysis of the 3 methods to do linear or polynomial regression, stability of these algorirths is  mentioned and 
]

#outline()

#pagebreak()

= Introduction

Given $D subset RR^2$, a dataset, approximating this set through a _continuous_ $f: RR -> RR$ is a very important problem in statistics, we will derive the 2 most important and most used methods to do this: linear and polynomial regression. Both are based on the least squares minimization problem. We will also discuss the conditioning number of the problems shown. A computational approach to regression is shown as well. We discuss how the condition number changes when the matrix is QR or SVD decomposed, and the algorithms for such decompositions are built.

= Condition of a Problem

A _problem_ is usually described as a function $f: X -> Y$ from a #text(weight: "bold")[normed] vector space $X$ of data (it has to be normed so qe can _quantify_ data) and a _normed_ vector space $Y$ of solutions, $f$ is not always a well-behaved continuous function, which is why we are interested in #text(weight: "bold")[well-conditioned] problems and not in  #text(weight: "bold")[ill-conditioned] problems, which we define:

#definition[(Well-Conditioned Problem)
  A problem $f:X -> Y$ is _well-conditioned_ at $x_0 in X$ $<=> forall epsilon > 0, exists delta > 0 |  ||x - x_0|| < delta => ||f(x) - f(x_0)|| < epsilon$.
] <definition_of_stable_problem>

This means that small perturbations in $x$ lead to small changes in $f(x)$, a problem is #text(weight: "bold")[ill-conditioned] if $f(x)$ can suffer huge changes with small changes in $x$.

We usually say $f$ is well-conditioned if it is well-conditioned $forall x in X$, if there is at least one $x_i$ in which the problem is ill-conditioned, then we can use that whole problem is ill-conditioned.


== The Condition number of a problem

Conditioning numbers are a tool to quantify how well/ill conditioned a problem is:

#definition[(Absolute Conditioning Number)
  Let $delta x$ be a small pertubation of $x$, so $delta f = f(x + delta x) - f(x)$. The #text(weight: "bold")[absolute] conditioning number of $f$ is:

  $
    hat(kappa) = lim_(delta -> 0) sup_(||delta x|| <= delta) (||delta f||) / (||delta x||)
  $
] <definition_absolute_conditioning_number>

The limit of the supremum can be seen as the supremum of all _infinitesimal_ perturbations, so this can be rewritten as:

$
  hat(kappa) = sup_(delta x) (||delta f||) / (||delta x||)
$

If $f$ is differentiable, we can evaluate the abs.conditioning number using its derivative, if $J$ is the matrix whose $i times j$ entry is the derivative $(diff f_i) / (diff x_j)$ (jacobian of $f$), then we know that $delta f approx J(x) delta x$, with equality in the limit $||delta x|| -> 0$. So the absolute conditioning number of $f$ becomes:

$
  hat(kappa) = ||J(x)||,
$ <absolute_conditioning_number__jacobian>

== The Relative Conditioning Number

When, instead of analyzing the whole set $X$ of data, we are interested in _relative_ changes, we use the #text(weight: "bold")[relative condition number:]

#definition[(Relative Condition Number)
  Given $f:X -> Y$ a problem, the _relative condition number_ $kappa(x)$ at $x in X$ is:

  $
    kappa(x) = lim_(delta -> 0) sup_(||delta x|| <= delta) ((||delta f||) / (||f(x)||)) dot ((||delta x||) / (||x||))^(-1)
  $
] <definition_relative_condition_number>

Or, as we did in @definition_absolute_conditioning_number, assuming that $delta f$ and $delta x$ are infinitesimal:

$
  kappa(x) = sup_(delta x) ((||delta f||) / (||f(x)||)) dot ((||delta x||) / (||x||))^(-1)
$

If $f$ is differentiable:

$
  kappa(x) = (||J(x)||) dot ((||f(x)||) / (||x||))^(-1)
$ <relative_condition_number_through_jacobian>

Relative condition numbers are more useful than absolute conditioning numbers because the #text(weight: "bold")[floating point arithmetic] used in many computers produces _relative_ errors, the latter is not a highlight of this discussion.

Here are some examples of conditioning:

#example[
  Consider the problem of obtaining the scalar $x/2$ from $x in RR$. The function $f(x) = x/2$ is differentiablle, so by @relative_condition_number_through_jacobian:

  $
    kappa(x) = (||J||) dot ((||f(x)||) / (||x||))^(-1) = (1/2) dot ((x/2) / x)^(-1) = 1.
  $
]

This problem is well-conditioned ($kappa$ is small).

#example[
  Consider the problem of computing the scalar $x_1 - x_2$ from $(x_1, x_2) in RR^2$ (Use the $infinity$-norm in $RR^2$ for simplicity). The function associated is differentiable and the jacobian is:

  $
    J = mat((diff f) / (diff x_1) , (diff f) / (diff x_2)) = mat(1,-1)
  $
  With $||J||_infinity$ = 2, so the condition number is:

  $
    kappa = (||J||_infinity) dot ((||f(x)||) / (||x||))^(-1) = 2 / (|x_1 - x_1| dot max{|x_1|, |x_2|})
  $

  This problem can be ill-conditioned if $|x_1 - x_2| approx 0$ ($kappa$ gets huge), and well-conditioned otherwise
].

== Condition Number of Matrices

We will deduce the conditioning number of a matrix from the conditioning number of _matrix-vector_ multiplication:

Consider the problem of obtaining $A x$ given $A in CC^(m times n)$. We will calculate the relative condition number with respect to perturbations on $x$. Directly from @definition_relative_condition_number, we have:

$
  kappa = sup_(delta x) (||A (x + delta x) - A x||) / (||A x||) dot ((||delta x||) / (||x||))^(-1) = sup_(delta x) (||A delta x||) / (||delta x||) dot ((||A x||) / (||x||))^(-1)
$ <beginning_proof>

Since $sup_(forall x) (||A delta x||) / (||delta x||) = ||A|| $, we have:

$
  kappa = ||A|| dot (||x||) / (||A x||)
$ <condition_number_of_matrix_vector_multiplication>

This is a precise formula as a function of $(A, x)$.

The following theorem will be useful in a near future:

#theorem[
  $forall x in CC^n, A in CC^(n times n), det(A) != 0$, the following holds:

  $
      (||x||) / (||A x||) <= ||A^(-1)||
  $
] <theorem_inequality_norms>

#proof[
  Since $forall A, B in CC^(n times n), ||A B || <= ||A|| ||B||$, we have:

  $
    ||A A^(-1) x || <= ||A x|| ||A^(-1)|| <=> (||x||) / (||A x||) <= ||A^(-1)||
  $
]

So using this in @condition_number_of_matrix_vector_multiplication , we can write:

$
  kappa <= ||A|| dot ||A^(-1)||
$

Or:

$
  kappa = alpha ||A|| dot ||A^(-1)||
$

With

$
  alpha = (||x||) / (||A x||) dot (||A^(-1)||)^(-1)
$ <end_proof>

From @theorem_inequality_norms, we can choose $x$ to make $alpha = 1$, and therefore $kappa = ||A|| dot ||A^(-1)||$.

Consider now the problem of calculating $A^(-1) b$ given $A in CC^(n times n)$. This is mathematically identical to the problem we just analyzed, so the following theorem has already been proven:

#theorem[
  Let $A in CC^(n times n), det(A) != 0$, and consider the problem of computing $b$, from $A x = b$, by perturbating $x$. Then the following holds:

  $
    kappa = ||A|| (||x||) / (||b||) <= ||A|| dot ||A^(-1)||
  $
  
  Where $kappa$ is the condition number of the problem.
]

#proof[
  Read from @beginning_proof to @end_proof.
]

Finally, $||A|| dot ||A^(-1)||$ is so useful it has a name: #text(weight: "bold")[the condition number of A] (relative to the norm $norm(dot)$)

If $A$ is singular, usually we write $kappa(A) = infinity$. Notice that if $||dot|| = ||dot||_2$, then $||A|| = sigma_1$ and $||A^(-1)|| = 1/sigma_m$, so:

$
  kappa(A) = sigma_1 / sigma_m
$

This is the condition number of $A$ with respect to the $2$-norm, which is the most used norm in practice. The condition number of a matrix is a measure of how sensitive the solution of a system of equations is to perturbations in the data. A large condition number indicates that the matrix is ill-conditioned, meaning that small changes in the input can lead to large changes in the output.

= Linear Regression (1a)
<section_simple_linear_regression>

Given a dataset of equally spaced points $D := {t_i = i/m}, i = 0, 1, dots, m in RR$, linear regression consists of finding the best _line_ $f(t) = alpha + beta t$ that approximates the points $(t_i, b_i) in RR^2$, where $b_i$ are arbitrary

Approximating $2$ points in $RR^2$ by a line is trivial, now approximating more points is a task that requires linear algebra. To see this, we will analyze the following example to build intuition for the general case:

#figure(
  image("least_squares_idea.png", width: 80%),
  caption: [
    A glimpse into what we want to see
  ],
)

Given the points $(1, 1), (2, 2), (3, 2) in RR^2$, we have $(t_1, b_1) = (1, 1), (t_2, b_2) = (2, 2), (t_3, b_3) = (3, 2) $ we would like a _line_ $f(t) = y(t) = alpha + beta t$ that best approximates $(t_i, b_i)$, in other words, since we know that the line does not pass through all 3 points, we would like to find the _closest_ line to #text(weight: "bold")[each point] of the dataset $D$, so the system:

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

Clearly has no solution, (the line does not cross the 3 points), but it has a _closest solution_, which we can find through #text(weight: "bold")[minimizing] the errors produced by this approximation.

Let $x^* != x$ be a solution to the system, let the error produced by approximating the points through a line be $e = A x - b$, we want the smaller error _square_ possible (that is why least squares). We square the error to avoid and detect outliers, so:

$
  e_1^2 + e_2^2 + e_3^2
$

Is what we want to minimize, where $e_i$ is the error (distance) from the ith point to the line:

#figure(
  image("least_squares_with_errors.png", width: 80%),
  caption: [
    The errors (distances)
  ],
)

So we will project $b$ into $C(A)$, giving us the closest solution, and the least squares solutions is when $hat(x)$ minimizes $||A x - b||^2$, this occurs when the residual $e = A x - b$ is orthogonal to $C(A)$, since $N(A^*) perp C(A)$ and the dimensions sum up the left dimension of the matrix, so by the well-known projection formula, we have:

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

So the system to find $hat(x) = mat(hat(alpha);hat(beta))$ becomes:

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

The system in @system_1 is _precisely_ what is obtained after using partial derivatives to minimize the erros sum as a function of $(alpha, beta)$:

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

If we have $n > 3$ points to approximate through a line, the reasoning is analogous:

Going back to $D$, we want to find the extended system as we did in @system_2, so let the best line be:
$
  f(t) = alpha + beta t
$

That best approximates the points $(0, b_0), (1/m, b_1), dots , (1, b_m)$. The system is:

$
  f(0) = b_0 = alpha,\
  f(1/m) = b_1 = alpha + beta/m,\
  f(2/m) = b_2 = alpha + 2/m beta\
  dots\
  f(1) = b_m = alpha + beta
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

And the least squares optimal solution is:

$
  mat(hat(alpha); hat(beta)) = mat(
    m + 1, sum_(i = 1)^m t_i;
    sum_(i = 1)^m t_i, sum_(i = 1)^m t_i^2;
  )^(-1) dot mat(
    sum_(i = 1)^m b_i;
    sum_(i = 1)^m i/m dot b_i
  )
$

= How the conditioning number of A changes (1b)
<condition_number_linear_regression_system>

We are interested in the condition number of linear regression, which is the condition number of the matrix $A$ in @matrix_system_least_squares. We will analyze how the condition number of $A$ changes with respect to perturbations $m$, the number of points in the dataset. A computational approach is appropriate.

Here is a python code that numerically calculates many values of $kappa(A) = f(m)$ as a function of $m$:

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

def main(): 
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
    plt.ylabel('Condition number κ₂(A)')
    plt.title('Condition number of A(m) over [0, M]')
    plt.grid(True)
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
```

Good plots are:

#figure(
  image("condition_number.png", width: 80%),
  caption: [
    Condition number of A(m) over [0, 100]
  ],
) <condition_number_plot>

#figure(
  image("condition_number_2.png", width: 80%),
  caption: [
    Condition number of A(m) over [0, 10000]
  ],
) <condition_number_plot_2>


@condition_number_plot and @condition_number_plot_2 show us that it looks like $f(m) = kappa(A_m)$ converges to a real number, we will evaluate this hypothesis below:

Using $||dot||_2$ the conditioning number of $hat(A) = A^* A$ in @matrix_system_least_squares is:

$
  kappa(hat(A)) = norm(hat(A))_2 dot norm(hat(A)^(-1))_2 = sigma_1 / sigma_m 
$

Singular Values are better explored in @section_SVD_decompositon. Now we will calculate the singular values of $hat(A)$, which are the square roots of the eigenvalues of $hat(A)$ (see @theorem_eigencalculation_of_svd). So we have:
$
  det(hat(A) - lambda I) = 0 <=> det(mat(
    m+1 - lambda , (m + 1) / 2;
    (m + 1) / 2 , ((m + 1) (2m + 1)) / 6 - lambda
  )) = 0\ 

  <=> (m + 1 - lambda) [((m + 1) (2m + 1)) / 6 - lambda] - ((m + 1)/2)^2 = 0\

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
  kappa(A) = sigma_1 / sigma_m = sqrt((m+1)/(12m) [(8m + 1) + sqrt(52m^2-  8m + 1)]) / sqrt((m+1)/(12m) [(8m + 1) - sqrt(52m^2-  8m + 1)])\

  = sqrt(((8m + 1) + sqrt(52m^2-  8m + 1)) / ((8m + 1) - sqrt(52m^2-  8m + 1)))
$ <condition_number_problem_b>

And the limit as $m$ grows is:

$
  lim_(m -> oo) kappa(A) = lim_(m -> oo) sqrt(((8m + 1) + sqrt(52m^2-  8m + 1)) / ((8m + 1) - sqrt(52m^2-  8m + 1)))\
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
  lim_(m -> oo) kappa(A) = (8 + sqrt(52)) / sqrt(12) = (4 + sqrt(13)) / sqrt(3)  
$

A very good visualization of this is:

#figure(
  image("desmos_condition_number.png", width: 60%),
  caption: [
    The purple line is the limit and the red is the function @condition_number_problem_b
  ],
) <condition_number_desmos_plot>


@condition_number_desmos_plot shows the function approaching the limit. One could say that this problem is well conditioned, for $kappa(A_m) < (4 + sqrt(13)) / sqrt(3), forall m > 0$, and $(4 + sqrt(13)) / sqrt(3)$ is not a very big number. We will not go deep into the discussion of how well-condition this problem is, but we can say that the condition number of $A$ is not a problem for the linear regression algorithm.

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

And the least squares _polynomial regression_ solution is:

$
  mat(hat(phi_0); hat(phi_1); hat(phi_2); dots.v; hat(phi_n)) = mat(
    m + 1, sum_(i = 1)^m i/m, dots, sum_(i = 1)^m (i/m)^n;
    sum_(i = 1)^m i/m, sum_(i = 1)^m (i/m)^2, dots, sum_(i = 1)^m (i/m)^(n+1);
    sum_(i = 1)^m (i/m)^2, sum_(i = 1)^m (i/m)^3, dots, sum_(i = 1)^m (i/m)^(n+2);
    dots.v, dots.v, dots.v, dots.v;
    sum_(i = 1)^m (i/m)^n, sum_(i = 1)^m (i/m)^(n+1), dots, sum_(i = 1)^m (i/m)^(2n)
 )^(-1) dot mat(
    sum_(i = 0)^m b_i; sum_(i = 0)^m (i b_i) / m; sum_(i = 0)^m (i/m)^2 m; dots.v; sum_(i = 0)^m (i/m)^n b_i
  )
$
= Computing the matrix A given (m,n) (1d)
<poly_ls_section>

Here is a python function that calculates the polynomial regression matrix from @matrix_system_polynomial_least_squares, given the dimensions $(m, n)$:

```python
import numpy as np

def poly_ls(m, n):

    """
    Build the (n+1) x (n+1) matrix A for least-squares polynomial fitting.
    
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
    for j in range(n+1):
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
  A(1, 1) = mat(
    2, 1;
    1, 1
  )\

  A(2, 2) = mat(
    3, 1.5, 1.25;
    1.5,   1.25,  1.125;
    1.25,  1.125, 1.062
  )\

  A(2, 3) = mat(
    3, 1.5, 1.25, 1.125;
    1.5,   1.25,  1.125, 1.062;
    1.25,  1.125, 1.062, 1.031;
    1.125, 1.062, 1.031, 1.016
  )
$

= How Perturbations Affect The Condition Number of A (1e)

Still on polynomial regression, in this section we analyze what happens to $kappa(A)$, when $A$ is perturbated with $m = 100$ and $n = 1, dots, 20$.

We will run _poly_ls(m, n)_ built in @poly_ls_section for $m = 100$ and $n = 1, dots, 20$ and then numerically calculate the condition number of the matrices, the following code is used:

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
  image("(m,n)condition_perturbation.png", width: 80%),
  caption: [
    Growth of the condition number of A(m) for polynomial regression
  ],
) <condition_number_perturbation_plot>

@condition_number_perturbation_plot Shows that _magic_ happens





= Condition Analysis of a Different Dataset
== A Different Dataset

If we change $S := {(t_i, b_i) | t_i = i / m, i = 0, 1, dots, m}$ to $hat(S) = {(t_i, b_i) | t_i = i/m - 1/2}$, the polynomial regression becomes:

$
  p(t_0 = -1/2) = phi_0 + phi_1 (-1/2) + dots + phi_n (-1/2)^n = b_0\
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
  ), hat(Phi))\

  = mat(
    n + 1, sum_(i = 0)^n (i/m - 1/2), sum_(i = 0)^n (i/m - 1/2)^2, dots,  sum_(i = 0)^n (i/m - 1/2)^n;
    sum_(i = 0)^n i/m - 1/2,sum_(i = 0)^n (i/m - 1/2)^2, sum_(i = 0)^n (i/m - 1/2)^3, dots, sum_(i = 0)^n (i/m - 1/2)^(n+1);
    dots.v, dots.v, dots.v, dots.v, dots.v;
    sum_(i = 0)^n (i/m - 1/2)^n, sum_(i = 0)^n (i/m - 1/2)^(n+1), sum_(i = 0)^n (i/m - 1/2)^(n+2), dots, sum_(i = 0)^n (i/m - 1/2)^(2n)
  ) dot mat(
    hat(phi_0);
    hat(phi_1);
    dots.v;
    hat(phi_n)
  )\

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

$therefore$ The least squares optimal solution $(hat(phi_0), hat(phi_1) , dots, hat(phi_n))$ is:

$
  mat(
    n + 1, sum_(i = 0)^n (i/m - 1/2), dots,  sum_(i = 0)^n (i/m - 1/2)^n;
    sum_(i = 0)^n i/m - 1/2, sum_(i = 0)^n (i/m - 1/2)^2,dots, sum_(i = 0)^n (i/m - 1/2)^(n+1);
    dots.v, dots.v, dots.v, dots.v;
    sum_(i = 0)^n (i/m - 1/2)^n, sum_(i = 0)^n (i/m - 1/2)^(n+1), dots, sum_(i = 0)^n (i/m - 1/2)^(2n)
  )^(-1) dot mat(
    sum_(i = 0)^n b_i;
    sum_(i = 0)^n (i/m - 1/2) b_i;
    sum_(i = 0)^n (i/m - 1/2)^2 b_i;
    dots.v;
    sum_(i = 0)^n (i/m - 1/2)^n b_i
  )
$ <new_dataset_polynomial_regression>

We provide numerical examples in the next section for a better visualization of @new_dataset_polynomial_regression.

== How Conditioning changes (1f)

BOTA AS PORRA DOS GRAFICO


= Least Squares with QR adn SVD decompositions

We have shown the solutions to the least squares problem $A x = b$, but this problem could be solved with factorizations of $A$, such as the QR and SVD, in the following sections we will define these factorizations and use them to solve the least squares problem.


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
  image("full_qr.png", width: 60%),
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
      1, 2, 3;
      4, 5, 6;
      7, 8, 9
    ) = 
      Q = mat(
        1/sqrt(66), -7/sqrt(246), -2/15;
        4/sqrt(66), 1/sqrt(246), -11/15;
        7/sqrt(66), -5/sqrt(246), 8/15
      ),

      R = mat(
        sqrt(66), 5 sqrt(66)/6, 4 sqrt(66)/3;
        0, sqrt(246)/6, 5 sqrt(246)/6;
        0, 0, 0
      )
    $

    This can be verified by computing $Q R$ and checking that it equals $A$. You can also verify that $Q^* Q = I$, which shows that $Q$ has orthonormal columns.
  ]

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
  image("hyperellipse.png", width: 80%),
  caption: [
    SVD of a $2 times 2$ matrix
  ],
)

So the independent directions $v_1, v_2$ have been mapped to another set of orthogonal directions $sigma_1 v_1, sigma_2 v_2$, so with $S:= {v in CC^n | ||v|| = 1}$ as the unit ball,  let's define:

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
  We proceed by fixing the largest image of $A$ and using induction on the dimension of $A$:

  Let $sigma_1 = ||A||_2$. There must exist unitary vectors $u_1, v_1 in CC^n$ such that $A v_1 = sigma_1 u_1$. PROVA ESSA POPRRA DIREITO
]

Is the SVD factorization of A. There are more about the SVD on computing $U, Sigma, V^*$, as we will show below:

#theorem[
  $forall A in CC^(m times n)$, the following holds:
  
  - The eigenvalues of $A^*A$ are the singular values _squared_ of $A$, and the column-eigenvectors of $A^*A$ form the matrix $V$.

  - The eigenvalues of $A A^*$ are the singular values _squared_ of $A$, and the column-eigenvectors of $A A^*$ form the matrix $U$.
]<theorem_eigencalculation_of_svd>

#proof[
  PROVA ESSA PORRA DIREITO
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

== A Python-Implementation of Least Squares with Decompositions (2a)
<section_python_qr_svd>

Here we will write code that solves the least squares problem usig the 2 factorizations shown in @section_QR_decomposition and @section_SVD_decompositon, as well as the ordinary approach to least squares shown in @section_simple_linear_regression.

== Examples (2b)

We will also use these algorithms to do linear regression on the simple functions $f, g, h: RR -> RR$ defined as:

$
  f(t) = sin(t)\
  g(t) = e^t\
  h(t) = cos(3t)
$ <functions_to_be_numerically_analysed>

BOTA OS CODIGO

== The polynomial approach: An efficiency analysis (2c)

Here we will analyse what happens when we do _polynomial_ regression with the tools shown in @section_python_qr_svd. The same functions of @functions_to_be_numerically_analysed will be used here, with polynomials of degree up to $n = 15$:

BOTA AS PORRA DOS CODIGO


