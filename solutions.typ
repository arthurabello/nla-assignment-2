#import "@preview/ctheorems:1.1.3": *
#import "@preview/plotst:0.2.0": *
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

#let example = thmplain("example", "Example").with(numbering: none)
#let proof = thmproof("proof", "Proof")

#align(right, text(13pt)[ //L
  FGV - EMAP
])

#align(center, text(22pt)[
  Assignment 2 - Numerical Linear Algebra
])
#align(center, text(18pt)[
  Prof.: Bernardo Freitas Paulo da Costa
])
#align(center, text(16pt)[
  TA: Beatriz Lúcia Teixeira de Souza
])
#align(center, text(16pt)[
  Student: Arthur Rabello Oliveira

  #datetime.today().display("[day]/[month]/[year]")
])

#pagebreak()

#outline()

#pagebreak()

= Problem 1: Traditional Least Squares ( a - f )

== Linear Regression (a)

We have a set of equally spaced points $S := {t_i = i/m}, i = 0, 1, dots, m$, we will find the best _line_ $f(t) = alpha + beta t$ that approximates the points $(t_i, b_i) in RR^2$

The system of equations to be solved is to be given as a function of $t_i, b_i, m$.

#text(weight: "bold")[Solution:]

Approximating 2 points in $RR^2$ by a line is trivial, now approximating more than 2 points is a task that requires linear algebra. To see this, we will analyze the following example to build intuition for the general case:

#figure(
  image("least_squares_idea.png", width: 80%),
  caption: [
    A glimpse into what we want to see
  ],
)

Given the points $(1, 1), (2, 2(, (3, 2) in RR^2$, we would like a _line_ $y = alpha + beta t$ that best approximates these 3 points, in other words, since we know that the line does not pass through all of the 3 points, we would like to find the _closest_ line to the line that would pass throught the 3 points, so the system:

$
  f(1) = alpha + beta = 1\
  f(2) = alpha + 2 beta = 2\
  f(3) = alpha + 3 beta = 2
$

Clearly has no solution, (the line does not cross the 3 points), but it has a _closest solution_, which we can find throught #text(weight: "bold")[projections], the system is:

$
  underbrace(mat(
    1, 1;
    1, 2;
    1, 3
  ), A) dot underbrace(mat(alpha;beta), x) = underbrace(mat(1;2;2), b)
$

Let $x^* != x$ be a solution to the system, we want to #text(weight: "bold")[minimize] the error produced by approximating the points through a line, so if the #text(weight: "bold")[error] is $e = A x - b$, we want the smaller error _square_ possible (that is why least squares). We square the error to avoid and detect outliers, so:

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

In this case, we will project this system into the column space of the matrix $A$, giving us the closest solution, and the least squares solutions is when $hat(x)$ minimizes $||A x - b||^2$, this occurs when the residual $e = A x - b$ is orthogonal to $C(A)$, since $N(A^T) perp C(A)$ and the dimensions sum up the left dimension of the matrix, so:

$
  A^T A hat(x) = A ^T b\

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
  e_1^2 = (alpha + beta - 1)^2\
  e_2^2 = (alpha + 2 beta - 2)^2\
  e_3^2 = (alpha + 3 beta - 2)^2
$

The system in @system_1 is _precisely_ what is obtained after using partial derivatives to minimize the erros sum as a function of $(alpha, beta)$:

$
  f(alpha, beta) = (alpha + beta - 1)^2 + (alpha + 2 beta - 2)^2 + (alpha + 3 beta - 2)^2\
  = 3 alpha^2 + 14 beta^2 + 12 alpha beta - 10 alpha- 22 beta + 9,\

  (diff f) / (diff alpha) = (diff f) / (diff beta) = 0 <=> 6 alpha + 12 d - 10 = 28 beta + 12 alpha - 22 = 0 <=> cases(
    3c + 6d = 11, 6c + 14d = 11
  )
$

This new system has a solution in $hat(alpha) = 2/3 , hat(beta) = 1/2$, so the equation of the optimal line, obtained through _linear regression_ (or least squares) is:

$
  y(t) = 2/3 + 1/2 t.
$

If we have $n > 3$ points to approximate throught a line, the reasoning is analogous:






