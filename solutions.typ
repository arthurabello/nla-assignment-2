#import "@preview/ctheorems:1.1.3": *
#import "@preview/plotst:0.2.0": *
#show: thmrules.with(qed-symbol: $square$)

#set heading(numbering: "1.1.")

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

#set page(numbering: "1")
#set heading(numbering: "1.")

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

Given the points $(1, 1), (2, 2(, (3, 2) in RR^2$, we would like a _line_ $y = c + d t$ that best approximates these 3 points, in other words, since we know that the line does not pass through all of the 3 points, we would like to find the _closest_ line to the line that would pass throught the 3 points, so the system:

$
  f(1) = c + d = 1\
  f(2) = c + 2d = 2\
  f(3) = c + 3d = 2
$

Clearly has no solution, but

#figure(
  image("Least_Squares_Idea.png", width: 80%),
  caption: [
    A glimpse into what we want to see
  ],
)

