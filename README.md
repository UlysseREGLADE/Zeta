# A geometrical summation method for the Riemann zêta function

Codes used for the ArXiv preprint:

[https://arxiv.org/pdf/1903.10853.pdf](https://arxiv.org/pdf/1903.10853.pdf)

These codes are written in Python3, and make use of the following libraries:

+ numpy 1.15.3
+ sympy 1.3
+ mpmath 1.1.0

## Overview of the code

### radial_convergence_plot.py

This code was used to create some plots described in the article:

![Alt text](images/Fig5.png?raw=true "Example of plot")

For instance, this plot shows <img src="https://latex.codecogs.com/gif.latex?\zeta(\frac{1}{2}-2i)" title="\zeta(\frac{1}{2}-2i)" /> in the complex plane (<font color="green">green</font>), the original deverging series <img src="https://latex.codecogs.com/gif.latex?\left\{&space;\zeta_n(\frac{1}{2}-2i)&space;\right\}_{n&space;\in&space;\mathbb{N}}" title="\left\{ \zeta_n(\frac{1}{2}-2i) \right\}_{n \in \mathbb{N}}" /> (<font color="red">red</font>), and the converging series <img src="https://latex.codecogs.com/gif.latex?\left\{&space;c_n(\frac{1}{2}-2i)&space;\right\}_{n&space;\in&space;\mathbb{N}}" title="\left\{ \zeta_n(\frac{1}{2}-2i) \right\}_{n \in \mathbb{N}}" /> (<font color="blue">blue</font>). In the preprint, this last series is shown to converge to the actual value of zêta.

### Uz_calculator.py

This code is meant to compute the approximated expression of *Uz* for first known non-trivial zeros of the Riemann zêta function. The table of these values is in the preprint, here are the first values: 8, 14, 18, 24, 18, 32, 38, 40, 46...

### domination_function.py

This code was used to generate an upper and lower bounds for the function <img src="https://latex.codecogs.com/gif.latex?d_{x,&space;y}(n)"> described in the preprint. The expression of these bounds can be written as:

<p align="center">
<img src="https://latex.codecogs.com/gif.latex?\frac{&space;\sum_{i=0}^{27}&space;n^ip_i(x,&space;y)}{&space;qn^{28}&space;&plus;&space;\sum_{i=0}^{27}&space;n^iq_i(x,&space;y)&space;}&space;<&space;d_{x,&space;y}(n)&space;<&space;\frac{&space;\sum_{i=0}^{26}&space;n^iP_i(x,&space;y)}{&space;Qn^{27}&space;&plus;&space;\sum_{i=0}^{26}&space;n^iQ_i(x,&space;y)&space;}" title="\frac{ \sum_{i=0}^{27} n^ip_i(x, y)}{ qn^{28} + \sum_{i=0}^{27} n^iq_i(x, y) } < d_{x, y}(n) < \frac{ \sum_{i=0}^{26} n^iP_i(x, y)}{ Qn^{27} + \sum_{i=0}^{26} n^iQ_i(x, y) }" />
</p>

Where *q* and *Q* are relative integers. Here is a plot of such a bounding.

![Alt text](images/Fig8.png?raw=true "Upper and lower bound for d")
