# Legendre Polynomial

This repository contains a Python script for projecting functions onto a set of normalized Legendre polynomials using numerical integration.

## Mathematical Overview

The script performs projections of a function onto a set of orthogonal polynomials, specifically **Legendre polynomials**, over a given interval. 
These polynomials are generated and normalized using the Gram-Schmidt orthogonalization process. 

The Legendre polynomials \( P_n(x) \) are a sequence of orthogonal polynomials that arise in various areas of numerical analysis, physics, and engineering. 
In this script, we perform the projection of a function \( f(x) \) onto these Legendre polynomials. 
The projection process involves calculating the coefficients \( \langle f, P_n \rangle \) using numerical integration, 
and then constructing the projection of the function onto the polynomial basis.

## Functions Implemented

1. **`integration_simpson_sampled(a, b, fx, N)`**:
   - Numerically integrates the function `fx` over the interval \([a, b]\) using **Simpson's rule** with \( N \) subintervals.
   
2. **`get_projection(BASE, f, a, b)`**:
   - Projects the function `f(x)` onto the polynomial basis `BASE` over the interval \([a, b]\).
   - The projection is computed using the previously defined numerical integration function.

3. **`legendre_polynomial_normalized(number_functions, N, a, b)`**:
   - Generates and normalizes a set of Legendre polynomials over the interval \([a, b]\).
   - The normalization is performed using the Gram-Schmidt process to ensure the polynomials are orthogonal and normalized.

4. **`main()`**:
   - The main function demonstrates the projection of two specific functions, \( \sin(x) \) and \( e^{\sqrt{\left| \frac{x}{7} \right|}} \), onto the Legendre polynomial basis. The projections are then visualized along with the original functions.


## **Plt**
<p align="center">
   <img src="https://github.com/user-attachments/assets/694667f8-5808-4c10-9df3-5531f440b1f8" width="30%" />
   <img src="https://github.com/user-attachments/assets/80e67e40-a79a-4c27-87ed-02af847d80b2" width="60%" />
</p>
