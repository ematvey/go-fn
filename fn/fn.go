// Copyright 2012 - 2013 The Fn Authors. All rights reserved. See the LICENSE file.

package fn

// Special math functions
/* To do:
Better implementation of Riemann Zeta
*/

import (
	"math"
)

var log func(float64) float64 = math.Log
var abs func(float64) float64 = math.Abs
var floor func(float64) float64 = math.Floor
var isNaN func(float64) bool = math.IsNaN
var exp func(float64) float64 = math.Exp
var negInf float64 = math.Inf(-1)
var Γ = math.Gamma
var GammaF = math.Gamma
var sqrt2pi = math.Sqrt(2 * math.Pi)
var logsqrt2pi = math.Log(math.Sqrt(2 * math.Pi))
var LnΓp = LnGammaP
var LnΓpRatio = LnGammaPRatio

var lanczos_coef []float64 = []float64{
	0.99999999999980993,
	676.5203681218851,
	-1259.1392167224028,
	771.32342877765313,
	-176.61502916214059,
	12.507343278686905,
	-0.13857109526572012,
	9.9843695780195716e-6,
	1.5056327351493116e-7}

func isOdd(k float64) bool {
	if k != 2*floor(k/2.0) {
		return true
	}
	return false
}

func isInt(x float64) bool {
	if abs((x)-floor((x)+0.5)) <= 1e-7 {
		return true
	}
	return false
}

// Round to nearest integer
func Round(x float64) float64 {
	var i float64
	f := math.Floor(x)
	c := math.Ceil(x)
	if x-f < c-x {
		i = f
	} else {
		i = c
	}
	return i
}

// Arithmetic mean
func ArithMean(data *Vector) float64 {
	n := data.L
	sum := 0.0
	for i := 0; i < n; i++ {
		sum += data.Get(i)
	}
	return sum / float64(n)
}

// Geometric mean
func GeomMean(data *Vector) float64 {
	n := data.L
	sum := 0.0
	for i := 0; i < n; i++ {
		sum += math.Log(data.Get(i))
	}
	return math.Exp(sum / float64(n))
}

// Harmonic mean
func HarmonicMean(data *Vector) float64 {
	n := data.L
	sum := 0.0
	for i := 0; i < n; i++ {
		sum += 1.0 / data.Get(i)
	}
	return float64(n) / sum
}

// Generalized mean
func GenMean(data *Vector, p float64) float64 {
	n := data.L
	sum := 0.0
	for i := 0; i < n; i++ {
		sum += math.Pow(data.Get(i), p)
	}
	return math.Pow(sum/float64(n), 1/p)
}

// Bernoulli number
// Akiyama–Tanigawa algorithm for Bn
func Bn(n int64) float64 {
	var m int64
	a := make([]float64, n)
	for m = 0; m <= n; m++ {
		a[m] = 1 / float64(m+1)
		for j := m; j >= 1; j-- {
			a[j-1] = float64(j) * (a[j-1] - a[j])
		}
	}
	return a[0] // (which is Bn)
}

// H returns the generalized harmonic number of order n of m. 
func H(n int64, m float64) float64 {
	var i int64
	h := 0.0
	for i = 1; i <= n; i++ {
		h += math.Pow(float64(i), m)
	}
	return h
}

// Generalized harmonic number
func H2(n int64, q, s float64) float64 {
	var i int64
	h := 0.0
	for i = 1; i <= n; i++ {
		h += math.Pow((float64(i) + q), -s)
	}
	return h
}
