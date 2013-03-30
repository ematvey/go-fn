// Copyright 2012 - 2013 The Fn Authors. All rights reserved. See the LICENSE file.

package fn

//The Gamma function and relatives.

import (
	"math"
)

//Natural logarithm of the Gamma function
func LnΓ(x float64) (res float64) {
	res = (x-0.5)*math.Log(x+4.5) - (x + 4.5)
	res += logsqrt2pi
	res += math.Log1p(
		76.1800917300/(x+0) - 86.5053203300/(x+1) +
			24.0140982200/(x+2) - 1.23173951600/(x+3) +
			0.00120858003/(x+4) - 0.00000536382/(x+5))

	return
}

/*
//Upper incomplete Gamma function	DOES NOT WORK FOR FLOAT, ONLY INT S, needs to be reimplemented
func IΓ(s, x float64) float64 { 
	if s < 0 {
		return 1
	}
	return (s-1) * IΓ(s-1, x) + math.Pow(x, s-1) * math.Exp(-x)
}
*/

//Upper incomplete Gamma function	// did not pass test for IΓ(1.45896, 3.315) == 0.0706743424609074192334
func IΓ(s, x float64) float64 {
	return IGamC(s, x) * Γ(s)
}

//Upper incomplete Gamma function for integer s only
func IΓint(s int64, x float64) float64 {
	if s < 0 {
		return 1
	}
	return float64(s-1)*IΓint(s-1, x) + math.Pow(x, float64(s-1))*math.Exp(-x)
}

/*
//Lower incomplete Gamma function   BUGGY!!!
func Iγ(s, x float64) float64 { 
	if s < 0 {
		return 1
	}
	return (s-1) * Iγ(s-1, x) - math.Pow(x, s-1) * math.Exp(-x)
}
*/

//Lower incomplete Gamma function
func Iγ(s, x float64) float64 {
	if s < 0 {
		return 1
	}
	return IGam(s, x) * Γ(s)
}

//Lower incomplete Gamma function for integer s only
func Iγint(s int64, x float64) float64 {
	if s < 0 {
		return 1
	}
	return Γ(float64(s)) - IΓint(s, x)
}

// Regularized Gamma function
func Γr(s, x float64) float64 {
	return Iγ(s, x) / Γ(s)
}

func GammaP(p int, x float64) (r float64) {
	pf := float64(p)
	r = math.Pow(math.Pi, 0.25*pf*(pf-1))
	for j := float64(1); j <= pf; j++ {
		r *= GammaF(x + .5*(1-j))
	}
	return
}

func LnGammaP(p int, x float64) (r float64) {
	pf := float64(p)
	r = pf * (pf - 1) * .25 * math.Log(math.Pi)
	for j := float64(1); j <= pf; j++ {
		r += LnΓ(x + .5*(1-j))
	}
	return
}

func GammaPRatio(p int, x, y float64) (r float64) {
	pf := float64(p)
	for j := float64(1); j <= pf; j++ {
		r *= GammaF(x + .5*(1-j))
		r /= GammaF(y + .5*(1-j))
	}
	return
}

//LnΓp(x)/LnΓp(y)
func LnGammaPRatio(p int, x, y float64) (r float64) {
	pf := float64(p)
	for j := float64(1); j <= pf; j++ {
		r += LnΓ(x + .5*(1-j))
		r -= LnΓ(y + .5*(1-j))
	}
	return
}
