// Copyright 2012 - 2013 The Fn Authors. All rights reserved. See the LICENSE file.

package fn

//The Beta function and relatives.

import (
	"math"
)

func betaContinuedFraction(α, β, x float64) float64 {

	var aa, del, res, qab, qap, qam, c, d, m2, m, acc float64
	var i int64
	const eps = 2.2204460492503131e-16
	const maxIter = 1000000000

	acc = 1e-16
	qab = α + β
	qap = α + 1.0
	qam = α - 1.0
	c = 1.0
	d = 1.0 - qab*x/qap

	if math.Abs(d) < eps {
		d = eps
	}
	d = 1.0 / d
	res = d

	for i = 1; i <= maxIter; i++ {
		m = (float64)(i)
		m2 = 2 * m
		aa = m * (β - m) * x / ((qam + m2) * (α + m2))
		d = 1.0 + aa*d
		if math.Abs(d) < eps {
			d = eps
		}
		c = 1.0 + aa/c
		if math.Abs(c) < eps {
			c = eps
		}
		d = 1.0 / d
		res *= d * c
		aa = -(α + m) * (qab + m) * x / ((α + m2) * (qap + m2))
		d = 1.0 + aa*d
		if math.Abs(d) < eps {
			d = eps
		}
		c = 1.0 + aa/c
		if math.Abs(c) < eps {
			c = eps
		}
		d = 1.0 / d
		del = d * c
		res *= del
		if math.Abs(del-1.0) < acc {
			return res
		}
	}

	panic("betaContinuedFraction(): α or β too big, or maxIter too small")
	return -1.00
}

//B returns the Beta function.
func B(x float64, y float64) float64 {
	return Γ(x) * Γ(y) / Γ(x+y)
}

//LogBeta function
func LnB(x float64, y float64) float64 {
	return LnΓ(x) + LnΓ(y) - LnΓ(x+y)
}

//BetaIncReg returns the Non-regularized incomplete Beta function.
func IB(a, b, x float64) float64 {
	return BetaIncReg(a, b, x) * math.Exp(LnΓ(a)+LnΓ(b)-LnΓ(a+b))
}

//BetaIncReg returns the Regularized incomplete Beta function.
func BetaIncReg(α, β, x float64) float64 {
	var y, res float64
	y = math.Exp(LnΓ(α+β) - LnΓ(α) - LnΓ(β) + α*math.Log(x) + β*math.Log(1.0-x))
	switch {
	case x == 0:
		res = 0.0
	case x == 1.0:
		res = 1.0
	case x < (α+1.0)/(α+β+2.0):
		res = y * betaContinuedFraction(α, β, x) / α
	default:
		res = 1.0 - y*betaContinuedFraction(β, α, 1.0-x)/β
	}
	return res
}
