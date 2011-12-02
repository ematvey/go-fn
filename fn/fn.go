package fn

import (
	"fmt"
	"math"
)

//Fact(n) = n*Fact(n-1)
func Fact(n int64) int64 {
	return PartialFact(n, 0)
}
//LnFact(n) = math.Log(n)+LnFact(n-1)
func LnFact(n int64) float64 {
	return LnPartialFact(n, 0)
}

//returns Fact(n)/Fact(m)
func PartialFact(n int64, m int64) int64 {
	if n == m {
		return 1
	}
	return n * PartialFact(n-1, m)
}

//returns LnFact(n)-LnFact(m)
func LnPartialFact(n int64, m int64) float64 {
	if n == m {
		return 0
	}
	return math.Log(float64(n)) + LnPartialFact(n-1, m)
}

func Choose(n int64, i int64) int64 {
	smaller := i
	if n-i < smaller {
		smaller = n - i
	}
	return PartialFact(n, smaller) / Fact(smaller)
}

func LnChoose(n int64, i int64) float64 {
	smaller := i
	if n-i < smaller {
		smaller = n - i
	}
	return LnPartialFact(n, smaller) - LnFact(smaller)
}

func ChooseMany(i []int64) int64 {
	return 0
}

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

//The Gamma function
var Γ = math.Gamma
var GammaF = math.Gamma

var sqrt2pi = math.Sqrt(2 * math.Pi)
var logsqrt2pi = math.Log(math.Sqrt(2 * math.Pi))

//Natural logarithm of the Gamma function
func LnΓ(x float64) (res float64) {
	res = (x - 0.5) * math.Log(x+4.5) - (x + 4.5)
	res += logsqrt2pi
	res += math.Log(1.0 +
		76.1800917300/(x+0) - 86.5053203300/(x+1) +
		24.0140982200/(x+2) - 1.23173951600/(x+3) +
		0.00120858003/(x+4) - 0.00000536382/(x+5))

	return
}

//Upper incomplete Gamma function
func IΓ(s, x float64) float64 { 
	if s < 0 {
		return 1
	}
	return (s-1) * IΓ(s-1, x) + math.Pow(x, s-1) * math.Exp(-x)
}

//Lower incomplete Gamma function
func Iγ(s, x float64) float64 { 
	if s < 0 {
		return 1
	}
	return (s-1) * Iγ(s-1, x) - math.Pow(x, s-1) * math.Exp(-x)
}

//Beta function
func B(x float64, y float64) float64 {
	return Γ(x) * Γ(y) / Γ(x+y)
}

//Non regularized incomplete Beta function
func IB(a, b, x float64) float64 { 
	return BetaIncReg(a, b, x) * math.Exp(LnΓ(a) + LnΓ(b) - LnΓ(a + b))
}

//Regularized incomplete Beta function
func BetaIncReg(α float64, β, x float64) float64 { 
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

//LogBeta function
func LnB(x float64, y float64) float64 {
	return LnΓ(x) + LnΓ(y) - LnΓ(x+y)
}


func GammaP(p int, x float64) (r float64) {
	pf := float64(p)
	r = math.Pow(math.Pi, 0.25 * pf * (pf - 1))
	for j := float64(1); j <= pf; j++ {
		r *= GammaF(x + .5*(1-j))
	}
	return
}

var LnΓp = LnGammaP
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
var LnΓpRatio = LnGammaPRatio
func LnGammaPRatio(p int, x, y float64) (r float64) {
	pf := float64(p)
	for j := float64(1); j <= pf; j++ {
		r += LnΓ(x + .5*(1-j))
		r -= LnΓ(y + .5*(1-j))
	}
	return
}


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

	panic(fmt.Sprintf("betaContinuedFraction(): α or β too big, or maxIter too small"))
	return -1.00
}


