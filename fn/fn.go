// Special math functions
/* To do:
Better implementation of Riemann Zeta
*/
package fn

import (
	"fmt"
	"math"
)

//Fact(n) = n*Fact(n-1)
func Fact(n int64) int64 {
	return PartialFact(n, 0)
}

/*
//LnFact(n) = math.Log(n)+LnFact(n-1)
func LnFact(n int64) float64 {
	var lnF float64
	if n == 0 {
		lnF = math.Log(1)
	} else {
		lnF = LnPartialFact(n, 0)
	}
	return lnF
}
*/

func LnFact(n int64) float64 {
	var (
		lnF float64
		i   int64
	)
	if n == 0 {
		lnF = 0
	} else {

		for i = 1; i <= n; i++ {
			lnF += math.Log(float64(i))
		}
	}
	return lnF
}

//LnFactBig(n) = Gamma(n+1)
func LnFactBig(n int64) float64 {
	return LnΓ(float64(n + 1))
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
	res = (x-0.5)*math.Log(x+4.5) - (x + 4.5)
	res += logsqrt2pi
	res += math.Log(1.0 +
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

//Beta function
func B(x float64, y float64) float64 {
	return Γ(x) * Γ(y) / Γ(x+y)
}

//Non regularized incomplete Beta function
func IB(a, b, x float64) float64 {
	return BetaIncReg(a, b, x) * math.Exp(LnΓ(a)+LnΓ(b)-LnΓ(a+b))
}

//Regularized incomplete Beta function
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

//LogBeta function
func LnB(x float64, y float64) float64 {
	return LnΓ(x) + LnΓ(y) - LnΓ(x+y)
}

func GammaP(p int, x float64) (r float64) {
	pf := float64(p)
	r = math.Pow(math.Pi, 0.25*pf*(pf-1))
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

// Binomial coefficient (in combinatorics, it gives the number of ways, disregarding order, 
// that k objects can be chosen from among n objects; more formally, the number of 
// k-element subsets (or k-combinations) of an n-element set)

func LnBinomCoeff(n, k int64) float64 {
	if k == 0 {
		return math.Log(1)
	}
	if n == 0 {
		panic("n == 0")
	}
	if n < 10 && k < 10 {
		return math.Log(BinomCoeff(n, k))
	}

	// else, use factorial formula
	return LnFactBig(n) - LnFactBig(k) - LnFactBig(n-k)
}

func BinomCoeff(n, k int64) float64 {
	if k == 0 {
		return 1
	}
	if n == 0 {
		return 0
	}
	// if n, k are small, use recursive formula
	if n < 10 && k < 10 {
		return BinomCoeff(n-1, k-1) + BinomCoeff(n-1, k)
	}

	// else, use factorial formula
	fmt.Println(LnFactBig(n), LnFactBig(k), LnFactBig(n-k))
	return Round(math.Exp(LnFactBig(n) - LnFactBig(k) - LnFactBig(n-k)))
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

// Harmonic mean
func HarmonicMean(data *Vector) float64 {
	n := data.L
	sum := 0.0
	for i := 0; i < n; i++ {
		sum += 1.0 / data.Get(i)
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

// Arithmetic mean
func ArithMean(data *Vector) float64 {
	n := data.L
	sum := 0.0
	for i := 0; i < n; i++ {
		sum += data.Get(i)
	}
	return sum / float64(n)
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

// Riemann zeta function ζ(s)
func RiemannZeta(s float64) float64 {
	var (
		x float64
		b []float64 = []float64{
			1.6666666666666666666666666666666666666666666666667e-1,
			-3.3333333333333333333333333333333333333333333333333e-2,
			2.3809523809523809523809523809523809523809523809524e-2,
			-3.3333333333333333333333333333333333333333333333333e-2,
			7.5757575757575757575757575757575757575757575757576e-2,
			-2.5311355311355311355311355311355311355311355311355e-1,
			1.1666666666666666666666666666666666666666666666667e+0,
			-7.0921568627450980392156862745098039215686274509804e+0}
	)

	switch s {
	case 0:
		x = -0.5
	case 1:
		x = math.NaN()
	case 1.5:
		x = 2.61237534868548834334856756792407163057080065240006340757332824881492776768827286099624386812631195238297
	case 2:
		x = 1.64493406684822643647241516664602518921894990120679843773555822937000747040320087383362890061975870
	case 3:
		x = 1.20205690315959428539973816151144999076498629234049888179227155534183820578631309018645587360933525814619915
	case 4:
		x = 1.08232323371113819151600369654116790277475095191872690768297621544412061618696884655690963594169991
	case 5:
		x = 1.03692775514336992633136548645703416805708091950191281197419267790380358978628148456004310655713333
	case 6:
		x = 1.01734306198444913971451792979092052790181749003285356184240866400433218290195789788277397793853517
	case 7:
		x = 1.00834927738192282683979754984979675959986356056523870641728313657160147831735573534609696891385132
	case 8:
		x = 1.00407735619794433937868523850865246525896079064985002032911020265258295257474881439528723037237197
	case 9:
		x = 1.00200839282608221441785276923241206048560585139488875654859661590978505339025839895039306912716958
	case 10:
		x = 1.00099457512781808533714595890031901700601953156447751725778899463629146515191295439704196861038565

	default:
		a := 12
		aa := float64(a)
		k := 8
		x = 0.0

		for i := 1; i < a; i++ {
			ii := float64(i)
			x += 1.0 / math.Pow(ii, s)
		}

		x += 1.0/((s-1)*math.Pow(aa, (s-1))) + 1/(2*math.Pow(aa, s))
		term := (s / 2) / math.Pow(aa, (s+1))
		x += term * b[0]

		for i := 1; i < k; i++ {
			ii := float64(i) + 1
			term *= (s + 2*ii - 2) * (s + 2*ii - 3) / (aa * aa * 2 * ii * (2*ii - 1))
			x += term * b[i]
		}
	}
	return x
}
