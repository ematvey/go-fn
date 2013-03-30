// Copyright 2012 - 2013 The Fn Authors. All rights reserved. See the LICENSE file.

package fn

// Factorials.

import (
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

//PartialFact returns Fact(n)/Fact(m)
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
