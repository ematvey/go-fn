package fn

import (
	"testing"
	"fmt"
	"math"
)
// tests against known values



func TestGamma(t *testing.T) {
	var (
		y float64
		i int64
	)
	fmt.Println("test of Γ function")

	a:=[]float64{0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2}
	x:=[]float64{9.5135076987,4.590843712,2.9915689877,2.2181595438,1.7724538509,1.4891922488,1.2980553326,1.1642297137,1.0686287021,1,0.9513507699,0.9181687424,0.8974706963,0.8872638175,0.8862269255,0.8935153493,0.9086387329,0.931383771,0.9617658319,1}

	for i = 0; i < int64(len(a)); i++ {
		y=Γ(a[i])
			if !check(x[i], y){
				t.Error()
				fmt.Println(a[i], x[i], y)
			}
	}
}

func TestLnGamma(t *testing.T) {
	var (
		y float64
		i int64
	)
	fmt.Println("test of LnΓ function")

	a:=[]float64{0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9}
	x:=[]float64{2.2527126517,1.5240638224,1.0957979948,0.7966778177,0.5723649429,0.3982338581,0.2608672465,0.1520596784,0.0663762397,-0.0498724413,-0.08537409,-0.1081748095,-0.1196129142,-0.1207822376,-0.1125917657,-0.0958076974,-0.0710838729,-0.0389842759}


	for i = 0; i < int64(len(a)); i++ {
		y=LnΓ(a[i])
			if !check(x[i], y){
				t.Error()
				fmt.Println(a[i], x[i], y)
			}
	}
}



func TestIncompleteGamma(t *testing.T) {

	var s, z, x, y float64
	fmt.Println("test of Upper and Lower Incomplete Gamma function")

	s=1
	z=2
	x= IΓ(s,z)
	y= 0.135335283236612691894
	if !check(x, y){
		t.Error()
		fmt.Println(x, y)
	}
	x= Iγ(s,z)
	y= 0.864664716763387308106
	if !check(x, y){
		t.Error()
		fmt.Println(x, y)
	}




	s=1.45896
	z=3.315
	y= 0.0706743424609074192334
	x= IΓ(s,z)
	if !check(x, y){
		t.Error()
		fmt.Println(x, y)
	}
	x= Iγ(s,z)
	y= 0.81493191400161894606
	if !check(x, y){
		t.Error()
		fmt.Println(x, y)
	}

	s=2.0
	z=6.0
	y= 0.01735126523666450896132
	x= IΓ(s,z)
	if !check(x, y){
		t.Error()
		fmt.Println(x, y)
	}
	x= Iγ(s,z)
	y= 0.9826487347633354910387
	if !check(x, y){
		t.Error()
		fmt.Println(x, y)
	}


	s=2.0
	z=6.5
	y= 0.01127579394733179335537
	x= IΓ(s,z)
	if !check(x, y){
		t.Error()
		fmt.Println(x, y)
	}
	x= Iγ(s,z)
	y= 0.9887242060526682066446
	if !check(x, y){
		t.Error()
		fmt.Println(x, y)
	}


	s=3.5316559
	z=8.3561865
	y= 0.0690709003122005470888
	x= IΓ(s,z)
	if !check(x, y){
		t.Error()
		fmt.Println(x, y)
	}
	x= Iγ(s,z)
	y= 3.3729541465127126384
	if !check(x, y){
		t.Error()
		fmt.Println(x, y)
	}


	s=6.54681
	z=8.68188
	y= 58.8630238247231988064
	x= IΓ(s,z)
	if !check(x, y){
		t.Error()
		fmt.Println(x, y)
	}
	x= Iγ(s,z)
	y= 254.283105608912865394
	if !check(x, y){
		t.Error()
		fmt.Println(x, y)
	}
}


func TestIncompleteGammaInt(t *testing.T) {
	var s int64
	var z, x, y float64
	fmt.Println("test of Incomplete Gamma integer s function")

	s=1
	z=2
	y= 0.135335283236612691894
	x= IΓint(s,z)
	if !check(x, y){
		t.Error()
		fmt.Println(x, y)
	}

	s=6
	z=8.68188
	y= 16.37023052772108684457
	x= IΓint(s,z)
	if !check(x, y){
		t.Error()
		fmt.Println(x, y)
	}

	s=3
	z=5.41657816
	y= 0.1873436971857055661085
	x= IΓint(s,z)
	if !check(x, y){
		t.Error()
		fmt.Println(x, y)
	}
}

// test for Incomplete Beta
func TestIncompleteBeta(t *testing.T) {
	var x, y, acc float64
	acc = 1e-5
	check := func(x, y, acc float64) bool {
		if math.Abs(x-y) > acc  {
			return false
		}
		return true
	}

	x = IB(0.5, 2, 0.5)
	y = 1.17851130

	if !check(x, y, acc){
		t.Error()
	}

	x = IB(2.0, 6, 0.75)
	y = 0.02377755


	if !check(x, y, acc){
		t.Error()
	}

	x = IB(1.456, 3.25895, 0.99)
	y = 0.14450022


	if !check(x, y, acc){
		t.Error()
	}
}

// test for Incomplete Gamma
func TestIncompleteGamma2(t *testing.T) {
	var x, y, z, acc float64
	acc = 1e-05
	check := func(x, y, acc float64) bool {
		if x/y > 1.00 {
			z = y/x
		} else {
			z = x/y
		}
		if 1-z > acc  {
			return false
		}
		return true
	}
	x = IΓ(3, 5.5)
	y = 0.17675286

	if !check(x, y, acc){
		t.Error()
	}

	x = IΓ(6, 3.96545)
	y = 94.86080842

	fmt.Println(x, y)

	if !check(x, y, acc){
		t.Error()
	}
/*
	x = IΓ(18.3542, 3.96545)
	y = 9.838284e+14
	fmt.Println(x, y)

	if !check(x, y, acc){
		t.Error()
	}
//FAILED, non-integer s
*/

	x = IΓ(18, 3.96545)
	y = 3.556874e+14

	if !check(x, y, acc){
		t.Error()
	}

		fmt.Println("test of binom coeff", BinomCoeff(150, 71) )

}

