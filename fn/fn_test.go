package fn

import (
	"testing"
	"fmt"
//	"math"
)

// test against known values

func TestUpperIncompleteGamma(t *testing.T) {
	x:= IΓ(1, 2.0)
	y:= 0.13533528  

	if !check(x, y){
		t.Error()
		fmt.Println(x, y)
	}


	x= IΓ(1.45896, 3.315)
	y= 0.07067434   

	if !check(x, y){
		t.Error()
		fmt.Println(x, y)
	}

	x= IΓ(2, 6)
	y= 0.01735127 

	if !check(x, y){
		t.Error()
		fmt.Println(x, y)
	}

	x= IΓ(2, 6.5)
	y= 0.01127579  

	if !check(x, y){
		t.Error()
		fmt.Println(x, y)
	}
}

func TestLowerIncompleteGamma(t *testing.T) {
	x:= Iγ(1, 2.0)
	y:= 0.86466472

	if !check(x, y){
		t.Error()
		fmt.Println(x, y)
	}

}
