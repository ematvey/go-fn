package fn

import (
	"fmt"
	"testing"
)

// Test against R:choose

func TestFChoose(t *testing.T) {
	const delta = 1e-4

	fmt.Println("Testing FChoose #1")
	n, k :=15.3, 4.0
	a := FChoose(n, k)
	x := 1491.327

	if abs(x-a) > delta {
		fmt.Println("failed: ",  a, x)
		t.Error()
	}

	fmt.Println("Testing FChoose #2")
	n, k =21.3234, 6.0
	a = FChoose(n, k)
	x = 60263.41

	if abs(x-a) > delta {
		fmt.Println("failed: ",  a, x)
		t.Error()
	}

	fmt.Println("Testing FChoose #3")
	n, k =21.3234, 34.0
	a = FChoose(n, k)
	x = 2.690386e-11

	if abs(x-a) > delta {
		fmt.Println("failed: ",  a, x)
		t.Error()
	}

	fmt.Println("Testing FChoose #4")
	n, k =0.77, 34.0
	a = FChoose(n, k)
	x = -0.0003862904

	if abs(x-a) > delta {
		fmt.Println("failed: ",  a, x)
		t.Error()
	}
}
