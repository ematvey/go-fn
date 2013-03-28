package fn

import (
	"fmt"
	"testing"
)

// Test against R:choose

func TestFChoose(t *testing.T) {
	const delta = 1e-4
n, k :=15.3, 4.0

	fmt.Println("Testing FChoose")
	a := FChoose(n, k)
	x := 1491.327

	if abs(x-a) > delta {
		fmt.Println("failed: ",  a, x)
		t.Error()
	}
}
