#include "L.h"

void mult_poly_taylor(Complex *poly1, Complex *poly2, Complex *poly_product, int precision_digits)
{
	int j,k;
	for (k=0;k<=precision_digits;k++)
		poly_product[k] = 0;
	
	for (k=0;k<=precision_digits;k++)
		for (j=0;j<=k;j++)
			poly_product[k] = poly_product[k] + poly1[j]*poly2[k-j];
}
