# polynomials
Polynomial Handler in Python 3.9+ (supporting matplotlib plotting).

## Definition and artithmetic of polynomials
Polynomials are instances of the `polynomial` class, with parameters `coeffs: list` - list of polynomial coefficients - and `domain: tuple | list`, and optional parameters `variable: char` and `label: str`. For example, the definition `P = polynomial([1, 3, -2], (-10, 10), variable="t", label="f(t)")` refers to the polynomial `f:[-10, 10] -> C` given by `f(t) = 1 + 3t - 2t^2`.

The set of polynomials over the complex numbers as defined in the code forms a group under addition and a monoid under multiplication. Hence, the `__add__`, `__sub__`, `__neg__`, `__mul__`, `__rmul__`, and `__pow__` are all defined for polynomials, where `__add__, __sub__, __mul__, __rmul__` all return polynomials over the shared domain of the two operands.
Division of polynomials, however, is not defined as it does not always result in a polynomial output. The `__call__` method is set to `evaluate`, i.e. you can call a polynomial object as a function of a single variable.

```
P = polynomial([1, 2], [-1, 1], label="f(t)")
print(P(1)) # this outputs the value of f(1) = 1 + 2(1) = 3
```

## Extras
Included in this class is function composition, differentiation, integration, root finding, and plotting (using matplotlib).
Moreover, the Chebyshev polynomials of the first kind `T(n)` and second kind `U(n)` are defined. Their definition is recursive (i.e. to generate `T(3)` you must first generate `T(2)` and `T(1)` &c.), but uses caching.
