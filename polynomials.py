# polynomial handler
# imports
from dataclasses import dataclass
from typing import Union
from functools import lru_cache
from matplotlib import pyplot as plt
import matplotlib as mpl
import re


def frange(start, stop, step):
	X = [start]
	while start < stop-step:
		start += step
		X += [start]
	return X

def prod(X): 
	x = X[0]
	for i in X[1:]: x*=i
	return x

def msum(X):
	x = X[0]
	for i in X[1:]: x+=i
	return x

SUP = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")
SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
superscript = lambda n: str(n).translate(SUP)
subscript = lambda n: str(n).translate(SUB)
REG = str.maketrans({superscript(n):("" if n==1 else "^"+str(n)) for n in range(10)} | {subscript(n):"_"+str(n) for n in range(10)})
regularscript = lambda n: str(n).translate(REG)

class polynomial:
	def __init__(self, coeffs : list, domain : Union[tuple, list], variable = "x", label = None):
		self.coeffs = coeffs
		self.degree = len(coeffs)-1
		self.variable = variable
		self.domain = domain
		self.label = label

	# shared domain
	def getSharedDomain(self, other):
		a, b = self.domain
		x, y = other.domain
		return (max(a, x), min(b, y))

	# str method
	def __str__(self):
		item = lambda k: str(self.coeffs[k])+f"{self.variable}"+superscript(k)
		ret = " + ".join(map(item, range(self.degree+1))).replace("+ -", "- ").replace(self.variable+superscript(0), "").replace(" 1"+self.variable, " "+self.variable)
		ret = re.sub(" [+-] 0x[⁰¹²³⁴⁵⁶⁷⁸⁹]*", "", ret).replace("0 - ", "-").replace("0 + ", "")
		return ret

	def __repr__(self): return str(self)

	def __add__(self, other):
		if isinstance(other, polynomial):
			residual = self.coeffs[other.degree+1:] if self.degree > other.degree else other.coeffs[self.degree+1:]
			return polynomial(list(map(sum, zip(self.coeffs, other.coeffs)))+residual, self.getSharedDomain(other), self.variable, label=self.label)

	def __sub__(self, other):
		return self + -other

	def __neg__(self):
		return (-1)*self

	def __mul__(self, other):
		if isinstance(other, polynomial):
			c = lambda k: sum([self.coeff(i)*other.coeff(k-i) for i in range(k+1)])
			#'''
			coeffs = [c(k) for k in range(self.degree + other.degree + 1)]
			'''
			coeffs = []
			for k in range(self.degree+other.degree+1):
				C = c(k)
				print(k, [self.coeff(i)*other.coeff(k-i) for i in range(k+1)], C)
				coeffs+=[C]
			#'''
			return polynomial(coeffs, self.getSharedDomain(other), self.variable, label=self.label)

		elif type(other) == Union[int, float, complex]:
			return polynomial([other*c for c in self.coeffs], self.domain, self.variable, label=self.label)

	def __rmul__(self, other):
		if isinstance(other, polynomial):
			coeffs = [sum([self.coeff(i)*other.coeff(k-i) for i in range(k+1)]) for k in range(self.degree + other.degree + 1)]
			return polynomial(coeffs, self.getSharedDomain(other), self.variable)
		elif type(other) in [int, float, complex]:
			return polynomial([other*c for c in self.coeffs], self.domain, self.variable, label=self.label)

	def __pow__(self, ind : int):
		if ind == 0: return polynomial([1], self.domain, self.variable, label=self.label)
		P = self
		for _ in range(ind - 1):
			P *= self
		return P

	
	# retrieve x^n coefficient
	def coeff(self, n : int) -> Union[int, float, complex]:
		return 0 if self.degree < n else self.coeffs[n]

	# evaluate P(x) at x
	def evaluate(self, x):
		a, b = self.domain
		if not (a <= x <= b): raise ValueError("input not in domain")
		return sum([self.coeffs[k]*(x**k) for k in range(len(self.coeffs))])

	# find ∂ₓP(x)
	def derivative(self):
		C = [self.coeffs[idx+1]*(idx+1) for idx, _ in enumerate(self.coeffs[:-1])]
		return polynomial(C, self.domain, self.variable, label=r"\frac{d}{d"+self.variable+"}"+self.label)

	# find ∫ P(x)dx
	def integral(self, c):
		C = [c if idx == 0 else self.coeffs[idx-1]/idx for idx, _ in enumerate(self.coeffs)]
		return polynomial(C, self.domain, self.variable, label=r"\int"+self.label+f"d{self.variable}")

	# newton-raphson method starting at x0
	def newton(self, x0, iterations=5):
		for _ in range(iterations):
			x0 = x0 - self.evaluate(x0)/self.derivative().evaluate(x0)
		return x0

	# durand-kerner algorithm
	def solve(self,iterations=20):
		roots = [complex(.4, .9)**k for k in range(len(self.coeffs)-1)]
		new_roots = roots

		for _ in range(iterations):
			for idx, root in enumerate(roots):
				new_roots[idx] -= self.evaluate(root)/prod([root - roots[k] if k!=idx else 1 for k in range(len(roots))])
			roots = new_roots

		return [root.real for root in roots]

	# finding (p o q)(x)
	def compose(self, other):
		S = polynomial([0], self.variable)
		for k in range(other.degree+1):
			print(a:=P*polynomial([0 if i < k else Q.coeff(k) for i in range(k+1)]))
			S += a
		return S

	# plotting
	def plot(self, ax=None, X : Union[list, type(None)] = None):
		if not ax : ax = plt.subplot()
		if X == None: X = frange(self.domain[0], self.domain[1], 0.01)
		Y = list(map(self, X))
		ax.plot(X, Y, label=str(self) if not self.label else self.label)
		title = regularscript(self.label if self.label else self)
		title = re.sub(r"_(?P<first>[0123456789])_(?P<second>[0123456789])", r"_{\g<first>\g<second>}", title)
		title = re.sub(r"\^(?P<first>[0123456789])\^(?P<second>[0123456789])", r"\^{\g<first>\g<second>}", title)
		ax.set_title(r"$f(x) = "+title+r"$")

	# making "f(x)" the same as "f.evaluate(x)"
	__call__ = evaluate

# chebyshev polynomial of the first kind
@lru_cache
def T(n): 
	if not n: return polynomial([1], (-1, 1), label="T"+subscript(0)+"(x)")
	if n == 1: return polynomial([0, 1], (-1, 1), label="T"+subscript(1)+"(x)")
	return polynomial([0, 2], (-1, 1), label="T"+subscript(n)+"(x)")*T(n-1) - T(n-2)

# chebyshev polynomial of the second kind
@lru_cache
def U(n): 
	if not n: return polynomial([1], (-1, 1), label="U"+subscript(0)+"(x)")
	if n == 1: return polynomial([0, 2], (-1, 1), label="U"+subscript(1)+"(x)")
	return polynomial([0, 2], (-1, 1), label="U"+subscript(n)+"(x)")*U(n-1) - U(n-2)

# runtime
if __name__ == "__main__":
	f = polynomial([4, -3, 0, 1], (-10, 10))
	g = polynomial([1, 0, 0, 1], (-10, 10))
	print(f*g)
