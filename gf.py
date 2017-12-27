import pprint
import numpy as np

class Polynom:
	def __init__(self, num):
		if num < 0 :
			raise ValueError("Got negative number in Polynom initialization")
		self._num = num
		self._bin  = [ int(i) for i in list( bin(num)[2:] ) ]
		# Не может быть отрицательной, т.к. числа больше или равны 0
		self._power = len( self._bin ) - 1

	@property
	def num(self):
		return self._num

	@property
	def power(self):
		return self._power

	def bin(self, index):
		return self._bin[self._power - index]

	def __add__(self, polynom):
		return Polynom( self._num ^ polynom.num )

	def __sub__(self, polynom):
		return self.__add__(polynom)

	def __mul__(self, polynom):
		result = Polynom(0)
		for i in range(polynom.power, -1, -1):
			if polynom.bin(i) != 0:
				result = result + Polynom(self._num << i)

		return result

	def __neg__(self):
		return Polynom(self.num)

	def __repr__(self):
		output = ""
		for digit, power  in zip( self._bin, range(self._power, 0, -1) ):
			if digit == 1:
				if power != 1:
					output += "x^{}+".format(power)
				else:
					output += "x+"

		output += "" if self._bin[-1] == 0 else "1+"
		
		if output == "":
			return "0"
		else:
			return output[:-1]


class PolynomF(Polynom):
	def __init__(self, num, table):
		super().__init__(num)
		self._table = table
		self._primpoly = None

	@property
	def primpoly(self):
		if self._primpoly is None:
			power = 1
			for row in table:
				if row[1] != 2**power:
					self._primpoly = 2**power + row[1]
					break
				power += 1

		return self._primpoly

	@property
	def table(self):
		return self._table

	def __add__(self, polynom):
		poly = super().__add__(polynom)
		return PolynomF(poly.num, self.table)

	def __sub__(self, polynom):
		return self.__add__(polynom)

	def __neg__(self):
		poly = super().__neg__()
		return PolynomF(poly.num, self.table)

	def __mul__(self, polyf):
		table = self.table
		poly1num = self.num
		poly2num = polyf.num

		if poly1num == 0 or poly2num == 0:
			return PolynomF(0, self.table)

		m, _ = table.shape
		power1, power2 = None, None
		for key, row in enumerate(table):
			if row[1] == poly1num:
				power1  = key + 1
			if row[1] == poly2num:
				power2  = key + 1

		mul_key = (power1 + power2) % m
		mul_poly = PolynomF(table[mul_key - 1][1], self.table)

		return mul_poly

	def __pow__(self, power):
		if power < 1: raise NotImplementedError("Polynomial negative power isn't implemented")
		
		one_poly = PolynomF(self.num, self.table)
		poly = one_poly
		for i in range(power-1):
			poly = poly * one_poly
		return poly

	def __truediv__(self, polyf):
		table = self.table
		poly2num = polyf.num

		if poly2num == 0:
			raise ZeroDivisionError("Division by zero in polynomial division")

		if self.num == 0:
			return PolynomF(0, self.table)

		m, _ = table.shape
		power = None
		for key, row in enumerate(table):
			if row[1] == poly2num:
				power = key + 1
				break

		div_key = (-power) % m
		mul_poly = PolynomF(table[div_key - 1][1], self.table)

		return self.__mul__(mul_poly)



def gen_pow_matrix(primpoly):

	def pow_in_field(polx, ppoly):
		# Т.к. многочлены примитивные, то примитивный элемент - это x
		polx = polx * Polynom(2)

		if polx.power == ppoly.power:
			polx = polx + Polynom(1 << polx.power)
			polx = polx + ( ppoly + Polynom(1 << ppoly.power) )
		
		return polx

	ppoly = Polynom( primpoly )
	m = 2**(ppoly.power)-1
	table = np.zeros( (m, 2), dtype="int64" )

	px = Polynom(1)
	for i in range(m):
		px = pow_in_field(px, ppoly)
		table[i][1] = px.num
	
	for i in range(m):	
		index = None
		for j in range(m):
			if table[j][1] == i+1:
				index = j+1
				break
		table[i][0] = index


	return table

def add(X, Y):
	z = np.zeros( X.shape, dtype="int64" )
	for i in range(len(X)):
		for j in range(len(Y)):
			z[i][j] = (Polynom(X[i][j]) + Polynom(Y[i][j])).num

	return z

def sum(X, axis=0):
	m, n = X.shape

	result = None
	if axis == 0:
		result = np.zeros(n, dtype="int64")

		for j in range(n):
			poly = Polynom(0)
			for i in range(m):
				poly = poly + Polynom(X[i][j])
			result[j] = poly.num
	else:
		result = np.zeros(m, dtype="int64")

		for i in range(m):
			poly = Polynom(0)
			for j in range(n):
				poly = poly + Polynom(X[i][j])
			result[i] = poly.num

	return result

def prod(X, Y, pm):
	m, n = X.shape
	z = np.zeros( X.shape, dtype="int64" )
	for i in range(m):
		for j in range(n):
			z[i][j] = PolynomF(X[i][j], pm) * PolynomF(Y[i][j], pm)

def divide(X, Y, pm):
	m, n = X.shape
	z = np.zeros( X.shape, dtype="int64" )
	for i in range(m):
		for j in range(n):
			z[i][j] = PolynomF(X[i][j], pm) / PolynomF(Y[i][j], pm)

def linsolve(A, b, pm):
	def pfsum(pair):
		s = PolynomF(0, pm)
		for i in pair:
			s += i
		return s

	b = b.reshape( (len(b), 1) )
	matrix = np.append(A, b, axis=1)
	m = matrix.tolist()

	for i in range(len(m)):
		for j in range(len(m[0])):
			m[i][j] = PolynomF(m[i][j], pm)

	try:
		for col in range(len(m[0])):
			for row in range(col+1, len(m)):
				r = [(rowValue * (-(m[row][col] / m[col][col]))) for rowValue in m[col]]
				m[row] = [pfsum(pair) for pair in zip(m[row], r)]
		ans = []
		m.reverse()
		for sol in range(len(m)):
				if sol == 0:
					ans.append(m[sol][-1] / m[sol][-2])
				else:
					inner = PolynomF(0, pm)
					for x in range(sol):
						inner = inner + (ans[x]*m[sol][-2-x])
					ans.append((m[sol][-1]-inner)/m[sol][-sol-2])
	except ZeroDivisionError:
		return np.nan
	else:
		ans.reverse()
		result = [ poly.num for poly in ans ]
		return result

def minpoly(x, pm):
	def minpoly_b(b, pm):
		orig_b = b
		x = PolynomF(2, pm) # x
		power = 1
		minimal = x - b
		roots = [b.num]
		while True:
			curr_b = b**(2**power)
			#print("{} степень {}".format(curr_b, 2**power))
			if curr_b.num == orig_b.num:
				minimal = Polynom(minimal.num)
				break

			if curr_b.num == x.num:
				minimal = Polynom(orig_b.primpoly)
				roots = [2**(2**i) for i in range(minimal.power)]
				break

			minimal = minimal * (x - curr_b)
			roots.append(curr_b.num)
			power += 1

		return minimal, roots

	roots = []
	minimal = Polynom(1)
	for b in x:
		minimal_pol, roots_pol = minpoly_b(PolynomF(b, pm), pm)
		minimal = minimal * minimal_pol
		roots = roots + roots_pol

	roots = sorted(list(set(roots)))
	return [minimal.bin(i) for i in range(minimal.power+1)][::-1], roots


def polyval(p, x, pm):
	result = []
	for val in x:
		single_result = PolynomF(0, pm)
		poly = PolynomF(val, pm)
		for coeff, power in zip(p, range(len(p))[::-1]):
			if coeff == 1:
				if power == 0:
					single_result = single_result + PolynomF(1, pm)
				else:
					single_result = single_result +  poly**power

		result.append(single_result)

	return result

table = gen_pow_matrix(11)
print( polyval(np.array([1, 1, 0]), np.array([3, 4]), table) )

# A = np.array([ [4, 6, 4], [6, 1, 7], [1, 6, 3] ], dtype="int64")
# b = np.array( [5, 3, 1], dtype="int64" )
# print( linsolve(A, b, table) )

#print( PolynomF(2, table) / PolynomF(5, table) )

# print( Polynom(21) * Polynom(12) )
# print( Polynom(21) )
# print( Polynom(12) )

#pprint.pprint( add( np.array([ [1, 2], [3, 4] ]), np.array([ [4, 3], [2, 1] ]) ) )
#pprint.pprint( sum(np.array([[1, 2], [3, 4]]), 1) )

