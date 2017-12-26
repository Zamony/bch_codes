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
	def __init__(self, num, primpoly, table):
		super().__init__(num)
		self._primpoly = primpoly
		self._table = table

	@property
	def primpoly(self):
		return self._primpoly

	@property
	def table(self):
		return self._table

	def __mul__(self, polyf):
		table = self.table
		poly1num = self.num
		poly2num = polyf.num

		m, _ = table.shape
		power1, power2 = None, None
		for key, row in enumerate(table):
			if row[1] == poly1num:
				power1  = key + 1
			if row[1] == poly2num:
				power2  = key + 1

		mul_key = (power1 + power2) % m

		return table[mul_key - 1][1]

	def __truediv__(self, polyf):
		table = self.table
		poly2num = polyf.num

		m, _ = table.shape
		power = None
		for key, row in enumerate(table):
			if row[1] == poly2num:
				power = key + 1
				break

		div_key = (-power) % m
		mul_poly = PolynomF(table[div_key - 1][1], self.primpoly, self.table)

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


def mulf(table, poly1num, poly2num):
	m, _ = table.shape
	power1, power2 = None, None
	for key, row in enumerate(table):
		if row[1] == poly1num:
			power1  = key + 1
		if row[1] == poly2num:
			power2  = key + 1

	mul_key = (power1 + power2) % m

	return table[mul_key - 1][1] 

def divf(table, poly1num, poly2num):
	m, _ = table.shape
	power = None
	for key, row in enumerate(table):
		if row[1] == poly2num:
			power = key + 1
			break

	div_key = (-power) % m

	return mulf(table, poly1num, table[div_key - 1][1])


def prod(X, Y, pm):
	m, n = X.shape
	z = np.zeros( X.shape, dtype="int64" )
	for i in range(m):
		for j in range(n):
			z[i][j] = mulf(pm, X[i][j], Y[i][j])

def divide(X, Y, pm):
	m, n = X.shape
	z = np.zeros( X.shape, dtype="int64" )
	for i in range(m):
		for j in range(n):
			z[i][j] = divf(pm, X[i][j], Y[i][j])


primpoly = 11
table = gen_pow_matrix(primpoly)
print( PolynomF(4, primpoly, table) / PolynomF(7, primpoly, table) )

# print( Polynom(21) * Polynom(12) )
# print( Polynom(21) )
# print( Polynom(12) )

#pprint.pprint( add( np.array([ [1, 2], [3, 4] ]), np.array([ [4, 3], [2, 1] ]) ) )
#pprint.pprint( sum(np.array([[1, 2], [3, 4]]), 1) )

