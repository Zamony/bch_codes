import pprint

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



def gen_pow_matrix(primpoly):

	def pow_in_field(polx, ppoly):
		# Т.к. многочлены примитивные, то примитивный элемент - это x
		polx = polx * Polynom(2)
		print(polx)

		if polx.power == ppoly.power:
			polx = polx + Polynom(1 << polx.power)
			polx = polx + ( ppoly + Polynom(1 << ppoly.power) )
		
		return polx

	ppoly = Polynom( primpoly )
	m = 2**(ppoly.power)-1
	table = [ [None, None] for _ in range(m) ]

	px = Polynom(1)
	for i in range(m):
		px = pow_in_field(px, ppoly)
		table[i][1] = px
		table[i][0] = px.num

	return table


# print( Polynom(21) * Polynom(12) )
# print( Polynom(21) )
# print( Polynom(12) )
