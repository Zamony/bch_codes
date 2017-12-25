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
		return self._bin[index]

	def __add__(self, polynom):
		return Polynom( self._num ^ polynom.num )

	def __sub__(self, polynom):
		return self.__add__(polynom)

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


print(Polynom(5)+Polynom(1))

