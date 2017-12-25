class Polynom:
	def __init__(self, num):
		if num < 0 :
			raise ValueError("Got negative number in Polynom initialization")
		self._num = num
		self._bin  = bin(num)[2:]

	@property
	def num(self):
		return self._num

	def __add__(self, polynom):
		return Polynom( self._num ^ polynom.num )

	def __sub__(self, polynom):
		return self.__add__(polynom)

	def __repr__(self):
		output = ""
		max_power = len(self._bin)
		for digit, power  in zip( self._bin, range(max_power-1, 0, -1) ):
			if digit == '1':
				if power != 1:
					output += "x^{}+".format(power)
				else:
					output += "x+"

		output += "" if self._bin.endswith('0') else "1+"
		
		if output == "":
			return "0"
		else:
			return output[:-1]


print(Polynom(5)+Polynom(2))

# with open("primpoly.txt") as file:
# 	content = file.read()
# 	separated = content.split(", ")
# 	nums = tuple(( int(num) for num in separated ))
# 	print( nums[0] )
