nums = []
with open("primpoly.txt") as file:
	content = file.read()
	separated = content.split(", ")
	nums = [ int(num) for num in separated ]
	nums.sort()

with open("prim.txt", 'w') as outfile:
	for num in nums:
		outfile.write(str(num)+"\n")