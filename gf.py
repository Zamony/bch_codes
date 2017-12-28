import pprint
import numpy as np

class Polynom:
    def __init__(self, num):
        if num < 0 :
            raise ValueError("Got negative number in Polynom initialization")
        self._num = num
        self._bin  = [ int(i) for i in list( bin(num)[2:] ) ]
        self._power = len( self._bin ) - 1

    @property
    def num(self):
        return self._num

    @property
    def power(self):
        return self._power

    @property
    def binary(self):
        return self._bin.copy()

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

    def __truediv__(self, polyf):
        p1 = np.array(self.binary, dtype="int64")
        p2 = np.array(polyf.binary, dtype="int64")
        
        q_deg = p1.size - p2.size
        if q_deg < 0:
            return Polynom(0), Polynom(l2_to_num(p1))

        q = np.zeros(q_deg + 1, dtype=bool)
        while p1.size >= p2.size:
            cur_q_pow = q_deg - (p1.size - p2.size)
            q[cur_q_pow] = 1
            sub_poly = np.asarray((Polynom(l2_to_num(q[cur_q_pow:])) * Polynom(l2_to_num(p2))).binary)
            p1 = np.asarray((Polynom(l2_to_num(p1)) + Polynom(l2_to_num(sub_poly))).binary)
        r = p1
        return Polynom(l2_to_num(q)), Polynom(l2_to_num(r))

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
            for row in self.table:
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


def l2_to_num(binary):
    num = 0
    for i, coeff in enumerate(reversed(binary)):
        num += coeff * (2**i)
    return num

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
    x, y = X.copy(), Y.copy()
    x = x.astype(int)
    y = y.astype(int)
    
    is_vector = len(x.shape) == 1
    if is_vector:
        m = x.size
        n = 1
    else:
        m, n = x.shape

    x = x.flatten()
    y = y.flatten()

    z = np.zeros( m*n, dtype="int64" )
    for i in range( m*n ):
        z[i] = (Polynom(x[i]) + Polynom(y[i])).num

    if is_vector:
        return z
    
    z = z.reshape( (m, n) )
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
    x, y = X.copy(), Y.copy()
    x = x.astype(int)
    y = y.astype(int)
    
    is_vector = len(x.shape) == 1
    if is_vector:
        m = x.size
        n = 1
    else:
        m, n = x.shape

    x = x.flatten()
    y = y.flatten()

    z = np.zeros( m*n, dtype="int64" )
    for i in range( m*n ):
        z[i] = (PolynomF(x[i], pm) * PolynomF(y[i], pm)).num

    if is_vector:
        return z
    
    z = z.reshape( (m, n) )
    return z

def divide(X, Y, pm):
    x, y = X.copy(), Y.copy()
    x = x.astype(int)
    y = y.astype(int)
    
    is_vector = len(x.shape) == 1
    if is_vector:
        m = x.size
        n = 1
    else:
        m, n = x.shape

    x = x.flatten()
    y = y.flatten()

    z = np.zeros( m*n, dtype="int64" )
    for i in range( m*n ):
        z[i] = (PolynomF(x[i], pm) / PolynomF(y[i], pm)).num

    if is_vector:
        return z
    
    z = z.reshape( (m, n) )
    return z

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
        return np.array(result, dtype="int64")

def minpoly(x, pm):
    def minpoly_b(b, pm):
        orig_b = b
        x = PolynomF(2, pm) # x
        power = 1
        roots = [b.num]
        primpoly_power = Polynom(b.primpoly).power

        while True:
            curr_b = b**(2**power)
            if curr_b.num == x.num:
                roots = [(x**(2**i)).num for i in range(primpoly_power)]
                break

            if curr_b.num == orig_b.num:
                break

            roots.append(curr_b.num)
            power += 1

        return roots

    roots = []
    minimal = np.array([1], dtype="int64")
    for b in x:
        roots_pol = minpoly_b(PolynomF(b, pm), pm)
        roots = roots + roots_pol

    roots = sorted(list(set(roots)))
    for root in roots:
        minimal = polyprod(minimal,  np.array([1, root]), pm ) 

    roots = np.array(roots, dtype="int64")

    return minimal, roots


def polyval(p, x, pm):
    result = []
    for val in x:
        single_result = PolynomF(0, pm)
        poly = PolynomF(val, pm)
        for coeff, power in zip(p, range(len(p))[::-1]):
            if coeff != 0:
                if power == 0:
                    single_result = single_result + PolynomF(coeff, pm)
                else:
                    single_result = single_result + PolynomF(coeff, pm) * (poly**power)

        result.append(single_result.num)

    return np.array(result, dtype="int64")


def polyprod(p1, p2, pm):
    if p1.size >= p2.size:
        f1, f2 = p2, p1
    else:
        f1, f2 = p1, p2
    conv = np.zeros(f1.size + f2.size - 1)
    for i in range(f1.size):
        conv[i:(i + f2.size)] = add(conv[i:(i + f2.size)], prod(f1[i] * np.ones(f2.size), f2, pm))
    return conv.astype(int)

def polysum(p1, p2):
    p1, p2 = p1.copy(), p2.copy()
    pow_diff = p1.size - p2.size

    if pow_diff > 0:
        p2 = np.append(np.zeros(pow_diff), p2)
    elif pow_diff < 0:
        p1 = np.append(np.zeros(-pow_diff), p1)

    poly_sum = add(p1, p2)
    nonzero_index = np.nonzero(poly_sum)[0]
    if nonzero_index.size == 0:
        return np.array([0])

    poly_sum = poly_sum[nonzero_index[0]:]

    return poly_sum

def polydiv(p1, p2, pm):
    if p2.size > p1.size:
        return (np.array([0], dtype=int), p1)

    p1 = p1.astype(int).copy()
    q = np.zeros(p1.size - p2.size + 1, dtype=int)

    for i in range(q.size):
        q[i] = divide(np.array([p1[0]]), np.array([p2[0]]), pm)[0]
        poly_prod = polyprod(q[i:], p2, pm)
        if poly_prod.size > p1.size:
            q[i] = 0
        else:
            p1 = polysum(p1, polyprod(q[i:], p2, pm))

    return (q, p1)

def euclid(p1, p2, pm, max_deg=0):
    p1 = p1[np.nonzero(p1)[0][0]:].copy()
    p2 = p2[np.nonzero(p2)[0][0]:].copy()

    if p2.size > p1.size:
        r_prev_prev, r_prev = p2, p1
    else:
        r_prev_prev, r_prev = p1, p2

    x_prev_prev, x_prev = np.array([1]), np.array([0])
    y_prev_prev, y_prev = np.array([0]), np.array([1])

    while True:
        (q, r) = polydiv(r_prev_prev, r_prev, pm)
        x = polysum(x_prev_prev, polyprod(q, x_prev, pm))
        y = polysum(y_prev_prev, polyprod(q, y_prev, pm))
        r_prev_prev, r_prev = r_prev, r
        x_prev_prev, x_prev = x_prev, x
        y_prev_prev, y_prev = y_prev, y

        if (max_deg > 0 and r.size - 1 <= max_deg) or (max_deg == 0 and r.size == 1):
            break

    return (r_prev, x_prev, y_prev)


# table = gen_pow_matrix(11)
# print(table)
# print( minpoly(np.array([2, 4]), table) )

# print( 
#     euclid(
#         np.array([7, 3, 5], dtype="int64"),
#         np.array([2, 0, 4], dtype="int64"),
#         table
#     )
# )

#print( Polynom(3) / Polynom(4) )

#print( minpoly(np.array([6, 4]), table) )

# print( polyval( np.array([2, 5, 6], dtype="int64"), np.array([3, 4], dtype="int64"), table ) )

# A = np.array([ [4, 6, 4], [6, 1, 7], [1, 6, 3] ], dtype="int64")
# b = np.array( [5, 3, 1], dtype="int64" )
# print( linsolve(A, b, table) )

#print( PolynomF(2, table) / PolynomF(5, table) )

# print( Polynom(21) * Polynom(12) )
# print( Polynom(21) )
# print( Polynom(12) )

#pprint.pprint( add( np.array([ [1, 2], [3, 4] ]), np.array([ [4, 3], [2, 1] ]) ) )
#pprint.pprint( sum(np.array([[1, 2], [3, 4]]), 1) )

