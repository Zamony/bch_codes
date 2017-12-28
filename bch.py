import numpy as np
import gf

def form_bch_matrix(fc, lc):
    matrix = np.zeros( (fc.size, fc.size), dtype="int64" )
    els = np.concatenate([fc, lc[1:]])

    for i in range(fc.size):
        for j in range(fc.size):
            matrix[i][j] = els[i+j]

    return matrix

class BCH:
    def __init__(self, n, t):
        q = np.log2(n + 1).astype(int)
        if q < 2 or q > 16:
            raise ValueError("log2(n + 1) not in [2, 16]")

        primpoly = None
        with open("primpoly2.txt") as file:
            content = file.read()
            separated = content.split(", ")
            for num in separated:
                num = int(num)
                if gf.Polynom(num).power == q:
                    primpoly = num
                    break

        self.n = n
        self.pm = gf.gen_pow_matrix(primpoly)
        self.R = self.pm[:(2 * t), 1]
        self.g = gf.minpoly(self.R, self.pm)[0]
        #print((self.g, self.R, self.pm ))

    def dist(self):
        k = self.n - self.g.size + 1
        U = np.eye(k, dtype="int64")
        V = self.encode(U)

        min_dist = self.n + 1
        for num in range(1, 2**k):
            binary = [int(i) for i in list(bin(num)[2:])]
            binary = [0] * ( k - len(binary) ) + binary
            binary = np.array(binary, dtype="int64").reshape(k, 1)

            dist = np.sum(np.sum(V * binary, axis=0) % 2)
            if dist < min_dist:
                min_dist = dist

        return min_dist

    def encode(self, U):
        k = U.shape[1]
        m = self.g.size - 1
        n_msgs = U.shape[0]
        V = np.zeros((n_msgs, k + m), dtype=int)

        x_m = gf.Polynom( 2**m )
        for i in range(n_msgs):
            s = x_m * gf.Polynom( gf.l2_to_num(U[i, :]) )
            r = ( s / gf.Polynom( gf.l2_to_num(self.g)) )[1]
            v = s + r
            v = np.array(v.binary, dtype="int64")
            V[i, -v.size:] = v.copy()

        return V

    def decode(self, W, method="euclid"):
        n_msgs = W.shape[0]
        n = W.shape[1]
        t = int(self.R.size / 2)
        V = np.zeros((n_msgs, n))

        for msg_idx in range(n_msgs):
            s = gf.polyval(W[msg_idx, :], self.R, self.pm)

            # no errors
            if np.count_nonzero(s) == 0:
                V[msg_idx, :] = W[msg_idx, :]
                continue

            # decode using PGZ
            if method == 'pgz':
                for nu in range(t, 0, -1):
                    A = form_bch_matrix( s[:nu], s[(nu-1) : (2*nu-1)] )
                    b = s[nu : (2*nu)]
                    x = gf.linsolve(A, b, self.pm)

                    if np.all(np.isnan(x) == False):
                        break
                else:
                    W[msg_idx, :] = np.nan
                    continue

                Lambda = np.append(x, [1])

            # else decode using euclid
            else:
                S = np.append(s[::-1], [1])
                z_d = np.zeros((2 * t + 2), dtype="int64")
                z_d[0] = 1
                Lambda = gf.euclid(z_d, S, self.pm, max_deg=t)[2]

            err_pos = np.where(gf.polyval(Lambda, self.pm[:, 1], self.pm) == 0)[0]
            v = W[msg_idx, :].copy()
            v[err_pos] = np.abs(v[err_pos] - 1)

            if np.count_nonzero(gf.polyval(v, self.R, self.pm)) == 0:
                V[msg_idx, :] = v.copy()
            else:
                V[msg_idx, :] = np.nan

        return V.astype(int)

bch = BCH(7, 1)
msg = np.array([[0,0,1,1]], dtype="int64")
print( bch.encode(msg) )
print( bch.decode(np.array([[0, 1, 1, 1, 1, 0, 1]], dtype="int64"), method='pgz') )
print( bch.dist() )

