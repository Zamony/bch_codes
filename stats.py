import random
import numpy as np
import bch

def stat_test(n, t, count_of_tests):
    bch_code = bch.BCH(n, t)
    m = bch_code.g.size - 1
    k = n - m
    msgs = []
    for i in range(2**k):
        msg = [int(bit) for bit in list(bin(i)[2:])]
        msg = [0]*( k - len(msg) ) + msg
        msgs.append(msg)

    msgs = np.array(msgs, dtype="int64")
    enc_msgs = bch_code.encode(msgs)

    wrongly_decoded = 0
    not_decoded = 0
    succ_decoded = 0
    

    for iteration in range(count_of_tests):
        err_enc_msgs = enc_msgs.copy()
        for msg in err_enc_msgs:

            errs_count = random.randint(0, t)
            positions = list( range(n-1) )
            err_poses = random.sample(positions, errs_count)

            for err_pos in err_poses:
                msg[err_pos] = (msg[err_pos] + 1) % 2

        dec_msgs = bch_code.decode(err_enc_msgs)
        for key, msg in enumerate(dec_msgs):
            if np.all( np.isnan(msg) ):
                not_decoded += 1
                continue

            msg = msg.astype(int)
            if np.all( msg == enc_msgs[key] ):
                succ_decoded += 1
            else:
                wrongly_decoded += 1


    return (succ_decoded, wrongly_decoded, not_decoded)


if __name__ == "__main__":
    print( stat_test(15, 3, 30) )