import math
import sys
import traceback
from gmpy2 import mpz
from gmpy2 import isqrt
from time  import time


class PiChudnovsky:
    # Constants used in the algorithm
    A = 13591409
    B = 545140134
    C = 640320
    D = 426880
    E = 10005
    C3_24  = C ** 3 // 24
    # Number of digits per term in the series
    DIGITS_PER_TERM = math.log(53360 ** 3) / math.log(10)  #=> 14.181647462725476
    # Name of the file to write the computed value of pi
    FILENAME = "pi.txt"

    def __init__(self, digits):
        """ Initialization
        :param int digits: digits of PI computation
        """
        self.digits = digits
        # Number of terms required for desired digits of pi
        self.n      = math.floor(self.digits / self.DIGITS_PER_TERM + 1)
        # Precision of the computation based on number of digits
        self.prec   = math.floor((self.digits + 1) * math.log2(10))

    def compute(self):
        """ Computation """
        try:
            tm_s = time() # Start timer
            # Compute P, Q and T using binary splitting algorithm
            p, q, t = self.__bsa(0, self.n)
            # Compute the square root of (E * 10^(2*digits))
            one_sq = mpz(10) ** (2 * self.digits)
            sqrt_c = isqrt(self.E * one_sq)
            # Compute pi using the formula and write it to the file
            pi = (q * self.D * sqrt_c) // t
            with open(self.FILENAME, "w") as f:
                f.write(str(pi))
            # Return elapsed time
            return time() - tm_s
        except Exception as e:
            raise

    def __bsa(self, a, b):
        """ PQT computation by BSA(= Binary Splitting Algorithm)
        :param int a: positive integer
        :param int b: positive integer
        :return list [int p_ab, int q_ab, int t_ab]
        """
        try:
            if a + 1 == b:
                if a == 0:
                    p_ab = q_ab = mpz(1)
                else:
                    # Compute P and Q for the current value of a
                    p_ab = mpz((6 * a -5) * (2 * a - 1) * (6 * a - 1))
                    q_ab = mpz(a * a * a * self.C3_24)
                # Compute T for the current value of a
                t_ab = p_ab * (self.A + self.B * a)
                if a & 1:
                    t_ab *= -1
            else:
                # Divide the interval [a,b] into two halves and compute P, Q and T for each half recursively
                m = (a + b) // 2
                p_am, q_am, t_am = self.__bsa(a, m)
                p_mb, q_mb, t_mb = self.__bsa(m, b)
                p_ab = p_am * p_mb
                q_ab = q_am * q_mb
                t_ab = q_mb * t_am + p_am * t_mb
            return [p_ab, q_ab, t_ab]
        except Exception as e:
            raise


if __name__ == '__main__':
    try:
        if len(sys.argv) < 2:
            digits = pow(10,8)
        else:
            digits = int(sys.argv[1])
        print("#### PI COMPUTATION ( {} digits )".format(digits))
        obj = PiChudnovsky(digits)
        tm = obj.compute()
        print("  Output  file:", "pi.txt")
        print("  Elapsed time: {} seconds".format(tm))
    except Exception as e:
        traceback.print_exc()
        sys.exit(1)