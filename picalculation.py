# # Initialize denominator
# k = 1

# # Initialize sum
# s = 0

# for i in range(10000000):

# 	# even index elements are positive
# 	if i % 2 == 0:
# 		s += 4/k
# 	else:

# 		# odd index elements are negative
# 		s -= 4/k

# 	# denominator is odd
# 	k += 2

# print(s)

#Program that prints the approximate value of pi
#4/1 - 4/3 + 4/5 - 4/7 + 4/9 - 4/11...
#nth term of the sequence is: \n", (4/(2*n-1)*(-1)**(1-n))) using simple math.
#Code begins

# import math
# def T(n):
#     Tn =(4/(2*n-1)*(-1)**(1-n)) #using some simple math for arithmetic progressions
#     if n==0:
#         return 0
#     if n==1:
#         return 4
#     else:
#         return Tn+T(n-1)
# n=eval(input("Enter the value of n: \n"))
# print ("The approximate value of pi is:", T(n))
# print ("The difference from the actual value of pi is:", T(n)-math.pi)

# #Code ends

# # monte carlo method

# import random

# INTERVAL = 1000

# circle_points = 0
# square_points = 0

# # Total Random numbers generated= possible x
# # values* possible y values
# for i in range(INTERVAL**2):

# 	# Randomly generated x and y values from a
# 	# uniform distribution
# 	# Range of x and y values is -1 to 1
# 	rand_x = random.uniform(-1, 1)
# 	rand_y = random.uniform(-1, 1)

# 	# Distance between (x, y) from the origin
# 	origin_dist = rand_x**2 + rand_y**2

# 	# Checking if (x, y) lies inside the circle
# 	if origin_dist <= 1:
# 		circle_points += 1

# 	square_points += 1

# 	# Estimating value of pi,
# 	# pi= 4*(no. of points generated inside the
# 	# circle)/ (no. of points generated inside the square)
# 	pi = 4 * circle_points / square_points

# ## print(rand_x, rand_y, circle_points, square_points, "-", pi)
# # print("\n")

# print("Final Estimation of Pi=", pi)

import math
import sys
import traceback
from gmpy2 import mpz
from gmpy2 import isqrt
from time  import time


class PiChudnovsky:
    A = 13591409
    B = 545140134
    C = 640320
    D = 426880
    E = 10005
    C3_24  = C ** 3 // 24
    DIGITS_PER_TERM = math.log(53360 ** 3) / math.log(10)  #=> 14.181647462725476
    FILENAME = "pi.txt"

    def __init__(self, digits):
        """ Initialization
        :param int digits: digits of PI computation
        """
        self.digits = digits
        self.n      = math.floor(self.digits / self.DIGITS_PER_TERM + 1)
        self.prec   = math.floor((self.digits + 1) * math.log2(10))

    def compute(self):
        """ Computation """
        try:
            tm_s = time()
            p, q, t = self.__bsa(0, self.n)
            one_sq = mpz(10) ** (2 * self.digits)
            sqrt_c = isqrt(self.E * one_sq)
            pi = (q * self.D * sqrt_c) // t
            with open(self.FILENAME, "w") as f:
                f.write(str(pi))
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
                    p_ab = mpz((6 * a -5) * (2 * a - 1) * (6 * a - 1))
                    q_ab = mpz(a * a * a * self.C3_24)
                t_ab = p_ab * (self.A + self.B * a)
                if a & 1:
                    t_ab *= -1
            else:
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