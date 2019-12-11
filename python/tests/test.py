from qif import channel
import numpy as np
from fractions import Fraction
import sys

print("running tests")


def f(i):
	return i

np.set_printoptions(formatter={'object':lambda x: str(x)})

q = Fraction(1,4)
# qq = channel.rat_test(q)
# print(qq)
# print(type(q))
# print(type(qq))
# print(q.denominator)


A = [[0,1,2],[3,4,5],[6,7,8.]];

pi = np.array([.25,.25,.25,.25])
pir = np.array([q,q,q,q])
Ad = np.array([pi, pi, pi, pi])
Ar = np.array([[q,q,q,q], [q,q,2*q,0*q]]);


d = channel.draw(Ar, pir, 5);
print(d)
print(type(d[0,0]));
# print(channel.draw(Ar, pir)[0]);


# print(channel.identity(4))
sys.exit()

C = np.array(A)
Cr = np.array(Ar)
print(channel.is_proper(C))
print(channel.assert_proper(Cr))
print(channel.normalize(C))
print(channel.normalize(Cr))
sys.exit()



Ar = np.array([[q+1,q+2], [q+3,q+4]]);
# print(Ar)
# print(Ar.dtype)

Ad = np.array(A, order='F')
print(Ad.flags)

print(channel.poly(Ad))
print(channel.poly(Ar))



print(channel.ff(4));
sys.exit()



def euclid(a,b):
	return abs(a-b)
	# return np.ones((3,3), order='F')

f = channel.metric(euclid)
print(f(35,120))



# B = [[0,1,2],[3,4,5],[6,Fraction(2,3),8]];
# An = np.array(A).astype('object') + Fraction()
# print(An*An)
# print(np.array(B).dtype)
# print(np.array(B))
# print(type(An.dtype))
# print(type(np.object))

# sys.exit()



a = np.ones((3,3), order='F')
a *= A;
# a[0,1] = 2;
pointer, read_only_flag = a.__array_interface__['data']
print("a", "0%x" % pointer, read_only_flag)
print(a.flags)
# while 1:
channel.test_mat(a)

# print(b)
# pointer, read_only_flag = b.__array_interface__['data']
# print("b", "0%x" % pointer, read_only_flag)

sys.exit()



x = channel.ret_mat()
print(x.flags)
x[0,0] = 9
pointer, read_only_flag = x.__array_interface__['data']
print("0%x" % pointer, read_only_flag)

# y = channel.ret_eigen()
# pointer, read_only_flag = y.__array_interface__['data']
# print("0%x" % pointer, read_only_flag)

# print(y)
# print(channel.sum(x))
print(x)