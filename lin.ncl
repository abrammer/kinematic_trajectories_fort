external EX01 "./lint.so"
 xx=0.0
 y = EX01::abcd(10., 20., 0.5, xx)
 print(y)


delete(y)
 xx=0.0
 y = EX01::lbcd(10., 20., 0.5, xx)
print(y)

