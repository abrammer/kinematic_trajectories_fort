 external EX02 "./test.so"
  
  x = (/2.,2.5,3.5,4/)
  fx = 1/x
  
  x1 = 3.
 print( "x  :  f(x)") 
  print(x+" : "+ fx)
 print("+++++++++++++++")
print("Matlab Interp") 
  print(x1+" : "+EX02::y_interp(x, fx, x1,dimsizes(x)))
print("Janiga Interp")
  print(x1+" : "+EX02::y2_interp(x,fx, x1))


 print("+++++++++++++++")
  x = (/1.,2.,3.,4./)
  fx = (/1.,4.,5.,8./)
  x1 = 2.7
 
  print( "x  :  f(x)")
  print(x+" : "+ fx)
  print("+++++++++++++++")
  print("Matlab Interp")
  print(x1+" : "+EX02::y_interp(x, fx, x1,dimsizes(x)))
  print("Janiga Interp")
  print(x1+" : "+EX02::y2_interp(x,fx, x1))

xx = 0.
 int = (EX02::linint(10., 100., 0.5, xx))
print(int)
print(xx)

