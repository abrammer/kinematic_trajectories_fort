 external EX02 "./test.so"
  
  x = (/2.,2.5,3.5,4/)
  fx = 1/x
  
  x1 = 3.
 print( "x  :  f(x)") 
  print(x+" : "+ fx)
 print("+++++++++++++++")
print("Nevillie Interp") 
  print(x1+" : "+EX02::neville_interp(x, fx, x1,dimsizes(x)))


 print("+++++++++++++++")
  x = (/1.,2.,3.,4./)
  fx = (/1.,4.,5.,8./)
  x1 = 2.7
 
  print( "x  :  f(x)")
  print(x+" : "+ fx)
  print("+++++++++++++++")
  print("Nevillie Interp")
  print(x1+" : "+EX02::neville_interp(x, fx, x1,dimsizes(x)))

  print("+++++++++++++++")
print("Linear interp")
print("0.52 between 10 : 100")
print( EX02::linear_interp(10., 100., 0.52))

