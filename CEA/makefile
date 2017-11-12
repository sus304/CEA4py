all:
ifeq (\$(OS), Windows_NT)
	# for Windows
	gfortran -o FCEA2.exe cea2.f
	gfortran -o b1b2b3.exe b1b2b3.f
	gfortran -o syntax.exe syntax.f
else
	# UNAME = \${shell uname}
	# ifeq (\(UNAME), Linux)
		# for Linux
		gfortran -o FCEA2 cea2.f
		gfortran -o b1b2b3 b1b2b3.f
		gfortran -o syntax syntax.f
	# endif

endif