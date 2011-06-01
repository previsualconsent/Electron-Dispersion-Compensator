#Makfe file for EDC
tuning:tuning.for
	fl32 /Op /4Yb /4Yf /nologo tuning.for
	tuning.exe
	
path:path.for
	fl32 /Op /4Yb /4Yf /nologo path.for
	path.exe

bgradient.tuning.3d:bgradient.tuning.3d.for
	fl32 /Op /4Yb /4Yf /nologo bgradient.tuning.3d.for
	bgradient.tuning.3d.exe
	
EvsPOS:EvsPOS.for
	fl32 /Op /4Yb /4Yf /nologo EvsPOS.for
	EvsPOS.exe

mBvsdx:mBvsdx.for
	fl32 /Op /4Yb /4Yf /nologo mBvsdx.for
	mBvsdx.exe

foo:foo.for
	fl32 /Op /4Yb /4Yf /nologo foo.for
	foo.exe

