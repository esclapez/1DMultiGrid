EXEname := Driver.exe

fCOMP := gfortran -g -fbounds-check -O3 -ffree-line-length-512
XTRALIBS :=  -lm -lgfortran -lblas -llapack

f90_sources += mod_prec_defs.f90 mod_utils.f90 class_LinOpBC.f90 class_LinOpLvl.f90 class_ABecLvl.f90 class_ABecCecLvl.f90 class_LinOp.f90 class_ABecLap.f90 class_ABecCecLap.f90 class_mg1d.f90 mgDriver.f90

f90_objects := $(f90_sources:%.f90=%.o)

${EXEname}: ${f90_objects}
	${fCOMP} -o ${EXEname} ${f90_objects} ${XTRALIBS}

clean:
	\rm -rf ${EXEname} ${f90_objects} *.mod

%.o: %.f90
	${fCOMP} -c $^ -o $*.o
