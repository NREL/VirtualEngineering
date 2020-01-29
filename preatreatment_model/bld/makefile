SOURCEDIR=../src
FORT = gfortran
F2PY = f2py

PRETREATFILES=${SOURCEDIR}/globalvars.o ${SOURCEDIR}/exprparser.o ${SOURCEDIR}/polynomial.o\
	      ${SOURCEDIR}/quadrature.o ${SOURCEDIR}/inputs.o ${SOURCEDIR}/transport.o ${SOURCEDIR}/element.o\
	      ${SOURCEDIR}/outputs.o ${SOURCEDIR}/gmres.o ${SOURCEDIR}/timestepper.o
#FLAGS2= -g -fbounds-check -ffpe-trap=invalid,zero,overflow
PRETREATOBJFILES=exprparser.o polynomial.o inputs.o quadrature.o globalvars.o element.o transport.o outputs.o gmres.o timestepper.o
FLAGS2= -O3
%.o: %.f90
	${FORT} -c ${FFLAGS} $< -o $@

pretreat: ${PRETREATFILES}
	mv ${SOURCEDIR}/*.o .
	${FORT} -shared -fPIC -o libptreat.so *.o
	${FORT} -L. -lptreat ${PRETREATOBJFILES} ${SOURCEDIR}/main.f90 -o pretreat.exe

ptpython: ${PRETREATFILES}
	mv ${SOURCEDIR}/*.o .
	${FORT} -shared -fPIC -o libptreat.so *.o
	${F2PY} -c -L. -lptreat -m pt ${SOURCEDIR}/f2pymain.f90

clean:
	rm -rf *.o *.mod *.exe *.so pt.*
