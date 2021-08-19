CXX = g++
FLAGS = -O3

EXE = run_two_phase_batch_example run_batch_lignocell_EH_standalone

all: ${EXE}
	echo "everything compiled"

run_two_phase_batch_example: two_phase_batch_model.cpp
	${CXX} ${FLAGS} two_phase_batch_model.cpp -o run_two_phase_batch_model

run_batch_lignocell_EH_standalone: batch_lignocell_EH_standalone.cpp
	${CXX} ${FLAGS} batch_lignocell_EH_standalone.cpp -o run_batch_lignocell_EH_standalone

clean:
	rm -f ${EXE}
