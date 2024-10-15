import pt
import numpy as np
import sys

solnvec = pt.main(sys.argv[1], sys.argv[2], sys.argv[3])
np.savetxt('pt_solnvec.csv', solnvec, delimiter=',')
