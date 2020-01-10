import numpy as np
import sys
import time

print('Hello')
if len(sys.argv) > 1:
	input_val = sys.argv[1]
	try:
		input_val = float(input_val)
		print('Input Number =', input_val)
		print('Double Input Number =', 2*input_val)
	except:
		print('Input should be a number', input_val)


N = 5

for k in range(5):
	print('Script waiting %d seconds...' % (N-k))
	time.sleep(1)

print('Script finished!')