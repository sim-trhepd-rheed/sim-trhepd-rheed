import numpy as np
import argparse

parser = argparse.ArgumentParser(description='test')
parser.add_argument('--filename', default='surf-bulkP.s')
parser.add_argument('--omega',    default=0.5)
parser.add_argument('--verbose',  action="store_true")
args = parser.parse_args()

if args.verbose:
   verbose_mode = True
else:
   verbose_mode = False

input_filename = args.filename
omega          = float(args.omega)

if verbose_mode:
   print('arg:filename=',args.filename)
   print('omega   =',args.omega)

sigma = 0.5 * omega / (np.sqrt(2.0*np.log(2.0)))

#def g(x):
#    g = (0.939437 / omega) * np.exp(-2.77259 * (x ** 2.0 / omega ** 2.0))
#    return g

def g(x):
    g = (1.0 / (sigma*np.sqrt(2.0*np.pi))) * np.exp(-0.5 * x ** 2.0 / sigma ** 2.0)
    return g

file_input = open(input_filename, "r")
Clines = file_input.readlines()
file_input.close()

number_of_lines        = int(len(Clines))
number_of_header_lines = 4
if verbose_mode:
   print("number_of_lines         =", number_of_lines)
   print("number_of_header_lines  =", number_of_header_lines)

#degree_list = []
#C_list = []

# Read the file header 
line = Clines[0]
#print("file header :", line, end="")
line = Clines[1]
#print("file header :", line, end="")
line = line.replace(",", "")
data = line.split()
#print(data)
number_of_azimuth_angles  = int(data[0])
number_of_glancing_angles = int(data[1])
number_of_beams           = int(data[2])

if verbose_mode:
   print("number of azimuth angles  = ", number_of_azimuth_angles)
   print("number of glancing angles = ", number_of_glancing_angles)
   print("number of beams           = ", number_of_beams)

# Define the array for the rocking curve data.
#   Note the components with (beam-index)=0 are the degree data
RC_data_org = np.zeros((number_of_glancing_angles, number_of_beams+1))
RC_data_cnv = np.zeros((number_of_glancing_angles, number_of_beams+1))

# Read the file header 
line = Clines[2]
#print("file header :", line, end="")
line = Clines[3]
#print("file header :", line, end="")
line = line.replace(",", "")
data = line.split()
#print(data)

if verbose_mode:
   print("beam index (p,q): i  p_i q_i")
   for beam_index in range(number_of_beams):
      print(beam_index, data[beam_index*2], data[beam_index*2+1])

for g_angle_index in range(number_of_glancing_angles):
    line_index = number_of_header_lines + g_angle_index
    line = Clines[ line_index ]
#   print("data line: ", line_index, g_angle_index, line)
    line = line.replace(",", "")
    data = line.split()
#   print(data)
    RC_data_org[g_angle_index,0]=float(data[0])
    RC_data_cnv[g_angle_index,0]=float(data[0])
    for beam_index in range(number_of_beams):
        RC_data_org[g_angle_index, beam_index+1] = data[beam_index+1]

angle_interval = RC_data_org[1,0] - RC_data_org[0,0]

if verbose_mode:
   print('angle_ interval=', angle_interval)

for beam_index in range(number_of_beams):
    for index in range(number_of_glancing_angles):
        integral = 0.0
        for index2 in range(number_of_glancing_angles):
            integral += RC_data_org[index2,beam_index+1] * g(RC_data_org[index,0] - RC_data_org[index2,0]) * angle_interval
        RC_data_cnv[index, beam_index+1]=integral

#np.savetxt("original.txt", RC_data_org, fmt="%.5e")
np.savetxt("convolution.txt", RC_data_cnv, fmt="%.5e")


