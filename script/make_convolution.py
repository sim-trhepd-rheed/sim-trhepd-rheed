import numpy as np

first_line = 5
last_line = 60
row_number = 8

omega = 0.5
def g(x):
    g = (0.939437 / omega) * np.exp(-2.77259 * (x ** 2.0 / omega ** 2.0))
    return g

file_input = open("surf-bulkP.s", "r")
Clines = file_input.readlines()
file_input.close()

degree_list = []
C_list = []
for line in Clines[first_line -1:last_line]:
    line = line.replace(",", "")
    data = line.split()
    degree_list.append(float(data[0]))
    C_list.append(float(data[row_number - 1]))

print("len(C_list):", len(C_list))

new_C_list = []
for index in range(len(C_list)):
    integral = 0.0
    for index2 in range(len(C_list)):
        integral += C_list[index2] * g(degree_list[index] - degree_list[index2]) * 0.1
    new_C_list.append(integral)

C_list = new_C_list

file_output = open("convolution.txt", "w")
for i in range(len(C_list)):
    file_output.write("%f %.9f\n"%(degree_list[i], C_list[i]))
file_output.close()
