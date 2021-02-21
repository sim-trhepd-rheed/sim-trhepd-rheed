# Calculation of Ge(001)-c4x2

The present sample presents the calculation of Ge(001)-c4x2 surface structure in the one-beam condition [1].

[1] Two-stage data-analysis method for total-reflection high-energy positron diffraction (TRHEPD),  Kazuyuki Tanaka, Izumi Mochizuki, Takashi Hanada, Ayahiko Ichimiya, Toshio Hyodo, Takeo Hoshi, JJAP Conf. Series, in press; Preprint:https://arxiv.org/abs/2002.12165 

The required input files are

bulk.txt

surf.txt

The output files are stored in the output/ directory.

## (1) bulk-part calculation

First, run the bulk-part calculation by typing

$ bulk.exe

Then the following output will appear 

0:electron 1:positron ?

P

input-filename (end=e) ? :

bulk.txt

output-filename :

bulkP.b

and the output file bulkP.b will be generated. 

## (2) surface-part calculation

Second, run the surface-part calculation by typing 

$ surf.ext

Then the following output will appear 

 bulk-filename (end=e) ? :

bulkP.b             

structure-filename (end=e) ? :

surf.txt            

output-filename :

surf-bulkP.s                      

and the output file surf-bulkP.s wil be generated.

## (3) Convolution procedure 

Third, the convolution procedure should be performed with the python3 script

script/make_convolution.py

in the present package. In general, one should may the source for the lines 

first_line = 5

last_line = 74

row_number = 2

omega = 0.5

The above values are valid for the present sample and one does not need to edit the source. 

$ python3 make_convolution.py 

Then the following output will appear 

len(C_list): 70

and the output file convolution.txt wil be generated and contains the rocking-curve data as

0.100000 0.002374995

0.200000 0.003614789

0.300000 0.005023215

.

.

.

6.800000 0.000155263

6.900000 0.000133134

7.000000 0.000101161

