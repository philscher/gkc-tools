Descritpion
________________

A code for solving the reduced MHD equations in two dimensional sheared slab
geometry. Parallelized using OpenMP.



Notes
__________________

Usually stack size is increased, however, as we use large static arrays, with a compile
time dependent stack frame (?). Thus enable (-heap-array) to declare automatic arrays
on the heat and not on the stack.

Use icc 13.0 or later
