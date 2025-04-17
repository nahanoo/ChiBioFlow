from sympy import preview
from os.path import join

def write_e(expr,f):
    preview(expr,viewer='file',filename=join('equations',f),dvioptions=['-D','600','-bg','Transparent'])

J_Ct = '$J_{Ct} = \mu_{Ct}\frac{C}{K_C + C}$'
write_e(J_Ct,'J_Ct.svg')
