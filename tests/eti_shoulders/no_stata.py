"""Exercise the real plugin ABI without Stata, including the comparison path."""
import ctypes as ct,sys,struct,csv
from pathlib import Path
import numpy as np
host,plugin,fixture=sys.argv[1:]
h=ct.CDLL(host)
h.host_open.argtypes=[ct.c_char_p];h.host_open.restype=ct.c_int
h.host_data.argtypes=[ct.c_int,ct.c_int];h.host_column.argtypes=[ct.c_int,ct.POINTER(ct.c_double)]
h.host_call.argtypes=[ct.c_int,ct.POINTER(ct.c_char_p)];h.host_macro.argtypes=[ct.c_char_p];h.host_macro.restype=ct.c_char_p
h.host_matrix.argtypes=[ct.c_int,ct.c_int]
assert h.host_open(plugin.encode())==0
z=(np.arange(400,dtype=float)-200)/256
h.host_data(len(z),1);h.host_column(0,z.ctypes.data_as(ct.POINTER(ct.c_double)))
def call(*args):
 a=(ct.c_char_p*len(args))(*(str(x).encode() for x in args));assert h.host_call(len(args),a)==0,args
def macro(k):return h.host_macro(('_'+k).encode()).decode()
def scalars():
 names=macro('scalar_names').split();h.host_matrix(1,len(names));call('fetch','scalars','out')
 return dict(zip(names,[struct.unpack('>d',bytes.fromhex(v))[0] for v in macro('bits').split()]))
for row in csv.DictReader(open(fixture)):
 call('npmle',.08,.03,0,0,1,4,25,1000,200,.01,10000,3,3,0,1e-9,1e-7,'_none',row['basis'])
 n=scalars();call('release')
 call('eti',.08,.03,0,0,1,12,45,1.34,1,200,n['latent_left'],n['latent_right'],n['lambda_at_grid_edge'],n['shoulders_weak'],n['shoulders_inconclusive'],row['basis'],'converged')
 e=scalars();call('release')
 for k in ('saez','continuous'):assert abs(e['npmle.'+k+'.e']-float(row[k]))<1e-8
 assert e['npmle.shoulder_left']==n['latent_left']
print('PASS: linear/cubic NPMLE and both ETI rules through no-Stata ABI')
