import numpy as np
import cmath
import copy

N_POINT=1024
BIT_TOT=10
BIT_MAT=3

def binary_print(x):
    for i in x:
        print(bin(i)[2:].zfill(BIT_TOT),end=',')
    print()

def generate_signal(samples = 2048, fs = 44_100):
    f0, f1, f2, f3, f4 = 1_000, 4_000, 6_000, 8_000, 17_357
    ts = 1 / fs
    time = np.linspace(0, samples*ts, samples)

    tone_0 = 1 * np.sin(2*np.pi*f0*time)
    tone_1 = 0.8 * np.sin(2*np.pi*f1*time)
    tone_2 = 2 * np.sin(2*np.pi*f2*time)
    tone_3 = 0.1 * np.sin(2*np.pi*f3*time)
    tone_4 = 0.3* np.sin(2*np.pi*f4*time)
    signal = tone_0 + tone_1 + tone_2*1j + tone_3*1j + tone_4*1j
    return signal

def shuffle(x_in):
    x=copy.deepcopy(x_in)
    l=N_POINT
    for i in range(BIT_TOT-BIT_MAT):
        y=[]
        for sta in range(0,N_POINT,l):
            z=[[],[]]
            for j in range(sta,sta+l):
                z[(j-sta)%2].append(x[j])
            y+=z[0]+z[1]
        x=y
        l=l//2
    return x    
    
def get_omg_mat_8():
    mat=[]
    for i in range(8):
        # mat.append([])
        for j in range(8):
            v=cmath.exp(-2j * cmath.pi * j*i / 8)
            if np.fabs(v.real)<1e-14:
                real=0
            else:
                real=v.real
            if np.fabs(v.imag)<1e-14:
                image=0
            else:
                image=v.imag
            # mat[i].append((float(real),float(image)))
            mat.append((int(real*32),int(image*32)))
            # mat[i].append(v)
            # mat.append(float(real))
            # mat.append(float(image))
    return mat

def get_swap_cycle(dst):
    vis=[0 for i in range(N_POINT)]
    res=[]
    for i in range(N_POINT):
        if vis[i]==1: continue
        sta=i
        f=i
        t=-1
        c=[]
        while t!=sta:
            c.append(f)
            vis[f]=1
            t=dst[f]
            f=dst[f]
        res.append(c)
    return res

def shuffle_by_swap_cycle(x,sc):
    x=copy.deepcopy(x)
    for c in sc:
        # print([x[i] for i in c])
        if len(c)==1: continue
        tmp=x[c[-1]]
        x[c[-1]]=x[c[0]]
        for i in range(1,len(c)-1):
            x[c[i-1]]=x[c[i]]
        x[c[-2]]=tmp
        # print([x[i] for i in c])
        # print('-'*32)
    return x


def print_cl():
    x=generate_signal(N_POINT, 44_100)
    # x=[i for i in range(N_POINT)]
    s=[i for i in range(N_POINT)]
    d=shuffle(s)
    sc=get_swap_cycle(d)

    cl=[]
    out_sc=[]
    for c in sc:
        cl.append(len(c))
        out_sc+=c
    print(cl)
    # print(out_sc)
    # print(len(sc))

    # res=shuffle_by_swap_cycle(x,sc)
    # ans=shuffle(x)
    # print(res[:5])
    # print(ans[:5])

    # out_sc={2:[],4:[]}
    # for c in sc:
    #     if len(c)==1: continue
    #     out_sc[len(c)]+=c
    # print(len(out_sc[2]))
    # print(out_sc[2])
    # print(len(out_sc[4]))
    # print(out_sc[4])


x=get_omg_mat_8()
# for line in x:
#     for r,i in line:
#         if np.fabs(r)<1e-14:
#             r=0
#         if np.fabs(i)<1e-14:
#             i=0
#         print(r,i)
print(x)



