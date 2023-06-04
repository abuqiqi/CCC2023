from fft_window import *
import json

#### VMC validators ####

# Note: For stream-based Window Fn, the low level IP multipies
# SSR*2 to get input and output port number. But for VMC
# users, SSR = no. of in or out ports. Hence we divide it by
# 2 here before passing on to low level IP.

def vmc_validate_coeff_type(args):
  data_type = args["data_type"]
  coeff_type = args["coeff_type"]
  return fn_validate_coeff_type(data_type, coeff_type)

def vmc_validate_point_size(args):
  point_size = args["point_size"]
  data_type = args["data_type"]
  return fn_validate_point_size(point_size, data_type)

def vmc_validate_shift_val(args):
  data_type = args["data_type"]
  shift_val = args["shift_val"]
  return fn_validate_shift(data_type, shift_val)

def vmc_validate_input_window_size(args):
  point_size = args["point_size"]
  input_window_size = args["input_window_size"]
  return fn_validate_window_vsize(point_size, input_window_size)

def vmc_validate_ssr(args):
  data_type = args["data_type"]
  point_size = args["point_size"]
  interface_type = 1
  ssr = args["ssr"]//2
  return fn_validate_ssr(data_type, point_size, interface_type, ssr)
  
def vmc_validate_is_dyn_pt_size(args):
  point_size = args["point_size"]
  ssr = args["ssr"]//2
  dyn_pt = 1 if args["is_dyn_pt_size"] else 0
  return fn_validate_dyn_pt_size(point_size, ssr, dyn_pt)

def vmc_validate_coeff(args):
  dyn_pt = 1 if args["is_dyn_pt_size"] else 0
  return fn_validate_weights(args["point_size"], dyn_pt, args["coeff"])

#### VMC graph generator ####
def vmc_generate_graph(name, args):
  tmpargs = {}
  tmpargs["TT_DATA"] = args["data_type"]
  tmpargs["TT_COEFF"] = args["coeff_type"]
  tmpargs["TP_POINT_SIZE"] = args["point_size"]
  tmpargs["TP_WINDOW_VSIZE"] = args["input_window_size"]
  tmpargs["TP_SHIFT"] = args["shift_val"]
  tmpargs["TP_API"] = 1
  tmpargs["TP_SSR"] = args["ssr"]
  tmpargs["TP_DYN_PT_SIZE"] = 1 if args["is_dyn_pt_size"] else 0
  tmpargs["weights"] = args["coeff"]
  return generate_graph(name, tmpargs)
