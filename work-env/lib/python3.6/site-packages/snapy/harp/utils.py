from numpy import array
from ..athena.athena_read import athinput
import re

def get_rt_bands(inpfile):
  inp = athinput(inpfile)

  # count number of bands
  num_bands = 0
  for key in inp['radiation'].keys():
    if re.match('b[0-9]+$', key):
      num_bands += 1

  # read frequency
  band_info = []
  for i in range(num_bands):
    band_info.append(list(map(float, inp['radiation']['b%d' % (i+1)].split())))
  band_info = array(band_info)

  return band_info

def get_ray_out(inpfile):
  inp = athinput(inpfile)

  # out directions
  outdir = inp['radiation']['outdir'].split(' ')
  num_dirs = len(outdir)

  # polar and azimuthal angle
  amu, aphi = [], []
  for i in range(num_dirs):
    m = re.search('\((.*),(.*)\)', outdir[i])
    if m.group(1) != '':
      amu.append(float(m.group(1)))
    else:
      amu.append(0.)
    if m.group(2) != '':
      aphi.append(float(m.group(2)))
    else:
      aphi.append(0.)
  return array(amu), array(aphi)

def get_sample_pressure(inpfile):
  inp = athinput(inpfile)
  return list(map(float, inp['inversion']['PrSample'].split()))
